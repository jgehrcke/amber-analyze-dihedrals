#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 Jan-Philip Gehrcke, TU Dresden
# http://gehrcke.de
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import sys
import StringIO
import argparse
import itertools
from fnmatch import fnmatchcase
from matplotlib import pyplot
import pandas as pd
import numpy as np
import logging
import brewer2mpl


logging.basicConfig(
    format='%(asctime)s,%(msecs)-6.1f %(levelname)s: %(message)s',
    datefmt='%H:%M:%S')
log = logging.getLogger()
log.setLevel(logging.DEBUG)


# Reserve name for global options object.
options = None
# Reserve name for global open window option, which is interpreted at the
# end of the runtime of the program, potentially instructing matplotlib to
# open/display all figures in their windows.
open_figure_windows = False

def main():
    global options
    parser = argparse.ArgumentParser(
        description=("Program for evaluating torsion angle data sets as "
            "created by cpptraj (from AmberTools 13)."))
    parser.add_argument('dihedraldatafile', action="store",
        help=("Path to data file as produced by cpptraj's 'dihedral' command. "
            "Might contain multiple columns (cpptraj 'datasets'). Column "
            "names (headers in the first line) are interpreted as angle data "
            "set names."
            ))
    parser.add_argument('-t', '--two-dimensional', action="append", nargs='+',
        metavar=("'nameX,nameY'", "'title'"),
        help=("Comma-separated list of two angle data set names to be plotted "
            "in a two-dimensional heat map histogram. Names can be used from "
            "original "
            "data or from merged data sets. Data corresponding to the first "
            "name will be placed on the horizontal axis. If a data set name "
            "contains one of the strings ['psi', 'phi'], the axis label will "
            "be set to 'greek_letter / degrees' correspondingly. A second "
            "argument can be provided to the option containing a descriptive "
            "title for the entire plot. Other arguments to this option are "
            "ignored. This option can be specified multiple times in order to "
            "plot various data set pairs."))
    parser.add_argument('--pdf', action="store_true", default=False,
        help=("Instead of displaying a figure, write a PDF file. The file is "
            "written to the current working directory. Its name is prepended "
            "with a prefix, if provided by the --imagefile-prefix option "
            "(this prefix might contain path separators). The rest of the "
            "name is generated from the name(s) of the data set(s) involved "
            "in the figure."))
    parser.add_argument('--png', action="store_true", default=False,
        help=("Instead of displaying a figure, write a PNG file. See --pdf "
            "for further details. When this option is specified, the "
            "--resolution option takes effect."))
    parser.add_argument('-i', '--imagefile-prefix', action="store", default='',
        metavar='filename-prefix',
        help=("A common prefix that is prepended to all image "
            "filenames before writing. Takes effect when --png and/or --pdf "
            "is specified."))
    parser.add_argument('-r', '--resolution', action="store",
        type=int, default=150, metavar="DPI",
        help=("Image resolution (DPI) when saving plot as PNG. Forwarded to "
            "matplotlib's 'dpi' argument of the 'savefig()'' function. "
            "Default: 150."))
    parser.add_argument('-m', '--merge', nargs=2, action="append",
        metavar=("name", "'wildcard'"),
        help=("This option consumes the following two arguments. The first "
            "defines the name of the merged data set. The second must be a "
            "Unix shell-style wildcard that will be matched against original "
            "data set names. The values "
            "of an original data set are added to the merged data set (of "
            "given name) if the name of the original data set matches the "
            "wildcard. The option can be "
            "used multiple times. Multiple occurrences of this option "
            "can share the merged data set name (leading to a single merged "
            "data set based on multiple wildcards). Multiple occurrences "
            "of this option must not share the same wildcard. If an "
            "angle name matches multiple wildcards, it is matched with "
            "the first one as specified on the command line. Wildcards "
            "are case-sensitive. Example: "
            "--merge typeA 'ARG*LYS*psi' --merge typeA 'ARG?CYX?psi'. This "
            "creates one merged data set with name 'typeA' and adds data of "
            "original angle data sets whose names match either the first or "
            "the second wildcard. Make sure to quote wildcards. Reminder: "
            "In a wildcard, * matches everything while ? matches a single "
            "character. [ABC] matches A,B, or C."
            ""
            ))
    parser.add_argument('-b', '--bins', action="store", type=int, default=20,
        help="Number of histogram bins for each dimension. Default: 20.")
    options = parser.parse_args()

    if options.merge:
        # Validate user-given input regarding data merging, provide useful
        # output.
        log.info("Merge option was specified %s time(s). Validate.",
            len(options.merge))
        merge_names, merge_wildcards = zip(*options.merge)
        # Make sure there are no wildcard repetitions.
        if not len(options.merge) == len(set(merge_wildcards)):
            sys.exit("Repetitions of the same merge wildcard are not allowed.")
        log.info("Grouping merge option(s) by name.")
        merge_groups = {n:[] for n in merge_names}
        for name, wc in zip(merge_names, merge_wildcards):
            merge_groups[name].append(wc)
        for n, wcs in merge_groups.iteritems():
            log.info("Merge group:\n  name: %r\n  wildcard(s): %s)", n, wcs)

    twodim_hist_infos = []
    for twodim_hist_info_list in options.two_dimensional:
        # Validate user-given input regarding 2D histogram plotting.
        log.info("2D plotting option was specified. Validate.")
        # One argument to this option is required. It must contains a comma-
        # separated list of two data set names.
        x_y_names_for_2d_plot = twodim_hist_info_list[0].split(',')
        if not len(x_y_names_for_2d_plot) == 2:
            sys.exit(("Exactly two comma-separated data set names must be "
                "provided as an argument to the 2D plotting option. "
                "%s found." % len(x_y_names_for_2d_plot)))
        # One more argument is optional, if given it contains the plot title.
        title_for_2d_plot = None
        if len(twodim_hist_info_list) > 1:
            title_for_2d_plot = twodim_hist_info_list[1]
        twodim_hist_infos.append({
            'title': title_for_2d_plot,
            'x_y_names': x_y_names_for_2d_plot
        })


    # Read raw data to pandas DataFrame.
    original_df = parse_dihed_datafile()

    # Merge data if applicable.
    merged_series = None
    if options.merge:
        merged_series = merge_dataseries_by_wildcards(
            original_df, merge_groups)

    # Plot 2D histogram if applicable.
    for twodim_hist_info in twodim_hist_infos:
        histogram_from_dataset_names(
            twodim_hist_info['x_y_names'],
            twodim_hist_info['title'],
            original_df,
            merged_series)

    log.info("Data processing and plotting finished.")
    if open_figure_windows:
        log.info(("There are figures ready for being displayed. Calling "
            "pyplot.show()."))
        pyplot.show()


def histogram_from_dataset_names(
        dataset_names,
        title,
        original_df,
        merged_series):
    """Currently, `dataset_names` must be of length 1 (1D hist) or 2 (2D hist).
    The data sets are first searched for in the list of merged data series,
    then in the original dataframe. If two names are provided and found (in one
    or the other dataframe), the datasets are validated to be of the same
    length before considering them for plotting. If two data set names are
    provided, the first data set ends up on the x-axis (horizontal axis) of the
    2D histogram.
    """
    global open_figure_windows
    log.info("Instructed to plot histogram from data set(s) with name(s) %s.",
        dataset_names)
    def get_series(name):
        if name in merged_series:
            log.info("Found set '%s' in merged data set. Use it.", name)
            return merged_series[name]
        if name in original_df:
            log.info("Found set '%s' in original data set. Use it.", name)
            return original_df[name]
        log.error(("Dataset name '%s' was specified for plotting but does not "
            "exist in merged or original data.", name))
        sys.exit(1)

    if len(dataset_names) == 1:
        raise NotImplementedError("1D histogram is yet to be implemented.")
    elif len(dataset_names) == 2:
        series_x, series_y = [get_series(n) for n in dataset_names]
        log.info(("Angle data set names to be used for 2D histogram: '%s' "
            "(horizontal axis (x)) and '%s' (vertical axis, (y))." % (
                series_x.name, series_y.name)))
        log.info("Length of x data set: %s", len(series_x))
        log.info("Length of y data set: %s", len(series_y))
        if not len(series_x) == len(series_y):
            log.error(("I don't want to build a 2D histogram from two data "
                "sets different in length."))
            sys.exit(1)
        log.info("Plot it (using matplotlib).")
        if title is None:
            # Build default title from data set names.
            t = "dihedral '%s' vs. dihedral '%s'" % (
                series_x.name, series_y.name)
        else:
            t = title
        filename_wo_ext = None
        if options.pdf or options.png:
            log.info("Don't open plot in window (PNG/PDF output specified).")
            log.info(("Building image filename w/o extension, using prefix "
                "'%s', and both data set names."), options.imagefile_prefix)
            filename_wo_ext = "%s%s-%s" % (
                options.imagefile_prefix, series_x.name, series_y.name)
            log.info("Filename w/o extension: '%s'", filename_wo_ext)
        else:
            open_figure_windows = True
            log.info(("Planning to open plot in window, since no image file "
                "output has been specified."))
        create_2d_hist(
            series_x=series_x,
            series_y=series_y,
            title=t,
            xlabel="%s / degrees" % util_greek_map(series_x.name),
            ylabel="%s / degrees" % util_greek_map(series_y.name),
            save_png=options.png,
            save_pdf=options.pdf,
            filename_wo_ext=filename_wo_ext,
            resolution=options.resolution,
            axis_range=[[-180, 180], [-180, 180]] # Default for the moment.
            )
    else:
        raise Exception("Don't you dare asking for a higher dimension.")


def util_greek_map(a):
    """If string `a` contains greek letter from `translate` mapping, replace
    `a` with the greek representation. Otherwise return `a` unchanged.
    """
    translate = {
        'phi': r'$\phi$',
        'psi': r'$\psi$'
        }
    for greek_letter in translate:
        if greek_letter in a:
            return translate[greek_letter]
    return a


def create_2d_hist(
        series_x,
        series_y,
        title,
        xlabel,
        ylabel,
        save_png,
        save_pdf,
        filename_wo_ext,
        resolution,
        axis_range=None):
    """Create a 2D histogram (heat map) from data in the two DataSeries objects
    `series_x` and `series_y`.
    """
    log.info("Creating new figure.")
    fig = pyplot.figure()
    log.info("Calling 'hist2d', using %s bins.", options.bins)
    pyplot.hist2d(
        series_x.values,
        series_y.values,
        bins=options.bins,
        range=axis_range,
        cmap=brewer2mpl.get_map('Greys', 'sequential', 9).mpl_colormap)
    pyplot.title(title)
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    pyplot.colorbar()
    if save_pdf:
        fn = "%s.pdf" % filename_wo_ext
        log.info("(Over)writing '%s'.", fn)
        pyplot.savefig(fn)
    if save_png:
        fn = "%s.png" % filename_wo_ext
        log.info("(Over)writing '%s'.", fn)
        pyplot.savefig(fn, dpi=resolution)
    if save_png or save_pdf:
        # Don't display this figure in any window later on.
        pyplot.close(fig)


def merge_dataseries_by_wildcards(df, merge_groups):
    """Create and return dictionary containing only merged data sets. Modify
    original dataframe `df` on the fly (remove those series (columns) that)
    have been merged.
    """
    # Create dictionary containing the names of the merge groups as keys.
    # For each key, build a list of matching pandas DataSeries as value. A
    # data series is added to a key when the wildcard corresponding to the
    # merge group matches the original data set (column) name.
    log.info("Merging data according to command line options provided.")
    merge_groups_data_series = {name:[] for name in merge_groups}
    for name, wildcards in merge_groups.iteritems():
        for original_name in df.columns:
            for wc in wildcards:
                if fnmatchcase(original_name, wc):
                    log.info("Column '%s' matches merge wildcard '%s'.",
                        original_name, wc)
                    log.debug("Append DataSeries to merge group '%s'.", name)
                    merge_groups_data_series[name].append(df[original_name])
                    log.debug("Delete column from original DataFrame.")
                    del df[original_name]

    log.info("Groups identified:")
    for name, series_list in merge_groups_data_series.iteritems():
        log.info("  Name '%s': %s columns", name, len(series_list))
        for i, s in enumerate(series_list):
            log.info("        Column %s: '%s' with %s values.",
                i, s.name, len(s))

    log.info("Merging data within groups.")
    merged_series_list = []
    for name, series_list in merge_groups_data_series.iteritems():
        merged_series = pd.concat(series_list)
        log.debug("Type of merged_series: %s", type(merged_series))
        merged_series.name = name
        log.info("Created series with name '%s' containing %s values.",
            merged_series.name, len(merged_series))
        merged_series_list.append(merged_series)
    log.debug("Creating dictionary from multiple DataSeries.")
    merged_series =  {s.name:s for s in merged_series_list}
    #log.info("Created DataFrame from merged series. Details:\n%s",
    #    df_merged_series)
    log.info("Created dictionary from merged series. Details:")
    for n, s in merged_series.iteritems():
        log.info("    series name: '%s', series length: %s", n, len(s))
    #log.info("Head:\n%s", df_merged_series.head())
    return merged_series


def parse_dihed_datafile():
    log.info("Reading '%s' to pandas DataFrame.", options.dihedraldatafile)
    # Transform output of cpptraj into classical CSV shape, use an in-memory
    # buffer for this.
    with open(options.dihedraldatafile) as f:
        lines = f.readlines()
    # Remove leading '#' in first line, strip leading and trailing white spaces
    # and replace whitespace delimiters with commas.
    firstline = ','.join(lines[0].strip().strip('#').split())
    otherlines_gen = (','.join(l.strip().split()) for l in lines[1:])
    csv_buffer = StringIO.StringIO()
    csv_buffer.write("%s\n" % firstline)
    csv_buffer.write("\n".join(otherlines_gen))
    csv_buffer.seek(0)
    df = pd.read_csv(csv_buffer)
    log.debug("Columns read:\n%s", "\n".join(c for c in df.columns))
    log.debug("Dataframe head:\n%s", df.head())
    log.info("Built pandas DataFrame with %s columns and %s rows.",
        len(df.columns), len(df))
    return df


if __name__ == "__main__":
    main()
