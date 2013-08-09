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

import StringIO
import argparse
import itertools
from matplotlib import pyplot
import pandas as pd
import numpy as np
import logging


logging.basicConfig(
    format='%(asctime)s,%(msecs)-6.1f %(levelname)s: %(message)s',
    datefmt='%H:%M:%S')
log = logging.getLogger()
log.setLevel(logging.DEBUG)


# Reserve name for global options object.
options = None

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
    parser.add_argument('-t', '--two-dimensional', action="store",
        metavar="name1,name2",
        help=("Comma-separated list of two angle data set names to be plotted "
            "in a two-dimensional heat map. Names can be used from original "
            "data or from merged data sets."))
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
            "character. Docs: http://docs.python.org/2/library/fnmatch.html"
            ""
            ))




    parser.add_argument('-b', '--bins', action="store", type=int, default=20,
        help="Number of histogram bins for each dimension. Default: 20.")
    options = parser.parse_args()
    df = parse_dihed_datafile()

    df_merged_series = None
    if options.merge:
        suffixes_for_merge = options.merge.split(',')
        log.info("Angle name suffixes for merge:\n%s",
            "\n".join(suffixes_for_merge))
        df_merged_series = merge_dataframe_by_suffix(df, suffixes_for_merge)

    if options.two_dimensional:
        sfx2d_1, sfx2d_2 = options.two_dimensional.split(',')
        log.info(("Angle name suffixes for 2D histogram for matching angle "
            "name prefixes:\n%s"), "\n".join([sfx2d_1, sfx2d_2]))
        if df_merged_series:
            if all(s in df_merged_series for s in [sfx2d_1, sfx2d_2]):
                log.info("Found 2D plot suffixes (%s,%s) in merged data set.",
                    sfx2d_1, sfx2d_2)
                create_2d_hist(df_merged_series, sfx2d_1, sfx2d_2)

def util_greek_map(a):
    """If string `a` is in `translate` mapping, replace it with greek
    representation.
    """
    translate = {
        'phi': r'$\phi$',
        'psi': r'$\psi$'
        }
    if a in translate:
        return translate[a]
    return a

def create_2d_hist(df, cname_x, cname_y):
    """Create a 2D histogram (heat map) from data in DataFrame `df` (columns
    `chead1` and `chead2`).
    """
    df_cnames = ",".join("%r" % c for c in df.columns)
    log.debug("Creating 2D histogram from DataFrame (columns '%s')",
        df_cnames)
    log.debug("Horizontal axis: '%s', vertical axis: '%s'", cname_x, cname_y)
    pyplot.hist2d(df[cname_x].values, df[cname_y].values, bins=options.bins)
    pyplot.xlabel("%s / degrees" % util_greek_map(cname_x))
    pyplot.ylabel("%s / degrees" % util_greek_map(cname_y))
    pyplot.colorbar()
    pyplot.show()


def merge_dataframe_by_suffix(df, suffixes_for_merge):
    column_groups = {sfx:[] for sfx in suffixes_for_merge}
    for sfx in suffixes_for_merge:
        for c in df.columns:
            if c.endswith(sfx):
                log.info("Column '%s' matches merge suffix '%s'.", c, sfx)
                log.debug("Append data series (column) to groups for merge.")
                column_groups[sfx].append(df[c])
                log.debug("Delete column from original dataframe.")
                del df[c]

    log.info("Groups identified:")
    for sfx, series_list in column_groups.iteritems():
        log.info("  suffix '%s': %s columns" , sfx, len(series_list))
        for i, s in enumerate(series_list):
            log.info("        Column %s: '%s' with %s values.",
                i, s.name, len(s))
        #    print " type: %s" % type(s)
        #    print s.head()

    log.info("Merging data within groups.")

    merged_series_list = []
    for sfx, series_list in column_groups.iteritems():
        merged_series = pd.concat(series_list)
        log.debug("Type of merged_series: %s", type(merged_series))
        merged_series.name = sfx
        log.info("Created series with name '%s' containing %s values.",
            merged_series.name, len(merged_series))
        merged_series_list.append(merged_series)
    df_merged_series = pd.DataFrame(
        {s.name:s.values for s in merged_series_list})
    log.info("Created dataframe from merges series. Details:\n%s",
        df_merged_series)
    log.info("Head:\n%s", df_merged_series.head())
    return df_merged_series




def parse_dihed_datafile():
    log.info("Reading '%s' to pandas dataframe.", options.dihedraldatafile)
    # Bring output of cpptraj into classical CSV shape.
    with open(options.dihedraldatafile) as f:
        lines = f.readlines()
    # Remove leading '#' in first line, strip white spaces and replace
    # whitespace delimiter with comma.
    firstline = ','.join(lines[0].strip().strip('#').split())
    otherlines = (','.join(l.strip().split()) for l in lines[1:])
    csv_buffer = StringIO.StringIO()
    csv_buffer.write("%s\n" % firstline)
    csv_buffer.write("\n".join(otherlines))
    csv_buffer.seek(0)
    df = pd.read_csv(csv_buffer)
    log.debug("Columns read:\n%s", "\n".join(c for c in df.columns))
    log.debug("Dataframe head:\n%s", df.head())
    return df











if __name__ == "__main__":
    main()


