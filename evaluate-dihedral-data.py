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

import os
import sys
import argparse
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
    parser = argparse.ArgumentParser(description='Dickes Tool.')
    parser.add_argument('dihedraldatafile', action="store",
        help="Datafile as produces by the 'dihedral' command of cpptraj.")
    parser.add_argument('-t', '--two-dimensional', action="store",
        help=("Comma-separated list of two angle name suffixes to be plotted "
            "in a two-dimensional heat map if prefix matches."))
    parser.add_argument('-m', '--merge', action="store",
        help=("Comma-separated list of angle name suffixes. For each suffix "
            "given, a merged data set will be created from all angles with "
            "this suffix in their names, the prefix is ignored."))
    options = parser.parse_args()
    df = parse_dihed_datafile()

    if options.merge:
        suffixes_for_merge = options.merge.split(',')
        log.info("Angle name suffixes for merge:\n%s",
            "\n".join(suffixes_for_merge))
        merge_dataframe_by_suffix(df, suffixes_for_merge)

    if options.two_dimensional:
        suffixes_for_2d = options.two_dimensional.split(',')
        log.info(("Angle name suffixes for 2D histogram for matching angle "
            "name prefixes:\n%s"), "\n".join(suffixes_for_merge))
        assert len(suffixes_for_2d) == 2, "Exactly two suffixes allowed."


def parse_dihed_datafile():
    log.info("Reading '%s' to pandas dataframe.", options.dihedraldatafile)
    # `index_col=False` is required, otherwise indices are NaN.
    df = pd.read_csv(
        options.dihedraldatafile,
        delim_whitespace=True,
        index_col=False)
    log.debug("Columns read:\n%s", "\n".join(c for c in df.columns))
    # We expect the first column name be prefixed with a '#' character.
    # Get rid of that (renaming the columns requires re-assigning the entire
    # `df.columns` list or using the convoluted `rename` method.)
    log.debug("Filtering column names.")
    columnnames = list(df.columns[:])
    for i, name in enumerate(columnnames):
        if name.startswith('#'):
            columnnames[i] = name[1:]
    df.columns = columnnames
    log.debug("Columns after filtering:\n%s" % "\n".join(c for c in df.columns))
    log.debug("Dataframe head:\n%s", df.head())
    return df











if __name__ == "__main__":
    main()


