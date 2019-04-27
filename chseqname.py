#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © %YEAR% %USER% <%MAIL%>
#
# Distributed under terms of the %LICENSE% license.

import sys;
import re;
import argparse;

"""
%HERE%
"""


parser = argparse.ArgumentParser();
parser.add_argument("-f", "--fm", type = int);
parser.add_argument("-t", "--to", type = int);
parser.add_argument("-p", "--pr", type = str);
parser.add_argument("-s", "--sb", type = str);
parser.add_argument("-o", "--qp", type = str, default = "chr");
parser.add_argument("-q", "--qs", type = str, default = ".fa");
args = parser.parse_args();

for i in range(args.fm, args.to):
	print("mv %s%02d%s %s%02d%s" % (args.pr, i, args.sb, args.qp, i-args.fm+1, args.qs));

