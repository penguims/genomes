#!/usr/local/bin/python
# -*- coding: UTF-8 -*-

from pyecharts import Bar, Grid, Scatter, Style;
from pyecharts.engine import create_default_environment;
import argparse;

parser = argparse.ArgumentParser();
parser.add_argument("-i", "--infile", help = "input file", required = True);
parser.add_argument("-x", "--xaxis", help = "X column", type = int, required = True);
parser.add_argument("-y", "--yaxis", help = "Y Column", type = int, required = True);
parser.add_argument("-c", "--category", help = "category column", type = int, required = True);
parser.add_argument("-t", "--title", help = "chart title", default = "eChart Scatter");
parser.add_argument("-m", "--comment", help = "comment symbot", action = "store_false");
parser.add_argument("-a", "--xmax", help = "max x axis value", type = int, default = 3500);
parser.add_argument("-b", "--ymax", help = "max y axis value", type = int, default = 50000);
parser.add_argument("-n", "--ignore", help = "ignore data the larger than x and y max value", action = "store_false");
parser.add_argument("-o", "--outfile", help = "output file", default = "echarts.html");
parser.add_argument("-f", "--type", help = "default output format", default = "html");
parser.add_argument("-w", "--width", help = "figure width", type = int, default = 1024);
parser.add_argument("-r", "--height", help = "figure height", type = int, default = 768);
parser.add_argument("-z", "--zoom", help = "data zoom", action = "store_false");
args = parser.parse_args();

xdata={};
ydata={};
xaxis_name = "";
yaxis_name = "";
title = "";
with open(args.infile, "r") as fh:
	for ln in fh:
		cells = ln.strip().split("\t");
		key=cells[args.category];
		if ln[0] == "#" and args.comment:
			xaxis_name = cells[args.xaxis];
			yaxis_name = cells[args.yaxis];
			title = cells[args.category]+": "+xaxis_name+" vs. "+yaxis_name;
			continue;
		if not (key in xdata):
			xdata[key]=[];
			ydata[key]=[];
		if args.ignore:
			if cells[args.yaxis] == "-":
				continue;
			if float(cells[args.xaxis]) > args.xmax or float(cells[args.yaxis]) > args.ymax:
				continue;
		xdata[key].append(float(cells[args.xaxis]));
		if cells[args.yaxis] == "-":
			cells[args.yaxis] = "0";
		ydata[key].append(float(cells[args.yaxis]));

grid = Grid(width = args.width, height = args.height);
style = Style(
	title_pos = "center", title_text_size = 24,
);
if args.title == "eChart Scatter":
	pass;
else:
	title = args.title;

scatter = Scatter(
	title, 
	**style.init_style
);


s_style = style.add(
	legend_pos = "center", legend_top = "90%",
	legend_text_size = 18, line_opacity = "80%",
	xaxis_label_textsize = 18, yaxis_label_textsize = 18,
	xaxis_name_size = 18, yaxis_name_size = 18,
	xaxis_name_gap = 30, yaxis_name_gap = 80,
	xaxis_name_pos = "center", yaxis_name_pos = "middle",
	xaxis_line_width = 2, yaxis_line_width = 2,
	xaxis_splitline_show = False, yaxis_splitline_show = False,
	is_splitline_show = False,
	symbol_size = 15,
	is_datazoom_show = args.zoom,
);
for key in xdata:
	scatter.add(
		key, xdata[key], ydata[key],
		xaxis_max = args.xmax, yaxis_max = args.ymax,
		xaxis_name = xaxis_name, yaxis_name = yaxis_name,
		**s_style
	);

grid.add(
	scatter, 
	grid_bottom = "20%"
);
env = create_default_environment(args.type);
env.render_chart_to_file(grid, path = args.outfile);
