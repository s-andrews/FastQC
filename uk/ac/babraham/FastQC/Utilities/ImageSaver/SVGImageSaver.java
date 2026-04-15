/**
 * Copyright 2009- 21 Simon Andrews
 *
 *    This file is part of SeqMonk.
 *
 *    SeqMonk is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    SeqMonk is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with SeqMonk; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
package uk.ac.babraham.FastQC.Utilities.ImageSaver;

import java.awt.Component;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * A utility class which acts as a wrapper for the SVG or PNG generating
 * code which can be used to save (almost) any component which uses the
 * standard Graphics interface to draw itself.
 */
public class SVGImageSaver {

	private static final Pattern LINE_PATTERN = Pattern.compile(
		"<line x1=\"(\\d+)\" y1=\"(\\d+)\" x2=\"(\\d+)\" y2=\"(\\d+)\" stroke=\"(rgb\\([^)]+\\))\" stroke-width=\"(\\d+)\"/>"
	);

	private static final Pattern RECT_PATTERN = Pattern.compile(
		"<rect width=\"(\\d+)\" height=\"(\\d+)\" x=\"(\\d+)\" y=\"(\\d+)\" fill=\"(rgb\\([^)]+\\))\"/>"
	);

	public static String saveImage (Component c, OutputStream os) {
		PrintWriter pr = new PrintWriter(os);
		String svgData = SVGGenerator.writeSVG(c);
		svgData = optimizeSvg(svgData);
		pr.write(svgData);
		pr.flush();
		return(svgData);
	}

	private static String optimizeSvg(String svg) {
		svg = mergeLinesToPolylines(svg);
		svg = mergeRects(svg);
		return svg;
	}

	/** Merge consecutive same-colour line elements into polyline elements. */
	private static String mergeLinesToPolylines(String svg) {
		String[] lines = svg.split("\n");
		StringBuilder result = new StringBuilder(svg.length());
		ArrayList<String[]> group = new ArrayList<>();
		String groupStroke = null;
		String groupWidth = null;

		for (String line : lines) {
			Matcher m = LINE_PATTERN.matcher(line.trim());
			if (m.matches()) {
				String stroke = m.group(5);
				String width = m.group(6);
				if ("rgb(0,0,0)".equals(stroke) || "rgb(180,180,180)".equals(stroke)) {
					flushLineGroup(result, group, groupStroke, groupWidth);
					group.clear();
					groupStroke = null;
					result.append(line).append("\n");
				} else if (stroke.equals(groupStroke) && width.equals(groupWidth)) {
					group.add(new String[]{m.group(1), m.group(2), m.group(3), m.group(4)});
				} else {
					flushLineGroup(result, group, groupStroke, groupWidth);
					group.clear();
					groupStroke = stroke;
					groupWidth = width;
					group.add(new String[]{m.group(1), m.group(2), m.group(3), m.group(4)});
				}
			} else {
				flushLineGroup(result, group, groupStroke, groupWidth);
				group.clear();
				groupStroke = null;
				result.append(line).append("\n");
			}
		}
		flushLineGroup(result, group, groupStroke, groupWidth);
		return result.toString();
	}

	private static void flushLineGroup(StringBuilder sb, ArrayList<String[]> group, String stroke, String width) {
		if (group.size() <= 2) {
			for (String[] seg : group) {
				sb.append("<line x1=\"").append(seg[0]).append("\" y1=\"").append(seg[1])
				  .append("\" x2=\"").append(seg[2]).append("\" y2=\"").append(seg[3])
				  .append("\" stroke=\"").append(stroke).append("\" stroke-width=\"").append(width)
				  .append("\"/>\n");
			}
			return;
		}
		StringBuilder points = new StringBuilder();
		points.append(group.get(0)[0]).append(",").append(group.get(0)[1]);
		for (String[] seg : group) {
			points.append(" ").append(seg[2]).append(",").append(seg[3]);
		}
		sb.append("<polyline points=\"").append(points)
		  .append("\" stroke=\"").append(stroke).append("\" stroke-width=\"").append(width)
		  .append("\" fill=\"none\"/>\n");
	}

	/** Merge consecutive same-colour filled rects on the same row into wider rects. */
	private static String mergeRects(String svg) {
		String[] lines = svg.split("\n");
		StringBuilder result = new StringBuilder(svg.length());
		ArrayList<int[]> group = new ArrayList<>();
		String groupFill = null;

		for (String line : lines) {
			Matcher m = RECT_PATTERN.matcher(line.trim());
			if (m.matches()) {
				int w = Integer.parseInt(m.group(1));
				int h = Integer.parseInt(m.group(2));
				int x = Integer.parseInt(m.group(3));
				int y = Integer.parseInt(m.group(4));
				String fill = m.group(5);
				if (w > 100 || h > 100) {
					flushRectGroup(result, group, groupFill);
					group.clear();
					groupFill = null;
					result.append(line).append("\n");
				} else if (fill.equals(groupFill) && !group.isEmpty()
						&& y == group.get(0)[3] && h == group.get(0)[1]
						&& x == group.get(0)[2] + group.get(0)[0] * group.size()) {
					group.add(new int[]{w, h, x, y});
				} else {
					flushRectGroup(result, group, groupFill);
					group.clear();
					groupFill = fill;
					group.add(new int[]{w, h, x, y});
				}
			} else {
				flushRectGroup(result, group, groupFill);
				group.clear();
				groupFill = null;
				result.append(line).append("\n");
			}
		}
		flushRectGroup(result, group, groupFill);
		return result.toString();
	}

	private static void flushRectGroup(StringBuilder sb, ArrayList<int[]> group, String fill) {
		if (group.isEmpty()) return;
		int mergedWidth = group.get(0)[0] * group.size();
		int[] first = group.get(0);
		sb.append("<rect width=\"").append(mergedWidth).append("\" height=\"").append(first[1])
		  .append("\" x=\"").append(first[2]).append("\" y=\"").append(first[3])
		  .append("\" fill=\"").append(fill).append("\"/>\n");
	}

}
