/**
 * Copyright Copyright 2013-17 Simon Andrews
 *
 *    This file is part of FastQC.
 *
 *    FastQC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    FastQC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with FastQC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

package uk.ac.babraham.FastQC.Graphs;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;

import javax.swing.JPanel;

import uk.ac.babraham.FastQC.Modules.ModuleConfig;
import uk.ac.babraham.FastQC.Utilities.HotColdColourGradient;

public class TileGraph extends JPanel {

	private String [] xLabels;
	private int [] tiles;
	private double [][]tileBaseMeans;
	private HotColdColourGradient gradient = new HotColdColourGradient();

	public TileGraph (String [] xLabels, int [] tiles, double [][] tileBaseMeans) {
		this.xLabels = xLabels;
		this.tiles = tiles;
		this.tileBaseMeans = tileBaseMeans;

	}
	

	private int getY(double y) {
		return (getHeight()-40) - (int)(((getHeight()-80)/(double)(tiles.length))*y);
	}


	public void paint (Graphics g) {
		super.paint(g);

		g.setColor(Color.WHITE);
		g.fillRect(0, 0, getWidth(), getHeight());
		g.setColor(Color.BLACK);

		int lastY = 0;

		int xOffset = 0;

		for (int i=0;i<tiles.length;i++) {
			String label = ""+tiles[i];
			int width = g.getFontMetrics().stringWidth(label);
			if (width > xOffset) {
				xOffset = width;
			}

			int thisY = getY(i);
			if (i>0 && thisY+g.getFontMetrics().getAscent() > lastY) continue;
			
			g.drawString(label, 2, getY(i));
			lastY = thisY;
		}

		// Give the x axis a bit of breathing space
		xOffset += 5;

		// Draw the graph title
		String graphTitle = "Quality per tile";
		int titleWidth = g.getFontMetrics().stringWidth(graphTitle);
		g.drawString(graphTitle, (xOffset + ((getWidth()-(xOffset+10))/2)) - (titleWidth/2), 30);


		// Now draw the axes
		g.drawLine(xOffset, getHeight()-40, getWidth()-10,getHeight()-40);
		g.drawLine(xOffset, getHeight()-40, xOffset, 40);

		// Draw the xLabel under the xAxis
		String xLabel = "Position in read (bp)";
		g.drawString(xLabel, (getWidth()/2) - (g.getFontMetrics().stringWidth(xLabel)/2), getHeight()-5);


		// Now draw the data points
		int baseWidth = (getWidth()-(xOffset+10))/xLabels.length;
		if (baseWidth<1) baseWidth=1;

		//		System.out.println("Base Width is "+baseWidth);

		// First draw faint boxes over alternating bases so you can see which is which

		// Let's find the longest label, and then work out how often we can draw labels

		int lastXLabelEnd = 0;
		g.setColor(Color.BLACK);

		for (int base=0;base<xLabels.length;base++) {

			String baseNumber = ""+xLabels[base];
			int baseNumberWidth = g.getFontMetrics().stringWidth(baseNumber);
			int baseNumberPosition =  (baseWidth/2)+xOffset+(baseWidth*base)-(baseNumberWidth/2);

			if (baseNumberPosition > lastXLabelEnd) {
				g.drawString(baseNumber,baseNumberPosition, getHeight()-25);
				lastXLabelEnd = baseNumberPosition+baseNumberWidth+5;
			}
		}

		// Now draw the datasets

		if (g instanceof Graphics2D) {
			((Graphics2D)g).setStroke(new BasicStroke(2));
			((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		}

		for (int tile=0;tile<tiles.length;tile++) {
			for (int base=0;base<xLabels.length;base++) {

				g.setColor(getColour(tile,base));

				int x=xOffset+(baseWidth*base);
				int y=getY(tile+1);
				g.fillRect(x, y, baseWidth, getY(tile)-getY(tile+1));

			}

		}


	}
	
	private Color getColour(int tile, int base) {		
		return gradient.getColor(0-tileBaseMeans[tile][base], 0, ModuleConfig.getParam("tile", "error"));
	}


}
