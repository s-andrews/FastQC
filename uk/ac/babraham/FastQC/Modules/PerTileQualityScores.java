/**
 * Copyright Copyright 2010-17 Simon Andrews
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
package uk.ac.babraham.FastQC.Modules;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

import javax.swing.JPanel;
import javax.xml.stream.XMLStreamException;

import uk.ac.babraham.FastQC.Graphs.BaseGroup;
import uk.ac.babraham.FastQC.Graphs.TileGraph;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;
import uk.ac.babraham.FastQC.Sequence.QualityEncoding.PhredEncoding;
import uk.ac.babraham.FastQC.Utilities.QualityCount;

public class PerTileQualityScores extends AbstractQCModule {


	public HashMap<Integer, QualityCount []> perTileQualityCounts = new HashMap<Integer, QualityCount[]>();
	private int currentLength = 0;
	private double [][] means = null;
	private String [] xLabels;
	private int [] tiles;
	private int high = 0;
	PhredEncoding encodingScheme;
	private boolean calculated = false;
	
	private long totalCount = 0;
	
	private int splitPosition = -1;
	
	private double maxDeviation = 0;

	private boolean ignoreInReport = false;

	public JPanel getResultsPanel() {

		if (!calculated) getPercentages();

		return new TileGraph(xLabels, tiles, means);

	}

	public boolean ignoreFilteredSequences() {
		return true;
	}

	public boolean ignoreInReport () {
		if (ignoreInReport || ModuleConfig.getParam("tile", "ignore") > 0  || currentLength == 0) {
			return true;
		}
		return false;
	}

	private synchronized void getPercentages () {

		char [] range = calculateOffsets();
		encodingScheme = PhredEncoding.getFastQEncodingOffset(range[0]);
		high = range[1] - encodingScheme.offset();
		if (high < 35) {
			high = 35;
		}

		BaseGroup [] groups = BaseGroup.makeBaseGroups(currentLength);

		Integer [] tileNumbers = perTileQualityCounts.keySet().toArray(new Integer[0]);
		
		Arrays.sort(tileNumbers);
		
		tiles = new int[tileNumbers.length];
		for (int i=0;i<tiles.length;i++) {
			tiles[i] = tileNumbers[i];
		}
		
		means = new double[tileNumbers.length][groups.length];
		xLabels = new String[groups.length];

		for (int t=0;t<tileNumbers.length;t++){
			for (int i=0;i<groups.length;i++) {
				if (t==0)
					xLabels[i] = groups[i].toString();

				int minBase = groups[i].lowerCount();
				int maxBase = groups[i].upperCount();
				means[t][i] = getMean(tileNumbers[t],minBase,maxBase,encodingScheme.offset());
			}
		}
		
		// Now we normalise across each column to see if there are any tiles with unusually
		// high or low quality.
		
		double maxDeviation = 0;
		
		double [] averageQualitiesPerGroup = new double[groups.length];
		
		for (int t=0;t<tileNumbers.length;t++) {		
			for (int i=0;i<groups.length;i++) {
				averageQualitiesPerGroup[i] += means[t][i];
			}
		}
		
		for (int i=0;i<averageQualitiesPerGroup.length;i++) {
			averageQualitiesPerGroup[i] /= tileNumbers.length;
		}

		for (int i=0;i<groups.length;i++) {
			for (int t=0;t<tileNumbers.length;t++) {
				means[t][i] -= averageQualitiesPerGroup[i];
				if (Math.abs(means[t][i])> maxDeviation) {
					maxDeviation = Math.abs(means[t][i]);
				}
			}
		}
		
		this.maxDeviation = maxDeviation;

		calculated = true;

	}

	private char [] calculateOffsets () {
		// Works out from the set of chars what is the most
		// likely encoding scale for this file.

		char minChar = 0;
		char maxChar = 0;

		// Iterate through the tiles to check them all in case
		// we're dealing with unrepresentative data in the first one.
		Iterator<QualityCount[]> qit = perTileQualityCounts.values().iterator();
		
		while (qit.hasNext()) { 
			
			QualityCount [] qualityCounts = qit.next();
	
			for (int q=0;q<qualityCounts.length;q++) {
				if (minChar == 0) {
					minChar = qualityCounts[q].getMinChar();
					maxChar = qualityCounts[q].getMaxChar();
				}
				else {
					if (qualityCounts[q].getMinChar() < minChar) {
						minChar = qualityCounts[q].getMinChar();
					}
					if (qualityCounts[q].getMaxChar() > maxChar) {
						maxChar = qualityCounts[q].getMaxChar();
					}
				}
			}
		}

		return new char[] {minChar,maxChar};
	}

	public void processSequence(Sequence sequence) {

		// Check if we can skip counting because the module is being ignored anyway
		if (totalCount == 0) {
			if (ModuleConfig.getParam("tile", "ignore") > 0) {
				ignoreInReport = true;
			}
		}

		
		// Don't waste time calculating this if we're not going to use it anyway
		if (ignoreInReport) return;
		
		// Don't bother with sequences with zero length as they don't have any 
		// quality information anyway.
		if (sequence.getQualityString().length() == 0) {
			return;
		}
				
		calculated = false;

		// Try to find the tile id.  This can come in one of two forms:
		//		@HWI-1KL136:211:D1LGAACXX:1:1101:18518:48851 3:N:0:ATGTCA
		//                                ^
		//		@HWUSI-EAS493_0001:2:1:1000:16900#0/1
		//                         ^

		// These would appear at sections 2 or 3 of an array split on :

		// This module does quite a lot of work and ends up being the limiting
		// step when calculating.  We'll therefore take only a sample of the 
		// sequences to try to get a representative selection.  We'll use the 
		// first 10k sequences in case we're dealing with a very small file
		// and then take 10% of the rest.
		
		++totalCount;
		if (totalCount > 10000 && totalCount % 10 != 0) return;
		
		// First try to split the id by :
		int tile = 0;
		
		String [] splitID = sequence.getID().split(":");


		// If there are 7 or more fields then it's a 1.8+ file
		try {

						
			if (splitPosition >=0) {
				// We've found a split position before so let's try to use it again
				
				
				if (splitID.length <= splitPosition) {
					// There isn't enough data in this header to split the way we did before
					throw new NumberFormatException("Can't extract a number - not enough data");
				}

				tile = Integer.parseInt(splitID[splitPosition]);
			}
			
			else if (splitID.length>=7) {
				splitPosition = 4;
				tile = Integer.parseInt(splitID[4]);
			}
			else if (splitID.length >=5) {
				splitPosition = 2;
				// We can try the older format
				tile = Integer.parseInt(splitID[2]);
			}
			else {
				// We're not going to get a tile out of this
				ignoreInReport = true;
				return;
			}
			
			
		}
		catch (NumberFormatException nfe) {
			// This doesn't conform
			ignoreInReport = true;
			return;
		}

		char [] qual = sequence.getQualityString().toCharArray();
		if (currentLength < qual.length) {

			Iterator<Integer> tiles = perTileQualityCounts.keySet().iterator();
			while (tiles.hasNext()) {
				int thisTile = tiles.next();

				QualityCount [] qualityCounts = perTileQualityCounts.get(thisTile);
				QualityCount [] qualityCountsNew = new QualityCount[qual.length];

				for (int i=0;i<qualityCounts.length;i++) {
					qualityCountsNew[i] = qualityCounts[i];
				}
				for (int i=qualityCounts.length;i<qualityCountsNew.length;i++) {
					qualityCountsNew[i] = new QualityCount();				
				}
				perTileQualityCounts.put(thisTile, qualityCountsNew);
			}

			currentLength = qual.length;

		}

		if (! perTileQualityCounts.containsKey(tile)) {
			
			if (perTileQualityCounts.size() > 1000) {
				// There are too many tiles, so we're probably parsing this wrong.
				// Let's give up
				System.err.println("Too many tiles (>1000) so giving up trying to do per-tile qualities since we're probably parsing the file wrongly");
				ignoreInReport = true;
				perTileQualityCounts.clear();
				return;
			}
			
			QualityCount [] qualityCounts = new QualityCount[currentLength];
			for (int i=0;i<currentLength;i++) {
				qualityCounts[i] = new QualityCount();
			}

			perTileQualityCounts.put(tile, qualityCounts);
		}

		QualityCount [] qualityCounts = perTileQualityCounts.get(tile);

		for (int i=0;i<qual.length;i++) {
			qualityCounts[i].addValue(qual[i]);
		}

	}

	public void reset () {
		totalCount = 0;
		perTileQualityCounts = new HashMap<Integer, QualityCount[]>();
	}

	public String description() {
		return "Shows the perl tile Quality scores of all bases at a given position in a sequencing run";
	}

	public String name() {
		return "Per tile sequence quality";
	}

	public boolean raisesError() {
		if (!calculated) getPercentages();

		if (maxDeviation > ModuleConfig.getParam("tile", "error")) return true;
		return false;
	}

	public boolean raisesWarning() {
		if (!calculated) getPercentages();
		
		if (maxDeviation > ModuleConfig.getParam("tile", "warn")) return true;
		return false;
	}

	public void makeReport(HTMLReportArchive report) throws IOException,XMLStreamException {
		if (!calculated) getPercentages();
		
		writeDefaultImage(report, "per_tile_quality.png", "Per tile quality graph", Math.max(800, xLabels.length*15), 600);

		StringBuffer sb = report.dataDocument();
		sb.append("#Tile\tBase\tMean\n");

		for (int t=0;t<tiles.length;t++) {
			for (int i=0;i<means[t].length;i++) {

				sb.append(tiles[t]);
				sb.append("\t");

				sb.append(xLabels[i]);
				sb.append("\t");

				sb.append(means[t][i]);

				sb.append("\n");
			}
		}
	}

	private double getMean (int tile, int minbp, int maxbp, int offset) {
		int count = 0;
		double total = 0;

		QualityCount [] qualityCounts = perTileQualityCounts.get(tile);

		for (int i=minbp-1;i<maxbp;i++) {
			if (qualityCounts[i].getTotalCount() > 0) {
				count++;
				total += qualityCounts[i].getMean(offset);
			}
		}

		if (count > 0) {
			return total/count;
		}
		return 0;

	}


}
