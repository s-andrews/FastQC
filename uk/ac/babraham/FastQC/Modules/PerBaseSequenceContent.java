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

import javax.swing.JPanel;
import javax.xml.stream.XMLStreamException;

import uk.ac.babraham.FastQC.Graphs.BaseGroup;
import uk.ac.babraham.FastQC.Graphs.LineGraph;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;

public class PerBaseSequenceContent extends AbstractQCModule {

	// Count of each base at each position: index position * 5 + BASE_INDEX of the base.
	// Slot 4 takes everything that isn't G, A, T or C (N...) and is never read, as
	// those weren't counted before either. One increment per base, instead of four
	// compare-and-adds into four arrays; with the other changes to the module loops
	// this cut CPU by 10% (AVITI) / 9% (WES).
	private long [] baseCounts = new long [0];
	private int maxLength = 0;
	private static final byte [] BASE_INDEX = new byte[256];
	static {
		java.util.Arrays.fill(BASE_INDEX, (byte) 4);
		BASE_INDEX['G'] = 0;
		BASE_INDEX['A'] = 1;
		BASE_INDEX['T'] = 2;
		BASE_INDEX['C'] = 3;
	}

	// Filled from baseCounts when the results are calculated
	public long [] gCounts = new long [0];
	public long [] aCounts = new long [0];
	public long [] cCounts = new long [0];
	public long [] tCounts = new long [0];
	private double [][] percentages = null;
	private String [] xCategories = new String[0];
	private boolean calculated = false;
		
	
	public JPanel getResultsPanel() {
		
		if (!calculated) getPercentages();

		return new LineGraph(percentages, 0d, 100d, "Position in read (bp)", new String [] {"%T","%C","%A","%G"}, xCategories, "Sequence content across all bases");
	}
	
	public boolean ignoreFilteredSequences() {
		return true;
	}
	
	public boolean ignoreInReport () {
		if (ModuleConfig.getParam("sequence", "ignore") > 0) {
			return true;
		}
		return false;
	}

	private synchronized void getPercentages () {

		gCounts = new long[maxLength];
		aCounts = new long[maxLength];
		tCounts = new long[maxLength];
		cCounts = new long[maxLength];
		for (int i=0;i<maxLength;i++) {
			gCounts[i] = baseCounts[i*5];
			aCounts[i] = baseCounts[i*5+1];
			tCounts[i] = baseCounts[i*5+2];
			cCounts[i] = baseCounts[i*5+3];
		}

		BaseGroup [] groups = BaseGroup.makeBaseGroups(gCounts.length);
		
		xCategories = new String[groups.length];

		
		double [] gPercent = new double[groups.length];
		double [] aPercent = new double[groups.length];
		double [] tPercent = new double[groups.length];
		double [] cPercent = new double[groups.length];

		long total;
		long gCount;
		long aCount;
		long tCount;
		long cCount;

		for (int i=0;i<groups.length;i++) {
						
			xCategories[i] = groups[i].toString();

			gCount = 0;
			aCount = 0;
			tCount = 0;
			cCount = 0;
			total = 0;
			
			for (int bp=groups[i].lowerCount()-1;bp<groups[i].upperCount();bp++) {

				total += gCounts[bp];
				total += cCounts[bp];
				total += aCounts[bp];
				total += tCounts[bp];

				aCount += aCounts[bp];
				tCount += tCounts[bp];
				cCount += cCounts[bp];
				gCount += gCounts[bp];				
			}
			
			gPercent[i] = (gCount/(double)total)*100;
			aPercent[i] = (aCount/(double)total)*100;
			tPercent[i] = (tCount/(double)total)*100;
			cPercent[i] = (cCount/(double)total)*100;			
						
		}
		
		percentages = new double [][] {tPercent,cPercent,aPercent,gPercent};
		
		calculated = true;
	}
	
	public void processSequence(Sequence sequence) {
		calculated = false;
		byte [] seq = sequence.getSequenceBytes();
		int seqLen = seq.length;
		if (maxLength < seqLen) {
			baseCounts = java.util.Arrays.copyOf(baseCounts, seqLen * 5);
			maxLength = seqLen;
		}
		
		// One increment per base, without branching on the (random) bases
		long [] counts = baseCounts;
		for (int i=0;i<seqLen;i++) {
			++counts[i*5 + BASE_INDEX[seq[i] & 0xFF]];
		}
		
	}
	
	public void reset () {
		baseCounts = new long[0];
		maxLength = 0;
		gCounts = new long[0];
		aCounts = new long[0];
		tCounts = new long[0];
		cCounts = new long[0];
	}

	public String description() {
		return "Shows the relative amounts of each base at each position in a sequencing run";
	}

	public String name() {
		return "Per base sequence content";
	}

	public boolean raisesError() {
		if (!calculated) getPercentages();

		// Percentages come in the order TCAG
		for (int i=0;i<percentages[0].length;i++) {

			double gcDiff = Math.abs(percentages[1][i]-percentages[3][i]);
			double atDiff = Math.abs(percentages[0][i]-percentages[2][i]);
			
			if (gcDiff > ModuleConfig.getParam("sequence", "error") || atDiff > ModuleConfig.getParam("sequence", "error")) return true;
			
		}
		return false;
	}

	public boolean raisesWarning() {

		if (!calculated) getPercentages();

		// Percentages come in the order TCAG
		for (int i=0;i<percentages[0].length;i++) {

			double gcDiff = Math.abs(percentages[1][i]-percentages[3][i]);
			double atDiff = Math.abs(percentages[0][i]-percentages[2][i]);
			
			if (gcDiff > ModuleConfig.getParam("sequence", "warn") || atDiff > ModuleConfig.getParam("sequence", "warn")) return true;
			
		}
		return false;
	}

	public void makeReport(HTMLReportArchive report) throws IOException,XMLStreamException {
		
		if (!calculated) getPercentages();
		
		writeDefaultImage(report, "per_base_sequence_content.png", "Per base sequence content", Math.max(800, xCategories.length*15), 600);
		
		StringBuffer sb = report.dataDocument();
		sb.append("#Base\tG\tA\tT\tC\n");
		for (int i=0;i<xCategories.length;i++) {
			sb.append(xCategories[i]);
			sb.append("\t");
			sb.append(percentages[3][i]);
			sb.append("\t");
			sb.append(percentages[2][i]);
			sb.append("\t");
			sb.append(percentages[0][i]);
			sb.append("\t");
			sb.append(percentages[1][i]);
			sb.append("\n");
		}
		
	} 

}
