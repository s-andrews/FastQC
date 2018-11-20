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

import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Graphs.LineGraph;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;

public class SequenceLengthDistribution extends AbstractQCModule {

	private long [] lengthCounts = new long[0];
	private double [] graphCounts = null;
	private String [] xCategories = new String[0];
	private double max = 0;
	private boolean calculated = false;
	
	public JPanel getResultsPanel() {
	
		if (!calculated) calculateDistribution();
				
		return new LineGraph(new double [][] {graphCounts}, 0d, max, "Sequence Length (bp)",new String [] {"Sequence Length"}, xCategories, "Distribution of sequence lengths over all sequences");
	}
	
	public boolean ignoreFilteredSequences() {
		return true;
	}
	
	public boolean ignoreInReport () {
		if (ModuleConfig.getParam("sequence_length", "ignore") > 0) {
			return true;
		}
		return false;
	}
	
	private synchronized void calculateDistribution () {
		int maxLen = 0;
		int minLen = -1;
		max = 0;
		
		// Find the min and max lengths
		for (int i=0;i<lengthCounts.length;i++) {
			if (lengthCounts[i]>0) {
				if (minLen < 0) {
					minLen = i;
				}
				maxLen = i;
			}
		}
		
		// We can get a -1 value for min if there aren't any valid sequences
		// at all.
		
		if (minLen < 0) minLen = 0;
		
		// We put one extra category either side of the actual size
		if (minLen>0) minLen--;
		maxLen++;
		
		int [] startAndInterval = getSizeDistribution(minLen, maxLen);
				
		// Work out how many categories we need
		int categories = 0;
		int currentValue = startAndInterval[0];
		while (currentValue<= maxLen) {
			++categories;
			currentValue+= startAndInterval[1];
		}
		
		graphCounts = new double[categories];
		xCategories = new String[categories];
		
		for (int i=0;i<graphCounts.length;i++) {
			
			int minValue = startAndInterval[0]+(startAndInterval[1]*i);
			int maxValue = (startAndInterval[0]+(startAndInterval[1]*(i+1)))-1;

			if (maxValue > maxLen) {
				maxValue = maxLen;
			}
			
			for (int bp=minValue;bp<=maxValue;bp++) {
				if (bp < lengthCounts.length) {
					graphCounts[i] += lengthCounts[bp];
				}
			}

			if (startAndInterval[1] == 1) {
				xCategories[i] = ""+minValue;
			}
			else {
				xCategories[i] = minValue+"-"+maxValue;
			}
			
			if (graphCounts[i] > max) max = graphCounts[i];
		}
		calculated = true;
	}

	public void processSequence(Sequence sequence) {
		int seqLen = sequence.getSequence().length();

		if (seqLen+2 > lengthCounts.length) {
			long [] newLengthCounts = new long[seqLen+2];
			for (int i=0;i<lengthCounts.length;i++) {
				newLengthCounts[i] = lengthCounts[i];
			}
			lengthCounts = newLengthCounts;
		}
		
		++lengthCounts[seqLen];
		
	}
	
	private int [] getSizeDistribution (int min, int max) {
		
		// We won't group if they've asked us not to
		if (FastQCConfig.getInstance().nogroup) {
			return(new int [] {min,1});
		}
		
		int base = 1;
		
		while (base > (max-min)) {
			base /= 10;
		}
				
		int interval;
		int starting;

		int [] divisions = new int [] {1,2,5};
		
		OUTER: while (true) {
			
			for (int d=0;d<divisions.length;d++) {
				int tester = base * divisions[d];
				if (((max-min) / tester) <= 50) {
					interval = tester;
					break OUTER;
				}
			}
		
			base *=10;
			
		}
		
		// Now we work out the first value to be plotted
		int basicDivision = (int)(min/interval);
				
		int testStart = basicDivision * interval;
				
		starting = testStart;
		
		return new int[] {starting,interval};
		
	}
	
	
	
	public void reset () {
		lengthCounts = new long[0];
	}

	public String description() {
		return "Shows the distribution of sequence length over all sequences";
	}

	public String name() {
		return "Sequence Length Distribution";
	}

	public boolean raisesError() {
		if (!calculated) calculateDistribution();
		
		// See if they've turned this test off
		if (ModuleConfig.getParam("sequence_length", "error") == 0) {
			return false;
		}

		// We might not have any sequences so only check if we do
		if (lengthCounts.length > 0) {
			// Empty sequences get us an error
			if (lengthCounts[0] > 0) {
				return true;
			}
		}
		return false;
	}

	public boolean raisesWarning() {
		if (!calculated) calculateDistribution();

		// See if they've turned this test off
		if (ModuleConfig.getParam("sequence_length", "warn") == 0) {
			return false;
		}
		
		// Warn if they're not all the same length
		boolean seenLength = false;
		for (int i=0;i<lengthCounts.length;i++) {
			if (lengthCounts[i] > 0) {
				if (seenLength) {
					return true;
				}
				else {
					seenLength = true;
				}
			}
		}
		return false;
	}

	public void makeReport(HTMLReportArchive report) throws IOException,XMLStreamException {
		if (!calculated) calculateDistribution();

		writeDefaultImage(report, "sequence_length_distribution.png", "Sequence length distribution",  Math.max(800, graphCounts.length*15), 600);
		
		StringBuffer sb = report.dataDocument();
		sb.append("#Length\tCount\n");
		for (int i=0;i<xCategories.length;i++) {
			// Remove any padding we added to make the graph look better
			if ((i==0 || i== xCategories.length-1) && graphCounts[i] == 0) {
				continue;
			}
			sb.append(xCategories[i]);
			sb.append("\t");
			sb.append(graphCounts[i]);
			sb.append("\n");
		}
	}

}
