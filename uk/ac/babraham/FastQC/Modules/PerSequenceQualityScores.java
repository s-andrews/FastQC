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

import javax.swing.JPanel;
import javax.xml.stream.XMLStreamException;

import uk.ac.babraham.FastQC.Graphs.LineGraph;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;
import uk.ac.babraham.FastQC.Sequence.QualityEncoding.PhredEncoding;

public class PerSequenceQualityScores extends AbstractQCModule {

	private HashMap<Integer, Long> averageScoreCounts = new HashMap<Integer, Long>();
	private double [] qualityDistribution = null;
	private int [] xCategories = new int[0];
	private char lowestChar = 126;
	private int maxCount = 0;
	private int mostFrequentScore;
	private boolean calculated = false;
	
	public JPanel getResultsPanel() {
	
		if (!calculated) calculateDistribution();
				
		return new LineGraph(new double [][] {qualityDistribution}, 0d, maxCount, "Mean Sequence Quality (Phred Score)",new String [] {"Average Quality per read"}, xCategories, "Quality score distribution over all sequences");
	}
	
	public boolean ignoreInReport () {
		// We don't show this if they didn't have any quality data.
		if (ModuleConfig.getParam("quality_sequence", "ignore") > 0  || averageScoreCounts.size() == 0) {
			return true;
		}
		return false;
	}
	
	private synchronized void calculateDistribution () {
		
		PhredEncoding encoding = PhredEncoding.getFastQEncodingOffset(lowestChar);
		
		Integer [] rawScores = averageScoreCounts.keySet().toArray(new Integer [0]);
		Arrays.sort(rawScores);
		
		// We'll run from the lowest to the highest
		qualityDistribution = new double [1+(rawScores[rawScores.length-1]-rawScores[0])] ;
		
		xCategories = new int[qualityDistribution.length];
		
		for (int i=0;i<qualityDistribution.length;i++) {
			xCategories[i] = (rawScores[0]+i)-encoding.offset();
			if (averageScoreCounts.containsKey(rawScores[0]+i)) {
				qualityDistribution[i] = averageScoreCounts.get(rawScores[0]+i);
			}
		}

		
		for (int i=0;i<qualityDistribution.length;i++) {
			if (qualityDistribution[i]>maxCount) {
				maxCount = (int)qualityDistribution[i];
				mostFrequentScore = xCategories[i];
			}
		}
				
		calculated = true;
	}

	public void processSequence(Sequence sequence) {
				
		char [] seq = sequence.getQualityString().toCharArray();
		int averageQuality = 0;
		
		for (int i=0;i<seq.length;i++) {
			if (seq[i] < lowestChar) {
				lowestChar = seq[i];
			}
			averageQuality += seq[i];
		}

		if (seq.length > 0) {
			averageQuality /= seq.length;
					
			if (averageScoreCounts.containsKey(averageQuality)) {
				long currentCount = averageScoreCounts.get(averageQuality);
				currentCount++;
				averageScoreCounts.put(averageQuality, currentCount);
			}
			else {
				averageScoreCounts.put(averageQuality, 1L);
			}
		}
	}
	
	public void reset () {
		averageScoreCounts.clear();
		lowestChar = 126;
		maxCount = 0;
		calculated = false;
	}

	public String description() {
		return "Shows the distribution of average quality scores for whole sequences";
	}

	public String name() {
		return "Per sequence quality scores";
	}

	public boolean raisesError() {
		if (!calculated) calculateDistribution();

		if (mostFrequentScore <= ModuleConfig.getParam("quality_sequence", "error")) return true;
		
		return false;
	}

	public boolean raisesWarning() {
		if (!calculated) calculateDistribution();

		if (mostFrequentScore <= ModuleConfig.getParam("quality_sequence", "warn")) return true;
		
		return false;
	}

	public void makeReport(HTMLReportArchive report) throws IOException,XMLStreamException {
		if (!calculated) calculateDistribution();
		
		writeDefaultImage(report, "per_sequence_quality.png", "Per Sequence quality graph", 800, 600);

		StringBuffer sb = report.dataDocument();
		sb.append("#Quality\tCount\n");
		for (int i=0;i<xCategories.length;i++) {
			sb.append(xCategories[i]);
			sb.append("\t");
			sb.append(qualityDistribution[i]);
			sb.append("\n");
		}
	}


	public boolean ignoreFilteredSequences() {
		return true;
	}

}
