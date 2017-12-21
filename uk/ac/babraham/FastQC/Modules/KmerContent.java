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

import java.awt.BorderLayout;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.apache.commons.math3.distribution.BinomialDistribution;


import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Graphs.BaseGroup;
import uk.ac.babraham.FastQC.Graphs.LineGraph;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;

public class KmerContent extends AbstractQCModule {

	private Hashtable<String, Kmer> kmers = new Hashtable<String, Kmer>((int)Math.pow(4, MAX_KMER_SIZE));
	
	private int longestSequence = 0;
	
	/* 2D array, first dimension is the position in the sequence, second is Kmer length */
	private long [][] totalKmerCounts = new long [0][0];
	
	private long skipCount = 0;
	
	private static int MIN_KMER_SIZE = 7;
	private static int MAX_KMER_SIZE = 7;
	
	public boolean calculated = false;
	
	// This is the full set of Kmers to be reported
	private Kmer [] enrichedKmers = null;
	
	// This is the data for the Kmers which are going to be placed on the graph
	private double [][] enrichments = null;
	
	// For the graph we also need to know the scale we need to use on the axes.
	private double minGraphValue = 0;
	private double maxGraphValue = 0;
	
	
	private String [] xCategories = new String[0];
	private String [] xLabels = new String[0];
	
	BaseGroup [] groups;
	public KmerContent () {
		if (FastQCConfig.getInstance().kmer_size != null) {
			int kmerSize = FastQCConfig.getInstance().kmer_size;
			MIN_KMER_SIZE = kmerSize;
			MAX_KMER_SIZE = kmerSize;
		}
	}
	
	public boolean ignoreFilteredSequences() {
		return true;
	}
	
	public boolean ignoreInReport () {
		if (ModuleConfig.getParam("kmer", "ignore") > 0) {
			return true;
		}
		return false;
	}
	
	public JPanel getResultsPanel() {
		
		if (!calculated) calculateEnrichment();
		JPanel returnPanel = new JPanel();
		returnPanel.setLayout(new BorderLayout());
		returnPanel.add(new JLabel("Overrepresented Kmers",JLabel.CENTER),BorderLayout.NORTH);
		
		JSplitPane splitPanel = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
		
		if (enrichedKmers.length > 0) {
			TableModel model = new ResultsTable(enrichedKmers);
			splitPanel.setBottomComponent(new JScrollPane(new JTable(model)));
			splitPanel.setTopComponent(new LineGraph(enrichments, minGraphValue, maxGraphValue, "Position in read (bp)", xLabels, xCategories, "Log2 Obs/Exp"));
			returnPanel.add(splitPanel,BorderLayout.CENTER);
		}
		else {
			returnPanel.add(new JLabel("There are no overrepresented Kmers",JLabel.CENTER),BorderLayout.CENTER);
		}
		
		return returnPanel;
	}
	
	/**
	 * This method simply keeps a count of the number of Kmers of a given size
	 * seen at each position within the run.  We can use this later on to calculate
	 * the enrichment of the Kmers we actually count.
	 * 
	 * We take in the Kmer sequence even though this isn't used in the total counts
	 * we do this because we don't want to count Kmers with Ns in them, but we do
	 * need to ensure that the data structure is expanded to the right size, and if
	 * we have libraries where later positions are Ns in all sequences then our
	 * data structure ends up too short and we crash. 
	 * 
	 * @param position Position within the read.  0 indexed
	 * @param kmerLength Actual length of the Kmer analysed
	 */
	private void addKmerCount (int position,int kmerLength, String kmer) {
	
		
		if (position >= totalKmerCounts.length) {
			// We need to expand the array
			long [][] newCounts = new long[position+1][];
			for (int i=0;i<totalKmerCounts.length;i++) {
				newCounts[i] = totalKmerCounts[i];
			}
			for (int i=totalKmerCounts.length;i<newCounts.length;i++) {
				newCounts[i] = new long[MAX_KMER_SIZE];
			}
			
			totalKmerCounts = newCounts;
		}
		
		if (kmer.indexOf("N") >=0) return;

		++totalKmerCounts[position][kmerLength-1];
		
	}

	private synchronized void calculateEnrichment () {
		
		/*
		 * For each Kmer we work out whether there is a statistically
		 * significant deviation in its coverage at any given position
		 * compared to its average coverage over all positions.
		 */
		
		
		// We'll be grouping together positions later so make up the groups now
		groups = BaseGroup.makeBaseGroups((longestSequence-MIN_KMER_SIZE)+1);

		Vector<Kmer>unevenKmers = new Vector<Kmer>();
				
		Iterator<Kmer> rawKmers = kmers.values().iterator();
		
		while (rawKmers.hasNext()) {
			Kmer k = rawKmers.next();
			char [] chars = k.sequence().toCharArray();

			
			long totalKmerCount = 0;

			
			// This gets us the total number of Kmers of this type in the whole
			// dataset.
			for (int i=0;i<totalKmerCounts.length;i++) {
				totalKmerCount += totalKmerCounts[i][k.sequence().length()-1];
			}

			// This is the expected proportion of all Kmers which have this
			// specific Kmer sequence.  We no longer make any attempt to judge
			// overall enrichment or depletion of this sequence since once you
			// get to longer lengths the distribution isn't flat anyway
			
			float expectedProportion = k.count/(float)totalKmerCount;

			// We now want to go through each of the positions looking for whether
			// this Kmer was seen an unexpected number of times compared to what we
			// expected from the global values
							
				
			float [] obsExpPositions = new float[groups.length];
			float [] binomialPValues = new float[groups.length];
				
			long [] positionCounts = k.getPositions();
				
			for (int g=0;g<groups.length;g++) {
				// This is a summation of the number of Kmers of this length which
				// fall into this base group
				long totalGroupCount = 0;
				
				// This is a summation of the number of hit Kmers which fall within
				// this base group.
				long totalGroupHits = 0;
				for (int p=groups[g].lowerCount()-1;p<groups[g].upperCount() && p < positionCounts.length ;p++) {
					totalGroupCount += totalKmerCounts[p][chars.length-1];
					totalGroupHits += positionCounts[p];
				}
							
				float predicted = expectedProportion * totalGroupCount;
//				obsExpPositions[g] = (float)(Math.log(totalGroupHits/predicted)/Math.log(2));
				obsExpPositions[g] = (float)(totalGroupHits/predicted);
				
				// Now we can run a binomial test to see if there is a significant
				// deviation from what we expect given the number of observations we've
				// made
				
				BinomialDistribution bd = new BinomialDistribution((int)totalGroupCount, expectedProportion);
				if (totalGroupHits > predicted) {
					binomialPValues[g] = (float)((1 - bd.cumulativeProbability((int)totalGroupHits)) * Math.pow(4,chars.length));
				}
				else {
					binomialPValues[g] = 1;
				}
				
			}
			
			k.setObsExpPositions(obsExpPositions);
			
			
			// To keep this we need a p-value below 0.01 and an obs/exp above 5 (actual values are log2 transformed)
			float lowestPValue = 1;
			for (int i=0;i<binomialPValues.length;i++) {
//				if (binomialPValues[i] < 0.01 && obsExpPositions[i] > (Math.log(5)/Math.log(2))) {
				if (binomialPValues[i] < 0.01 && obsExpPositions[i] > 5) {
					if (binomialPValues[i]<lowestPValue) {
						lowestPValue = binomialPValues[i];
					}
				}
			}
			
			if (lowestPValue < 0.01) {
				k.setLowestPValue(lowestPValue);
				unevenKmers.add(k);
			}
			
			
		}
		
		Kmer [] finalKMers = unevenKmers.toArray(new Kmer[0]);
		
		// We sort by the highest degree of enrichment over the average		
		Arrays.sort(finalKMers);
		
		// So we don't end up with stupidly long lists of Kmers in the
		// report we'll only report the top 20
		if (finalKMers.length > 20) {
			Kmer [] shortenedKmers = new Kmer [20];
			for (int i=0;i<shortenedKmers.length;i++) {
				shortenedKmers[i] = finalKMers[i];
			}
			
			finalKMers = shortenedKmers;
		}
		
		// Now we take the enrichment positions for the top 6 hits and
		// record these so we can plot them on a line graph
		enrichments = new double [Math.min(6, finalKMers.length)][];
		xLabels = new String[enrichments.length];
		
		xCategories = new String [groups.length];
		
		for (int i=0;i<xCategories.length;i++) {
			xCategories[i] = groups[i].toString();
		}
		
		for (int k=0;k<enrichments.length;k++) {
			enrichments[k] = new double[groups.length];
			
			float [] obsExpPos = finalKMers[k].getObsExpPositions();
						
			for (int g=0;g<groups.length;g++) {				
				enrichments[k][g] = obsExpPos[g];
				if (obsExpPos[g] > maxGraphValue) maxGraphValue = obsExpPos[g];
				if (obsExpPos[g] < minGraphValue) minGraphValue = obsExpPos[g];
			}
			
			xLabels[k] = finalKMers[k].sequence();
			
		}
		
		minGraphValue = 0;
		
//		System.err.println("Max value="+maxGraphValue+" min value="+minGraphValue);
		
		this.enrichedKmers = finalKMers;		
		
		// Delete the initial data structure so we don't suck up more memory
		// than we have to.
		kmers.clear();
		
		calculated = true;
	}

	
	public void processSequence(Sequence sequence) {
		calculated = false;
		
		/*
		 * The processing done by this module is quite intensive so to speed things
		 * up we don't look at every sequence.  Instead we take only 2% of the 
		 * submitted sequences and extrapolate from these to the full set in the file.
		 */
		++skipCount;
		if (skipCount % 50 != 0) return;
		
		/*
		 * This module uses horrible amounts of memory if allowed to store the full
		 * Kmer content for all positions in really long reads (pacbio reads were the
		 * ones which really broke this).  We'll therefore limit our read lengths to
		 * 500bp since specific Kmer positions beyond that are not likely to be useful
		 */

		String seq;
		
		if (sequence.getSequence().length() > 500) {
			seq = sequence.getSequence().substring(0, 500);		
		}
		else {
			seq = sequence.getSequence();
		}

		if (seq.length() > longestSequence) {
			longestSequence = seq.length();
		}
						
		// Now we go through all of the Kmers to count these
		for (int kmerSize=MIN_KMER_SIZE;kmerSize<=MAX_KMER_SIZE;kmerSize++) {
			for (int i=0;i<=seq.length()-kmerSize;i++) {
				
				String kmer = seq.substring(i, i+kmerSize);
				
				if (kmer.length() != kmerSize) {
					throw new IllegalStateException("String length "+kmer.length()+" wasn't the same as the kmer length "+kmerSize);
				}
				
				// Add to the counts before skipping Kmers containing Ns (see
				// explanation in addKmerCount for the reasoning).
				addKmerCount(i, kmerSize, kmer);
				
				// Skip Kmers containing N
				if (kmer.indexOf("N") >=0) continue;

				if (kmers.containsKey(kmer)) {
					kmers.get(kmer).incrementCount(i);
				}
				else {
					kmers.put(new String(kmer), new Kmer(kmer,i,(seq.length()-kmerSize)+1));
				}

			}
		}
	}
	
	public void reset () {
		calculated = false;
		totalKmerCounts = new long[0][0];
		longestSequence = 0;
		skipCount = 0;
		enrichedKmers = null;
		kmers.clear();
	}

	public String description() {
		return "Identifies short sequences which have uneven representation";
	}

	public String name() {
		return "Kmer Content";
	}

	public boolean raisesError() {
		if (!calculated) calculateEnrichment();
		
		// We raise an error if the most enriched kmer is seen more than 100 times
		// more frequently than we expect.
		
		if (enrichedKmers.length > 0 && 0-Math.log10(enrichedKmers[0].pValue()) > ModuleConfig.getParam("kmer", "error")) return true;
		return false;
	}

	public boolean raisesWarning() {
		if (!calculated) calculateEnrichment();
		
		// We raise a warning if there are any enriched kmers
		if (enrichedKmers.length > 0 && 0-Math.log10(enrichedKmers[0].pValue()) > ModuleConfig.getParam("kmer", "warn")) return true;
		return false;
	}

	public void makeReport(HTMLReportArchive report) throws IOException,XMLStreamException {
		if (!calculated) calculateEnrichment();
		
		if (enrichedKmers.length > 0) {			
			writeSpecificImage(report, new LineGraph(enrichments, minGraphValue, maxGraphValue, "Position in read (bp)", xLabels, xCategories, "Log2 Obs/Exp"),"kmer_profiles.png", "Kmer graph", Math.max(800, groups.length*15), 600);
		}
		
		
		ResultsTable table = new ResultsTable(enrichedKmers);
		
		XMLStreamWriter xhtml = report.xhtmlStream();
		
		if (enrichedKmers.length == 0)
			{
			xhtml.writeStartElement("p");
			xhtml.writeCharacters("No overrepresented Kmers");
			xhtml.writeEndElement();
			}
		
		else
			{
			super.writeTable(report, table);
			}
	}
	
	private class Kmer implements Comparable<Kmer>{
		
		private String sequence;
		private long count = 0;
		private float lowestPValue = 0;
		private float [] obsExpPositions = null;
		private long [] positions = new long[0];
		
		public Kmer (String sequence, int position, int seqLength) {

			// Do this slightly convoluted dance to try to avoid
			// keeping the whole original sequence in memory
			char [] chars = sequence.toCharArray();
			this.sequence = new String(chars);
			count = 1;
			positions = new long[seqLength];
			++positions[position];
		}
		
		public void incrementCount (int position) {
			++count;
			
			if (position >= positions.length) {
				long [] newPositions = new long[position+1];
				for (int i=0;i<positions.length;i++) {
					newPositions[i] = positions[i];
				}
				positions = newPositions;
			}
			
			++positions[position];
			
		}
		
		public long [] getPositions () {
			return positions;
		}
		
		public String sequence () {
			return sequence;
		}
		
		public long count () {
			return count;
		}
		
		public void setLowestPValue (float p) {
			this.lowestPValue = p;
		}
		
		public void setObsExpPositions (float [] oePositions) {
			this.obsExpPositions = oePositions;
		}
		
		public float [] getObsExpPositions () {
			return obsExpPositions;
		}
		
		public float pValue () {
			return lowestPValue;
		}
				
		public float maxObsExp () {
			float max = 0;
			for (int i=0;i<obsExpPositions.length;i++) {
				if (obsExpPositions[i]>max) max = obsExpPositions[i];
			}
			return max;
		}

		public int maxPosition () {
			float max = 0;
			int position = 0;
			for (int i=0;i<obsExpPositions.length;i++) {
				if (obsExpPositions[i]>max) {
					max = obsExpPositions[i];
					position = i+1;
				}
			}
			
			if (position == 0) {
				System.err.println("No value > 0 for "+sequence);
				position = 1;
			}
			
			return position;
		}

		public int compareTo(Kmer o) {
			return Float.compare(o.maxObsExp(), maxObsExp());
		}
	}
	
	
	private class ResultsTable extends AbstractTableModel
		{
		private static final long serialVersionUID = 1L;
		private Kmer [] kmers;
		
		public ResultsTable (Kmer [] kmers) {
			this.kmers = kmers;
		}
		
		
		// Sequence - Count - Obs/Exp
		public int getColumnCount() {
			return 5;
		}

		public int getRowCount() {
			return kmers.length;
		}

		public Object getValueAt(int rowIndex, int columnIndex) {
			switch (columnIndex) {
				case 0: return kmers[rowIndex].sequence();
				case 1: return kmers[rowIndex].count()*5;
				case 2: return kmers[rowIndex].pValue();
				case 3: return kmers[rowIndex].maxObsExp();
				case 4: return groups[kmers[rowIndex].maxPosition()-1].toString();
						
			}
			return null;
		}
		
		public String getColumnName (int columnIndex) {
			switch (columnIndex) {
				case 0: return "Sequence";
				case 1: return "Count";
				case 2: return "PValue";
				case 3: return "Obs/Exp Max";
				case 4: return "Max Obs/Exp Position";
			}
			return null;
		}
		
		public Class<?> getColumnClass (int columnIndex) {
			switch (columnIndex) {
			case 0: return String.class;
			case 1: return Integer.class;
			case 2: return Float.class;
			case 3: return Float.class;
			case 4: return String.class;
		}
		return null;
			
		}
	}

}
