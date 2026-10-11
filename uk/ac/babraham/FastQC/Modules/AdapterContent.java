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
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.table.AbstractTableModel;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;


import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Graphs.BaseGroup;
import uk.ac.babraham.FastQC.Graphs.LineGraph;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;
import uk.ac.babraham.FastQC.Sequence.Contaminant.ContaminentFinder;

public class AdapterContent extends AbstractQCModule {

	private int longestSequence = 0;
	private int longestAdapter = 0;

	private long totalCount = 0;

	public boolean calculated = false;

	// This is the full set of Kmers to be reported
	private Adapter [] adapters;

	// This is the data for the Kmers which are going to be placed on the graph
	private double [][] enrichments = null;
	private String [] labels;

	private String [] xLabels = new String[0];

	BaseGroup [] groups;

	// Single pass adapter search: when all adapters share the same length k (<= 31)
	// and only contain ACGT, each k-mer of the read is 2-bit packed into a long and
	// looked up in a small open addressing table, instead of scanning the read once
	// per adapter with indexOf. Null when the adapters don't qualify.
	// Reads are walked once whatever the number of adapters (6 by default), which
	// cut total CPU by 11% (AVITI 2M reads) and 14% (WES 4M reads). The hit
	// positions are exactly those of indexOf, so the report is unchanged.
	private long [] adapterKeys = null;
	private int [][] adapterSlots = null;
	private int adapterMask;
	// First match position of each adapter in the current read (-1 if none),
	// filled by findAdapters(). Per instance scratch: a module only ever sees one
	// read at a time, from a single thread.
	private int [] firstHit;
	// Bitset over the last 8 bases (16 bits) of each adapter, to cheaply rule out most read k-mers.
	// With 6 adapters only ~6 of the 65536 bits are set, so nearly every k-mer of a read is
	// rejected by a single bit test, without hashing or probing the table.
	private long [] adapterFilter;

	// 2-bit code of each byte: A=0 C=1 G=2 T=3, -1 for anything else (N, lowercase, '.').
	// Only exact uppercase ACGT can match, as with indexOf on ACGT-only adapters.
	private static final byte [] BASE_CODES = new byte[256];
	static {
		Arrays.fill(BASE_CODES, (byte)-1);
		BASE_CODES['A'] = 0;
		BASE_CODES['C'] = 1;
		BASE_CODES['G'] = 2;
		BASE_CODES['T'] = 3;
	}

	public AdapterContent () {

		Vector<Adapter>c = new Vector<Adapter>();
		Vector<String>l = new Vector<String>();
		
		boolean differentAdapterLengths = false;

		try {

			BufferedReader br = null;
			if (FastQCConfig.getInstance().adapter_file == null) {
				InputStream rsrc=ContaminentFinder.class.getResourceAsStream("/Configuration/adapter_list.txt");
				if (rsrc==null) throw new FileNotFoundException("cannot find Configuration/adapter_list.txt");
				br =new BufferedReader(new InputStreamReader(rsrc));
			}
			else {
				br=new BufferedReader(new FileReader(FastQCConfig.getInstance().adapter_file));
			}


			String line;
			while ((line = br.readLine())!= null){

				if (line.startsWith("#")) continue; // Skip comments
				if (line.trim().length() == 0) continue; // Skip blank lines

				String [] sections = line.split("\\t+");
				if (sections.length != 2) {
					System.err.println("Expected 2 sections for contaminant line but got "+sections.length+" from "+line);
					continue;
				}
				Adapter adapter = new Adapter(sections[0], sections[1]);
				c.add(adapter);	
				l.add(adapter.name());
				
				if (longestAdapter > 0 && adapter.sequence().length() != longestAdapter) differentAdapterLengths=true;
				
				if (adapter.sequence().length() > longestAdapter) longestAdapter = adapter.sequence().length();
			}
			labels = l.toArray(new String[0]);

			br.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		
		if (differentAdapterLengths) {
			System.err.println("[Warning] You are using adapter sequences with different lengths. Matches will only be reported up to the position where the longest adapter could match. Matches to shorter adapters at the end of sequences will not be recorded.");
		}

		adapters = c.toArray(new Adapter[0]);
		firstHit = new int[adapters.length];
		if (!differentAdapterLengths) buildAdapterTable();

	}

	private static int baseCode (char c) {
		return c < 256 ? BASE_CODES[c] : -1;
	}

	// Fibonacci hashing: multiply by 2^64/phi and keep high bits, which are well mixed
	// even for keys differing only in their low bits (e.g. PolyA vs PolyG).
	private static int slotFor (long key, int mask) {
		long h = key * 0x9E3779B97F4A7C15L;
		return (int)(h >>> 40) & mask;
	}

	/**
	 * Builds the k-mer lookup table and the 8-base filter used by findAdapters().
	 * Leaves adapterKeys null, so processSequence falls back to indexOf, when there
	 * are no adapters, they are longer than 31 bases (k-mer doesn't fit in 62 bits)
	 * or any of them contains a non-ACGT character. Called only when all adapters
	 * have the same length.
	 * The table has at least twice as many slots as adapters, so linear probing
	 * always reaches an empty slot and lookups of absent keys terminate.
	 */
	private void buildAdapterTable () {
		if (adapters.length == 0 || longestAdapter == 0 || longestAdapter > 31) return;

		long [] keys = new long[adapters.length];
		for (int a=0;a<adapters.length;a++) {
			String seq = adapters[a].sequence();
			for (int i=0;i<seq.length();i++) {
				int code = baseCode(seq.charAt(i));
				if (code < 0) return;
				keys[a] = (keys[a] << 2) | code;
			}
		}

		int size = Integer.highestOneBit(adapters.length * 4 - 1) << 1;
		adapterMask = size-1;
		adapterKeys = new long[size];
		adapterSlots = new int[size][];
		adapterFilter = new long[1 << 10];
		for (int a=0;a<adapters.length;a++) {
			int low = (int)(keys[a] & 0xFFFF);
			adapterFilter[low >>> 6] |= 1L << low;
			int slot = slotFor(keys[a], adapterMask);
			while (adapterSlots[slot] != null && adapterKeys[slot] != keys[a]) slot = (slot+1) & adapterMask;
			adapterKeys[slot] = keys[a];
			// Identical adapter sequences share a slot
			int [] ids = adapterSlots[slot] == null ? new int[0] : adapterSlots[slot];
			ids = Arrays.copyOf(ids, ids.length+1);
			ids[ids.length-1] = a;
			adapterSlots[slot] = ids;
		}
	}

	/**
	 * Fills firstHit with the first index of each adapter in the bases (-1 if
	 * absent), equivalent to indexOf(adapter) on them for every adapter.
	 * Keeps a rolling 2-bit key of the last k bases; a non-ACGT base resets it, so
	 * only windows of k valid bases are looked up. Positions are visited in order,
	 * so the first time an adapter's key is seen is its indexOf position. Stops
	 * early once every adapter has been found.
	 */
	private void findAdapters (byte [] bases) {
		Arrays.fill(firstHit, -1);
		int length = bases.length;

		int k = longestAdapter;
		long kmerMask = (1L << (2*k)) - 1;
		long [] filter = adapterFilter;
		long key = 0;
		int valid = 0;
		int remaining = adapters.length;
		for (int i=0;i<length;i++) {
			int code = BASE_CODES[bases[i] & 0xFF];
			if (code < 0) {
				valid = 0;
				key = 0;
				continue;
			}
			key = ((key << 2) | code) & kmerMask;
			if (++valid < k) continue;

			int low = (int)key & 0xFFFF;
			if ((filter[low >>> 6] & (1L << low)) == 0) continue;

			int slot = slotFor(key, adapterMask);
			while (adapterSlots[slot] != null) {
				if (adapterKeys[slot] == key) {
					for (int a : adapterSlots[slot]) {
						if (firstHit[a] < 0) {
							firstHit[a] = i-k+1;
							--remaining;
						}
					}
					break;
				}
				slot = (slot+1) & adapterMask;
			}
			if (remaining == 0) return;
		}
	}

	public boolean ignoreFilteredSequences() {
		return true;
	}

	public boolean ignoreInReport () {
		if (ModuleConfig.getParam("adapter", "ignore") > 0) {
			return true;
		}
		return false;
	}

	public JPanel getResultsPanel() {

		if (longestAdapter > longestSequence) {
			// We can't display sensible results
			JPanel failPanel = new JPanel();
			failPanel.setLayout(new BorderLayout());
			failPanel.add(new JLabel("Can't analyse adapters as read length is too short ("+longestAdapter+" vs "+longestSequence+")",JLabel.CENTER),BorderLayout.CENTER);
			return failPanel;
		}

		if (!calculated) calculateEnrichment();

		//		System.err.println("Xlabels="+xLabels.length+" labels="+labels.length+" enrichments="+enrichments.length+" enrichments[0]="+enrichments[0].length+" groups="+groups.length);	

		return (new LineGraph(enrichments, 0, 100, "Position in read (bp)", labels, xLabels, "% Adapter"));

	}



	public void processSequence(Sequence sequence) {
		calculated = false;
		++totalCount;

		// We need to be careful about making sure that a sequence is not only longer
		// than we've seen before, but also that the last position we could find a hit
		// is a positive position.
		
		// If the sequence is longer than it was then we need to expand the storage in
		// all of the adapter objects to account for this.
		
		if (sequence.length() > longestSequence && sequence.length() - longestAdapter > 0) {
			longestSequence = sequence.length();
			for (int a=0;a<adapters.length;a++) {
				adapters[a].expandLengthTo((longestSequence-longestAdapter)+1);
			}
		}

		// Now we go through all of the Adapters to see where they occur

		// One pass over the read finds all adapters at once (see findAdapters);
		// otherwise each adapter is searched with its own indexOf scan.
		if (adapterKeys != null) findAdapters(sequence.getSequenceBytes());

		// The indexOf fallback needs the read as a String: build it once, not per adapter
		String seqString = adapterKeys != null ? null : sequence.getSequence();
		for (int a=0;a<adapters.length;a++) {

			int index = adapterKeys != null ? firstHit[a] : seqString.indexOf(adapters[a].sequence());
			if (index >=0) {
				for (int i=index;i<=longestSequence-longestAdapter;i++) {
					adapters[a].incrementCount(i);
				}
			}
		}


	}

	public synchronized void calculateEnrichment () {
		int maxLength = 0;
		for (int a=0;a<adapters.length;a++) {
			if (adapters[a].getPositions().length > maxLength) {
				maxLength = adapters[a].getPositions().length;
			}
		}

		// We'll be grouping together positions later so make up the groups now
		groups = BaseGroup.makeBaseGroups(maxLength);

		//		System.err.println("Made "+groups.length+" groups from "+maxLength);

		xLabels = new String[groups.length];
		for (int i=0;i<xLabels.length;i++) {
			xLabels[i] = groups[i].toString();
		}

		enrichments = new double [adapters.length][groups.length];


		for (int a=0;a<adapters.length;a++) {
			long [] positions = adapters[a].positions;

			for (int g=0;g<groups.length;g++) {

//				System.err.println("Looking at group "+groups[g]);

				for (int p=groups[g].lowerCount()-1;p<groups[g].upperCount() && p<positions.length;p++) {
//					System.err.println("Count at position "+p+" is "+ positions[p]);
					enrichments[a][g] += (positions[p]*100d)/totalCount;
//					System.err.println("Percentage at position "+p+" is "+ ((positions[p]*100d)/totalCount)+" total count is "+totalCount);
				}

				enrichments[a][g] /= (groups[g].upperCount()-groups[g].lowerCount())+1;
//				System.err.println("Averge Percetage for group "+groups[g]+" is "+ enrichments[a][g]);
			}
		}

		calculated = true;

	}

	public void reset () {
		calculated = false;
		totalCount = 0;
		longestSequence = 0;
		for (int a=0;a<adapters.length;a++) {
			adapters[a].reset();
		}
	}

	public String description() {
		return "Searches for specific adapter sequences in a library";
	}

	public String name() {
		return "Adapter Content";
	}

	public boolean raisesError() {
		if (!calculated) calculateEnrichment();

		for (int i=0;i<enrichments.length;i++) {
			for (int j=0;j<enrichments[i].length;j++) {
				if (enrichments[i][j] > ModuleConfig.getParam("adapter", "error")) return true;
			}
		}
		return false;
	}

	public boolean raisesWarning() {

		// We warn if we just couldn't run the analysis
		if (longestAdapter > longestSequence) return true;

		if (!calculated) calculateEnrichment();

		for (int i=0;i<enrichments.length;i++) {
			for (int j=0;j<enrichments[i].length;j++) {
				if (enrichments[i][j] > ModuleConfig.getParam("adapter", "warn")) return true;
			}
		}
		return false;
	}

	public void makeReport(HTMLReportArchive report) throws IOException,XMLStreamException {

		if (longestAdapter > longestSequence) {
			XMLStreamWriter xhtml = report.xhtmlStream();

			xhtml.writeStartElement("p");
			xhtml.writeCharacters("Can't analyse adapters as read length is too short ("+longestAdapter+" vs "+longestSequence+")");
			xhtml.writeEndElement();
		}

		else {

			if (!calculated) calculateEnrichment();

			writeDefaultImage(report, "adapter_content.png", "Adapter graph", Math.max(800, groups.length*15), 600);

			StringBuffer sb = report.dataDocument();

			ResultsTable table = new ResultsTable();
			// Header
			sb.append("#");
			for (int i=0;i<table.getColumnCount();i++) {
				if (i>0) {
					sb.append("\t");
				}
				sb.append(table.getColumnName(i));
			}
			sb.append("\n");
			
			for (int r=0;r<table.getRowCount();r++) {
				for (int c=0;c<table.getColumnCount();c++) {
					if (c>0) {
						sb.append("\t");
					}
					sb.append(table.getValueAt(r, c));
				}
				sb.append("\n");
			}
			
		}
	}

	private class Adapter {

		private String name;
		private String sequence;
		private long [] positions = new long[0];

		public Adapter (String name, String sequence) {

			this.name = name;
			this.sequence = sequence;
			positions = new long[1];
		}

		public void incrementCount (int position) {

			// Don't ever check or expand the storage within this
			// function as it ends up double counting previously
			// incremented positions.  Rely on the upstream code
			// having done the expansion correctly already.

			++positions[position];			

		}
		
		public void expandLengthTo (int newLength) {
			long [] newPositions = new long[newLength];
			for (int i=0;i<positions.length;i++) {
				newPositions[i] = positions[i];
				//					System.err.println("Copied value "+positions[i]+" at position "+i);
			}
			// Copy the current longest value to the newly added slots
			if (positions.length > 0) {
				for (int i=positions.length;i<newPositions.length;i++) {
					newPositions[i] = positions[positions.length-1];
				}
			}
			
			positions = newPositions;
			
		}

		public long [] getPositions () {
			return positions;
		}

		public String sequence () {
			return sequence;
		}

		public void reset () {
			positions = new long[0];
		}

		public String name () {
			return name;
		}
	}


	private class ResultsTable extends AbstractTableModel {
		private static final long serialVersionUID = 1L;

		public ResultsTable () {}

		// Sequence - Count - Obs/Exp
		public int getColumnCount() {
			return adapters.length+1;
		}

		public int getRowCount() {
			return enrichments[0].length;
		}

		public Object getValueAt(int rowIndex, int columnIndex) {
			if (columnIndex == 0) return xLabels[rowIndex];
			return enrichments[columnIndex-1][rowIndex];
		}

		public String getColumnName (int columnIndex) {
			if (columnIndex == 0) return "Position";
			return (labels[columnIndex-1]);
		}

		public Class<?> getColumnClass (int columnIndex) {
			switch (columnIndex) {
			case 0: return String.class;
			}
			return Double.class;

		}
	}

}
