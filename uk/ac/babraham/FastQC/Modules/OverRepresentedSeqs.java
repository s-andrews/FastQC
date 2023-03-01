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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;
import uk.ac.babraham.FastQC.Sequence.Contaminant.ContaminantHit;
import uk.ac.babraham.FastQC.Sequence.Contaminant.ContaminentFinder;

public class OverRepresentedSeqs extends AbstractQCModule {

	protected HashMap<String, Integer>sequences = new HashMap<String, Integer>();
	protected long count = 0;
	private OverrepresentedSeq [] overrepresntedSeqs = null;
	private boolean calculated = false;
	private boolean frozen = false;
	private DuplicationLevel duplicationModule;
	
	// This is the number of different sequences we want to track
	private final int OBSERVATION_CUTOFF = 100000;
	// This is a count of how many unique sequences we've seen so far
	// so we know when to stop adding them.
	private int uniqueSequenceCount = 0;
	// This was the total count at the point at which we saw our total
	// number of unique sequences, so we know what to correct by when
	// extrapolating to the whole file
	protected long countAtUniqueLimit = 0;
	
	public OverRepresentedSeqs () {
		duplicationModule = new DuplicationLevel(this);
	}
	
	public boolean ignoreInReport () {
		if (ModuleConfig.getParam("overrepresented", "ignore") > 0) {
			return true;
		}
		return false;
	}
	
	public String description() {
		return "Identifies sequences which are overrepresented in the set";
	}
	
	public boolean ignoreFilteredSequences() {
		return true;
	}
	
	public DuplicationLevel duplicationLevelModule () {
		return duplicationModule;
	}

	public JPanel getResultsPanel() {
		JPanel returnPanel = new JPanel();
		returnPanel.setLayout(new BorderLayout());
		returnPanel.add(new JLabel("Overrepresented sequences",JLabel.CENTER),BorderLayout.NORTH);
		
		if (!calculated) getOverrepresentedSeqs();
		
		if (overrepresntedSeqs.length > 0) {
			TableModel model = new ResultsTable(overrepresntedSeqs);
			JTable table = new JTable(model);
			table.setCellSelectionEnabled(true);
			returnPanel.add(new JScrollPane(table),BorderLayout.CENTER);
		}
		else {
			returnPanel.add(new JLabel("There are no overrepresented sequences",JLabel.CENTER),BorderLayout.CENTER);
		}
		
		return returnPanel;
	
	}
	
	public DuplicationLevel getDuplicationLevelModule () {
		return duplicationModule;
	}
	private synchronized void getOverrepresentedSeqs () {

		// If the duplication module hasn't already done
		// its calculation it needs to do it now before
		// we stomp all over the data
		duplicationModule.calculateLevels();
		
		Iterator<String> s = sequences.keySet().iterator();
		List<OverrepresentedSeq>keepers = new ArrayList<OverrepresentedSeq>();
		
		while (s.hasNext()) {
			String seq = s.next();
			double percentage = ((double)sequences.get(seq)/count)*100;
			if (percentage > ModuleConfig.getParam("overrepresented", "warn")) {
				OverrepresentedSeq os = new OverrepresentedSeq(seq, sequences.get(seq), percentage);
				keepers.add(os);
			}
		}
		
		overrepresntedSeqs = keepers.toArray(new OverrepresentedSeq[0]);
		Arrays.sort(overrepresntedSeqs);
		calculated  = true;
		sequences.clear();
		
	}
	
	public void reset () {
		count = 0;
		sequences.clear();
	}

	public String name() {
		return "Overrepresented sequences";
	}

	public void processSequence(Sequence sequence) {
		
		calculated = false;
		
		++count;
		
		// Since we rely on identity to match sequences we can't trust really long
		// sequences, so anything over 75bp gets truncated to 50bp.
		String seq = sequence.getSequence();

		// For long sequences (above 50bp) we truncate to 50bp so that random
		// base miscalls don't mess things up.  We also allow the user to 
		// specify a shorter truncation length on the command line in the 
		// dup_length option.
		
		if (FastQCConfig.getInstance().dupLength != 0) {
			seq = new String(seq.substring(0, FastQCConfig.getInstance().dupLength));			
		}
		else if (seq.length() > 50) {
			seq = new String(seq.substring(0, 50));
		}
				
		if (sequences.containsKey(seq)) {
			sequences.put(seq, sequences.get(seq)+1);
			
			// We need to increment the count at unique limit just in case
			// we never hit the unique sequence limit, so we need to know 
			// that we'd actually seen all of the sequences.
			if (! frozen) {
				countAtUniqueLimit = count;
			}
		}
		else {
			if (! frozen) {
				sequences.put(seq, 1);
				++uniqueSequenceCount;
				countAtUniqueLimit = count;
				if (uniqueSequenceCount == OBSERVATION_CUTOFF) {
					frozen = true;
				}

			}
		}		
	}
	
	@SuppressWarnings("serial")
	private class ResultsTable extends AbstractTableModel {
		
		private OverrepresentedSeq [] seqs;
		
		public ResultsTable (OverrepresentedSeq [] seqs) {
			this.seqs = seqs;
		}
		
		
		// Sequence - Count - Percentage
		public int getColumnCount() {
			return 4;
		}

		public int getRowCount() {
			return seqs.length;
		}

		public Object getValueAt(int rowIndex, int columnIndex) {
			switch (columnIndex) {
				case 0: return seqs[rowIndex].seq();
				case 1: return seqs[rowIndex].count();
				case 2: return seqs[rowIndex].percentage();
				case 3: return seqs[rowIndex].contaminantHit();
					
			}
			return null;
		}
		
		public String getColumnName (int columnIndex) {
			switch (columnIndex) {
				case 0: return "Sequence";
				case 1: return "Count";
				case 2: return "Percentage";
				case 3: return "Possible Source";
			}
			return null;
		}
		
		public Class<?> getColumnClass (int columnIndex) {
			switch (columnIndex) {
			case 0: return String.class;
			case 1: return Integer.class;
			case 2: return Double.class;
			case 3: return String.class;
		}
		return null;
			
		}
	}
	
	private class OverrepresentedSeq implements Comparable<OverrepresentedSeq>{
		
		private String seq;
		private int count;
		private double percentage;
		private ContaminantHit contaminantHit;
		
		public OverrepresentedSeq (String seq, int count, double percentage) {
			this.seq = seq;
			this.count = count;
			this.percentage = percentage;
			this.contaminantHit = ContaminentFinder.findContaminantHit(seq);
		}
		
		public String seq () {
			return seq;
		}
		
		public int count () {
			return count;
		}
		
		public double percentage () {
			return percentage;
		}
		
		public String contaminantHit () {
			if (contaminantHit == null) {
				return "No Hit";
			}
			else {
				return contaminantHit.toString();
			}
		}

		public int compareTo(OverrepresentedSeq o) {
			return o.count-count;
		}
	}

	public boolean raisesError() {
		if (!calculated) getOverrepresentedSeqs();
		if (overrepresntedSeqs.length>0) {
			if (overrepresntedSeqs[0].percentage > ModuleConfig.getParam("overrepresented", "error")) {
				return true;
			}
		}
		return false;
	}

	public boolean raisesWarning() {
		if (!calculated) getOverrepresentedSeqs();

		if (overrepresntedSeqs.length > 0) return true;
		return false;
	}

	public void makeReport(HTMLReportArchive report) throws IOException,XMLStreamException {
		if (!calculated) getOverrepresentedSeqs();
		ResultsTable table = new ResultsTable(overrepresntedSeqs);
				
		if (overrepresntedSeqs.length == 0)
			{
			XMLStreamWriter w=report.xhtmlStream();
			w.writeStartElement("p");
			w.writeCharacters("No overrepresented sequences");
			w.writeEndElement();
			}
		
		else {
			super.writeTable(report, table);
			}	
		}

}
