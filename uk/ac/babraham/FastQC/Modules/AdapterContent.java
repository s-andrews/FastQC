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
	public AdapterContent () {

		Vector<Adapter>c = new Vector<Adapter>();
		Vector<String>l = new Vector<String>();

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
				if (adapter.sequence().length() > longestAdapter) longestAdapter = adapter.sequence().length();
			}
			labels = l.toArray(new String[0]);

			br.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		adapters = c.toArray(new Adapter[0]);

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
		
		if (sequence.getSequence().length() > longestSequence && sequence.getSequence().length() - longestAdapter > 0) {
			longestSequence = sequence.getSequence().length();
			for (int a=0;a<adapters.length;a++) {
				adapters[a].expandLengthTo((longestSequence-longestAdapter)+1);
			}
		}

		// Now we go through all of the Adapters to see where they occur

		for (int a=0;a<adapters.length;a++) {

			int index = sequence.getSequence().indexOf(adapters[a].sequence());
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
