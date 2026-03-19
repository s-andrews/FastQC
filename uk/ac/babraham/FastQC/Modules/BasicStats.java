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

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.xml.stream.XMLStreamException;

import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;
import uk.ac.babraham.FastQC.Sequence.QualityEncoding.PhredEncoding;

public class BasicStats extends AbstractQCModule {

	private String name = null;
	private long actualCount = 0;
	private long filteredCount = 0;
	private int minLength = 0;
	private int maxLength = 0;
	private long totalBases = 0;
	private long gCount = 0;
	private long cCount = 0;
	private long aCount = 0;
	private long tCount = 0;
	@SuppressWarnings("unused")
	private long nCount = 0;
	private char lowestChar = 126;
	private String fileType = null;
	
	public String description() {
		return "Calculates some basic statistics about the file";
	}
	
	public boolean ignoreFilteredSequences() {
		return false;
	}

	public JPanel getResultsPanel() {
		JPanel returnPanel = new JPanel();
		returnPanel.setLayout(new BorderLayout());
		returnPanel.add(new JLabel("Basic sequence stats",JLabel.CENTER),BorderLayout.NORTH);
		
		TableModel model = new ResultsTable();
		returnPanel.add(new JScrollPane(new JTable(model)),BorderLayout.CENTER);
		
		return returnPanel;
	
	}
	
	public void reset () {
		minLength = 0;
		maxLength = 0;
		gCount = 0;
		cCount = 0;
		aCount = 0;
		tCount = 0;
		nCount = 0;
	}

	public String name() {
		return "Basic Statistics";
	}
	
	public void setFileName (String name) {
		this.name = name;
		
		this.name = this.name.replaceFirst("stdin:", "");
	}

	public void processSequence(Sequence sequence) {

		if (name == null) setFileName(sequence.file().name());
		
		// If this is a filtered sequence we simply count it and move on.
		if (sequence.isFiltered()) {
			filteredCount++;
			return;
		}
		
		actualCount++;
		
		totalBases += sequence.getSequence().length();
		
		if (fileType == null) {
			if (sequence.getColorspace() != null) {
				fileType = "Colorspace converted to bases";
			}
			else {
				fileType = "Conventional base calls";
			}
		}
		
		if (actualCount == 1) {
			minLength = sequence.getSequence().length();
			maxLength = sequence.getSequence().length();
		}
		else {
			if (sequence.getSequence().length() < minLength) minLength = sequence.getSequence().length();
			if (sequence.getSequence().length() > maxLength) maxLength = sequence.getSequence().length();
		}

		char [] chars = sequence.getSequence().toCharArray();
		for (int c=0;c<chars.length;c++) {			
			switch (chars[c]) {
				case 'G': ++gCount;break;
				case 'A': ++aCount;break;
				case 'T': ++tCount;break;
				case 'C': ++cCount;break;
				case 'N': ++nCount;break;			
			}
		}
		
		chars = sequence.getQualityString().toCharArray();
		for (int c=0;c<chars.length;c++) {
			if (chars[c] < lowestChar) {
				lowestChar = chars[c];
			}
		}
	}
	
	public boolean raisesError() {
		return false;
	}

	public boolean raisesWarning() {
		return false;
	}
	
	public boolean ignoreInReport () {
		return false;
	}
	
	public static String formatLength (long originalLength) {
		
		double length = originalLength;
		
		String unit = " bp";

		if (length >= 1000000000) {
			length /= 1000000000;
			unit = " Gbp";
		}

		else if (length >= 1000000) {
			length /= 1000000;
			unit = " Mbp";
		}
		else if (length >= 1000) {
			length /=1000;
			unit = " kbp";
		}

		
		String rawLength = ""+length;
		char [] chars = rawLength.toCharArray();
		
		int lastIndex = 0;
		
		// Go through until we find a dot (if there is one)
		for (int i=0;i<chars.length;i++) {
			lastIndex = i;
			if (chars[i] == '.') break;
		}
		
		// We keep the next char as well if they are non
		// zero numbers
		
		if (lastIndex+1 < chars.length && chars[lastIndex+1] != '0') {
			lastIndex+=1;
		}
		else if (lastIndex > 0 && chars[lastIndex] == '.') {
			lastIndex -= 1; // Lose the dot if its the last character
		}

		char [] finalChars = new char[lastIndex+1];
		for (int i=0;i<=lastIndex;i++) {
			finalChars[i] = chars[i];
		}
		
		return new String(finalChars)+unit;
	}


	public void makeReport(HTMLReportArchive report) throws XMLStreamException,IOException {
		super.writeTable(report, new ResultsTable());
	}

	@SuppressWarnings("serial")
	private class ResultsTable extends AbstractTableModel {
				
		private String [] rowNames = new String [] {
				"Filename",
				"File type",
				"Encoding",
				"Total Sequences",
				"Total Bases",
				"Sequences flagged as poor quality",
				"Sequence length",
				"%GC",
		};		
		
		// Sequence - Count - Percentage
		public int getColumnCount() {
			return 2;
		}
	
		public int getRowCount() {
			return rowNames.length;
		}
	
		public Object getValueAt(int rowIndex, int columnIndex) {
			switch (columnIndex) {
				case 0: return rowNames[rowIndex];
				case 1:
					switch (rowIndex) {
					case 0 : return name;
					case 1 : return fileType;
					case 2 : return PhredEncoding.getFastQEncodingOffset(lowestChar);
					case 3 : return ""+actualCount;
					case 4 : return BasicStats.formatLength(totalBases);
					case 5 : return ""+filteredCount;
					case 6 :
						if (minLength == maxLength) {
							return ""+minLength;
						}
						else {
							return minLength+"-"+maxLength;
						}
						
						
					case 7 : 
						if (aCount+tCount+gCount+cCount > 0) {
							return ""+(((gCount+cCount)*100)/(aCount+tCount+gCount+cCount));
						}
						else {
							return 0;
						}
					
					}
			}
			return null;
		}
		
		public String getColumnName (int columnIndex) {
			switch (columnIndex) {
				case 0: return "Measure";
				case 1: return "Value";
			}
			return null;
		}
		
		public Class<?> getColumnClass (int columnIndex) {
			switch (columnIndex) {
			case 0: return String.class;
			case 1: return String.class;
		}
		return null;
			
		}
	}
	

}
