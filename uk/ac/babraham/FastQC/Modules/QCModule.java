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

import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.Sequence;

public interface QCModule {

	public void processSequence(Sequence sequence);

	public JPanel getResultsPanel();
	
	public String name ();
	
	public String description ();
	
	public void reset ();
	
	public boolean raisesError();
	
	public boolean raisesWarning();
	
	public boolean ignoreFilteredSequences();
	
	/**
	 * Allows you to say that this module shouldn't be included in the final report.
	 * Useful for modules which have a use under some circumstances but not others.
	 * @return
	 */
	public boolean ignoreInReport();

	public void makeReport(HTMLReportArchive report) throws XMLStreamException,IOException;
	
	
}
