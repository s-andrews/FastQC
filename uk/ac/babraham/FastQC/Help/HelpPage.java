/**
 * Copyright 2009-15 Simon Andrews
 *
 *    This file is part of SeqMonk.
 *
 *    SeqMonk is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    SeqMonk is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with SeqMonk; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
package uk.ac.babraham.FastQC.Help;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Vector;

import javax.swing.tree.DefaultMutableTreeNode;

/**
 * The Class HelpPage represents a single page of information in 
 * the help system
 */
public class HelpPage extends DefaultMutableTreeNode {
	
	/** The file. */
	private File file;
	
	/** The name. */
	private String name;
	
	/**
	 * Instantiates a new help page.
	 * 
	 * @param file the file
	 */
	public HelpPage (File file) {
		this.file = file;
		name = file.getName();
		name = this.name.replaceFirst("\\.[hH][tT][mM][lL]?$", "");
		
		String [] nameSections = name.split(" ");
		if (nameSections.length > 1) {
			// We have two sections so check if the first is just integers
			// separated by dots.  If it is then we can lose it.
			String [] numbers = nameSections[0].split("\\.");
			for (int n=0;n<numbers.length;n++) {
				try {
					Integer.parseInt(numbers[n]);
				}
				catch(NumberFormatException nfe) {
					return; // We don't want to chop this off.
				}
			}
			
			// If we get here then we want to chop the first part
			// of the name off
			StringBuffer sb = new StringBuffer(nameSections[1]);
			for (int s=2;s<nameSections.length;s++) {
				sb.append(" ");
				sb.append(nameSections[s]);
			}
			name = sb.toString();
			
		}
		
	}
	
	/**
	 * Contains string.
	 * 
	 * @param searchTerm the search term
	 * @param hits the hits
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	@SuppressWarnings("unchecked")
	public void containsString (String searchTerm, Vector<HelpPage>hits) throws IOException {
				
		// Since this will be part of a search thread then take a quick
		// break in case we're trying to do anything else.
		try {
			Thread.sleep(10);
		} 
		catch (InterruptedException e) {}
		
		// We don't want to be trying to open directories
		if (isLeaf()) {
			BufferedReader br = new BufferedReader(new FileReader(file));
			searchTerm = searchTerm.toLowerCase();
			String line;
			while ((line = br.readLine()) != null) {
//				System.out.println("Read line "+line);
				if (line.toLowerCase().indexOf(searchTerm)!=-1) {
					hits.add(this);
					break;
				}
			}
		}

		// Extend the search to our children
		Enumeration kids = children();
		while (kids.hasMoreElements()) {
			Object node = kids.nextElement();
			if (node instanceof HelpPage) {
				((HelpPage)node).containsString(searchTerm, hits);
			}
		}
	}
	
	/* (non-Javadoc)
	 * @see javax.swing.tree.DefaultMutableTreeNode#toString()
	 */
	public String toString () {
		return name;
	}
	
	/* (non-Javadoc)
	 * @see javax.swing.tree.DefaultMutableTreeNode#isLeaf()
	 */
	public boolean isLeaf() {
		if (file.isDirectory()) return false;
		else return true;
	}
	
	/**
	 * Gets the file.
	 * 
	 * @return the file
	 */
	public File getFile () {
		return file;
	}

}
