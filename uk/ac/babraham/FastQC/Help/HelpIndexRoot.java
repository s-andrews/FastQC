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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Vector;

import javax.swing.tree.DefaultMutableTreeNode;


/**
 * The Class HelpIndexRoot is the root node of the tree of help files.
 */
public class HelpIndexRoot extends DefaultMutableTreeNode {
	
	/** The fs. */
	public FileSorter fs = new FileSorter();
	
	/**
	 * Instantiates a new help index root.
	 * 
	 * @param startingLocation the starting location
	 */
	public HelpIndexRoot (File startingLocation) {
		super("Help Contents");
		
		if (!startingLocation.exists() || !startingLocation.isDirectory()) {
			throw new IllegalArgumentException("Couldn't find help file directory at '"+startingLocation.getAbsolutePath()+"'");
		}
		
		addSubfiles(startingLocation, this);
	}
	
	/**
	 * Adds the subfiles.
	 * 
	 * @param directory the directory
	 * @param node the node
	 */
	private void addSubfiles (File directory, DefaultMutableTreeNode node) {
		File [] files = directory.listFiles();
		
		Arrays.sort(files,fs);
		for (int f=0;f<files.length;f++) {
			if (files[f].isDirectory()) {
				HelpPage h = new HelpPage(files[f]);
				node.add(h);
				addSubfiles(files[f], h);
			}
			else if (files[f].getName().toLowerCase().endsWith(".html") || files[f].getName().toLowerCase().endsWith(".htm")){
				HelpPage h = new HelpPage(files[f]);
				node.add(h);
			}
			// Skip files which aren't html (eg images etc)
		}
	}
	
	/**
	 * Find pages for term.
	 * 
	 * @param searchTerm the search term
	 * @return the help page[]
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	@SuppressWarnings("unchecked")
	public HelpPage [] findPagesForTerm (String searchTerm) throws IOException {
		Vector<HelpPage>hits = new Vector<HelpPage>();
				
		Enumeration kids = children();
		while (kids.hasMoreElements()) {
			Object node = kids.nextElement();
			if (node instanceof HelpPage) {
				((HelpPage)node).containsString(searchTerm, hits);
			}
		}
		
		return hits.toArray(new HelpPage[0]);
	}
	

	/**
	 * The Class FileSorter.
	 */
	private class FileSorter implements Comparator<File> {

		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(File f1, File f2) {
			
			// The file names should be preceeded by a series of
			// integers separated by dots (eg 1.2.1).  We therefore
			// split these out to compare the individual sections

			int [] numbers1;
			int [] numbers2;

			
			try {
				numbers1 = getNumberArray(f1);
				numbers2 = getNumberArray(f2);
			}
			catch (NumberFormatException nfe) {
				return f1.getName().compareTo(f2.getName());
			}
			
			int shortest = numbers1.length;
			if (numbers2.length < shortest) shortest = numbers2.length;
			
			for (int i=0;i<shortest;i++) {
				if (numbers1[i] != numbers2[i]) {
					return numbers1[i]-numbers2[i];
				}
			}
			
			// If we get here then the shortest number string wins
			return numbers1.length - numbers2.length;
			
		}
		
		/**
		 * Gets the number array.
		 * 
		 * @param f the f
		 * @return the number array
		 * @throws NumberFormatException the number format exception
		 */
		private int [] getNumberArray (File f) throws NumberFormatException {
			String [] numberStrings = f.getName().split(" ")[0].split("\\.");
			int [] ints = new int [numberStrings.length];
			for (int i=0;i<numberStrings.length;i++) {
				ints[i] = Integer.parseInt(numberStrings[i]);
			}
			
			return ints;
		}
		
	}

	
}
