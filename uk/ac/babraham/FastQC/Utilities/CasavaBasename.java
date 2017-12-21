/**
 * Copyright Copyright 2011-17 Simon Andrews
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

package uk.ac.babraham.FastQC.Utilities;

import java.io.File;
import java.util.Hashtable;
import java.util.Vector;

public class CasavaBasename {

	/**
	 * This method finds the core name from a CASAVA 1.8 fastq file.  It strips off the
	 * part which indicates that this file is one of a set and returns the base name with
	 * this part removed.
	 * 
	 * If the filename does not conform to standard CASAVA naming then a NameFormatException
	 * is thrown.
	 * 
	 * @param originalName
	 * @return
	 * @throws NameFormatException
	 */
	
	public static String getCasavaBasename (String originalName) throws NameFormatException {
		
		// Find the base name.  We need to remove the 123 numbers from
		// files of the form:
		//
		// anyold_text_123.fastq.gz
		//
		// where the base name will be
		//
		// anyoldtext.fastq.gz
		
		// The file must usually end with .fastq.gz, but you can tell cassava to not
		// compress, in which case the .gz is missing.

		
		if (originalName.endsWith(".fastq.gz")) {
			
			// They must have an _ 13 chars before the end
			if (originalName.substring(originalName.length()-13, originalName.length()-12).equals("_")) {
				
				// They must have numbers for the 3 positions before .fastq
				try {
					Integer.parseInt(originalName.substring(originalName.length()-12, originalName.length()-9));
					
					// If we get here then everything is OK to use the base name from this file
					String baseName = originalName.substring(0,originalName.length()-13)+".fastq.gz";
					return baseName;
				}
				catch (NumberFormatException nfe) {}
			}
		}

		
		else if (originalName.endsWith(".fastq")) {
			
			// They must have an _ 10 chars before the end
			if (originalName.substring(originalName.length()-10, originalName.length()-9).equals("_")) {
				
				// They must have numbers for the 3 positions before .fastq
				try {
					Integer.parseInt(originalName.substring(originalName.length()-9, originalName.length()-6));
					
					// If we get here then everything is OK to use the base name from this file
					String baseName = originalName.substring(0,originalName.length()-10)+".fastq";
					return baseName;
				}
				catch (NumberFormatException nfe) {}
			}
		}

		
		throw new NameFormatException();
	}
	
	public static File [][] getCasavaGroups (File [] files) {
		Hashtable<String, Vector<File>> fileBases = new Hashtable<String, Vector<File>>();
		
		for (int f=0;f<files.length;f++) {

			// If a file forms part of a CASAVA group then put it into that
			// group.
			try {
				String baseName = CasavaBasename.getCasavaBasename(files[f].getName());
				if (! fileBases.containsKey(baseName)) {
					fileBases.put(baseName,new Vector<File>());
				}
				fileBases.get(baseName).add(files[f]);

			}
			
			// If the file name doesn't appear to be part of a CASAVA group
			// then add it as a singleton
			catch (NameFormatException nfe) {
				
				System.err.println("File '"+files[f].getName()+"' didn't look like part of a CASAVA group");
				Vector<File> newVector = new Vector<File>();
				newVector.add(files[f]);
				fileBases.put(files[f].getName(), newVector);				
			}
			
		}
		
		String [] baseNames = fileBases.keySet().toArray(new String [0]);
		
		File [][] fileGroups = new File[baseNames.length][];
		
		for (int i=0;i<baseNames.length;i++) {
			fileGroups[i] = fileBases.get(baseNames[i]).toArray(new File[0]);
		}
		
		return fileGroups;		
	}
	
}
