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

public class NanoporeBasename {

	/**
	 * This method finds the core name from an ONT fast5 file.  It strips off the
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
	
	public static String getNanoporeBasename (String originalName) throws NameFormatException {
		
		// Files from nanopores look like: Computer_Samplename_number_chXXX_fileXXX_strand.fast5
		// We need to reduce this to Computer_Samplename_number
		
		// Some more recent files have names which are just: Computer_Samplename_number.fast5 so 
		// we need to account for those too.
				
		String [] subNames = originalName.replaceAll(".fast5$", "").split("_");
		
		if (subNames.length < 3) {
			throw new NameFormatException();
		}

		String basename = subNames[0]+"_"+subNames[1]+"_"+subNames[2];
		
		System.err.println("Basename is "+basename);
		
		return basename;
		
	}
	
	public static File [][] getNanoporeGroups (File [] files) {
		Hashtable<String, Vector<File>> fileBases = new Hashtable<String, Vector<File>>();
		
		for (int f=0;f<files.length;f++) {

			
			if (files[f].getName().contains("muxscan")) continue; // Control files not containing real data.
			
			// If a file forms part of a nanopore group then put it into that
			// group.
			try {
				String baseName = NanoporeBasename.getNanoporeBasename(files[f].getName());
				if (! fileBases.containsKey(baseName)) {
					fileBases.put(baseName,new Vector<File>());
				}
				fileBases.get(baseName).add(files[f]);

			}
			
			// If the file name doesn't appear to be part of a nanopore group
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
