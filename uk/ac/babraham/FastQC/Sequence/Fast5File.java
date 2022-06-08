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
package uk.ac.babraham.FastQC.Sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleReader;

public class Fast5File implements SequenceFile {

	private Sequence nextSequence = null;
	private File file;
	private String name;
	private IHDF5SimpleReader reader;
	private String [] readPaths = new String[] {""};
	
	private int readPathsIndexPosition = 0;
	
	private String [] rdfPaths = new String [] {
			"Analyses/Basecall_2D_000/BaseCalled_template/Fastq",
			"Analyses/Basecall_2D_000/BaseCalled_2D/Fastq",
			"Analyses/Basecall_1D_000/BaseCalled_template/Fastq",
			"Analyses/Basecall_1D_000/BaseCalled_1D/Fastq"
	};



	protected Fast5File(File file) throws SequenceFormatException, IOException {
		this.file = file;
		name = file.getName();

		reader = HDF5Factory.openForReading(file);
		
		// These files have changed structure over time.  Originally they contained
		// a single read per file where the base of the heirarchy was the read 
		// itself.
		//
		// Later the files moved to having multiple reads per file.  Now there is
		// an additional top level folder per read with the sub-structure being the
		// same as it used to be for the individual reads.
		//
		// We need to account for both of these structures.
		
		
		// See if we can see a bunch of paths starting with "read_" at the top of the
		// heirarchy.  If we can then we substitute the read paths for these.
		
		List<String> topLevelFolders = reader.getGroupMembers("/");
		
		List<String> readFolders = new ArrayList<String>();
		
		for (String folder : topLevelFolders) {
			System.err.println("Looking at "+folder);
			
			if (folder.startsWith("read_")) {
				readFolders.add(folder+"/");
			}
		}
		
		if (readFolders.size() > 0) {
			// We have read folders so we'll replace the default readPaths with
			// the list we made
			
			readPaths = readFolders.toArray(new String[0]);
		}
		
	}

	public String name() {
		return name;
	}

	public int getPercentComplete() {
		return (readPathsIndexPosition*100) / readPaths.length;		
	}

	public boolean isColorspace() {
		return false;
	}

	public boolean hasNext() {
		return readPathsIndexPosition < readPaths.length;
	}

	public Sequence next() throws SequenceFormatException {
		
		for (int r=0;r<rdfPaths.length;r++) {
			
				if (reader.exists(readPaths[readPathsIndexPosition]+rdfPaths[r])) {
		
					String fastq = reader.readString(readPaths[readPathsIndexPosition]+rdfPaths[r]);
		
					String [] sections = fastq.split("\\n");
			
					if (sections.length != 4) {
						throw new SequenceFormatException("Didn't get 4 sections from "+fastq);
					}
			
					Sequence seq = new Sequence(this, sections[1].toUpperCase(),sections[3], sections[0]);
					++readPathsIndexPosition;
					
					if(readPathsIndexPosition >= readPaths.length) {
						reader.close();
					}
					
					return(seq);
				}
		}
		
		throw new SequenceFormatException("No valid fastq paths found in "+file);
	}

	public void remove() {
		// No action here
	}

	public File getFile() {
		return file;
	}

}
