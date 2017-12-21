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

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5SimpleReader;

public class Fast5File implements SequenceFile {

	private Sequence nextSequence = null;
	private File file;

	private String name;

	protected Fast5File(File file) throws SequenceFormatException, IOException {
		this.file = file;
		name = file.getName();

		IHDF5SimpleReader reader = HDF5Factory.openForReading(file);
		
		String [] rdfPaths = new String [] {
				"Analyses/Basecall_2D_000/BaseCalled_template/Fastq",
				"Analyses/Basecall_2D_000/BaseCalled_2D/Fastq",
				"Analyses/Basecall_1D_000/BaseCalled_template/Fastq",
				"Analyses/Basecall_1D_000/BaseCalled_1D/Fastq"
		};
		
		boolean foundPath = false;
		for (int r=0;r<rdfPaths.length;r++) {
		
				if (reader.exists(rdfPaths[r])) {
		
					foundPath = true;
					String fastq = reader.readString(rdfPaths[r]);
		
					String [] sections = fastq.split("\\n");
			
					if (sections.length != 4) {
						throw new SequenceFormatException("Didn't get 4 sections from "+fastq);
					}
			
					nextSequence = new Sequence(this, sections[1].toUpperCase(),sections[3], sections[0]);
					break;
				}
		}

		reader.close();

		if (!foundPath) {
			throw new SequenceFormatException("No valid fastq paths found in "+file);
		}
		
	}

	public String name() {
		return name;
	}

	public int getPercentComplete() {
		if (! hasNext()) return 100;

		return 0;		
	}

	public boolean isColorspace() {
		return false;
	}

	public boolean hasNext() {
		return nextSequence != null;
	}

	public Sequence next() throws SequenceFormatException {
		Sequence seq = nextSequence;
		nextSequence = null;
		return seq;
	}

	public void remove() {
		// No action here
	}

	public File getFile() {
		return file;
	}

}
