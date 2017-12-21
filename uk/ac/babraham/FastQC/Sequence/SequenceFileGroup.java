/**
 * Copyright Copyright 2013-17 Simon Andrews
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

import uk.ac.babraham.FastQC.Utilities.CasavaBasename;
import uk.ac.babraham.FastQC.Utilities.NameFormatException;

public class SequenceFileGroup implements SequenceFile {

	private File [] files;
	private SequenceFile sequenceFile;
	private File groupFile;
	private int currentIndex = 0;

	public SequenceFileGroup( File [] files) throws IOException, SequenceFormatException {
		this.files = files;
		sequenceFile = SequenceFactory.getSequenceFile(files[0]);
				
		try {
			String baseName = CasavaBasename.getCasavaBasename(sequenceFile.name());
			if (sequenceFile.getFile().getParent() == null) {
				groupFile = new File(baseName);
			} else {
				groupFile = new File(sequenceFile.getFile().getParent() + "/"
						+ baseName);
			}
		} catch (NameFormatException nfe) {
			groupFile = sequenceFile.getFile();
		}
	}

	public File getFile() {
		return groupFile;
	}

	public int getPercentComplete() {
		return ((100 * currentIndex) / files.length)
				+ (sequenceFile.getPercentComplete() / files.length);
	}

	public boolean hasNext() {
		if (sequenceFile.hasNext()) {
			return true;
		} 
		else {
			while (currentIndex < files.length - 1) {
				++currentIndex;
				try {
					sequenceFile = SequenceFactory.getSequenceFile(files[currentIndex]);
				}
				catch (Exception e) {
					e.printStackTrace();
					return false;
				}
				if (sequenceFile.hasNext()) break;
			}
			return sequenceFile.hasNext();
		}
	}

	public boolean isColorspace() {
		return sequenceFile.isColorspace();
	}

	public String name() {
		return groupFile.getName();
	}

	public Sequence next() throws SequenceFormatException {
		return sequenceFile.next();
	}

}
