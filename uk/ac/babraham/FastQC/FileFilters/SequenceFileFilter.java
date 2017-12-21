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
package uk.ac.babraham.FastQC.FileFilters;

import java.io.File;

import javax.swing.filechooser.FileFilter;

public class SequenceFileFilter extends FileFilter {

	public boolean accept(File f) {
		if (f.isDirectory() 
				|| f.getName().toLowerCase().endsWith(".txt.gz") 
				|| f.getName().toLowerCase().endsWith(".fastq.gz") 
				|| f.getName().toLowerCase().endsWith(".fq.gz") 
				|| f.getName().toLowerCase().endsWith(".fq") 
				|| f.getName().toLowerCase().endsWith(".txt.bz2") 
				|| f.getName().toLowerCase().endsWith(".fastq.bz2")
				|| f.getName().toLowerCase().endsWith(".txt") 
				|| f.getName().toLowerCase().endsWith(".fastq") 
				|| f.getName().toLowerCase().endsWith(".bam") 
				|| f.getName().toLowerCase().endsWith(".sam")
				|| f.getName().toLowerCase().endsWith(".compact-reads")
				|| f.getName().toLowerCase().endsWith(".goby")
				
		) {
			return true;
		}
		else {
			return false;
		}
	}

	public String getDescription() {
		return "Sequence Files";
	}

}
