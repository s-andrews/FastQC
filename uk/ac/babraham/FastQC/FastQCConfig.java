/**
 * Copyright Copyright 2012-17 Simon Andrews
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
package uk.ac.babraham.FastQC;

import java.io.File;

public class FastQCConfig {
	
	private static FastQCConfig instance = new FastQCConfig();
	public boolean nogroup = false;
	public boolean expgroup = false;
	public boolean quiet = false;
	public boolean show_version = false;
	public Integer kmer_size = null;
	public Integer threads = null;
	public boolean showUpdates = true;
	public File output_dir = null;
	public boolean casava = false;
	public boolean nano = false;
	public boolean nofilter = false;
	public Boolean do_unzip = null;
	public boolean delete_after_unzip = false;
	public String lineSeparator = System.getProperty("line.separator");
	public String sequence_format = null;
	public File contaminant_file = null;
	public File adapter_file = null;
	public File limits_file = null;
	public int minLength = 0;
	public int dupLength = 0;
	public boolean svg_output = false;

	private FastQCConfig () {
		
		// Output dir
		if (System.getProperty("fastqc.output_dir") != null) {
			output_dir = new File(System.getProperty("fastqc.output_dir"));
			if (!(output_dir.exists() && output_dir.canWrite())) {
				throw new IllegalArgumentException("Output dir "+output_dir+" doesn't exist or isn't writeable");
			}
		}
		
		// Contaminant file
		if (System.getProperty("fastqc.contaminant_file") != null) {
			contaminant_file = new File(System.getProperty("fastqc.contaminant_file"));
			if (!(contaminant_file.exists() && contaminant_file.canRead())) {
				throw new IllegalArgumentException("Contaminant file "+contaminant_file+" doesn't exist or can't be read");
			}
		}

		// Adapter file
		if (System.getProperty("fastqc.adapter_file") != null) {
			adapter_file = new File(System.getProperty("fastqc.adapter_file"));
			if (!(adapter_file.exists() && adapter_file.canRead())) {
				throw new IllegalArgumentException("Adapter file "+adapter_file+" doesn't exist or can't be read");
			}
		}

		// Limits file
		if (System.getProperty("fastqc.limits_file") != null) {
			limits_file = new File(System.getProperty("fastqc.limits_file"));
			if (!(limits_file.exists() && limits_file.canRead())) {
				throw new IllegalArgumentException("Limits file "+limits_file+" doesn't exist or can't be read");
			}
		}
		
		// Threads
		if (System.getProperty("fastqc.threads") != null) {
			threads = Integer.parseInt(System.getProperty("fastqc.threads"));
			if (threads < 1) {
				throw new IllegalArgumentException("Number of threads must be >= 1");
			}
		}
		
		// Kmer size
		if (System.getProperty("fastqc.kmer_size") != null) {
			kmer_size = Integer.parseInt(System.getProperty("fastqc.kmer_size"));
		}

		
		// Min length
		if (System.getProperty("fastqc.min_length") != null) {
			minLength = Integer.parseInt(System.getProperty("fastqc.min_length"));
		}

		// Dup length
		if (System.getProperty("fastqc.dup_length") != null) {
			dupLength = Integer.parseInt(System.getProperty("fastqc.dup_length"));
		}

		
		// Quiet
		if (System.getProperty("fastqc.quiet") != null && System.getProperty("fastqc.quiet").equals("true")) {
			quiet = true;
		}
		
		// Casava
		if (System.getProperty("fastqc.casava") != null && System.getProperty("fastqc.casava").equals("true")) {
			casava = true;
		}

		// Nanopore
		if (System.getProperty("fastqc.nano") != null && System.getProperty("fastqc.nano").equals("true")) {
			nano = true;
		}

		// SVG
		if (System.getProperty("fastqc.svg") != null && System.getProperty("fastqc.svg").equals("true")) {
			svg_output = true;
		}

		
		// Nofilter
		if (System.getProperty("fastqc.nofilter") != null && System.getProperty("fastqc.nofilter").equals("true")) {
			nofilter = true;
		}

		
		// No group
		if (System.getProperty("fastqc.nogroup") != null && System.getProperty("fastqc.nogroup").equals("true")) {
			nogroup = true;
		}

		// Exponential group
		if (System.getProperty("fastqc.expgroup") != null && System.getProperty("fastqc.expgroup").equals("true")) {
			expgroup = true;
		}

		// Unzip
		if (System.getProperty("fastqc.unzip") != null && System.getProperty("fastqc.unzip").equals("true")) {
			do_unzip = true;
			
			if (System.getProperty("fastqc.delete") != null && System.getProperty("fastqc.delete").equals("true")) {
				delete_after_unzip = true;
			}
		}
			
		
		// Sequence Format
		if (System.getProperty("fastqc.sequence_format") != null) {
			setSequenceFormat(System.getProperty("fastqc.sequence_format"));
		}
		
	};
	
	public void setSequenceFormat (String sequenceFormat) {
		if (sequenceFormat.equals("fastq") || sequenceFormat.equals("sam") || sequenceFormat.equals("bam") || sequenceFormat.equals("sam_mapped") || sequenceFormat.equals("bam_mapped")) {
			sequence_format = sequenceFormat;
		}
		else {
			throw new IllegalArgumentException("Sequence format '"+sequenceFormat+"' wasn't recognised");
		}
	}
	
	public void setCasavaMode (boolean casava) {
		this.casava = casava;
	}

	public static FastQCConfig getInstance() {
		return instance;
	}


}
