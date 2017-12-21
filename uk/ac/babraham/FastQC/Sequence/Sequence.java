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

public class Sequence {

	private String sequence;
	private String quality;
	private String id;
	private SequenceFile file;
	private String colorspace;
	private boolean isFiltered;
	
	public Sequence (SequenceFile file,String sequence, String quality, String id) {
		this.id = id;
		this.file = file;
		this.sequence = sequence.toUpperCase();
		this.quality = quality;
		this.colorspace = null;
		this.isFiltered = false;
	}
	
	public Sequence (SequenceFile file,String sequence, String colorspace, String quality, String id) {
		this.id = id;
		this.file = file;
		this.sequence = sequence;
		this.quality = quality;
		this.colorspace = colorspace;
	}
	
	public void setIsFiltered (boolean isFiltered) {
		this.isFiltered = isFiltered;
	}
	
	public boolean isFiltered () {
		return isFiltered;
	}
	
	public SequenceFile file () {
		return file;
	}
	
	public String getSequence () {
		return sequence;
	}
	
	public String getColorspace () {
		return colorspace;
	}
	
	public String getQualityString () {
		return quality;
	}
	
	public String getID () {
		return id;
	}
	
}
