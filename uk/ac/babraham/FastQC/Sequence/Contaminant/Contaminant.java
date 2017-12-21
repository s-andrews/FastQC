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
package uk.ac.babraham.FastQC.Sequence.Contaminant;

public class Contaminant {

	private String name;
	private char [] forward;
	private char [] reverse;
	
	public Contaminant (String name, String sequence) {
		this.name = name;
		
		sequence = sequence.toUpperCase();
		forward = sequence.toCharArray();
		reverse = new char[forward.length];
		for (int c=0;c<forward.length;c++) {
			int revPos = (reverse.length-1)-c;
			switch (forward[c]) {
				case 'G':
					reverse[revPos] = 'C';
					break;
				case 'A':
					reverse[revPos] = 'T';
					break;
				case 'T':
					reverse[revPos] = 'A';
					break;
				case 'C':
					reverse[revPos] = 'G';
					break;
				default:
					throw new IllegalArgumentException("Contaminant contained the illegal character '"+forward[c]+"'");
			}
		}
		
	}
	
	public ContaminantHit findMatch (String query) {
		query = query.toUpperCase();
		
		// We have a special case for queries between 8 - 20bp where we will allow a hit
		// if it's an exact substring of this contaminant
		if (query.length()<20 && query.length()>=8) {
			if ((new String(forward)).contains(query)) {
				return new ContaminantHit(this, ContaminantHit.FORWARD, query.length(), 100);
			}
			if ((new String(reverse)).contains(query)) {
				return new ContaminantHit(this, ContaminantHit.REVERSE, query.length(), 100);
			}
			
		}
		
		
		char [] q = query.toCharArray();
		
		ContaminantHit bestHit = null;
		
		// We're going to allow only one mismatch and will require 
		// a match of at least 20bp to consider this a match at all
		
		for (int offset=0-(forward.length-20);offset<q.length-20;offset++) {
			ContaminantHit thisHit = findMatch(forward,q,offset,ContaminantHit.FORWARD);
//			System.out.println("Best match from offset "+offset+" was "+thisHit);
			if (thisHit == null) continue;
			if (bestHit == null || thisHit.length()>bestHit.length()) {
				bestHit = thisHit;
			}
		}

		for (int offset=0-(forward.length-20);offset<q.length-20;offset++) {
			ContaminantHit thisHit = findMatch(reverse,q,offset,ContaminantHit.REVERSE);
			if (thisHit == null) continue;
			if (bestHit == null || thisHit.length()>bestHit.length()) {
				bestHit = thisHit;
			}
		}
		
		return bestHit;
		
	}
	
	private ContaminantHit findMatch (char [] ca, char [] cb, int offset, int direction) {
		
		ContaminantHit bestHit = null;
		
		int mismatchCount = 0;
		int start = 0;
		int end = 0;
				
		for (int i=0;i<ca.length;i++) {
			if (i+offset < 0) {
				start=i+1;
				continue;
			}
			if (i+offset >= cb.length) break;
			
			if (ca[i] == cb[i+offset]) {
				end = i;
			}
			else {
				++mismatchCount;
				if (mismatchCount>1) {
					// That's the end of this match, see if it's worth recording
					if (1+(end-start) > 20) {
						int id = (((1+(end-start))-(mismatchCount-1))*100)/(1+(end-start));
						if (bestHit == null || bestHit.length()< 1+(end-start) || (bestHit.length() == 1+(end-start) && bestHit.percentID()<id)) {
//							System.out.println("New best hit from "+start+"-"+end);
							bestHit = new ContaminantHit(this, direction, 1+(end-start), id);
						}
					}
					start = i+1;
					end = i+1;
					mismatchCount = 0;
				}
			}		
		}
		
		// See if we ended with a match.
		if (1+(end-start) > 20) {
			int id = (((1+(end-start))-mismatchCount)*100)/(1+(end-start));
			if (bestHit == null || bestHit.length()< 1+(end-start) || (bestHit.length() == 1+(end-start) && bestHit.percentID()<id)) {
//				System.out.println("New best hit from "+start+"-"+end);
				bestHit = new ContaminantHit(this, direction, 1+(end-start), id);
			}
		}
		
		return bestHit;		
		
	}
	
	public String name () {
		return name;
	}
	
	
}
