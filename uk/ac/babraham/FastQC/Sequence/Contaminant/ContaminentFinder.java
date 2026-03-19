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

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Vector;

import uk.ac.babraham.FastQC.FastQCConfig;

public class ContaminentFinder {

	private static Contaminant [] contaminants;
	
	public static ContaminantHit findContaminantHit (String sequence) {
		if (contaminants == null) {
			contaminants = makeContaminantList();
		}
		
		ContaminantHit bestHit = null;
		
		for (int c=0;c<contaminants.length;c++) {
			ContaminantHit thisHit = contaminants[c].findMatch(sequence);
			
//			System.out.println("Best hit from "+c+" is "+thisHit);
			
			if (thisHit == null) continue; // No hit
			
			if (bestHit == null || thisHit.length() > bestHit.length() || (thisHit.length() == bestHit.length() && thisHit.percentID() > bestHit.percentID())) {
				bestHit = thisHit;
			}
			
		}
		
		return bestHit;
		
	}
	
	
	private static Contaminant [] makeContaminantList () {
		Vector<Contaminant>c = new Vector<Contaminant>();
		
		try {
			
			BufferedReader br = null;
			if (FastQCConfig.getInstance().contaminant_file == null) {
				InputStream rsrc=ContaminentFinder.class.getResourceAsStream("/Configuration/contaminant_list.txt");
				if (rsrc==null) throw new FileNotFoundException("cannot find Confituration/contaminant_list.txt");
				br =new BufferedReader(new InputStreamReader(rsrc));
			}
			else {
				br=new BufferedReader(new FileReader(FastQCConfig.getInstance().contaminant_file));
			}
			
			

			String line;
			while ((line = br.readLine())!= null){
				
				if (line.startsWith("#")) continue; // Skip comments
				if (line.trim().length() == 0) continue; // Skip blank lines
				
				String [] sections = line.split("\\t+");
				if (sections.length != 2) {
					System.err.println("Expected 2 sections for contaminant line but got "+sections.length+" from "+line);
					continue;
				}
				Contaminant con = new Contaminant(sections[0], sections[1]);
				c.add(con);	
			}
			
			br.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		
		return c.toArray(new Contaminant[0]);
	}
	/*
	public static void main (String [] args) {
		
		Config cfg=new Config();
		String query = "agagtgtagatctccgtggtcgccgtatca";
		
		ContaminantHit c = findContaminantHit(cfg,query);
		
		System.out.println("Query was "+query.length()+"bp Found hit "+c);
		
	}*/
	
}
