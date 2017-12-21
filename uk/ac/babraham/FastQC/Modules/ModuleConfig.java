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
package uk.ac.babraham.FastQC.Modules;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Sequence.Contaminant.ContaminentFinder;

public class ModuleConfig {

	private static HashMap<String, Double>parameters = readParams();
	
	
	private static HashMap<String, Double> readParams () {
		
		HashMap<String, Double>params = new HashMap<String, Double>();
		
		// Set the defaults to use if we don't have any overrides		
		params.put("duplication:warn",70d);
		params.put("duplication:error",	50d);
		params.put("kmer:warn",	2d);
		params.put("kmer:error",5d);
		params.put("n_content:warn",5d);
		params.put("n_content:error",20d);
		params.put("overrepresented:warn",	0.1);
		params.put("overrepresented:error",1d);
		params.put("quality_base_lower:warn",10d);
		params.put("quality_base_lower:error",	5d);
		params.put("quality_base_median:warn",	25d);
		params.put("quality_base_median:error",	20d);
		params.put("sequence:warn",	10d);
		params.put("sequence:error",20d);
		params.put("gc_sequence:warn",15d);
		params.put("gc_sequence:error",30d);
		params.put("quality_sequence:warn",	20d);
		params.put("quality_sequence:error",27d);
		params.put("tile:warn",	5d);
		params.put("tile:error",10d);
		params.put("sequence_length:warn",1d);
		params.put("sequence_length:error",	1d);
		params.put("adapter:warn", 5d);
		params.put("adapter:error", 10d);
		
		params.put("duplication:ignore",0d);
		params.put("kmer:ignore",0d);
		params.put("n_content:ignore",0d);
		params.put("overrepresented:ignore",0d);
		params.put("quality_base:ignore",0d);
		params.put("sequence:ignore",0d);
		params.put("gc_quality:ignore",0d);
		params.put("quality_sequence:ignore",0d);
		params.put("tile:ignore",0d);
		params.put("sequence_length:ignore",0d);
		params.put("adapter:ignore",0d);
				
		// Now read the config file to see if there are updated values for any of these.
		
		BufferedReader br;
		
		try {
			if (FastQCConfig.getInstance().limits_file == null) {
				InputStream rsrc=ContaminentFinder.class.getResourceAsStream("/Configuration/limits.txt");
				if (rsrc==null) throw new FileNotFoundException("cannot find Configuration/limits.txt");
				br =new BufferedReader(new InputStreamReader(rsrc));
			}
			else {
				br=new BufferedReader(new FileReader(FastQCConfig.getInstance().limits_file));
			}
			
			String line;
			while ((line = br.readLine()) != null) {
				
				if (line.startsWith("#")) continue;
				
				if (line.trim().length() == 0) continue;
				
				String [] sections = line.split("\\s+");
				if (sections.length != 3) {	
					System.err.println("Config line '"+line+"' didn't contain the 3 required sections");
				}
				
				if (!(sections[1].equals("warn") || sections[1].equals("error") || sections[1].equals("ignore"))) {
					System.err.println("Second config field must be error, warn or ignore, not '"+sections[1]+"'");
					continue;
				}
				
				double value;
				try {
					value = Double.parseDouble(sections[2]);
				}
				catch (NumberFormatException nfe) {
					System.err.println("Value "+sections[2]+" didn't look like a number");
					continue;
				}
				
				String key = sections[0]+":"+sections[1];
				params.put(key, value);
				
				
			}
		
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		

		
		
		return 	params;
		
	}
	
	public static Double getParam(String module, String level) {
		
		if (!(level.equals("warn") || level.equals("error") || level.equals("ignore"))) {
			throw new IllegalArgumentException("Level must be warn, error or ignore");
		}
		
		String key = module+":"+level;
		
		if (! parameters.containsKey(key)) {
			throw new IllegalArgumentException("No key called "+key+" in the config data");
		}
		
		return parameters.get(key);
		
	}
	
	
}
