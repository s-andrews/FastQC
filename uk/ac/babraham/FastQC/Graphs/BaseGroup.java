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
package uk.ac.babraham.FastQC.Graphs;

import java.util.Vector;

import uk.ac.babraham.FastQC.FastQCConfig;

/**
 * The base group class is a simple way to create a set of bins into
 * which positions within a read can be put such that early positions
 * get a group to themselves and later positions get averaged so that
 * general trends can still be observed.
 *
 */
public class BaseGroup {

	private int lowerCount;
	private int upperCount;

	public static BaseGroup [] makeBaseGroups (int maxLength) {
		

		// They might have set a fixed max length.  If the observed
		// length is longer than this then tough - they'll have to deal
		// with it, but if not then we'll use the global value instead 
		// of theirs
		
		if (FastQCConfig.getInstance().minLength > maxLength) {
			maxLength = FastQCConfig.getInstance().minLength;
		}
		
		if (FastQCConfig.getInstance().nogroup) {
			return(makeUngroupedGroups(maxLength));
		}
		

		if (FastQCConfig.getInstance().nogroup) {
			return(makeUngroupedGroups(maxLength));
		}
		else if (FastQCConfig.getInstance().expgroup) {
			return(makeExponentialBaseGroups(maxLength));
		}
		else {
			return(makeLinearBaseGroups(maxLength));
		}
	}

	public static BaseGroup [] makeUngroupedGroups (int maxLength) {

		int startingBase = 1;
		int interval = 1;

		Vector<BaseGroup> groups = new Vector<BaseGroup>();

		while (startingBase <= maxLength) {

			int endBase = startingBase+(interval-1);
			if (endBase > maxLength) endBase = maxLength;

			BaseGroup bg = new BaseGroup(startingBase, endBase);
			groups.add(bg);

			startingBase += interval;
		}

		return groups.toArray(new BaseGroup[0]);

	}

	public static BaseGroup [] makeExponentialBaseGroups (int maxLength) {

		int startingBase = 1;
		int interval = 1;

		Vector<BaseGroup> groups = new Vector<BaseGroup>();

		while (startingBase <= maxLength) {

			int endBase = startingBase+(interval-1);
			if (endBase > maxLength) endBase = maxLength;

			BaseGroup bg = new BaseGroup(startingBase, endBase);
			groups.add(bg);

			startingBase += interval;

			// See if we need to increase the interval
			if (startingBase == 10 && maxLength > 75) {
				interval = 5;
			}
			if (startingBase == 50 && maxLength > 200) {
				interval = 10;
			}
			if (startingBase == 100 && maxLength > 300) {
				interval = 50;
			}
			if (startingBase == 500 && maxLength > 1000) {
				interval = 100;
			}
			if (startingBase == 1000 && maxLength > 2000) {
				interval = 500;
			}


		}

		return groups.toArray(new BaseGroup[0]);

	}

	private static int getLinearInterval (int length) {
		// The the first 9bp as individual residues since odd stuff 
		// can happen there, then we find a grouping value which gives 
		// us a total set of groups below 75.  We limit the intervals
		// we try to sensible whole numbers.
		
		int [] baseValues = new int [] {2,5,10};
		int multiplier = 1;
		
		while (true) {
			for (int b=0;b<baseValues.length;b++) {
				int interval = baseValues[b] * multiplier;
				int groupCount = 9 + ((length-9)/interval);
				if ((length-9) % interval != 0) {
					groupCount += 1;
				}
				
				if (groupCount < 75) return interval;
			}
			
			multiplier *= 10;
			
			if (multiplier == 10000000) {
				throw new IllegalStateException("Couldn't find a sensible interval grouping for length '"+length+"'");
			}
		}
		
	}
	
	public static BaseGroup [] makeLinearBaseGroups (int maxLength) {

		// For lengths below 75bp we just return everything.
		if (maxLength <= 75) return makeUngroupedGroups(maxLength);
		
		// We need to work out what interval we're going to use.
		
		int interval = getLinearInterval(maxLength);
		
		
		int startingBase = 1;

		Vector<BaseGroup> groups = new Vector<BaseGroup>();

		while (startingBase <= maxLength) {

			int endBase = startingBase+(interval-1);

			if (startingBase < 10) endBase = startingBase;

			if (startingBase == 10 && interval > 10) {
				endBase = interval-1;
			}
			
			if (endBase > maxLength) endBase = maxLength;

			BaseGroup bg = new BaseGroup(startingBase, endBase);
			groups.add(bg);

			if (startingBase < 10) {
				startingBase +=1;
			}
			else if (startingBase == 10 && interval > 10) {
				startingBase = interval;
			}
			else {
				startingBase += interval;
			}

		}

		return groups.toArray(new BaseGroup[0]);

	}


	/**
	 * 
	 * @param lowerCount
	 * @param upperCount
	 */

	private BaseGroup (int lowerCount, int upperCount) {
		this.lowerCount = lowerCount;
		this.upperCount = upperCount;
	}

	public int lowerCount () {
		return lowerCount;
	}

	public int upperCount () {
		return upperCount;
	}

	public boolean containsValue (int value) {
		return value>=lowerCount && value<=upperCount;
	}

	public String toString () {
		if (lowerCount == upperCount) {
			return ""+lowerCount;
		}
		else {
			return ""+lowerCount+"-"+upperCount;
		}
	}

}
