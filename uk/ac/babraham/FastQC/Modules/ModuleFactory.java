/**
 * Copyright Copyright 2014-17 Simon Andrews
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

import java.util.ArrayList;
import java.util.List;

public class ModuleFactory {

	public static QCModule [] getStandardModuleList () {

		OverRepresentedSeqs os = new OverRepresentedSeqs();
		SequenceLengthDistribution ld = new SequenceLengthDistribution();
		BasicStats bs = new BasicStats();
		bs.addLengthDistribution(ld);
		
		// Modules disabled in the limits file are never shown in the report, so we
		// don't create them at all rather than paying to process every sequence.
		// Length distribution, overrepresented seqs and duplication are always kept
		// as other modules depend on their data.
		List<QCModule> modules = new ArrayList<QCModule>();
		modules.add(bs);
		if (enabled("quality_base")) modules.add(new PerBaseQualityScores());
		if (enabled("tile")) modules.add(new PerTileQualityScores());
		if (enabled("quality_sequence")) modules.add(new PerSequenceQualityScores());
		if (enabled("sequence")) modules.add(new PerBaseSequenceContent());
		if (enabled("gc_sequence")) modules.add(new PerSequenceGCContent());
		if (enabled("n_content")) modules.add(new NContent());
		modules.add(ld);
		modules.add(os.duplicationLevelModule());
		modules.add(os);
		if (enabled("adapter")) modules.add(new AdapterContent());
		if (enabled("kmer")) modules.add(new KmerContent());

		QCModule [] module_list = modules.toArray(new QCModule[0]);
	
		return (module_list);
	}

	private static boolean enabled (String module) {
		return ModuleConfig.getParam(module, "ignore") == 0;
	}
	
}
