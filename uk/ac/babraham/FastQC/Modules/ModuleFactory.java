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

import uk.ac.babraham.FastQC.FastQCConfig;

public class ModuleFactory {

	private FastQCConfig config;
	
	public ModuleFactory(FastQCConfig config) {
		this.config = config;
	}

	public QCModule [] getStandardModuleList () {

		OverRepresentedSeqs os = new OverRepresentedSeqs(config);
		
		QCModule [] module_list = new QCModule [] {
				new BasicStats(config),
				new PerBaseQualityScores(config),
				new PerTileQualityScores(config),
				new PerSequenceQualityScores(config),
				new PerBaseSequenceContent(config),
				new PerSequenceGCContent(config),
				new NContent(config),
				new SequenceLengthDistribution(config),
				os.duplicationLevelModule(),
				os,
				new AdapterContent(config),
				new KmerContent(config),
			};
	
		return (module_list);
	}
	
}
