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
package uk.ac.babraham.FastQC.Analysis;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import uk.ac.babraham.FastQC.Modules.BasicStats;
import uk.ac.babraham.FastQC.Modules.QCModule;
import uk.ac.babraham.FastQC.Sequence.Sequence;
import uk.ac.babraham.FastQC.Sequence.SequenceFile;
import uk.ac.babraham.FastQC.Sequence.SequenceFormatException;

public class AnalysisRunner implements Runnable {

	private SequenceFile file;
	private QCModule [] modules;
	private List<AnalysisListener> listeners = new ArrayList<AnalysisListener>();
	private int percentComplete = 0;
	
	public AnalysisRunner (SequenceFile file) {
		this.file = file;
	}
	
	public void addAnalysisListener (AnalysisListener l) {
		if (l != null && !listeners.contains(l)) {
			listeners.add(l);
		}
	}

	public void removeAnalysisListener (AnalysisListener l) {
		if (l != null && listeners.contains(l)) {
			listeners.remove(l);
		}
	}

	
	public void startAnalysis (QCModule [] modules) {
		this.modules = modules;
		for (int i=0;i<modules.length;i++) {
			modules[i].reset();
		}
		AnalysisQueue.getInstance().addToQueue(this);
	}

	public void run() {

		Iterator<AnalysisListener> i = listeners.iterator();
		while (i.hasNext()) {
			i.next().analysisStarted(file);
		}

		
		int seqCount = 0;
		while (file.hasNext()) {
			++seqCount;
			Sequence seq;
			try {
				seq = file.next();
			}
			catch (SequenceFormatException e) {
				i = listeners.iterator();
				while (i.hasNext()) {
					i.next().analysisExceptionReceived(file,e);
				}
				return;
			}
			
			for (int m=0;m<modules.length;m++) {
				if (seq.isFiltered() && modules[m].ignoreFilteredSequences()) continue;
				modules[m].processSequence(seq);
			}
			
			if (seqCount % 1000 == 0) {
			if (file.getPercentComplete() >= percentComplete+5) {
			
				percentComplete = (((int)file.getPercentComplete())/5)*5;
				
				i = listeners.iterator();
					while (i.hasNext()) {
						i.next().analysisUpdated(file,seqCount,percentComplete);
					}
					try {
						Thread.sleep(10);
					} 
					catch (InterruptedException e) {}
			}
			}
		}
		
		// We need to account for their potentially being no sequences
		// in the file.  In this case the BasicStats module never gets 
		// the file name so we need to explicitly pass it.
		
		if (seqCount == 0) {
			for (int m=0; m<modules.length; m++) {
				if (modules[m] instanceof BasicStats) {
					((BasicStats)modules[m]).setFileName(file.name());
				}
			}
		}
		
		i = listeners.iterator();
		while (i.hasNext()) {
			i.next().analysisComplete(file,modules);
		}

	}
	
}
