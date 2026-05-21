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
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.CountDownLatch;

import uk.ac.babraham.FastQC.Modules.BasicStats;
import uk.ac.babraham.FastQC.Modules.QCModule;
import uk.ac.babraham.FastQC.Sequence.Sequence;
import uk.ac.babraham.FastQC.Sequence.SequenceFile;
import uk.ac.babraham.FastQC.Sequence.SequenceFormatException;

public class AnalysisRunner implements Runnable {

	private SequenceFile file;
	private QCModule [] modules;
	private List<AnalysisListener> listeners = new ArrayList<AnalysisListener>();
	private volatile int percentComplete = 0;
	
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

		// Reader thread reads sequences in batches; NUM_PROCESSORS threads
		// drain those batches and each runs a disjoint subset of modules.
		final int BATCH_SIZE = 1024;
		final int QUEUE_CAPACITY = 32;
		final int NUM_PROCESSORS = 3;

		@SuppressWarnings("unchecked")
		ArrayBlockingQueue<Sequence[]>[] procQueues = new ArrayBlockingQueue[NUM_PROCESSORS];
		for (int p = 0; p < NUM_PROCESSORS; p++) {
			procQueues[p] = new ArrayBlockingQueue<>(QUEUE_CAPACITY);
		}
		final Sequence[] POISON = new Sequence[0];
		final SequenceFormatException[] readerError = { null };
		final Throwable[] processorError = { null };

		QCModule[][] moduleSplits = new QCModule[NUM_PROCESSORS][];
		int splitSize = (modules.length + NUM_PROCESSORS - 1) / NUM_PROCESSORS;
		for (int p = 0; p < NUM_PROCESSORS; p++) {
			int from = p * splitSize;
			int to = Math.min(from + splitSize, modules.length);
			moduleSplits[p] = Arrays.copyOfRange(modules, from, to);
		}

		// The same Sequence[] reference is published to every processor
		// queue. This is safe because Sequence is populated by file.next()
		// before put() and modules only read from it.
		final int[] seqCountHolder = { 0 };

		Thread reader = new Thread(() -> {
			int readerSeqCount = 0;
			try {
				Sequence[] batch = new Sequence[BATCH_SIZE];
				int idx = 0;
				while (file.hasNext()) {
					batch[idx++] = file.next();
					++readerSeqCount;
					if (idx == BATCH_SIZE) {
						for (int p = 0; p < NUM_PROCESSORS; p++) {
							procQueues[p].put(batch);
						}
						batch = new Sequence[BATCH_SIZE];
						idx = 0;

						// Same cadence as the single-threaded version: every
						// 1000 reads (rounded up to 1024 by BATCH_SIZE), fire
						// analysisUpdated only when we've crossed another 5%.
						if (file.getPercentComplete() >= percentComplete + 5) {
							percentComplete = (((int) file.getPercentComplete()) / 5) * 5;
							for (int li = 0; li < listeners.size(); li++) {
								listeners.get(li).analysisUpdated(file, readerSeqCount, percentComplete);
							}
						}
					}
				}
				if (idx > 0) {
					Sequence[] partial = Arrays.copyOf(batch, idx);
					for (int p = 0; p < NUM_PROCESSORS; p++) {
						procQueues[p].put(partial);
					}
				}
			}
			catch (SequenceFormatException e) {
				readerError[0] = e;
			}
			catch (InterruptedException e) {
				Thread.currentThread().interrupt();
			}
			finally {
				seqCountHolder[0] = readerSeqCount;
				// POISON must always be delivered, otherwise processor threads
				// block forever on take() and the run never completes.
				for (int p = 0; p < NUM_PROCESSORS; p++) {
					try { procQueues[p].put(POISON); } catch (InterruptedException e) {
						Thread.currentThread().interrupt();
					}
				}
			}
		}, "fastqc-reader");
		reader.setDaemon(true);
		reader.start();

		CountDownLatch processorsDone = new CountDownLatch(NUM_PROCESSORS);

		for (int p = 0; p < NUM_PROCESSORS; p++) {
			final QCModule[] myModules = moduleSplits[p];
			final ArrayBlockingQueue<Sequence[]> myQueue = procQueues[p];

			Thread processor = new Thread(() -> {
				try {
					while (true) {
						Sequence[] batch = myQueue.take();
						if (batch == POISON) break;

						for (int b = 0; b < batch.length; b++) {
							Sequence seq = batch[b];
							for (int m = 0; m < myModules.length; m++) {
								if (seq.isFiltered() && myModules[m].ignoreFilteredSequences()) continue;
								myModules[m].processSequence(seq);
							}
						}
					}
				}
				catch (InterruptedException e) {
					Thread.currentThread().interrupt();
				}
				catch (RuntimeException e) {
					// Don't let a misbehaving module deadlock the pipeline:
					// drain remaining batches so the reader doesn't block on
					// put(), and surface the failure after the join.
					processorError[0] = e;
					try {
						while (myQueue.take() != POISON) { /* drain */ }
					}
					catch (InterruptedException ie) {
						Thread.currentThread().interrupt();
					}
				}
				finally {
					processorsDone.countDown();
				}
			}, "fastqc-proc-" + p);
			processor.setDaemon(true);
			processor.start();
		}

		try {
			processorsDone.await();
			reader.join();
		}
		catch (InterruptedException e) {
			Thread.currentThread().interrupt();
		}

		int seqCount = seqCountHolder[0];

		if (readerError[0] != null) {
			i = listeners.iterator();
			while (i.hasNext()) {
				i.next().analysisExceptionReceived(file, readerError[0]);
			}
			return;
		}
		if (processorError[0] != null) {
			SequenceFormatException sfe = new SequenceFormatException(
					"Module processing failed: " + processorError[0].getMessage());
			sfe.initCause(processorError[0]);
			i = listeners.iterator();
			while (i.hasNext()) {
				i.next().analysisExceptionReceived(file, sfe);
			}
			return;
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
