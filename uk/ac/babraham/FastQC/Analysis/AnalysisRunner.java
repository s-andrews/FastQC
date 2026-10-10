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
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.Semaphore;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import uk.ac.babraham.FastQC.FastQCConfig;
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

	// We might have read length filters. If we do then we want to 
	// check for these quickly so we'll make some booleans to speed
	// this up.
	private boolean filterLengths = false;
	private boolean minLengthFilter = false;
	private int minLength;
	private boolean maxLengthFilter = false;
	private int maxLength;
	
	public AnalysisRunner (SequenceFile file) {
		this.file = file;
		
		if (FastQCConfig.getInstance().minLength != 0) {
			filterLengths = true;
			minLengthFilter = true;
			minLength = FastQCConfig.getInstance().minLength;
		}
		
		if (FastQCConfig.getInstance().maxLength != 0) {
			filterLengths = true;
			maxLengthFilter = true;
			maxLength = FastQCConfig.getInstance().maxLength;
		}
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
		try {
			int numProcessors = AnalysisQueue.getInstance().getProcessorsPerFile();
			if (numProcessors == 0) {
				runSequential();
				return;
			}
			runParallel(numProcessors);
		}
		catch (Throwable t) {
			// Anything escaping the run (a module failing on the -t 1 path, an Error
			// such as OutOfMemoryError) would end this thread without telling the
			// listeners: the slot would never be freed and OfflineRunner would wait
			// for this file forever. Report it as a failure of the file instead.
			SequenceFormatException sfe = new SequenceFormatException("Analysis failed: " + t);
			sfe.initCause(t);
			Iterator<AnalysisListener> i = listeners.iterator();
			while (i.hasNext()) {
				i.next().analysisExceptionReceived(file, sfe);
			}
		}
	}

	private void runSequential() {

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
				
				if (filterLengths) {
					if (minLengthFilter) {
						if (seq.getSequence().length() < minLength) {
							continue;
						}
					}
					if (maxLengthFilter) {
						if (seq.getSequence().length() > maxLength) {
							continue;
						}
					}
				}
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

	private void runParallel(final int numProcessors) {

		Iterator<AnalysisListener> i = listeners.iterator();
		while (i.hasNext()) {
			i.next().analysisStarted(file);
		}

		// Reader fills batches of sequences, which are queued for each module.
		// A module processes its batches in order, one thread at a time, but
		// any of the numProcessors threads picks up whichever module has work
		// pending, so the load balances itself whatever each module costs.
		final int BATCH_SIZE = 1024;
		final int BATCHES_IN_FLIGHT = 32;

		final Semaphore inFlight = new Semaphore(BATCHES_IN_FLIGHT);
		final LinkedBlockingQueue<ModuleWork> ready = new LinkedBlockingQueue<ModuleWork>();
		final ModuleWork POISON = new ModuleWork(null);
		final ModuleWork[] work = new ModuleWork[modules.length];
		for (int m = 0; m < modules.length; m++) {
			work[m] = new ModuleWork(modules[m]);
		}
		// One count per batch not yet processed by every module, plus one held
		// by the reader until it's done; whoever takes it to zero stops the processors.
		final AtomicInteger pending = new AtomicInteger(1);
		final Runnable finishIfDone = () -> {
			if (pending.decrementAndGet() == 0) {
				for (int p = 0; p < numProcessors; p++) ready.add(POISON);
			}
		};
		final SequenceFormatException[] readerError = { null };
		final AtomicReference<Throwable> processorError = new AtomicReference<Throwable>();
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
						dispatch(new Batch(batch, work.length), work, ready, inFlight, pending);
						batch = new Sequence[BATCH_SIZE];
						idx = 0;

						// Fire analysisUpdated only when the file advances another 5%.
						if (file.getPercentComplete() >= percentComplete + 5) {
							percentComplete = (((int) file.getPercentComplete()) / 5) * 5;
							for (int li = 0; li < listeners.size(); li++) {
								listeners.get(li).analysisUpdated(file, readerSeqCount, percentComplete);
							}
						}
					}
				}
				if (idx > 0) {
					dispatch(new Batch(Arrays.copyOf(batch, idx), work.length), work, ready, inFlight, pending);
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
				// Must always release the reader's count or processors block on take() forever.
				finishIfDone.run();
			}
		}, "fastqc-reader");
		reader.setDaemon(true);
		reader.start();

		CountDownLatch processorsDone = new CountDownLatch(numProcessors);

		for (int p = 0; p < numProcessors; p++) {
			Thread processor = new Thread(() -> {
				try {
					while (true) {
						ModuleWork w = ready.take();
						if (w == POISON) break;

						// Only the thread holding the module polls its queue. It can be
						// empty: a thread finishing the module may find a batch queued, and
						// schedule it again after another thread has processed that batch.
						Batch batch = w.batches.poll();
						if (batch == null) {
							w.scheduled.set(false);
							w.schedule(ready);
							continue;
						}
						if (processorError.get() == null) {
							try {
								processBatch(w.module, batch.sequences);
							}
							catch (RuntimeException e) {
								// Keep consuming batches so the reader doesn't block, and
								// surface the failure to listeners after join().
								processorError.compareAndSet(null, e);
							}
						}
						w.scheduled.set(false);
						w.schedule(ready);

						if (batch.modulesRemaining.decrementAndGet() == 0) {
							inFlight.release();
							finishIfDone.run();
						}
					}
				}
				catch (InterruptedException e) {
					Thread.currentThread().interrupt();
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
		if (processorError.get() != null) {
			SequenceFormatException sfe = new SequenceFormatException(
					"Module processing failed: " + processorError.get().getMessage());
			sfe.initCause(processorError.get());
			i = listeners.iterator();
			while (i.hasNext()) {
				i.next().analysisExceptionReceived(file, sfe);
			}
			return;
		}

		// We need to account for their potentially being no sequences
		// in the file.  In this case the BasicStats module never gets
		// the file name so we need to explicitly pass it.  The sequence
		// counts aren't reliable here as they may have been reduced due
		// to sequence length filtering, ie sequences read != sequences used.
		// There isn't really a downside to setting the file name again
		// so let's do that.

		for (int m=0; m<modules.length; m++) {
			if (modules[m] instanceof BasicStats) {
				((BasicStats)modules[m]).setFileName(file.name());
			}
		}

		i = listeners.iterator();
		while (i.hasNext()) {
			i.next().analysisComplete(file,modules);
		}

	}
	
	private void processBatch (QCModule module, Sequence[] batch) {
		for (int b = 0; b < batch.length; b++) {
			Sequence seq = batch[b];

			// We need to apply length filters here if we have them.
			if (filterLengths) {
				if (minLengthFilter) {
					if (seq.getSequence().length() < minLength) {
						continue;
					}
				}
				if (maxLengthFilter) {
					if (seq.getSequence().length() > maxLength) {
						continue;
					}
				}
			}

			if (seq.isFiltered() && module.ignoreFilteredSequences()) continue;
			module.processSequence(seq);
		}
	}

	private static void dispatch (Batch batch, ModuleWork[] work, BlockingQueue<ModuleWork> ready, Semaphore inFlight, AtomicInteger pending) throws InterruptedException {
		inFlight.acquire();
		pending.incrementAndGet();
		for (int m = 0; m < work.length; m++) {
			work[m].batches.add(batch);
			work[m].schedule(ready);
		}
	}

	private static class Batch {
		final Sequence[] sequences;
		final AtomicInteger modulesRemaining;

		Batch (Sequence[] sequences, int modules) {
			this.sequences = sequences;
			this.modulesRemaining = new AtomicInteger(modules);
		}
	}

	/**
	 * The batches waiting for a module. The module is on the ready queue, or
	 * being processed, while scheduled is set, so only one thread runs it at a time.
	 */
	private static class ModuleWork {
		final QCModule module;
		final ConcurrentLinkedQueue<Batch> batches = new ConcurrentLinkedQueue<Batch>();
		final AtomicBoolean scheduled = new AtomicBoolean(false);

		ModuleWork (QCModule module) {
			this.module = module;
		}

		void schedule (BlockingQueue<ModuleWork> ready) {
			if (!batches.isEmpty() && scheduled.compareAndSet(false, true)) {
				ready.add(this);
			}
		}
	}

}
