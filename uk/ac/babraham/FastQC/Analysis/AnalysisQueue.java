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

import java.util.concurrent.LinkedBlockingDeque;
import java.util.concurrent.atomic.AtomicInteger;

import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Modules.QCModule;
import uk.ac.babraham.FastQC.Sequence.SequenceFile;

public class AnalysisQueue implements Runnable, AnalysisListener{

	private static AnalysisQueue instance = new AnalysisQueue();
	
	private LinkedBlockingDeque<AnalysisRunner>queue = new LinkedBlockingDeque<AnalysisRunner>();
	
	private AtomicInteger availableSlots = new AtomicInteger(1);
	private AtomicInteger usedSlots = new AtomicInteger(0);
	private volatile int processorsPerFile = 0;

	// Per-file pipeline = 1 reader + up to MAX_PROCESSORS_PER_FILE processor
	// threads. Past 3 processors gains plateau as the reader becomes the bottleneck.
	private static final int MAX_PROCESSORS_PER_FILE = 3;
	private static final int THREADS_PER_FILE = 1 + MAX_PROCESSORS_PER_FILE;

	public static AnalysisQueue getInstance () {
		return instance;
	}

	private AnalysisQueue () {
		// Sized for a single file; OfflineRunner overrides via configure().
		configure(1);

		Thread t = new Thread(this);
		t.start();
	}

	/**
	 * Size the CPU budget. -t/-Dfastqc.threads is a total-thread budget
	 * split between outer concurrency (files in parallel) and inner
	 * concurrency (per-file pipeline). Unset defaults to
	 * <code>min(THREADS_PER_FILE * max(1, expectedFiles), availableProcessors)</code>.
	 * A budget of 1 makes AnalysisRunner take its single-threaded path.
	 * Call before submitting any AnalysisRunner.
	 */
	void configure (int expectedFiles) {
		int totalThreads;
		if (FastQCConfig.getInstance().threads != null) {
			totalThreads = FastQCConfig.getInstance().threads;
		}
		else {
			int filesForBudget = Math.max(1, expectedFiles);
			long requested = (long) THREADS_PER_FILE * filesForBudget; // long to avoid overflow
			totalThreads = (int) Math.min(requested, Runtime.getRuntime().availableProcessors());
		}
		processorsPerFile = totalThreads <= 1 ? 0 : Math.min(MAX_PROCESSORS_PER_FILE, totalThreads - 1);
		int actualThreadsPerFile = processorsPerFile == 0 ? 1 : 1 + processorsPerFile;
		availableSlots.set(Math.max(1, totalThreads / actualThreadsPerFile));
	}

	int getProcessorsPerFile() {
		return processorsPerFile;
	}
	
	public void addToQueue (AnalysisRunner runner) {
		queue.add(runner);
	}

	public void run() {

		while (true) {
//			System.err.println("Status available="+availableSlots+" used="+usedSlots+" queue="+queue.size());
			if (availableSlots.intValue() > usedSlots.intValue() && queue.size() > 0) {
				usedSlots.incrementAndGet();
				AnalysisRunner currentRun = queue.getFirst();
				queue.removeFirst();
				currentRun.addAnalysisListener(this);
				Thread t = new Thread(currentRun);
				t.start();
			}
			
			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {}
		}
	}

	public void analysisComplete(SequenceFile file, QCModule[] results) {
		usedSlots.decrementAndGet();
	}

	public void analysisUpdated(SequenceFile file, int sequencesProcessed, int percentComplete) {}

	public void analysisExceptionReceived(SequenceFile file, Exception e) {
		usedSlots.decrementAndGet();
	}

	public void analysisStarted(SequenceFile file) {}
	
	
}
