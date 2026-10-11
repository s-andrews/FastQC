# FastQC Java performance work

Goal: see whether the Java FastQC can match the performance of the Rust rewrite
([FastQC-Rust](https://github.com/ewels/FastQC-Rust), [performance report](https://ewels.github.io/FastQC-Rust/about/performance/report/))
while producing byte-identical results.

## Achieved so far

Through the `fastqc` launcher at its defaults (`report.sh`, Apple M5 Pro; optimised Java and
Rust median of 3 runs per file, Java 0.12.1 and master one run per file),
with `fastqc_data.txt` identical to FastQC-Rust v1.1 and Java 0.13.0 on every file:

| Data set | Java 0.12.1 | Java 0.13.0 (master) | Optimised Java | Rust v1.1 |
|---|---|---|---|---|
| WES, 4 files (29.3 GB) | 19m 47s | 7m 00s | **1m 18s** | 1m 32s |
| AVITI, 4 files (2.3 GB) | 1m 21s | 33.3 s | **10.6 s** | 9.0 s |

- **WES**: 15.2× faster than 0.12.1 and 5.4× faster than master; 15% faster than Rust, and faster
  on every file. Java uses ~2.1× Rust's CPU and ~10× its memory (~570 MB vs ~50 MB).
- **AVITI**: 7.6× faster than 0.12.1; 4–10% behind Rust on the larger files, 24–39% behind on the
  smaller ones (18% overall), mostly the ~0.7 s fixed cost of a JVM run and report generation.
- **lib-rapidgzip**: a pure-Java parallel gzip decompressor, as a separate module with its own
  build and tests, verified identical to `gzip -dc` on the AVITI and WES files.
- Every change is a branch of a stack (see below), opened as stacked PRs, each benchmarked step
  by step, with output checked identical to Rust's.

### Against the published report

The [FastQC-Rust report](https://ewels.github.io/FastQC-Rust/about/performance/report/) ran on an
M1 Pro (8+2 cores) with Java 21. Speed-ups over Java 0.12.1 carry over between the machines: the
report has 2.7× (WES) / 2.5× (AVITI) for Java 0.13.0 and 12.4× / 9.7× for Rust v1.1; here they
are 2.8× / 2.4× and 12.9× / 9.0×. Optimised Java is 15.2× / 7.6×.

Caveats: Java 0.12.1 and master one run per file; Java uses more threads than Rust (up to 7 per
file plus decompression, vs 6), so on a machine with fewer cores its relative position would
likely be worse. Not run: the report's 57.6 GB ONT file and its BAM (source not published).

## Benchmark setup

Everything lives in `benchmark/` (data, tools, build output and results are git-ignored).

| Script | What it does |
|---|---|
| `bench.sh [dataset...]` | Compiles the Java sources with `javac`, downloads FastQC-Rust v1.1.0, runs each tool under `/usr/bin/time`, appends wall/CPU/RSS to `results/results.tsv` and checks that `fastqc_data.txt` is identical across tools |
| `compare.sh [base] [tip] [dataset]` | In a git worktree, benchmarks `base` then each commit of the stack in order and prints a median summary; `MODE=alone` instead applies each commit alone on `base` |
| `report.sh file.fastq.gz...` | Mimics the published FastQC-Rust report: Java 0.12.1 (`fastqc` on the PATH), master and the optimised branch (`OPT_REF`) through their own launchers, `optbest` (`java` directly) and Rust, one run per file; results in `results/report/report.tsv`, output checked against Rust |

Options (env): `TOOLS=java,rust`, `REPS`, `THREADS`, `MEMORY`, `LABEL`, `SRC_DIR`, `PROFILE=1` (async-profiler 4.5 CPU flamegraph + collapsed stacks in `results/profile/`).

Datasets (same files as the Rust report, from ENA):

| Name | File | Size | Use |
|---|---|---|---|
| `small` | first 2M reads of ERR16944282_1 | 227 MB | quick iteration (default) |
| `aviti` | ERR16944282_1 (Element AVITI) | 510 MB | confirmation |
| `wes` | SRR7890918_1 (Illumina WES, HCC1395 tumor) | 7.9 GB | final validation |

Machine: Apple Silicon, 15 cores, Corretto 25 (code compiled for Java 21). Each tool runs with its default thread count
(Java run directly: 1 reader + up to 6 processors per file, plus decompression threads; Rust: up to 6). Prefer CPU time for comparing
changes, wall time is noisier.

## Baseline (master, 496d60e)

Profile on `aviti` (CPU share of ~24 s total):

| Where | CPU |
|---|---|
| `fastqc-reader` thread: zlib inflate 15%, `readLine` 3%, `toUpperCase` 1.4%, CRC32 1% | 26% |
| AdapterContent: 6 × `String.indexOf` per read | 17.5% |
| BasicStats, PerBaseSequenceContent, PerBaseQualityScores, GC, N content, ... | ~40% |
| KmerContent (disabled in `limits.txt`, still computed) | 7.6% |

The completion and queue checks also polled with `Thread.sleep(1000)` / `sleep(500)`, adding up to 1.5 s
of idle wall time per run (wall times were quantized to whole seconds).

## Stacked branches

```
master
 └─ perf/skip-disabled-modules
     └─ perf/single-pass-adapters
         └─ perf/no-sleep-polling
             └─ perf/java21
                 └─ perf/threaded-decompression
                     └─ perf/gzip-buffer
                         └─ perf/module-scheduling
                             └─ perf/lib-rapidgzip
                                 └─ perf/parallel-gzip
                                     └─ perf/more-processors
                                         └─ perf/branchless-counting
                                             └─ perf/byte-sequences
                                                 └─ perf/rapidgzip-speedups
                                                     └─ perf/rapidgzip-decoder
                                                         └─ perf/module-loops
                                                             └─ perf/memory
                                                                 └─ perf/font-registration
                                                                     └─ perf/launcher-threads
```

| Branch | Change |
|---|---|
| `perf/skip-disabled-modules` | Don't create modules with `ignore 1` in the limits file (Kmer by default); they were fed every read but never reported. Basic stats, length distribution, overrepresented and duplication are always kept as they depend on each other. |
| `perf/single-pass-adapters` | Find all adapters in one pass: 2-bit packed k-mers looked up in a small hash table, screened by a bitset over the last 8 bases. Falls back to `indexOf` for adapters of mixed length, > 31 bp or non-ACGT. Unit test checks it against `indexOf` on 100k random reads. |
| `perf/no-sleep-polling` | `OfflineRunner` awaits a `CountDownLatch`, `AnalysisQueue` waits on its monitor and is notified when a run is queued or a slot frees. |
| `perf/java21` | Compile with `release 21` (was source/target 1.8), CI on Java 21, docs updated. |
| `perf/threaded-decompression` | `ThreadedInputStream`: gzip/bzip2 decompression runs on its own thread, handing 256 KB chunks to the parser through a bounded queue (not used with `-t 1`). |
| `perf/gzip-buffer` | 64 KB input buffer for `GZIPInputStream` (default 512 bytes, one native read per 512 compressed bytes). |
| `perf/module-scheduling` | Modules were split statically into equal sized groups per processor thread (one thread had 5.2 s of work, another 1.8 s). Now each module has its own queue of batches, processed in order by one thread at a time, and any free processor thread picks up whichever module has work. A second commit fixes a race (found on WES) where a module could be scheduled with an empty queue, failing with a NullPointerException in ~1 of 2 WES runs; results were never affected. |
| `perf/lib-rapidgzip` | New module `lib-rapidgzip/` (own Gradle build and tests, run in CI): a pure-Java port of [rapidgzip](https://github.com/mxmlnkn/rapidgzip) parallel gzip decompression. See `lib-rapidgzip/README.md`. Output verified identical to `gzip -dc` on the aviti and WES files. |
| `perf/parallel-gzip` | FastQC reads regular `.gz` files with `RapidGzipInputStream` (4 threads, or half of `-t`); named pipes and `-t 1` keep `GZIPInputStream`. Ant build compiles the library sources. |
| `perf/more-processors` | Up to 6 processor threads per file (was 3): processors became the critical path. |
| `perf/branchless-counting` | Base counting in BasicStats, PerBaseSequenceContent and PerSequenceGCContent adds each comparison's 0/1 instead of branching on random bases (mostly mispredicted). These modules bounded the wall time on WES, each running on one thread at a time: 4M WES reads 3.7 s → 1.95 s. |
| `perf/byte-sequences` | `Sequence` holds bases (upper case) and qualities as bytes; FASTQ lines are split as bytes (8 bytes at a time, `\n`/`\r`/`\r\n` like `readLine`), only the id line is decoded, and per-base modules read bytes instead of `String.charAt`. Identical results for ASCII FASTQ (non-ASCII bytes, invalid FASTQ, were UTF-8 decoded before). |
| `perf/rapidgzip-speedups` | lib-rapidgzip: resolve markers through a lookup table (no mispredicted branch per value), `System.arraycopy` for non-overlapping copies, and check the window for markers every 64 KB instead of tracking them per copy. Decompression alone on 4 threads: 1.89 s → 1.53 s. |
| `perf/rapidgzip-decoder` | lib-rapidgzip: bit buffer kept in locals in the decoding loops (refill only below 48 bits, so several literals per refill), short back references copied with a loop, marker checks backed off exponentially, and chunk buffers reused through a pool. Decompression alone on 4 threads: 1.53 s → 1.31 s. |
| `perf/module-loops` | One increment per base: flat position × base / position × quality tables (filled into the arrays and `QualityCount` objects the reports use when calculating), NContent counts only Ns, BasicStats counts bytes 4-way interleaved. AdapterContent looks up the 8 bases at every (k-7)th position in a bitset of the adapters' 8-mers and checks candidate starts on a hit, instead of a lookup at every base. |
| `perf/memory` | Less allocation per read (skip the `+` line without copying, id kept as bytes and decoded on demand, overrepresented sequences built from the first 50 bytes), 16 read batches in flight instead of 32, `QualityCount` stores only the range of characters seen (per-tile qualities 26 MB → 6.8 MB), lib-rapidgzip reuses its I/O buffers and Huffman tables. |
| `perf/font-registration` | Don't register the bundled fonts with the `GraphicsEnvironment`, which enumerated every installed font (~150 ms per run). |
| `perf/launcher-threads` | The launchers no longer pass `-t 4` by default, so a file gets 7 threads capped at the CPUs; GC limited to 4 threads per file processed in parallel instead of 1; `-t` help and warnings updated from 4 to 7 threads per file. WES tumor_1 through the launcher: 42 s → 20 s (Rust 25 s). |

### Step by step (`compare.sh master perf/launcher-threads`, median of 3)

Each row is a branch tip, measured with `bench.sh` (`java` run directly, all threads). `small` is
the first 2M reads of ERR16944282_1 (AVITI), `wes-small` the first 4M reads of SRR7890918_1.

| Branch | small wall (s) | small CPU (s) | wes-small wall (s) | wes-small CPU (s) |
|---|---|---|---|---|
| master | 4.41 | 10.75 | 6.40 | 17.01 |
| skip-disabled-modules | 4.39 | 9.76 | 6.39 | 15.10 |
| single-pass-adapters | 4.42 | 8.68 | 6.39 | 13.04 |
| no-sleep-polling | 3.65 | 8.73 | 5.09 | 12.88 |
| java21 | 3.61 | 8.59 | 5.29 | 13.38 |
| threaded-decompression | 3.25 | 9.07 | 5.12 | 13.03 |
| gzip-buffer | 3.18 | 9.03 | 5.17 | 12.78 |
| module-scheduling | 2.90 | 8.68 | 3.77 | 12.75 |
| lib-rapidgzip (not used yet) | 2.90 | 8.82 | 3.77 | 12.78 |
| parallel-gzip | 2.76 | 11.53 | 3.88 | 15.25 |
| more-processors | 2.35 | 12.86 | 3.69 | 16.07 |
| branchless-counting | 2.32 | 12.80 | 2.00 | 12.05 |
| byte-sequences | 2.25 | 11.38 | 1.94 | 9.99 |
| rapidgzip-speedups | – | 10.92 | 1.72 | 9.26 |
| rapidgzip-decoder | 2.10 | 10.00 | 1.68 | 9.07 |
| module-loops | – | 9.66 | 1.67 | 8.24 |
| memory | 1.91 | 8.72 | 1.68 | 7.94 |
| font-registration | 1.72 | 8.16 | 1.51 | 7.76 |
| launcher-threads (launcher only) | 1.66 | 8.06 | 1.53 | 7.83 |

–: two of the three runs stalled 1.5–2.5 s idle (machine noise, CPU unaffected); a 5-run repeat
of those steps was clean and in line with the neighbouring rows. Wall times before
no-sleep-polling are quantized to whole seconds by the polling. `launcher-threads` only changes
the launchers, which `bench.sh` doesn't use: through the launcher it takes WES tumor_1 from
42.2 s to 19.4 s (Rust 24.8 s); a later A/B of the stack before and after the review fixes gave
medians of 20.3 s and 20.2 s on the same file.

## lib-rapidgzip

Decompression only (aviti, 510 MB → 2.94 GB, `-Xmx512m`):

| Decoder | Wall (s) | CPU (s) |
|---|---|---|
| `GZIPInputStream` | 3.58 | 3.55 |
| rapidgzip, 4 threads | 1.89 | 7.6 |
| rapidgzip, 8 threads | 1.42 | 8.3 |

In FASTQ about 55% of the output of a speculatively decoded chunk refers back to the unknown
window, so chunks stay in 16-bit marker mode from start to end. That makes the per-chunk memory
large (~20× the compressed chunk size), so chunks in flight are capped to a quarter of the heap,
and costs CPU: the Java decoder uses ~2.3× zlib's CPU.

## Optimised Java vs Rust v1.1.0 (default threads)

| Dataset | Version | Java wall (s) | Java CPU (s) | Rust wall (s) | Rust CPU (s) | Wall ratio |
|---|---|---|---|---|---|---|
| aviti | master | 7.4–9.4 | 22–24 | 2.2 | 6.1 | ~3.8× |
| aviti | `perf/no-sleep-polling` | 7.65 | 19.0 | 2.10 | 5.9 | 3.6× |
| aviti | `perf/module-scheduling` | 5.26 | 18.0 | 2.06 | 5.7 | 2.6× |
| aviti | `perf/more-processors` | 3.87 | 26.3 | 2.17 | 6.0 | 1.8× |
| aviti | `perf/rapidgzip-speedups` | 3.18 | 19.7 | 2.16 | 6.0 | 1.5× |
| aviti | `perf/module-loops` | 2.72–2.94 | 15.6 | 2.15 | 6.1 | 1.3× |
| wes-small (4M reads) | `perf/rapidgzip-speedups` | 1.90 | 9.4 | 1.00 | 3.1 | 1.9× |
| wes-small (4M reads) | `perf/module-loops` | 1.68 | 8.2 | 1.04 | 3.2 | 1.6× |
| wes | `perf/no-sleep-polling` | 115.8 | 316 | 27.0 | 81 | 4.3× |
| wes | `perf/module-scheduling` | 79.3 | 313 | 25.5 | 81 | 3.1× |
| wes | `perf/more-processors` | 81.4 | 408 | 25.1 | 79 | 3.2× |
| wes | `perf/branchless-counting` | 30.4 | 270 | 26.0 | 81 | 1.17× |
| wes | `perf/byte-sequences` | 31.6 | 230 | 25.5 | 78 | 1.24× |
| wes | `perf/rapidgzip-speedups` | 23.9 | 205 | 25.8 | 81 | 0.93× |
| wes | `perf/module-loops` | **23.0–24.5** | 180–190 | 25.6–26.1 | 79–81 | **0.92×** |

Output identical to Rust in every run. RSS: Java ~500–650 MB, Rust 50–80 MB.

## Fixed cost

A FastQC run on a 5-read file took 0.81 s in Java vs 0.05 s for Rust: JVM startup, class loading
and JIT, and report generation (PNG/SVG charts, fonts, HTML, zip: ~370 ms of CPU). It's negligible
on WES but a large part of the remaining gap on the smaller files.

**Font registration** (`perf/font-registration`): `GraphicsEnvironment.registerFont` made the JDK
enumerate every installed font, ~150 ms of every run. FastQC only uses the `Font` objects it loads,
so registering them isn't needed: 5 reads 0.82 s → 0.70 s, reports byte-identical.

### AOT cache (JDK 25) — evaluated, not adopted

Tried the JDK 25 AOT cache (JEP 483 class loading and linking, JEP 515 method profiles, created in
one step with `-XX:AOTCacheOutput` from a training run on 200k reads) and a CDS dynamic archive,
on Corretto 25 (median of 3–7 runs):

| File | Classes dir | jar | jar + CDS | jar + AOT cache |
|---|---|---|---|---|
| 5 reads (before the font fix) | 0.82 s | 0.81 s | 0.77 s | 0.73 s |
| 5 reads (after the font fix) | 0.67 s | – | – | 0.60 s |
| `small` (after the font fix) | 1.66 s | – | – | 1.56 s |
| AVITI (after the font fix) | 2.67 s | – | – | 2.66 s |

- The cache only applies to classes loaded from jars: with FastQC's class directory on the
  classpath the JVM refuses to create it (`non-empty directory`), so FastQC would have to ship as
  a `fastqc.jar`.
- A cache is only valid for the exact JDK build and classpath that created it, so it would have to
  be created on each user's machine (first run with `-XX:AOTCacheOutput`, per user and JDK
  version, with concurrent first runs to handle).

**Decision:** not adopted for now. ~0.07–0.1 s on small files and nothing on large ones doesn't
justify changing the distribution to a jar (build, release workflow, three launchers) plus cache
management in the launchers. Worth revisiting if FastQC moves to a jar anyway, or with JDKs that
add AOT compiled code (Leyden), when it becomes mostly launcher work (opt-in on Java 25+, e.g. a
`--aot-cache` option).

### GraalVM native image — evaluated, not adopted

Built FastQC (`perf/font-registration` code, as a jar) with Oracle GraalVM 25 `native-image`
(`-march=native`, reachability metadata from the tracing agent). AWT isn't supported by native
images on macOS, so this ran in Linux arm64 containers on this Mac (Docker, 15 CPUs, 8 GB), against
the same jar on Temurin 25 in the same VM. `-Xmx512m`, median of 3; aviti runs alternated
JVM/native as the Docker VM is noisy. The XSL-FO output was skipped in both (see below).

| File | JVM (Temurin 25) | Native | Native + PGO (trained on `small`) |
|---|---|---|---|
| tiny (5 reads) | 0.33 s, 0.51 s CPU, 176 MB | 0.29 s, 0.10 s CPU, 67 MB | 0.31 s, 0.10 s CPU, 54 MB |
| `small` | 1.42 s, 8.2 s CPU, 463 MB | 1.64 s, 9.1 s, 283 MB | 1.44 s, 7.3 s, 262 MB |
| `wes-small` | 1.35 s, 8.6 s CPU, 399 MB | 1.37 s, 9.0 s, 251 MB | 1.21 s, 7.5 s, 247 MB |
| aviti | 2.65–3.34 s, 16.4–19.7 s CPU, 450–500 MB | 3.0 s, 20.1 s, 368 MB | 2.97–4.01 s, 16.4–21.0 s, 250–340 MB |

- **Memory**: 30–70% less in all cases.
- **Speed**: no overall gain. Startup CPU is ~5× lower, but on Linux the JVM already starts in
  ~0.33 s. With PGO native is even on small files, ~10% faster on 4M WES reads and 10–20% slower
  on AVITI (the JIT's run-time profiles win on the long hot loops); without PGO it's 10–25% slower
  everywhere. G1 (`--gc=G1`) instead of the serial GC didn't help.
- Report data identical; chart PNGs differ byte-wise but look the same (font rasterisation).
- On Linux the JVM's fixed cost is much lower than on macOS (5 reads: 0.33 s vs 0.70 s).

Obstacles to shipping it:
- macOS: no AWT in native images, so no charts; Mac users would stay on the JVM.
- Linux: needs `freetype`, `fontconfig` and a system font, and `-Djava.home` pointing to a JDK so
  the font code finds its configuration (else "Fontconfig head is null").
- The XSL-FO output compiles `fastqc2fo.xsl` into a class at run time (XSLTC), which native images
  can't do; the agent's experimental class-define support crashed the build. Fix would be to
  precompile the translet at build time.
- Reflection/resource/JNI metadata to maintain, per-OS/arch builds (58 MB executable plus AWT
  libraries), a PGO training run in the release pipeline, 1–2 min builds.

Also found: when report generation throws an `Error` (not an `Exception`), the analysis thread dies
without counting down the files remaining and FastQC waits forever. It happens on master too; worth
fixing regardless.

**Decision:** not adopted, as CPU and speed matter more than memory here and native gives no speed
gain. Worth revisiting only if memory becomes the priority (many small concurrent jobs).

## Plan

Remaining gap (`perf/memory`, fixed cost from `perf/font-registration`, vs Rust v1.1):

| | Java | Rust | Gap |
|---|---|---|---|
| Wall, WES (7.9 GB) | 23.0–24.5 s | 25.6–26.1 s | none (0.92×) |
| Wall, AVITI (510 MB) | 2.75–2.85 s | 2.15 s | 1.3× |
| Wall, 4M WES reads | 1.69 s | 1.04 s | 1.6× |
| CPU, WES | 180–190 s | 80 s | ~2.3× |
| CPU, AVITI | 15.4–16.1 s | 6.1 s | ~2.6× |
| Fixed cost per run (5 reads) | 0.70 s | 0.05 s | ~0.65 s |
| Peak RSS, default launcher (`-Xmx512m`) | ~640 MB | ~50–75 MB | ~10× |
| Peak RSS, `-Xms32m -Xmx320m` | ~380–400 MB (+3% time) | ~50–75 MB | ~6× |
| Live heap at peak (AVITI) | ~240 MB | | |
| Busy threads (default) | ~11 | ≤ 6 | |

- **Wall time** is on par on large files, but behind on small and medium ones, mostly the fixed
  cost of a JVM run and report generation.
- **CPU**: the gzip decoder (16-bit marker mode for all of a FASTQ chunk, ~1.4× zlib's CPU on one
  thread) and the modules. Java matches Rust's wall time by using more cores, so on a machine with
  fewer cores, or at an equal thread budget, it would be slower.
- **Memory**: at the default the JVM commits the whole 512 MB heap up front (initial heap = 1/64
  of RAM, capped at `-Xmx`), whatever is live. Live data is mostly lib-rapidgzip chunks in flight
  (16 bit values, ~20 bytes per compressed byte each, capped to a quarter of the heap) and reads in
  flight.

1. **Fixed cost** — font registration done (`perf/font-registration`); AOT cache / CDS evaluated
   and not adopted (see *Fixed cost*). Next: render the module charts in parallel (~370 ms of
   report CPU, sequential now). Matters for small and medium files.
2. **Memory** — lower priority than speed and CPU: only changes that cost no time. Done:
   per-read allocations, quality counts and decompressor buffers (`perf/memory`). Not adopted as
   they cost time: `-Xms64m` in the launchers, 512 KB gzip chunks, decoding chunks twice, native
   image (see *Memory* and *Fixed cost*).
3. **Thread budget** — launcher defaults done (`perf/launcher-threads`: no `-t 4`, 4 GC threads
   per file; GC threads 1/2/4/default are within noise on WES, 4 is ~4% faster than 1 on AVITI).
   Next: make `-t` cap all threads (decompression included), and compare with Rust at an equal
   thread count.

## Memory

Peak RSS is ~640 MB at the launcher's default `-Xmx512m`, but that's mostly heap the JVM commits
up front: on machines with more than 32 GB of RAM the default initial heap (1/64 of RAM) is the
whole `-Xmx`. At the end of an AVITI run the JVM has committed 665 MB: 537 MB of heap, 70 MB of G1
structures, ~30 MB of code and metadata.

Live data peaks at ~240 MB (`perf/memory`, AVITI): lib-rapidgzip chunks in flight (~83 MB of 16 bit
values plus resolved bytes; chunks in flight are capped to a quarter of the heap), reads in flight,
and the modules' counts.

AVITI with `-Xms32m` (grow the heap as needed), by `-Xmx`:

| `-Xmx` | Wall (s) | Peak RSS (MB) |
|---|---|---|
| 512 (default, no `-Xms`) | 2.74–2.78 | ~640 |
| 512 | 2.80–2.86 | 430–460 |
| 320 | 2.90–2.94 | 380–400 |
| 256 | 3.24–3.28 | 330–350 |
| 192 | 3.32–3.42 | ~300 |

Below ~320 MB the decompressor gets fewer chunks in flight and slows down. 512 KB chunks save a
further ~60 MB for 3–5% more time. Tried and dropped: resolving markers in place in a byte array
(holds the 16 bit array until the chunk is read, so more memory, and slower copies).

### JVM heap tuning

AVITI, Corretto 25 on macOS, medians of 3–5 alternated runs (the machine was busy, load ~7–10,
so differences of a few % are noise):

| Option | Wall (s) | CPU (s) | Peak RSS (MB) |
|---|---|---|---|
| default (`-Xmx512m`) | 2.63–3.11 | 15.5–17.8 | ~600 |
| `-Xms32m` / `-Xms64m` / `-Xms128m` | within noise of default | within noise | 390–405 |
| `-XX:+UseCompactObjectHeaders` | +10–15% | +10–15% | no change (~400 with `-Xms`) |
| `-XX:+UseParallelGC`, `-XX:+UseSerialGC` (with `-Xms32m`) | +15–25% | +15–20% | 565–595 |
| ZGC, `-XX:SoftMaxHeapSize=256m` | +55% | +75% | 466 |
| `-XX:ReservedCodeCacheSize=48m -XX:MaxMetaspaceSize=64m` | no gain | | −20 |

Small and 4M WES reads, default vs `-Xms64m`: +1–2% wall (noise), 600 → 400 MB and 554 → 389 MB.

Only the initial heap size helps: G1 then grows the heap to what's needed (~350–400 MB here)
instead of committing the whole `-Xmx` at start.

**The launchers' own settings**: `fastqc` defaults to `-t 4` and adds `-XX:ParallelGCThreads=1`
(to limit GC CPU on big shared machines), while the benchmarks above run `java` directly, with all
threads. AVITI through the launcher's settings:

| Settings | Wall (s) | CPU (s) | Peak RSS (MB) |
|---|---|---|---|
| `-t 4 -XX:ParallelGCThreads=1` (launcher default) | 3.50–3.55 | 12.9–13.2 | 560–600 |
| same + `-Xms64m` | 3.58–3.61 | 13.3–13.4 | 343–405 |
| no `-t`, default GC threads (`bench.sh`) | 2.50–2.61 | 14.9–15.6 | ~600 |

So users of the `fastqc` command get ~3.5 s on AVITI, not the ~2.6 s measured by `bench.sh`, and the
comparison with Rust (which uses up to 6 threads by default) should be redone with the launcher, as
part of the thread budget item of the plan.

**Decision:** `-Xms64m` not adopted. With the launcher's single GC thread it saves ~220 MB but costs
~1.5% wall time and ~2.6% CPU, and speed and CPU take priority over memory. (It was tried as a
launcher change, `-Xms64m` in `fastqc` and `fastqc.pl`, and dropped.)
