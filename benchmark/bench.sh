#!/usr/bin/env bash
#
# Self-contained FastQC benchmark: Java (this source tree) vs FastQC-Rust.
# Datasets mirror https://ewels.github.io/FastQC-Rust/about/performance/report/
#
# Usage: benchmark/bench.sh [dataset...]
#   datasets: small (first 2M reads of aviti, for quick iteration, default)
#             wes-small (first 4M reads of wes)
#             aviti (ERR16944282_1, 510 MB) | wes (SRR7890918_1, 7.9 GB, final validation)
#
# Env:
#   TOOLS=java,rust   tools to run
#   REPS=1            repetitions per tool/dataset
#   THREADS=          passed as -t to both tools (unset = each tool's default)
#   MEMORY=512        Java heap in MB (same as the fastqc launcher default)
#   JAVA_OPTS=        extra JVM options for the Java runs
#   TIMEOUT=600       seconds before a run is killed and counted as failed
#   SRC_DIR=path      Java source tree to compile (default: this repo), e.g. a git worktree
#   LABEL=name        tag recorded with the commit in results.tsv
#   PROFILE=1         profile the Java run with async-profiler (CPU, 1ms) and write
#                     results/profile/<dataset>.{html,collapsed} (timings are then skewed)
#
# Everything (data, Rust binary, compiled classes, results) stays under benchmark/
# and is ignored by git via benchmark/.gitignore.
set -euo pipefail

BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$BENCH_DIR")"
SRC_DIR="${SRC_DIR:-$REPO_DIR}"
DATA_DIR="$BENCH_DIR/data"
TOOLS_DIR="$BENCH_DIR/tools"
CLASSES_DIR="$BENCH_DIR/build/classes"
RESULTS_DIR="$BENCH_DIR/results"

TOOLS="${TOOLS:-java,rust}"
REPS="${REPS:-1}"
THREADS="${THREADS:-}"
MEMORY="${MEMORY:-512}"
PROFILE="${PROFILE:-}"
RUST_VERSION="v1.1.0"
ASPROF_VERSION="4.5"

# name -> "url md5"
dataset_spec() {
  case "$1" in
    aviti) echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR169/082/ERR16944282/ERR16944282_1.fastq.gz 67270051a5f8e11931066e47c50f6cf7" ;;
    wes)   echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789/008/SRR7890918/SRR7890918_1.fastq.gz 8d8f82e4f35e1b4babc4bddb5150d50b" ;;
    *) echo "Unknown dataset: $1 (use: aviti, wes)" >&2; exit 1 ;;
  esac
}

log() { echo "[bench] $*" >&2; }

fetch_dataset() {
  local spec url md5 file
  # Subsets of a dataset for quick iteration: name -> "source dataset, reads, file"
  local subset=""
  case "$1" in
    small) subset="aviti 2000000 ERR16944282_1.2M.fastq.gz" ;;
    wes-small) subset="wes 4000000 SRR7890918_1.4M.fastq.gz" ;;
  esac
  if [[ -n "$subset" ]]; then
    read -r from reads name <<< "$subset"
    file="$DATA_DIR/$name"
    if [[ ! -f "$file" ]]; then
      local src; src=$(fetch_dataset "$from")
      log "Creating $file (first $reads reads of $from)"
      # head closes the pipe early; ignore gunzip's SIGPIPE
      { gunzip -c "$src" || true; } | head -n $((reads * 4)) | gzip -n -6 > "$file.tmp" && mv "$file.tmp" "$file"
    fi
    echo "$file"; return
  fi
  spec=$(dataset_spec "$1"); url=${spec% *}; md5=${spec#* }
  file="$DATA_DIR/$(basename "$url")"
  mkdir -p "$DATA_DIR"
  if [[ ! -f "$file.ok" ]]; then
    log "Downloading $url"
    curl -fL -C - -o "$file" "$url"
    log "Verifying md5"
    [[ "$(md5 -q "$file" 2>/dev/null || md5sum "$file" | cut -d' ' -f1)" == "$md5" ]] \
      || { echo "md5 mismatch for $file" >&2; exit 1; }
    touch "$file.ok"
  fi
  echo "$file"
}

fetch_rust() {
  local os arch asset bin
  case "$(uname -s)" in Darwin) os=macos ;; Linux) os=linux ;; *) echo "unsupported OS" >&2; exit 1 ;; esac
  case "$(uname -m)" in arm64|aarch64) arch=aarch64 ;; x86_64) arch=x86_64 ;; esac
  asset="fastqc-$os-$arch"
  bin="$TOOLS_DIR/rust-$RUST_VERSION/$asset/fastqc"
  if [[ ! -x "$bin" ]]; then
    log "Downloading FastQC-Rust $RUST_VERSION"
    mkdir -p "$TOOLS_DIR/rust-$RUST_VERSION"
    curl -fsSL "https://github.com/ewels/FastQC-Rust/releases/download/$RUST_VERSION/$asset.tar.gz" \
      | tar xz -C "$TOOLS_DIR/rust-$RUST_VERSION"
    xattr -dr com.apple.quarantine "$TOOLS_DIR/rust-$RUST_VERSION" 2>/dev/null || true
  fi
  echo "$bin"
}

fetch_asprof() {
  local os lib
  case "$(uname -s)" in Darwin) os=macos ;; Linux) os=linux-$( [[ "$(uname -m)" == x86_64 ]] && echo x64 || echo arm64) ;; esac
  lib="$TOOLS_DIR/async-profiler-$ASPROF_VERSION-$os/lib/libasyncProfiler.$( [[ $os == macos ]] && echo dylib || echo so)"
  if [[ ! -f "$lib" ]]; then
    log "Downloading async-profiler $ASPROF_VERSION"
    local url="https://github.com/async-profiler/async-profiler/releases/download/v$ASPROF_VERSION/async-profiler-$ASPROF_VERSION-$os"
    if [[ $os == macos ]]; then
      curl -fsSL -o "$TOOLS_DIR/asprof.zip" "$url.zip" && unzip -qo "$TOOLS_DIR/asprof.zip" -d "$TOOLS_DIR" && rm "$TOOLS_DIR/asprof.zip"
      xattr -dr com.apple.quarantine "$TOOLS_DIR/async-profiler-$ASPROF_VERSION-$os" 2>/dev/null || true
    else
      curl -fsSL "$url.tar.gz" | tar xz -C "$TOOLS_DIR"
    fi
  fi
  echo "$lib"
}

# Runs the Java command twice more under async-profiler: thread-split collapsed stacks + HTML flamegraph
profile_java() {
  local input=$1 name=$2 out="$RESULTS_DIR/profile" fmt
  mkdir -p "$out/tmp"
  for fmt in "collapsed,threads,file=$out/$name.collapsed" "flamegraph,file=$out/$name.html"; do
    java -Xmx${MEMORY}m -Djava.awt.headless=true ${JAVA_OPTS:-} \
      -agentpath:"$ASPROF_LIB"=start,event=cpu,interval=1ms,$fmt \
      -cp "$(java_classpath)" -Dfastqc.output_dir="$out/tmp" -Dfastqc.quiet=true \
      ${THREADS:+-Dfastqc.threads=$THREADS} \
      uk.ac.babraham.FastQC.FastQCApplication "$input" >/dev/null 2>&1
  done
  rm -rf "$out/tmp"
  log "Profile: $out/$name.html (flamegraph), $out/$name.collapsed"
  # Top self-time frames
  awk '{n=$NF; t+=n; sub(/ [0-9]+$/,""); k=split($0,a,";"); s[a[k]]+=n}
       END{for(x in s) printf "  %5.1f%%  %s\n",100*s[x]/t,x}' "$out/$name.collapsed" | sort -rn | head -15 >&2
}

java_classpath() {
  local jars
  jars=$(ls "$REPO_DIR"/*.jar | tr '\n' ':')
  # classes first, then the repo root for Templates/Configuration resources
  echo "$CLASSES_DIR:$REPO_DIR:$jars"
}

compile_java() {
  log "Compiling Java sources -> $CLASSES_DIR"
  rm -rf "$CLASSES_DIR"; mkdir -p "$CLASSES_DIR"
  (cd "$SRC_DIR" && find uk net org $( [[ -d lib-rapidgzip/src/main/java ]] && echo lib-rapidgzip/src/main/java ) -name '*.java' > "$CLASSES_DIR/sources.txt")
  (cd "$SRC_DIR" && javac -nowarn -Xlint:none --release 21 -encoding UTF-8 -d "$CLASSES_DIR" \
      -cp "$(java_classpath)" @"$CLASSES_DIR/sources.txt" 2>&1 | grep -v '^Note:' >&2 || true)
  [[ -f "$CLASSES_DIR/uk/ac/babraham/FastQC/FastQCApplication.class" ]] \
    || { echo "Java compile failed" >&2; exit 1; }
}

# Prints "wall_s cpu_s rss_mb" for the given command, using /usr/bin/time
measure() {
  local tlog; tlog=$(mktemp)
  if [[ "$(uname -s)" == Darwin ]]; then
    /usr/bin/time -l -o "$tlog" timeout "${TIMEOUT:-600}" "$@" >"$tlog.out" 2>&1 || { echo "[bench] FAILED: $*" >&2; tail -20 "$tlog.out" >&2; return 1; }
    awk '/ real /{w=$1; u=$3; s=$5} /maximum resident set size/{r=$1/1048576}
         END{printf "%.2f %.2f %.0f\n", w, u+s, r}' "$tlog"
  else
    /usr/bin/time -f '%e %U %S %M' -o "$tlog" timeout "${TIMEOUT:-600}" "$@" >"$tlog.out" 2>&1 || { echo "[bench] FAILED: $*" >&2; tail -20 "$tlog.out" >&2; return 1; }
    awk '{printf "%.2f %.2f %.0f\n", $1, $2+$3, $4/1024}' "$tlog"
  fi
  rm -f "$tlog" "$tlog.out"
}

run_tool() {
  local tool=$1 input=$2 outdir=$3
  mkdir -p "$outdir"
  case "$tool" in
    java)
      measure java -Xmx${MEMORY}m -XX:+ExitOnOutOfMemoryError -Djava.awt.headless=true ${JAVA_OPTS:-} \
        -cp "$(java_classpath)" -Dfastqc.output_dir="$outdir" -Dfastqc.quiet=true \
        ${THREADS:+-Dfastqc.threads=$THREADS} \
        uk.ac.babraham.FastQC.FastQCApplication "$input" ;;
    rust)
      measure "$RUST_BIN" --quiet -o "$outdir" ${THREADS:+-t $THREADS} "$input" ;;
  esac
}

# Compare fastqc_data.txt of each tool against the first one (ignoring the version header)
check_outputs() {
  local name=$1 ref="" tool f
  for tool in ${TOOLS//,/ }; do
    f=$(find "$RESULTS_DIR/out/$tool/$name" -name fastqc_data.txt -o -name '*_fastqc.zip' | head -1)
    if [[ "$f" == *.zip ]]; then
      unzip -p "$f" '*/fastqc_data.txt' | grep -v '^##FastQC' > "$RESULTS_DIR/out/$tool/$name.data.txt"
    else
      grep -v '^##FastQC' "$f" > "$RESULTS_DIR/out/$tool/$name.data.txt"
    fi
    if [[ -z "$ref" ]]; then ref=$tool; continue; fi
    if cmp -s "$RESULTS_DIR/out/$ref/$name.data.txt" "$RESULTS_DIR/out/$tool/$name.data.txt"; then
      log "$name: $tool output identical to $ref"
    else
      log "$name: WARNING $tool output differs from $ref (diff $RESULTS_DIR/out/{$ref,$tool}/$name.data.txt)"
    fi
  done
}

main() {
  local datasets=("${@:-small}")
  [[ ",$TOOLS," == *",java,"* ]] && compile_java
  [[ ",$TOOLS," == *",rust,"* ]] && RUST_BIN=$(fetch_rust)
  [[ -n "$PROFILE" ]] && ASPROF_LIB=$(fetch_asprof)

  mkdir -p "$RESULTS_DIR"
  local csv="$RESULTS_DIR/results.tsv"
  [[ -f "$csv" ]] || printf 'timestamp\tcommit\tdataset\ttool\tthreads\trep\twall_s\tcpu_s\trss_mb\n' > "$csv"
  local commit; commit=$(git -C "$SRC_DIR" rev-parse --short HEAD)$(git -C "$SRC_DIR" diff --quiet || echo "-dirty")${LABEL:+:$LABEL}

  printf '\n%-8s %-6s %4s %10s %10s %9s\n' dataset tool rep wall_s cpu_s rss_mb
  for ds in "${datasets[@]}"; do
    local input; input=$(fetch_dataset "$ds")
    for tool in ${TOOLS//,/ }; do
      for rep in $(seq 1 "$REPS"); do
        local out="$RESULTS_DIR/out/$tool/$ds"
        rm -rf "$out"
        read -r wall cpu rss < <(run_tool "$tool" "$input" "$out")
        printf '%-8s %-6s %4s %10s %10s %9s\n' "$ds" "$tool" "$rep" "$wall" "$cpu" "$rss"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date +%FT%T)" "$commit" "$ds" "$tool" \
          "${THREADS:-default}" "$rep" "$wall" "$cpu" "$rss" >> "$csv"
      done
    done
    check_outputs "$ds"
    [[ -n "$PROFILE" ]] && profile_java "$input" "$ds"
  done
  log "Results appended to $csv"
}

main "$@"
