#!/usr/bin/env bash
#
# Compares FastQC versions the way the FastQC-Rust performance report does
# (https://ewels.github.io/FastQC-Rust/about/performance/report/): each tool with its default
# settings, through its own launcher, one file at a time, one run per file.
#
# Usage: benchmark/report.sh file.fastq.gz...
#
# Tools (TOOLS=java012,master,opt,rust):
#   java012  FastQC 0.12.1 as installed (`fastqc` on the PATH, eg Homebrew)
#   master   FastQC 0.13.0, the master branch, through its `fastqc` launcher
#   opt      the optimised branch (OPT_REF, default HEAD), through its `fastqc` launcher
#   optbest  the same build run with java directly, without the launcher's -t 4 and
#            -XX:ParallelGCThreads=1, so it uses as many threads as it wants
#   rust     FastQC-Rust v1.1.0 (downloaded by bench.sh)
#
# Env: REPS (default 1), TIMEOUT (seconds per run, default 900), OPT_REF.
set -euo pipefail

BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$BENCH_DIR")"
OUT_DIR="$BENCH_DIR/results/report"
DIST_DIR="$BENCH_DIR/build/dist"
TOOLS="${TOOLS:-java012,master,opt,rust}"
REPS="${REPS:-1}"
TIMEOUT="${TIMEOUT:-900}"
OPT_REF="${OPT_REF:-HEAD}"
RUST_BIN="$BENCH_DIR/tools/rust-v1.1.0/fastqc-macos-aarch64/fastqc"

log() { echo "[report] $*" >&2; }

# An installed-style FastQC: compiled classes and resources next to the jars and the launcher
make_dist() {
  local name=$1 ref=$2 dir="$DIST_DIR/$1" wt="$BENCH_DIR/build/worktree-$1"
  local commit; commit=$(git -C "$REPO_DIR" rev-parse "$ref")
  if [[ -f "$dir/.commit" && "$(cat "$dir/.commit")" == "$commit" ]]; then echo "$dir"; return; fi
  log "Building $name from $ref ($commit)"
  rm -rf "$dir"; mkdir -p "$dir"
  git -C "$REPO_DIR" worktree remove --force "$wt" 2>/dev/null || true
  git -C "$REPO_DIR" worktree add -q --detach "$wt" "$commit"
  (cd "$wt" && find uk net org $( [[ -d lib-rapidgzip/src/main/java ]] && echo lib-rapidgzip/src/main/java ) -name '*.java' > "$dir/.sources" \
    && javac -nowarn -Xlint:none -encoding UTF-8 -d "$dir" -cp "$(ls "$REPO_DIR"/*.jar | tr '\n' ':')" @"$dir/.sources" 2>&1 | grep -v '^Note:' >&2 || true)
  (cd "$wt" && find uk net org -type f ! -name '*.java' | while read -r f; do mkdir -p "$dir/$(dirname "$f")"; cp "$f" "$dir/$f"; done)
  cp -R "$wt/Configuration" "$wt/Templates" "$wt/fastqc" "$dir/"
  cp "$REPO_DIR"/*.jar "$dir/"
  rm -f "$dir/.sources"
  git -C "$REPO_DIR" worktree remove --force "$wt"
  [[ -f "$dir/uk/ac/babraham/FastQC/FastQCApplication.class" ]] || { echo "build of $name failed" >&2; exit 1; }
  echo "$commit" > "$dir/.commit"
  echo "$dir"
}

command_for() {
  local tool=$1 out=$2 file=$3
  case "$tool" in
    java012) echo "fastqc --quiet --outdir $out $file" ;;
    master|opt) echo "python3 $DIST_DIR/$tool/fastqc --quiet --outdir $out $file" ;;
    optbest) echo "java -Xmx512m -XX:+ExitOnOutOfMemoryError -Djava.awt.headless=true -Dfastqc.quiet=true -Dfastqc.output_dir=$out -cp $DIST_DIR/opt:$(ls "$REPO_DIR"/*.jar | tr '\n' ':') uk.ac.babraham.FastQC.FastQCApplication $file" ;;
    rust) echo "$RUST_BIN --quiet -o $out $file" ;;
  esac
}

main() {
  [[ $# -gt 0 ]] || { echo "usage: $0 file.fastq.gz..." >&2; exit 1; }
  [[ ",$TOOLS," == *",master,"* ]] && make_dist master master >/dev/null
  [[ ",$TOOLS," == *",opt,"* || ",$TOOLS," == *",optbest,"* ]] && make_dist opt "$OPT_REF" >/dev/null
  [[ ",$TOOLS," == *",java012,"* ]] && { fastqc --version | grep -q "0.12.1" || { echo "fastqc 0.12.1 not on the PATH" >&2; exit 1; }; }
  [[ ",$TOOLS," == *",rust,"* && ! -x "$RUST_BIN" ]] && { echo "run benchmark/bench.sh once to download FastQC-Rust" >&2; exit 1; }

  mkdir -p "$OUT_DIR"
  local tsv="$OUT_DIR/report.tsv"
  [[ -f "$tsv" ]] || printf 'timestamp\tfile\ttool\trep\twall_s\tcpu_s\trss_mb\n' > "$tsv"
  printf '%-40s %-8s %3s %9s %9s %7s\n' file tool rep wall_s cpu_s rss_mb
  for file in "$@"; do
    local name; name=$(basename "$file")
    for tool in ${TOOLS//,/ }; do
      for rep in $(seq 1 "$REPS"); do
        local out="$OUT_DIR/out/$tool/${name%%.*}"; rm -rf "$out"; mkdir -p "$out"
        local tlog; tlog=$(mktemp)
        # shellcheck disable=SC2046
        if ! /usr/bin/time -l -o "$tlog" timeout "$TIMEOUT" $(command_for "$tool" "$out" "$file") > "$tlog.out" 2>&1; then
          echo "FAILED $tool $name" >&2; tail -5 "$tlog.out" >&2; rm -f "$tlog" "$tlog.out"; continue
        fi
        read -r wall cpu rss < <(awk '/ real /{w=$1; u=$3; s=$5} /maximum resident set size/{r=$1/1048576} END{printf "%.2f %.2f %.0f\n", w, u+s, r}' "$tlog")
        rm -f "$tlog" "$tlog.out"
        printf '%-40s %-8s %3s %9s %9s %7s\n' "$name" "$tool" "$rep" "$wall" "$cpu" "$rss"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date +%FT%T)" "$name" "$tool" "$rep" "$wall" "$cpu" "$rss" >> "$tsv"
      done
    done
    # Results must be identical to Rust's (ignoring the version line)
    local ref="$OUT_DIR/out/rust/${name%%.*}"
    for tool in ${TOOLS//,/ }; do
      [[ $tool == rust || ! -d "$ref" ]] && continue
      local a b
      a=$(unzip -p "$ref"/*_fastqc.zip '*/fastqc_data.txt' 2>/dev/null | grep -v '^##FastQC' | md5 -q 2>/dev/null || true)
      b=$(unzip -p "$OUT_DIR/out/$tool/${name%%.*}"/*_fastqc.zip '*/fastqc_data.txt' 2>/dev/null | grep -v '^##FastQC' | md5 -q 2>/dev/null || true)
      [[ -n "$a" && "$a" == "$b" ]] && log "$name: $tool output identical to rust" || log "$name: $tool output DIFFERS from rust"
    done
  done
  log "Results appended to $tsv"
}

main "$@"
