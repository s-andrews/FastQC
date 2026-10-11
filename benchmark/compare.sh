#!/usr/bin/env bash
#
# Benchmark each commit of a branch stack:
#   MODE=step (default)  the base, then each commit of the stack in order, so each
#                        result compares with the previous one (what each stacked PR adds)
#   MODE=alone           the base, base + only that commit (cherry-picked), and the tip
#
# Usage: benchmark/compare.sh [base] [tip] [dataset]
#   defaults: base=master tip=HEAD dataset=small; REPS (default 5) and TOOLS (default java)
set -euo pipefail
BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"; REPO_DIR="$(dirname "$BENCH_DIR")"
DATASET="${3:-small}"
# resolve in the main repo, as HEAD inside the worktree is the worktree's own
BASE=$(git -C "$REPO_DIR" rev-parse "${1:-master}"); TIP=$(git -C "$REPO_DIR" rev-parse "${2:-HEAD}")
export REPS="${REPS:-5}" TOOLS="${TOOLS:-java}"
WT="$BENCH_DIR/build/worktree"

cleanup() { git -C "$REPO_DIR" worktree remove --force "$WT" 2>/dev/null || true; }
trap cleanup EXIT
cleanup

SUMMARY=$(mktemp)
# Prints the runs, and records the median wall/CPU/RSS of each variant for the summary
bench() {
  echo "== $1"
  SRC_DIR="$WT" LABEL="$1" "$BENCH_DIR/bench.sh" "$DATASET" 2>&1 | grep -E '^[a-z-]+ +[a-z]+ +[0-9]|identical|differs|FAILED' \
    | tee >(awk -v name="$1" '$2=="java"{w[n]=$4; c[n]=$5; r[n]=$6; n++}
        function med(a,  i,j,t,b){for(i=0;i<n;i++)b[i]=a[i]; for(i=0;i<n;i++)for(j=i+1;j<n;j++)if(b[j]<b[i]){t=b[i];b[i]=b[j];b[j]=t}; return b[int(n/2)]}
        END{if(n) printf "%-50s %7.2f %7.2f %6d\n", name, med(w), med(c), med(r)}' >> "$SUMMARY")
}
summary() {
  sleep 1
  printf '\n%-50s %7s %7s %6s\n' "median ($REPS runs, $DATASET)" wall_s cpu_s rss_mb
  cat "$SUMMARY"; rm -f "$SUMMARY"
}

git -C "$REPO_DIR" worktree add -q --detach "$WT" "$BASE"
bench baseline
if [[ "${MODE:-step}" == step ]]; then
  for c in $(git -C "$REPO_DIR" rev-list --reverse "$BASE..$TIP"); do
    git -C "$WT" checkout -q --detach "$c"
    bench "step:$(git -C "$REPO_DIR" log -1 --format=%f "$c" | cut -c1-40)"
  done
  summary
  exit 0
fi
for c in $(git -C "$REPO_DIR" rev-list --reverse "$BASE..$TIP"); do
  git -C "$WT" checkout -q --detach "$BASE"
  git -C "$WT" cherry-pick -n "$c" >/dev/null
  bench "alone:$(git -C "$REPO_DIR" log -1 --format=%f "$c" | cut -c1-40)"
  git -C "$WT" reset -q --hard
done
git -C "$WT" checkout -q --detach "$TIP"
bench stack
summary
