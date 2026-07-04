#!/bin/bash
# gpu_validation_sweep.sh -- consolidated correctness sweep for the RG_0S GPU path.
#
# Runs one fixed-basis eigenvalue (EXPC_VALS G 201 on a Carbon 201-function basis)
# under every GPU configuration and records the energy. Because the basis is FIXED,
# all configs solve the identical eigenproblem, so every energy must agree within
# the CUDA atomicAdd noise floor (~1e-8 per eigenvalue). CPU matrix elements are
# exact (no atomicAdd), so the CPU reference is the ground truth for H,S.
#
# It grabs the "Energy:" line (printed right after the eigensolve) and kills the
# run before the slow expectation-value tail, so each config takes seconds.
#
# Submit as an sbatch job with 2 GPUs (for the 2-rank config). Writes results.txt.
#
# Usage (from a run dir containing the binary `ecg_bin` and base deck `base.txt`):
#   sbatch --account=hpcnc --partition=NVIDIA --nodes=1 --ntasks=2 --gres=gpu:2 \
#          --cpus-per-task=4 --time=00:40:00 --output=sweep-%j.out gpu_validation_sweep.sh
set -u

source /etc/profile.d/lmod.sh 2>/dev/null
module purge >/dev/null 2>&1
module load NVHPC/25.9-CUDA-12.9.1

HERE=$(pwd)
BIN=$HERE/ecg_bin
BASE=$HERE/base.txt          # 201-basis deck whose single step is "EXPC_VALS G 201" (with USE_GPU present)
RES=$HERE/results.txt
: > "$RES"

# run_cfg NAME NP USEGPU [ENV=VAL ...]
run_cfg() {
  local name=$1 np=$2 usegpu=$3; shift 3
  local d=$HERE/cfg_$name
  rm -rf "$d"; mkdir -p "$d"; cd "$d"
  cp "$BIN" ./ecg_bin
  cp "$BASE" ./inout.txt
  [ "$usegpu" = no ] && sed -i '/^ USE_GPU/d' inout.txt
  env "$@" mpirun --bind-to none -np "$np" ./ecg_bin > out.log 2>&1 &
  local pid=$!
  local i
  for i in $(seq 1 180); do
    grep -q 'Energy:' out.log 2>/dev/null && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 1
  done
  local E
  E=$(grep -m1 'Energy:' out.log 2>/dev/null | grep -oE '[-]37\.[0-9]+')
  local eig
  eig=$(grep -m1 'GPU eigensolver' out.log 2>/dev/null | wc -l)
  local chunk
  chunk=$(grep -m1 'GPU ME' out.log 2>/dev/null | sed 's/^ *//')
  # tear down this config's run so the GPUs are free for the next one
  kill "$pid" 2>/dev/null
  pkill -P "$pid" 2>/dev/null
  sleep 2
  pkill -f "$d/ecg_bin" 2>/dev/null
  sleep 2
  printf '%-14s np=%s  E=%-22s eig_gpu=%s  %s\n' "$name" "$np" "${E:-FAILED}" "$eig" "$chunk" | tee -a "$RES"
  cd "$HERE"
}

echo "=== RG_0S GPU validation sweep (fixed-basis EXPC_VALS G 201, Carbon K=201) ===" | tee -a "$RES"
nvidia-smi --query-gpu=index,name --format=csv,noheader | tee -a "$RES"
echo | tee -a "$RES"

run_cfg cpu_ref     1 no                              # CPU matrix elements + CPU DSYGVX (ground truth)
run_cfg gpu_me_a    1 yes                             # GPU ME + CPU eig (noise-floor run 1)
run_cfg gpu_me_b    1 yes                             # GPU ME + CPU eig (noise-floor run 2)
run_cfg gpu_eig     1 yes ECG_GPU_EIG=1               # GPU ME + GPU cuSOLVER eig
run_cfg persist_off 1 yes ECG_GPU_PERSIST=0           # Phase 1 off (legacy malloc/free)
run_cfg batch_64    1 yes ECG_GPU_BATCH=64            # Phase 2 chunking forced (318 chunks)
run_cfg rank2       2 yes                             # 2 ranks / 2 GPUs

echo | tee -a "$RES"
echo "Done. All energies should agree within ~1e-8 (atomicAdd noise floor)." | tee -a "$RES"
