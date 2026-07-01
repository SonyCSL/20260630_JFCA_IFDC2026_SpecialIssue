#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${HOME}/work/nonaka/Chao1_Intensity}"
cd "${ROOT_DIR}"

Q1_TAG="${Q1_TAG:-stage2_q1_formula_brite_pathway_raw_gpupc2_20260629_000132}"
SSC_TAG="${SSC_TAG:-stage2_q2_nonformula9_ssc_raw_fourthroot_20260626_w4_171430}"
SSC_FIRST_CASE="${SSC_FIRST_CASE:-nonformula9_raw_ssc0935}"
EXPECTED_SHARDS="${EXPECTED_SHARDS:-34}"
POLL_SECONDS="${POLL_SECONDS:-300}"
MIN_AVAILABLE_GIB="${MIN_AVAILABLE_GIB:-40}"
WATCH_TAG="${WATCH_TAG:-watch_q1_then_resume_nonformula9_ssc_$(date '+%Y%m%d_%H%M%S')}"
LOG="logs/${WATCH_TAG}.log"
STATUS="logs/${WATCH_TAG}.status"

log_msg() {
  mkdir -p "$(dirname "${LOG}")"
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "${LOG}"
}

write_status() {
  printf '%s\n' "$*" > "${STATUS}"
}

available_gib() {
  free -g | awk '/^Mem:/ {print $7}'
}

q1_complete() {
  [[ -f "out_sensitivity/${Q1_TAG}/RUN_COMPLETE.txt" ]] || return 1
  if [[ -d "out_heatmap_cutoff_quantification_from_csv/${Q1_TAG}" ]]; then
    [[ -f "out_heatmap_cutoff_quantification_from_csv/${Q1_TAG}/RUN_COMPLETE.txt" ]] || return 1
  fi
  return 0
}

heavy_processes() {
  pgrep -af '[v]36_9_tomato_long_q012_patch_4gpu_v4_stage2_inext_q2fast.py|[T]omato_0324_v6_4gpu_stage2_inext.py|[s]tage2_q1_formula_brite_pathway_raw_gpupc2' || true
}

successful_shard_count() {
  local run_dir="out_sensitivity/${SSC_TAG}/upstream_${SSC_FIRST_CASE}"
  grep -Rsl '\[SUCCESS\] v36.9 output written to:' "${run_dir}/_shards"/worker_*/worker.log 2>/dev/null | wc -l
}

merge_first_case_if_needed() {
  local run_dir="out_sensitivity/${SSC_TAG}/upstream_${SSC_FIRST_CASE}"
  if [[ -s "${run_dir}/comparison_outer_v36_9_q012.csv" && -s "${run_dir}/comparison_cells_long_all_q012_v36_9.csv" && -s "${run_dir}/group_outer_v36_9_q012.csv" && -f "${run_dir}/RUN_COMPLETE.txt" ]]; then
    log_msg "first upstream case already merged: ${run_dir}"
    return 0
  fi
  local count
  count="$(successful_shard_count)"
  if [[ "${count}" != "${EXPECTED_SHARDS}" ]]; then
    log_msg "ERROR expected ${EXPECTED_SHARDS} successful shards for ${run_dir}, found ${count}"
    return 1
  fi
  log_msg "merge first upstream case from existing shards: ${run_dir}"
  .venv/bin/python scripts/merge_stage2_q2fast_shards_safe_20260625.py \
    --run-dir "${run_dir}" \
    --expected-shards "${EXPECTED_SHARDS}" \
    2>&1 | tee -a "${LOG}"
  log_msg "merge complete: ${run_dir}"
}

resume_ssc_wrapper() {
  local session="${SSC_SESSION:-stage2_q2_nonformula9_ssc_resume_after_q1_$(date '+%Y%m%d_%H%M%S')}"
  if tmux has-session -t "${session}" 2>/dev/null; then
    log_msg "resume session already exists: ${session}"
  else
    log_msg "launch resume SSC wrapper session=${session} run_tag=${SSC_TAG}"
    tmux new-session -d -s "${session}" \
      "cd '${ROOT_DIR}' && RUN_TAG='${SSC_TAG}' UPSTREAM_WORKERS=4 DOWNSTREAM_WORKERS=4 CUTOFF_CHUNK_SIZE=6 GENERATE_OVERLAYS=1 GENERATE_PPT=0 bash scripts/run_stage2_q2_nonformula9_ssc_raw_fourthroot_20260625.sh"
  fi
  write_status "LAUNCHED session=${session} run_tag=${SSC_TAG}"
}

main() {
  log_msg "START watcher tag=${WATCH_TAG} q1_tag=${Q1_TAG} ssc_tag=${SSC_TAG}"
  while true; do
    local avail heavy
    avail="$(available_gib)"
    heavy="$(heavy_processes)"
    if q1_complete && [[ -z "${heavy}" && -n "${avail}" && "${avail}" -ge "${MIN_AVAILABLE_GIB}" ]]; then
      log_msg "q1 complete and resources ok: avail_gib=${avail}"
      break
    fi
    log_msg "wait q1_complete=$(q1_complete && echo yes || echo no) avail_gib=${avail:-unknown} heavy_present=$([[ -n "${heavy}" ]] && echo yes || echo no)"
    write_status "WAITING q1_complete=$(q1_complete && echo yes || echo no) avail_gib=${avail:-unknown} heavy_present=$([[ -n "${heavy}" ]] && echo yes || echo no)"
    sleep "${POLL_SECONDS}"
  done

  merge_first_case_if_needed
  resume_ssc_wrapper
}

main "$@"
