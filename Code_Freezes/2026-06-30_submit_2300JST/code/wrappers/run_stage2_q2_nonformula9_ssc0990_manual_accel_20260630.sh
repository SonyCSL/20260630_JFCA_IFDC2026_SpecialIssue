#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${HOME}/work/nonaka/Chao1_Intensity}"
cd "${ROOT_DIR}"

PY="${PY:-.venv/bin/python}"
RUN_TAG="${RUN_TAG:-stage2_q2_nonformula9_ssc_raw_fourthroot_20260626_w4_171430}"
CASE_LABEL="${CASE_LABEL:-nonformula9_raw_ssc0990}"
UP_DIR="out_sensitivity/${RUN_TAG}/upstream_${CASE_LABEL}"
DOWN_DIR="out_heatmap_cutoff_quantification_from_csv/${RUN_TAG}/${CASE_LABEL}"
LOG="logs/${RUN_TAG}_${CASE_LABEL}_manual_accel_20260630.log"

DOMAINS="${DOMAINS:-Brite,Pathway,ATC_L1,ATC_L2,ATC_L3,Network,Disease_NE,Disease_ICD11,Disease_PathCL}"
SUBSETS="${SUBSETS:-All,PM,SM}"
UPSTREAM_MODES="${UPSTREAM_MODES:-Annual,First3,Last3,Combo}"
DOWNSTREAM_MODES="${DOWNSTREAM_MODES:-Period1_2015_2017,Period2_2018_2020,Combo6yr}"
GPU_IDS_CSV="${TOMATO_4GPU_GPU_IDS:-0,1,2,3}"
CUTOFF_CHUNK_SIZE="${CUTOFF_CHUNK_SIZE:-6}"
START_CUTOFF="${START_CUTOFF:-25}"
END_CUTOFF="${END_CUTOFF:-199}"
START_WORKER_IDX="${START_WORKER_IDX:-5}"
OUTER_REP="${TOMATO_OUTER_REP:-24}"
BOOT_MATRIX="${TOMATO_BOOT_MATRIX:-12}"
BASE_SEED_DEFAULT="${TOMATO_BASE_SEED:-12345}"
EXPECTED_SHARDS="${EXPECTED_SHARDS:-34}"
ADOPT_INCOMPLETE="${ADOPT_INCOMPLETE:-0}"
SKIP_MERGE_AFTER_QUEUE="${SKIP_MERGE_AFTER_QUEUE:-0}"
LAUNCH_REC=""

log_msg() {
  mkdir -p "$(dirname "${LOG}")"
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "${LOG}" >&2
}

required_worker_files_present() {
  local dir="$1"
  test -s "${dir}/comparison_outer_v36_9_q012.csv" &&
    test -s "${dir}/comparison_cells_long_all_q012_v36_9.csv" &&
    test -s "${dir}/group_outer_v36_9_q012.csv" &&
    test -s "${dir}/group_summary_v36_9_q012.csv" &&
    test -s "${dir}/comparison_agg_v36_9_q012.csv" &&
    grep -q '\[SUCCESS\]' "${dir}/worker.log" 2>/dev/null
}

launch_one() {
  local worker_idx="$1"
  local gpu_id="$2"
  local cutoff_start="$3"
  local cutoff_end="$4"
  local shard_dir="${UP_DIR}/_shards/worker_$(printf '%03d' "${worker_idx}")_gpu${gpu_id}_cutoffs_$(printf '%03d' "${cutoff_start}")_$(printf '%03d' "${cutoff_end}")"

  if required_worker_files_present "${shard_dir}"; then
    log_msg "skip worker=${worker_idx}; complete shard exists: ${shard_dir}"
    LAUNCH_REC=""
    return 0
  fi
  if test -e "${shard_dir}"; then
    if [[ "${ADOPT_INCOMPLETE}" == "1" ]]; then
      log_msg "adopt existing incomplete shard worker=${worker_idx} gpu=${gpu_id} cutoffs=${cutoff_start}-${cutoff_end}: ${shard_dir}"
      LAUNCH_REC="0:${worker_idx}:${gpu_id}:${cutoff_start}:${cutoff_end}:${shard_dir}"
      return 0
    fi
    local archive="${UP_DIR}/_archive_manual_accel_20260630/$(basename "${shard_dir}")_$(date '+%H%M%S')"
    mkdir -p "$(dirname "${archive}")"
    log_msg "archive incomplete shard ${shard_dir} -> ${archive}"
    mv "${shard_dir}" "${archive}"
  fi
  mkdir -p "${shard_dir}"
  log_msg "launch worker=${worker_idx} gpu=${gpu_id} cutoffs=${cutoff_start}-${cutoff_end}"
  env \
    TOMATO_4GPU_WORKER=1 \
    TOMATO_4GPU_GPU_ID="${gpu_id}" \
    TOMATO_4GPU_CUTOFF_START="${cutoff_start}" \
    TOMATO_4GPU_CUTOFF_END="${cutoff_end}" \
    TOMATO_4GPU_NUM_WORKERS=4 \
    TOMATO_BASE_SEED="$((BASE_SEED_DEFAULT + (cutoff_start - 1) * 1000))" \
    TOMATO_POC_MODE=0 \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${UPSTREAM_MODES}" \
    TOMATO_Q_FILTER=q2 \
    TOMATO_OUTER_REP="${OUTER_REP}" \
    TOMATO_BOOT_MATRIX="${BOOT_MATRIX}" \
    TOMATO_REFERENCE_TARGET_SC=0.990 \
    TOMATO_INTENSITY_TRANSFORM=raw \
    TOMATO_OUT_DIR="${shard_dir}" \
    "${PY}" v36_9_tomato_long_q012_patch_4gpu_v4_stage2_inext_q2fast.py \
    > "${shard_dir}/worker.log" 2>&1 &
  LAUNCH_REC="$!:${worker_idx}:${gpu_id}:${cutoff_start}:${cutoff_end}:${shard_dir}"
}

run_manual_queue() {
  IFS=',' read -r -a gpu_ids <<< "${GPU_IDS_CSV}"
  local -a ranges=()
  local start="${START_CUTOFF}"
  local worker_idx="${START_WORKER_IDX}"
  while [[ "${start}" -le "${END_CUTOFF}" ]]; do
    local end=$((start + CUTOFF_CHUNK_SIZE - 1))
    if [[ "${end}" -gt "${END_CUTOFF}" ]]; then
      end="${END_CUTOFF}"
    fi
    ranges+=("${worker_idx}:${start}:${end}")
    worker_idx=$((worker_idx + 1))
    start=$((end + 1))
  done

  local next=0
  local -a active=()
  log_msg "manual queue start: ranges=${#ranges[@]} gpu_ids=${GPU_IDS_CSV}"
  while [[ "${next}" -lt "${#ranges[@]}" || "${#active[@]}" -gt 0 ]]; do
    while [[ "${next}" -lt "${#ranges[@]}" && "${#active[@]}" -lt "${#gpu_ids[@]}" ]]; do
      IFS=':' read -r wid cstart cend <<< "${ranges[$next]}"
      local gpu="${gpu_ids[$((next % ${#gpu_ids[@]}))]}"
      LAUNCH_REC=""
      launch_one "${wid}" "${gpu}" "${cstart}" "${cend}" || exit 1
      if [[ "${LAUNCH_REC}" == *:* ]]; then
        active+=("${LAUNCH_REC}")
      fi
      next=$((next + 1))
    done

    local -a still=()
    local rec pid wid gpu cstart cend shard_dir
    for rec in "${active[@]}"; do
      IFS=':' read -r pid wid gpu cstart cend shard_dir <<< "${rec}"
      if required_worker_files_present "${shard_dir}"; then
        log_msg "finish worker=${wid} gpu=${gpu} cutoffs=${cstart}-${cend}"
      elif [[ "${pid}" == "0" ]] || kill -0 "${pid}" 2>/dev/null; then
        still+=("${rec}")
      else
        log_msg "ERROR worker=${wid} process ended without complete outputs: ${shard_dir}"
        tail -n 40 "${shard_dir}/worker.log" | tee -a "${LOG}" || true
        exit 1
      fi
    done
    active=("${still[@]}")
    sleep 20
  done
  log_msg "manual queue complete"
}

wait_all_expected_shards() {
  log_msg "wait for expected shard completeness: expected=${EXPECTED_SHARDS}"
  local complete
  while true; do
    complete=0
    local dir
    for dir in "${UP_DIR}"/_shards/worker_*_cutoffs_*; do
      test -d "${dir}" || continue
      if required_worker_files_present "${dir}"; then
        complete=$((complete + 1))
      fi
    done
    log_msg "complete_shards=${complete}/${EXPECTED_SHARDS}"
    [[ "${complete}" -ge "${EXPECTED_SHARDS}" ]] && break
    sleep 60
  done
}

run_merge_and_downstream() {
  log_msg "start safe merge ${UP_DIR}"
  wait_all_expected_shards
  "${PY}" scripts/merge_stage2_q2fast_shards_safe_20260625.py \
    --run-dir "${UP_DIR}" \
    --expected-shards "${EXPECTED_SHARDS}" 2>&1 | tee -a "${LOG}"

  if test -e "${DOWN_DIR}" && ! test -f "${DOWN_DIR}/RUN_COMPLETE.txt"; then
    local archive="out_heatmap_cutoff_quantification_from_csv/${RUN_TAG}/_archive_manual_accel_20260630/$(basename "${DOWN_DIR}")_$(date '+%H%M%S')"
    mkdir -p "$(dirname "${archive}")"
    log_msg "archive incomplete downstream ${DOWN_DIR} -> ${archive}"
    mv "${DOWN_DIR}" "${archive}"
  fi
  log_msg "start downstream ${CASE_LABEL}"
  env \
    TOMATO_RUN_DIR="${UP_DIR}" \
    TOMATO_OUT_DIR="${DOWN_DIR}" \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${DOWNSTREAM_MODES}" \
    TOMATO_Q_FILTER=q2 \
    TOMATO_ESTIMATE_FILTER=Stage2_iNEXT_TD_m_est \
    TOMATO_METRIC_FILTER=CliffsDelta \
    TOMATO_FIGURE_FILTER=01,06,06_2 \
    TOMATO_EXPORT_FORMATS=png \
    TOMATO_4GPU_NUM_WORKERS=4 \
    TOMATO_4GPU_GPU_IDS="${GPU_IDS_CSV}" \
    "${PY}" Tomato_0324_v6_4gpu_stage2_inext.py 2>&1 | tee -a "${LOG}"
  touch "${DOWN_DIR}/RUN_COMPLETE.txt"
  log_msg "finish downstream ${CASE_LABEL}: ${DOWN_DIR}"
}

main() {
  mkdir -p "${UP_DIR}/_shards" "${DOWN_DIR}" logs
  log_msg "START manual accel ${CASE_LABEL} run_tag=${RUN_TAG}"
  run_manual_queue
  if [[ "${SKIP_MERGE_AFTER_QUEUE}" == "1" ]]; then
    log_msg "SKIP_MERGE_AFTER_QUEUE=1; leaving shards complete check/merge/downstream for a separate controlled step"
    return 0
  fi
  run_merge_and_downstream
  log_msg "DONE manual accel ${CASE_LABEL}"
}

main "$@"
