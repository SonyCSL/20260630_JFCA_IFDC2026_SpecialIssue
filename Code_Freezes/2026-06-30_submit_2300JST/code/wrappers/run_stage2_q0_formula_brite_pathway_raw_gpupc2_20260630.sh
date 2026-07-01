#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/home/nonaka/work/nonaka/Chao1_Intensity}"
PY="${PY:-${ROOT_DIR}/.venv/bin/python}"
RUN_TAG="${RUN_TAG:-stage2_q0_formula_brite_pathway_raw_gpupc2_$(date +%Y%m%d_%H%M%S)}"
OUT_PARENT="${OUT_PARENT:-${ROOT_DIR}/out_sensitivity/${RUN_TAG}}"
DOWN_PARENT="${DOWN_PARENT:-${ROOT_DIR}/out_heatmap_cutoff_quantification_from_csv/${RUN_TAG}}"
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/logs}"
LOG="${LOG_DIR}/${RUN_TAG}.log"

UPSTREAM_SCRIPT="${UPSTREAM_SCRIPT:-${ROOT_DIR}/v36_9_tomato_long_q012_patch_4gpu_v4_stage2_inext_q2fast.py}"
DOWNSTREAM_SCRIPT="${DOWNSTREAM_SCRIPT:-${ROOT_DIR}/Tomato_0324_v6_4gpu_stage2_inext.py}"
SAFE_MERGE_SCRIPT="${SAFE_MERGE_SCRIPT:-${ROOT_DIR}/scripts/merge_stage2_q2fast_shards_safe_20260625.py}"

DOMAINS="${DOMAINS:-Formula,Brite,Pathway}"
SUBSETS="${SUBSETS:-All,PM,SM}"
UPSTREAM_MODES="${UPSTREAM_MODES:-Annual,First3,Last3,Combo}"
DOWNSTREAM_MODES="${DOWNSTREAM_MODES:-2015,2016,2017,2018,2019,2020,Period1_2015_2017,Period2_2018_2020,Combo6yr}"
FIGURE_FILTER="${FIGURE_FILTER:-01,06,06_2}"

OUTER_REP="${TOMATO_OUTER_REP:-24}"
BOOT_MATRIX="${TOMATO_BOOT_MATRIX:-12}"
REFERENCE_TARGET_SC="${TOMATO_REFERENCE_TARGET_SC:-0.965}"
WORKERS="${WORKERS:-8}"
CUTOFF_CHUNK_SIZE="${CUTOFF_CHUNK_SIZE:-1}"
GPU_IDS="${GPU_IDS:-0,1,2,3}"

RUN_FULL="${RUN_FULL:-1}"
RUN_DOWNSTREAM="${RUN_DOWNSTREAM:-1}"
# Existing overlay generators are tuned for the q1/q2 main figure family and
# failed on the paired q1 run; keep q0 overlays as a separate explicit step.
GENERATE_OVERLAYS="${GENERATE_OVERLAYS:-0}"

mkdir -p "${OUT_PARENT}" "${DOWN_PARENT}" "${LOG_DIR}"

log_msg() {
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" | tee -a "${LOG}"
}

require_file() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    log_msg "missing required file: ${path}"
    return 1
  fi
}

has_upstream_csvs() {
  local dir="$1"
  [[ -s "${dir}/comparison_cells_long_all_q012_v36_9.csv" ]] &&
    [[ -s "${dir}/comparison_outer_v36_9_q012.csv" ]] &&
    [[ -s "${dir}/group_outer_v36_9_q012.csv" ]]
}

try_safe_merge() {
  local dir="$1"
  if [[ ! -f "${SAFE_MERGE_SCRIPT}" ]]; then
    log_msg "safe merge unavailable: ${SAFE_MERGE_SCRIPT}"
    return 1
  fi
  if [[ ! -d "${dir}/_shards" ]]; then
    log_msg "safe merge skipped; no shard root: ${dir}/_shards"
    return 1
  fi
  local shard_count
  shard_count=$(find "${dir}/_shards" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
  log_msg "safe merge check: dir=${dir}, shards=${shard_count}"
  if [[ "${shard_count}" -lt 199 ]]; then
    return 1
  fi
  "${PY}" "${SAFE_MERGE_SCRIPT}" --run-dir "${dir}" --expected-shards 199 2>&1 | tee -a "${LOG}"
}

run_upstream() {
  local out_dir="${OUT_PARENT}/upstream_raw_q0"
  if [[ "${RUN_FULL}" != "1" ]]; then
    log_msg "skip full upstream; RUN_FULL=${RUN_FULL}"
    return 0
  fi
  if has_upstream_csvs "${out_dir}" && [[ -f "${out_dir}/RUN_COMPLETE.txt" ]]; then
    log_msg "skip full upstream; existing complete output: ${out_dir}"
    return 0
  fi
  if try_safe_merge "${out_dir}"; then
    touch "${out_dir}/RUN_COMPLETE.txt"
    log_msg "finish full upstream via safe merge: ${out_dir}"
    return 0
  fi

  log_msg "start full q0 upstream: domains=${DOMAINS}, workers=${WORKERS}, chunk=${CUTOFF_CHUNK_SIZE}, outer=${OUTER_REP}, boot=${BOOT_MATRIX}"
  rm -rf "${out_dir}"
  mkdir -p "${out_dir}"
  set +e
  env \
    PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    TOMATO_POC_MODE=0 \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${UPSTREAM_MODES}" \
    TOMATO_Q_FILTER=q0 \
    TOMATO_OUTER_REP="${OUTER_REP}" \
    TOMATO_BOOT_MATRIX="${BOOT_MATRIX}" \
    TOMATO_REFERENCE_TARGET_SC="${REFERENCE_TARGET_SC}" \
    TOMATO_INTENSITY_TRANSFORM=raw \
    TOMATO_4GPU_NUM_WORKERS="${WORKERS}" \
    TOMATO_4GPU_GPU_IDS="${GPU_IDS}" \
    TOMATO_4GPU_CUTOFF_CHUNK_SIZE="${CUTOFF_CHUNK_SIZE}" \
    TOMATO_OUT_DIR="${out_dir}" \
    "${PY}" "${UPSTREAM_SCRIPT}" 2>&1 | tee -a "${LOG}"
  local rc=${PIPESTATUS[0]}
  set -e
  if [[ "${rc}" -ne 0 ]]; then
    log_msg "full upstream returned rc=${rc}; attempting safe merge before failing"
    if try_safe_merge "${out_dir}"; then
      touch "${out_dir}/RUN_COMPLETE.txt"
      log_msg "finish full upstream via safe merge after rc=${rc}: ${out_dir}"
      return 0
    fi
    return "${rc}"
  fi
  touch "${out_dir}/RUN_COMPLETE.txt"
  log_msg "finish full q0 upstream: ${out_dir}"
}

run_downstream() {
  local run_dir="${OUT_PARENT}/upstream_raw_q0"
  local out_dir="${DOWN_PARENT}/raw_q0_figures"
  if [[ "${RUN_DOWNSTREAM}" != "1" ]]; then
    log_msg "skip downstream; RUN_DOWNSTREAM=${RUN_DOWNSTREAM}"
    return 0
  fi
  if ! has_upstream_csvs "${run_dir}"; then
    log_msg "downstream blocked; missing upstream CSVs: ${run_dir}"
    return 1
  fi
  if [[ -f "${out_dir}/RUN_COMPLETE.txt" ]] &&
    [[ -n "$(find "${out_dir}" -name '01_reference_bridge_boxplots_paired.png' -print -quit 2>/dev/null)" ]] &&
    [[ -n "$(find "${out_dir}" -name '06_sweep_sc_heatmap.png' -print -quit 2>/dev/null)" ]] &&
    [[ -n "$(find "${out_dir}" -name '06_2_sweep_sc_large_effect_mask.png' -print -quit 2>/dev/null)" ]]; then
    log_msg "skip downstream; existing complete figures: ${out_dir}"
    return 0
  fi
  log_msg "start downstream q0: figures=${FIGURE_FILTER}"
  rm -rf "${out_dir}"
  mkdir -p "${out_dir}"
  env \
    TOMATO_RUN_DIR="${run_dir}" \
    TOMATO_OUT_DIR="${out_dir}" \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${DOWNSTREAM_MODES}" \
    TOMATO_Q_FILTER=q0 \
    TOMATO_METRIC_FILTER=CliffsDelta \
    TOMATO_FIGURE_FILTER="${FIGURE_FILTER}" \
    TOMATO_EXPORT_FORMATS=png \
    TOMATO_4GPU_NUM_WORKERS="${WORKERS}" \
    "${PY}" "${DOWNSTREAM_SCRIPT}" 2>&1 | tee -a "${LOG}"
  touch "${out_dir}/RUN_COMPLETE.txt"
  log_msg "finish downstream q0: ${out_dir}"
}

generate_overlays() {
  if [[ "${GENERATE_OVERLAYS}" != "1" ]]; then
    log_msg "skip q0 overlays; GENERATE_OVERLAYS=${GENERATE_OVERLAYS}"
    return 0
  fi
  log_msg "q0 overlays requested, but no q0-specific overlay step is defined in this wrapper"
  return 1
}

write_manifest() {
  {
    echo "RUN_TAG=${RUN_TAG}"
    echo "OUT_PARENT=${OUT_PARENT}"
    echo "DOWN_PARENT=${DOWN_PARENT}"
    echo "DOMAINS=${DOMAINS}"
    echo "SUBSETS=${SUBSETS}"
    echo "UPSTREAM_MODES=${UPSTREAM_MODES}"
    echo "DOWNSTREAM_MODES=${DOWNSTREAM_MODES}"
    echo "Q=q0"
    echo "TRANSFORM=raw"
    echo "OUTER_REP=${OUTER_REP}"
    echo "BOOT_MATRIX=${BOOT_MATRIX}"
    echo "REFERENCE_TARGET_SC=${REFERENCE_TARGET_SC}"
    echo "WORKERS=${WORKERS}"
    echo "CUTOFF_CHUNK_SIZE=${CUTOFF_CHUNK_SIZE}"
    echo "GPU_IDS=${GPU_IDS}"
    echo "LOG=${LOG}"
    date '+COMPLETED_AT=%F %T %Z'
  } > "${OUT_PARENT}/RUN_MANIFEST.txt"
}

main() {
  log_msg "run start: ${RUN_TAG}"
  require_file "${UPSTREAM_SCRIPT}"
  require_file "${DOWNSTREAM_SCRIPT}"
  run_upstream
  run_downstream
  generate_overlays
  write_manifest
  touch "${OUT_PARENT}/RUN_COMPLETE.txt" "${DOWN_PARENT}/RUN_COMPLETE.txt"
  log_msg "run complete: OUT_PARENT=${OUT_PARENT}; DOWN_PARENT=${DOWN_PARENT}"
}

main "$@"
