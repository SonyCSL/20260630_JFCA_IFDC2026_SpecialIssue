#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${HOME}/work/nonaka/Chao1_Intensity}"
cd "${ROOT_DIR}"

PY="${PY:-.venv/bin/python}"
RUN_TAG="${RUN_TAG:-stage2_q2_nonformula9_ssc_raw_fourthroot_20260626_w4_171430}"
OUT_PARENT="${OUT_PARENT:-out_sensitivity/${RUN_TAG}}"
DOWN_PARENT="${DOWN_PARENT:-out_heatmap_cutoff_quantification_from_csv/${RUN_TAG}}"
LOG="logs/${RUN_TAG}_raw_urgent_20260630.log"

DOMAINS="${DOMAINS:-Brite,Pathway,ATC_L1,ATC_L2,ATC_L3,Network,Disease_NE,Disease_ICD11,Disease_PathCL}"
SUBSETS="${SUBSETS:-All,PM,SM}"
UPSTREAM_MODES="${UPSTREAM_MODES:-Annual,First3,Last3,Combo}"
DOWNSTREAM_MODES="${DOWNSTREAM_MODES:-Period1_2015_2017,Period2_2018_2020,Combo6yr}"
UPSTREAM_WORKERS="${UPSTREAM_WORKERS:-4}"
DOWNSTREAM_WORKERS="${DOWNSTREAM_WORKERS:-4}"
GPU_IDS="${TOMATO_4GPU_GPU_IDS:-0,1,2,3}"
CUTOFF_CHUNK_SIZE="${CUTOFF_CHUNK_SIZE:-6}"
OUTER_REP="${TOMATO_OUTER_REP:-24}"
BOOT_MATRIX="${TOMATO_BOOT_MATRIX:-12}"
FIGURE_FILTER="${FIGURE_FILTER:-01,06,06_2}"
EXPECTED_SHARDS="${EXPECTED_SHARDS:-34}"

log_msg() {
  mkdir -p "$(dirname "${LOG}")"
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "${LOG}"
}

upstream_complete() {
  local dir="$1"
  test -s "${dir}/comparison_outer_v36_9_q012.csv" &&
    test -s "${dir}/comparison_cells_long_all_q012_v36_9.csv" &&
    test -s "${dir}/group_outer_v36_9_q012.csv" &&
    test -f "${dir}/RUN_COMPLETE.txt"
}

downstream_complete() {
  local dir="$1"
  test -f "${dir}/RUN_COMPLETE.txt" &&
    test -n "$(find "${dir}" -name '01_reference_bridge_boxplots_paired.png' -print -quit 2>/dev/null)" &&
    test -n "$(find "${dir}" -name '06_sweep_sc_heatmap.png' -print -quit 2>/dev/null)" &&
    test -n "$(find "${dir}" -name '06_2_sweep_sc_large_effect_mask.png' -print -quit 2>/dev/null)"
}

shard_count() {
  local dir="$1"
  find "${dir}/_shards" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' '
}

safe_merge_if_possible() {
  local label="$1"
  local dir="${OUT_PARENT}/upstream_${label}"
  if upstream_complete "${dir}"; then
    log_msg "skip safe merge ${label}; complete upstream exists"
    return 0
  fi
  local n
  n="$(shard_count "${dir}")"
  if [[ "${n}" != "${EXPECTED_SHARDS}" ]]; then
    log_msg "safe merge not applicable for ${label}; shard_count=${n}, expected=${EXPECTED_SHARDS}"
    return 1
  fi
  log_msg "start safe merge ${label}: ${dir}"
  "${PY}" scripts/merge_stage2_q2fast_shards_safe_20260625.py \
    --run-dir "${dir}" \
    --expected-shards "${EXPECTED_SHARDS}" \
    2>&1 | tee -a "${LOG}"
  log_msg "finish safe merge ${label}: ${dir}"
}

run_upstream() {
  local label="$1"
  local ssc="$2"
  local out_dir="${OUT_PARENT}/upstream_${label}"

  if upstream_complete "${out_dir}"; then
    log_msg "skip upstream ${label}; complete output exists"
    return 0
  fi
  if test -d "${out_dir}/_shards"; then
    safe_merge_if_possible "${label}" && return 0
  fi
  if test -e "${out_dir}"; then
    log_msg "ERROR incomplete upstream exists and cannot be safely reused: ${out_dir}"
    return 2
  fi

  log_msg "start upstream ${label}: raw, reference_sc=${ssc}, workers=${UPSTREAM_WORKERS}, gpu_ids=${GPU_IDS}"
  env \
    TOMATO_POC_MODE=0 \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${UPSTREAM_MODES}" \
    TOMATO_Q_FILTER=q2 \
    TOMATO_OUTER_REP="${OUTER_REP}" \
    TOMATO_BOOT_MATRIX="${BOOT_MATRIX}" \
    TOMATO_REFERENCE_TARGET_SC="${ssc}" \
    TOMATO_INTENSITY_TRANSFORM=raw \
    TOMATO_4GPU_NUM_WORKERS="${UPSTREAM_WORKERS}" \
    TOMATO_4GPU_GPU_IDS="${GPU_IDS}" \
    TOMATO_4GPU_CUTOFF_CHUNK_SIZE="${CUTOFF_CHUNK_SIZE}" \
    TOMATO_OUT_DIR="${out_dir}" \
    "${PY}" v36_9_tomato_long_q012_patch_4gpu_v4_stage2_inext_q2fast.py \
    2>&1 | tee -a "${LOG}"
  touch "${out_dir}/RUN_COMPLETE.txt"
  log_msg "finish upstream ${label}: ${out_dir}"
}

run_downstream() {
  local label="$1"
  local run_dir="${OUT_PARENT}/upstream_${label}"
  local out_dir="${DOWN_PARENT}/${label}"

  if downstream_complete "${out_dir}"; then
    log_msg "skip downstream ${label}; complete output exists"
    return 0
  fi
  if test -e "${out_dir}" && ! downstream_complete "${out_dir}"; then
    local archive="${DOWN_PARENT}/_archive_${RUN_TAG}/$(basename "${out_dir}")_$(date '+%H%M%S')"
    mkdir -p "$(dirname "${archive}")"
    log_msg "archive incomplete downstream ${label}: ${out_dir} -> ${archive}"
    mv "${out_dir}" "${archive}"
  fi

  log_msg "start downstream ${label}: periods=${DOWNSTREAM_MODES}, figures=${FIGURE_FILTER}"
  env \
    TOMATO_RUN_DIR="${run_dir}" \
    TOMATO_OUT_DIR="${out_dir}" \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${DOWNSTREAM_MODES}" \
    TOMATO_Q_FILTER=q2 \
    TOMATO_ESTIMATE_FILTER=Stage2_iNEXT_TD_m_est \
    TOMATO_METRIC_FILTER=CliffsDelta \
    TOMATO_FIGURE_FILTER="${FIGURE_FILTER}" \
    TOMATO_EXPORT_FORMATS=png \
    TOMATO_4GPU_NUM_WORKERS="${DOWNSTREAM_WORKERS}" \
    TOMATO_4GPU_GPU_IDS="${GPU_IDS}" \
    "${PY}" Tomato_0324_v6_4gpu_stage2_inext.py \
    2>&1 | tee -a "${LOG}"
  touch "${out_dir}/RUN_COMPLETE.txt"
  log_msg "finish downstream ${label}: ${out_dir}"
}

run_case() {
  local label="$1"
  local ssc="$2"
  run_upstream "${label}" "${ssc}"
  run_downstream "${label}"
}

main() {
  mkdir -p "${OUT_PARENT}" "${DOWN_PARENT}" logs
  log_msg "START raw urgent SSC ${RUN_TAG}"
  log_msg "scope domains=${DOMAINS}; subsets=${SUBSETS}; downstream_modes=${DOWNSTREAM_MODES}; upstream_modes=${UPSTREAM_MODES}"
  nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv,noheader | tee -a "${LOG}" || true
  run_case "nonformula9_raw_ssc0935" "0.935"
  run_case "nonformula9_raw_ssc0990" "0.990"
  find "${DOWN_PARENT}" -name '*.png' | sort > "${DOWN_PARENT}/png_file_list_raw_urgent_20260630.txt"
  touch "${OUT_PARENT}/RUN_COMPLETE_RAW_URGENT_20260630.txt" "${DOWN_PARENT}/RUN_COMPLETE_RAW_URGENT_20260630.txt"
  log_msg "DONE raw urgent SSC ${RUN_TAG} png_count=$(wc -l < "${DOWN_PARENT}/png_file_list_raw_urgent_20260630.txt")"
}

main "$@"
