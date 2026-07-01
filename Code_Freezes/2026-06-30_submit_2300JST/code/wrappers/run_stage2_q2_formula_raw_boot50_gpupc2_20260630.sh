#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${ROOT_DIR:-${HOME}/work/nonaka/Chao1_Intensity}"
cd "${ROOT_DIR}"

PY="${PY:-.venv/bin/python}"
RUN_TAG="${RUN_TAG:-stage2_q2_formula_raw_boot50_gpupc2_20260630_$(date +%H%M%S)}"
OUT_PARENT="${OUT_PARENT:-out_sensitivity/${RUN_TAG}}"
DOWN_PARENT="${DOWN_PARENT:-out_heatmap_cutoff_quantification_from_csv/${RUN_TAG}}"
LOG="${LOG:-logs/${RUN_TAG}.log}"

UPSTREAM_SCRIPT="${UPSTREAM_SCRIPT:-v36_9_tomato_long_q012_patch_4gpu_v4_stage2_inext_q2fast.py}"
DOWNSTREAM_SCRIPT="${DOWNSTREAM_SCRIPT:-Tomato_0324_v6_4gpu_stage2_inext.py}"
SAFE_MERGE_SCRIPT="${SAFE_MERGE_SCRIPT:-scripts/merge_stage2_q2fast_shards_safe_20260625.py}"

DOMAINS="${DOMAINS:-Formula}"
SUBSETS="${SUBSETS:-All,PM,SM}"
UPSTREAM_MODES="${UPSTREAM_MODES:-First3,Last3,Combo}"
DOWNSTREAM_MODES="${DOWNSTREAM_MODES:-Period1_2015_2017,Period2_2018_2020,Combo6yr}"
Q_FILTER="${Q_FILTER:-q2}"
TRANSFORM="${TRANSFORM:-raw}"
OUTER_REP="${TOMATO_OUTER_REP:-24}"
BOOT_MATRIX="${TOMATO_BOOT_MATRIX:-50}"
REFERENCE_TARGET_SC="${TOMATO_REFERENCE_TARGET_SC:-0.965}"
SWEEP_SC_GRID="${TOMATO_SWEEP_SC_GRID:-}"

WORKERS="${WORKERS:-8}"
DOWNSTREAM_WORKERS="${DOWNSTREAM_WORKERS:-8}"
GPU_IDS="${TOMATO_4GPU_GPU_IDS:-0,1,2,3}"
CUTOFF_CHUNK_SIZE="${CUTOFF_CHUNK_SIZE:-1}"
EXPECTED_SHARDS="${EXPECTED_SHARDS:-199}"
FIGURE_FILTER="${FIGURE_FILTER:-01,06,06_2}"

RUN_SMOKE="${RUN_SMOKE:-1}"
RUN_FULL="${RUN_FULL:-1}"
RUN_DOWNSTREAM="${RUN_DOWNSTREAM:-1}"
SMOKE_OUTER_REP="${SMOKE_OUTER_REP:-1}"
SMOKE_BOOT_MATRIX="${SMOKE_BOOT_MATRIX:-50}"
SMOKE_SWEEP_SC_GRID="${SMOKE_SWEEP_SC_GRID:-0.5,0.99,1.0}"
SMOKE_CUTOFF_START="${SMOKE_CUTOFF_START:-151}"
SMOKE_CUTOFF_END="${SMOKE_CUTOFF_END:-151}"

mkdir -p "${OUT_PARENT}" "${DOWN_PARENT}" "$(dirname "${LOG}")"

log_msg() {
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
  local dir="$1"
  if upstream_complete "${dir}"; then
    log_msg "skip safe merge; upstream complete: ${dir}"
    return 0
  fi
  local n
  n="$(shard_count "${dir}")"
  if [[ "${n}" != "${EXPECTED_SHARDS}" ]]; then
    log_msg "safe merge not applicable; shard_count=${n}, expected=${EXPECTED_SHARDS}"
    return 1
  fi
  log_msg "start safe merge: ${dir}"
  "${PY}" "${SAFE_MERGE_SCRIPT}" --run-dir "${dir}" --expected-shards "${EXPECTED_SHARDS}" 2>&1 | tee -a "${LOG}"
  touch "${dir}/RUN_COMPLETE.txt"
  log_msg "finish safe merge: ${dir}"
}

run_smoke() {
  local out_dir="${OUT_PARENT}/smoke_formula_raw_boot50_cutoff${SMOKE_CUTOFF_START}_outer${SMOKE_OUTER_REP}"
  if [[ "${RUN_SMOKE}" != "1" ]]; then
    log_msg "skip smoke; RUN_SMOKE=${RUN_SMOKE}"
    return 0
  fi
  if upstream_complete "${out_dir}"; then
    log_msg "skip smoke; complete output exists: ${out_dir}"
    return 0
  fi
  if [[ -e "${out_dir}" ]]; then
    log_msg "ERROR incomplete smoke output exists: ${out_dir}"
    return 2
  fi
  log_msg "start smoke: cutoff=${SMOKE_CUTOFF_START}-${SMOKE_CUTOFF_END}, outer=${SMOKE_OUTER_REP}, boot=${SMOKE_BOOT_MATRIX}, sweep=${SMOKE_SWEEP_SC_GRID}"
  env \
    PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    TOMATO_POC_MODE=0 \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${UPSTREAM_MODES}" \
    TOMATO_Q_FILTER="${Q_FILTER}" \
    TOMATO_OUTER_REP="${SMOKE_OUTER_REP}" \
    TOMATO_BOOT_MATRIX="${SMOKE_BOOT_MATRIX}" \
    TOMATO_REFERENCE_TARGET_SC="${REFERENCE_TARGET_SC}" \
    TOMATO_SWEEP_SC_GRID="${SMOKE_SWEEP_SC_GRID}" \
    TOMATO_INTENSITY_TRANSFORM="${TRANSFORM}" \
    TOMATO_4GPU_WORKER=1 \
    TOMATO_4GPU_GPU_ID=0 \
    TOMATO_4GPU_CUTOFF_START="${SMOKE_CUTOFF_START}" \
    TOMATO_4GPU_CUTOFF_END="${SMOKE_CUTOFF_END}" \
    TOMATO_OUT_DIR="${out_dir}" \
    "${PY}" "${UPSTREAM_SCRIPT}" 2>&1 | tee -a "${LOG}"
  touch "${out_dir}/RUN_COMPLETE.txt"
  log_msg "finish smoke: ${out_dir}"
}

run_upstream() {
  local out_dir="${OUT_PARENT}/upstream_formula_raw_boot50"
  if [[ "${RUN_FULL}" != "1" ]]; then
    log_msg "skip full upstream; RUN_FULL=${RUN_FULL}"
    return 0
  fi
  if upstream_complete "${out_dir}"; then
    log_msg "skip full upstream; complete output exists: ${out_dir}"
    return 0
  fi
  if [[ -d "${out_dir}/_shards" ]]; then
    safe_merge_if_possible "${out_dir}" && return 0
    log_msg "ERROR incomplete upstream shards cannot be safely reused: ${out_dir}"
    return 2
  fi
  if [[ -e "${out_dir}" ]]; then
    log_msg "ERROR incomplete upstream output exists without reusable shards: ${out_dir}"
    return 2
  fi

  log_msg "start full upstream: domains=${DOMAINS}, subsets=${SUBSETS}, modes=${UPSTREAM_MODES}, outer=${OUTER_REP}, boot=${BOOT_MATRIX}, workers=${WORKERS}, chunk=${CUTOFF_CHUNK_SIZE}"
  local sweep_env=()
  if [[ -n "${SWEEP_SC_GRID}" ]]; then
    sweep_env=(TOMATO_SWEEP_SC_GRID="${SWEEP_SC_GRID}")
  fi
  env \
    PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    TOMATO_POC_MODE=0 \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${UPSTREAM_MODES}" \
    TOMATO_Q_FILTER="${Q_FILTER}" \
    TOMATO_OUTER_REP="${OUTER_REP}" \
    TOMATO_BOOT_MATRIX="${BOOT_MATRIX}" \
    TOMATO_REFERENCE_TARGET_SC="${REFERENCE_TARGET_SC}" \
    TOMATO_INTENSITY_TRANSFORM="${TRANSFORM}" \
    TOMATO_4GPU_NUM_WORKERS="${WORKERS}" \
    TOMATO_4GPU_GPU_IDS="${GPU_IDS}" \
    TOMATO_4GPU_CUTOFF_CHUNK_SIZE="${CUTOFF_CHUNK_SIZE}" \
    TOMATO_OUT_DIR="${out_dir}" \
    "${sweep_env[@]}" \
    "${PY}" "${UPSTREAM_SCRIPT}" 2>&1 | tee -a "${LOG}"
  touch "${out_dir}/RUN_COMPLETE.txt"
  log_msg "finish full upstream: ${out_dir}"
}

run_downstream() {
  local run_dir="${OUT_PARENT}/upstream_formula_raw_boot50"
  local out_dir="${DOWN_PARENT}/formula_raw_boot50_figures"
  if [[ "${RUN_DOWNSTREAM}" != "1" ]]; then
    log_msg "skip downstream; RUN_DOWNSTREAM=${RUN_DOWNSTREAM}"
    return 0
  fi
  if ! upstream_complete "${run_dir}"; then
    log_msg "downstream blocked; missing complete upstream: ${run_dir}"
    return 1
  fi
  if downstream_complete "${out_dir}"; then
    log_msg "skip downstream; complete output exists: ${out_dir}"
    return 0
  fi
  if [[ -e "${out_dir}" ]]; then
    local archive="${DOWN_PARENT}/_archive_$(basename "${out_dir}")_$(date '+%H%M%S')"
    log_msg "archive incomplete downstream: ${out_dir} -> ${archive}"
    mkdir -p "$(dirname "${archive}")"
    mv "${out_dir}" "${archive}"
  fi

  log_msg "start downstream: modes=${DOWNSTREAM_MODES}, figures=${FIGURE_FILTER}"
  env \
    TOMATO_RUN_DIR="${run_dir}" \
    TOMATO_OUT_DIR="${out_dir}" \
    TOMATO_DOMAIN_FILTER="${DOMAINS}" \
    TOMATO_SUBSET_FILTER="${SUBSETS}" \
    TOMATO_MODE_FILTER="${DOWNSTREAM_MODES}" \
    TOMATO_Q_FILTER="${Q_FILTER}" \
    TOMATO_ESTIMATE_FILTER=Stage2_iNEXT_TD_m_est \
    TOMATO_METRIC_FILTER=CliffsDelta \
    TOMATO_FIGURE_FILTER="${FIGURE_FILTER}" \
    TOMATO_EXPORT_FORMATS=png \
    TOMATO_4GPU_NUM_WORKERS="${DOWNSTREAM_WORKERS}" \
    TOMATO_4GPU_GPU_IDS="${GPU_IDS}" \
    "${PY}" "${DOWNSTREAM_SCRIPT}" 2>&1 | tee -a "${LOG}"
  touch "${out_dir}/RUN_COMPLETE.txt"
  find "${out_dir}" -name '*.png' | sort > "${DOWN_PARENT}/png_file_list_formula_raw_boot50.txt"
  log_msg "finish downstream: ${out_dir}; png_count=$(wc -l < "${DOWN_PARENT}/png_file_list_formula_raw_boot50.txt")"
}

write_manifest() {
  cat > "${OUT_PARENT}/RUN_MANIFEST.txt" <<MANIFEST
RUN_TAG=${RUN_TAG}
HOST=$(hostname)
STARTED_AT=$(date '+%Y-%m-%d %H:%M:%S %Z')
ROOT_DIR=${ROOT_DIR}
OUT_PARENT=${OUT_PARENT}
DOWN_PARENT=${DOWN_PARENT}
UPSTREAM=${OUT_PARENT}/upstream_formula_raw_boot50
DOWNSTREAM=${DOWN_PARENT}/formula_raw_boot50_figures
DOMAINS=${DOMAINS}
SUBSETS=${SUBSETS}
UPSTREAM_MODES=${UPSTREAM_MODES}
DOWNSTREAM_MODES=${DOWNSTREAM_MODES}
Q_FILTER=${Q_FILTER}
TRANSFORM=${TRANSFORM}
OUTER_REP=${OUTER_REP}
BOOT_MATRIX=${BOOT_MATRIX}
REFERENCE_TARGET_SC=${REFERENCE_TARGET_SC}
SWEEP_SC_GRID=${SWEEP_SC_GRID:-default_0.01_to_0.99_plus_asymptotic}
WORKERS=${WORKERS}
DOWNSTREAM_WORKERS=${DOWNSTREAM_WORKERS}
GPU_IDS=${GPU_IDS}
CUTOFF_CHUNK_SIZE=${CUTOFF_CHUNK_SIZE}
FIGURE_FILTER=${FIGURE_FILTER}
BOOT12_BASELINE_CANDIDATE=/home/nonaka/work/nonaka/Chao1_Intensity/out_sensitivity/Full_v36_9_stage2_inext_q2fast_boot12_16w_20260622_004302
MANIFEST
}

main() {
  log_msg "START ${RUN_TAG}"
  write_manifest
  nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv,noheader | tee -a "${LOG}" || true
  run_smoke
  run_upstream
  run_downstream
  touch "${OUT_PARENT}/RUN_COMPLETE.txt" "${DOWN_PARENT}/RUN_COMPLETE.txt"
  nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv,noheader | tee -a "${LOG}" || true
  log_msg "DONE ${RUN_TAG}"
}

main "$@"
