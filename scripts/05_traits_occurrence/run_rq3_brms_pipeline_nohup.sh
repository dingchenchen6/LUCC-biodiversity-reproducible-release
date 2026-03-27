#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LOG_DIR="${RQ3_BRMS_LOG_DIR:-${REPO_ROOT}/results/traits_fg/rq3_brms_server_runner}"
RUN_STAMP="$(date +%Y%m%d_%H%M%S)"

mkdir -p "${LOG_DIR}"
cd "${REPO_ROOT}"

MASTER_LOG="${LOG_DIR}/${RUN_STAMP}_rq3_brms_pipeline_master.log"
PID_FILE="${LOG_DIR}/${RUN_STAMP}_rq3_brms_pipeline.pid"

nohup Rscript "${SCRIPT_DIR}/00_run_rq3_brms_server.R" > "${MASTER_LOG}" 2>&1 &
RUN_PID=$!
echo "${RUN_PID}" > "${PID_FILE}"

echo "RQ3 brms pipeline started in background."
echo "PID: ${RUN_PID}"
echo "Master log: ${MASTER_LOG}"
echo "PID file: ${PID_FILE}"
