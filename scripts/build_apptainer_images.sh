#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

exec bash "${ROOT_SCRIPTS_DIR}/build_apptainer_images.sh" "$@"

