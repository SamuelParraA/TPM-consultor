#!/usr/bin/env bash
set -Eeuo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
python "${SCRIPT_DIR}/verify_bundle.py" --bundle "${BUNDLE_DIR}"
