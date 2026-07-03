#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")"
exec python3 app.py --host 0.0.0.0 --port 8765
