#!/usr/bin/env bash
# Linux/macOS Bash entry point. The Python standard-library helper handles
# resumable downloads, process accounting, locking and atomic completion files.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$SCRIPT_DIR/copsac_batch.py" "$@"
