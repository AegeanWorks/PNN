#!/usr/bin/env bash
# ==============================================================================
# PNN Discovery Pipeline - Turnkey One-Click Launcher (Linux / macOS)
# ==============================================================================
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

echo "[PNN] Starting PNN Perineuronal Net Drug Discovery Pipeline..."
if command -v python3 &>/dev/null; then
    python3 run.py "$@"
else
    python run.py "$@"
fi

echo "[PNN] Turnkey execution succeeded! Dashboard opened in your web browser."
