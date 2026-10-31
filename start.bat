@echo off
REM ==============================================================================
REM PNN Discovery Pipeline - Turnkey One-Click Launcher (Windows)
REM ==============================================================================
title PNN Discovery Pipeline
cd /d "%~dp0"

echo [PNN] Starting PNN Perineuronal Net Drug Discovery Pipeline...
python run.py %*

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] Pipeline encountered an error (exit code: %ERRORLEVEL%).
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo [PNN] Turnkey execution succeeded! Dashboard opened in your web browser.
pause
