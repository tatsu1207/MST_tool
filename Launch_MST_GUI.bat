@echo off
title SourceTracker2 MST GUI
echo ========================================
echo   SourceTracker2 MST Tool - Starting...
echo ========================================
echo.

echo Starting Streamlit GUI...
echo.
echo The app will open at: http://localhost:8501
echo.
echo Press Ctrl+C to stop the server when done.
echo ========================================
echo.

REM Start the browser after a short delay (in background)
start "" cmd /c "timeout /t 4 /nobreak >nul && start http://localhost:8501"

REM Convert Windows path to WSL path and run the Streamlit app
for /f "usebackq tokens=*" %%i in (`wsl wslpath -u "%~dp0"`) do set "WSL_DIR=%%i"
wsl bash -l -c "cd '%WSL_DIR%' && ./run_app.sh"

pause
