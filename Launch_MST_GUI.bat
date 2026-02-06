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

REM Run the Streamlit app in WSL using the run_app.sh script
wsl bash -l -c "cd /mnt/c/Users/User/Desktop/MST_tool && ./run_app.sh"

pause
