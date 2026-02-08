@echo off
REM This script creates a desktop shortcut for the MST GUI

set "SCRIPT_DIR=%~dp0"
REM Remove trailing backslash to avoid escaping issues in PowerShell strings
if "%SCRIPT_DIR:~-1%"=="\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
set "SHORTCUT_NAME=MST Tool GUI"
set "DESKTOP=%USERPROFILE%\Desktop"

echo Creating desktop shortcut...

powershell -Command "$WS = New-Object -ComObject WScript.Shell; $SC = $WS.CreateShortcut('%DESKTOP%\%SHORTCUT_NAME%.lnk'); $SC.TargetPath = '%SCRIPT_DIR%\Launch_MST_GUI.bat'; $SC.WorkingDirectory = '%SCRIPT_DIR%'; $SC.Description = 'Launch SourceTracker2 MST Tool GUI'; $SC.Save()"

if exist "%DESKTOP%\%SHORTCUT_NAME%.lnk" (
    echo.
    echo Success! Shortcut created on your desktop: "%SHORTCUT_NAME%"
    echo.
    echo To change the icon:
    echo   1. Right-click the shortcut on your desktop
    echo   2. Select Properties
    echo   3. Click "Change Icon..."
    echo   4. Browse to choose an icon file
) else (
    echo.
    echo Failed to create shortcut. You can manually create one:
    echo   1. Right-click on Launch_MST_GUI.bat
    echo   2. Select "Create shortcut"
    echo   3. Move the shortcut to your desktop
)

echo.
pause
