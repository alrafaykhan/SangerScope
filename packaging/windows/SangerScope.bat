@echo off
rem Starts SangerScope and opens it in the default browser.
title SangerScope
cd /d "%~dp0"
"%~dp0python\python.exe" -m sangerscope app
if errorlevel 1 (
    echo.
    echo SangerScope stopped with an error. The message above says why.
    pause
)
