@echo off
cd /d "%~dp0"
git push -u origin main
echo.
echo Done! Press any key to close.
pause >nul
