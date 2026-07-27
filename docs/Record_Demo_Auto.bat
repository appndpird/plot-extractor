@echo off
rem ============================================================
rem  Plot Extractor - Steps Recorder (hands-free start/stop)
rem  Double-click this to START recording + open the app.
rem  Do the demo, then double-click  Stop_Demo.bat  to save.
rem  Result: %USERPROFILE%\Desktop\PlotExtractor_Demo.zip
rem ============================================================
set "OUT=%USERPROFILE%\Desktop\PlotExtractor_Demo.zip"
echo Starting Steps Recorder (background)...
start "" psr.exe /start /output "%OUT%" /sc 1 /maxsc 300 /gui 0
timeout /t 2 >nul
echo Launching Plot Extractor...
start "" "%~dp0..\PlotExtractor.bat"
echo.
echo  Recording. Perform the demo (see docs\DEMO_SCRIPT.md).
echo  When finished, double-click  Stop_Demo.bat  to save the recording.
echo  It will be saved to: %OUT%
