@echo off
rem ============================================================
rem  Plot Extractor - Steps Recorder demo helper
rem  Opens Steps Recorder AND the app, side by side.
rem  In Steps Recorder:  click  Start Record  ->  do the demo
rem  (follow DEMO_SCRIPT.md)  ->  click  Stop Record  ->  Save.
rem  Output is a .zip containing an .mht slideshow of every click.
rem ============================================================
echo.
echo  Opening Steps Recorder and Plot Extractor...
echo.
echo  1) In Steps Recorder, click  START RECORD
echo  2) Click through the app  (see docs\DEMO_SCRIPT.md)
echo  3) Click  STOP RECORD  in Steps Recorder, then Save the .zip
echo.
start "" psr.exe
timeout /t 1 >nul
start "" "%~dp0..\PlotExtractor.bat"
