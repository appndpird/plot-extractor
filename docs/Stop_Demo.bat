@echo off
rem ============================================================
rem  Stops the Steps Recorder session started by
rem  Record_Demo_Auto.bat and saves the .zip to your Desktop.
rem ============================================================
echo Stopping Steps Recorder and saving...
psr.exe /stop
timeout /t 2 >nul
echo.
echo  Saved to:  %USERPROFILE%\Desktop\PlotExtractor_Demo.zip
echo  Open the .zip and view the .mht file to see every step.
pause
