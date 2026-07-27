@echo off
rem PlotExtractor launcher - uses Python 3.11 which has the Metashape module installed
start "" "%LOCALAPPDATA%\Programs\Python\Python311\pythonw.exe" "%~dp0plot_extractor.py"
