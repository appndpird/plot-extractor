@echo off
rem Build p1tool.exe against a Phase One ImageSDK for Windows (C++).
rem   build_p1tool.bat "<ImageSDK folder with include\ and lib\>"
rem The SDK is licensed by Phase One A/S and is not part of this repository.
rem Needs Visual Studio 2022 (Build Tools) with the C++ workload.
setlocal
if "%~1"=="" goto usage
set "SDK=%~1"
set "VSW=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
set "VSDEV="
for /f "usebackq delims=" %%i in (`call "%VSW%" -latest -products * -property installationPath`) do set "VSDEV=%%i"
if not defined VSDEV goto novs
call "%VSDEV%\VC\Auxiliary\Build\vcvars64.bat" >nul
cd /d "%~dp0"
if not exist obj mkdir obj
cl /nologo /EHsc /std:c++17 /O2 /MD /utf-8 /W3 /wd4251 /D_CRT_SECURE_NO_WARNINGS /Foobj\ /I "%SDK%\include" src\p1tool.cpp /Fe:p1tool.exe /link /LIBPATH:"%SDK%\lib" ImageSdkCpp.lib
exit /b %ERRORLEVEL%
:usage
echo Usage: build_p1tool.bat "ImageSDK folder"
exit /b 2
:novs
echo Visual Studio C++ tools not found
exit /b 2
