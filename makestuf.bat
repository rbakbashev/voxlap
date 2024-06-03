cls

@call "C:\Program Files (x86)\Microsoft Visual Studio 8\VC\bin\vcvars32.bat"
@set include=%include%;C:\mssdk\include
@set lib=%lib%;C:\mssdk\lib

cl /w /c /J /TP voxlap5.c   /Ox /Ob2 /Gs /MD
REM ml /w /c /coff v5.asm
cl /w /c /J /TP sdlmain.c   /Ox /Ob2 /Gs /MD /DUSEKZ /DZOOM_TEST /DNOSOUND
link voxlap5 v5 sdlmain lib\x86\SDL2main.lib lib\x86\SDL2.lib ddraw.lib dinput.lib ole32.lib dxguid.lib user32.lib gdi32.lib /opt:nowin98

voxlap5.exe
