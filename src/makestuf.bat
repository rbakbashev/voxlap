cls

@call "C:\Program Files (x86)\Microsoft Visual Studio 8\VC\bin\vcvars32.bat"
@set include=%include%;C:\mssdk\include
@set lib=%lib%;C:\mssdk\lib

nmake game.c

game.exe
