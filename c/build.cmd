call "C:\Program Files\Microsoft Visual Studio 8\VC\vcvarsall.bat" x86
rem "C:\Program Files\Microsoft Visual Studio 8\Common7\IDE\VCExpress.exe"  /build "Debug|Win32" solver.vcproj /Project solver
cl.exe /Od /I "C:\myprojects\expat2.0.1\Source\lib\\" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "XML_STATIC" /D "USE_EXPAT" /D "_CRT_SECURE_NO_WARNINGS" /D "_UNICODE" /D "UNICODE" /Gm /RTC1 /MT /Za /Fo"Debug\\" /Fd"Debug\vc80.pdb" /W3 /c /Wp64 /ZI /TC ".\fea_solver.c"
link.exe Debug/fea_solver.obj /OUT:"C:\myprojects\finite-strains\fea-solver\c\Debug\solver.exe" /INCREMENTAL /LIBPATH:"C:\Program Files\Microsoft SDKs\Windows\v6.1\Lib\;C:\myprojects\expat2.0.1\Source\win32\bin\debug\" /MANIFEST /MANIFESTFILE:"Debug\solver.exe.intermediate.manifest" /DEBUG /PDB:"c:\myprojects\finite-strains\fea-solver\c\Debug\solver.pdb" /SUBSYSTEM:CONSOLE /MACHINE:X86 libexpatMT.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib 

