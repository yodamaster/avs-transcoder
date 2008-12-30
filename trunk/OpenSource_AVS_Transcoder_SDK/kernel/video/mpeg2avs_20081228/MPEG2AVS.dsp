# Microsoft Developer Studio Project File - Name="MPEG2AVS" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=MPEG2AVS - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "MPEG2AVS.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "MPEG2AVS.mak" CFG="MPEG2AVS - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "MPEG2AVS - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "MPEG2AVS - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=xicl6.exe
RSC=rc.exe

!IF  "$(CFG)" == "MPEG2AVS - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir ".\bin"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I ".\lencod" /I ".\lencod\inc" /I ".\mpeg2dec" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "TRANSCODING" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=xilink6.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "MPEG2AVS - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ".\bin"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I ".\lencod" /I ".\lencod\inc" /I ".\mpeg2dec" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=xilink6.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "MPEG2AVS - Win32 Release"
# Name "MPEG2AVS - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "AVS ENCODE SRC"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\lencod\src\bitstream.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\block.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\block_const.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\configfile.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\fast_me.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\golomb.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\header.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\image.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\lencod.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\loopfilter.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\macroblock.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\mbuffer.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\memalloc.c
# End Source File
# Begin Source File

SOURCE=".\lencod\src\mv-search.c"
# End Source File
# Begin Source File

SOURCE=.\lencod\src\quadtree.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\ratectl.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\rdopt.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\rdopt_coding_state.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\refbuf.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\slice.c
# End Source File
# Begin Source File

SOURCE=.\lencod\src\vlc.c
# End Source File
# End Group
# Begin Group "MPEG2DEC SRC"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\mpeg2dec\display.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\getbits.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\getblk.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\gethdr.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\getpic.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\getvlc.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\idct.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\idctref.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\motion.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\mpeg2dec.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\recon.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\spatscal.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\store.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\subspic.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\systems.c
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\verify.c
# End Source File
# End Group
# Begin Source File

SOURCE=.\MPEG2AVS.C
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Group "AVS ENCODE HEADER"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\lencod\AVS_ENC_LIB.H
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\bitstream.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\block.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\configfile.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\contributors.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\defines.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\elements.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\fast_me.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\global.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\golomb.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\header.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\image.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\loopfilter.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\macroblock.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\mbuffer.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\memalloc.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\minmax.h
# End Source File
# Begin Source File

SOURCE=".\lencod\inc\mv-search.h"
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\ratectl.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\rdopt_coding_state.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\refbuf.h
# End Source File
# Begin Source File

SOURCE=.\lencod\inc\vlc.h
# End Source File
# End Group
# Begin Group "MPEG2 DEC HEADER"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\mpeg2dec\config.h
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\getvlc.h
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\mpeg2dec.h
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\mpeg2defines.h
# End Source File
# Begin Source File

SOURCE=.\mpeg2dec\mpeg2global.h
# End Source File
# End Group
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\bin\encoder.cfg
# End Source File
# End Target
# End Project
