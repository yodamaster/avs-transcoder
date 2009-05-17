tavs.exe ..\..\..\..\..\..\sequences\test.wmv test.avs -l 0 -t 2 -o 0 -b 2
AudioTransCoder.exe ..\..\..\..\..\..\sequences\test.wmv test.mp3
AvsCreator.exe test.avs test.asm -f
AvsCreator.exe test.mp3 test.asm
pause
