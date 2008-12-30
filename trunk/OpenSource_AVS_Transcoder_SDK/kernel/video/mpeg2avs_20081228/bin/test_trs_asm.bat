VideoTransCoder.exe test.avi test.avs -l 0 -t 2
AudioTransCoder.exe test.avi test.mp3
AvsCreator.exe test.avs test.asm -f
AvsCreator.exe test.mp3 test.asm
pause
