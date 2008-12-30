#ifndef _MPEG2_AUDIO_ENCODER_H_
#define _MPEG2_AUDIO_ENCODER_H_

#ifdef __cplusplus
extern "C" {
#endif 

typedef struct
{
	int sample_rate;    //sampling frequency,support 16/22.5/24/32/44.1/48kHz
	int num_channel;        //number of channels, mono/stereo
	int bit_rate;      //target bit rate, default is 128kbps (32-384kbps)
	/* allowable bitrate for 16,22.5,24kHz are 8,16,24,32,40,48,56,64,80,96,112,128,144,160kbps
	   for 32,44.1,48kHz are 32,48,56,64,80,96,112,128,160,192,224,256,320,384kbps*/

} MPEG2AUDIO_ENCODER_PARAM;

typedef struct
{
	short * insamp;   //pointer to input sample buffer, size is 2304
	unsigned char* Bitstream;   //pointer to output bitstream 
	int Length;     // number of encoded bits per frame

} MPEG2AUDIO_ENCODER_FRAME;

#define DllImport __declspec(dllimport)
#define DllExport __declspec(dllexport)

DllExport int OpenMpeg2AudioEncoder(MPEG2AUDIO_ENCODER_PARAM* Info);
DllExport void CloseMpeg2AudioEncoder(MPEG2AUDIO_ENCODER_FRAME* Frame);
DllExport int EncodeOneMpeg2AudioFrame(MPEG2AUDIO_ENCODER_FRAME* Frame);

#ifdef __cplusplus
}
#endif

#endif // _MPEG2_AUDIO_ENCODER_H_