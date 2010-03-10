// AudioTransCoder.cpp : Defines the entry point for the console application.
//

#include <windows.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

/*
#include <libtsp.h>
#include <libtsp/AFpar.h>
*/

#include <atlbase.h>
#include <streams.h>
#include <qedit.h>         // for Null Renderer

#include <stdio.h>
#include <math.h>

#include "grabber.h"
#include "BladeMP3EncDLL.h"
#include "mpeg4ip_getopt.h"

#include "wavfmt.h"

#include "sam_encode.h"
#include "av3enc.h"
#include "getcmdarg.h"    


double max_e=0.,min_e=100.;


BEINITSTREAM		beInitStream=NULL;
BEENCODECHUNK		beEncodeChunk=NULL;
BEDEINITSTREAM		beDeinitStream=NULL;
BECLOSESTREAM		beCloseStream=NULL;
BEVERSION			beVersion=NULL;
BEWRITEVBRHEADER	beWriteVBRHeader=NULL;
BEWRITEINFOTAG		beWriteInfoTag=NULL;

#define BUFSIZE_StrInfo (80) 

HINSTANCE	hDLL			=NULL;
PBYTE		pMP3Buffer		=NULL;
PBYTE		pWAVBuffer		=NULL;

BE_CONFIG	beConfig		={0,};

DWORD		dwSamples		=0;
DWORD		dwMP3Buffer		=0;
HBE_STREAM	hbeStream		=0;
BE_ERR		err				=0;
DWORD		dwWrite=0;
int			bit_rate;
int			BUFFERLENGTH;
int			buffer_start, buffer_end;
CHAR		pFileInput[500];		//输入输出文件
CHAR		pFileOutputAudio[500];
FILE*		inputFile;
FILE*		m2aFile;
FILE*		fpErrorInfo;
double      length;
bool		flag = true;
int			SampleRate;
///////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
int			Vbit_rate;
int			ReSampleRate;
double		multiple;
int			mode;
int			Orimode;

	char* usageString = 
		" <options> <file>\n"
		"  Options:\n"
		"  -bitrate=<bitrate>			CBR bitrate, VBR min bitrate\n"
		"  -Vbitrate=<Vbitrate>			VBR Max bitrate\n"
		"  -samplerate=<samplerate>		ReSampleRate\n"
		"  -mode=<mode>					STEREO-0,JSTEREO-1,DUALCHANNEL-2,MONO-3\n"
		"  -help						Display the Usage\n"
		;

////////////////////////////////////////////////////////////////////////////////////////////////////////*///

/************************************************************************/
/*  Added By YueLei Xie                 
    Global Values for AV3Encoder
*/
/************************************************************************/

typedef enum
{
	OUT_MP3 = 0,
	OUT_AV3
}OUT_FILE_TYPE;

    OUT_FILE_TYPE out_file_type = OUT_AV3;
	int frames, currentFrame,FrameCnt,samplesRead = 0,cutOff = -1,bitRate = 0;
	unsigned long samplesInput, maxBytesOutput, totalBytesWritten=0;
	unsigned int AVSVer = AVS1,useSquarePolar = 1;
	AV3EncFramePtr hEncoder;
	AV3EncCfgPtr configPtr;
	float *pcmbuf;
	unsigned char *bytebuf;
	unsigned char *bitbuf;
	int shortctl = SHORTCTL_NORMAL;
	//	enum stream_format stream;
	char *audioFileName = NULL;
	char *av3FileName = NULL;
	char *av3FileExt = NULL;
	FILE *outfile = NULL;
	pcmfile_t *infile = NULL;

	// wlei [20060223]
	unsigned char chRawStreamLength[4];

	// xun 2006-4-16 buffer control
	int  iBufferSize;
	float fltStartRatio;

	int         narg;
	int         fn;
	cmd_params  param;
	cmd_option* option;

	char* format_string;
  
	/* functios for AV3Enc */
	int av3enc_init(WAVEFORMATEX * pWavFmt);
	int av3enc_close();

	extern void sam_cbc_init(int fsidx, int bitrate, FILE *outFile);

/************************************************************************/
/*  Ended By YueLei Xie                                                 */
/************************************************************************/

HRESULT GetPin(IBaseFilter * pFilter, PIN_DIRECTION dirrequired,  int iNum, IPin **ppPin);
IPin *  GetInPin ( IBaseFilter *pFilter, int Num );
IPin *  GetOutPin( IBaseFilter *pFilter, int Num );
int GetVideoInfo(char* FileName, double* fps, double* length);

/* Added By YueLei Xie */
FILE * fpWavFile;

RIFF_HEADER chunk_riff;
FMT_BLOCK chunk_wave;
FACT_BLOCK  chunk_fact;
DATA_BLOCK  chunk_data;

size_t  tot_write;

/* Ended By YueLei Xie */


HRESULT Callback(IMediaSample* pSample, REFERENCE_TIME* StartTime,
  REFERENCE_TIME* StopTime, BOOL TypeChanged)
{
	int s;
	int iSampleSize;
	int bytesWritten;
	/*
	if (flag == true && *StartTime > 0)   //在音频前面插入0  用于TS流时应去掉
	{
		iSampleSize = SampleRate * (*StartTime) / 10000000 * 4;
		int sample_read = 0;
		
		while(iSampleSize>=(BUFFERLENGTH-buffer_start))
		{
			memset(&pWAVBuffer[buffer_start], 0, BUFFERLENGTH-buffer_start);
			sample_read += BUFFERLENGTH-buffer_start;
			iSampleSize -= BUFFERLENGTH-buffer_start;
			buffer_start = 0;
			
			for ( s = 0; s < BUFFERLENGTH; s += 2 )
			{
				pcmbuf[s / 2] = *((short *)(pWAVBuffer + s));
			}
			bytesWritten = AV3EncEncode(hEncoder,
				(int *)pcmbuf,
				BUFFERLENGTH / 2,
				bitbuf);
		}
		
		if(iSampleSize!=0)
		{
			memset(&pWAVBuffer[buffer_start], 0, iSampleSize);
			buffer_start += iSampleSize;
		}	
		flag = false;
	}	//在音频前面插入0  用于TS流时应去掉
	*/
	iSampleSize = pSample->GetActualDataLength();

	unsigned char *pBuffer;
	pSample->GetPointer(&pBuffer);

	if (pBuffer != NULL)
	{
		int sample_read = 0;
		
		while(iSampleSize>=(BUFFERLENGTH-buffer_start))
		{
///////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
			int i, NumOfSamples = BUFFERLENGTH-buffer_start;
			short *temp =(short *) pWAVBuffer;
			if( mode == 3 && Orimode != 1)
			{
				for( i = 0; i < NumOfSamples; i += 2)
				{
					pWAVBuffer[buffer_start + i] = (unsigned char) *(pBuffer + sample_read +  i*2 );
					pWAVBuffer[buffer_start + i+1] = (unsigned char) *(pBuffer + sample_read + i*2 + 1);
				}
				for( i = 0; i < NumOfSamples/2; i += 2)
					temp[ i/2 ]=(short)( (double)temp[i/2] * multiple );
			}
			else
			{
				memcpy(&pWAVBuffer[buffer_start], (unsigned char *)pBuffer + sample_read, BUFFERLENGTH-buffer_start);
				for( i = 0; i < NumOfSamples; i += 2)
					temp[ i/2 ]=(short)( (double)temp[i/2] * multiple );
			}
/////////////////////////////////////////////////////////////////////////////////////////////////////////*///
			sample_read += BUFFERLENGTH-buffer_start;
			iSampleSize -= BUFFERLENGTH-buffer_start;
			buffer_start = 0;

			if ( out_file_type == OUT_MP3 )
			{
				///////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
				if( mode == 3 && Orimode != 1)
					err = beEncodeChunk(hbeStream, BUFFERLENGTH/4, (short *) pWAVBuffer, pMP3Buffer, &dwWrite);
				else
					err = beEncodeChunk(hbeStream, BUFFERLENGTH/2, (short *) pWAVBuffer, pMP3Buffer, &dwWrite);
				/////////////////////////////////////////////////////////////////////////////////////////////////////////*///

				if(err != BE_ERR_SUCCESSFUL)
				{
					beCloseStream(hbeStream);
					fprintf(stdout,"beEncodeChunk() failed (%lu)", err);
					return -1;
				}

				fwrite(pMP3Buffer,1,dwWrite, m2aFile);
			}
			else
			{
				/* Added By YueLei */
				for ( s = 0;  s < BUFFERLENGTH; s += 2 )
				{
					pcmbuf[s / 2] = *((short *) (pWAVBuffer + s));
				}

				bytesWritten = AV3EncEncode(hEncoder,(int * ) pcmbuf, BUFFERLENGTH / 2, bitbuf);
			}
		}
		
		if(iSampleSize!=0)
		{
////////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
			int i, NumOfSamples = iSampleSize;
			short *temp =(short *) pWAVBuffer;
			if( mode == 3 && Orimode != 1)
			{
				for( i = 0; i < NumOfSamples; i += 2)
				{
					pWAVBuffer[buffer_start + i] = (unsigned char) *(pBuffer + sample_read +  i*2 );
					pWAVBuffer[buffer_start + i+1] = (unsigned char) *(pBuffer + sample_read + i*2 + 1);
				}
				for( i = 0; i < NumOfSamples/2; i += 2)
					temp[ i/2 ]=(short)( (double)temp[i/2] * multiple );
			}
			else
			{
				memcpy(&pWAVBuffer[buffer_start], (unsigned char *)pBuffer + sample_read, iSampleSize);
				for( i = 0; i < NumOfSamples; i += 2)
					temp[ i/2 ]=(short)( (double)temp[i/2] * multiple );
			}
/////////////////////////////////////////////////////////////////////////////////////////////////////////*///
			buffer_start += iSampleSize;
		}	
	}

	int iPercent;
	iPercent = (*StartTime) / (length * 100000);
	/*
	char StrInfoBuffer[BUFSIZE_StrInfo];

	memset(StrInfoBuffer, 0, BUFSIZE_StrInfo);
	sprintf(StrInfoBuffer, "\rProgress:\t  %d%%", iPercent);
	
	WriteFile( stdout, StrInfoBuffer, BUFSIZE_StrInfo, &dwWrite, 0);
	FlushFileBuffers(stdout);
	*/
	fprintf(stdout, "\rProgress:\t%d%%", iPercent);
	fflush(stdout);
	
	return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
void Initial(int argc, char* argv[])
{

	while (true) 
	{
		int c = -1;
		int option_index = 0;
		static struct option long_options[] = {
			{ "bitrate", 1, 0, 'b' },
			{ "help", 0, 0, '?' },
			{ "samplerate", 1, 0, 's' },
			{ "mode", 1, 0, 'm' },
			{ "Vbitrate", 1, 0, 'v' },
			{ NULL, 0, 0, 0 }
		};
		
		c = getopt_long_only(argc, argv, "b:s:m:v:?",
			long_options, &option_index);												
		
		if (c == -1)
			break;
		
		switch (c) {
		case 'b':
			if (optarg == NULL) {
				fprintf(stderr, "No bitrate specified \n");
				exit(EXIT_FAILURE);
			}
			if (sscanf(optarg, "%u", &bit_rate) != 1) {
				fprintf(stderr,	"Bad bitrate specified: %s\n",	optarg);
				exit(EXIT_FAILURE);
			}
			break;

		case 's':
			if (optarg == NULL) {
				fprintf(stderr, "No samplerate specified \n");
				exit(EXIT_FAILURE);
			}
			if (sscanf(optarg, "%u", &ReSampleRate) != 1) {
				fprintf(stderr,	"Bad samplerate specified: %s\n",	optarg);
				exit(EXIT_FAILURE);
			}
			break;
			
		case 'm':
			if (optarg == NULL) {
				fprintf(stderr, "No mode specified \n");
				exit(EXIT_FAILURE);
			}
			if (sscanf(optarg, "%u", &mode) != 1) {
				fprintf(stderr,	"Bad mode specified: %s\n",	optarg);
				exit(EXIT_FAILURE);
			}
			break;
			
		case 'v':
			if (optarg == NULL) {
				fprintf(stderr, "No VBR bitrate specified \n");
				exit(EXIT_FAILURE);
			}
			if (sscanf(optarg, "%u", &Vbit_rate) != 1) {
				fprintf(stderr,	"Bad VBR bitrate specified: %s\n",	optarg);
				exit(EXIT_FAILURE);
			}
			break;
			
		case '?':
			fprintf(stderr, "usage: %s", usageString);
			exit(EXIT_SUCCESS);

		default:
			fprintf(stderr, "Unknown option specified, ignoring: %c\n", c);
    }
  }
  // check that we have at least one non-option argument
  if ((argc - optind) < 1) {
	  fprintf(stderr, "usage: %s", usageString);                         
	  exit(EXIT_FAILURE);
  }

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////*///


void init_wav_file(WAVEFORMATEX * pWF)
{
	
	memset(&chunk_riff, 0, sizeof(chunk_riff));
	memset(&chunk_wave, 0, sizeof(chunk_wave));
	memset(&chunk_data, 0, sizeof(chunk_data));

	/* RIFF CHUNK */
	strcpy(chunk_riff.strRiffID, "RIFF");
	chunk_riff.dwSize = 0;
	strcpy(chunk_riff.strWavFmt, "WAVE");

	/* FOMAT CHUNK */
	strcpy(chunk_wave.szFmtID, "fmt ");
	chunk_wave.dwFmtSize = pWF->cbSize;
	chunk_wave.dwFmtSize = 16;
	chunk_wave.wavFormat.wFormatTag       = pWF->wFormatTag;
	chunk_wave.wavFormat.wChannels        = pWF->nChannels;
	chunk_wave.wavFormat.dwSamplesPerSec  = pWF->nSamplesPerSec;
	chunk_wave.wavFormat.dwAvgBytesPerSec = pWF->nAvgBytesPerSec;
	chunk_wave.wavFormat.wBitsPerSample   = pWF->wBitsPerSample;
	chunk_wave.wavFormat.wBlockAlign      = pWF->nBlockAlign;

	/* DATA CHUNK */
	strcpy(chunk_data.szDataID, "data");

	
	fwrite((char *)&chunk_riff, 1, sizeof(chunk_riff), fpWavFile);
	fwrite((char *)&chunk_wave, 1, sizeof(chunk_wave), fpWavFile);
	fwrite((char *)&chunk_data, 1, sizeof(chunk_data), fpWavFile);
	
}

int main(int argc, char* argv[])
{
    CoInitialize( NULL );
	//Get input params
    if (argc <= 1)
	{
		fprintf(stderr, "usage: %s", usageString);                         
		exit(EXIT_FAILURE);
	}
	
	memcpy(pFileInput, argv[1], 500);
////////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
	memcpy(pFileOutputAudio, argv[2], 500);
//	strcat(pFileOutputAudio, ".mp3");

	out_file_type = OUT_AV3;

	char * pStr = strstr(pFileOutputAudio, ".mp3");

	if ( pStr && strlen(pStr) == 4 )
	{
		out_file_type = OUT_MP3;
	}

	bit_rate = 64;
	Vbit_rate = 0;///Lirh 2008.10
	ReSampleRate = 0;
	multiple = 1;
	mode = -1;///Lirh 2008.10
	Initial(argc, argv);///Lirh 2008.10
///////////////////////////////////////////////////////////////////////////////////////////////////////////*/

	
	double fps;
	GetVideoInfo(pFileInput, &fps, &length);	//Get fps, length of the video from file

	buffer_start = 0;

    // The sample grabber is not in the registry, so create it with 'new'.
    HRESULT hr = S_OK;
	CComPtr< IGrabberSample > pGrab;
    
    pGrab.CoCreateInstance( CLSID_GrabberSample );
    //pGrab->AddRef();

    // Set the callback function of the filter.
    pGrab->SetCallback(&Callback);

    // Set up a partially specified media type.
    CMediaType mt;

 	mt.SetType( &MEDIATYPE_Audio );
 	mt.SetSubtype( &MEDIASUBTYPE_PCM );

    hr = pGrab->SetAcceptedMediaType(&mt);

    // Create the filter graph manager.
    CComPtr<IFilterGraph> pGraph;
    hr = pGraph.CoCreateInstance( CLSID_FilterGraph );

    // Query for other useful interfaces.
    CComQIPtr<IGraphBuilder, &IID_IGraphBuilder> pBuilder(pGraph);
    CComQIPtr<IMediaSeeking, &IID_IMediaSeeking> pSeeking(pGraph);
    CComQIPtr<IMediaControl, &IID_IMediaControl> pControl(pGraph);
    CComQIPtr<IMediaFilter, &IID_IMediaFilter> pMediaFilter(pGraph);
    CComQIPtr<IMediaEvent, &IID_IMediaEvent> pEvent(pGraph);

    // Add a source filter to the graph.
    CComPtr<IBaseFilter> pSource;
    hr = pBuilder->AddSourceFilter( A2WBSTR(pFileInput), L"Source", &pSource);

	if( FAILED( hr ) )
    {
        _tprintf( TEXT("Could not find the input file\r\n") );
        return -1;
    }

    // Add the sample grabber to the graph.
	CComQIPtr< IBaseFilter, &IID_IBaseFilter > pGrabberBase( pGrab );
	
    hr = pBuilder->AddFilter(pGrabberBase, L"Grabber");

	//Enumerating Pins of Source filter, and find number of pins
	IEnumPins  *pEnum;
    IPin       *pPin;
	int			i, num = 0;

    hr = pSource->EnumPins(&pEnum);
    if (FAILED(hr))
    {
        _tprintf( TEXT("Could not Enumerate Pins of Source filter\r\n") );
		return -1;
    }
    while(pEnum->Next(1, &pPin, 0) == S_OK)
    {
		num++;
		pPin->Release();
	}
	pEnum->Release();

    // Find the input and output pins, and connect them.
    IPin *pGrabIn = GetInPin(pGrabberBase, 0);
	for (i=0; i<num; i++)
	{
		IPin *pSourceOut = GetOutPin(pSource, i);
		hr = pBuilder->Connect(pSourceOut, pGrabIn);
		if (hr == S_OK)		//some times hr = VFW_S_PARTIAL_RENDER ?
		{
			break;
		}
	}
	if( FAILED( hr ) )
	{
		_tprintf( TEXT("Could not connect source filter to grabber\r\n") );
		fpErrorInfo = fopen("AudioErrorInfo.txt", "a");
		SYSTEMTIME time;
		TCHAR szDate[64], szTime[64];
		GetLocalTime(&time);
		GetDateFormat (LOCALE_USER_DEFAULT, LOCALE_NOUSEROVERRIDE | DATE_SHORTDATE,
                    &time, NULL, szDate, sizeof (szDate)) ;
		GetTimeFormat ( LOCALE_USER_DEFAULT, LOCALE_NOUSEROVERRIDE | 
                        TIME_NOTIMEMARKER | TIME_FORCE24HOURFORMAT,
                    &time, NULL, szTime, sizeof (szTime)) ;

		fwrite("*****", 5, 1, fpErrorInfo);
		fwrite(szDate, strlen(szDate), 1, fpErrorInfo);
		fwrite(" ", 1, 1, fpErrorInfo);
		fwrite(szTime, strlen(szTime), 1, fpErrorInfo);
		fwrite("*****\r\n", 7, 1, fpErrorInfo);
		fwrite("Could not connect source filter to grabber\r\n", 44, 1, fpErrorInfo);
		fwrite("InputFile: ", 12, 1, fpErrorInfo);
		fwrite(pFileInput, strlen(pFileInput), 1, fpErrorInfo);
		fwrite("\r\n\r\n", 4, 1, fpErrorInfo);
		fclose(fpErrorInfo);
		return -1;
	}
    
    // Create the Null Renderer filter and add it to the graph.
    CComPtr<IBaseFilter> pNull;
    hr = pNull.CoCreateInstance(CLSID_NullRenderer);
    hr = pBuilder->AddFilter(pNull, L"Renderer");

    // Get the other input and output pins, and connect them.
    IPin *pGrabOut = GetOutPin(pGrabberBase, 0);
    IPin *pNullIn = GetInPin(pNull, 0);
    hr = pBuilder->Connect(pGrabOut, pNullIn);

    CMediaType  mt1;

    hr = pGrab->GetConnectedMediaType( &mt1 );

//FORMAT_WAVEFORMATEX
	
	WAVEFORMATEX * pWF = (WAVEFORMATEX*) mt1.pbFormat;

	if ( out_file_type == OUT_MP3 )
	{
		hDLL = LoadLibrary("lame_enc.dll");
		if ( NULL == hDLL )
		{
			fprintf(stdout,"Error loading lame_enc.DLL");
			return -1;
		}

		// Get Interface functions from the DLL
		beInitStream	= (BEINITSTREAM) GetProcAddress(hDLL, TEXT_BEINITSTREAM);
		beEncodeChunk	= (BEENCODECHUNK) GetProcAddress(hDLL, TEXT_BEENCODECHUNK);
		beDeinitStream	= (BEDEINITSTREAM) GetProcAddress(hDLL, TEXT_BEDEINITSTREAM);
		beCloseStream	= (BECLOSESTREAM) GetProcAddress(hDLL, TEXT_BECLOSESTREAM);
		beVersion		= (BEVERSION) GetProcAddress(hDLL, TEXT_BEVERSION);
		beWriteVBRHeader= (BEWRITEVBRHEADER) GetProcAddress(hDLL,TEXT_BEWRITEVBRHEADER);
		beWriteInfoTag  = (BEWRITEINFOTAG) GetProcAddress(hDLL,TEXT_BEWRITEINFOTAG);

		// Check if all interfaces are present
		if (!beInitStream || !beEncodeChunk || !beDeinitStream || !beCloseStream || !beVersion || !beWriteVBRHeader)
		{
			printf("Unable to get LAME interfaces");
			return -1;
		}

		memset(&beConfig,0,sizeof(beConfig));					// clear all fields

		// use the LAME config structure
		beConfig.dwConfig = BE_CONFIG_LAME;

		// this are the default settings for testcase.wav
		beConfig.format.LHV1.dwStructVersion	= 1;
		beConfig.format.LHV1.dwStructSize		= sizeof(beConfig);		
		beConfig.format.LHV1.dwSampleRate		= pWF->nSamplesPerSec;		//Sample rate
		beConfig.format.LHV1.dwReSampleRate		= pWF->nSamplesPerSec;					// DON"T RESAMPLE
		SampleRate = pWF->nSamplesPerSec;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
		if(ReSampleRate != 0)
		{
			beConfig.format.LHV1.dwReSampleRate	= ReSampleRate;					
			SampleRate = ReSampleRate;
			multiple = sqrt( (float)pWF->nSamplesPerSec / ReSampleRate );
		}
		if(mode == -1)
		{
			if (pWF->nChannels == 1)									//Channels
			{
				beConfig.format.LHV1.nMode			= BE_MP3_MODE_MONO;	// OUTPUT IN MONO
			}
			else
			{
				beConfig.format.LHV1.nMode			= BE_MP3_MODE_JSTEREO;	// OUTPUT IN STREO
			}
		}
		else
			beConfig.format.LHV1.nMode			= mode;
		Orimode = pWF->nChannels;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		beConfig.format.LHV1.dwBitrate			= bit_rate;					// MINIMUM BIT RATE
		beConfig.format.LHV1.nPreset			= LQP_HIFI;		// QUALITY PRESET SETTING
		beConfig.format.LHV1.dwMpegVersion		= MPEG2;				// MPEG VERSION (I or II)
		beConfig.format.LHV1.dwPsyModel			= 0;					// USE DEFAULT PSYCHOACOUSTIC MODEL 
		beConfig.format.LHV1.dwEmphasis			= 0;					// NO EMPHASIS TURNED ON
		beConfig.format.LHV1.bOriginal			= TRUE;					// SET ORIGINAL FLAG
		beConfig.format.LHV1.bWriteVBRHeader	= TRUE;					// Write INFO tag

		////////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10

		if(Vbit_rate >= bit_rate)
		{
			beConfig.format.LHV1.dwMaxBitrate		= Vbit_rate;					// MAXIMUM BIT RATE
			beConfig.format.LHV1.bWriteVBRHeader	= TRUE;					// YES, WRITE THE XING VBR HEADER
			beConfig.format.LHV1.bEnableVBR			= TRUE;					// USE VBR
			beConfig.format.LHV1.nVbrMethod			= VBR_METHOD_DEFAULT;	// VBR METHOD
			beConfig.format.LHV1.nVBRQuality		= 6;					// SET VBR QUALITY
		}
		else
		{
			beConfig.format.LHV1.dwMaxBitrate		= bit_rate;					// MAXIMUM BIT RATE
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//	beConfig.format.LHV1.bCRC				= TRUE;					// INSERT CRC
		//	beConfig.format.LHV1.bCopyright			= TRUE;					// SET COPYRIGHT FLAG	
		//	beConfig.format.LHV1.bPrivate			= TRUE;					// SET PRIVATE FLAG
		beConfig.format.LHV1.bNoRes				= false;					// No Bit resorvoir

		// Init the MP3 Stream
		err = beInitStream(&beConfig, &dwSamples, &dwMP3Buffer, &hbeStream);

		// Check result
		if (err != BE_ERR_SUCCESSFUL)
		{
			fprintf(stdout,"Error opening encoding stream (%lu)", err);
			return -1;
		}


		// Allocate MP3 buffer
		pMP3Buffer = new BYTE[dwMP3Buffer];

		// Allocate WAV buffer
		if (pWF->nChannels == 1)
		{
			pWAVBuffer = new BYTE[dwSamples];

			BUFFERLENGTH = dwSamples;

		}
		else	
		{
			pWAVBuffer = new BYTE[dwSamples*2];

			BUFFERLENGTH = dwSamples*2;
		}

		m2aFile = fopen(pFileOutputAudio, "wb");
	}
	else
	{
		SampleRate = pWF->nSamplesPerSec;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
		if(ReSampleRate != 0)
		{
			beConfig.format.LHV1.dwReSampleRate	= ReSampleRate;					
			SampleRate = ReSampleRate;
			multiple = sqrt( (float)pWF->nSamplesPerSec / ReSampleRate );
		}

		Orimode = pWF->nChannels;


		/*set default value*/
		param.bitrate = bit_rate;
		param.cutoff = -1;
		param.stereo = 1;
		param.wincontrol = 0;
		param.outFile = NULL;
		param.inFile = NULL;
		param.showHelp = 0;
		param.format = 0;
		param.LFE_3SFB = 1;
		param.cbr	   = 1;

		m2aFile = fopen(pFileOutputAudio, "wb");	//Open output file

		av3enc_init(pWF);

		BUFFERLENGTH = 1024 * pWF->nBlockAlign;

		// Allocate WAV buffer
		pWAVBuffer = new BYTE[BUFFERLENGTH];
		
		/* Added By YueLei Xie */
		  
		/* */
	}

	FreeMediaType( mt1 );
	
	


	fprintf(stdout, "InputFile: %s\r\n", pFileInput);
	fprintf(stdout, "OutputFile: %s\r\n", pFileOutputAudio);
	fprintf(stdout, "BitRate: %d\r\n", bit_rate);
	fprintf(stdout, "Encoding, please wait...\r\n"); 
	fflush(stdout);

	hr = pMediaFilter->SetSyncSource(NULL);

	if ( FAILED( hr ) )
    {
        _tprintf( TEXT("Could not SetSyncSource pMediaFilter\r\n") );
        return -1;
    }

    hr = pControl->Run();
	
	long EvCode = 0;

	hr = pEvent->WaitForCompletion( INFINITE, &EvCode );
	

	if ( out_file_type == OUT_MP3 )
	{
		if (dwWrite)
		{
			fwrite( pMP3Buffer, 1, dwWrite, m2aFile);
		}

		beCloseStream( hbeStream );

		// Delete WAV buffer
		delete [] pWAVBuffer;
	}
	else
	{
		av3enc_close();

		delete [] pWAVBuffer;
	}
	
	fprintf(stdout, "\rProgress:\t100%%");
	fflush(stdout);
	printf("\nMission completed!\r\n");
	
//  pGrab.Release();
//	pGrab.~CComPtr;
    CoUninitialize();
	pGrab.p->Release();
    //pSeeking.p->Release();
	pControl.p->Release();
	pMediaFilter.p->Release();
	pEvent.p->Release();
	pSource.p->Release();
	pGrabberBase.p->Release();
	pNull.p->Release();
	pBuilder.p->Release();
	//pGraph.p->Release();
	exit(0);
    return 0;
}

HANDLE hf = NULL;


HRESULT GetPin( IBaseFilter * pFilter, PIN_DIRECTION dirrequired, int iNum, IPin **ppPin)
{
    CComPtr< IEnumPins > pEnum;
    *ppPin = NULL;

    HRESULT hr = pFilter->EnumPins(&pEnum);
    if(FAILED(hr)) 
        return hr;

    ULONG ulFound;
    IPin *pPin;
    hr = E_FAIL;

    while(S_OK == pEnum->Next(1, &pPin, &ulFound))
    {
        PIN_DIRECTION pindir = (PIN_DIRECTION)3;

        pPin->QueryDirection(&pindir);
        if(pindir == dirrequired)
        {
            if(iNum == 0)
            {
                *ppPin = pPin;  // Return the pin's interface
                hr = S_OK;      // Found requested pin, so clear error
                break;
            }
            iNum--;
        } 

        pPin->Release();
    } 

    return hr;
}


IPin * GetInPin( IBaseFilter * pFilter, int nPin )
{
    CComPtr<IPin> pComPin=0;
    GetPin(pFilter, PINDIR_INPUT, nPin, &pComPin);
    return pComPin;
}


IPin * GetOutPin( IBaseFilter * pFilter, int nPin )
{
    CComPtr<IPin> pComPin=0;
    GetPin(pFilter, PINDIR_OUTPUT, nPin, &pComPin);
    return pComPin;
}

int GetVideoInfo(char* FileName, double* fps, double* length)
{
	//////////////////////////////////////////////////////////////////////////
	//获取源文件帧率等信息 hhan 2006-08-20
	//**Obtain the authored frame rate of a video file without rendering it
	
	HRESULT res;
	IMediaDet *pMediaDet;
	
	res = CoCreateInstance(CLSID_MediaDet, NULL, CLSCTX_INPROC, IID_IMediaDet,
		(void **)&pMediaDet);
	if(FAILED(res))
	{
		printf("CoCreateInstance(): res=0x%08x\n", res);
		return -1;
	}
	
//	OLECHAR bstrFilename[MAX_PATH];
//	MultiByteToWideChar(CP_ACP, 0, argv[1], -1, bstrFilename, sizeof(bstrFilename) / sizeof(OLECHAR));
//	memcpy(bstrFilename, argv[1], MAX_PATH);
	res = pMediaDet->put_Filename(A2WBSTR(FileName));
	if(FAILED(res))
	{
		printf("put_Filename(): res=0x%08x\n", res);
		return -1;
	}
	
	long lStreams, lStreamSelected;
	
	res = pMediaDet->get_OutputStreams(&lStreams);
	if(FAILED(res))
	{
		printf("get_OutputStreams(): res=0x%08x\n", res);
		return -1;
	}
	
	AM_MEDIA_TYPE amMediaType;
	BITMAPINFOHEADER *pbih;
	
	lStreamSelected = -1;
	for(long lStream = 0; lStream < lStreams; lStream++)
	{
		res = pMediaDet->put_CurrentStream(lStream);
		if(FAILED(res))
		{
			printf("put_CurrentStream(): res=0x%08x\n", res);
			return -1;
		}
		
		res = pMediaDet->get_StreamMediaType(&amMediaType);
		if(FAILED(res))
		{
			printf("get_StreamMediaType(): res=0x%08x\n", res);
			return -1;
		}
		
		char cSelected = ' ';
		if(IsEqualGUID(amMediaType.majortype, MEDIATYPE_Video))
		{
			lStreamSelected = lStream;
			cSelected = '*';
		}
		
		OLECHAR szGUID[256];
		
		StringFromGUID2(amMediaType.majortype, szGUID, sizeof(szGUID));
		printf("[%2d] majortype: %S %c\n", lStream,
			IsEqualGUID(amMediaType.majortype, MEDIATYPE_Video) ? L"MEDIATYPE_Video" :
		IsEqualGUID(amMediaType.majortype, MEDIATYPE_Audio) ? L"MEDIATYPE_Audio" :
		szGUID, cSelected);
		
		StringFromGUID2(amMediaType.formattype, szGUID, sizeof(szGUID));
		printf("     formattype: %S\n",
			IsEqualGUID(amMediaType.formattype, FORMAT_VideoInfo) ? L"FORMAT_VideoInfo" :
		IsEqualGUID(amMediaType.formattype, FORMAT_WaveFormatEx) ? L"FORMAT_WaveFormatEx" :
		szGUID);
		
		if(IsEqualGUID(amMediaType.formattype, FORMAT_VideoInfo))
		{
			pbih = &((VIDEOINFOHEADER *)amMediaType.pbFormat)->bmiHeader;
			printf("     (%dx%dx%d '%4.4s')\n", pbih->biWidth, pbih->biHeight,
				pbih->biBitCount, &pbih->biCompression);
		}
	}
	
	if(lStreamSelected == -1)
	{
		printf("No video stream found!\n");
		return -1;
	}
	
	res = pMediaDet->put_CurrentStream(lStreamSelected);
	if(FAILED(res))
	{
		printf("put_CurrentStream(): res=0x%08x\n", res);
		return -1;
	}
	
	if(FAILED(pMediaDet->get_FrameRate(fps)) || FAILED(pMediaDet->get_StreamLength(length)))
	{
		printf("get_FrameRate/StreamLength(): res=0x%08x\n", res);
		return -1;
	}
	
	UINT nFrames = (UINT)(*length * *fps);
	printf("Length of Video: %gs\r\nAuthored Frame Rate: %gfps\r\n(%d frames of %dx%dx%d)\n", *length, *fps, nFrames,
		pbih->biWidth, pbih->biHeight, pbih->biBitCount);
	
	//End of 获取源文件帧率等信息 hhan 2006-08-20
	//**Obtain the authored frame rate of a video file without rendering it
	
	//////////////////////////////////////////////////////////////////////////
	return 0;
}

static const char *logo =
"\n  AVS Audio encoder version 1.0 \n";
static const char *help =
"\n%s\t[-b[bitrate]][-c[cutoff]][-s<+,->][-lfe<0,1>][-o[outfilename]]\n"
"    \t[-t<0,1,2>][-f<0,1,2>][-h/-help] [infilename]\n\n"
"  -b[bitrate]      :\t set bitrate in kbps\n"
"  -c[cutoff]       :\t set cutoff frequency in Hz\n"
"  -s<+,->          :\t enable('+') or disable('-') stereo coding\n"
"  -o[outfilename]  :\t specify output file name\n"
"  -t<0,1,2>        :\t long/short window control\n"
"                   :\t 0: normal; 1: only long; 2: only short\n"
"  -f<0,1,2>        :\t set output format\n"
"                   :\t 0: raw; 1: storage; 2: transport\n"
"  -h/-help         :\t show this help info\n"
"  -lfe<0,1>        :\t set lfe 3sfb coding\n"
"                   :\t 0:normal coding; 1: 3sfb coding"
"  -vbr<+,->		:\t enable('+') or disable ('-') vbr mode, default enable\n"
"\n  for example: %s -b128 -c16000 -s+ -o outfile.av3 infile.wav\n"
"               %s /b128 /c16000 /s- /o outfile.av3 infile.wav\n"
"\n  IF ANY COMMENTS AND SUGGESTIONS, PLEASE CONTACT AVSaudio@jdl.ac.cn\n";

char* binSet[] = {"+", "-", "\0"};
char* finSet[] = {"0", "1", "2", "\0"};

const cmd_switch  swtArr[] = 
{                                   
	{"b",     1, NULL,    0, 0},  /*bitrate*/ 
	{"c",     1, NULL,    0, 1},  /*cutoff frequency*/
	{"s",     3, binSet,  0, 2},  /*stereo enable or disable*/
	{"o",     2, NULL,    0, 3},  /*output file name*/
	{"t",     3, finSet,  0, 4},  /*long/short window control*/
	{"f",     3, finSet,  0, 5},  /*bitstream format option*/  
	{"h",     0, NULL,    0, 7},  /*show help*/
	{"help",  0, NULL,    0, 6},  /*show help*/
	{"lfe",   1, NULL,    0, 8},  /* lfe 3sfb coding */
	{"vbr",   3, binSet,  0, 9},  /* cbr/vbr mode */ 
	{"\0",    2, NULL,    0, 10}   /*necessary for ending*/
};

static int getParam(cmd_params* param, cmd_option* option, int narg)
{
	int i;
	int f;

	i = 0;
	f = 0;
	while(i < narg){
		char* opt;
		opt = option[i].opt;
		switch(option[i].swIdx)
		{
		case 0:
			if(opt)
				param->bitrate = atof(opt);
			param->bitrate=128;
			break;
		case 1:
			if(opt)
				param->cutoff = atoi(opt);
			break;
		case 2:
			if(opt && !strcmp(opt, "-"))
				param->stereo = 0;
			break;        
		case 3:
			param->outFile = opt;
			break;                        
		case 4:
			if(opt && !strcmp(opt, "1"))
				param->wincontrol = 1;
			else if(opt && !strcmp(opt, "2"))
				param->wincontrol = 2;    
			break;
		case 5:
			if(opt && !strcmp(opt, "1"))
				param->format = 1;
			else if(opt && !strcmp(opt, "2"))
				param->format = 2;    
			break;
		case 6:
		case 7:
			param->showHelp = 1;
			break;
		case 8:
			param->LFE_3SFB = atoi(opt);
			break;
		case 9:			
			if(opt && !strcmp(opt, "-"))
				param->cbr = 1;
			break;     
		case 10:
			param->inFile = opt;
			f++;
			break;
		default:
			break;             
		}
		i++;
	}
	return f;
}

/* channel remapping, especially for LFE */
void ReChannelMap(int numChannels, unsigned char* SpkrConfig, AV3EncFramePtr hEncoder)
{
	int i;
	int isLFE=0;

	/* if LFE */
	i = 0;
	while(SpkrConfig[i] != 4 && i<numChannels)
		i ++;
	if(i<numChannels)
		isLFE++;

	/* ReChannelMapping, not use default config */
	if(isLFE)
	{
		hEncoder->config.isLFE = 1;
		for(i=4; i<numChannels; i++)
			hEncoder->config.channel_map[i-1] = i;
		/* LFE channel */
		hEncoder->config.channel_map[numChannels-1] = 3;
	}
	else {
		for(i=4; i<=numChannels; i++)
			hEncoder->config.channel_map[i-1] = i;
	}
}

int av3enc_init(WAVEFORMATEX * pWavFmt)
{
	useSquarePolar = param.stereo;
	bitRate = param.bitrate*1000;
	cutOff = param.cutoff;
	shortctl = param.wincontrol;

	infile = (pcmfile_t *)malloc(sizeof(*infile));
	memset(infile, 0, sizeof(*infile));

	infile->channels = pWavFmt->nChannels;
	infile->samplerate = pWavFmt->nSamplesPerSec;


	/* open the encoder library */
	hEncoder = AV3EncOpen((int)infile->samplerate, infile->channels,
		&samplesInput, &maxBytesOutput);

	hEncoder->config.LFE_3SFB = param.LFE_3SFB;
	hEncoder->config.isLFE = 0;
	if(infile->channels>6) {
		printf("more than 6 channels, not supported yet\n");
		return 0;
	}
	if(infile->channels>3)
		ReChannelMap(infile->channels, infile->f->SpkrConfig, hEncoder);

	pcmbuf = (float *)malloc(samplesInput*sizeof(float));
	bitbuf = (unsigned char*)malloc(maxBytesOutput*sizeof(unsigned char));
	bytebuf = (unsigned char *)malloc(samplesInput * 4);

	if (cutOff <= 0)
	{
		if (cutOff < 0) 
			cutOff = 0;
		else 
			cutOff = infile->samplerate / 2;
	}
	if (cutOff > (infile->samplerate / 2))
		cutOff = infile->samplerate / 2;

	/* put the options in the configuration struct */
	configPtr = &(hEncoder->config);
	configPtr->AVSVersion = AVSVer;
	configPtr->allowSPSC = useSquarePolar;
	configPtr->bandWidth = cutOff;

	switch (shortctl){
		case SHORTCTL_NOSHORT:
			fprintf(stderr, "disabling short blocks\n");
			configPtr->shortctl = shortctl;
			break;
		case SHORTCTL_NOLONG:
			fprintf(stderr, "disabling long blocks\n");
			configPtr->shortctl = shortctl;
			break;
	}

	if(bitRate<=0)
		bitRate = 128000*infile->channels;
	if(hEncoder->config.isLFE && hEncoder->config.LFE_3SFB)
		configPtr->bitRate = bitRate /(infile->channels-1);
	else
		configPtr->bitRate = bitRate / infile->channels;

	configPtr->outputFormat = param.format;

	// bit reservior buffer setting - xun 2006-4-16
	if (param.cbr)
	{
		iBufferSize = 6144;
		fltStartRatio = 1.0;
	}
	else
	{
		iBufferSize = 614400;	// set a larget value to simulate infinite buffer
		fltStartRatio = 0.5;
	}
	if(hEncoder->config.isLFE)	// LFE bits provided by total bit reservior buffer
		hEncoder->iBufferSize = iBufferSize * (infile->channels - 1);
	else
		hEncoder->iBufferSize = iBufferSize * infile->channels;
	hEncoder->fltBufferFullness = (float)hEncoder->iBufferSize * fltStartRatio;
	hEncoder->fltAvgBitPerFrm = (float)bitRate * (float)BLOCK_LEN_LONG / (float)infile->samplerate;

	if (!AV3EncSetConfiguration(hEncoder, configPtr)) {
		fprintf(stderr, "Unsupported output format!\n");
		return 1;
	}

	/* open the av3 output file */
	//outfile = fopen(av3FileName, "wb");

	outfile = m2aFile;

	cutOff = configPtr->bandWidth;
	bitRate = configPtr->bitRate;

	/* entropy coding initialization */
	/* per channel bitrate */
	sam_cbc_init(hEncoder->config.sampleRateIdx, bitRate, outfile);
	
	/*
	BitRate = bitrate;

	sam_scale_bits_init(hEncoder->config.sampleRateIdx); // zhanjie 0708


	BSHCInit(outfile);
	*/

	/* bitRate -> per channel */
	hEncoder->config.top_layer = (bitRate / 1000) - 16;
	if (hEncoder->config.top_layer > 63) hEncoder->config.top_layer = 63;

		switch(configPtr->outputFormat) {
	case 0:
		format_string = "RAW";
		break;
	case 1:
		format_string = "AASF";		
		break;
	case 2:
		format_string = "AATF";
		break;
	default:
		fprintf(stderr,"unsupported format");
		}

    	/* encoding loop */
		/*
		for ( ; ;)
		{
			int bytesWritten;

			int s;

			samplesRead = fread(bytebuf, 1, samplesInput * 2, raw_infile);

			for (s = 0; s < samplesRead; s += 2 )
			{
				pcmbuf[s / 2] = *((short * )(bytebuf+s));
			}
			//printf("%d ", samplesRead);

			samplesRead /= 2;


			bytesWritten = AV3EncEncode(hEncoder,
				(int *)pcmbuf,
				samplesRead,
				bitbuf);				 

		}
		*/
	return 0;
}



int av3enc_close()
{
	// write AASF header (renew raw_stream_length) wlei [20060223]
	if (configPtr->outputFormat == 1)
	{
		chRawStreamLength[0] = (totalBytesWritten>>24)&0xFF;
		chRawStreamLength[1] = (totalBytesWritten>>16)&0xFF;
		chRawStreamLength[2] = (totalBytesWritten>>8)&0xFF;
		chRawStreamLength[3] = totalBytesWritten&0xFF;
		//	fseek(outfile, 0L, SEEK_SET);
		fseek(outfile, 9L, SEEK_SET);

		fwrite(chRawStreamLength, 1, 4, outfile);
	}

	fclose(outfile);

	AV3EncClose(hEncoder);



	if(option)	free(option);
	if(infile)	free(infile);
	if (pcmbuf) free(pcmbuf);
	if (bitbuf) free(bitbuf);

	return 0;
}