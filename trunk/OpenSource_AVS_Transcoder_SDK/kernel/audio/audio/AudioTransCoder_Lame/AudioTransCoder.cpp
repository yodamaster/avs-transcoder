// AudioTransCoder.cpp : Defines the entry point for the console application.
//
#include <atlbase.h>
#include <streams.h>
#include <qedit.h>         // for Null Renderer

#include <stdio.h>
#include <math.h>

#include "grabber.h"
#include "BladeMP3EncDLL.h"
#include "mpeg4ip_getopt.h"

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
CHAR		pFileInput[100];		//输入输出文件
CHAR		pFileOutputAudio[100];
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

HRESULT GetPin(IBaseFilter * pFilter, PIN_DIRECTION dirrequired,  int iNum, IPin **ppPin);
IPin *  GetInPin ( IBaseFilter *pFilter, int Num );
IPin *  GetOutPin( IBaseFilter *pFilter, int Num );
int GetVideoInfo(char* FileName, double* fps, double* length);


HRESULT Callback(IMediaSample* pSample, REFERENCE_TIME* StartTime,
  REFERENCE_TIME* StopTime, BOOL TypeChanged)
{
	int iSampleSize;
/*	if (flag == true && *StartTime > 0)   //在音频前面插入0  用于TS流时应去掉
	{
		iSampleSize = SampleRate * (*StartTime) / 10000000 * 4;
		int sample_read = 0;
		
		while(iSampleSize>=(BUFFERLENGTH-buffer_start))
		{
			memset(&pWAVBuffer[buffer_start], 0, BUFFERLENGTH-buffer_start);
			sample_read += BUFFERLENGTH-buffer_start;
			iSampleSize -= BUFFERLENGTH-buffer_start;
			buffer_start = 0;
			
			err = beEncodeChunk(hbeStream, BUFFERLENGTH/2, (short *) pWAVBuffer, pMP3Buffer, &dwWrite);
			
			if(err != BE_ERR_SUCCESSFUL)
			{
				beCloseStream(hbeStream);
				fprintf(stdout,"beEncodeChunk() failed (%lu)", err);
				return -1;
			}
			
			fwrite(pMP3Buffer,1,dwWrite, m2aFile);		
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
	char StrInfoBuffer[BUFSIZE_StrInfo];

	memset(StrInfoBuffer, 0, BUFSIZE_StrInfo);
	sprintf(StrInfoBuffer, "\rProgress:\t  %d%%", iPercent);
	
	WriteFile( stdout, StrInfoBuffer, BUFSIZE_StrInfo, &dwWrite, 0);
	FlushFileBuffers(stdout);
	
	//fprintf(stdout, "\rProgress:\t%d%%", iPercent);
	//fflush(stdout);
	
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
int main(int argc, char* argv[])
{
    CoInitialize( NULL );
	//Get input params
    if (argc <= 1)
	{
		fprintf(stderr, "usage: %s", usageString);                         
		exit(EXIT_FAILURE);
	}
	
	memcpy(pFileInput, argv[1], 100);
////////////////////////////////////////////////////////////////////////////////////////////////////////////Lirh 2008.10
	memcpy(pFileOutputAudio, argv[2], 100);
//	strcat(pFileOutputAudio, ".mp3");
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
		multiple = sqrt( pWF->nSamplesPerSec / ReSampleRate );
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
	
	FreeMediaType( mt1 );
	
	m2aFile = fopen(pFileOutputAudio, "wb");	//Open output file

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
	
	err = beDeinitStream(hbeStream, pMP3Buffer, &dwWrite);
	
	if (dwWrite)
	{
		fwrite( pMP3Buffer, 1, dwWrite, m2aFile);
	}
	
	beCloseStream( hbeStream );
	
	// Delete WAV buffer
	delete [] pWAVBuffer;
	
	// Delete MP3 Buffer
	delete [] pMP3Buffer;
	
	fclose(m2aFile);

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