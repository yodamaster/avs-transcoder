/*
*****************************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2003, Advanced Audio Video Coding Standard, Part II
*
* DISCLAIMER OF WARRANTY
*
* The contents of this file are subject to the Mozilla Public License
* Version 1.1 (the "License"); you may not use this file except in
* compliance with the License. You may obtain a copy of the License at
* http://www.mozilla.org/MPL/
*
* Software distributed under the License is distributed on an "AS IS"
* basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
* License for the specific language governing rights and limitations under
* the License.
*                     
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
* The AVS Working Group doesn't represent or warrant that the programs
* furnished here under are free of infringement of any third-party patents.
* Commercial implementations of AVS, including shareware, may be
* subject to royalty fees to patent holders. Information regarding
* the AVS patent policy for standardization procedure is available at 
* AVS Web site http://www.avs.org.cn. Patent Licensing is outside
* of AVS Working Group.
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
************************************************************************
*/

/*
*************************************************************************************
* File name: 
* Function: 
*
*************************************************************************************
*/
#include "transcoder.h"
#include "multithead.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <memory.h>
#include <string.h>
#include <BaseTsd.h>
#include <Windows.h>
#include <WinBase.h>
#include <WinNT.h>
#include <assert.h>

#include <atlbase.h>
#include <streams.h>
#include <qedit.h>         // for Null Renderer
#include <Dvdmedia.h>


#include "grabber.h"
#include "Convert.h"

#include "wmsdkidl.h"
#include "wmsdk.h"
#include "rostream.h"

#define TLS         __declspec(thread)
#define MAX_THREAD_NUM 4
#define MAX_IMG_WIDTH     352
#define MAX_IMG_HEIGHT    288

#pragma warning (disable:4005)
#pragma warning (disable:4018)
#pragma warning (disable:4244)


PTHREAD   p_thread;
PWTHREAD  p_wthread;
HANDLE thread_handle[MAX_THREAD_NUM];
HANDLE pEvtReady[MAX_THREAD_NUM];
HANDLE write_handle;
HANDLE write_evt;
HANDLE hMutex;
unsigned char *write_buf[MAX_THREAD_NUM];
int write_buf_flag[MAX_THREAD_NUM];
int iThreadnum = 4;
int tmp[MAX_THREAD_NUM];
int iGOPlen = 15;
int B_num = 2;
int last_thread_num = 0;
int last_frm_num = 0;

int ttfrms=0,cdfrms=0;

__declspec( thread ) c_avs_enc    **p_c_avs_enc;
unsigned char ***pYUVbuf;
InputParameters   inputs, *input;

ColorSpaceConversions conv;
double        length;
bool      bUDF_ASFBitrate = false;      //Flag, True if user specified the Bitrate for converting an asf video, otherwise, False.
bool          bFramerateFlag = false;
int           iOrgWidth;
int           iOrgHeight;
unsigned char*      pYUVBuffer;
unsigned char*      tmpBuffer;
FILE* fpFrameTime;

HRESULT GetPin(IBaseFilter * pFilter, PIN_DIRECTION dirrequired,  int iNum, IPin **ppPin);
IPin *  GetInPin ( IBaseFilter *pFilter, int Num );
IPin *  GetOutPin( IBaseFilter *pFilter, int Num );
int     GetEncodeParams(int argc, char* argv[], InputParameters* input);
void  YUV_Scale_v1(unsigned char* srcBuffer, unsigned char *dstBuffer, int wsrc, int hsrc, int wdst, int hdst);
void  YUV_Scale_v2(unsigned char* srcBuffer, unsigned char *dstBuffer, int wsrc, int hsrc, int wdst, int hdst); //wangyue20081107
int    GetVideoInfo(char* FileName, double* fps, double* length);
int    GetASFVidoInfo(char* FileName);
HRESULT GetStreamNumbers(IWMProfile* pProfile);

int *         dsmptmp1;
int *         dsmptmp2;
int *         dsmptmp3;
int *         dsmptmp4;

const int filter16[8][16][12] = {   // sine, N = 3
  { // D = 1
    {0,0,0,0,0,128,0,0,0,0,0,0},
    {0,0,0,2,-6,127,7,-2,0,0,0,0},
    {0,0,0,3,-12,125,16,-5,1,0,0,0},
    {0,0,0,4,-16,120,26,-7,1,0,0,0},
    {0,0,0,5,-18,114,36,-10,1,0,0,0},
    {0,0,0,5,-20,107,46,-12,2,0,0,0},
    {0,0,0,5,-21,99,57,-15,3,0,0,0},
    {0,0,0,5,-20,89,68,-18,4,0,0,0},
    {0,0,0,4,-19,79,79,-19,4,0,0,0},
    {0,0,0,4,-18,68,89,-20,5,0,0,0},
    {0,0,0,3,-15,57,99,-21,5,0,0,0},
    {0,0,0,2,-12,46,107,-20,5,0,0,0},
    {0,0,0,1,-10,36,114,-18,5,0,0,0},
    {0,0,0,1,-7,26,120,-16,4,0,0,0},
    {0,0,0,1,-5,16,125,-12,3,0,0,0},
    {0,0,0,0,-2,7,127,-6,2,0,0,0}
  },
  { // D = 1.5
    {0,2,0,-14,33,86,33,-14,0,2,0,0},
    {0,1,1,-14,29,85,38,-13,-1,2,0,0},
    {0,1,2,-14,24,84,43,-12,-2,2,0,0},
    {0,1,2,-13,19,83,48,-11,-3,2,0,0},
    {0,0,3,-13,15,81,53,-10,-4,3,0,0},
    {0,0,3,-12,11,79,57,-8,-5,3,0,0},
    {0,0,3,-11,7,76,62,-5,-7,3,0,0},
    {0,0,3,-10,3,73,65,-2,-7,3,0,0},
    {0,0,3,-9,0,70,70,0,-9,3,0,0},
    {0,0,3,-7,-2,65,73,3,-10,3,0,0},
    {0,0,3,-7,-5,62,76,7,-11,3,0,0},
    {0,0,3,-5,-8,57,79,11,-12,3,0,0},
    {0,0,3,-4,-10,53,81,15,-13,3,0,0},
    {0,0,2,-3,-11,48,83,19,-13,2,1,0},
    {0,0,2,-2,-12,43,84,24,-14,2,1,0},
    {0,0,2,-1,-13,38,85,29,-14,1,1,0}
  },

  { // D = 2
    {2,0,-10,0,40,64,40,0,-10,0,2,0},
    {2,1,-9,-2,37,64,42,2,-10,-1,2,0},
    {2,1,-9,-3,34,64,44,4,-10,-1,2,0},
    {2,1,-8,-5,31,63,47,6,-10,-2,3,0},
    {1,2,-8,-6,29,62,49,8,-10,-2,3,0},
    {1,2,-7,-7,26,61,52,10,-10,-3,3,0},
    {1,2,-6,-8,23,60,54,13,-10,-4,3,0},
    {1,2,-6,-9,20,59,56,15,-10,-4,3,1},
    {1,2,-5,-9,18,57,57,18,-9,-5,2,1},
    {1,3,-4,-10,15,56,59,20,-9,-6,2,1},
    {0,3,-4,-10,13,54,60,23,-8,-6,2,1},
    {0,3,-3,-10,10,52,61,26,-7,-7,2,1},
    {0,3,-2,-10,8,49,62,29,-6,-8,2,1},
    {0,3,-2,-10,6,47,63,31,-5,-8,1,2},
    {0,2,-1,-10,4,44,64,34,-3,-9,1,2},
    {0,2,-1,-10,2,42,64,37,-2,-9,1,2}
  },

  { // D = 2.5

    {0,-4,-7,11,38,52,38,11,-7,-4,0,0},
    {0,-4,-7,9,37,51,40,13,-6,-7,2,0},
    {0,-3,-7,8,35,51,41,14,-5,-7,1,0},
    {0,-2,-8,6,33,51,42,16,-5,-7,2,0},
    {0,-2,-8,5,32,50,43,18,-4,-8,2,0},
    {0,-2,-8,4,30,50,45,19,-3,-8,1,0},
    {0,-1,-8,2,28,49,46,21,-2,-8,1,0},
    {0,-1,-8,1,26,49,47,23,-1,-8,0,0},
    {0,0,-8,0,24,48,48,24,0,-8,0,0},
    {0,0,-8,-1,23,47,49,26,1,-8,-1,0},
    {0,1,-8,-2,21,46,49,28,2,-8,-1,0},
    {0,1,-8,-3,19,45,50,30,4,-8,-2,0},
    {0,2,-8,-4,18,43,50,32,5,-8,-2,0},
    {0,2,-7,-5,16,42,51,33,6,-8,-2,0},
    {0,1,-7,-5,14,41,51,35,8,-7,-3,0},
    {0,2,-7,-6,13,40,51,37,9,-7,-4,0}


  },
  { // D = 3
    {-2,-7,0,17,35,43,35,17,0,-7,-5,2},
    {-2,-7,-1,16,34,43,36,18,1,-7,-5,2},
    {-1,-7,-1,14,33,43,36,19,1,-6,-5,2},
    {-1,-7,-2,13,32,42,37,20,3,-6,-5,2},
    {0,-7,-3,12,31,42,38,21,3,-6,-5,2},
    {0,-7,-3,11,30,42,39,23,4,-6,-6,1},
    {0,-7,-4,10,29,42,40,24,5,-6,-6,1},
    {1,-7,-4,9,27,41,40,25,6,-5,-6,1},
    {1,-6,-5,7,26,41,41,26,7,-5,-6,1},
    {1,-6,-5,6,25,40,41,27,9,-4,-7,1},
    {1,-6,-6,5,24,40,42,29,10,-4,-7,0},
    {1,-6,-6,4,23,39,42,30,11,-3,-7,0},
    {2,-5,-6,3,21,38,42,31,12,-3,-7,0},
    {2,-5,-6,3,20,37,42,32,13,-2,-7,-1},
    {2,-5,-6,1,19,36,43,33,14,-1,-7,-1},
    {2,-5,-7,1,18,36,43,34,16,-1,-7,-2}
  },
  { // D = 3.5
    {-6,-3,5,19,31,36,31,19,5,-3,-6,0},
    {-6,-4,4,18,31,37,32,20,6,-3,-6,-1},
    {-6,-4,4,17,30,36,33,21,7,-3,-6,-1},
    {-5,-5,3,16,30,36,33,22,8,-2,-6,-2},
    {-5,-5,2,15,29,36,34,23,9,-2,-6,-2},
    {-5,-5,2,15,28,36,34,24,10,-2,-6,-3},
    {-4,-5,1,14,27,36,35,24,10,-1,-6,-3},
    {-4,-5,0,13,26,35,35,25,11,0,-5,-3},
    {-4,-6,0,12,26,36,36,26,12,0,-6,-4},
    {-3,-5,0,11,25,35,35,26,13,0,-5,-4},
    {-3,-6,-1,10,24,35,36,27,14,1,-5,-4},
    {-3,-6,-2,10,24,34,36,28,15,2,-5,-5},
    {-2,-6,-2,9,23,34,36,29,15,2,-5,-5},
    {-2,-6,-2,8,22,33,36,30,16,3,-5,-5},
    {-1,-6,-3,7,21,33,36,30,17,4,-4,-6},
    {-1,-6,-3,6,20,32,37,31,18,4,-4,-6}
  },
  { // D = 4
    {-9,0,9,20,28,32,28,20,9,0,-9,0},
    {-9,0,8,19,28,32,29,20,10,0,-4,-5},
    {-9,-1,8,18,28,32,29,21,10,1,-4,-5},
    {-9,-1,7,18,27,32,30,22,11,1,-4,-6},
    {-8,-2,6,17,27,32,30,22,12,2,-4,-6},
    {-8,-2,6,16,26,32,31,23,12,2,-4,-6},
    {-8,-2,5,16,26,31,31,23,13,3,-3,-7},
    {-8,-3,5,15,25,31,31,24,14,4,-3,-7},
    {-7,-3,4,14,25,31,31,25,14,4,-3,-7},
    {-7,-3,4,14,24,31,31,25,15,5,-3,-8},
    {-7,-3,3,13,23,31,31,26,16,5,-2,-8},
    {-6,-4,2,12,23,31,32,26,16,6,-2,-8},
    {-6,-4,2,12,22,30,32,27,17,6,-2,-8},
    {-6,-4,1,11,22,30,32,27,18,7,-1,-9},
    {-5,-4,1,10,21,29,32,28,18,8,-1,-9},
    {-5,-4,0,10,20,29,32,28,19,8,0,-9}
  },
  { // D = 6
    {-6,8,13,18,20,22,20,18,13,8,4,-10},
    {-6,8,13,17,20,21,20,18,13,9,4,-9},
    {-6,8,12,17,20,21,20,18,14,9,4,-9},
    {-7,7,12,17,20,21,21,18,14,9,5,-9},
    {-7,7,12,16,20,21,21,18,14,10,5,-9},
    {-7,7,12,16,20,21,21,18,14,10,5,-9},
    {-8,7,11,16,20,21,21,19,15,10,5,-9},
    {-8,6,11,16,19,21,21,19,15,11,6,-9},
    {-8,6,11,15,19,21,21,19,15,11,6,-8},
    {-9,6,11,15,19,21,21,19,16,11,6,-8},
    {-9,5,10,15,19,21,21,20,16,11,7,-8},
    {-9,5,10,14,18,21,21,20,16,12,7,-7},
    {-9,5,10,14,18,21,21,20,16,12,7,-7},
    {-9,5,9,14,18,21,21,20,17,12,7,-7},
    {-9,4,9,14,18,20,21,20,17,12,8,-6},
    {-9,4,9,13,18,20,21,20,17,13,8,-6}
  }
};
// ----------------------------------------------------------------------------
// get current time in milli-second
// ----------------------------------------------------------------------------
static __inline unsigned int
GetTime()
{
  struct _timeb time_ms;

  _ftime(&time_ms);            // start time ms
  return ((unsigned int)(time_ms.time * 1000 + time_ms.millitm));
}

HRESULT Callback(IMediaSample* pSample, REFERENCE_TIME* StartTime, REFERENCE_TIME* StopTime, BOOL TypeChanged)
{
  int i;
  int iPercent;
  unsigned char *pBuffer;
  static DWORD iRet = 0;
  static int iCount = 0;
  static int Frm_num = 0;
  REFERENCE_TIME iLastFrameTime;

  int iFlag=0;
  static REFERENCE_TIME lastst=*StartTime;
  const int iPeriod = (int)(10000000/input->fr);

  if(!input->fr)
  {
    iFlag=1;
  }
  else
  {
    REFERENCE_TIME st=((REFERENCE_TIME)(*StartTime/iPeriod))*iPeriod;
    REFERENCE_TIME ed=st+(REFERENCE_TIME)iPeriod;
    if(abs(*StartTime-st)<=(*StartTime-lastst)/2 || abs(*StartTime-ed)<(*StartTime-lastst)/2)
    {
      iFlag=1;
      cdfrms++;
    }
    lastst=*StartTime;
  }

  ttfrms++;

  if(iFlag)
  {
    pSample->GetPointer(&pBuffer);
    const int iSampleSize = pSample->GetActualDataLength();
    if (iCount == 0)
    {
      if (Frm_num < iThreadnum * iGOPlen)
        iRet = Frm_num / iGOPlen;
      else
        iRet = WaitForMultipleObjects(iThreadnum, pEvtReady, FALSE, INFINITE) - WAIT_OBJECT_0;
      tmp[iRet] = Frm_num;
    }
    ////Trans RGB24 to IYUV
    //conv.RGB24_to_YV12(pBuffer, pYUVBuffer, iOrgWidth, iOrgHeight);
    conv.YUY2_to_YV12(pBuffer, pYUVBuffer, iOrgWidth, iOrgHeight);
    ////xzhao YVU==>YUV?
    memcpy(tmpBuffer, pYUVBuffer + iOrgWidth * iOrgHeight, iOrgWidth * iOrgHeight / 4);
    memcpy(pYUVBuffer + iOrgWidth * iOrgHeight, pYUVBuffer + iOrgWidth * iOrgHeight * 5 / 4, iOrgWidth * iOrgHeight / 4);
    memcpy(pYUVBuffer + iOrgWidth * iOrgHeight * 5 / 4, tmpBuffer, iOrgWidth * iOrgHeight / 4);

    if (iCount % (B_num + 1) == 0 && iCount!=0)
      YUV_Scale_v2(pYUVBuffer, pYUVbuf[iRet][iCount - B_num], iOrgWidth, iOrgHeight, input->img_width, input->img_height);
    else if (iCount == 0 || ((iCount > iGOPlen - B_num - 1) && iGOPlen%(B_num+1)!=1))         //xzhao 20081107
      YUV_Scale_v2(pYUVBuffer, pYUVbuf[iRet][iCount], iOrgWidth, iOrgHeight, input->img_width, input->img_height);
    else
      YUV_Scale_v2(pYUVBuffer, pYUVbuf[iRet][iCount + 1], iOrgWidth, iOrgHeight, input->img_width, input->img_height);


    Frm_num++;
    iCount++;
    p_thread[iRet].p_enc->dec_frm_num = iCount;
    last_thread_num = iRet;
    last_frm_num = iCount;

    if (iCount == iGOPlen)
    {
      //p_thread[iRet].m_bRunning = TRUE;
      WaitForSingleObject(hMutex, INFINITE);
      for (i=0; i<iThreadnum*2; i++)
      {
        if (p_wthread->order[i][0] == -1)
        {
          p_wthread->order[i][0] = iRet;
          break;
        }
      }
      ReleaseMutex(hMutex);
      ResumeThread(thread_handle[iRet]);
      iCount = 0;
      last_frm_num = 0;
      iLastFrameTime = (*StartTime) / 10000;
      iPercent = (int)((*StartTime) / (length * 100000));
      fprintf(stderr, "\rFrame: %d\t\t", Frm_num);
      fprintf(stderr, "Start Time(ms): %d\t\t", iLastFrameTime);
      fprintf(stderr, "%d%%", iPercent);
      fflush(stderr);
    }
    iLastFrameTime = (*StartTime) / 10000;
    fwrite(&iLastFrameTime, sizeof(int), 1, fpFrameTime);  //写时间戳
    return S_OK;
  }
  else
  {
    return S_OK;
  }
}

int main(int argc, char* argv[])
{
  unsigned int iTime;
  int i, j;
  char time_stamp_file[100];

  c_avs_enc    **p_c_avs_enc;
  CoInitialize( NULL );

  input = &inputs;
  if (GetEncodeParams(argc, argv, input) == -1)
  {
    return 0;
  }

  double fRateTable[] = {24000.0/1001, 24, 25, 30000.0/1001, 30, 50, 60000.0/1001, 60};
  double fFrameRate = fRateTable[input->frame_rate_code - 1];    //FrameRateCode to FrameRate use Table-Driven Method

  double fps = 25;
  GetVideoInfo(input->infile, &fps, &length);  //Get fps, length of the video from file
  fps = (fps == 0 ? 25 : fps);

  if(!bFramerateFlag)
  {
    input->fr= (float)fps;
  }


  //读取asf,wmv文件的码率信息
  int iBitrate; //bit per second
  int iThreshold = 300000;  //BitRate Threshold

  //如果GetEncodeParams()已经指定了码率，则以用户指定的为准
  //否则，赋予码率值
  if (input->RCEnable == 1)
  {
    bUDF_ASFBitrate = TRUE;
  }
  if(!bUDF_ASFBitrate)
  {
    iBitrate = GetASFVidoInfo(input->infile);
    if (iBitrate <= 0)
    {
      printf("Failed in GetASFVidoInfo() Function!\r\n");
    }
    else
    {
      input->bit_rate = (int_32_t)(iBitrate > iThreshold ? iBitrate / fps * fFrameRate : iThreshold / fps * fFrameRate);  // convert to kbps
      bUDF_ASFBitrate = true;
    }
  }

  strcpy(time_stamp_file, input->outfile);
  char* str = strrchr(time_stamp_file, '.');
  *str = '\0';
  strcat(time_stamp_file, ".dat");
  if ((fpFrameTime = fopen(time_stamp_file, "wb")) == NULL)
  {
    printf("Could not open time stamp file!\r\n");
    return -1;
  }

  // The sample grabber is not in the registry, so create it with 'new'.
  HRESULT hr = S_OK;
  CComPtr< IGrabberSample > pGrab;

  pGrab.CoCreateInstance( CLSID_GrabberSample );

  // Set the callback function of the filter.
  pGrab->SetCallback(&Callback);


  // Set up a partially specified media type.
  CMediaType mt;

  mt.SetType(&MEDIATYPE_Video);
  // xzhao
  mt.SetSubtype(&MEDIASUBTYPE_YUY2);
  //mt.SetSubtype(&MEDIASUBTYPE_IYUV);
  //mt.SetSubtype(&MEDIASUBTYPE_RGB24);

  hr = pGrab->SetAcceptedMediaType(&mt);

  // Create the filter graph manager.
  CComPtr<IFilterGraph> pGraph;
  hr = pGraph.CoCreateInstance( CLSID_FilterGraph );

  // Query for other useful interfaces.
  CComQIPtr<IGraphBuilder, &IID_IGraphBuilder> pBuilder(pGraph);
  //    CComQIPtr<IMediaSeeking, &IID_IMediaSeeking> pSeeking(pGraph);
  CComQIPtr<IMediaControl, &IID_IMediaControl> pControl(pGraph);
  CComQIPtr<IMediaFilter, &IID_IMediaFilter> pMediaFilter(pGraph);
  CComQIPtr<IMediaEvent, &IID_IMediaEvent> pEvent(pGraph);

  // Add a source filter to the graph.
  CComPtr<IBaseFilter> pSource;
  hr = pBuilder->AddSourceFilter( A2WBSTR(input->infile), L"Source", &pSource);

  if( FAILED( hr ) )
  {
    printf( "Could not find the input file\r\n");
    return -1;
  }


  // Add the sample grabber to the graph.
  CComQIPtr< IBaseFilter, &IID_IBaseFilter > pGrabberBase( pGrab );

  hr = pBuilder->AddFilter(pGrabberBase, L"Grabber");

  //Enumerating Pins of Source filter, and find number of pins
  IEnumPins  *pEnum;
  IPin       *pPin;
  int  num = 0;

  hr = pSource->EnumPins(&pEnum);
  if (FAILED(hr))
  {
    printf("Could not Enumerate Pins of Source filter\r\n");
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
    if (hr == S_OK)    //some times hr = VFW_S_PARTIAL_RENDER ?
    {
      break;
    }
  }
  if( FAILED( hr ) )
  {
    printf("Could not connect source filter to grabber\r\n");
    //    fpErrorInfo = fopen("VideoErrorInfo.txt", "a");
    SYSTEMTIME time;
    TCHAR szDate[64], szTime[64];
    GetLocalTime(&time);
    GetDateFormat (LOCALE_USER_DEFAULT, LOCALE_NOUSEROVERRIDE | DATE_SHORTDATE,
      &time, NULL, szDate, sizeof (szDate)) ;
    GetTimeFormat ( LOCALE_USER_DEFAULT, LOCALE_NOUSEROVERRIDE |
      TIME_NOTIMEMARKER | TIME_FORCE24HOURFORMAT,
      &time, NULL, szTime, sizeof (szTime)) ;

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
  //FORMAT_VideoInfo
  //Get video info(weight,height,bitrate)
  if(mt1.formattype == FORMAT_VideoInfo)
  {
    VIDEOINFOHEADER * vih = (VIDEOINFOHEADER*) mt1.pbFormat;
    if (!bUDF_ASFBitrate && vih->dwBitRate > iThreshold && vih->dwBitRate < 10000000)
    {
      input->bit_rate = (int)(vih->dwBitRate / fps * fFrameRate);
    }
    else if (!bUDF_ASFBitrate && vih->dwBitRate <= iThreshold && vih->dwBitRate != 0)
    {
      input->bit_rate = (int)(iThreshold / fps * fFrameRate);
    }
    else if (!bUDF_ASFBitrate && (vih->dwBitRate > 10000000 || vih->dwBitRate == 0))
    {
      printf("Wrong with Bitrate! Please set a Bitrate!\r\n");
    }

    if (input->img_width == 0 || input->img_height == 0)
    {
      // The width must be a multiples of 4!
      if((vih->bmiHeader.biWidth%4)!=0)
        iOrgWidth = input->img_width = 4*((int)(vih->bmiHeader.biWidth/4)+1);
      else
        iOrgWidth = input->img_width = vih->bmiHeader.biWidth;
      iOrgHeight = input->img_height = vih->bmiHeader.biHeight>0 ?vih->bmiHeader.biHeight:-vih->bmiHeader.biHeight;
      if((input->img_width%16)!=0)
        input->img_width = 16*((int)(iOrgWidth/16)+1);
    }
    else
    {
      iOrgWidth = vih->bmiHeader.biWidth;
      iOrgHeight = vih->bmiHeader.biHeight>0 ?vih->bmiHeader.biHeight:-vih->bmiHeader.biHeight;
    }
  }
  else if(mt1.formattype == FORMAT_VideoInfo2)
  {
    VIDEOINFOHEADER2 * vih = (VIDEOINFOHEADER2*) mt1.pbFormat;
    if (!bUDF_ASFBitrate && vih->dwBitRate > iThreshold && vih->dwBitRate < 10000000)
    {
      input->bit_rate = (int)(vih->dwBitRate / fps * fFrameRate);
    }
    else if (!bUDF_ASFBitrate && vih->dwBitRate <= iThreshold && vih->dwBitRate != 0)
    {
      input->bit_rate = (int)(iThreshold / fps * fFrameRate);
    }
    else if (!bUDF_ASFBitrate && (vih->dwBitRate > 10000000 || vih->dwBitRate == 0))
    {
      printf("Wrong with Bitrate! Please set a Bitrate!\r\n");
    }

    if (input->img_width == 0 || input->img_height == 0)
    {
      // The width must be a multiples of 4!
      if((vih->bmiHeader.biWidth%4)!=0)
        iOrgWidth = input->img_width = 4*((int)(vih->bmiHeader.biWidth/4)+1);
      else
        iOrgWidth = input->img_width = vih->bmiHeader.biWidth;
      iOrgHeight = input->img_height = vih->bmiHeader.biHeight>0 ?vih->bmiHeader.biHeight:-vih->bmiHeader.biHeight;


      if((input->img_width%16)!=0)
        input->img_width = 16*((int)(iOrgWidth/16)+1);
    }
    else
    {
      iOrgWidth = vih->bmiHeader.biWidth;
      iOrgHeight = vih->bmiHeader.biHeight>0 ?vih->bmiHeader.biHeight:-vih->bmiHeader.biHeight;
    }
  }

  FreeMediaType( mt1 );

  iThreadnum = input->thread_num;
  iGOPlen = input->GopLength;
  B_num = input->successive_Bframe;
  printf("Encoding Thread Number: %d\n", iThreadnum);

  pYUVbuf = new unsigned char** [iThreadnum];

  p_c_avs_enc   = new c_avs_enc* [iThreadnum];

  //为转码器需要的对象准备必要的内存
  /*avs_enc*/
  for (i=0; i<iThreadnum; i++)
  {
    p_c_avs_enc[i] = new c_avs_enc;
    p_c_avs_enc[i]->input = &(p_c_avs_enc[i]->inputs);
    memcpy(p_c_avs_enc[i]->input, input, sizeof(InputParameters));
    p_c_avs_enc[i]->p_avs_enc_frame = (avs_enc_frame_t*)malloc(sizeof(avs_enc_frame_t));
    p_c_avs_enc[i]->p_avs_enc_frame->bitstream = (unsigned char *)malloc(input->img_width*input->img_height*4*sizeof(char) * iGOPlen);
    if (p_c_avs_enc[i]->p_avs_enc_frame->bitstream == NULL)
    {
      printf("error in alloc memory\n");
      exit(0);
    }
    write_buf[i] = (unsigned char *)malloc(input->img_width*input->img_height*4*sizeof(char) * iGOPlen);
    pYUVbuf[i] = new unsigned char* [iGOPlen];
    for (j=0; j<iGOPlen; j++)
    {
      pYUVbuf[i][j] = (unsigned char *) malloc(input->img_width * input->img_height * 3 /2);
    }
    p_c_avs_enc[i]->p_avs_enc_frame->inputfrm   = pYUVbuf[i];
    write_buf_flag[i] = 0;
  }

  pYUVBuffer = (unsigned char *) malloc(iOrgWidth * iOrgHeight * 3 /2);
  tmpBuffer = (unsigned char *) malloc(iOrgWidth * iOrgHeight / 4);

  p_thread = (PTHREAD)malloc(iThreadnum * sizeof(THREAD));
  p_wthread = (PWTHREAD)malloc(sizeof(WTHREAD));
  //创建转码线程
  create_multithread_transcoder(iThreadnum, p_c_avs_enc);

  hr = pMediaFilter->SetSyncSource(NULL);

  if( FAILED( hr ) )
  {
    printf("Could not SetSyncSource pMediaFilter\r\n");
    return -1;
  }

  iTime = GetTime();

  hr = pControl->Run( );

  if( FAILED( hr ) )
  {
    printf("Could not Run\r\n");
    return -1;
  }

  long EvCode = 0;

  hr = pEvent->WaitForCompletion(INFINITE, &EvCode);

  if (last_frm_num > 0)
    ResumeThread(thread_handle[last_thread_num]);

  WaitForMultipleObjects(iThreadnum, pEvtReady, TRUE, INFINITE);
  p_wthread->m_bRunning = FALSE;
  ResumeThread(write_handle);
  for (i=0; i<iThreadnum; i++)
  {
    p_thread[i].m_bRunning = FALSE;
    ResumeThread(thread_handle[i]);
  }

  WaitForMultipleObjects(iThreadnum, thread_handle, TRUE, 3000);
  WaitForSingleObject(write_handle, INFINITE);

  iTime = GetTime() - iTime;
  fprintf(stderr,"\n:::::: TOTAL ENCODING TIME(ms): %u ::::::\n", iTime);
  fflush(stderr);

  for (i=0; i<iThreadnum; i++)
  {
    //free(p_c_avs_enc[i]->p_stats);
    free(write_buf[i]);
    free(p_c_avs_enc[i]->p_avs_enc_frame->bitstream);
    free(p_c_avs_enc[i]->p_avs_enc_frame);
    delete p_c_avs_enc[i];
    for (j=0; j<iGOPlen; j++)
    {
      free(pYUVbuf[i][j]);
      //free(pMBinfo[i][j]);
    }
    delete [] pYUVbuf[i];
    //delete [] pMBinfo[i];
    CloseHandle(thread_handle[i]);
    CloseHandle(pEvtReady[i]);
  }

  printf("TOTAL FRAMES: %d, CODED FRAMES: %d",ttfrms,cdfrms);

  CloseHandle(write_evt);
  CloseHandle(write_handle);
  delete [] p_c_avs_enc;
  delete [] pYUVbuf;
  //delete [] pMBinfo;
  free(pYUVBuffer);
  free(tmpBuffer);
  free(p_thread);
  free(p_wthread);

  fclose(fpFrameTime);

  CoUninitialize();
  pGrab.p->Release();
  pControl.p->Release();
  pMediaFilter.p->Release();
  pEvent.p->Release();
  pSource.p->Release();
  pGrabberBase.p->Release();
  pNull.p->Release();
  pBuilder.p->Release();
  exit(0);
  return 0;
}

DWORD WINAPI run_thread(LPVOID pArg)
{
  int n;
  avs_encoder_create(((PTHREAD)pArg)->p_enc);
  while (1)
  {
    n = ((PTHREAD)pArg)->iThreadOrder;
    ((PTHREAD)pArg)->p_enc->current_encoded_frame = tmp[n];
    avs_encoder_encode(((PTHREAD)pArg)->p_enc);
    memcpy(write_buf[n], ((PTHREAD)pArg)->p_enc->p_avs_enc_frame->bitstream, ((PTHREAD)pArg)->p_enc->p_avs_enc_frame->length);
    write_buf_flag[n] = 1;
    p_wthread->nbit[n] = ((PTHREAD)pArg)->p_enc->p_avs_enc_frame->length;

	ResumeThread(write_handle);

    SetEvent(((PTHREAD)pArg)->m_hEvent);
    SuspendThread(((PTHREAD)pArg)->m_hThread);
    if(!((PTHREAD)pArg)->m_bRunning)
    {
      avs_encoder_destroy(((PTHREAD)pArg)->p_enc);
      ExitThread(0);
    }
  }
  return 1;
}

DWORD WINAPI write_thread(LPVOID pArg)
{
  int n, i;
  int end_code = 0xb1010000;
  FILE* f = fopen(((PWTHREAD)pArg)->output_file_name, "wb");
  if(f==NULL)
  {
    printf ("\nCan't open file %s",((PWTHREAD)pArg)->output_file_name);
    exit(-1);
  }

  while (1)
  {
    n = ((PWTHREAD)pArg)->order[0][0];
    if (n != -1 && write_buf_flag[n] == 1)
    {
      fwrite(write_buf[n], 1, ((PWTHREAD)pArg)->nbit[n], f);
      WaitForSingleObject(hMutex, INFINITE);
      for (i=0; i<iThreadnum*2-1; i++)
      {
        ((PWTHREAD)pArg)->order[i][0] = ((PWTHREAD)pArg)->order[i+1][0];
      }
      ((PWTHREAD)pArg)->order[iThreadnum*2-1][0] = -1;
      write_buf_flag[n] = 0;
      ReleaseMutex(hMutex);
      SetEvent(((PWTHREAD)pArg)->m_hEvent);
    }
    else
      SuspendThread(((PTHREAD)pArg)->m_hThread);
    if(!((PWTHREAD)pArg)->m_bRunning)
    {
      for (i=0; i<iThreadnum; i++)
      {
        n = ((PWTHREAD)pArg)->order[i][0];
        if (write_buf_flag[n] == 1)
        {
          fwrite(write_buf[n], 1, ((PWTHREAD)pArg)->nbit[n], f);
        }
      }
      if (last_frm_num > 0)
        fwrite(write_buf[last_thread_num], 1, ((PWTHREAD)pArg)->nbit[last_thread_num], f);
      fwrite(&end_code, 1, sizeof(int), f);
      fclose(f);
      SetEvent(((PWTHREAD)pArg)->m_hEvent);
      ExitThread(0);
    }
  }
  return 1;
}

void create_multithread_transcoder(int thread_num, c_avs_enc **p_c_avs_enc)
{
  int i;

  for (i=0; i<thread_num; i++)
  {
    p_thread[i].m_bRunning = FALSE;
    pEvtReady[i]   = CreateEvent(NULL, FALSE, FALSE, NULL);
  }
  write_evt = CreateEvent(NULL, FALSE, FALSE, NULL);
  hMutex = CreateMutex(NULL, 0, NULL);
  for (i=0; i<thread_num; i++)
  {
    thread_handle[i] = CreateThread(
      NULL,                       // [in ] pointer to SECURITY_ATTRIBUTES structure
      0,                          // [in ] initial size of the stack, in bytes
      run_thread,               // [in ] thread function
      &p_thread[i],               // [in ] thread argument
      CREATE_SUSPENDED,      // [in ] additional creation flags
      &p_thread[i].nThreadId);    // [out] thread identifier
    SetThreadPriority(thread_handle[i], THREAD_PRIORITY_HIGHEST); //THREAD_PRIORITY_ABOVE_NORMAL
    p_thread[i].m_hThread = thread_handle[i];
    p_thread[i].m_hEvent = pEvtReady[i];
    p_thread[i].iThreadOrder = i;
    p_thread[i].p_enc = p_c_avs_enc[i];
    p_thread[i].m_bRunning = TRUE;
    p_wthread->order[i][0] = -1;
    p_wthread->order[i+iThreadnum][0] = -1;
  }
  write_handle = CreateThread(
    NULL,                       // [in ] pointer to SECURITY_ATTRIBUTES structure
    0,                          // [in ] initial size of the stack, in bytes
    write_thread,               // [in ] thread function
    p_wthread,               // [in ] thread argument
    CREATE_SUSPENDED,      // [in ] additional creation flags
    &p_wthread->nThreadId);    // [out] thread identifier
  SetThreadPriority(write_handle, THREAD_PRIORITY_HIGHEST); //THREAD_PRIORITY_ABOVE_NORMAL
  p_wthread->m_hThread = write_handle;
  p_wthread->m_bRunning = TRUE;
  p_wthread->m_hEvent = write_evt;
  strcpy(p_wthread->output_file_name, p_c_avs_enc[0]->inputs.outfile);
}

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

/***********************************************************************
Function: YUV_Scale
Input:  unsigned char* srcBuffer    Source File
unsigned char *dstBuffer    Destiny File
int wsrc, int hsrc        Original Size
int wdst, int hdst        New Size
Return: void
***********************************************************************/

void YUV_Scale_v2(unsigned char* srcBuffer, unsigned char *dstBuffer, 
                  int wsrc, int hsrc, int wdst, int hdst)
{
  if (wsrc == wdst && hsrc == hdst)
  {
    memcpy(dstBuffer, srcBuffer, wsrc * hsrc * 3 / 2);
    return;
  }

  double xscale = wsrc * 1.0 / wdst;
  double yscale = hsrc * 1.0 / hdst;
  int crop_w = wsrc;
  int crop_h = hsrc;
  int filterw, filterh;
  int x16,y16;
  int x, y;
  int i, j;
  int *px, *py;
  int m;

  dsmptmp1 = (int *) calloc(hsrc * wdst, sizeof(int));
  dsmptmp2 = (int *) calloc(max(hsrc,hdst) * wdst, sizeof(int));
  dsmptmp3 = (int *) calloc(hdst * wdst, sizeof(int));
  dsmptmp4 = (int *) calloc(hdst * wdst, sizeof(int));

  if(crop_w*7 > 20*wdst) filterw = 6;
  else if(crop_w*2 > 5*wdst) filterw = 5;
  else if(crop_w*1 > 2*wdst) filterw = 4;
  else if(crop_w*3 > 5*wdst) filterw = 3;
  else if(crop_w*4 > 5*wdst) filterw = 2;
  else if(crop_w*19 > 20*wdst) filterw = 1;
  else filterw = 0;


  if(crop_h*7 > 20*hdst) filterh = 6;
  else if(crop_h*2 > 5*hdst) filterh = 5;
  else if(crop_h*1 > 2*hdst) filterh = 4;
  else if(crop_h*3 > 5*hdst) filterh = 3;
  else if(crop_h*4 > 5*hdst) filterh = 2;
  else if(crop_h*19 > 20*hdst) filterh = 1;
  else filterh = 0;

  px = new int[wdst];
  py = new int[hdst];
  for( i = 0; i < wdst; i++ )
  {
    px[i] =  (i*crop_w*16 + 4*2*crop_w - 4*2*wdst + wdst/2) / wdst;
  }
  for( j = 0; j < hdst; j++ )
  {
    py[j] = ( j*crop_h*16 + 4*2*crop_h - 4*2*hdst + hdst/2) / hdst;
  }

  for( j = 0; j < hsrc; j++ ) 
  {
    //----- down sample row -----
    for(  i = 0; i < wdst; i++ )
    {
      x16 = px[i]&0x0f;
      x = px[i]>>4;
      dsmptmp1[j * wdst + i] = 0;
      for( int k = 0; k < 12; k++ )
      {
        m = x - 5 + k;
        if( m<0 ) m = 0;
        else if( m>(wsrc-1) ) m=wsrc-1;
        dsmptmp1[j * wdst + i] += filter16[filterw][x16][k] * ( *(srcBuffer+(j * wsrc + m)));
      }
    }
  }


  for( i = 0; i < wdst; i++ )
  {
    //----- down sample column -----
    for( j = 0; j < hdst; j++ )
    {
      y16 = py[j]&0x0f;
      y = py[j]>>4;
      dsmptmp2[j * wdst + i] = 0;
      for( int k = 0; k < 12; k++ )
      {
        m = y - 5 + k;
        if( m<0 ) m = 0;
        else if( m>(hsrc-1) ) m=hsrc-1;
        dsmptmp2[j * wdst + i] += filter16[filterh][y16][k] * dsmptmp1[m * wdst + i];
      }
    }

  }
  for (y = 0; y < hdst; y++)
  {
    for (x = 0; x < wdst; x++)
    {
      *(dstBuffer + (y * wdst + x)) = unsigned char(Clip3( 0, 255,(dsmptmp2[y * wdst + x]+(1<<13) ) / (1<<14))) ;
    }
  }
  delete [] px;
  delete [] py;


  // chroma
  wsrc >>= 1;
  hsrc >>= 1;
  wdst >>= 1;
  hdst >>= 1;
  crop_w = wsrc;
  crop_h = hsrc;

  px = new int[wdst];
  py = new int[hdst];
  for( i = 0; i < wdst; i++ )
  {
    px[i] =  (i*crop_w*16 + 4*crop_w - 4*wdst + wdst/2) / wdst;
  }
  for( j = 0; j < hdst; j++ )
  {
    py[j] = ( j*crop_h*16 + 4*2*crop_h - 4*2*hdst + hdst/2) / hdst;
  }

  for( j = 0; j < hsrc; j++ ) 
  {

    //----- down sample row -----
    for(  i = 0; i < wdst; i++ )
    {
      x16 = px[i]&0x0f;
      x = px[i]>>4;
      dsmptmp1[j * wdst + i] = 0;
      dsmptmp2[j * wdst + i] = 0;
      for( int k = 0; k < 12; k++ )
      {
        m = x - 5 + k;
        if( m<0 ) m = 0;
        else if( m>(wsrc-1) ) m=wsrc-1;
        dsmptmp1[j * wdst + i] += filter16[filterw][x16][k] * ( *(srcBuffer+(j * wsrc + (wsrc * hsrc * 4) +m)));
        dsmptmp2[j * wdst + i] += filter16[filterw][x16][k] * ( *(srcBuffer+(j * wsrc + (wsrc * hsrc * 5) +m)));
      }
    }
  }


  for( i = 0; i < wdst; i++ )
  {
    //----- down sample column -----
    for( j = 0; j < hdst; j++ )
    {
      y16 = py[j]&0x0f;
      y = py[j]>>4;
      dsmptmp3[j * wdst + i] = 0;
      dsmptmp4[j * wdst + i] = 0;

      for( int k = 0; k < 12; k++ )
      {
        m = y - 5 + k;
        if( m<0 ) m = 0;
        else if( m>(hsrc-1) ) m=hsrc-1;
        dsmptmp3[j * wdst + i] += filter16[filterh][y16][k] * dsmptmp1[m * wdst + i];
        dsmptmp4[j * wdst + i] += filter16[filterh][y16][k] * dsmptmp2[m * wdst + i];
      }
    }

  }


  for (y = 0; y < hdst; y++)
  {
    for (x = 0; x < wdst; x++)
    {
      *(dstBuffer + (wdst * hdst * 4) + (y * wdst + x)) = unsigned char(Clip3( 0, 255,(dsmptmp3[y * wdst + x]+(1<<13) ) / (1<<14))) ;
      *(dstBuffer + (wdst * hdst * 5) + (y * wdst + x)) = unsigned char(Clip3( 0, 255,(dsmptmp4[y * wdst + x]+(1<<13) ) / (1<<14))) ;
    }
  }
  delete [] px;
  delete [] py;

  // xzhao
  free(dsmptmp1);
  free(dsmptmp2);
  free(dsmptmp3);
  free(dsmptmp4);  
}

void YUV_Scale_v1(unsigned char* srcBuffer, unsigned char *dstBuffer,
                  int wsrc, int hsrc, int wdst, int hdst)
{
  if (wsrc == wdst && hsrc == hdst)
  {
    memcpy(dstBuffer, srcBuffer, wsrc * hsrc * 3 / 2);
    return;
  }

  double xscale = wsrc * 1.0 / wdst;
  double yscale = hsrc * 1.0 / hdst;


  int x, y;
  //int xsrc, ysrc;

  int m;


  int aiFilter[13]  = { 2, 0, -4, -3, 5, 19, 26, 19, 5, -3, -4, 0, 2 };

  dsmptmp1 = (int *) calloc(hsrc * wdst, sizeof(int));
  dsmptmp2 = (int *) calloc(max(hsrc,hdst) * wdst, sizeof(int));
  dsmptmp3 = (int *) calloc(hdst * wdst, sizeof(int));
  dsmptmp4 = (int *) calloc(hdst * wdst, sizeof(int));


  for( int y = 0; y < hsrc; y++ )
  {

    //----- down sample row -----
    for( int x = 0; x < wdst; x++ )
    {
      dsmptmp1[y * wdst + x] = 0;
      for( int k = 0; k < 13; k++ )
      {
        m = Clip3(0, wsrc - 1 ,int(x * xscale) + k - 6 );
        dsmptmp1[y * wdst + x] += aiFilter[k] * ( *(srcBuffer+(y * wsrc + m)));
      }
    }
  }


  for( int x = 0; x < wdst; x++ )
  {
    //----- down sample column -----
    for( int y = 0; y < hdst; y++ )
    {
      dsmptmp2[y * wdst + x] = 0;
      for( int k = 0; k < 13; k++ )
      {
        int m = Clip3( 0, hsrc - 1 ,int(y * yscale) + k - 6 );
        dsmptmp2[y * wdst + x] += aiFilter[k] * dsmptmp1[m * wdst + x];
      }
    }

  }
  for (y = 0; y < hdst; y++)
  {

    for (x = 0; x < wdst; x++)
    {


      *(dstBuffer + (y * wdst + x)) = unsigned char(Clip3( 0, 255,(dsmptmp2[y * wdst + x]+2048)>>12)) ;

    }
  }


  // chroma
  wsrc >>= 1;
  hsrc >>= 1;
  wdst >>= 1;
  hdst >>= 1;
  for( int y = 0; y < hsrc; y++ )
  {

    //----- down sample row -----
    for( int x = 0; x < wdst; x++ )
    {
      dsmptmp1[y * wdst + x] = 0;
      dsmptmp2[y * wdst + x] = 0;
      for( int k = 0; k < 13; k++ )
      {
        m = Clip3(0, wsrc - 1 ,int(x * xscale) + k - 6 );
        dsmptmp1[y * wdst + x] += aiFilter[k] * ( *(srcBuffer+(y * wsrc + (wsrc * hsrc * 4) +m)));
        dsmptmp2[y * wdst + x] += aiFilter[k] * ( *(srcBuffer+(y * wsrc + (wsrc * hsrc * 5) +m)));
      }
    }
  }


  for( int x = 0; x < wdst; x++ )
  {
    //----- down sample column -----
    for( int y = 0; y < hdst; y++ )
    {
      dsmptmp3[y * wdst + x] = 0;
      dsmptmp4[y * wdst + x] = 0;
      for( int k = 0; k < 13; k++ )
      {
        int m = Clip3( 0, hsrc - 1 ,int(y * yscale) + k - 6 );
        dsmptmp3[y * wdst + x] += aiFilter[k] * dsmptmp1[m * wdst + x];
        dsmptmp4[y * wdst + x] += aiFilter[k] * dsmptmp2[m * wdst + x];
      }
    }

  }

  for (y = 0; y < hdst; y++)
  {

    for (x = 0; x < wdst; x++)
    {


      *(dstBuffer + (wdst * hdst * 4) + (y * wdst + x)) = unsigned char(Clip3( 0, 255,(dsmptmp3[y * wdst + x]+2048)>>12)) ;
      *(dstBuffer + (wdst * hdst * 5) + (y * wdst + x)) = unsigned char(Clip3( 0, 255,(dsmptmp4[y * wdst + x]+2048)>>12)) ;

    }
  }


  // xzhao
  free(dsmptmp1);
  free(dsmptmp2);
  free(dsmptmp3);
  free(dsmptmp4);
}

int GetVideoInfo(char* FileName, double* fps, double* length)
{
  //////////////////////////////////////////////////////////////////////////
  /*获取源文件帧率等信息 hhan 2006-08-20
  **Obtain the authored frame rate of a video file without rendering it
  */
  HRESULT res;
  IMediaDet *pMediaDet;

  res = CoCreateInstance(CLSID_MediaDet, NULL, CLSCTX_INPROC, IID_IMediaDet,
    (void **)&pMediaDet);
  if(FAILED(res))
  {
    printf("CoCreateInstance(): res=0x%08x\n", res);
    return -1;
  }

  //  OLECHAR bstrFilename[MAX_PATH];
  //  MultiByteToWideChar(CP_ACP, 0, argv[1], -1, bstrFilename, sizeof(bstrFilename) / sizeof(OLECHAR));
  //  memcpy(bstrFilename, argv[1], MAX_PATH);
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
  printf("Length of Video: %gs\r\nAuthored Frame Rate: %gfps\r\n(%d frames of %dx%dx%d)\n", *length, *fps, nFrames, pbih->biWidth, pbih->biHeight, pbih->biBitCount);

  pMediaDet->Release();
  return 0;
}

/*
**获取ASF文件的码率信息
*/
int GetASFVidoInfo(char* FileName)
{
  HRESULT hr;

  IWMSyncReader*  m_pReader;
  CROStream*      m_pStream;
  //  IWMReader* pReader;
  IWMProfile* pProfile = NULL;
  //   IWMStreamConfig* ppConfig;

  //  IWMMetadataEditor* pEditor;
  IWMHeaderInfo* pInfo = NULL;

  //   hr = WMCreateEditor(&pEditor);
  //  hr = WMCreateReader(NULL, WMT_RIGHT_PLAYBACK OR WMT_RIGHT_COPY, &pReader);

  //  OLECHAR bstrFilename[MAX_PATH];
  //  MultiByteToWideChar(CP_ACP, 0, FileName, -1, bstrFilename, sizeof(bstrFilename) / sizeof(OLECHAR));

  hr = WMCreateSyncReader(  NULL, 0, &m_pReader );
  if ( FAILED( hr ) )
  {
    _tprintf( _T( "Could not create reader (hr=0x%08x).\n" ), hr );
    return( hr );
  }

  m_pStream = new CROStream;
  if( NULL == m_pStream )
  {
    hr = E_OUTOFMEMORY;
    _tprintf( _T( "Could not open file (hr=0x%08x).\n" ), hr );
    return( hr );
  }

  hr = m_pStream->Open( (LPCTSTR)FileName );

  if( FAILED( hr ) )
  {
    _tprintf( _T( "Could not open file (hr=0x%08x).\n" ) ,hr );
    return( hr );
  }

  hr = m_pReader->OpenStream( m_pStream );
  if ( FAILED( hr ) )
  {
    _tprintf( _T( "Could not open file (hr=0x%08x).\n" ), hr );
    return( hr );
  }


  hr = m_pReader->QueryInterface( IID_IWMProfile, ( VOID ** )&pProfile );
  if ( FAILED( hr ) )
  {
    _tprintf( _T(  "Could not QI for IWMProfile (hr=0x%08x).\n" ), hr );
    return( hr );
  }

  DWORD dwStreamCount = 0;
  hr = pProfile->GetStreamCount( &dwStreamCount );
  if ( FAILED( hr ) )
  {
    _tprintf( _T( "Could not get stream count: (hr=0x%08x)\n" ), hr );
    return( hr );
  }

  _tprintf( _T( "This Windows Media file has %d stream(s)\n" ), dwStreamCount );
  _tprintf( _T( "\n" ) );

  for ( DWORD dwIndex = 0; dwIndex < dwStreamCount; dwIndex++ )
  {
    IWMStreamConfig *pConfig = NULL;
    hr = pProfile->GetStream( dwIndex, &pConfig );
    if ( FAILED( hr ) )
    {
      _tprintf( _T( "Could not get the stream: (hr=0x%08x)\n" ), hr );
      return( hr );
    }

    GUID guid = GUID_NULL;
    hr = pConfig->GetStreamType( &guid );
    if ( FAILED( hr ) )
    {
      _tprintf( _T( "Could not get the stream type: (hr=0x%08x)\n" ), hr );
      return( hr );
    }
    else
    {
      if ( WMMEDIATYPE_Video == guid )
      {
        WORD wStreamNum = 0;
        hr = pConfig->GetStreamNumber( &wStreamNum );
        if ( FAILED( hr ) )
        {
          _tprintf( _T( "Could not get stream number: (hr=0x%08x)\n" ), hr );
          return( hr );
        }
        DWORD dwBitrate = 0;
        hr = pConfig->GetBitrate( &dwBitrate );
        if ( FAILED( hr ) )
        {
          _tprintf( _T( "Could not get bit rate: (hr=0x%08x)\n" ), hr );
          return( hr );
        }

        _tprintf( _T( "Video Stream properties:\n" ) );
        _tprintf( _T( "Stream number: %d\n" ), wStreamNum );
        _tprintf( _T( "Bitrate: %d bps \n" ), dwBitrate );
        //         hr = PrintCodecName( pConfig );
        _tprintf( _T( "\n" ) );
        return dwBitrate;
      }
      else if ( WMMEDIATYPE_Audio == guid )
      {
        WORD wStreamNum = 0;
        hr = pConfig->GetStreamNumber( &wStreamNum );
        if ( FAILED( hr ) )
        {
          _tprintf( _T( "Could not get stream number: (hr=0x%08x)\n" ), hr );
          return( hr );
        }
        DWORD dwBitrate = 0;
        hr = pConfig->GetBitrate( &dwBitrate );
        if ( FAILED( hr ) )
        {
          _tprintf( _T( "Could not get bit rate: (hr=0x%08x)\n" ), hr );
          return( hr );
        }
        _tprintf( _T( "Audio Stream properties:\n" ) );
        _tprintf( _T( "Stream number: %d\n" ), wStreamNum );
        _tprintf( _T( "Bitrate: %d bps \n" ), dwBitrate );
        //        hr = PrintCodecName( pConfig );
        _tprintf( _T( "\n" ) );
      }
      else if ( WMMEDIATYPE_Script == guid )
      {
        WORD wStreamNum = 0;
        hr = pConfig->GetStreamNumber( &wStreamNum );
        if ( FAILED( hr ) )
        {
          _tprintf( _T( "Could not get stream number: (hr=0x%08x)\n" ), hr );
          return( hr );
        }
        DWORD dwBitrate = 0;
        hr = pConfig->GetBitrate( &dwBitrate );
        if ( FAILED( hr ) )
        {
          _tprintf( _T( "Could not get bit rate: (hr=0x%08x)\n" ), hr );
          return( hr );
        }
        _tprintf( _T( "Script Stream properties:\n" ) );
        _tprintf( _T( "Stream number: %d\n" ), wStreamNum );
        _tprintf( _T( "Bitrate: %d bps \n" ), dwBitrate );
        _tprintf( _T( "\n" ) );
      }
    }
    pConfig->Release();
  }

  return hr;
}
/************************************************************************
Function: GetEncodeParams
Input: int argc, char* argv[]     arguments
InputParameters* input    pointer to encode param

Output: int     0 right, other wrong
************************************************************************/
int GetEncodeParams(int argc, char* argv[], InputParameters* input)
{
  int i;
  int err;

  //##########################################################################################
  //# Files
  //##########################################################################################
  //input->infile
  if (argc <= 1)
  {
    //output the help command
    printf("usage\n");
    printf("VideoTransCoder.exe input.m2v out.avs\n");
    printf("command options\n");
    printf("\t -w width\n");
    printf("\t -h height\n");
    printf("\t -c chroma_format\n");
    printf("\t -p progressive_frame (1(default):progressive, 0:interlace)\n");
    printf("\t -g gop_length (18:default)\n");
    printf("\t -r reference_frame_numer(2:default)\n");
    printf("\t -b b_frame_num(2:default)\n");
    printf("\t -q default_qp_for_i_b_p_frame(26:default)\n");
    printf("\t -a rate control enable flag(0:default, disable, 1:enable)\n");
    printf("\t -f frame_rate_code(2:default) Frame rate code [1)24000/1001, 2)24, 3)25, 4)30000/1001, 5)30, 6)50, 7)60000/1001, 8)60]\n");
    printf("\t -s bit_rate(1000000:default) in Kbps, s means speed\n");
    printf("\t -o rdo flag(1:default, enable, 0:disable)\n");
    printf("\t -m me search range(16:default)\n");
    printf("\t -n fps(default value is the decoded frame rate)\n");
    printf("\t -t thread number(1:default)\n");
    return -1;
  }

  if (*argv[1] != '-')
  {
    memcpy(input->infile, argv[1], 256);
  }
  else
  {
    printf("No input file!\r\n");
    return -1;
  }
  input->infile_header = 0;
  input->no_frames = 2;          //number of frames to be encoded (if -1, encoding forever until EOF)
  input->img_width = 0;          //Image width  in pels (must be a multiple of 16)
  input->img_height = 0;          //Image height in pels (must be a multiple of 8 )
  strcpy(input->TraceFile, "trace_enc.txt");
  strcpy(input->ReconFile, "test_rec.yuv");
#ifdef _OUTPUT_DEC_IMG_
  strcpy(input->DecRecFile, "org_dec.yuv");
#endif
  //input->outfile
  strcpy(input->outfile, input->infile);
  strcat(input->outfile, ".avs");

  //##########################################################################################
  //# Encoder Control
  //##########################################################################################
  input->thread_num = 1;      //Number of thread(s)
  input->GopLength  = 18;     //Length of one GOP, must be a multiple of (B frame number + 1)
  input->qp0 = 36;            //QP of frame I (0  ~ 63)
  input->qpN = 36;            //QP of frame P (0  ~ 63)
  input->qpB = 36;            //QP of frame B (0  ~ 63)
  input->hadamard = 0;
  input->search_range = 16;        //Search range
  input->no_multpred  = 2;          //Number of reference frames (1, 2)
  input->InterSearch16x16 = 1;
  input->InterSearch16x8 = 1;
  input->InterSearch8x16 = 1;
  input->InterSearch8x8 = 1;
  input->fr = 0;
  //input->intra_period=1;
  //##########################################################################################
  //# B Frames
  //##########################################################################################
  input->successive_Bframe = 1;      //Number of B frames (0, 1, 2)
  //##########################################################################################
  //# RD Optimization
  //##########################################################################################
  input->rdopt = 1;            //RDO switch (0: OFF, 1: ON)

  //##########################################################################################
  //# Additional Stuff
  //#########################################################################################
  input->progressive_frame = 1;      //is progressive sequence? (1: progressive, 0: interlace)
  input->InterlaceCodingOption = 0;

  //##########################################################################################
  //# Loop filter parameters
  //##########################################################################################
  input->loop_filter_disable = 1;
  input->loop_filter_parameter_flag = 1;
  input->alpha_c_offset = 6;
  input->beta_offset    = 6;

  //##########################################################################################
  //# Slice parameters
  //##########################################################################################
  input->slice_row_nr = 0;

  //##########################################################################################
  //# Weighting Prediction parameters
  //##########################################################################################
  input->picture_weighting_flag = 0;

  //##########################################################################################
  //#frame rate
  //###########################################################################################
  input->frame_rate_code = 2;        //Frame rate code [1)24000/1001, 2)24, 3)25, 4)30000/1001, 5)30, 6)50, 7)60000/1001, 8)60]

  //###########################################################################################
  //#chroma format parameter
  //###########################################################################################
  input->chroma_format = 1;        //Chroma format (1: 4:2:0,     2: 4:2:2)

  //########################################################################################
  //#Rate control
  //########################################################################################
  input->RCEnable = 0;          //Flag of rate control (0: off, 1: one-pass, 2: multi-pass?)
  input->bit_rate = 1000000;        //Bitrate (in bps)
  input->SeinitialQP = 36;        //Initial QP for the first I frame(available with one-pass rate control)
  input->basicunit = 1;
  input->channel_type = 0;        //   # type of channel( 1=time varying channel; 0=Constant channel)


  if (argc >= 4)
  {
    if (*argv[2] != '-')
    {
      strcpy(input->outfile, argv[2]);
    }

    for ( i = 0, err = 0; ++i < argc  &&  !err; )
    {
      char*  token;
      char*  nextArg;
      int    argUsed;
      token = argv[i];
      if ( *token++ == '-' )
      {
        argUsed = 0;
        nextArg = i+1 < argc  ?  argv[i+1]  :  "";

        switch (*token++)
        {
        case 'w':
          if (*token == 0)
          {
            input->img_width = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong Width param!\r\n");
            err = 1;
            return -1;
          }

        case 'h':
          if (*token == 0)
          {
            input->img_height = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong Height param!\r\n");
            err = 1;
            return -1;
          }

        case 'c':
          if (*token == 0)
          {
            input->chroma_format = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong ChromaFormat param!\r\n");
            err = 1;
            return -1;
          }

        case 'p':
          if (*token == 0)
          {
            input->progressive_frame = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong ProgressiveSequence param!\r\n");
            err = 1;
            return -1;
          }

        case 'g':
          if (*token == 0)
          {
            input->GopLength = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong GOPLength param!\r\n");
            err = 1;
            return -1;
          }

        case 'r':
          if (*token == 0)
          {
            input->no_multpred = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong RefNumber param!\r\n");
            err = 1;
            return -1;
          }

        case 'b':
          if (*token == 0)
          {
            input->successive_Bframe = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong BFrameNumber param!\r\n");
            err = 1;
            return -1;
          }
        case 'q':
          if (*token == 0)
          {
            input->qp0 = input->qpN = input->qpB = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong Flag of QP param!\r\n");
            err = 1;
            return -1;
          }

        case 'a':
          if (*token == 0)
          {
            input->RCEnable = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong Flag of Rate Control param!\r\n");
            err = 1;
            return -1;
          }

        case 'f':
          if (*token == 0)
          {
            input->frame_rate_code = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong FrameRate param!\r\n");
            err = 1;
            return -1;
          }

        case 's':
          if (*token == 0)
          {
            input->bit_rate = atoi(nextArg) * 1000;
            bUDF_ASFBitrate = true;
            break;
          }
          else
          {
            printf("Wrong Bitrate param!\r\n");
            err = 1;
            return -1;
          }

        case 'o':
          if (*token == 0)
          {
            input->rdopt = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong RDO param!\r\n");
            err = 1;
            return -1;
          }

        case 'm':
          if (*token == 0)
          {
            input->search_range = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong SearchRange param!\r\n");
            err = 1;
            return -1;
          }

        case 'i':
          if (*token == 0)
          {
            //pEncParam->iSeqHeader = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong Insert Sequence Header param!\r\n");
            err = 1;
            return -1;
          }

        case 't':
          if (*token == 0)
          {
            input->thread_num = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong ThreadNumber param!\r\n");
            err = 1;
            return -1;
          }

        case 'l':
          if (*token == 0)
          {
            input->loop_filter_disable = atoi(nextArg);
            break;
          }
          else
          {
            printf("Wrong ThreadNumber param!\r\n");
            err = 1;
            return -1;
          }

        case 'n':
          if (*token == 0)
          {
            input->fr = (float)atof(nextArg);
            bFramerateFlag = true;
            break;
          }
          else
          {
            printf("Wrong FrameRate!\r\n");
            err = 1;
            return -1;
          }
        default:
          printf("Wrong params!\r\n");
          err = 1;
          return -1;

        }
      }
    }
  }
  else
  {
    //output the help command
    printf("VideoTransCoder.exe\n");
    printf("will output the help message\n");
    printf("usage\n");
    printf("VideoTransCoder.exe input.m2v out.avs\n");
    printf("command options\n");
    printf("\t -w width\n");
    printf("\t -h height\n");
    printf("\t -c chroma_format\n");
    printf("\t -p progressive_frame (1(default):progressive, 0:interlace)\n");
    printf("\t -g gop_length (18:default)\n");
    printf("\t -r reference_frame_numer(2:default)\n");
    printf("\t -b b_frame_num(2:default)\n");
    printf("\t -q default_qp_for_i_b_p_frame(26:default)\n");
    printf("\t -a rate control enable flag(0:default, disable, 1:enable)\n");
    printf("\t -f frame_rate_code(2:default) Frame rate code [1)24000/1001, 2)24, 3)25, 4)30000/1001, 5)30, 6)50, 7)60000/1001, 8)60]\n");
    printf("\t -s bit_rate(1000000:default) in Kbps, s means speed\n");
    printf("\t -o rdo flag(1:default, enable, 0:disable)\n");
    printf("\t -m me search range(16:default)\n");
    printf("\t -n fps(default value is the decoded frame rate)\n");
    printf("\t -t thread number(1:default)\n");
    return -1;
  }
  return 0;
}