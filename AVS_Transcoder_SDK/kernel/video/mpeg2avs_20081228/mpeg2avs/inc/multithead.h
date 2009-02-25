#ifndef _MULTI_THREAD
#define _MULTI_THREAD

#include <Windows.h>

#ifndef FALSE
#define FALSE                   0
#endif

#ifndef TRUE
#define TRUE                    1
#endif

typedef void* PTRANSCODER_CREATE;
typedef struct node
  {
  struct node* pNext;
  void* pData;
  } NODE, *PNODE, *PLIST; //定义线程队列

// struct of encode parameters
typedef struct param_t
  {
  // version
  int       iProfileID;         // profile ID: 0x20
  int       iLevelID;           // level   ID: 0x20, 0x22, 0x40, 0x42

  // properties of input sequence
  int       iImageWidth;        // image width  in pels (input param, must be multiple of 16)
  int       iInputHeight;       // image height in pels (input param, must be multiple of 8 )
  int       iImageHeight;       // image height which is fixed to be a multiple of 16 pels [not input from config file]
  int       iChromaFormat;      // chroma format (1: 4:2:0, 2: 4:2:2)
  int       iProgressiveSeq;    // is progressive sequence? (1: progressive, 0: interlace)
  int       iTotalFrameNumber;  // number of frames to be encoded  (if -1, encoding forever until EOF)

  // GOP structure
  int       iGOPLength;         // length of one GOP, must be a multiple of (B frame number +1)
  int       iRefNumber;         // number of reference frames      (   1, 2)
  int       iBFrameNumber;      // number of B frames              (0, 1, 2)
  int       iCloseGOP;          // is close GOP?                   (1: yes, 0: no)

  // rate control & quantization parameters
  int       iTwoPassFlag;       // flag of rate control, 0: OFF,  1: one pass,  2: multi-pass
  int       iFrameRateCode;     // frame rate code  [1)24000/1001, 2)24, 3)25, 4)30000/1001, 5)30, 6)50, 7)60000/1001, 8)60]
  int       iBitrate;           // bitrate (in kbps)
  int       iInitialQP;         // initial QP for the first I frame  (available if iTwoPassFlag = 1)
  int       iMinQP;             // MIN value for QP (available if iTwoPassFlag = 1)
  int       iMaxQP;             // MAX value for QP (available if iTwoPassFlag = 1)
  int       iQPI;               // QP of frame I    (0  ~ 63)
  int       iQPP;               // QP of frame P    (0  ~ 63)
  int       iQPB;               // QP of frame B    (0  ~ 63)

  // other encoder controls
  int       iRDO;               // RDO switch       (0: OFF, 1: ON)
  int       iSearchRange;       // search range of motion search
  int       iBbvBufferSize;     // BBV buffer size  (in bytes)
  int       iSeqHeader;         // insert sequence header before each GOP? (1: yes, 0: no, only first GOP)
  int       iThreadNum;         // number of thread(s)

  // input & output file name
  char      Mpeg2FileName[256];// input original source file name (YUV)
  char      AVSRecFileName[256];// output reconstructed  file name (YUV)
  char      AVSBSFileName [256];// output bitstream file name
  } ENCODERPARAM, *PENCODERPARAM;
typedef struct tag_thread
  {
    BOOL      m_bRunning;         // this thread is still running?
    HANDLE    m_hThread;          // handle to the current thread
    HANDLE    m_hMutex;
    HANDLE    m_hEvent;
    ULONG     nThreadId;
    int       iThreadOrder;
    c_avs_enc *p_enc;
    //线程的参数
    char rec_file_name[100];
    char output_file_name[100];
    char config_file_name[100];
  }THREAD, *PTHREAD;
typedef struct tag_wthread
{
  BOOL      m_bRunning;         // this thread is still running?
  HANDLE    m_hThread;          // handle to the current thread
  HANDLE    m_hMutex;
  HANDLE    m_hEvent;
  ULONG     nThreadId;
  int       iThreadOrder;
  //c_avs_enc *p_enc;
  //void*     bitstream[4];
  int      nbit[4];
  int       order[8][2];
  //线程的参数
  char rec_file_name[100];
  char output_file_name[100];
  char config_file_name[100];
}WTHREAD, *PWTHREAD;
void mutithreads_test(int thread_num);
void create_multithread_transcoder(int thread_num, c_avs_enc **p_c_avs_enc);
DWORD WINAPI run_thread(LPVOID pArg);

#endif