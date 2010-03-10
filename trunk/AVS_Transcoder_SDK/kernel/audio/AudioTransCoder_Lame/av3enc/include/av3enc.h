/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2004,  Audio Video Coding Standard, Part III
*
* This software module was originally developed by
* edited by
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the users without any
* license fee or royalty on an "as is" basis. The AVS disclaims
* any and all warranties, whether express, implied, or statutory,
* including any implied warranties of merchantability or of fitness
* for a particular purpose. In no event shall the contributors or 
* the AVS be liable for any incidental, punitive, or consequential
* damages of any kind whatsoever arising from the use of this program.
*
* This disclaimer of warranty extends to the user of this program
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The AVS does not represent or warrant that the program furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of AVS, including shareware, may be
* subject to royalty fees to patent holders. Information regarding
* the AVS patent policy is available from the AVS Web site at
* http://www.avs.org.cn
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
************************************************************************
*/

#ifndef _av3enc_H_
#define _av3enc_H_

#include <libtsp.h>
#include <libtsp/AFpar.h>


/* AVS ID's */
#define AVS1 1

#define SHORTCTL_NORMAL    0
#define SHORTCTL_NOSHORT   1
#define SHORTCTL_NOLONG    2

#define min(x,y) ((x) > (y) ? (y) : (x))
#define max(x,y) ((x) > (y) ? (x) : (y))

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef float psyfloat;

#define MAX_CHANNELS 8
#define MAX_BANDS	128
#define FRAME_LEN 1024

#define SINE_WINDOW 0
#define NSFB_LONG  51
#define NSFB_SHORT 15
#define MAX_SHORT_WINDOWS 8
#define BLOCK_LEN_LONG 1024
#define BLOCK_LEN_SHORT 128
#define MAX_SCFAC_BANDS ((NSFB_SHORT+1)*MAX_SHORT_WINDOWS)

#define FLP_MAX_ORDER			    20
#define FLP_GAIN_THRESH				2.7
#define MAX_SAMPLING_RATES		    12

#define MAX_FILL_BITS				2160	// xun
#define FILL_DFT					0		// xun
////////////////////////////////////////////////

enum WINDOW_TYPE {
    ONLY_LONG_WINDOW,
    LONG_SHORT_WINDOW,
    ONLY_SHORT_WINDOW,
    SHORT_LONG_WINDOW
};


enum stream_format {
  RAW_STREAM = 0,
  SAVE_STREAM = 1,
  TRANSFORM_STREAM = 2,
};

typedef struct
{
   AFILE *f;
   int channels;
   float samplerate;
   int samples;
} pcmfile_t;

typedef struct {                           
    int     bitrate;                   
    int     cutoff;
    int     stereo;                    
    int     wincontrol;
	int     format;
    char*   outFile;                 
    char*   inFile;                      
    int     showHelp;
	/* LFE coding */
	int		LFE_3SFB;
	int		cbr;
} cmd_params;

typedef struct
{
	int		nAllUsedFLP;						/* use FLPVQ in all frame */
	int		nFlpPresent;						/* use FLPVQ in a frame  */
	int		nPredOrder;				 		    /* predictor order	      */
	double	dPredGain;							/* predictor gain  		*/
	int		nPredDirection;						/* predictor direction		*/
	int		nPredBandLength;							/* the length of predictor scalefactor band */
	double	adPredCoeffs[FLP_MAX_ORDER + 1];			/* predictor coefficient */
	double	adReflectCoeffs[FLP_MAX_ORDER + 1];			/* reflection coefficient*/
	int		anQuantedReflectIndex[FLP_MAX_ORDER + 1];	/* reflection coefficient after quantization */
	
	int		nFlpMinBandNum;				/* the index of begin scalefactor band for smooth-change signal*/
	int		nFlpMaxBandNum;				/* the index of end scalefactor band for smooth-change signal*/
	int		nFlpMaxOrderNum;			/* the maximum predictor order*/
	
	double	adLineSpectrumFreq[FLP_MAX_ORDER];			/* save line spectrum frequency*/
	int		nVQCode1Index;				/* code book 1 for vector quantization */
	int		nVQCode2Index;				/* code book 2 for vector quantization */
	int		nVQCode3Index;				/* code book 3 for vector quantization */
	int		nVQCode4Index;				/* code book 4 for vector quantization */
} FLPVQInfo;

typedef struct
{
	int nAllUsedMR;						
	int nMRSwitchOn;
	int anLayoutOfSubband;
    int nCBNum;               /////// int
    int *anCBWidth;
} MRInfo;

typedef struct {
    int present;
    int ps_used[MAX_SCFAC_BANDS];
}PSInfo;

typedef struct {
    int present;
    int ch_is_left;
    int paired_ch;
    int common_window;	
	FLPVQInfo	flpInfo;			/* FLPVQ structure */
	MRInfo mrInfo;					/* MR structure */
    PSInfo psInfo;
} ChannelInfo;

typedef struct {
    int window_shape;
    int block_type;
    int desired_block_type;
	enum SIGNAL_TYPE signal_type;

    int global_gain;
    int scale_factor[MAX_SCFAC_BANDS];

    int num_window_groups;
    int window_group_length[8];
    int max_sfb;
    int nr_of_sfb;

    int swb;

    int sfb_offset[128];
    int lastx;
    double avgenrg;

	double *requantFreq;
    int *quantFreq;
} CoderInfo;

typedef struct {
  unsigned long sampling_rate;  /* the following entries are for this sampling rate */
  int num_cb_long;
  int num_cb_short;
  int cb_width_long[NSFB_LONG];
  int cb_width_short[NSFB_SHORT];
} SR_INFO;

typedef struct AV3EncCfg
{
    /* config version */
    int version;

    /* copyright string */
    char *copyright;

    /* AVS version */
    unsigned int AVSVersion;

	/* samplerate of AV3 file */
	unsigned int sampleRate;                  ////////// long
	unsigned int sampleRateIdx;               ////////// int

	/* number of channels in AV3 file */
    unsigned int numChannels;                 ////////// int 

    /* Allow Square Polar Stereo Coding */
    unsigned int allowSPSC;

	/* Allow flpvq */
	unsigned int allowFLPVQ;

    /* bitrate / channel of AV3 file */
    unsigned long bitRate;

    /* AV3 file frequency bandwidth */
    unsigned int bandWidth;

    /* Bitstream output format (0 = Raw; 1 = ADTS) */
    unsigned int outputFormat;

    /* block type enforcing (SHORTCTL_NORMAL/SHORTCTL_NOSHORT/SHORTCTL_NOLONG) */
    int shortctl;

	int top_layer;

	int desbits;
	
	int channel_map[8];

	/* LFE coding */
	int isLFE;
	int LFE_3SFB;
} AV3EncCfg, *AV3EncCfgPtr;

enum SIGNAL_TYPE { 				
	STATIONARY_TYPE, 		
	TRANSIENT_TYPE			
};
///////////////////////////////////////////////
typedef struct {
    unsigned int usedBytes;

    /* frame number */
    unsigned int frameNum;
    unsigned int flushFrame;

	/* output bits difference in average bitrate mode */
    int bitDiff;

	/* Psydchoacoustics data */
	int psy_block_type[MAX_CHANNELS];

    /* sample buffers of current next and next next frame*/
    double *sampleBuff[MAX_CHANNELS];
    double *nextSampleBuff[MAX_CHANNELS];
   
    /* Filterbank buffers */
    double *freqBuff[MAX_CHANNELS];
    double *overlapBuff[MAX_CHANNELS];

	/* Scalefactorband data */
    SR_INFO *srInfo;

    /* Channel and Coder data for all channels */
    CoderInfo coderInfo[MAX_CHANNELS];
    ChannelInfo channelInfo[MAX_CHANNELS];

    /* Configuration data */
    AV3EncCfg config;

    /* quantizer specific config */
    double quality;

	/* the number of empty bits in the buffer */
    float fltBufferFullness;
	int	  iBufferSize;
	float fltAvgBitPerFrm;
} AV3EncFrame, *AV3EncFramePtr;

int  AV3EncSetConfiguration (AV3EncFramePtr hEncoder, AV3EncCfgPtr config);

AV3EncFramePtr AV3EncOpen(unsigned long sampleRate,
                                  unsigned int numChannels,
                                  unsigned long *inputSamples,
                                  unsigned long *maxOutputBytes);

int AV3EncEncode(AV3EncFramePtr hEncoder,
                          int *inputBuffer,
                          unsigned int samplesInput,
                          unsigned char *outputBuffer);

int AV3EncClose(AV3EncFramePtr hEncoder);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif