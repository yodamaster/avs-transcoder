/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2004,  Audio Video Coding Standard, Part III
*
* This software module was originally developed by
*
* JungHoe Kim (kjh94@samsung.com), Samsung AIT
*
* edited by
*
* Lei Miao (win.miaolei@samsung.com), Samsung AIT
* Lei Miao, CBC Multi-channel extension, 2005-09-19
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

#include <stdio.h>
#include <math.h>
#include "sam_encode.h"
#include "av3enc.h"


int BitRate;

extern char	mesg[];
extern int WriteAASFHeader(int ch,int freq,int bitrate,int aasf_flag);		//for AASF
extern int WriteAATFHeader(int nr_of_blocks, int freqIdx, int nch,int bitrate, 
						   int aatf_buffer_fullness,  int aatf_flag);		//for AATF

extern int BitstreamOpen(char *fname);
extern void BitstreamClose(void);
extern void output_byte(long byte,int len);
extern void start_outputing_bits();
extern void done_outputing_bits();
extern void FlushBuffer(void);
extern void FlushBufferWithLength(void);
extern int GetBitstreamSize(void);
extern int ByteAlign(void);
//~bs
extern void BSHCInit(FILE *name);
extern void BSHCClose(void);
extern int EncodeBSHCQuant(int model,int upper,int cur);
extern int EncodeBSHCSi(int siVal);
extern int EncodeBSHCSf(int sfVal);
extern int EncodeBSHCBin(int ms);
extern void EncodeBSHCStart(void);
extern void EncodeBSHCEnd(void);
extern void EncodeBSHCPutBits(int val,int len);
extern int EncodeBSHCGetSize();
extern void EncodeBSHCHeader(int ch,int freq);
extern int EncodeBSHCByteAlign(void);
extern int BSHCModelSelect(int bshcModel,int bpl);
extern void EncodeBSHCFlush(void);
extern int WriteAASFHeader(int ch,int freq,int bitrate,int aasf_flag);		//for AASF




extern void output_byte(long byte,int len);
extern void start_outputing_bits();
extern void done_outputing_bits();
extern void FlushBuffer(void);
extern void FlushBufferWithLength(void);
extern int BitstreamOpen(FILE *fname);
extern void BitstreamClose(void);
extern int GetBitstreamSize(void);
extern int ByteAlign(void);


extern void FlushBufferInit();
extern void FlushBufferWithLength_ext(void);
extern void FlushFrame(int fill_bits);

extern void FlushAATFMCFrame();
extern void FlushBufferInit();
// 2006-01-20 xusen
extern void ScanAATFFrame(void);

//for AASF
extern void FlushAASFHeaderBuffer();

//for AATF        
extern void Register_buffer(void);
extern int	EncodeBSHCGetAATFSize();
extern void Restore_buffer(void);
extern void FlushAATFHeaderBuffer();

extern void sam_scale_bits_init(int fsidx);

extern void BSHCFindModel(
						  int nch,
						  int	signal_type,
						  int	num_window_groups,
						  int window_group_length[],
						  int	*sample[][8],
						  int	*scalefactors[][8],
						  int	maxSfb,
						  short	*swb_offset[],
						  int top_layer,
						  int     *bshc_model[2][8]
);

extern int sam_encode_cbc(
						  ChannelInfo channelInfo[],				   
						  int outputformat,
						  int target_bitrate,
						  int	signal_type,
						  int	windowShape,
						  int	*sample[][8],
						  int *scalefac[][8],
						  int	maxSfb,
						  int minSfb,
						  int	num_window_groups,
						  int window_group_length[],
						  short	*swb_offset[],
						  int	*model_index[][8],
						  int	stereo_mode,
						  int	*stereo_info[],
						  int abits,
						  int nch,
						  int wflag,
						  int isExtended,
						  int channelIdx,
						  int mc_present,
						  int extended_bytes,
						  int *fill_bits);

extern int aatf_raw_data_block_error_check();

void init_cbc(int fsidx, int bitrate, FILE *outFile)
{

	BitRate = bitrate;

	sam_scale_bits_init(fsidx); // zhanjie 0708

	/* initialize huffman table */
    BSHCInit(outFile);
}

int sam_encode(ChannelInfo channelInfo[],int outputformat, int bitrate, int signal_type, int windowShape,
				short *swb_offset[],int *scalefac[][8],	int	num_window_groups,
				int	window_group_length[], int *quant[][8], int maxSfb,	int	stereo_mode,
				int	*stereo_info[],	int abits, int nch,
				int wflag,
				int isExtended,
				int channelIdx,
				int mc_present,
				int extended_bytes,
				int *fill_bits)
{
	int		ch;
	int		g, s;
	int		frame_length;
	int		frame_length1;
	int		minSfb=0;
	int		end_cband[2][8];
	int		ModelIndex[2][128];
	int     *bshc_model[2][8];
	int		band_snf[2][128];
	int		*cband_qbit[2][8];
	double  cband_codelen[2][128];
	double  *cband_spec_codelen[2][8];

	bitrate *= nch;
	if(maxSfb==3)
		/* LFE coding bitrate */
		bitrate = 16000;

	for(ch = 0; ch < nch; ch++) {
		for(g = 0; g < 128; g++)
			ModelIndex[ch][g] = 0;
	}

	for(ch = 0; ch < nch; ch++) {
		if(signal_type == TRANSIENT_TYPE) {
			s = 0;
			for(g = 0; g < num_window_groups; g++) { 
				end_cband[ch][g] = (swb_offset[g][maxSfb]+31)/32;
			}
		} else {
			end_cband[ch][0] = (swb_offset[0][maxSfb]+31)/32;
		}

		bshc_model[ch][0] = &(ModelIndex[ch][0]);
		for (g = 1; g < num_window_groups; g++) {
			bshc_model[ch][g] = bshc_model[ch][g-1] + end_cband[ch][g-1];
		}
	}
	/* 2005-12-1 xuhengyu */
//	BSHCFindModel(nch,signal_type,num_window_groups,
//		window_group_length,quant,scalefac,maxSfb,swb_offset,
//					bitrate/nch/1000-16,bshc_model);
	BSHCFindModel(nch,signal_type,num_window_groups,
		window_group_length,quant,scalefac,maxSfb,swb_offset,
		(bitrate/nch/1000-16) > 63 ? 63 : (bitrate/nch/1000-16),
		bshc_model);

	/********** B I T S T R E A M   M A K I N G *********/

	frame_length1 = sam_encode_cbc(channelInfo,
		outputformat, bitrate, signal_type, windowShape, quant,
			scalefac, maxSfb, minSfb, num_window_groups, 
			window_group_length, swb_offset,
			bshc_model, stereo_mode, stereo_info, 
			abits, nch, 0, isExtended, channelIdx, 
			mc_present, extended_bytes, fill_bits);


	frame_length = sam_encode_cbc(channelInfo,
		outputformat, bitrate, signal_type, windowShape, quant,
			scalefac, maxSfb, minSfb, num_window_groups,
			window_group_length, swb_offset,
			bshc_model, stereo_mode, stereo_info,
			frame_length1, nch, 1, isExtended, channelIdx, 
			mc_present, extended_bytes, fill_bits);

	if(abs(frame_length1-frame_length)>0)
	{
		fprintf(stderr,"\n\n\t\t\terror bitstream %d\t%d",frame_length1,frame_length);
	}
	return frame_length;
}

#define NR_OF_BLOCKS 1
#define AATF_FLAG    0
#define AATF_BUFFER_FULNESS 0x7FF

int cbc(
		ChannelInfo channelInfo[],
        int	nch,
		int sampleRateIdx,	
		int outputformat,
		int signal_type,
		int	windowShape,
		int num_window_groups,
		int	window_group_length[],
		short *swb_offset[],
		int	*scalefac[][8],
		int	*quant[][8],
		int	maxSfb,
		int stereo_mode,
		int	*stereo_info[],
		int frame_length,
		int wflag,
		int isExtended,
		int channelIdx,
		int numChannels,
		int mc_present,
		int extended_bytes,		
		AV3EncFramePtr hEncoder)
{
	int	frameSize;
	static int frameStart=1;

	int fill_bits_cnt = 0;	
	static int iTotalFrameSize = 0;
	float tmp = 0.0;
	
	//for AASFHeader
	int aasf_flag = 0;

	//for AATFHeader
	static int nr_of_blocks_remained = NR_OF_BLOCKS - 1;
	static int bit_count = 0;	
	int 	aatf_flag = AATF_FLAG;
	int		aatf_buffer_fullness = AATF_BUFFER_FULNESS;

	if(outputformat != 2 || nr_of_blocks_remained == (NR_OF_BLOCKS -1)) 
		EncodeBSHCStart();


/*******************************************************/
/*         C B C   F R A M E   H E A D E R           */
/*******************************************************/
	
	if(outputformat == 0 && frameStart) {
		frameStart = 0;
		EncodeBSHCHeader(numChannels,sampleRateIdx);
	}
	
	if(outputformat == 1 && frameStart) {
		frameStart = 0;
		if(!WriteAASFHeader(numChannels,sampleRateIdx,BitRate,aasf_flag))
		{
			fprintf(stderr,"\nerror in write av3 audio header.\n");
			return 0;
		}
	}

	if( (numChannels<=2) || !isExtended) {
		if(outputformat == 2 && nr_of_blocks_remained == NR_OF_BLOCKS - 1)
			bit_count += WriteAATFHeader( nr_of_blocks_remained, sampleRateIdx,	numChannels,
									  BitRate, aatf_buffer_fullness, aatf_flag);
		else
			nr_of_blocks_remained--;
	}

	// buffer control - xun 2006-4-16
	do {
		frameSize = sam_encode(channelInfo, outputformat,
								BitRate, signal_type, 
								windowShape, swb_offset, 
								scalefac, num_window_groups,
								window_group_length, quant,
								maxSfb, stereo_mode, 
								stereo_info, frame_length, 
								nch, wflag, 
								isExtended, channelIdx, 
								mc_present, extended_bytes,
								&fill_bits_cnt);
	
		// 20060116 Miao
		frameSize += fill_bits_cnt;
		
		iTotalFrameSize += frameSize;
		if (!isExtended)			
		{		
			float fltBufferStatus = hEncoder->fltBufferFullness;
			float fltAvgBitPerFrm = hEncoder->fltAvgBitPerFrm;

			tmp = fltBufferStatus - (float)iTotalFrameSize + fltAvgBitPerFrm;
			tmp -= (float)hEncoder->iBufferSize;
			if (tmp > 0.0)	// underflow
			{
				fill_bits_cnt = (int)tmp + 1;
				iTotalFrameSize -= frameSize;
				//printf("\n buffer underflow\n");
			}
		}
		else
			tmp = 0.0;
	}while(tmp > 0.0);

	// buffer control - xun 2006-4-16	
	if (hEncoder->fltBufferFullness - (float)iTotalFrameSize< 0.0)
	{
		iTotalFrameSize = 0;
		//printf("\n buffer overflow \n");
		return -1;		// overflow
	}

	if (!isExtended)
		iTotalFrameSize = 0;

	if (outputformat == 2) {
		if(numChannels<=2) {
			bit_count += frameSize;
			bit_count += aatf_raw_data_block_error_check();	//->not implemented
			
			if(!nr_of_blocks_remained){

				// write the real frame_length
				FlushAATFHeaderBuffer();		// xusen: write one aat frame in file

				bit_count =	0;
				nr_of_blocks_remained = NR_OF_BLOCKS - 1;
 			}
		}
		else if(isExtended)
			FlushBufferWithLength_ext();		// xusen: move side-channel data in OutputBuffer to OutputBuffer_ext
		else
		{
			FlushAATFMCFrame();
		}
	}
	else {
		if(numChannels<=2)
			EncodeBSHCFlush();
		else if(isExtended)
			FlushBufferWithLength_ext();
		else {
			FlushFrame(fill_bits_cnt);
		}
	}

	return frameSize;
}