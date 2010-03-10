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
*
* Lei Miao, fill data embedded, 2005-03-16
*
* Lei Miao, coding band bug fixed, 2005-06-15
* Lei Miao, flexible frame length/MC_present, 2005-06-15
* Lei Miao, CBC Multi-channel extension, 2005-09-19
* Lei Miao, CBC codebook reduction, 2005-12-22
* Xun Li  , fill bits for buffer control, 2006-4-16
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

#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "av3enc.h"

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


extern void Register_buffer(void);

extern void Restore_buffer(void);

extern int EncodeBSHCGetAATFSize(void);

int encode_cband_si(
					int	*model_index[][8],
					int	start_cband[8],
					int end_cband[8],
					int g,
					int nch);

int encode_scfband_si(
					  int	*scf[][8],
					  int	start_sfb[8],
					  int end_sfb[8],
					  int stereo_mode,
					  int *stereo_info[8],
					  int stereo_si_coded[8][MAX_SCFAC_BANDS],
					  int g,
					  int nch);

int encode_spectra(
				   int	*sample[][8],
				   int s_reg,
				   int e_reg,
				   int s_freq[8],
				   int e_freq[8],
				   int	min_bpl,
				   int available_len,
				   int *model_index[][8],
				   int *cur_bpl[][8],
				   int *coded_samp_bit[][8],
				   int *sign_coded[][8],
				   int *enc_vec[][8],
				   int	nch);

#define min(x,y) ((x) > (y) ? (y) : (x))
#define max(x,y) ((x) > (y) ? (x) : (y))

double  fs_tbl[12]=
{96000, 88200, 64000, /* zhanjie 0708 */ 48000., 44100., 32000., 24000., 22050, 16000, 12000, 11025, 8000};

int channel_mapping_51[3][2]=
{{2, 0}, {4, 2}, {3, 3}};

static int sampling_rate;

void init_layer_variables (
    int signal_type,
    int num_window_groups,
    int window_group_length[],
    int maxSfb,
    int base_band,
    short	*swb_offset[],
    int top_layer,
    int *slayer_size,
    int layer_reg[],
    int layer_max_freq[],
    int layer_max_cband[],
    int layer_max_qband[]    
);

/* function prototypes */
void sam_scale_bits_init(int fsidx)
{
  sampling_rate = (int)fs_tbl[fsidx];
}

int get_base_band() 
{
	int base_band;
	switch(sampling_rate) {
		case 96000:
		case 88200:
			base_band = 256;			//12000Hz zhanjie 0708
			break;
		case 64000:
			base_band = 320;			//10000Hz zhanjie 0708
			break;
		case 48000:
		case 44100:
			base_band = 288 + 32+64;	//9000Hz
			break;
		case 32000:
			base_band = 416;		// 6500 Hz
			break;
		case 24000:
		case 22050:
			base_band = 512;	//6000Hz
			break;
		case 16000:
			base_band = 640;	//5000Hz
			break;
		case 12000:
		case 11025:
			base_band = 768;	//4500Hz
			break;
		case  8000:
			base_band = 832;	//3250Hz
			break;
	}
	return base_band;
}

int enc_wflag = 0;

int sam_encode_cbc(
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
	int *fill_bits)
{
	int		i, ch, g, m, k;
	int     top_layer;
	int		sfb, qband, cband;
	int		base_snf;
	int		base_band;
	int		layer_reg[100];       /* Array Size : 64 -> 100  shpark 2000.04.20 */  
	int		layer_max_freq[100];  /* Array Size : 64 -> 100  shpark 2000.04.20 */  
	int		layer_max_cband[100]; /* Array Size : 64 -> 100  shpark 2000.04.20 */  
	int		layer_max_qband[100]; /* Array Size : 64 -> 100  shpark 2000.04.20 */  
	//int		layer_bit_flush[100];    /* Array Size : 64 -> 100  shpark 2000.04.20 */  
	int		layer_buf_offset[100];   /* Array Size : 64 -> 100  shpark 2000.04.20 */  
	int		layer_si_maxlen[100];    /* Array Size : 64 -> 100  shpark 2000.04.20 */  
	int		layer_extra_len[100];    /* Array Size : 64 -> 100  shpark 2000.04.20 */  
	int		layer_cw_offset[100];    /* Array Size : 64 -> 100  shpark 2000.04.20 */  

	int		start_freq[8];    /* Array Size : 64 -> 8  shpark 2000.04.20 */
	int		start_qband[8];   /* Array Size : 64 -> 8  shpark 2000.04.20 */
	int		start_cband[8];   /* Array Size : 64 -> 8  shpark 2000.04.20 */
	int		end_freq[8];      /* Array Size : 64 -> 8  shpark 2000.04.20 */
	int		end_cband[8];     /* Array Size : 64 -> 8  shpark 2000.04.20 */
	int		end_qband[8];     /* Array Size : 64 -> 8  shpark 2000.04.20 */
	int		s_freq[8];        /* Array Size : 64 -> 8  shpark 2000.04.20 */
	int		e_freq[8];        /* Array Size : 64 -> 8  shpark 2000.04.20 */

	int		cw_len;
	int		layer;
	int		slayer_size;
	int		si_offset;
	//int		est_total_len;
	int		available_len;
	//int		used_bits;
	int		cur_snf_buf[2][1024];
	int		*cur_snf[2][8];
	int		band_snf[2][8][32];
	int		sign_is_coded_buf[2][1024];
	int		*sign_is_coded[2][8];
	int		*coded_samp_bit[2][8];
	int		coded_samp_bit_buf[2][1024];
	int		*enc_vec[2][8];
	int		enc_vec_buf[2][256];
    int     stereo_si_coded[8][MAX_SCFAC_BANDS];
	int		code_len;
	
	int layerBit;
	int headerLen;
	int specBits=0;

    //for AATF 
	int bit_count0;
	int	bit_count1;

	/* fill element variables */
	int fill_enable=1;
	int fill_length=0;
	int fill_ele_number=0;	
	int fill_cnt[128];
	/* default fill */
	int fill_type=FILL_DFT;
	
	int enhance_FLCnt;
	int FL_fieldBits;
	int enhance_enable;
	int FL_Byte;
	int FL_Bit;
	int headerAlignBits, headerAlignBytes;

	int channel_configuration_index;
	int test=0;

	// 20060116 Miao
	int scalable_header_reuse=0;

	enc_wflag = wflag; 
	
	// 2006-4-16 xun
	*fill_bits = ((*fill_bits + 7) >> 3) << 3;

	/* ***************************************************** */
	/*                  Initialize variables                 */
	/* ***************************************************** */
	for(ch = 0; ch < nch; ch++) {
		int s;

		s = 0;
		for(g = 0; g < num_window_groups; g++) { 
			sign_is_coded[ch][g] = &(sign_is_coded_buf[ch][1024*s/8]);
			coded_samp_bit[ch][g] = &(coded_samp_bit_buf[ch][1024*s/8]);
			enc_vec[ch][g] = &(enc_vec_buf[ch][(1024*s/8)/4]);
			cur_snf[ch][g] = &(cur_snf_buf	[ch][1024*s/8]);
			s += window_group_length[g];

           for (i=0; i<MAX_SCFAC_BANDS; i++)
              stereo_si_coded[g][i] = 0;
		}

		for (i=0; i<1024; i++) {
			sign_is_coded_buf[ch][i] = 0;
			coded_samp_bit_buf[ch][i] = 0;
			cur_snf_buf[ch][i] = 0;
		}
		for (i=0; i<256; i++) {
			enc_vec_buf[ch][i] = 0;
		}
	}

	top_layer = ((target_bitrate/nch) / 1000) - 16;
	if (top_layer > 63) top_layer = 63; /* top_layer limit inserted by shpark 2000. 09. 25 */


	/******************************/
	/* Adjust necessary variables */
	/******************************/
	base_snf = 1;
    base_band = get_base_band();

	/* ************************************************
		Initialize
			layer_max_freq
			layer_max_cband
			layer_max_qband
			layer_bit_flush
	************************************************ */

   init_layer_variables (signal_type, num_window_groups, 
	   window_group_length, maxSfb,
       base_band, swb_offset, top_layer, &slayer_size, layer_reg,
       layer_max_freq,  layer_max_cband, layer_max_qband);

   /* maxSfb limitation */
   if (layer_max_qband[top_layer+slayer_size] < maxSfb)
		maxSfb = layer_max_qband[top_layer+slayer_size-1];

   for(g = 0; g < num_window_groups; g++) 
	  end_cband[g] = (swb_offset[g][maxSfb]+31)  / 32;

   if(signal_type == TRANSIENT_TYPE) {
	  for (layer = 0; layer < slayer_size; layer++) {
	 	 g = layer_reg[layer];
		 end_qband[g] = layer_max_qband[layer];
	  }
   }
   else
	  end_qband[0] = layer_max_qband[slayer_size-1];

	/* ***************************************************** */
	/*                  syntax of cbc_header()              */
	/* ***************************************************** */
	
	//HUFF start
   if(outputformat == 2)		
		Register_buffer();
   else
		EncodeBSHCStart();

   bit_count0 = EncodeBSHCGetSize();

	/* frame length */
	enhance_enable = 0;
	/* base frame length */
	if(abits == 20000)  // lixun 2008-7-28, max bits 6144 * 2=12288, so use 20000 as esmitated bits number
	{
		EncodeBSHCPutBits(0,8);
	}
	else
	{	
		enhance_enable = 1;
		enhance_FLCnt = 0;
		FL_Byte = abits >> 3;
		FL_Byte += extended_bytes;
		//bug fix: 20006-4-16 Xun
		if (!isExtended)
			FL_Byte += (*fill_bits>>3);

		/* 4 bytes: minimum header length */
		while( (FL_Byte-4) > ((1<<(7+enhance_FLCnt*3))-1) )
			enhance_FLCnt ++;
		FL_fieldBits = 7 + enhance_FLCnt*3;
		/* base frame length */
		EncodeBSHCPutBits((FL_Byte-4)>>(enhance_FLCnt*3),7);

		while(enhance_FLCnt) {
			EncodeBSHCPutBits(1, 1);
			EncodeBSHCPutBits(((FL_Byte-4)>>((enhance_FLCnt-1)*3))&0x7,3);
			enhance_FLCnt --;
		}
		/* frame length end flag */
		EncodeBSHCPutBits(0, 1);
	}

	/* MC_present */
	if(!isExtended) {
		EncodeBSHCPutBits(mc_present,1);
		if(mc_present) {
			/* FL/FR length */
			/* base frame length */
			if(abits == 20000)
			{
				EncodeBSHCPutBits(0,8);
			}
			else {
				enhance_enable = 1;
				enhance_FLCnt = 0;
				FL_Byte = abits >> 3;
				/* 4 bytes: minimum header length */
				while( (FL_Byte-4) > ((1<<(7+enhance_FLCnt*3))-1) )
					enhance_FLCnt ++;
				FL_fieldBits = 7 + enhance_FLCnt*3;
				/* base frame length */
				EncodeBSHCPutBits((FL_Byte-4)>>(enhance_FLCnt*3),7);
				while(enhance_FLCnt) {
					EncodeBSHCPutBits(1, 1);
					EncodeBSHCPutBits(((FL_Byte-4)>>((enhance_FLCnt-1)*3))&0x7,3);
					enhance_FLCnt --;
				}
				/* frame length end flag */
				EncodeBSHCPutBits(0, 1);
			}
		}
	}
	else {
		ch = 0;
		while(channel_mapping_51[ch][0] != channelIdx)
			ch ++;

		channel_configuration_index = channel_mapping_51[ch][1];

		/* channel configuration index: 4 */
		EncodeBSHCPutBits(channel_configuration_index,4);

		// 20060116 Miao
		/* scalable_header_reuse: 1 bit */
		EncodeBSHCPutBits(scalable_header_reuse, 1);

		/* reserved_bit : 1 */
		EncodeBSHCPutBits(0, 1);
	}
	
   	/* top layer index */
	EncodeBSHCPutBits(top_layer,6);

   	/* base snf */
	EncodeBSHCPutBits(base_snf-1,2);

   	/* base band */
	EncodeBSHCPutBits(base_band/32,5);

	/* ***************************************************** */
	/*               syntax of general_header()              */
	/* ***************************************************** */

	/* reserved_bit : 1 */
	EncodeBSHCPutBits(0,1);

	for(ch = 0; ch < nch; ch++) {
		EncodeBSHCPutBits(scalefac[ch][0][0],8);
	}

	EncodeBSHCPutBits((int)signal_type,1);

	EncodeBSHCPutBits((int)windowShape,1);

	if(signal_type == TRANSIENT_TYPE) {
		extern int siCodeMode;
		EncodeBSHCPutBits(maxSfb,4);

		for (i=1; i<window_group_length[0]; i++)
		{
			EncodeBSHCPutBits(1,1);
		}
		for(g = 1; g < num_window_groups; g++) {
			EncodeBSHCPutBits(0,1);
			for (i=1; i<window_group_length[g]; i++)
			{
				EncodeBSHCPutBits(1,1);
			}			
		}
		EncodeBSHCPutBits(siCodeMode,1);			
		
	} else {
		EncodeBSHCPutBits(maxSfb,6);
	}
	
	if(nch == 2) {
		EncodeBSHCPutBits(stereo_mode,2);
	}
	
	/* fill enable, 1 bit */
	// xun 2006-4-16 fill bits
	if (!isExtended) {
		if (*fill_bits == 0)
			fill_enable = 0;
		fill_length = *fill_bits >> 3;
		
		EncodeBSHCPutBits(fill_enable,1);

		if (*fill_bits > 0)
			fill_ele_number = (*fill_bits + MAX_FILL_BITS - 1)/MAX_FILL_BITS - 1;
		else
			fill_ele_number = 0;
		
		if(fill_enable)
		{			
			int tmp = *fill_bits;

			EncodeBSHCPutBits(fill_ele_number, 7);
			
			for (i = 0; i < fill_ele_number+1; i++)
			{
				fill_cnt[i] = min(tmp, MAX_FILL_BITS);
				tmp -= fill_cnt[i];
				fill_cnt[i] >>= 3;
								
				/* fill length field */
				if(fill_cnt[i] < 15)
				{
					EncodeBSHCPutBits(fill_cnt[i],4);
				}
				else
				{
					EncodeBSHCPutBits(15,4);
					EncodeBSHCPutBits(fill_cnt[i]-15,8);
				}
			}
		}
	}

	/* ***************************************************** */
	/*                  extension_coding_element()           */
	/* ***************************************************** */
	{		
		int chIndex = 0;
		FLPVQInfo *flpVQInfo;
		chIndex  = channelIdx;			
		for (ch = 0; ch < nch; ch++) {	
			flpVQInfo = &channelInfo[chIndex + ch].flpInfo; // wlei modified [2005-12-22]
			EncodeBSHCPutBits(flpVQInfo->nFlpPresent, 1);			
			if (flpVQInfo->nFlpPresent) {
				EncodeBSHCPutBits(flpVQInfo->nPredDirection, 1);
				EncodeBSHCPutBits(flpVQInfo->nVQCode1Index, 10);
				EncodeBSHCPutBits(flpVQInfo->nVQCode2Index, 10);
				EncodeBSHCPutBits(flpVQInfo->nVQCode3Index, 9);
				EncodeBSHCPutBits(flpVQInfo->nVQCode4Index, 8);
				chIndex ++;
			}			
		}	
	}

	

	if(!enhance_enable) {
		headerAlignBits = EncodeBSHCGetSize() & 0x7;
		if(headerAlignBits)
			headerAlignBits = 8 - headerAlignBits;
	}

	/*********** byte align *********/
 	EncodeBSHCByteAlign();

	if (outputformat == 2) 
		headerLen = EncodeBSHCGetAATFSize();
	else
		headerLen = EncodeBSHCGetSize();


	if(maxSfb == 0) {
		if(nch == 1) return headerLen;
		if(nch == 2 && maxSfb == 0) return headerLen;
	}
	
	for(ch = 0; ch < nch; ch++) {
		for(g = 0; g < num_window_groups; g++) {
			for(cband = 0; cband < end_cband[g]; cband++) {
				if(model_index[ch][g][cband]==0)
				{
					band_snf[ch][g][cband]=0;
					continue;
				}
				if(model_index[ch][g][cband] >= 9)
					band_snf[ch][g][cband] = model_index[ch][g][cband] -4;
				else
				{
					band_snf[ch][g][cband] = (model_index[ch][g][cband]+1)/2;
				}
			}
		}
	}

	/* ##################################################### */
	/* 	                CBC MAIN ROUTINE                     */
	/* ##################################################### */

	/* ***************************************************** */
	/*                  CBC_layer_stream()                   */
	/* ***************************************************** */
	for(g = 0; g < num_window_groups; g++)
		end_cband[g] = end_qband[g] = 0;

	for (layer=0; layer<(top_layer+slayer_size); layer++) {
		g = layer_reg[layer];
 	 	start_qband[g] = end_qband[g];
		start_cband[g] = end_cband[g];
		end_qband[g] = layer_max_qband[layer];
		end_cband[g] = layer_max_cband[layer];
	}

	/* to calculate core CBC data size */
	layer_buf_offset[0] = (double)16.*(abits-headerLen)
		/(double)(top_layer+16)/(double)slayer_size;
	for (layer = 1; layer < slayer_size-1; layer++)
	{
		layer_buf_offset[layer] = layer_buf_offset[layer-1] + 
			(double)16.*(abits-headerLen)/(top_layer+16)/slayer_size;
	}
	layer_buf_offset[slayer_size-1] = 
		(double)16.*(abits-headerLen)/(top_layer+16);
	layerBit = (abits-headerLen) - layer_buf_offset[slayer_size-1];
	for (layer = slayer_size; layer < top_layer+slayer_size-1; layer++)
	{
		layer_buf_offset[layer] = layer_buf_offset[layer-1] +
			layerBit/(top_layer);
	}
	layer_buf_offset[top_layer+slayer_size-1] = (abits-headerLen);

	for (layer = top_layer+slayer_size; layer < 100; layer++)
		layer_buf_offset[layer] = (abits-headerLen);

	for(i = 0; i < 100; i++) {
		layer_extra_len[i] = 0;
		layer_cw_offset[i] = 0;
	}

	//est_total_len = used_bits;

	for(g = 0; g < num_window_groups; g++) {
	   	start_qband[g] = 0;
		start_cband[g] = 0;
		start_freq[g] = 0;
	   	end_qband[g] = 0;
		end_cband[g] = 0;
		end_freq[g] = 0;
	}

	for(g = 0; g < num_window_groups; g++)
		s_freq[g] = 0;

	if (signal_type == TRANSIENT_TYPE) {
		for (layer = 0; layer < slayer_size; layer++) {
			g = layer_reg[layer];
			e_freq[g]  = layer_max_freq[layer];
		}
	}

	code_len = 0;
	
	available_len = 0;
	/* ##################################################### */
	/* 	                CBC MAIN LOOP                        */
	/* ##################################################### */
	for (layer = 0; layer < top_layer+slayer_size; layer++) {
		int min_snf;
		
		g = layer_reg[layer];

		start_freq[g]  = end_freq[g];
	   	start_qband[g] = end_qband[g];
		start_cband[g] = end_cband[g];

		end_qband[g] = layer_max_qband[layer];
		end_cband[g] = layer_max_cband[layer];
		end_freq[g]  = layer_max_freq[layer];

		for(ch = 0; ch < nch; ch++) {
			for(i = start_freq[g]; i < end_freq[g]; i++)  
				cur_snf[ch][g][i] = band_snf[ch][g][i/32];
		}

		if (signal_type == TRANSIENT_TYPE) {
			if(layer >= slayer_size) {
				g = layer_reg[layer];
				e_freq[g]  = layer_max_freq[layer];
			}
		} else {
			if (layer >= slayer_size) 
				e_freq[0] = layer_max_freq[layer];
			else
   				e_freq[0] = layer_max_freq[slayer_size-1];
		}
		if(layer==0)
			available_len = layer_buf_offset[0];
		else
			available_len += layer_buf_offset[layer]-
			layer_buf_offset[layer-1];
		
		//est_total_len += available_len;

		cw_len =  encode_cband_si(model_index, start_cband, end_cband, g, nch);
		code_len += cw_len;
		
		available_len -= cw_len;
		g = layer_reg[layer];
		
		/* side infomation : scalefactor */
		cw_len =  encode_scfband_si(scalefac, start_qband, 
			end_qband, stereo_mode, stereo_info,
						stereo_si_coded, g, nch);

		code_len += cw_len;
		available_len -= cw_len;
		
		min_snf = layer < slayer_size ? base_snf : 1;

		/* Bit-Sliced Spectral Data Coding */
		cw_len = encode_spectra(sample, g, g+1, start_freq, end_freq, 
				min_snf, available_len, model_index, cur_snf, coded_samp_bit, 
				sign_is_coded, enc_vec, nch);

		test += cw_len;
		
		code_len += cw_len;
		available_len -= cw_len;
		specBits+= cw_len;

		if(available_len > 0) {
			cw_len = encode_spectra(sample, 0, num_window_groups, s_freq, e_freq,
				1, available_len, model_index, cur_snf, coded_samp_bit, 
				sign_is_coded, enc_vec, nch);
			code_len += cw_len;
			available_len -= cw_len;
			specBits+= cw_len;

			test += cw_len;
		}
		//est_total_len -= available_len;

	}

	if(!enhance_enable) {
		int byte_incr = 0;
		int FL_Bit_cur = 0;
		int enhance_FLCnt_last = 0;
		// 20060116 Miao
		if(outputformat != 2)
			FL_Bit = EncodeBSHCGetSize();
		else
			FL_Bit = EncodeBSHCGetSize() - bit_count0;

		FL_Byte = (FL_Bit+7) >> 3;

		FL_Bit_cur = FL_Bit;
		while (1)
		{
			enhance_FLCnt = 0;
			FL_Byte = (FL_Bit_cur+7) >> 3;
			// 2006-4-16 Xun bug fixed
			if (!isExtended)
			{
				FL_Byte += extended_bytes;
				FL_Byte += fill_enable * fill_length;
			}
			/* 4 bytes: minimum header length */
			while((FL_Byte-4) > ((1<<(7+enhance_FLCnt*3))-1))
				enhance_FLCnt ++;
			
			// 2006-4-16 Xun
			if(!isExtended && mc_present) {
				int enhance_FLCnt2;
				FL_Byte -= extended_bytes;
				FL_Byte -= fill_enable * fill_length;

				enhance_FLCnt2 = 0;
				while((FL_Byte-4) > ((1<<(7+enhance_FLCnt2*3))-1))
					enhance_FLCnt2 ++;
				enhance_FLCnt += enhance_FLCnt2;
			}	
			if(enhance_FLCnt > enhance_FLCnt_last) {
				if(headerAlignBits < (enhance_FLCnt<<2)) {
					byte_incr = (enhance_FLCnt<<2) - headerAlignBits;					
					byte_incr = (byte_incr+7)>>3;
					// 2006-4-16 xun
					FL_Bit_cur = FL_Bit + (byte_incr<<3);
				}
				else
					break;
			}
			else
				break;
			enhance_FLCnt_last = enhance_FLCnt;
		}
		FL_Bit = FL_Bit_cur;
		FL_Byte = (FL_Bit+7) >> 3;	
		// 20060116 li xun
		for(i = 0; i < byte_incr; i ++)
			EncodeBSHCPutBits(0, 8);
	}

	EncodeBSHCByteAlign();

	/* fill data */
	// xun 2006-4-16 fill buffer
	if(!isExtended && fill_enable)
	{
		for (k = 0; k < fill_ele_number + 1; k++)
		{		
			EncodeBSHCPutBits(fill_type,4);
			switch(fill_type)
			{			
			case FILL_DFT:
				for(i = 0; i < 8*(fill_cnt[k]-1)+4; i ++)
					EncodeBSHCPutBits(1,1);
				break;
			default: 
				;
			}
		}
	}

	//EncodeBSHCEnd();
	// 20060116 Miao
	bit_count1 = EncodeBSHCGetSize()-*fill_bits;

	if(outputformat == 2 && !wflag)
		Restore_buffer();		
 	return (bit_count1 - bit_count0);

}

void init_layer_variables (
int signal_type,
int num_window_groups,
int window_group_length[],
int maxSfb,
int base_band,
short	*swb_offset[],
int top_layer,
int *slayer_size,
int layer_reg[],
int layer_max_freq[],
int layer_max_cband[],
int layer_max_qband[]
)
{
	int i, reg;
    int layer, g;
	int sfb, cband;
	int init_max_freq[8];
	int end_cband[8];
	int end_freq[8];
	int end_layer;

	/************************************************/
	/* calculate slayer_size                        */
	/************************************************/
	if(signal_type == TRANSIENT_TYPE) {
		*slayer_size = 0;
		for(g = 0; g < num_window_groups; g++) {
			init_max_freq[g] = (base_band * window_group_length[g]) / 32 * 4;
			switch(sampling_rate) {
			case 96000:
			case 88200:
			case 64000: /* zhanjie 0708 */
			case 48000:	case 44100:
				if((init_max_freq[g] % 32) >= 16)
					init_max_freq[g] = (init_max_freq[g]/32) * 32 + 20;
				else if((init_max_freq[g] % 32) >= 4)
					init_max_freq[g] = (init_max_freq[g]/32) * 32 + 8;
				break;
			case 32000:	case 24000:	case 22050:
				init_max_freq[g] = (init_max_freq[g]/16) * 16;
				break;
			case 16000:	case 12000:	case 11025:
				init_max_freq[g] = (init_max_freq[g]/32) * 32;
				break;
			case 8000:
				init_max_freq[g] = (init_max_freq[g]/64) * 64;
				break;
			}

			init_max_freq[g] = min(init_max_freq[g], swb_offset[g][maxSfb]);

			end_freq[g] = init_max_freq[g];
			end_cband[g] = (end_freq[g] + 31) / 32;
			*slayer_size += (init_max_freq[g]+31)/32;
		}
	} else { /* if(signal_type == TRANSIENT_TYPE) */
		g = 0;
		init_max_freq[g] = base_band;

		if (init_max_freq[g] > swb_offset[g][maxSfb])
			init_max_freq[g] = swb_offset[g][maxSfb];
		*slayer_size = (init_max_freq[g]+31)/32;
		end_freq[g] = init_max_freq[g];
		end_cband[g] = (end_freq[g] + 31) / 32;
	} /* if(signal_type == TRANSIENT_TYPE) ... else  */

	/************************************************/
	/* calculate layer_max_cband and layer_max_freq */
	/************************************************/
	for(g = 0, layer=0; g < num_window_groups; g++)
	for (cband = 1; cband <= end_cband[g]; layer++, cband++) {
		layer_reg[layer] = g;
		layer_max_freq[layer] = cband * 32;
		if (layer_max_freq[layer] > init_max_freq[g])
			layer_max_freq[layer] = init_max_freq[g];
		layer_max_cband[layer] = cband;
	}
	if(signal_type == TRANSIENT_TYPE) {	
	   layer = *slayer_size;
	   for(g = 0; g < num_window_groups; g++) {
	  	   for (i=0; i<window_group_length[g]; i++) { 
		  	   layer_reg[layer] = g;
			   layer++;
		   }
	   }
	   for (layer=*slayer_size+8; layer<100; layer++) 
	      layer_reg[layer] = layer_reg[layer-8];
    }
	else {
	   for (layer=*slayer_size; layer<100; layer++) 
	      layer_reg[layer] = 0;
	}
   
	for (layer = *slayer_size; layer <100; layer++) {
		reg = layer_reg[layer];
		switch(sampling_rate) {
			/* zhanjie 0708 begin */
		case 96000:
		case 88200:
			end_freq[reg] += 16;			
			break;
		case 64000:
			end_freq[reg] += 12;
			break;
			/* zhanjie 0708 ends */

	    case 48000: case 44100:
		    end_freq[reg] += 8;
		    if ((end_freq[reg]%32) != 0)
				end_freq[reg] += 4;
		    break;
		case 32000:	case 24000:	case 22050:
			end_freq[reg] += 16;
			break;
		case 16000:	case 12000:	case 11025:
			end_freq[reg] += 32;
		    break;
		default:
			end_freq[reg] += 64;
			break;
		}
		if (end_freq[reg] > swb_offset[reg][maxSfb] )
			end_freq[reg] = swb_offset[reg][maxSfb];
		//if( (top_layer+*slayer_size) <= layer)//
		//	end_freq[reg] = swb_offset[reg][maxSfb];


		layer_max_freq[layer] = end_freq[reg];
		layer_max_cband[layer] = (layer_max_freq[layer]+31)/32;
	}

	/*****************************/
	/* calculate layer_bit_flush */
	/*****************************/
#if 0
	for (layer = 0; layer < (top_layer+*slayer_size-1); layer++) {
		if (layer_max_cband[layer] != layer_max_cband[layer+1] ||
			layer_reg[layer] != layer_reg[layer+1] )
			layer_bit_flush[layer] = 1;
		else
			layer_bit_flush[layer] = 0;
	}
	for ( ; layer < 100; layer++)
		layer_bit_flush[layer] = 1;
#endif	// this part is conflit with standard spec. -- lixun 08/07/21

	/*****************************/
	/* calculate layer_max_qband */
	/*****************************/
	for (layer = 0; layer < 100; layer++) {
		end_layer = layer;
#if 0
		while ( !layer_bit_flush[end_layer] ) 
			end_layer++;
#endif		// this part is conflict with spec. -- lixun 08/07/21
		reg = layer_reg[layer];
		for (sfb=0; sfb<maxSfb; sfb++) {
			if (layer_max_freq[end_layer] <= swb_offset[reg][sfb+1]) {
				layer_max_qband[layer] = sfb+1;
				break;
			}
		}
		if (layer_max_qband[layer] > maxSfb)
			layer_max_qband[layer] = maxSfb;
	}
}

void get_end_band (
int signal_type,
int num_window_groups,
int window_group_length[],
int maxSfb,
short	*swb_offset[],
int top_layer,
int base_maxSfb[],
int end_sfb[],
int end_freq[],
int end_cband[]
)
{
	int i, reg;
    int layer, g;
	int sfb, cband;
	int init_max_freq[8];
    int slayer_size;
    int layer_reg[100];
    int base_band = get_base_band();

	/************************************************/
	/* calculate slayer_size                        */
	/************************************************/
	if(signal_type == TRANSIENT_TYPE) {
		slayer_size = 0;
		for(g = 0; g < num_window_groups; g++) {
			init_max_freq[g] = (base_band * window_group_length[g]) / 32 * 4;
			switch(sampling_rate) {
			case 96000: 
			case 88200:
			case 64000: // zhanjie 0708
			case 48000:	case 44100:
				if((init_max_freq[g] % 32) >= 16)
					init_max_freq[g] = (init_max_freq[g]/32) * 32 + 20;
				else if((init_max_freq[g] % 32) >= 4)
					init_max_freq[g] = (init_max_freq[g]/32) * 32 + 8;
				break;
			case 32000:	case 24000:	case 22050:
				init_max_freq[g] = (init_max_freq[g]/16) * 16;
				break;
			case 16000:	case 12000:	case 11025:
				init_max_freq[g] = (init_max_freq[g]/32) * 32;
				break;
			case 8000:
				init_max_freq[g] = (init_max_freq[g]/64) * 64;
				break;
			}

			end_freq[g] = init_max_freq[g];
			end_cband[g] = (end_freq[g] + 31) / 32;
			slayer_size += (init_max_freq[g]+31)/32;
		}
	} else { /* if(signal_type == TRANSIENT_TYPE) */
		g = 0;
		init_max_freq[g] = base_band;

		if (init_max_freq[g] > swb_offset[g][maxSfb])
			init_max_freq[g] = swb_offset[g][maxSfb];
		slayer_size = (init_max_freq[g]+31)/32;
		end_freq[g] = init_max_freq[g];
		end_cband[g] = (end_freq[g] + 31) / 32;
	} /* if(signal_type == TRANSIENT_TYPE) ... else  */

	/*****************************/
	/* calculate base_maxSfb */
	/*****************************/
	for (g = 0; g < num_window_groups; g++) {
		for (sfb=0; sfb<maxSfb; sfb++) {
			if (init_max_freq[g] <= swb_offset[g][sfb+1]) {
				base_maxSfb[g] = sfb+1;
				break;
			}
		}
		if (base_maxSfb[g] > maxSfb) base_maxSfb[g] = maxSfb;
	}

	/************************************************/
	/* calculate layer_max_cband and layer_max_freq */
	/************************************************/
	for(g = 0, layer=0; g < num_window_groups; g++)
	for (cband = 1; cband <= end_cband[g]; layer++, cband++) {
		layer_reg[layer] = g;
	}
	if(signal_type == TRANSIENT_TYPE) {
	   layer = slayer_size;
	   for(g = 0; g < num_window_groups; g++) {
	  	   for (i=0; i<window_group_length[g]; i++) { 
		  	   layer_reg[layer] = g;
			   layer++;
		   }
	   }
	   for (layer=slayer_size+8; layer<(top_layer+slayer_size); layer++) 
	      layer_reg[layer] = layer_reg[layer-8];
    }
	else {
	   for (layer=slayer_size; layer<(top_layer+slayer_size); layer++) 
	      layer_reg[layer] = 0;
	}
   
	for (layer = slayer_size; layer <(top_layer+slayer_size); layer++) {
		reg = layer_reg[layer];
		switch(sampling_rate) {
			/* zhanjie 0708 begin */
		case 96000:
		case 88200:
			end_freq[reg] += 16;			
			break;
		case 64000:
			end_freq[reg] += 12;
			break;
			/* zhanjie 0708 end */

	    case 48000: case 44100:
		    end_freq[reg] += 8;
		    if ((end_freq[reg]%32) != 0) end_freq[reg] += 4;
		    break;
		case 32000:	case 24000:	case 22050:
			end_freq[reg] += 16;
			break;
		case 16000:	case 12000:	case 11025:
			end_freq[reg] += 32;
		    break;
		default:
			end_freq[reg] += 64;
			break;
		}
	}


	/*****************************/
	/* calculate layer_max_qband */
	/*****************************/
	for (reg=0; reg < num_window_groups; reg++) {
	    if (end_freq[reg] > swb_offset[reg][maxSfb] )
		    end_freq[reg] = swb_offset[reg][maxSfb];
		end_cband[reg] = (end_freq[reg]+31)/32;
		for (sfb=0; sfb<maxSfb; sfb++) {
			if (end_freq[reg] <= swb_offset[reg][sfb+1]) {
				end_sfb[reg] = sfb+1;
				break;
			}
		}
		if (end_sfb[reg] > maxSfb)
			end_sfb[reg] = maxSfb;
	}
}

void get_base_maxSfb (
int signal_type,
int num_window_groups,
int window_group_length[],
int maxSfb,
short	*swb_offset[],
int base_maxSfb[]
)
{
    int g;
	int sfb;
    int base_band = get_base_band();
	int init_max_freq[8];

	/************************************************/
	/* calculate slayer_size                        */
	/************************************************/
	if(signal_type == TRANSIENT_TYPE) {
		for(g = 0; g < num_window_groups; g++) {
			init_max_freq[g] = (base_band * window_group_length[g]) / 32 * 4;
			switch(sampling_rate) {
			case 96000:
			case 88200:
			case 64000: /* zhanjie 0708 */
			case 48000:	case 44100:
				if((init_max_freq[g] % 32) >= 16)
					init_max_freq[g] = (init_max_freq[g]/32) * 32 + 20;
				else if((init_max_freq[g] % 32) >= 4)
					init_max_freq[g] = (init_max_freq[g]/32) * 32 + 8;
				break;
			case 32000:	case 24000:	case 22050:
				init_max_freq[g] = (init_max_freq[g]/16) * 16;
				break;
			case 16000:	case 12000:	case 11025:
				init_max_freq[g] = (init_max_freq[g]/32) * 32;
				break;
			case 8000:
				init_max_freq[g] = (init_max_freq[g]/64) * 64;
				break;
			}
		}
	} else { /* if(signal_type == TRANSIENT_TYPE) */
		init_max_freq[0] = base_band;
	} /* if(signal_type == TRANSIENT_TYPE) ... else  */

	/*****************************/
	/* calculate layer_max_qband */
	/*****************************/
	for (g = 0; g < num_window_groups; g++) {
		for (sfb=0; sfb<maxSfb; sfb++) {
			if (init_max_freq[g] <= swb_offset[g][sfb+1]) {
				base_maxSfb[g] = sfb+1;
				break;
			}
		}
		if (base_maxSfb[g] > maxSfb) base_maxSfb[g] = maxSfb;
	}
}

int cband_si_esc[9] = 
{-4, -3, -2, -1, 30, 31, 32, 33, 34};

/*--------------------------------------------------------------*/
/***************   coding band side infomation   ****************/
/*--------------------------------------------------------------*/
int encode_cband_si(
	int	*model_index[][8],
	int	start_cband[8],
	int end_cband[8],
	int g,
	int nch)
{
	int ch;
	int cband;
	int	m;
	int si_cw_len;
	int cband_model;
	
	static int prevSi[2][8];
	int t1,t2;
	extern int siCodeMode;

	si_cw_len = 0;
	for(ch = 0; ch < nch; ch++) {
		if(siCodeMode==0)
		{
			if(start_cband[g]==0 && g==0)
				prevSi[ch][0] = 0;
			for (cband = start_cband[g]; cband < end_cband[g]; cband++) {
				t2 = model_index[ch][g][cband];
				if(t2<-10000)
					t2 = t2;
		
				if((prevSi[ch][0]-t2+15)<0 || (prevSi[ch][0]-t2+15)>29 )
				{
					/* coding escape code */
					si_cw_len += EncodeBSHCSi(30);
					/* encoding 9 escape value with fixed 4 bit */
					t1 = prevSi[ch][0]-t2+15;
					m = 0;
					while(cband_si_esc[m] != t1)
						m ++;
					if(m & 0x8)
						si_cw_len += EncodeBSHCBin(1);
					else
						si_cw_len += EncodeBSHCBin(0);
					if(m & 0x4)
						si_cw_len += EncodeBSHCBin(1);
					else
						si_cw_len += EncodeBSHCBin(0);
					if(m & 0x2)
						si_cw_len += EncodeBSHCBin(1);
					else
						si_cw_len += EncodeBSHCBin(0);
					if(m & 0x1)
						si_cw_len += EncodeBSHCBin(1);
					else
						si_cw_len += EncodeBSHCBin(0);
				}
				else
					si_cw_len += EncodeBSHCSi(prevSi[ch][0]-t2+15);
				
				prevSi[ch][0] = t2;
			}
		}
		else
		{
			if(start_cband[g]==0)
				prevSi[ch][g] = 0;
			for (cband = start_cband[g]; cband < end_cband[g]; cband++) {
				t2 = model_index[ch][g][cband];
				
				if((prevSi[ch][g]-t2+15)<0 || (prevSi[ch][g]-t2+15)>29 )
				{
					/* coding escape code */
					si_cw_len += EncodeBSHCSi(30);
					/* encoding 9 escape value with fixed 4 bit */
					t1 = prevSi[ch][g]-t2+15;
					m = 0;
					while(cband_si_esc[m] != t1)
						m ++;
					if(m & 0x8)
						si_cw_len += EncodeBSHCBin(1);
					else
						si_cw_len += EncodeBSHCBin(0);
					if(m & 0x4)
						si_cw_len += EncodeBSHCBin(1);
					else
						si_cw_len += EncodeBSHCBin(0);
					if(m & 0x2)
						si_cw_len += EncodeBSHCBin(1);
					else
						si_cw_len += EncodeBSHCBin(0);
					if(m & 0x1)
						si_cw_len += EncodeBSHCBin(1);
					else
						si_cw_len += EncodeBSHCBin(0);
				}
				else
					si_cw_len += EncodeBSHCSi(prevSi[ch][g]-t2+15);
				
				prevSi[ch][g] = t2;
			}		
		}
	}

	return si_cw_len;
}
	
/*--------------------------------------------------------------*/
/***********  scalefactor band side infomation   ****************/
/*--------------------------------------------------------------*/
int encode_scfband_si(
	int	*scf[][8],
	int	start_sfb[8],
	int end_sfb[8],
	int stereo_mode,
	int *stereo_info[8],
	int stereo_si_coded[8][MAX_SCFAC_BANDS],
	int g,
	int nch)
{
	int ch;
	int	m;
	int sfb;
	int si_cw_len;
	
	static int prevSf[2][8];
	if(start_sfb[g]==0)
		for(ch = 0; ch < nch; ch++)
			prevSf[ch][g] = scf[ch][0][0];
	si_cw_len = 0;
	for(ch = 0; ch < nch; ch++) {
		
		for(sfb = start_sfb[g]; sfb < end_sfb[g]; sfb++) {
			if(stereo_si_coded[g][sfb] == 0) {
				if(stereo_mode != 2) {
					m = stereo_info[g][sfb];
					if(stereo_mode == 1) {						
						EncodeBSHCBin(m);
						si_cw_len++;
					}
				}
				stereo_si_coded[g][sfb] = 1;
			}
			
			if(!(g==0 && sfb==0))
			{
                                          
				if((scf[ch][g][sfb]-prevSf[ch][g]+64)<0 || (scf[ch][g][sfb]-prevSf[ch][g]+64)>127)
				{
					si_cw_len += EncodeBSHCSf(64);
				}
				else
                                           
				{
					if ((stereo_mode == 2 || (stereo_mode == 1 && stereo_info[g][sfb] == 1)) && ch!=0)
						si_cw_len += EncodeBSHCSf(scf[ch][g][sfb] - scf[0][g][sfb] + 64);
					else
						si_cw_len += EncodeBSHCSf(scf[ch][g][sfb]-prevSf[ch][g]+64);
					prevSf[ch][g] = scf[ch][g][sfb];
				}			
			}			
		}
	}

	return si_cw_len;
}

int cb_tab[16] = {0,1,2,3,4,3,3,5,6,3,3,5,3,5,5,5};

/*--------------------------------------------------------------*/
/***********     Bit-Sliced Spectral Data        ****************/
/*--------------------------------------------------------------*/
int encode_spectra(
	int	*sample[][8],
	int s_reg,
	int e_reg,
	int s_freq[8],
	int e_freq[8],
	int	min_bpl,
	int available_len,
	int *model_index[][8],
	int *cur_bpl[][8],
	int *coded_samp_bit[][8],
	int *sign_coded[][8],
	int *enc_vec[][8],
	int	nch)
{
	int ch, i, m, k, i4;
	int maxbpl;
	int	bpl;
	int	cband;
	int	maxhalf;
	int	cw_len;
	int	reg;

	
	int modelIdx;
	int j;

	int bit_mask[16] = { 0x0000, 0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040, 
						 0x0080, 0x0100, 0x0200, 0x0400, 0x0800, 0x1000, 0x2000, 0x4000 };  

	maxbpl = 0;
	for(ch = 0; ch < nch; ch++)
	{
		for(reg=s_reg; reg<e_reg; reg++)
		for(i = s_freq[reg]; i < e_freq[reg]; i++)
			if (maxbpl < cur_bpl[ch][reg][i]) 
				maxbpl = cur_bpl[ch][reg][i];
	}

	cw_len = available_len;
	for (bpl = maxbpl; bpl >= min_bpl; bpl--) {

		if (cw_len <= 0 ) 
		{
//			printf("\nbpl:%d, maxbpl:%d	", bpl, maxbpl);
			return (available_len-cw_len);
		}
		maxhalf = bit_mask[bpl]; 
		
		for(reg = s_reg; reg < e_reg; reg++)
	    for (i = s_freq[reg]; i < e_freq[reg]; i+=4) {
		for(ch = 0; ch < nch; ch++) {

			if(model_index[ch][reg][i/32]==0)
				continue;
	 		if ( cur_bpl[ch][reg][i] < bpl )
				continue; 

			if (coded_samp_bit[ch][reg][i]==0 || sign_coded[ch][reg][i]==1) {

				cband = i/32;
				i4 = i / 4;
				if (i%4==0) {
					enc_vec[ch][reg][i4] = 0;
				}
			
				for(j=0;j<4;j++)
				{
					if( abs(sample[ch][reg][i+j]) & maxhalf) 	m = 1;
					else 								    m = 0;

					enc_vec[ch][reg][i4] = (enc_vec[ch][reg][i4]<<1) | m;
					coded_samp_bit[ch][reg][i+j] = 
						(coded_samp_bit[ch][reg][i+j]<<1) | m;
					cur_bpl[ch][reg][i+j]-=1;
				}
				/* model selection */
				{
					int t[4];
					t[0] = (coded_samp_bit[ch][reg][i]&0x2)==0x2 ? 1:0;
					t[1] = (coded_samp_bit[ch][reg][i+1]&0x2)==0x2 ? 1:0;
					t[2] = (coded_samp_bit[ch][reg][i+2]&0x2)==0x2 ? 1:0;
					t[3] = (coded_samp_bit[ch][reg][i+3]&0x2)==0x2 ? 1:0;
					t[0] = (t[0]<<3) + (t[1]<<2) + (t[2]<<1) + (t[3]);

					modelIdx = BSHCModelSelect(model_index[ch][reg][cband],bpl);
					/* upper vector mapping, codebook reduction */
					t[0] = cb_tab[t[0]];
					k = EncodeBSHCQuant(modelIdx,t[0],enc_vec[ch][reg][i4]);
					
					cw_len -= k;
					if (coded_samp_bit[ch][reg][i] && sign_coded[ch][reg][i]==0)
					{
						if(sample[ch][reg][i] < 0)
							EncodeBSHCBin(1);
						else
							EncodeBSHCBin(0);
						cw_len--;
						sign_coded[ch][reg][i] = 0x01;
					}
					if (coded_samp_bit[ch][reg][i+1] && sign_coded[ch][reg][i+1]==0)
					{
						if(sample[ch][reg][i+1] < 0)
							EncodeBSHCBin(1);
						else
							EncodeBSHCBin(0);
						cw_len--;
						sign_coded[ch][reg][i+1] = 0x01;
					}
					if (coded_samp_bit[ch][reg][i+2] && sign_coded[ch][reg][i+2]==0)
					{
						if(sample[ch][reg][i+2] < 0)
							EncodeBSHCBin(1);
						else
							EncodeBSHCBin(0);
						cw_len--;
						sign_coded[ch][reg][i+2] = 0x01;
					}
					if (coded_samp_bit[ch][reg][i+3] && sign_coded[ch][reg][i+3]==0)
					{
						if(sample[ch][reg][i+3] < 0)
							EncodeBSHCBin(1);
						else
							EncodeBSHCBin(0);
						cw_len--;
						sign_coded[ch][reg][i+3] = 0x01;
					}
				}
			}
			if (cw_len <= 0 ) 
				return (available_len-cw_len);
		}
	} 
}
	return (available_len-cw_len);
}
