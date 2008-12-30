/*$T mpeg2avs.cpp GC 1.140 10/28/07 10:20:49 */

/* mpeg2avs.c */
//cl  mpeg2avs.cpp /I lencod /I mpeg2dec /c /D "WIN32"
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/timeb.h>
#include <memory.h>
#include <string.h>

#include "avs_enc_lib.h"
#include "mpeg2dec.h"

#include "mpeg2avs.h"

/*
 =======================================================================================================================
 =======================================================================================================================
 */

int_32_t trans_create(trans_create_t *create, char *encoder_cfg_file_name, c_avs_enc *p_avs_enc)
{
	memset(&create->decCreate, 0, sizeof(mpeg2_dec_create_t));
	create->decCreate.pBitstream = create->pBitstream;
	mpeg2_decoder_create(&create->decCreate);
	create->pYUVData = (unsigned char *) malloc(create->decCreate.width * create->decCreate.height * 3 / 2);

	/* init avs encoder */
	strcpy(create->avsCreate.strConfigFile, encoder_cfg_file_name);
	create->avsCreate.frame_rate_code = 3;

	/* Ö¡ÂÊ(1: 24000/1001,2: 24,3: 25,4: 30000/1001,5: 30,6: 50,7: 60000/1001,8: 60) */
	create->avsCreate.height = create->decCreate.height;
	create->avsCreate.width = create->decCreate.width;
	p_avs_enc->avs_enc_create(&create->avsCreate);
	create->pMbInfo = (MB_INFO *) malloc(sizeof(MB_INFO) * create->decCreate.width * create->decCreate.height / 256);
	create->decCreate.pBitstream = NULL;
	return 0;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int_32_t trans_transcode(trans_create_t *create, int_32_t transcoding_type, c_avs_enc *p_avs_enc)
{
	create->decFrame.bitstream = create->frame->pMpeg2Stream;
	create->decFrame.output = create->pYUVData;
	create->decFrame.MbInfo = create->pMbInfo;

	create->avsFrame.bitstream = create->frame->pAvsStream;

	mpeg2_decoder_decode(&create->decFrame, &create->decStats, transcoding_type);

	create->avsFrame.lastFrame = 0;
	create->avsFrame.type = create->decStats.type;
	create->avsFrame.version = AVS_RM_52C;
	create->avsFrame.input = create->decFrame.output;
	create->avsFrame.pEncMBInfo = create->decFrame.MbInfo;

	p_avs_enc->avs_enc_encode(&create->avsFrame, &create->avsStats);

	create->frame->nAvsPtr = create->avsFrame.length;
	create->frame->nMpeg2Ptr = create->decFrame.length;

	return 0;
}

/*
 =======================================================================================================================
 =======================================================================================================================
 */
int_32_t trans_destroy(trans_create_t *create, c_avs_enc *p_avs_enc)
{
	/*
	 * mpeg2_decoder_destroy( &decFrame, &decStats);
	 */
	create->avsFrame.bitstream = create->frame->pAvsStream;
	p_avs_enc->avs_enc_destroy(&create->avsFrame);
	create->frame->nAvsPtr = create->avsFrame.length;

	free(create->pYUVData);
	free(create->pMbInfo);
	return 0;
}
