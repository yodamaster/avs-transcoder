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
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <loopfilter.h>
#include "global.h"
void c_avs_enc:: stuffing_byte(int_32_t n)
  {
  int_32_t i;
  Bitstream *currStream;

  currStream = currBitStream;

  for(i=0; i<n; i++)
    {
    currStream->streamBuffer[currStream->byte_pos++] = 0x80;
    currStream->bits_to_go = 8;
    currStream->byte_buf   = 0;
    }
  }

int_32_t c_avs_enc:: start_slice()
{
  Bitstream *currStream;

  currStream = currBitStream;

  if (currStream->bits_to_go !=8 )
  {
    currStream->byte_buf <<= currStream->bits_to_go;
    currStream->byte_buf |= ( 1 << (currStream->bits_to_go - 1) );
    currStream->streamBuffer[currStream->byte_pos++] = currStream->byte_buf;
    currStream->bits_to_go = 8;
    currStream->byte_buf = 0;
  }
   else		// cjw 20060321 
  {
	  currStream->streamBuffer[currStream->byte_pos++] = 0x80;
	  currStream->bits_to_go = 8;
	  currStream->byte_buf = 0; 
  }

  return 0;
}

/*
*************************************************************************
* Function:This function terminates a picture
* Input:
* Output:
* Return: 0 if OK,                                                         \n
      1 in case of error
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::terminate_picture()
{
  Bitstream *currStream;

  currStream = currBitStream;
  currStream->byte_buf <<= currStream->bits_to_go;
  currStream->byte_buf |= (1 << (currStream->bits_to_go - 1) );
  currStream->streamBuffer[currStream->byte_pos++] = currStream->byte_buf;
  currStream->bits_to_go = 8;
  currStream->byte_buf = 0;

  return 0;
}

void c_avs_enc:: picture_data( )
  {
  Boolean end_of_picture = FALSE;
  int_32_t CurrentMbNumber=0;
  int_32_t MBRowSize = img->img_width_in_mb;
  int_32_t slice_nr = 0;
  int_32_t slice_qp = img->qp;
  int_32_t len, i, j;
  //init the intra pred mode
  for(i=0; i<img->width/B8_SIZE+100; i++)
  {
    for(j=0; j<img->height/B8_SIZE+100; j++)
    {
      img->ipredmode[i][j] = -1;
    }
  }

  for(i=0; i<img->width*img->height/256; i++)
  {
    img->mb_data[i].slice_nr = -1;
  }
  if (input->rdopt)
    {
    switch(img->type)
      {
      case INTRA_IMG:
        encode_one_macroblock = &c_avs_enc::encode_one_intra_macroblock_rdo;
        break;
      case INTER_IMG:
        encode_one_macroblock = &c_avs_enc::encode_one_inter_macroblock_rdo;
        break;
      case B_IMG:
        encode_one_macroblock = &c_avs_enc::encode_one_b_frame_macroblock_rdo;
        break;
      }
    }
  else
    {
      // xzhao
      //img->type=INTRA_IMG;
    switch(img->type)
      {
      case INTRA_IMG:
          encode_one_macroblock = &c_avs_enc::encode_one_intra_macroblock_not_rdo;
        break;
      case INTER_IMG:
          encode_one_macroblock = &c_avs_enc::encode_one_inter_macroblock_not_rdo;
        break;
      case B_IMG:
          encode_one_macroblock = &c_avs_enc::encode_one_b_frame_macroblock_not_rdo;
        break;
      }
    }

  while (end_of_picture == FALSE) // loop over macroblocks
    {
    set_MB_parameters(CurrentMbNumber);
    if (input->slice_row_nr && (img->current_mb_nr ==0 ||(img->current_mb_nr>0 && img->mb_data[img->current_mb_nr].slice_nr != img->mb_data[img->current_mb_nr-1].slice_nr)))
      {
      start_slice ();
      img->current_slice_qp = img->qp;
      img->current_slice_start_mb = img->current_mb_nr;

      len = SliceHeader(slice_nr, slice_qp);

      img->current_slice_nr = slice_nr;
      stat->bit_slice += len;
      slice_nr++;
      }
    start_macroblock();
    (this->*encode_one_macroblock)();

    write_one_macroblock(1);
    terminate_macroblock (&end_of_picture);

    proceed2nextMacroblock ();
    CurrentMbNumber++;
    }

  terminate_picture ();

  if(!input->loop_filter_disable)     //xzhao 20081227
	  DeblockFrame (img, imgY, imgUV);  //wangyue

}


void c_avs_enc:: top_field(Picture *pic)
{
  Boolean end_of_picture = FALSE;
  int_32_t CurrentMbNumber=0;
  int_32_t MBRowSize = img->width / MB_BLOCK_SIZE;
  int_32_t slice_nr =0;
  int_32_t slice_qp = img->qp;
  int_32_t len;

  img->top_bot = 0;  // Yulj 2004.07.20
  while (end_of_picture == FALSE) // loop over macroblocks
  {
    set_MB_parameters (CurrentMbNumber);
    if (input->slice_row_nr && (img->current_mb_nr ==0
      ||(img->current_mb_nr>0 && img->mb_data[img->current_mb_nr].slice_nr != img->mb_data[img->current_mb_nr-1].slice_nr)))
    {
      // slice header start jlzheng 7.1
      start_slice ();
      img->current_slice_qp = img->qp;
      img->current_slice_start_mb = img->current_mb_nr;
      len = SliceHeader(slice_nr, slice_qp);
      stat->bit_slice += len;
      slice_nr++;
      // slice header end
    }

    img->current_mb_nr_fld = img->current_mb_nr;
    start_macroblock ();
    (this->*encode_one_macroblock) ();
    write_one_macroblock (1);
    terminate_macroblock (&end_of_picture);
    proceed2nextMacroblock ();
    CurrentMbNumber++;
  }

  if(!input->loop_filter_disable)     // xzhao 20081227
	  DeblockFrame (img, imgY, imgUV);    // wangyue

  //rate control
  pic->bits_per_picture = 8 * (currBitStream->byte_pos);
}

void c_avs_enc:: bot_field(Picture *pic)
{
  Boolean end_of_picture = FALSE;
  int_32_t CurrentMbNumber=0;
  int_32_t MBRowSize = img->width / MB_BLOCK_SIZE;
  int_32_t slice_nr =0;
  int_32_t slice_qp = img->qp;
  int_32_t len;

  img->top_bot = 1; //Yulj 2004.07.20
  while (end_of_picture == FALSE) // loop over macroblocks
  {
    set_MB_parameters (CurrentMbNumber);
    if (input->slice_row_nr && (img->current_mb_nr ==0
      ||(img->current_mb_nr>0 && img->mb_data[img->current_mb_nr].slice_nr != img->mb_data[img->current_mb_nr-1].slice_nr)))
    {
      // slice header start  jlzheng 7.11
      start_slice ();
      img->current_slice_qp = img->qp;
      img->current_slice_start_mb = img->current_mb_nr;
      len = SliceHeader(slice_nr, slice_qp);
      stat->bit_slice += len;
      slice_nr++;
      // slice header end
    }

    img->current_mb_nr_fld = img->current_mb_nr+img->total_number_mb;
    start_macroblock ();
    (this->*encode_one_macroblock) ();
    write_one_macroblock (1);
    terminate_macroblock (&end_of_picture);
    proceed2nextMacroblock ();
    CurrentMbNumber++;
  }

  terminate_picture ();

  if(!input->loop_filter_disable)     //xzhao 20081227
	  DeblockFrame (img, imgY, imgUV); //wangyue

  pic->bits_per_picture = 8 * (currBitStream->byte_pos);
}

void c_avs_enc::store_field_MV ()
{
    int_32_t i, j;

    if (img->type != B_IMG)     //all I- and P-frames
  {
        if (!img->picture_structure)
    {
            for (i = 0; i < img->width / 8 + 4; i++)
      {
                for (j = 0; j < img->height / 16; j++)
        {
                    tmp_mv_frm[0][2 * j][i] = tmp_mv_frm[0][2 * j + 1][i] =
            tmp_mv_top[0][j][i];
                    tmp_mv_frm[0][2 * j][i] = tmp_mv_frm[0][2 * j + 1][i] = tmp_mv_top[0][j][i];        // ??
                    tmp_mv_frm[1][2 * j][i] = tmp_mv_frm[1][2 * j + 1][i] =
            tmp_mv_top[1][j][i] * 2;
                    tmp_mv_frm[1][2 * j][i] = tmp_mv_frm[1][2 * j + 1][i] = tmp_mv_top[1][j][i] * 2;    // ??

                    if (refFrArr_top[j][i] == -1)
          {
            refFrArr_frm[2 * j][i] =
              refFrArr_frm[2 * j + 1][i] = -1;
          }
                    else
          {
            refFrArr_frm[2 * j][i] =
              refFrArr_frm[2 * j + 1][i] =
              (int_32_t) (refFrArr_top[j][i] / 2);
          }
        }
      }
    }
        else
    {
            for (i = 0; i < img->width / 8 + 4; i++)
      {
                for (j = 0; j < img->height / 16; j++)
        {
                    tmp_mv_top[0][j][i] = tmp_mv_bot[0][j][i] =
            (int_32_t) (tmp_mv_frm[0][2 * j][i]);
                    tmp_mv_top[1][j][i] = tmp_mv_bot[1][j][i] =
            (int_32_t) ((tmp_mv_frm[1][2 * j][i]) / 2);


          if (refFrArr_frm[2 * j][i] == -1)
          {
            refFrArr_top[j][i] = refFrArr_bot[j][i] = -1;
          }
          else
          {
            refFrArr_top[j][i] = refFrArr_bot[j][i] =
              refFrArr_frm[2 * j][i] * 2;
          }
        }
      }
          }
      }
}

/*
*************************************************************************
* Function: allocates the memory for the coded picture data
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
void c_avs_enc:: AllocateBitstream()
{
  const int_32_t buffer_size = (bytes_y << 2);

  if ((currBitStream = (Bitstream *) calloc(1, sizeof(Bitstream))) == NULL)
    no_mem_exit ("malloc_slice: Bitstream");
  if ((currBitStream->streamBuffer = (byte *) calloc(buffer_size, sizeof(byte))) == NULL)
    no_mem_exit ("malloc_slice: StreamBuffer");

  currBitStream->bits_to_go = 8;
}
/*
*************************************************************************
* Function:free the allocated memory for the coded picture data
* Input:
* Output:
* Return:
*************************************************************************
*/
void c_avs_enc:: FreeBitstream()
{
  const int_32_t buffer_size = (img->width * img->height * 4);

  if (currBitStream->streamBuffer)
    free(currBitStream->streamBuffer);
  if (currBitStream)
    free(currBitStream);
}