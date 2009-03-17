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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <string.h>

#include "global.h"

////////////////////////////////////////////////////////////////////////////////
#ifdef SVA_START_CODE_EMULATION

TLS unsigned char bit[8] = {0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01};

void c_avs_enc::OpenORABS(OutputStream *p, char *fname)
{
  p->uPreBytes      = 0xffffffff;
  p->iBytePosition    = 0;
  p->iBitOffset      = 0;
  p->iNumOfStuffBits    = 0;
  p->iBitsCount      = 0;
}
void c_avs_enc::CloseORABS(OutputStream *p)
{
  if(p->iBitOffset)
  {
    memcpy(pOutBuffer + nOutBufPtr, p->buf,p->iBytePosition+1);
    nOutBufPtr += p->iBytePosition+1;
  }
  else
  {
    memcpy(pOutBuffer, p->buf,p->iBytePosition);
    nOutBufPtr += p->iBytePosition;
  }
}
void c_avs_enc::FlushORABS(OutputStream *p)
{
  fflush(p->f);
}
_inline void c_avs_enc::write_1_bit(OutputStream *p,int_32_t b)
{
  int_32_t i;
  if(p->iBitOffset==6)
  {
    i = p->uPreBytes & 0x003fffff;
    if(i == 0)
    {
      p->buf[p->iBytePosition] = 0x02;
      p->iBytePosition++;
      p->iBitOffset    = 0;
      p->uPreBytes    = 0x00000002;
    }
  }
  if(p->iBytePosition == SVA_STREAM_BUF_SIZE)
  {
    memcpy((unsigned char*)p_avs_enc_frame->bitstream + p_avs_enc_frame->length, p->buf, SVA_STREAM_BUF_SIZE);
    p_avs_enc_frame->length += SVA_STREAM_BUF_SIZE;
    p->iBytePosition  = 0;
    p->iBitOffset    = 0;
  }
  p->uPreBytes <<= 1;
  if(b)
  {
    p->buf[p->iBytePosition] |= bit[p->iBitOffset];
    p->uPreBytes |= 1;
  }
  else
  {
    p->buf[p->iBytePosition] &= (~bit[p->iBitOffset]);
  }
  p->iBitOffset++;
  if(p->iBitOffset==8)
  {
    p->iBitOffset = 0;
    p->iBytePosition++;
  }
}

int_32_t c_avs_enc::write_n_bit(OutputStream *p,int_32_t b,int_32_t n)
{
  if(n>30) return 1;
  while(n>0)
  {
    write_1_bit(p,b&(0x01<<(n-1)));
    n--;
  }
  return 0;
}

int_32_t c_avs_enc::write_a_byte(OutputStream *p,int_32_t b)
{
  int_32_t i;
  int_32_t j;
  i = p->uPreBytes & 0x0000ffff;
  j=b & 0x000000fc;
  if ((p->iBitOffset==0)&&((j != 0)||(i !=0)))
  {

    if(p->iBytePosition == SVA_STREAM_BUF_SIZE)
    {
      memcpy((unsigned char*)p_avs_enc_frame->bitstream + p_avs_enc_frame->length, p->buf, SVA_STREAM_BUF_SIZE);
      p_avs_enc_frame->length += SVA_STREAM_BUF_SIZE;
      p->iBytePosition  = 0;
      p->iBitOffset    = 0;
    }
    p->buf[p->iBytePosition] = b ;
    p->iBytePosition++;
    p->uPreBytes <<= 8;
    p->uPreBytes +=b;

  }
  else
  {
    for(i=8;i>0;i--)
    {
      write_1_bit(p, b&(0x01<<(i-1)));
    }
  }
  return 0;
}
// one bit "1" is added to the end of stream, then some bits "0" are added to byte aligned position.
int_32_t c_avs_enc::write_align_stuff(OutputStream *p)
{
  unsigned char c;

  c = 0xff << ( 8 - p->iBitOffset );
  p->buf[p->iBytePosition] = ( c & p->buf[p->iBytePosition] ) | (0x80>>(p->iBitOffset));
  p->iBitsCount += 8 - p->iBitOffset;
  p->uPreBytes  = (p->uPreBytes << (8 - p->iBitOffset)) & c;
  p->iNumOfStuffBits  += 8 - p->iBitOffset;
  p->iBitOffset = 0;
  p->iBytePosition++;
  return 0;
}
//---end
int_32_t c_avs_enc::write_start_code(OutputStream *p, unsigned char code)
{
  if(p->iBitOffset)  write_align_stuff(p);

  if(p->iBytePosition >= SVA_STREAM_BUF_SIZE-4 && p->iBytePosition >0 )
  {
    memcpy((unsigned char*)p_avs_enc_frame->bitstream + p_avs_enc_frame->length, p->buf, p->iBytePosition);
    p_avs_enc_frame->length += p->iBytePosition;
    //wangyue
    p->iBytePosition  = 0;
    p->iBitOffset    = 0;
  }
  p->buf[p->iBytePosition  ] = 0;
  p->buf[p->iBytePosition+1] = 0;
  p->buf[p->iBytePosition+2] = 1;
  p->buf[p->iBytePosition+3] = code;
  p->iBytePosition += 4;

  p->uPreBytes  = (uint_32_t)code + 256;
  return 0;
}

/*
*************************************************************************
* Function:Open the output file for the bytestream
* Input: The filename of the file to be opened
* Output:
* Return: none.Function terminates the program in case of an error
* Attention:
*************************************************************************
*/


void c_avs_enc::OpenBitStreamFile(char *Filename)
{
  OpenORABS(pORABS, Filename);
}
void c_avs_enc::CloseBitStreamFile()
{
  CloseORABS(pORABS);
}

/*
*************************************************************************
* Function:Write sequence header information
* Input:
* Output:
* Return: sequence header length, including stuffing bits
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::WriteSequenceHeader()
{
  Bitstream *bitstream;
  byte SequenceHeader[MAXHEADERSIZE];
  int_32_t  bitscount=0;
  int_32_t  stuffbits;
  int_32_t  i,j,k;
  if ((bitstream= (Bitstream*)calloc(1, sizeof(Bitstream)))==NULL)
    no_mem_exit("Seuqence Header: bitstream");

  bitstream->streamBuffer = SequenceHeader;
  bitstream->bits_to_go = 8;

  input->profile_id = 0x20;
  input->level_id   = 0x42;

  input->display_horizontal_size = input->img_width;
  input->display_vertical_size   = input->img_height;
  input->sample_precision        = 1;
  input->bbv_buffer_size         = 4096;//here we random give a value,but in fact it is not true.

  //xzhao 20081122
  input->aspect_ratio_information = 0x2;
  input->bit_rate_lower = (input->bit_rate/400) & 0x3FFFF;
  input->bit_rate_upper = (input->bit_rate/400 - input->bit_rate_lower)>>18;



  bitscount += u_v(32, "seqence start code",      0x1B0,                           bitstream);
  bitscount += u_v(8,  "profile_id",              input->profile_id,               bitstream);
  bitscount += u_v(8,  "level_id"  ,              input->level_id,                 bitstream);

  bitscount += u_v(1,  "progressive_sequence",    input->progressive_sequence,     bitstream);
  bitscount += u_v(14, "picture width",           input->img_width,                bitstream);
  bitscount += u_v(14, "picture height",          input->img_height,               bitstream);
  bitscount += u_v(2,  "chroma foramt",           input->chroma_format,            bitstream);
  bitscount += u_v(3,  "sample precision",        input->sample_precision,         bitstream);
  bitscount += u_v(4,  "aspect ratio information",input->aspect_ratio_information, bitstream);
  bitscount += u_v(4,  "frame rate code",         input->frame_rate_code,          bitstream);

  bitscount += u_v(18, "bit rate lower",          input->bit_rate_lower,           bitstream);
  bitscount += u_v(1,  "marker bit",              1,                               bitstream);
  bitscount += u_v(12, "bit rate upper",          input->bit_rate_upper,           bitstream);
  bitscount += u_v(1,  "low delay",               input->low_delay,                bitstream);
  bitscount += u_v(1,  "marker bit",              1,                               bitstream);
  bitscount += u_v(18, "bbv buffer size",         input->bbv_buffer_size,          bitstream);
  bitscount += u_v(3,  "reserved bits",           0,                               bitstream);

  k = bitscount >> 3;
  j = bitscount % 8;

  stuffbits = 8-(bitscount%8);

  if (stuffbits<8)
    bitscount+=u_v(stuffbits,"stuff bits for byte align",0,bitstream);

  write_start_code(pORABS, 0xB0);


  for(i=4;i<k;i++)
    write_n_bit(pORABS,SequenceHeader[i],8);
  //write_a_byte(pORABS,SequenceHeader[i]);
  //wangyue
  write_n_bit(pORABS,SequenceHeader[k],j);
  //wangyue
  write_align_stuff(pORABS);

  free(bitstream);

  return bitscount;
}

int_32_t c_avs_enc::WriteSequenceEnd()
{
  write_start_code(pORABS, 0xb1);
  return 32;
}

/*
*************************************************************************
* Function:Write sequence display extension information
* Input:
* Output:
* Return: sequence display extension information lenght
* Attention:
*************************************************************************
*/


int_32_t c_avs_enc::WriteSequenceDisplayExtension()
{
  Bitstream *bitstream;
  byte SequenceDisplayExtension[MAXHEADERSIZE];
  int_32_t  bitscount=0;
  int_32_t  stuffbits;
  int_32_t  i,j,k;

  memset(SequenceDisplayExtension, 0, MAXHEADERSIZE*sizeof(byte));
  if ((bitstream = (Bitstream*)calloc(1, sizeof(Bitstream)))==NULL)
    no_mem_exit("Sequence Display Extension: bitstream");
  input->video_format = 1;
  input->video_range  = 1;
  input->display_horizontal_size = img->width;
  input->display_vertical_size   = img->height;

  bitstream->streamBuffer = SequenceDisplayExtension;
  bitstream->bits_to_go = 8;

  bitscount += u_v(32,"sequence display extension start code",0x1B5,                    bitstream);
  bitscount += u_v(4, "extension id",                         2,                        bitstream);
  bitscount += u_v(3, "video format",                         input->video_format,      bitstream);
  bitscount += u_v(1, "video range",                          input->video_range,       bitstream);
  bitscount += u_v(1, "color description",                    input->color_description, bitstream);

  if(input->color_description)
  {
    bitscount += u_v(8,"color primaries",          input->color_primaries,          bitstream);
    bitscount += u_v(8,"transfer characteristics", input->transfer_characteristics, bitstream);
    bitscount += u_v(8,"matrix coefficients",      input->matrix_coefficients,      bitstream);
  }

  bitscount += u_v(14, "display horizontal size",input->display_horizontal_size, bitstream);
  bitscount += u_v(1,  "marker bit",             1,                              bitstream);
  bitscount += u_v(14, "display vertical size",  input->display_vertical_size,   bitstream);
  bitscount += u_v(2,  "reserved bits",          0,                              bitstream);

  k = bitscount / 3;
  j = bitscount % 8;

  stuffbits = 8-(bitscount%8);

  if (stuffbits<8)
  {
    bitscount += u_v(stuffbits,"stuff bits for byte align",0,bitstream);
  }

  write_start_code(pORABS, 0xb5);
  for(i=4; i<k; i++)
    write_n_bit(pORABS, SequenceDisplayExtension[i], 8);

  //write_a_byte(pORABS, SequenceDisplayExtension[i]);
  //wangyue
  write_n_bit(pORABS, SequenceDisplayExtension[k], j);
  write_align_stuff(pORABS);

  free(bitstream);

  return bitscount;
}

/*
*************************************************************************
* Function:Write copyright extension information
* Input:
* Output:
* Return: copyright extension information lenght
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc::WriteCopyrightExtension()
{
  Bitstream *bitstream;
  byte CopyrightExtension[MAXHEADERSIZE];
  int_32_t  bitscount=0;
  int_32_t  stuffbits;
  int_32_t  i,j,k;

  if ((bitstream = (Bitstream*) calloc(1, sizeof(Bitstream)))==NULL)
    no_mem_exit("Copyright Extension: bitstream");

  bitstream->streamBuffer = CopyrightExtension;
  bitstream->bits_to_go = 8;

  bitscount+=u_v(32,"copyright extension start code",0x1b5,bitstream);
  bitscount+=u_v(4,"extension id",4,bitstream);
  bitscount+=u_v(3,"copyright flag",cp->copyright_flag,bitstream);
  bitscount+=u_v(1,"copyright id",cp->copyright_id,bitstream);
  bitscount+=u_v(1,"original or copy",cp->original_or_copy,bitstream);
  bitscount+=u_v(64,"copyright number",cp->copyright_number,bitstream);

  k = bitscount >> 3;
  j = bitscount % 8;

  stuffbits = 8-(bitscount%8);

  if (stuffbits<8)
  {
    bitscount+=u_v(stuffbits,"stuff bits for byte align",0,bitstream);
  }

  write_start_code(pORABS, 0xb5);
  for(i=4;i<k;i++)
    write_n_bit(pORABS,CopyrightExtension[i],8);
  //write_a_byte(pORABS,CopyrightExtension[i]);
  //wangyue
  write_n_bit(pORABS,CopyrightExtension[k],j);
  write_align_stuff(pORABS);
  //  fwrite(CopyrightExtension,1,bitscount/8,f);

  free(bitstream);

  return bitscount;
}


/*
*************************************************************************
* Function:Write camera parameter extension information
* Input:
* Output:
* Return: camera parameter  extension information lenght
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc::WriteCameraParametersExtension()
{
  Bitstream *bitstream;
  byte CameraParametersExtension[MAXHEADERSIZE];
  int_32_t  bitscount=0;
  int_32_t  stuffbits;
  int_32_t  i,j,k;

  if ((bitstream = (Bitstream*)calloc(1, sizeof(Bitstream)))==NULL)
    no_mem_exit("Camera Parameters Extension: bitstream");

  bitstream->streamBuffer = CameraParametersExtension;
  bitstream->bits_to_go = 8;

  bitscount+=u_v(32,"camera parameters extension start code",0x1b5,bitstream);
  bitscount+=u_v(4,"extension id",11,bitstream);
  bitscount+=u_v(1,"reserved",0,bitstream);
  bitscount+=u_v(7,"camera id",camera->camera_id,bitstream);

  bitscount+=u_v(22,"height_of_image_device",camera->height_of_image_device,bitstream);
  bitscount+=u_v(22,"focal_length",camera->focal_length,bitstream);
  bitscount+=u_v(22,"f_number",camera->f_number,bitstream);
  bitscount+=u_v(22,"vertical_angle_of_view",camera->vertical_angle_of_view,bitstream);
  bitscount+=u_v(32,"camera_position_x",camera->camera_direction_x,bitstream);
  bitscount+=u_v(32,"camera_position_y",camera->camera_direction_y,bitstream);
  bitscount+=u_v(22,"camera_position_z",camera->camera_direction_z,bitstream);
  bitscount+=u_v(22,"camera_direction_x",camera->camera_direction_x,bitstream);
  bitscount+=u_v(22,"camera_direction_y",camera->camera_direction_y,bitstream);
  bitscount+=u_v(22,"camera_direction_z",camera->camera_direction_z,bitstream);
  bitscount+=u_v(22,"image_plane_vertical_x",camera->image_plane_vertical_x,bitstream);
  bitscount+=u_v(22,"image_plane_vertical_y",camera->image_plane_vertical_y,bitstream);
  bitscount+=u_v(32,"image_plane_vertical_z",camera->image_plane_vertical_z,bitstream);

  k = bitscount >> 3;
  j = bitscount % 8;

  stuffbits = 8-(bitscount%8);

  if (stuffbits<8)
  {
    bitscount+=u_v(stuffbits,"stuff bits for byte align",0,bitstream);
  }

  write_start_code(pORABS, 0xb5);
  for(i=4;i<k;i++)
    write_n_bit(pORABS,CameraParametersExtension[i],8);
  //write_a_byte(pORABS,CameraParametersExtension[i]);
  //wangyue
  //write_n_bit(pORABS,CameraParametersExtension[k],j);
  write_align_stuff(pORABS);
  //  fwrite(CameraParametersExtension,1,bitscount/8,f);

  free(bitstream);

  return bitscount;
}

/*
*************************************************************************
* Function:Write user data
* Input:
* Output:
* Return: user data length
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::WriteUserData(char *userdata)
{
  Bitstream *bitstream;
  byte UserData[MAXHEADERSIZE];
  int_32_t  bitscount=0;

  if ((bitstream = (Bitstream*) calloc(1, sizeof(Bitstream)))==NULL) no_mem_exit("User data: bitstream");
  bitstream->streamBuffer = UserData;
  bitstream->bits_to_go = 8;

  bitscount += u_v(32,"user data start code", 0x1B2, bitstream);
  write_start_code(pORABS, 0xB2);
  while (*userdata)
  {
    write_n_bit(pORABS, *userdata, 8);
    //write_a_byte(pORABS, *userdata);
    //wangyue
    bitscount += u_v(8, "user data", *userdata++, bitstream);
  }
  write_align_stuff(pORABS);
  free(bitstream);
  return bitscount;
}

/*
*************************************************************************
* Function:Write bit steam of one slice to file
* Input:
* Output:
* Return: none
* Attention:
*************************************************************************
*/



void c_avs_enc::WriteSlicetoFile()
{
  int_32_t n, i;
  n = currBitStream->byte_pos;


  write_start_code(pORABS, currBitStream->streamBuffer[3]);
  for(i=4;i<n;i++)
    write_n_bit(pORABS, currBitStream->streamBuffer[i],8);
  //   write_a_byte(pORABS, currBitStream->streamBuffer[i]);
  //wangyue

  write_align_stuff(pORABS);

  stat->bit_ctr += (uint_32_t)8*n;

  currBitStream->byte_pos = 0;
}

/*
*************************************************************************
* Function:Write bit steam to file
* Input:
* Output:
* Return: none
* Attention:
*************************************************************************
*/

void c_avs_enc::WriteBitstreamtoFile()
{
  int_32_t i;
  for(i=0;i<currBitStream->byte_pos;i++)
  {
    if(currBitStream->streamBuffer[i]==0 && currBitStream->streamBuffer[i+1]==0 && currBitStream->streamBuffer[i+2]==1)
    {
      write_start_code(pORABS, currBitStream->streamBuffer[i+3]);
      i=i+4;
    }
    write_n_bit(pORABS, currBitStream->streamBuffer[i],8);
  }
  write_align_stuff(pORABS);
  stat->bit_ctr += (uint_32_t)(8*currBitStream->byte_pos);
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////

/*
*************************************************************************
* Function:
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
void c_avs_enc::error(char *text, int_32_t code)
{
  fprintf(stderr, "%s\n", text);
  exit(code);
}
/*
*************************************************************************
* Function:
* Input:
* Output:
* Return:
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc::start_sequence()
{
  int_32_t len = 0;
  char id_string[255] = "OPEN_SOURCE_TRANSCODER";
  OpenBitStreamFile(input->outfile);
  len = WriteSequenceHeader();
  if (strlen(id_string) > 1)
    len += WriteUserData(id_string);
  return len;
}
/*
*************************************************************************
* Function:
* Input:
* Output:
* Return:
* Attention:Mainly flushing of everything Add termination symbol, etc.
*************************************************************************
*/
int_32_t c_avs_enc::terminate_sequence()
{
  int_32_t len;
  len = WriteSequenceEnd();
  CloseBitStreamFile();
  return len;   // make lint happy
}
