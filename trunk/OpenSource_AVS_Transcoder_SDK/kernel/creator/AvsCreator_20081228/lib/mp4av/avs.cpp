/*
 * The contents of this file are subject to the Mozilla Public
 * License Version 1.1 (the "License"); you may not use this file
 * except in compliance with the License. You may obtain a copy of
 * the License at http://www.mozilla.org/MPL/
 * 
 * Software distributed under the License is distributed on an "AS
 * IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * rights and limitations under the License.
 * 
 * The Original Code is MPEG4IP.
 * 
 * The Initial Developer of the Original Code is Cisco Systems Inc.
 * Portions created by Cisco Systems Inc. are
 * Copyright (C) Cisco Systems Inc. 2004.  All Rights Reserved.
 * 
 * Contributor(s): 
 *		Bill May wmay@cisco.com
 */

#include "mpeg4ip.h"
#include "mp4av_avs.h"
#include "mpeg4ip_bitstream.h"
//#define BOUND_VERBOSE 1

static uint8_t exp_golomb_bits[256] = {
8, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 
3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 
};

uint32_t avs_ue (CBitstream *bs)
{
  uint32_t bits, read;
  int bits_left;
  uint8_t coded;
  bool done = false;
  bits = 0;
  // we want to read 8 bits at a time - if we don't have 8 bits, 
  // read what's left, and shift.  The exp_golomb_bits calc remains the
  // same.
  while (done == false) {
    bits_left = bs->bits_remain();
    if (bits_left < 8) {
      read = bs->PeekBits(bits_left) << (8 - bits_left);
      done = true;
    } else {
      read = bs->PeekBits(8);
      if (read == 0) {
	bs->GetBits(8);
	bits += 8;
      } else {
	done = true;
      }
    }
  }
  coded = exp_golomb_bits[read];
  bs->GetBits(coded);
  bits += coded;

  //  printf("ue - bits %d\n", bits);
  return bs->GetBits(bits + 1) - 1;
}

int32_t avs_se (CBitstream *bs) 
{
  uint32_t ret;
  ret = avs_ue(bs);
  if ((ret & 0x1) == 0) {
    ret >>= 1;
    int32_t temp = 0 - ret;
    return temp;
  } 
  return (ret + 1) >> 1;
}

static void avs_decode_annexb( uint8_t *dst, int *dstlen,
                                const uint8_t *src, const int srclen )
{
  uint8_t *dst_sav = dst;
  const uint8_t *end = &src[srclen];

  while (src < end)
  {
    if (src < end - 3 && src[0] == 0x00 && src[1] == 0x00 &&
        src[2] == 0x03)
    {
      *dst++ = 0x00;
      *dst++ = 0x00;

      src += 3;
      continue;
    }
    *dst++ = *src++;
  }

  *dstlen = dst - dst_sav;
}

extern "C" bool avs_is_start_code (const uint8_t *pBuf) 
{
  if (pBuf[0] == 0 && 
      pBuf[1] == 0 && 
      (pBuf[2] == 1)) 
  {
    return true;
  }
  return false;
}

extern "C" uint32_t avs_find_next_start_code (const uint8_t *pBuf, 
					       uint32_t bufLen)
{
  uint32_t val, temp;
  uint32_t offset;

  offset = 0;

  if (pBuf[0] == 0 && 
      pBuf[1] == 0 && 
      pBuf[2] == 1) {
    pBuf += 3;
    offset = 3;
  }

  val = 0xffffffff;
  while (offset < bufLen - 3) {
    val <<= 8;
    temp = val & 0xff000000;
    val &= 0x00ffffff;
    val |= *pBuf++;
    offset++;
    if (val == AVS_START_CODE) {
      return offset - 3;
    }
  }
  return 0;
}

extern "C" uint8_t avs_nal_unit_type (const uint8_t *buffer)
{
  if (buffer[3] == 0xB0) {
	  return AVS_NAL_TYPE_SEQ_HEADER;
  }else if (buffer[3] == 0xB5) {
	  return AVS_NAL_TYPE_VIDEO_EXTENSION;
  }else if (buffer[3] == 0xB2) {
	  return AVS_NAL_TYPE_USER_DATA;
  }else if (buffer[3] == 0xB7) {
	  return AVS_NAL_TYPE_VIDEO_EDIT;
  }else if (buffer[3] == 0xB3) {
	  return AVS_NAL_TYPE_I_PIC_HEADER;
  }else if ((buffer[3] == 0xB6) && ((buffer[6] >> 6 ) == 0x01)) {
	  return AVS_NAL_TYPE_P_PIC_HEADER;
  }else if ((buffer[3] == 0xB6) && ((buffer[6] >> 6 ) == 0x02)) {
	  return AVS_NAL_TYPE_B_PIC_HEADER;
  }else if ((buffer[3] >= 0x00) && (buffer[3] <= 0xAF)) {
	  return AVS_NAL_TYPE_I_SLICE;
  }
}

extern "C" int avs_nal_unit_type_is_slice (uint8_t type)
{
  if ( (type >= AVS_NAL_TYPE_I_SLICE) && (type <= AVS_NAL_TYPE_B_SLICE) ) {
    return true;
  }
  return false;
}

/*
 * determine if the slice we decoded is a sync point
 */
extern "C" bool avs_picture_is_idr (avs_decode_t *dec) 
{
  if (!dec->pic_is_idr)
    return false;
  if (dec->picture_coding_type == AVS_TYPE_I) return true;
  return false;
}

extern "C" uint8_t avs_nal_ref_idc (const uint8_t *buffer)
{
  return 0;
}

int avs_read_seq_header (const uint8_t *buffer, 
			uint32_t buflen, 
			avs_decode_t *dec)
{
  CBitstream bs;
  uint32_t header;
  uint8_t tmp[2048];
  int tmp_len;

  header = 4;
  
  avs_decode_annexb( tmp, &tmp_len, buffer + header, MIN(buflen-header,2048) );
  bs.init(tmp, tmp_len * 8);

  //bs.set_verbose(true);
  try {
    dec->profile = bs.GetBits(8);
    dec->level = bs.GetBits(8);
    bs.GetBits(1); // progressive_sequence
    dec->horizontal_size = bs.GetBits(14);
	dec->vertical_size = bs.GetBits(14);
	bs.GetBits(2); // chroma_format
	bs.GetBits(3); // sample_precision
	bs.GetBits(4); // aspect_ratio

	uint32_t temp = bs.GetBits(4);
	switch( temp ) {
	case 1 :
		dec->frame_rate = (double)24000 / 1001;	   
		break;
	case 2 :
		dec->frame_rate = 24;
		break;
	case 3 :
		dec->frame_rate = 25;	   
		break;
	case 4 :
		dec->frame_rate = (double)30000/1001;
		break;
	case 5 :
		dec->frame_rate = 30;	   
		break;
	case 6 :
		dec->frame_rate = 50;
		break;
	case 7 :
		dec->frame_rate = (double)60000 / 1001;	   
		break;
	case 8 :
		dec->frame_rate = 60;
		break;
	default:
		exit (0);
	}
	
  } catch (...) {
    return -1;
  }
  return 0;
}

int avs_read_i_pic_hdr_info (const uint8_t *buffer, 
					   uint32_t buflen, 
					   avs_decode_t *dec)
{
	CBitstream bs;
	uint32_t header;
	uint8_t tmp[10]; 
	int tmp_len;
	
	header = 4;
	
	avs_decode_annexb( tmp, &tmp_len, buffer + header, 10);
	bs.init(tmp, tmp_len * 8);
	
	//bs.set_verbose(true);
	try {
		bs.GetBits(16);
		uint32_t temp = bs.GetBits(1);
		if (temp) {
			bs.GetBits(24);
		}
		dec->picture_distance = bs.GetBits(8);
		
	} catch (...) {
		return -1;
	}
	return 0;
}

int avs_read_pb_pic_hdr_info (const uint8_t *buffer, 
						   uint32_t buflen, 
						   avs_decode_t *dec)
{
	CBitstream bs;
	uint32_t header;
	uint8_t tmp[10];
	int tmp_len;
	
	header = 4;
	
	avs_decode_annexb( tmp, &tmp_len, buffer + header, 10);
	bs.init(tmp, tmp_len * 8);
	
	//bs.set_verbose(true);
	try {
		bs.GetBits(16);
		dec->picture_coding_type = bs.GetBits(2);
		dec->picture_distance = bs.GetBits(8);
		
	} catch (...) {
		return -1;
	}
	return 0;
}

extern "C" int avs_find_slice_type (const uint8_t *buffer, 
				     uint32_t buflen,
				     uint8_t *slice_type, 
				     bool noheader)
{
  uint32_t header;
  if (noheader) header = 1;
  else {
    header = 4;
  }
  CBitstream bs;
  bs.init(buffer + header, (buflen - header) * 8);
  try {
    avs_ue(&bs); // first_mb_in_slice
    *slice_type = avs_ue(&bs); // slice type
  } catch (...) {
    return -1;
  }
  return 0;
}

int avs_read_slice_info (const uint8_t *buffer, 
			  uint32_t buflen, 
			  avs_decode_t *dec)
{
  uint32_t header;
  uint8_t tmp[512]; /* Enough for the begining of the slice header */
  int tmp_len;

  header = 4;
  CBitstream bs;

  avs_decode_annexb( tmp, &tmp_len, buffer + header, MIN(buflen-header,512) );
  bs.init(tmp, tmp_len * 8);
  try {
    avs_ue(&bs); // first_mb_in_slice
  } catch (...) {
    return -1;
  }
  return 0;
}

static void avs_compute_poc( avs_decode_t *dec ) {
  const int max_frame_num = 1 << (dec->log2_max_frame_num_minus4 + 4);
  int field_poc[2] = {0,0};
  enum {
    AVS_PICTURE_FRAME,
    AVS_PICTURE_FIELD_TOP,
    AVS_PICTURE_FIELD_BOTTOM,
  } pic_type;

  /* FIXME FIXME it doesn't handle the case where there is a MMCO == 5
   * (MMCO 5 "emulates" an idr) */
  
  /* picture type */
  if (dec->frame_mbs_only_flag || !dec->field_pic_flag)
    pic_type = AVS_PICTURE_FRAME;
  else if (dec->bottom_field_flag)
    pic_type = AVS_PICTURE_FIELD_BOTTOM;
  else
    pic_type = AVS_PICTURE_FIELD_TOP;

  /* frame_num_offset */
  if (dec->nal_unit_type == AVS_NAL_TYPE_I_PIC_HEADER) {
    dec->pic_order_cnt_lsb_prev = 0;
    dec->pic_order_cnt_msb_prev = 0;
    dec->frame_num_offset = 0;
  } else {
    if (dec->frame_num < dec->frame_num_prev)
      dec->frame_num_offset = dec->frame_num_offset_prev + max_frame_num;
    else
      dec->frame_num_offset = dec->frame_num_offset_prev;
  }

  /* */
  if(dec->pic_order_cnt_type == 0) {
    const unsigned int max_poc_lsb = 1 << (dec->log2_max_pic_order_cnt_lsb_minus4 + 4);

    if (dec->pic_order_cnt_lsb < dec->pic_order_cnt_lsb_prev &&
        dec->pic_order_cnt_lsb_prev - dec->pic_order_cnt_lsb >= max_poc_lsb / 2)
      dec->pic_order_cnt_msb = dec->pic_order_cnt_msb_prev + max_poc_lsb;
    else if (dec->pic_order_cnt_lsb > dec->pic_order_cnt_lsb_prev &&
             dec->pic_order_cnt_lsb - dec->pic_order_cnt_lsb_prev > max_poc_lsb / 2)
      dec->pic_order_cnt_msb = dec->pic_order_cnt_msb_prev - max_poc_lsb;
    else
      dec->pic_order_cnt_msb = dec->pic_order_cnt_msb_prev;

    field_poc[0] = dec->pic_order_cnt_msb + dec->pic_order_cnt_lsb;
    field_poc[1] = field_poc[0];
    if (pic_type == AVS_PICTURE_FRAME)
      field_poc[1] += dec->delta_pic_order_cnt_bottom;

  } else if (dec->pic_order_cnt_type == 1) {
    int abs_frame_num, expected_delta_per_poc_cycle, expected_poc;

    if (dec->pic_order_cnt_cycle_length != 0)
      abs_frame_num = dec->frame_num_offset + dec->frame_num;
    else
      abs_frame_num = 0;

    if (dec->nal_ref_idc == 0 && abs_frame_num > 0)
      abs_frame_num--;

    expected_delta_per_poc_cycle = 0;
    for (int i = 0; i < (int)dec->pic_order_cnt_cycle_length; i++ )
      expected_delta_per_poc_cycle += dec->offset_for_ref_frame[i];

    if (abs_frame_num > 0) {
      const int poc_cycle_cnt = ( abs_frame_num - 1 ) / dec->pic_order_cnt_cycle_length;
      const int frame_num_in_poc_cycle = ( abs_frame_num - 1 ) % dec->pic_order_cnt_cycle_length;

      expected_poc = poc_cycle_cnt * expected_delta_per_poc_cycle;
      for (int i = 0; i <= frame_num_in_poc_cycle; i++)
        expected_poc += dec->offset_for_ref_frame[i];
    } else {
      expected_poc = 0;
    }

    if (dec->nal_ref_idc == 0)
      expected_poc += dec->offset_for_non_ref_pic;

    field_poc[0] = expected_poc + dec->delta_pic_order_cnt[0];
    field_poc[1] = field_poc[0] + dec->offset_for_top_to_bottom_field;

    if (pic_type == AVS_PICTURE_FRAME)
      field_poc[1] += dec->delta_pic_order_cnt[1];

  } else if (dec->pic_order_cnt_type == 2) {
    int poc;
    if (dec->nal_unit_type == AVS_NAL_TYPE_I_PIC_HEADER) {
      poc = 0;
    } else {
      const int abs_frame_num = dec->frame_num_offset + dec->frame_num;
      if (dec->nal_ref_idc != 0)
        poc = 2 * abs_frame_num;
      else
        poc = 2 * abs_frame_num - 1;
    }
    field_poc[0] = poc;
    field_poc[1] = poc;
  }

  /* */
  if (pic_type == AVS_PICTURE_FRAME)
    dec->pic_order_cnt = MIN(field_poc[0], field_poc[1] );
  else if (pic_type == AVS_PICTURE_FIELD_TOP)
    dec->pic_order_cnt = field_poc[0];
  else
    dec->pic_order_cnt = field_poc[1];
}


extern "C" int avs_detect_boundary (const uint8_t *buffer, 
									uint32_t buflen, 
									avs_decode_t *decode)
{
	uint8_t temp;
	avs_decode_t new_decode;
	int ret;
	int slice = 0;
	
	memcpy(&new_decode, decode, sizeof(new_decode));
	
	temp = new_decode.nal_unit_type = avs_nal_unit_type(buffer);
	new_decode.nal_ref_idc = avs_nal_ref_idc(buffer);
	ret = 0;
	switch (temp) {
	case AVS_NAL_TYPE_VIDEO_EXTENSION:
	case AVS_NAL_TYPE_USER_DATA:
	case AVS_NAL_TYPE_VIDEO_EDIT:
#ifdef BOUND_VERBOSE
		printf("nal type %d\n", temp);
#endif
		ret = 0;
		break;
	case AVS_NAL_TYPE_I_SLICE:
		slice = 1;
		// slice buffer - read the info into the new_decode, and compare.
		if (avs_read_slice_info(buffer, buflen, &new_decode) < 0) {
			// need more memory
			return -1;
		}
		if (decode->nal_unit_type != AVS_NAL_TYPE_I_SLICE ) {
			break;
		}
		if (decode->frame_num != new_decode.frame_num) {
#ifdef BOUND_VERBOSE
			printf("frame num values different\n");
#endif
			ret = 1;
			break;
		}
		if (decode->field_pic_flag != new_decode.field_pic_flag) {
			ret = 1;
#ifdef BOUND_VERBOSE
			printf("field pic values different\n");
#endif
			break;
		}
		if (decode->nal_ref_idc != new_decode.nal_ref_idc &&
			(decode->nal_ref_idc == 0 ||
			new_decode.nal_ref_idc == 0)) {
#ifdef BOUND_VERBOSE
			printf("nal ref idc values differ\n");
#endif
			ret = 1;
			break;
		}
		if (decode->frame_num == new_decode.frame_num &&
			decode->pic_order_cnt_type == new_decode.pic_order_cnt_type) {
			if (decode->pic_order_cnt_type == 0) {
				if (decode->pic_order_cnt_lsb != new_decode.pic_order_cnt_lsb) {
#ifdef BOUND_VERBOSE
					printf("pic order 1\n");
#endif
					ret = 1;
					break;
				}
				if (decode->delta_pic_order_cnt_bottom != new_decode.delta_pic_order_cnt_bottom) {
					ret = 1;
#ifdef BOUND_VERBOSE
					printf("delta pic order cnt bottom 1\n");
#endif
					break;
				}
			} else if (decode->pic_order_cnt_type == 1) {
				if (decode->delta_pic_order_cnt[0] != new_decode.delta_pic_order_cnt[0]) {
					ret =1;
#ifdef BOUND_VERBOSE
					printf("delta pic order cnt [0]\n");
#endif
					break;
				}
				if (decode->delta_pic_order_cnt[1] != new_decode.delta_pic_order_cnt[1]) {
					ret = 1;
#ifdef BOUND_VERBOSE
					printf("delta pic order cnt [1]\n");
#endif
					break;
					
				}
			}
		}
		if (decode->nal_unit_type == AVS_NAL_TYPE_I_PIC_HEADER &&
			new_decode.nal_unit_type == AVS_NAL_TYPE_I_PIC_HEADER) {
			if (decode->idr_pic_id != new_decode.idr_pic_id) {
#ifdef BOUND_VERBOSE
				printf("idr_pic id\n");
#endif
				
				ret = 1;
				break;
			}
		}
		break;
	case AVS_NAL_TYPE_SEQ_HEADER:
		if (avs_read_seq_header(buffer, buflen, &new_decode) < 0) {
			return -1;
		}
		ret = 0; //HuoLongshe_20071104
		break;
	case AVS_NAL_TYPE_I_PIC_HEADER:
		new_decode.pic_is_idr = 1;
		if (avs_read_i_pic_hdr_info(buffer, buflen, &new_decode) < 0) {
			return -1;
		}
		ret = 1;  //HuoLongshe_20071104
		break;
	case AVS_NAL_TYPE_P_PIC_HEADER:
	case AVS_NAL_TYPE_B_PIC_HEADER:
		new_decode.pic_is_idr = 0;
		if (avs_read_pb_pic_hdr_info(buffer, buflen, &new_decode) < 0) {
			return -1;
		}
		ret = 1;  //HuoLongshe_20071104
		break;
		
  } 
  
  if (new_decode.nal_unit_type != AVS_NAL_TYPE_I_SLICE) {
	  ret = 1;
  }else{
	  ret = 0;
  }
  
  if( slice ) {  // XXX we compute poc for every slice in a picture (but it's not needed)
	  avs_compute_poc( &new_decode );
  }
  
  
  // other types (6, 7, 8, 
#ifdef BOUND_VERBOSE
  if (ret == 0) {
	  printf("no change\n");
  }
#endif
  
  memcpy(decode, &new_decode, sizeof(*decode));
  
  return ret;
}























