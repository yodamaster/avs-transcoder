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
#include "mp4av_avsm.h"
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

uint32_t avsm_ue (CBitstream *bs)
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

int32_t avsm_se (CBitstream *bs) 
{
  uint32_t ret;
  ret = avsm_ue(bs);
  if ((ret & 0x1) == 0) {
    ret >>= 1;
    int32_t temp = 0 - ret;
    return temp;
  } 
  return (ret + 1) >> 1;
}

static void avsm_decode_annexb( uint8_t *dst, int *dstlen,
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

extern "C" bool avsm_is_start_code (const uint8_t *pBuf) 
{
  if (pBuf[0] == 0 && 
      pBuf[1] == 0 && 
      ((pBuf[2] == 1) ||
       ((pBuf[2] == 0) && pBuf[3] == 1))) {
    return true;
  }
  return false;
}

extern "C" uint32_t avsm_find_next_start_code (const uint8_t *pBuf, 
					       uint32_t bufLen)
{
  uint32_t val, temp;
  uint32_t offset;

  offset = 0;
  if (pBuf[0] == 0 && 
      pBuf[1] == 0 && 
      ((pBuf[2] == 1) ||
       ((pBuf[2] == 0) && pBuf[3] == 1))) {
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
    if (val == AVSM_START_CODE) {
      if (temp == 0) return offset - 4;
      return offset - 3;
    }
  }
  return 0;
}

extern "C" uint8_t avsm_nal_unit_type (const uint8_t *buffer)
{
  uint32_t offset;
  if (buffer[2] == 1) offset = 3;
  else offset = 4;
  return buffer[offset] & 0x1f;
}

extern "C" int avsm_nal_unit_type_is_slice (uint8_t type)
{
  if (type == AVSM_NAL_TYPE_SLICE ) {
    return true;
  }
  return false;
}

/*
 * determine if the slice we decoded is a sync point
 */
extern "C" bool avsm_picture_is_idr (avsm_decode_t *dec) 
{
  if (!dec->pic_is_idr)
    return false;
  if (dec->picture_coding_type == AVSM_TYPE_I) return true;
  return false;
}

extern "C" uint8_t avsm_nal_ref_idc (const uint8_t *buffer)
{
  uint32_t offset;
  if (buffer[2] == 1) offset = 3;
  else offset = 4;
  return (buffer[offset] >> 5) & 0x3;
}

int avsm_read_seq_info (const uint8_t *buffer, 
			uint32_t buflen, 
			avsm_decode_t *dec)
{
  CBitstream bs;
  uint32_t header;
  uint8_t tmp[2048]; /* Should be enough for all SPS (we have at worst 13 bytes and 496 se/ue in frext) */
  int tmp_len;

  if (buffer[2] == 1) header = 4;
  else header = 5;

  avsm_decode_annexb( tmp, &tmp_len, buffer + header, MIN(buflen-header,2048) );
  bs.init(tmp, tmp_len * 8);

  //bs.set_verbose(true);
  try {
    dec->profile = bs.GetBits(8);
    dec->level = bs.GetBits(8);
    avsm_ue(&bs); // seq_parameter_set_id
    dec->delta_time_picture_distance_1 = bs.GetBits(16);
    avsm_ue(&bs); // num_ref_frames
    uint32_t PicWidthInMbs = avsm_ue(&bs) + 1;
    dec->pic_width = PicWidthInMbs * 16;
    uint32_t PicHeightInMapUnits = avsm_ue(&bs) + 1;
    dec->pic_height = PicHeightInMapUnits * 16;
#if 0
    
    printf("   aspect_ratio: %u\n", bs->GetBits(4));
    temp = bs->GetBits(1);
    printf("   frame_cropping_flag: %u\n", temp);
    if (temp) {
      printf("     frame_crop_left_offset: %u\n", avsm_ue(bs));
      printf("     frame_crop_right_offset: %u\n", avsm_ue(bs));
      printf("     frame_crop_top_offset: %u\n", avsm_ue(bs));
      printf("     frame_crop_bottom_offset: %u\n", avsm_ue(bs));
    }
    temp = bs->GetBits(1);
    printf("   hrd_parameters_present_flag: %u\n", temp);
    if (temp) {
      avsm_hrd_parameters(bs);
    }
#endif
  } catch (...) {
    return -1;
  }
  return 0;
}

int avsm_read_pic_hdr_info (const uint8_t *buffer, 
					   uint32_t buflen, 
					   avsm_decode_t *dec)
{
	CBitstream bs;
	uint32_t header;
	uint8_t tmp[1024]; /* Should be enough for all SPS (we have at worst 13 bytes and 496 se/ue in frext) */
	int tmp_len;
	
	if (buffer[2] == 1) header = 4;
	else header = 5;
	
	avsm_decode_annexb( tmp, &tmp_len, buffer + header, MIN(buflen-header,2048) );
	bs.init(tmp, tmp_len * 8);
	
	//bs.set_verbose(true);
	try {
		dec->picture_coding_type = bs.GetBits(2);
		dec->picture_distance = bs.GetBits(8);
		avsm_ue(&bs); // pic_parameter_set_id
		if (dec->nal_unit_type == AVSM_NAL_TYPE_IDR_PIC_HEADER)
			dec->picture_distance_gap_minus1 = avsm_ue(&bs);
		dec->frame_num = bs.GetBits(5);
		
	} catch (...) {
		return -1;
	}
	return 0;
}

extern "C" int avsm_find_slice_type (const uint8_t *buffer, 
				     uint32_t buflen,
				     uint8_t *slice_type, 
				     bool noheader)
{
  uint32_t header;
  if (noheader) header = 1;
  else {
    if (buffer[2] == 1) header = 4;
    else header = 5;
  }
  CBitstream bs;
  bs.init(buffer + header, (buflen - header) * 8);
  try {
    avsm_ue(&bs); // first_mb_in_slice
    *slice_type = avsm_ue(&bs); // slice type
  } catch (...) {
    return -1;
  }
  return 0;
}

int avsm_read_slice_info (const uint8_t *buffer, 
			  uint32_t buflen, 
			  avsm_decode_t *dec)
{
  uint32_t header;
  uint8_t tmp[512]; /* Enough for the begining of the slice header */
  int tmp_len;

  if (buffer[2] == 1) header = 4;
  else header = 5;
  CBitstream bs;

  avsm_decode_annexb( tmp, &tmp_len, buffer + header, MIN(buflen-header,512) );
  bs.init(tmp, tmp_len * 8);
  try {
    avsm_ue(&bs); // first_mb_in_slice
  } catch (...) {
    return -1;
  }
  return 0;
}

static void avsm_compute_poc( avsm_decode_t *dec ) {
  const int max_frame_num = 1 << (dec->log2_max_frame_num_minus4 + 4);
  int field_poc[2] = {0,0};
  enum {
    AVSM_PICTURE_FRAME,
    AVSM_PICTURE_FIELD_TOP,
    AVSM_PICTURE_FIELD_BOTTOM,
  } pic_type;

  /* FIXME FIXME it doesn't handle the case where there is a MMCO == 5
   * (MMCO 5 "emulates" an idr) */
  
  /* picture type */
  if (dec->frame_mbs_only_flag || !dec->field_pic_flag)
    pic_type = AVSM_PICTURE_FRAME;
  else if (dec->bottom_field_flag)
    pic_type = AVSM_PICTURE_FIELD_BOTTOM;
  else
    pic_type = AVSM_PICTURE_FIELD_TOP;

  /* frame_num_offset */
  if (dec->nal_unit_type == AVSM_NAL_TYPE_IDR_PIC_HEADER) {
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
    if (pic_type == AVSM_PICTURE_FRAME)
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

    if (pic_type == AVSM_PICTURE_FRAME)
      field_poc[1] += dec->delta_pic_order_cnt[1];

  } else if (dec->pic_order_cnt_type == 2) {
    int poc;
    if (dec->nal_unit_type == AVSM_NAL_TYPE_IDR_PIC_HEADER) {
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
  if (pic_type == AVSM_PICTURE_FRAME)
    dec->pic_order_cnt = MIN(field_poc[0], field_poc[1] );
  else if (pic_type == AVSM_PICTURE_FIELD_TOP)
    dec->pic_order_cnt = field_poc[0];
  else
    dec->pic_order_cnt = field_poc[1];
}


extern "C" int avsm_detect_boundary (const uint8_t *buffer, 
				     uint32_t buflen, 
				     avsm_decode_t *decode)
{
  uint8_t temp;
  avsm_decode_t new_decode;
  int ret;
  int slice = 0;

  memcpy(&new_decode, decode, sizeof(new_decode));

  temp = new_decode.nal_unit_type = avsm_nal_unit_type(buffer);
  new_decode.nal_ref_idc = avsm_nal_ref_idc(buffer);
  ret = 0;
  switch (temp) {
  case AVSM_NAL_TYPE_RANDOM_ACCESS_POING:
  case AVSM_NAL_TYPE_SEI:
#ifdef BOUND_VERBOSE
    printf("nal type %d\n", temp);
#endif
    ret = 0;
    break;
  case AVSM_NAL_TYPE_SLICE:
    slice = 1;
    // slice buffer - read the info into the new_decode, and compare.
    if (avsm_read_slice_info(buffer, buflen, &new_decode) < 0) {
      // need more memory
      return -1;
    }
    if (decode->nal_unit_type != AVSM_NAL_TYPE_SLICE ) {
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
    if (decode->nal_unit_type == AVSM_NAL_TYPE_IDR_PIC_HEADER &&
	new_decode.nal_unit_type == AVSM_NAL_TYPE_IDR_PIC_HEADER) {
      if (decode->idr_pic_id != new_decode.idr_pic_id) {
#ifdef BOUND_VERBOSE
	printf("idr_pic id\n");
#endif
	
	ret = 1;
	break;
      }
    }
    break;
  case AVSM_NAL_TYPE_SEQ_PARAM:
    if (avsm_read_seq_info(buffer, buflen, &new_decode) < 0) {
      return -1;
    }
	break;
  case AVSM_NAL_TYPE_IDR_PIC_HEADER:
	new_decode.pic_is_idr = 1;
	if (avsm_read_pic_hdr_info(buffer, buflen, &new_decode) < 0) {
		return -1;
	}
	break;
  case AVSM_NAL_TYPE_NON_IDR_PIC_HEADER:
	new_decode.pic_is_idr = 0;
	if (avsm_read_pic_hdr_info(buffer, buflen, &new_decode) < 0) {
		return -1;
	}
	break;
    
  } 

  if (new_decode.nal_unit_type != AVSM_NAL_TYPE_SLICE) {
	  ret = 1;
  }else{
	  ret = 0;
  }

  

  /* save _prev values */
  if (ret)
  {
    /*new_decode.frame_num_offset_prev = decode->frame_num_offset;
    if (decode->pic_order_cnt_type != 2 || decode->nal_ref_idc != 0)
      new_decode.frame_num_prev = decode->frame_num;
    if (decode->nal_ref_idc != 0)
    {
      new_decode.pic_order_cnt_lsb_prev = decode->pic_order_cnt_lsb;
      new_decode.pic_order_cnt_msb_prev = decode->pic_order_cnt_msb;
    }*/
  }

  if( slice ) {  // XXX we compute poc for every slice in a picture (but it's not needed)
    avsm_compute_poc( &new_decode );
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

uint32_t avsm_read_sei_value (const uint8_t *buffer, uint32_t *size) 
{
  uint32_t ret = 0;
  *size = 1;
  while (buffer[*size] == 0xff) {
    ret += 255;
    *size = *size + 1;
  }
  ret += *buffer;
  return ret;
}

/*extern "C" const char *avsm_get_slice_name (uint8_t slice_type)
{
  if (AVSM_TYPE_IS_P(slice_type)) return "P";  
  if (AVSM_TYPE_IS_B(slice_type)) return "B";  
  if (AVSM_TYPE_IS_I(slice_type)) return "I";
  if (AVSM_TYPE_IS_SI(slice_type)) return "SI";
  if (AVSM_TYPE_IS_SP(slice_type)) return "SP";
  return "UNK";
}*/





















