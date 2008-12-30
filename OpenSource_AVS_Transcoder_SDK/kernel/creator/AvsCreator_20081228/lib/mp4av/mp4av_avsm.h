
#ifndef __MP4AV_AVSM_H__
#define __MP4AV_AVSM_H__ 1

#define AVSM_START_CODE 0x000001
#define AVSM_PREVENT_3_BYTE 0x000003

#define AVSM_PROFILE_BASELINE 66
#define AVSM_PROFILE_MAIN 77
#define AVSM_PROFILE_EXTENDED 88

#define AVSM_NAL_TYPE_NON_IDR_PIC_HEADER 0x1
#define AVSM_NAL_TYPE_IDR_PIC_HEADER 0x2
#define AVSM_NAL_TYPE_SLICE 0x3
#define AVSM_NAL_TYPE_SEQ_PARAM 0x4
#define AVSM_NAL_TYPE_PIC_PARAM 0x5
#define AVSM_NAL_TYPE_SEI 0x6
#define AVSM_NAL_TYPE_RANDOM_ACCESS_POING 0x7

#define AVSM_TYPE_I 0
#define AVSM_TYPE_P 1

typedef struct avsm_decode_t {
  uint8_t profile;
  uint8_t level;
  uint16_t delta_time_picture_distance_1;
  uint32_t log2_max_frame_num_minus4;
  uint32_t log2_max_pic_order_cnt_lsb_minus4;
  uint32_t pic_order_cnt_type;
  uint8_t frame_mbs_only_flag;
  uint8_t pic_order_present_flag;
  uint8_t delta_pic_order_always_zero_flag;
  int32_t offset_for_non_ref_pic;
  int32_t offset_for_top_to_bottom_field;
  uint32_t pic_order_cnt_cycle_length;
  int16_t offset_for_ref_frame[256];  

  uint8_t picture_coding_type;
  uint8_t picture_distance;
  uint8_t picture_distance_gap_minus1;
  uint32_t frame_num;
  uint8_t pic_is_idr;

  uint8_t nal_ref_idc;
  uint8_t nal_unit_type;

  uint8_t field_pic_flag;
  uint8_t bottom_field_flag;
  uint32_t idr_pic_id;
  uint32_t pic_order_cnt_lsb;
  int32_t delta_pic_order_cnt_bottom;
  int32_t delta_pic_order_cnt[2];

  uint32_t pic_width, pic_height;
  uint32_t slice_type;

  /* POC state */
  int32_t  pic_order_cnt;        /* can be < 0 */

  uint32_t  pic_order_cnt_msb;
  uint32_t  pic_order_cnt_msb_prev;
  uint32_t  pic_order_cnt_lsb_prev;
  uint32_t  frame_num_prev;
  int32_t  frame_num_offset;
  int32_t  frame_num_offset_prev;

} avsm_decode_t;

#ifdef __cplusplus
extern "C" {
#endif

  bool avsm_is_start_code(const uint8_t *pBuf);

uint32_t avsm_find_next_start_code(const uint8_t *pBuf, 
				   uint32_t bufLen);

  uint8_t avsm_nal_unit_type(const uint8_t *buffer);

  int avsm_nal_unit_type_is_slice(uint8_t nal_type);
  int avsm_find_slice_type(const uint8_t *buffer, 
			   uint32_t buflen, 
			   uint8_t *slice_type, 
			   bool noheader);
  uint8_t avsm_nal_ref_idc(const uint8_t *buffer);
  int avsm_detect_boundary(const uint8_t *buffer, 
			   uint32_t buflen, 
			   avsm_decode_t *decode);

  int avsm_read_slice_info(const uint8_t *buffer, 
			   uint32_t buflen, 
			   avsm_decode_t *dec);
  int avsm_read_pic_hdr_info(const uint8_t *buffer, 
			   uint32_t buflen, 
			   avsm_decode_t *dec);
  int avsm_read_seq_info(const uint8_t *buffer, 
			 uint32_t buflen, 
			 avsm_decode_t *dec);
  uint32_t avsm_read_sei_value (const uint8_t *buffer, uint32_t *size);

  bool avsm_picture_is_idr(avsm_decode_t *dec);

  const char *avsm_get_slice_name(uint8_t slice_type);

#ifdef __cplusplus
 }
#endif

#endif
