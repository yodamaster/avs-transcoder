
#ifndef __MP4AV_AVS_H__
#define __MP4AV_AVS_H__ 1

#define AVS_START_CODE 0x000001
#define AVS_PREVENT_3_BYTE 0x000003

#define AVS_PROFILE_BASELINE 32

#define AVS_NAL_TYPE_SEQ_HEADER 1
#define AVS_NAL_TYPE_VIDEO_EXTENSION 2
#define AVS_NAL_TYPE_USER_DATA 3
#define AVS_NAL_TYPE_VIDEO_EDIT 4
#define AVS_NAL_TYPE_I_PIC_HEADER 5
#define AVS_NAL_TYPE_P_PIC_HEADER 6
#define AVS_NAL_TYPE_B_PIC_HEADER 7
#define AVS_NAL_TYPE_I_SLICE 8
#define AVS_NAL_TYPE_P_SLICE 9
#define AVS_NAL_TYPE_B_SLICE 10

#define AVS_TYPE_I 0
#define AVS_TYPE_P 1
#define AVS_TYPE_B 2

typedef struct avs_decode_t {
	uint8_t profile;
	uint8_t level;
	uint8_t progressive_sequence;
	uint32_t horizontal_size, vertical_size;
	double frame_rate;
	
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
	uint8_t last_nal_unit_type;
	
	uint8_t field_pic_flag;
	uint8_t bottom_field_flag;
	uint32_t idr_pic_id;
	uint32_t pic_order_cnt_lsb;
	int32_t delta_pic_order_cnt_bottom;
	int32_t delta_pic_order_cnt[2];
	
	
	uint32_t slice_type;
	
	/* POC state */
	int32_t  pic_order_cnt;        /* can be < 0 */
	
	uint32_t  pic_order_cnt_msb;
	uint32_t  pic_order_cnt_msb_prev;
	uint32_t  pic_order_cnt_lsb_prev;
	uint32_t  frame_num_prev;
	int32_t  frame_num_offset;
	int32_t  frame_num_offset_prev;
	
} avs_decode_t;

#ifdef __cplusplus
extern "C" {
#endif
	
	bool avs_is_start_code(const uint8_t *pBuf);
	
	uint32_t avs_find_next_start_code(const uint8_t *pBuf, 
		uint32_t bufLen);
	
	uint8_t avs_nal_unit_type(const uint8_t *buffer);
	
	int avs_nal_unit_type_is_slice(uint8_t nal_type);
	int avs_find_slice_type(const uint8_t *buffer, 
		uint32_t buflen, 
		uint8_t *slice_type, 
		bool noheader);
	uint8_t avs_nal_ref_idc(const uint8_t *buffer);
	int avs_detect_boundary(const uint8_t *buffer, 
		uint32_t buflen, 
		avs_decode_t *decode);
	
	int avs_read_slice_info(const uint8_t *buffer, 
		uint32_t buflen, 
		avs_decode_t *dec);
	int avs_read_i_pic_hdr_info(const uint8_t *buffer, 
		uint32_t buflen, 
		avs_decode_t *dec);
	int avs_read_pb_pic_hdr_info(const uint8_t *buffer, 
		uint32_t buflen, 
		avs_decode_t *dec);
	int avs_read_seq_header(const uint8_t *buffer, 
		uint32_t buflen, 
		avs_decode_t *dec);
	uint32_t avs_read_sei_value (const uint8_t *buffer, uint32_t *size);
	
	bool avs_picture_is_idr(avs_decode_t *dec);
	
	const char *avs_get_slice_name(uint8_t slice_type);
	
#ifdef __cplusplus
}
#endif

#endif
