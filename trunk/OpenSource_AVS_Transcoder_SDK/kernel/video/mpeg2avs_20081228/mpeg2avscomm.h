#ifndef _MPEG2AVSCOMM_H_
#define _MPEG2AVSCOMM_H_
#include "typedef.h"
typedef struct
	{
	int_32_t mb_type;				// 0: skip, 1: intra, 2:intra

	int_32_t pdir;					// 0: foreward, 1: backward, 2: bidirect, -1: intra

	int_32_t mc_type;			// must be 	MC_FRAME = 2;
	int_32_t mv[2][2][2];		// mpeg2mv[mba][r][s][t]
	}MB_INFO;
	float global_frame_encoding_time[4];//0 is null, 1 is I, 2 is P, 3 is B
#endif
