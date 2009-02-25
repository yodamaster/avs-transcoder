#ifndef _SSE_HEADER_H_
#define _SSE_HEADER_H_
/* asm functions */
extern  "C" void  avs_quant_sse(int_32_t qp, int_32_t mode, int_16_t data[8][8]);
extern  "C" void  avs_dequant_sse(int_32_t qp, int_16_t data[8][8]);
extern  "C" void  avs_dct_sse (int_16_t curr_blk[8][8]);
extern  "C" void  avs_idct_sse(int_16_t curr_blk[8][8]);

#endif