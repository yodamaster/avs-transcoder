
#ifndef _DEFINES_H_
#define _DEFINES_H_
//#define ROI_ENABLE
#define AVS
#define TRACE  0  /* !< 0:Trace off 1:Trace on */
#define snprintf  _snprintf
#define _THREE_STEP_MOTION_SEARCH_
//#define _OUTPUT_DEC_IMG_
//#define _ME_FOR_RATE_CONTROL_
/*
 * define FastME ;
 * #define FIELDINTE
 */
#define AVS_OUT_BUFFER_SIZE  (1024 * 1024 * 4)

#define ET_SIZE        300  /* !< size of error text buffer */
#define NUM_PIC_TYPE      5

#define max(a, b)      (((a) > (b)) ? (a) : (b))
#define min(a, b)      (((a) < (b)) ? (a) : (b))
#define clamp(a, b, c)    ((a) < (b) ? (b) : ((a) > (c) ? (c) : (a)))  /* !< clamp a to the range of [b; c] */
#define Clip3(min, max, val)  (((val) < (min)) ? (min) : (((val) > (max)) ? (max) : (val)))

#define MAXMODE        13

#define B8_SIZE 8
#define MB_BLOCK_SIZE 16
#define BLOCK_SIZE 4
#define BLOCK_MULTIPLE (MB_BLOCK_SIZE / (2 * BLOCK_SIZE))

#define INIT_FRAME_RATE    30
#define EOSEOS      1    /* !< End Of Sequence */

#define MAX_SYMBOLS_PER_MB  1200    /* !< Maximum number of different syntax elements for one MB */
#define SVA_STREAM_BUF_SIZE 1024 
#endif
