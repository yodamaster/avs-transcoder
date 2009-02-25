/*$T typedef.h GC 1.140 10/28/07 15:50:20 */


/*$6
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */


#ifndef _TYPEDEF_H_
#define _TYPEDEF_H_
#include "disablewarning.h"
typedef unsigned char  byte;  /* !< byte type definition */
#define pel_t  byte
typedef char    int_8_t;
typedef short int  int_16_t;
typedef int        int_32_t;
typedef unsigned short  uint_16_t;
typedef unsigned int  uint_32_t;
#define TLS         __declspec(thread)
#endif
