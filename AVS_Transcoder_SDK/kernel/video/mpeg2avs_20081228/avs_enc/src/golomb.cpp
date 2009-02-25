/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2003, Advanced Audio Video Coding Standard, Part II
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the users without any
* license fee or royalty on an "as is" basis. The AVS disclaims
* any and all warranties, whether express, implied, or statutory,
* including any implied warranties of merchantability or of fitness
* for a particular purpose. In no event shall the contributors or 
* the AVS be liable for any incidental, punitive, or consequential
* damages of any kind whatsoever arising from the use of this program.
*
* This disclaimer of warranty extends to the user of this program
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The AVS does not represent or warrant that the program furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of AVS, including shareware, may be
* subject to royalty fees to patent holders. Information regarding
* the AVS patent policy is available from the AVS Web site at
* http://www.avs.org.cn
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
************************************************************************
*/

/*
*************************************************************************************
* File name: golomb.c
* Function: Description
*
*************************************************************************************
*/


#include <assert.h>
#include "global.h"


/*
*************************************************************************
* Function:
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/

void c_avs_enc::encode_golomb_word(uint_32_t symbol,uint_32_t grad0,uint_32_t max_levels,uint_32_t *res_bits,uint_32_t *res_len)
{
  uint_32_t level,res,numbits;
  
  res=1UL<<grad0;
  level=1UL;numbits=1UL+grad0;

  //find golomb level
  while( symbol>=res && level<max_levels )
  {
    symbol-=res;
    res=res<<1;
    level++;
    numbits+=2UL;
  }

  if(level>=max_levels)
  {
    if(symbol>=res)
      symbol=res-1UL;  //crop if too large.
  }

  //set data bits
  *res_bits=res|symbol;
  *res_len=numbits;
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

void c_avs_enc::encode_multilayer_golomb_word(uint_32_t symbol,const uint_32_t *grad,const uint_32_t *max_levels,uint_32_t *res_bits,uint_32_t *res_len)
{
  unsigned accbits,acclen,bits,len,tmp;

  accbits=acclen=0UL;
  
  while(1)
  {
    encode_golomb_word(symbol,*grad,*max_levels,&bits,&len);
    accbits=(accbits<<len)|bits;
    acclen+=len;
    assert(acclen<=32UL);  //we'l be getting problems if this gets longer than 32 bits.
    tmp=*max_levels-1UL;

    if(!(( len == (tmp<<1)+(*grad) )&&( bits == (1UL<<(tmp+*grad))-1UL )))  //is not last possible codeword? (Escape symbol?)
      break;

    tmp=*max_levels;
    symbol-=(((1UL<<tmp)-1UL)<<(*grad))-1UL;
    grad++;max_levels++;
  }
  *res_bits=accbits;
  *res_len=acclen;
}

int_32_t c_avs_enc::writeSyntaxElement_GOLOMB(SyntaxElement *se, Bitstream *bitstream)
  {
  uint_32_t bits,len,i;
  uint_32_t grad[4],max_lev[4];

  if(!( se->golomb_maxlevels & ~0xFF ))    //only bits 0-7 used? This means normal Golomb word.
    {
    encode_golomb_word(se->value1,se->golomb_grad,se->golomb_maxlevels,&bits,&len);
    }
  else
    {
    for(i=0UL;i<4UL;i++)
      {
      grad[i]=(se->golomb_grad>>(i<<3))&0xFFUL;
      max_lev[i]=(se->golomb_maxlevels>>(i<<3))&0xFFUL;
      }
    encode_multilayer_golomb_word(se->value1,grad,max_lev,&bits,&len);
    }

  se->len=len;
  se->bitpattern=bits;

  writeUVLC2buffer(se, bitstream);

#if TRACE
  snprintf(se->tracestring, TRACESTRING_SIZE, "Coefficients");
  if(se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}