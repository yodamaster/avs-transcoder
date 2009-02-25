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
* File name: rdopt_coding_state.c
* Function:  
    Storing/restoring coding state for
    Rate-Distortion optimized mode decision
*
*************************************************************************************
*/


#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "global.h"

void c_avs_enc:: delete_coding_state (CSptr cs)
{
  if (cs != NULL)
  {
    if (cs->bitstream != NULL)
      free (cs->bitstream);

    //=== coding state structure ===
    free (cs);
    cs=NULL;
  }

}

CSptr c_avs_enc:: create_coding_state ()
{
  CSptr cs;

  //=== coding state structure ===
  if ((cs = (CSptr) calloc (1, sizeof(CSobj))) == NULL)
    no_mem_exit("init_coding_state: cs");

  //=== important variables of data partition array ===

  if ((cs->bitstream = (Bitstream*) calloc (1, sizeof(Bitstream))) == NULL)
    no_mem_exit("init_coding_state: cs->bitstream");

  return cs;
}

void c_avs_enc::store_coding_state (CSptr cs)
{
  Bitstream            *bs_src, *bs_dest;
  Macroblock           *currMB  = &(img->mb_data [img->current_mb_nr]);

  if (!input->rdopt)
    return;

  //=== important variables of data partition array ===
  bs_src  =  currBitStream;
  bs_dest = cs->bitstream;

  memcpy (bs_dest, bs_src, sizeof(Bitstream));

  //=== syntax element number and bitcounters ===
  cs->currSEnr = currMB->currSEnr;
  memcpy (cs->bitcounter, currMB->bitcounter, MAX_BITCOUNTER_MB*sizeof(int_32_t));

  //=== elements of current macroblock ===
  memcpy (cs->mvd, currMB->mvd, 2*2*BLOCK_MULTIPLE*BLOCK_MULTIPLE*sizeof(int_32_t));
  cs->cbp_bits = currMB->cbp_bits;
}

void c_avs_enc:: reset_coding_state (CSptr cs)
{
  Bitstream            *bs_src, *bs_dest;
  Macroblock           *currMB  = &(img->mb_data [img->current_mb_nr]);

  if (!input->rdopt)
    return;

  bs_dest =   currBitStream;
  bs_src  =   cs->bitstream;

  memcpy (bs_dest, bs_src, sizeof(Bitstream));

  //=== syntax element number and bitcounters ===
  currMB->currSEnr = cs->currSEnr;
  memcpy (currMB->bitcounter, cs->bitcounter, MAX_BITCOUNTER_MB*sizeof(int_32_t));

  //=== elements of current macroblock ===
  memcpy (currMB->mvd, cs->mvd, 2*2*BLOCK_MULTIPLE*BLOCK_MULTIPLE*sizeof(int_32_t));
  currMB->cbp_bits = cs->cbp_bits;

}

