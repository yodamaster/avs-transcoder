/*
*****************************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2003, Advanced Audio Video Coding Standard, Part II
*
* DISCLAIMER OF WARRANTY
*
* The contents of this file are subject to the Mozilla Public License
* Version 1.1 (the "License"); you may not use this file except in
* compliance with the License. You may obtain a copy of the License at
* http://www.mozilla.org/MPL/
*
* Software distributed under the License is distributed on an "AS IS"
* basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
* License for the specific language governing rights and limitations under
* the License.
*                     
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
* The AVS Working Group doesn't represent or warrant that the programs
* furnished here under are free of infringement of any third-party patents.
* Commercial implementations of AVS, including shareware, may be
* subject to royalty fees to patent holders. Information regarding
* the AVS patent policy for standardization procedure is available at 
* AVS Web site http://www.avs.org.cn. Patent Licensing is outside
* of AVS Working Group.
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE AVS PATENT POLICY.
************************************************************************
*/

/*
*************************************************************************************
* File name: 
* Function: 
*
*************************************************************************************
*/
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "global.h"


/*
*************************************************************************
* Function:Reference buffer read, Full pel
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/

pel_t* c_avs_enc:: FastLineX (int_32_t dummy, pel_t* Pic, int_32_t y, int_32_t x)
{
  return Pic + y*img->width + x;
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
pel_t* c_avs_enc:: UMVLineX (int_32_t size, pel_t* Pic, int_32_t y, int_32_t x)
{
  int_32_t i, maxx;
  pel_t *Picy;

  Picy = Pic + max(0,min(img->height-1,y)) * img->width;

  if (x < 0)                            // Left edge
  {
    maxx = min(0,x+size);

    for (i = x; i < maxx; i++)
    {
      line[i-x] = Picy [0];             // Replicate left edge pixel
    }

    maxx = x+size;

    for (i = 0; i < maxx; i++)          // Copy non-edge pixels
      line[i-x] = Picy [i];
  }
  else if (x > img->width-size)         // Right edge
  {
    maxx = img->width;

    for (i = x; i < maxx; i++)
    {
      line[i-x] = Picy [i];             // Copy non-edge pixels
    }

    maxx = x+size;

    for (i = max(img->width,x); i < maxx; i++)
    {
      line[i-x] = Picy [img->width-1];  // Replicate right edge pixel
    }
  }
  else                                  // No edge
  {
    return Picy + x;
  }

  return line;
}

/*
*************************************************************************
* Function:Reference buffer, 1/4 pel
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/

pel_t c_avs_enc:: UMVPelY_14 (pel_t **Pic, int_32_t y, int_32_t x)
{
  int_32_t width4  = ((img->width+2*IMG_PAD_SIZE-1)<<2);
  int_32_t height4 = ((img->height+2*IMG_PAD_SIZE-1)<<2);
  x += IMG_PAD_SIZE << 2;
  y += IMG_PAD_SIZE << 2;
  
  if (x < 0)
    x &= 3;
  if (x > width4)
  {
    x &= 3;
    x += width4;
  }
  if (y < 0)
    y &= 3;
  if (y > height4)
  {
    y &= 3;
    y += height4;
  }

  return Pic [y][x];
}

pel_t c_avs_enc:: FastPelY_14 (pel_t **Pic, int_32_t y, int_32_t x)
{
  return Pic [IMG_PAD_SIZE*4+y][IMG_PAD_SIZE*4+x];
}


