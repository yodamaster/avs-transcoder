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


