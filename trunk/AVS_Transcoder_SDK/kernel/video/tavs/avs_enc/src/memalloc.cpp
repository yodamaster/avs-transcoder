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
#include <stdlib.h>
#include "global.h"
/*
*************************************************************************
* Function:Allocate 2D memory array -> unsigned char array2D[rows][columns]
* Input:
* Output:memory size in bytes
* Return: 
* Attention:
*************************************************************************
*/
void c_avs_enc::no_mem_exit(char *where)
  {
  snprintf(errortext, ET_SIZE, "Could not allocate memory: %s",where);
  error (errortext, 100);
  }
int_32_t c_avs_enc::get_mem2D(byte ***array2D, int_32_t rows, int_32_t columns)
{
  int_32_t i;

  //if((*array2D      = (byte**)calloc(rows,        sizeof(byte*))) == NULL)
  if(((*array2D) = (byte** )_aligned_malloc(rows*sizeof(byte*), 16)) == NULL)
    no_mem_exit("get_mem2D: array2D");
  
  //if(((*array2D)[0] = (byte* )calloc(columns*rows,sizeof(byte ))) == NULL)
  if(((*array2D)[0] = (byte* )_aligned_malloc(rows*columns*sizeof(byte ),16)) == NULL)
    no_mem_exit("get_mem2D: array2D");

  for(i=1;i<rows;i++)
    (*array2D)[i] = (*array2D)[i-1] + columns ;

  return rows*columns;
}

/*
*************************************************************************
* Function:Allocate 2D memory array -> int_32_t array2D[rows][columns]
* Input:
* Output:memory size in bytes
* Return: 
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::get_mem2Dint(int_32_t ***array2D, int_32_t rows, int_32_t columns)
{
  int_32_t i;

  //if((*array2D      = (int_32_t**)calloc(rows,        sizeof(int_32_t*))) == NULL)
  if(((*array2D) = (int_32_t** )_aligned_malloc(rows*sizeof(int_32_t*),16)) == NULL)
    no_mem_exit("get_mem2Dint: array2D");
 // if(((*array2D)[0] = (int_32_t* )calloc(rows*columns,sizeof(int_32_t ))) == NULL)
     if(((*array2D)[0] = (int_32_t* )_aligned_malloc(rows*columns*sizeof(int_32_t ),16)) == NULL)  
    no_mem_exit("get_mem2Dint: array2D");

  for(i=1 ; i<rows ; i++)
    (*array2D)[i] =  (*array2D)[i-1] + columns  ;

  return rows*columns*sizeof(int_32_t);
}
/*
*************************************************************************
* Function:Allocate 3D memory array -> unsigned char array3D[frames][rows][columns]
* Input:
* Output:memory size in bytes
* Return: 
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::get_mem3D(byte ****array3D, int_32_t frames, int_32_t rows, int_32_t columns)
{
  int_32_t  j;
  //if(((*array3D) = (byte***)calloc(frames,sizeof(byte**))) == NULL)
  if(((*array3D) = (byte*** )_aligned_malloc(frames*sizeof(byte**),16)) == NULL)
    no_mem_exit("get_mem3D: array3D");

  for(j=0;j<frames;j++)
    get_mem2D( (*array3D)+j, rows, columns ) ;

  return frames*rows*columns;
}

/*
*************************************************************************
* Function:Allocate 3D memory array -> int_32_t array3D[frames][rows][columns]
* Input:
* Output:memory size in bytes
* Return: 
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::get_mem3Dint(int_32_t ****array3D, int_32_t frames, int_32_t rows, int_32_t columns)
{
  int_32_t  j;

  //if(((*array3D) = (int_32_t***)calloc(frames,sizeof(int_32_t**))) == NULL)
  if(((*array3D) = (int_32_t*** )_aligned_malloc(frames*sizeof(int_32_t**),16)) == NULL)
    no_mem_exit("get_mem3Dint: array3D");

  for(j=0;j<frames;j++)
    get_mem2Dint( (*array3D)+j, rows, columns ) ;

  return frames*rows*columns*sizeof(int_32_t);
}

/*
*************************************************************************
* Function:Allocate 4D memory array -> int_32_t array3D[frames][rows][columns][component]
* Input:
* Output:memory size in bytes
* Return: 
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc::get_mem4Dint(int_32_t *****array4D, int_32_t idx, int_32_t frames, int_32_t rows, int_32_t columns )
{
  int_32_t  j;

  //if(((*array4D) = (int_32_t****)calloc(idx,sizeof(int_32_t**))) == NULL)
  if(((*array4D) = (int_32_t****)_aligned_malloc(idx*sizeof(int_32_t**),16)) == NULL)
    no_mem_exit("get_mem4Dint: array4D");

  for(j=0;j<idx;j++)
    get_mem3Dint( (*array4D)+j, frames, rows, columns ) ;

  return idx*frames*rows*columns*sizeof(int_32_t);
}

/*
*************************************************************************
* Function:free 2D memory array which was alocated with get_mem2D()
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/

void c_avs_enc::free_mem2D(byte **array2D)
{
  if (array2D)
  {
    if (array2D[0])
      _aligned_free (array2D[0]);
    else
      error ("free_mem2D: trying to free unused memory",100);

    _aligned_free (array2D);
  //free (array2D);
  } 
  else
  {
    error ("free_mem2D: trying to free unused memory",100);
  }

}

/*
*************************************************************************
* Function:free 2D memory array
      which was alocated with get_mem2Dint()
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/


void c_avs_enc::free_mem2Dint(int_32_t **array2D)
{
  if (array2D)
  {
    if (array2D[0]) 
      _aligned_free (array2D[0]);
    else
      error ("free_mem2D: trying to free unused memory",100);

     _aligned_free (array2D);

  }
  else
  {
    error ("free_mem2D: trying to free unused memory",100);
  }

}
/*
*************************************************************************
* Function:free 3D memory array
      which was alocated with get_mem3D()
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/

void c_avs_enc::free_mem3D(byte ***array3D, int_32_t frames)
{
  int_32_t i;

  if (array3D)
  {
    for (i=0;i<frames;i++)
    { 
      free_mem2D(array3D[i]);
    }
   _aligned_free (array3D);
  } 
  else
  {
    error ("free_mem3D: trying to free unused memory",100);
  }
}

/*
*************************************************************************
* Function:free 3D memory array 
      which was alocated with get_mem3Dint()
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/

void c_avs_enc::free_mem3Dint(int_32_t ***array3D, int_32_t frames)
{
  int_32_t i;

  if (array3D)
  {
    for (i=0;i<frames;i++)
    { 
      free_mem2Dint(array3D[i]);
    }
   _aligned_free (array3D);
  } 
  else
  {
    error ("free_mem3D: trying to free unused memory",100);
  }
}

/*
*************************************************************************
* Function:free 4D memory array 
      which was alocated with get_mem4Dint()
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/

void c_avs_enc::free_mem4Dint(int_32_t ****array4D, int_32_t idx, int_32_t frames )
{
  int_32_t  j;

  if (array4D)
  {
    for(j=0;j<idx;j++)
      free_mem3Dint( array4D[j], frames) ;
   _aligned_free (array4D);
  }
  else
  {
    error ("free_mem4D: trying to free unused memory",100);
  }

}
/*
*************************************************************************
* Function:Exit program if memory allocation failed (using error())
* Input: where
      string indicating which memory allocation failed
* Output:
* Return: 
* Attention:
*************************************************************************
*/
int_32_t c_avs_enc::get_mem2Dshort_int(int_16_t ***array2D, int_16_t rows, int_16_t columns)
  {
  int_16_t i;

  //if((*array2D      = (int_32_t**)calloc(rows,        sizeof(int_32_t*))) == NULL)
  if(((*array2D) = (int_16_t** )_aligned_malloc(rows*sizeof(int_16_t *),16)) == NULL)
    no_mem_exit("get_mem2Dint: array2D");
  // if(((*array2D)[0] = (int_32_t* )calloc(rows*columns,sizeof(int_32_t ))) == NULL)
  if(((*array2D)[0] = (int_16_t* )_aligned_malloc(rows*columns*sizeof(int_16_t ),16)) == NULL)  
    no_mem_exit("get_mem2Dint: array2D");

  for(i=1 ; i<rows ; i++)
    (*array2D)[i] =  (*array2D)[i-1] + columns  ;

  return rows*columns*sizeof(int_16_t);
  }
void c_avs_enc::free_mem2Dshort_int(int_16_t **array2D)
  {
  if (array2D)
    {
    if (array2D[0]) 
      _aligned_free (array2D[0]);
    else
      error ("free_mem2D: trying to free unused memory",100);

    _aligned_free (array2D);

    }
  else
    {
    error ("free_mem2D: trying to free unused memory",100);
    }
  }
