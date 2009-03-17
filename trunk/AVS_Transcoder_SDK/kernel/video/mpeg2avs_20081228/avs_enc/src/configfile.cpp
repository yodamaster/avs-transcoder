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
#include <string.h>
#include <assert.h>

#include "global.h"
TLS InputParameters configinput;
typedef struct tag_mapping
  {
  char    *TokenName;
  void    *Place;
  int_32_t  Type;
  } Mapping;

TLS Mapping Map[39];
/*
*************************************************************************
* Function:Parse the command line parameters and read the config files.
* Input: ac
         number of command line parameters
      av
        command line parameters
* Output:
* Return: 
* Attention:
*************************************************************************
*/
void c_avs_enc::Configure (char *av)
{
  char *content;  
  Map[0].TokenName  =    "GOPLength";                
  Map[1].TokenName  =    "FramesToBeEncoded";        
  Map[2].TokenName  =    "QPFirstFrame";             
  Map[3].TokenName  =    "QPRemainingFrame";         
  Map[4].TokenName  =    "UseHadamard";              
  Map[5].TokenName  =    "SearchRange";              
  Map[6].TokenName  =    "NumberReferenceFrames";    
  Map[7].TokenName  =    "SourceWidth";              
  Map[8].TokenName  =    "SourceHeight";             
  Map[9].TokenName  =    "InputFile";                
  Map[10].TokenName =    "InputHeaderLength";        
  Map[11].TokenName =    "OutputFile";               
  Map[12].TokenName =    "ReconFile";                
  Map[13].TokenName =    "TraceFile";                
  Map[14].TokenName =    "NumberBFrames";            
  Map[15].TokenName =    "QPBPicture";               
  Map[16].TokenName =    "InterSearch16x16";         
  Map[17].TokenName =    "InterSearch16x8";          
  Map[18].TokenName =    "InterSearch8x16";          
  Map[19].TokenName =    "InterSearch8x8";           
  Map[20].TokenName =    "RDOptimization";           
  Map[21].TokenName =    "InterlaceCodingOption";    
  Map[22].TokenName =    "LoopFilterDisable";        
  Map[23].TokenName =    "LoopFilterParameter";      
  Map[24].TokenName =    "LoopFilterAlphaOffset";    
  Map[25].TokenName =    "LoopFilterBetaOffset";     
  Map[26].TokenName =    "Progressive_frame";        
  Map[27].TokenName =    "Dct_Adaptive_Flag";        
  Map[28].TokenName =    "NumberOfRowsInSlice";      
  Map[29].TokenName =    "SliceParameter";           
  Map[30].TokenName =    "WeightEnable";             
  Map[31].TokenName =    "FrameRate";                
  Map[32].TokenName =    "ChromaFormat";             
  Map[33].TokenName =    "RateControlEnable";        
  Map[34].TokenName =    "Bitrate";                  
  Map[35].TokenName =    "InitialQP";                
  Map[36].TokenName =    "BasicUnit";                
  Map[37].TokenName =    "ChannelType";              
  Map[38].TokenName =    NULL;                       
  Map[0].Place  = &configinput.GopLength;                 
  Map[1].Place  = &configinput.no_frames;                 
  Map[2].Place  = &configinput.qp0;                       
  Map[3].Place  = &configinput.qpN;                       
  Map[4].Place  = &configinput.hadamard;                  
  Map[5].Place  = &configinput.search_range;              
  Map[6].Place  = &configinput.no_multpred;               
  Map[7].Place  = &configinput.img_width;                 
  Map[8].Place  = &configinput.img_height;                
  Map[9].Place  = &configinput.infile;                    
  Map[10].Place = &configinput.infile_header;             
  Map[11].Place = &configinput.outfile;                   
  Map[12].Place = &configinput.ReconFile;                 
  Map[13].Place = &configinput.TraceFile;                 
  Map[14].Place = &configinput.successive_Bframe;         
  Map[15].Place = &configinput.qpB;                       
  Map[16].Place = &configinput.InterSearch16x16;          
  Map[17].Place = &configinput.InterSearch16x8 ;          
  Map[18].Place = &configinput.InterSearch8x16;           
  Map[19].Place = &configinput.InterSearch8x8 ;           
  Map[20].Place = &configinput.rdopt;                     
  Map[21].Place = &configinput.InterlaceCodingOption;     
  Map[22].Place = &configinput.loop_filter_disable;       
  Map[23].Place = &configinput.loop_filter_parameter_flag;
  Map[24].Place = &configinput.alpha_c_offset;            
  Map[25].Place = &configinput.beta_offset;               
  Map[26].Place = &configinput.progressive_frame;         
  Map[27].Place = &configinput.dct_adaptive_flag;         
  Map[28].Place = &configinput.slice_row_nr;              
  Map[29].Place = &configinput.slice_parameter;           
  Map[30].Place = &configinput.picture_weighting_flag;    
  Map[31].Place = &configinput.frame_rate_code;           
  Map[32].Place = &configinput.chroma_format;             
  Map[33].Place = &configinput.RCEnable;                  
  Map[34].Place = &configinput.bit_rate;                  
  Map[35].Place = &configinput.SeinitialQP;               
  Map[36].Place = &configinput.basicunit;                 
  Map[37].Place = &configinput.channel_type;              
  Map[38].Place = NULL;
  Map[0].Type  =  0;
  Map[1].Type  =  0;
  Map[2].Type  =  0;
  Map[3].Type  =  0;
  Map[4].Type  =  0;
  Map[5].Type  =  0;
  Map[6].Type  =  0;
  Map[7].Type  =  0;
  Map[8].Type  =  0;
  Map[9].Type  =  1;
  Map[10].Type =  0;
  Map[11].Type =  1;
  Map[12].Type =  1;
  Map[13].Type =  1;
  Map[14].Type =  0;
  Map[15].Type =  0;
  Map[16].Type =  0;
  Map[17].Type =  0;
  Map[18].Type =  0;
  Map[19].Type =  0;
  Map[20].Type =  0;
  Map[21].Type =  0;
  Map[22].Type =  0;
  Map[23].Type =  0;
  Map[24].Type =  0;
  Map[25].Type =  0;
  Map[26].Type =  0;
  Map[27].Type =  0;
  Map[28].Type =  0;
  Map[29].Type =  0;
  Map[30].Type =  0;
  Map[31].Type =  0;
  Map[32].Type =  0;
  Map[33].Type =  0;
  Map[34].Type =  0;
  Map[35].Type =  0;
  Map[36].Type =  0;
  Map[37].Type =  0;
  Map[38].Type = -1;
  
  memset (&configinput, 0, sizeof (InputParameters));
  
  // Process default config file
  // Parse the command line
  
      content = GetConfigFileContent (av);
      printf ("Parsing Configfile %s", av);
      ParseContent (content, (int_32_t)strlen (content));
      printf ("\n");
      free (content);  
  printf ("\n");
}

/*
*************************************************************************
* Function: Alocate memory buf, opens file Filename in f, reads contents into
        buf and returns buf
* Input:name of config file
* Output:
* Return: 
* Attention:
*************************************************************************
*/


char* c_avs_enc::GetConfigFileContent (char *Filename)
{
  unsigned FileSize;
  FILE *f;
  char *buf;

  if (NULL == (f = fopen (Filename, "r")))
  {
    snprintf (errortext, ET_SIZE, "Cannot open configuration file %s.\n", Filename);
    error (errortext, 300);
  }

  if (0 != fseek (f, 0, SEEK_END))
  {
    snprintf (errortext, ET_SIZE, "Cannot fseek in configuration file %s.\n", Filename);
    error (errortext, 300);
  }

  FileSize = ftell (f);

  if (FileSize < 0 || FileSize > 60000)
  {
    snprintf (errortext, ET_SIZE, "Unreasonable Filesize %d reported by ftell for configuration file %s.\n", FileSize, Filename);
    error (errortext, 300);
  }

  if (0 != fseek (f, 0, SEEK_SET))
  {
    snprintf (errortext, ET_SIZE, "Cannot fseek in configuration file %s.\n", Filename);
    error (errortext, 300);
  }

  if ((buf = (char*) malloc (FileSize + 1))==NULL)
    no_mem_exit("GetConfigFileContent: buf");

  // Note that ftell() gives us the file size as the file system sees it.  The actual file size,
  // as reported by fread() below will be often smaller due to CR/LF to CR conversion and/or
  // control characters after the dos EOF marker in the file.

  FileSize = (unsigned int)fread (buf, 1, FileSize, f);
  buf[FileSize] = '\0';

  fclose (f);

  return buf;
}

/*
*************************************************************************
* Function: Parses the character array buf and writes global variable input, which is defined in
        configfile.h.  This hack will continue to be necessary to facilitate the addition of
       new parameters through the Map[] mechanism (Need compiler-generated addresses in map[]).
* Input:  buf
       buffer to be parsed
     bufsize
       buffer size of buffer
* Output:
* Return: 
* Attention:
*************************************************************************
*/

void c_avs_enc::ParseContent (char *buf, int_32_t bufsize)
{
  char *items[MAX_ITEMS_TO_PARSE];
  int_32_t MapIdx;
  int_32_t item     = 0;
  int_32_t InString = 0;
  int_32_t InItem   = 0;
  char *p      = buf;
  char *bufend = &buf[bufsize];
  int_32_t IntContent;
  int_32_t i;

  // Stage one: Generate an argc/argv-type list in items[], without comments and whitespace.
  // This is context insensitive and could be done most easily with lex(1).

  while (p < bufend)
  {
    switch (*p)
    {
      case 13:
        p++;
        break;
      case '#':                 // Found comment
        *p = '\0';              // Replace '#' with '\0' in case of comment immediately following integer or string
        while (*p != '\n' && p < bufend)  // Skip till EOL or EOF, whichever comes first
          p++;
        InString = 0;
        InItem = 0;
        break;
      case '\n':
        InItem = 0;
        InString = 0;
        *p++='\0';
        break;
      case ' ':
      case '\t':              // Skip whitespace, leave state unchanged
        if (InString)
          p++;
        else
        {                     // Terminate non-strings once whitespace is found
          *p++ = '\0';
          InItem = 0;
        }
        break;
      case '"':               // Begin/End of String
        *p++ = '\0';
        if (!InString)
        {
          items[item++] = p;
          InItem = ~InItem;
        }
        else
          InItem = 0;
        InString = ~InString; // Toggle
        break;
      default:
        if (!InItem)
        {
          items[item++] = p;
          InItem = ~InItem;
        }
        p++;
    }
  }

  item--;

  for (i=0; i<item; i+= 3)
  {
    if (0 > (MapIdx = ParameterNameToMapIndex (items[i])))
    {
      snprintf (errortext, ET_SIZE, " Parsing error in config file: Parameter Name '%s' not recognized.", items[i]);
      error (errortext, 300);
    }
    if (strcmp ("=", items[i+1]))
    {
      snprintf (errortext, ET_SIZE, " Parsing error in config file: '=' expected as the second token in each line.");
      error (errortext, 300);
    }

    // Now interprete the Value, context sensitive...
    switch (Map[MapIdx].Type)
    {
      case 0:           // Numerical
        if (1 != sscanf (items[i+2], "%d", &IntContent))
        {
          snprintf (errortext, ET_SIZE, " Parsing error: Expected numerical value for Parameter of %s, found '%s'.", items[i], items[i+2]);
          error (errortext, 300);
        }
        * (int_32_t *) (Map[MapIdx].Place) = IntContent;
        printf (".");
        break;
      case 1:
        strcpy ((char *) Map[MapIdx].Place, items [i+2]);
        printf (".");
        break;
      default:
        assert ("Unknown value type in the map definition of configfile.h");
    }
  }
  input = &inputs;
  memcpy (input, &configinput, sizeof (InputParameters));

}

/*
*************************************************************************
* Function:Return the index number from Map[] for a given parameter name.
* Input:parameter name string
* Output:
* Return: the index number if the string is a valid parameter name,         \n
          -1 for error
* Attention:
*************************************************************************
*/

int_32_t c_avs_enc:: ParameterNameToMapIndex (char *s)
{
  int_32_t i = 0;

  while (Map[i].TokenName != NULL)
    if (0==strcmp (Map[i].TokenName, s))
      return i;
    else
      i++;
    
  return -1;
};

/*
*************************************************************************
* Function:Checks the input parameters for consistency.
* Input:
* Output:
* Return: 
* Attention:
*************************************************************************
*/

void c_avs_enc::PatchInp ()
  {
  // consistency check of QPs
  input->fixed_picture_qp = 1;
  if (input->qp0 > MAX_QP || input->qp0 < MIN_QP)
    {
    snprintf(errortext, ET_SIZE, "Error input parameter quant_0,check configuration file");
    error (errortext, 400);
    }

  if (input->qpN > MAX_QP || input->qpN < MIN_QP)
    {
    snprintf(errortext, ET_SIZE, "Error input parameter quant_n,check configuration file");
    error (errortext, 400);
    }

  if (input->qpB > MAX_QP || input->qpB < MIN_QP)
    {
    snprintf(errortext, ET_SIZE, "Error input parameter quant_B,check configuration file");
    error (errortext, 400);
    }

  // consistency check no_multpred
  if (input->no_multpred<1) input->no_multpred=1;
  // consistency check size information
  if (input->img_height % 16 != 0 || input->img_width % 16 != 0)
    {
      input->stuff_height = input->img_height;
      input->stuff_width = input->img_width;
    if (input->img_height %16 != 0 )
      {
      input->stuff_height = input->img_height;
      input->img_height = input->img_height + 16 - (input->img_height%16);
      configinput.img_height = input->img_height;
      configinput.stuff_height = input->stuff_height;
      }
    //xzhao 20080709
    if (input->img_width %16 != 0 )
      {
      input->stuff_width = input->img_width;
      input->img_width = input->img_width + 16 - (input->img_width%16);
      configinput.img_width = input->img_width;
      configinput.stuff_width = input->stuff_width;
      }
    }
  else
    {
    //if(input->img_height-input->stuff_height>16)
      input->stuff_height = input->img_height;
    //if(input->img_width-input->stuff_width>16)
      input->stuff_width = input->img_width;
    }
  // check range of filter offsets
  if (input->alpha_c_offset > 6 || input->alpha_c_offset < -6)
    {
    snprintf(errortext, ET_SIZE, "Error input parameter LFAlphaC0Offset, check configuration file");
    error (errortext, 400);
    }

  if (input->beta_offset > 6 || input->beta_offset < -6)
    {
    snprintf(errortext, ET_SIZE, "Error input parameter LFBetaOffset, check configuration file");
    error (errortext, 400);
    }
  // Open Files
  if (strlen (input->infile) > 0 && (p_in=fopen(input->infile,"rb"))==NULL)
    {
    snprintf(errortext, ET_SIZE, "Input file %s does not exist",input->infile);    
    }

  if (strlen (input->ReconFile) > 0 && (p_rec=fopen(input->ReconFile, "wb"))==NULL)
    {
    snprintf(errortext, ET_SIZE, "Error open file %s", input->ReconFile);    
    }
#ifdef _OUTPUT_DEC_IMG_
  if (strlen (input->DecRecFile) > 0 && (p_org_dec=fopen(input->DecRecFile, "wb"))==NULL)
    {
    snprintf(errortext, ET_SIZE, "Error open file %s", input->DecRecFile);    
    }
#endif
  if (strlen (input->TraceFile) > 0 && (p_trace=fopen(input->TraceFile,"w"))==NULL)
    {
    snprintf(errortext, ET_SIZE, "Error open file %s", input->TraceFile);    
    }

  if (input->slice_row_nr==0)
    {
    input->slice_row_nr=input->img_height/16;
    }  
  // Set block sizes
  input->blc_size[0][0]=16;
  input->blc_size[0][1]=16;
  input->blc_size[1][0]=16;
  input->blc_size[1][1]=16;
  input->blc_size[2][0]=16;
  input->blc_size[2][1]= 8;
  input->blc_size[3][0]= 8;
  input->blc_size[3][1]=16;
  input->blc_size[4][0]= 8;
  input->blc_size[4][1]= 8;
  input->blc_size[8][0]= 8;
  input->blc_size[8][1]= 8;  
  if (input->slice_row_nr==0)
    {
    input->slice_row_nr=input->img_height/16;
    }
  }

