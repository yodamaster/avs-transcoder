/*
* The contents of this file are subject to the Mozilla Public
* License Version 1.1 (the "License"); you may not use this file
* except in compliance with the License. You may obtain a copy of
* the License at http://www.mozilla.org/MPL/
* df
* Software distributed under the License is distributed on an "AS
* IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
* implied. See the License for the specific language governing
* rights and limitations under the License.
* 
* The Original Code is MPEG4IP.
* 
* The Initial Developer of the Original Code is Cisco Systems Inc.
* Portions created by Cisco Systems Inc. are
* Copyright (C) Cisco Systems Inc. 2001-2004.  All Rights Reserved.
* 
* Portions created by Ximpo Group Ltd. are
* Copyright (C) Ximpo Group Ltd. 2003, 2004.  All Rights Reserved.
*
* Contributor(s): 
*		Dave Mackie			dmackie@cisco.com
*		Alix Marchandise-Franquet	alix@cisco.com
*		Ximpo Group Ltd.		mp4v2@ximpo.com
*/

#define MP4CREATOR_GLOBALS
#include "mp4creator.h"
#include "mpeg4ip_getopt.h"

#include "stdio.h"

// forward declarations
// AMR defines
#define AMR_TYPE_NONE 0
#define AMR_TYPE_AMR 1
#define AMR_TYPE_AMRWB 2

#define AMR_MAGIC_LEN_AMR 6
#define AMR_MAGIC_AMR "#!AMR\n"

#define AMR_MAGIC_LEN_AMRWB 9
#define AMR_MAGIC_AMRWB "#!AMR-WB\n"

MP4TrackId* CreateMediaTracks(
							  MP4FileHandle mp4File, 
							  const char* inputFileName,
							  bool doEncrypt,
							  bool doForcetimestamp=false);//Lirh20071113

void CreateHintTrack(
					 MP4FileHandle mp4File, 
					 MP4TrackId mediaTrackId,
					 const char* payloadName, 
					 bool interleave, 
					 u_int16_t maxPayloadSize,
					 bool doEncrypt);

void ExtractTrack(
				  MP4FileHandle mp4File, 
				  MP4TrackId trackId, 
				  const char* outputFileName);

void ExtractTimestamp( const char* mp4FileName );

bool ExtractAvsSlice( u_int8_t *source , u_int32_t sourceSize , u_int8_t*dest , u_int32_t*destSize ) ;
void load_avs_start_code( u_int8_t *dest , u_int32_t index ) ;
void adjust_Len_Pos( u_int8_t *pSample ) ;
// external declarations

// track creators
MP4TrackId Mp3Creator(MP4FileHandle mp4File, FILE* inFile, bool doEncrypt);

MP4TrackId AVSMCreator(MP4FileHandle mp4File, FILE *inFile);
MP4TrackId AVSCreator(MP4FileHandle mp4File, FILE *inFile, FILE* datFile=NULL);//Lirh20071113



// main routine
int main(int argc, char** argv)
{
	char* usageString = 
		" <options> <asm-file>\n"
		"  Options:\n"
		"  -create=<input-file>    Create track from <input-file>\n"
		"  -encrypt[=<track-id>]   Encrypt a track, also -E\n"
		"  -extract=<track-id>     Extract a track\n"
		"  -delete=<track-id>      Delete a track\n"
		"  -forcetimestamp         Use the timestamp file when creating video track\n"
		"  -gettimestamp           Get the timestamp from the asm file\n"
		"  -hint[=<track-id>]      Create hint track, also -H\n"
		"  -list                   List tracks in asm file\n"
		"  -mtu=<size>             Maximum Payload size for RTP packets in hint track\n"
		"  -optimize               Optimize asm file layout to interleaving audio and video, also -O\n"
		"  -rate=<fps>             Video frame rate, e.g. 30 or 29.97\n"
		"  -timescale=<ticks>      Time scale (ticks per second)\n"
		"  -use64bits              Use for large files\n"
		"  -use64bitstime          Use for 64 Bit times (not QT player compatible)\n"
		"  -verbose[=[1-5]]        Enable debug messages\n"
		"  -version                Display version information\n"
		"\n"
		"  北京大学数字媒体研究所 V1.8\n\n"
		;
	
	bool doCreate = false;
	bool doEncrypt = false;
	bool doExtract = false;
	bool doDelete = false;
	bool doForcetimestamp = false;//Lirh20071113
	bool doGettimestamp = false;//Lirh20080528
	bool doHint = false;
	bool doList = false;
	bool doOptimize = false;
	bool doIsma = false;
	uint64_t createFlags = 0;
	char* mp4FileName = NULL;
	char* inputFileName = NULL;
	char* outputFileName = NULL;
	char* payloadName = NULL;
	MP4TrackId hintTrackId = MP4_INVALID_TRACK_ID;
	MP4TrackId encryptTrackId = MP4_INVALID_TRACK_ID;
	MP4TrackId extractTrackId = MP4_INVALID_TRACK_ID;
	MP4TrackId deleteTrackId = MP4_INVALID_TRACK_ID;
	//u_int16_t maxPayloadSize = 1460;  //For UDP
	u_int16_t maxPayloadSize = 1440;  //For TCP, 1500-20(IP)-20(TCP)-12(RTP)-8(other)  --HuoLongshe

	Verbosity = MP4_DETAILS_ERROR;
	VideoFrameRate = 0;		// determine from input file
	TimeScaleSpecified = false;
	Mp4TimeScale = 90000;
	VideoProfileLevelSpecified = false;
	
	// begin processing command line
	ProgName = argv[0];					                               	
	
	while (true) {
		int c = -1;
		int option_index = 0;
		static struct option long_options[] = {
			{ "create", 1, 0, 'c' },
			{ "delete", 1, 0, 'd' },
			{ "extract", 2, 0, 'e' },
			{ "encrypt", 2, 0, 'E' },
			{ "forcetimestamp", 0, 0, 'f' },//Lirh20071113
			{ "gettimestamp", 0, 0, 'g' },//Lirh20080528
			{ "help", 0, 0, '?' },
			{ "hint", 2, 0, 'H' },
			{ "list", 0, 0, 'l' },
			{ "mtu", 1, 0, 'm' },
			{ "optimize", 0, 0, 'O' },
			{ "rate", 1, 0, 'r' },
			{ "timescale", 1, 0, 't' },
			{ "use64bits", 0, 0, 'u' },
			{ "use64bitstime", 0, 0, 'U' },
			{ "verbose", 2, 0, 'v' },
			{ "version", 0, 0, 'V' },
			{ NULL, 0, 0, 0 }
		};
		
		c = getopt_long_only(argc, argv, "aBc:Cd:e:E::fgGH::iIlL:m:Op:P:r:t:T:uUv::VZ",//Lirh20071113//Lirh20080528
			long_options, &option_index);												
		
		if (c == -1)
			break;
		
		switch (c) {
		case 'c':
			doCreate = true;
			inputFileName = optarg;
			break;
		case 'd':
			if (optarg == NULL) {
				fprintf(stderr, "%s:no track-id specified for delete\n", ProgName);
				exit(EXIT_COMMAND_LINE);
			}
			if (sscanf(optarg, "%u", &deleteTrackId) != 1) {
				fprintf(stderr, 
					"%s: bad track-id specified: %s\n",
					ProgName, optarg);
				exit(EXIT_COMMAND_LINE);
			}
			doDelete = true;
			break;
		case 'e':
			if (optarg == NULL) {
				fprintf(stderr, "%s:no track-id specified for extract\n", ProgName);
				exit(EXIT_COMMAND_LINE);
			}
			
			if (sscanf(optarg, "%u", &extractTrackId) != 1) {
				fprintf(stderr, 
					"%s: bad track-id specified: %s\n",
					ProgName, optarg);
				exit(EXIT_COMMAND_LINE);
			}
			doExtract = true;
			break;
		case 'E':
			doEncrypt = true;
			if (optarg) {
				// if the short version of option is given, optarg has 
				// an = at the beginning. this causes sscanf to fail. 
				// if the long version of the option is given, there 
				// is no =
				if ( optarg[0] == '=' ) optarg[0] = ' ';
				if (sscanf(optarg, "%d", &encryptTrackId) != 1) {
					fprintf(stderr, 
						"%s: bad track-id specified: %s\n",
						ProgName, optarg);
					exit(EXIT_COMMAND_LINE);
				}
			}	
			break;
	//---------------------------------------------------------
		case 'f':
			doForcetimestamp = true;
			break;
	//---------------------------------------------------------//Lirh20071113
	//---------------------------------------------------------
		case 'g':
			doGettimestamp = true;
			break;
	//---------------------------------------------------------//Lirh20071113
		case 'H':
			doHint = true;
			if (optarg) {
				// if the short version of option is given, optarg has 
				// an = at the beginning. this causes sscanf to fail. 
				// if the long version of the option is given, there 
				// is no =
				if ( optarg[0] == '=' ) optarg[0] = ' ';
				if (sscanf(optarg, "%u", &hintTrackId) != 1) {
					fprintf(stderr, 
						"%s: bad track-id specified: %s\n",
						ProgName, optarg);
					exit(EXIT_COMMAND_LINE);
				}
			}
			break;
		case 'l':
			doList = true;
			break;
		case 'm':
			u_int32_t mtu;
			if (optarg == NULL) {
				fprintf(stderr, "%s:no mtu specified\n", ProgName);
				exit(EXIT_COMMAND_LINE);
			}
			if (sscanf(optarg, "%u", &mtu) != 1 || mtu < 64) {
				fprintf(stderr, 
					"%s: bad mtu specified: %s\n",
					ProgName, optarg);
				exit(EXIT_COMMAND_LINE);
			}
			maxPayloadSize = mtu - 40;	// subtract IP, UDP, and RTP hdrs
			break;
		case 'O':
			doOptimize = true;
			break;
		case 'r':
			if (optarg == NULL) {
				fprintf(stderr, "%s:no rate specifed\n", ProgName);
				exit(EXIT_COMMAND_LINE);
			}
			if (sscanf(optarg, "%lf", &VideoFrameRate) != 1) {
				fprintf(stderr, 
					"%s: bad rate specified: %s\n",
					ProgName, optarg);
				exit(EXIT_COMMAND_LINE);
			}
			break;
		case 't':
			if (optarg == NULL) {
				fprintf(stderr, "%s:no time scale specifed\n", ProgName);
				exit(EXIT_COMMAND_LINE);
			}
			if (sscanf(optarg, "%u", &Mp4TimeScale) != 1) {
				fprintf(stderr, 
					"%s: bad timescale specified: %s\n",
					ProgName, optarg);
				exit(EXIT_COMMAND_LINE);
			}
			TimeScaleSpecified = true;
			break;
		case 'u':
			createFlags |= MP4_CREATE_64BIT_DATA;
			break;
		case 'U':
			createFlags |= MP4_CREATE_64BIT_TIME;
			break;  
		case 'v':
			Verbosity |= (MP4_DETAILS_READ | MP4_DETAILS_WRITE);
			if (optarg) {
				u_int32_t level;
				if (sscanf(optarg, "%u", &level) == 1) {
					if (level >= 2) {
						Verbosity |= MP4_DETAILS_TABLE;
					} 
					if (level >= 3) {
						Verbosity |= MP4_DETAILS_SAMPLE;
					} 
					if (level >= 4) {
						Verbosity |= MP4_DETAILS_HINT;
					} 
					if (level >= 5) {
						Verbosity = MP4_DETAILS_ALL;
					} 
				}
			}
			break;
		case '?':
			fprintf(stderr, "usage: %s %s", ProgName, usageString);
			exit(EXIT_SUCCESS);
		case 'V':
			fprintf(stderr, "%s - %s version %s\n", 
				ProgName, MPEG4IP_PACKAGE, MPEG4IP_VERSION);
			exit(EXIT_SUCCESS);
		default:
			fprintf(stderr, "%s: unknown option specified, ignoring: %c\n", 
				ProgName, c);
    }
  }
  
  // check that we have at least one non-option argument
  if ((argc - optind) < 1) {
	  fprintf(stderr, "usage: %s %s", ProgName, usageString);                         
	  exit(EXIT_COMMAND_LINE);
  }
  
  if ((argc - optind) == 1) {														
	  mp4FileName = argv[optind++];		
  } else {
	  // it appears we have two file names											
	  if (doExtract) {
		  mp4FileName = argv[optind++];
		  outputFileName = argv[optind++];												
	  } else {																		
		  if (inputFileName == NULL) {
			  // then assume -c for the first file name
			  doCreate = true;
			  inputFileName = argv[optind++];													
		  }
		  mp4FileName = argv[optind++];
	  }
  }																					
  
  // warn about extraneous non-option arguments
  if (optind < argc) {
	  fprintf(stderr, "%s: unknown options specified, ignoring: ", ProgName);
	  while (optind < argc) {
		  fprintf(stderr, "%s ", argv[optind++]);
	  }
	  fprintf(stderr, "\n");
  }
  
  // operations consistency checks
  
  if (!doList && !doCreate && !doHint && !doEncrypt  
      && !doOptimize && !doExtract && !doDelete && !doIsma && !doGettimestamp) {		//Lirh20080528
	  fprintf(stderr, 
		  "%s: no operation specified\n",
		  ProgName);
	  exit(EXIT_COMMAND_LINE);
  }
  if ((doCreate || doHint || doEncrypt || doIsma) && (doExtract || doGettimestamp)) {	//Lirh20080528
	  fprintf(stderr, 
		  "%s: extract operation must be done separately\n",
		  ProgName);
	  exit(EXIT_COMMAND_LINE);
  }
  if ((doCreate || doHint || doEncrypt || doIsma) && doDelete) {
	  fprintf(stderr, 
		  "%s: delete operation must be done separately\n",
		  ProgName);
	  exit(EXIT_COMMAND_LINE);
  }
  if ((doExtract || doGettimestamp) && doDelete) {										//Lirh20080528
	  fprintf(stderr, 
		  "%s: extract and delete operations must be done separately\n",
		  ProgName);
	  exit(EXIT_COMMAND_LINE);
  }
  
  // end processing of command line
  
  if (doList) {
	  // just want the track listing
	  char* info = MP4FileInfo(mp4FileName);
	  
	  if (!info) {
		  fprintf(stderr, 
			  "%s: can't open %s\n", 
			  ProgName, mp4FileName);
		  exit(EXIT_INFO);
	  }
	  
	  fputs(info, stdout);
	  free(info);
	  exit(EXIT_SUCCESS);
  } 
  
  // test if mp4 file exists
  bool mp4FileExists = (access(mp4FileName, F_OK) == 0);

  //---------------------------------------------------------
  if(doForcetimestamp)
  {
	  	char datFileName[MAX_PATH];
		strcpy(datFileName, inputFileName);
		char* extension1 = strrchr(datFileName, '.');
		strcpy(extension1, ".dat");
		if(GetFileAttributes( datFileName ) == -1)
		{
				  fprintf(stderr,
					  "%s: Cannot find the dat file: %s\n", ProgName, datFileName);
				  exit(EXIT_COMMAND_LINE);
		}

  }
  //---------------------------------------------------------//Lirh20071113
  //---------------------------------------------------------
  if(doGettimestamp)
  {
	  ExtractTimestamp(mp4FileName);
  }
  //---------------------------------------------------------//Lirh20080528
 
  MP4FileHandle mp4File;
  
  if (doCreate || doHint) {
	  if (!mp4FileExists) {																				
		  if (doCreate) {
			  const char* extension = strrchr(inputFileName, '.');
			  
			  if (extension == NULL) {
				  fprintf(stderr,
					  "%s: unknown file type: %s\n", ProgName, inputFileName);
				  exit(EXIT_COMMAND_LINE);
			  }			  
		      mp4File = MP4Create(mp4FileName, Verbosity, createFlags);
			  if (mp4File) {
				  MP4SetTimeScale(mp4File, Mp4TimeScale);										
			  }
		  } else if (doEncrypt) {
			  fprintf(stderr,
				  "%s: can't encrypt track in file that doesn't exist\n", 
				  ProgName);
			  exit(EXIT_CREATE_FILE);
		  } else {
			  fprintf(stderr,
				  "%s: can't hint track in file that doesn't exist\n", 
				  ProgName);
			  exit(EXIT_CREATE_FILE);
		  }
	  } else {
		  if (createFlags != 0) {
			  fprintf(stderr, "Must specify 64 bits on new file only");
			  exit(EXIT_CREATE_FILE);
		  }
		  mp4File = MP4Modify(mp4FileName, Verbosity);					
	  }
	  
	  if (!mp4File) {
		  // mp4 library should have printed a message
		  exit(EXIT_CREATE_FILE);
	  }
	  
	  bool allMpeg4Streams = true;									
	  
	  if (doCreate) {																
		  MP4TrackId* pCreatedTrackIds = 
			  CreateMediaTracks(mp4File, inputFileName,									
			  doEncrypt, doForcetimestamp);//Lirh20071113
		  
		  if (pCreatedTrackIds == NULL) {
			  MP4Close(mp4File);
			  exit(EXIT_CREATE_MEDIA);
		  }
		  
		  MP4TrackId* pTrackId = pCreatedTrackIds;
		  
		  
		  if (doHint) {
			  MP4Close(mp4File);		 
			  mp4File = MP4Modify(mp4FileName, Verbosity);					
			  pTrackId = pCreatedTrackIds;
			  
			  while (*pTrackId != MP4_INVALID_TRACK_ID) {							
				  CreateHintTrack(mp4File, *pTrackId, payloadName, 
					  0, maxPayloadSize, doEncrypt);					
				  pTrackId++;
			  }
		  }
	  } else if (doHint) {
		  // in this case, we only hint the track specified in the command line"-H/-hint"
		  
		  CreateHintTrack(mp4File, hintTrackId, payloadName, 				
			  0, maxPayloadSize, doEncrypt);
		  
		  uint32_t trackNum;
		  
		  trackNum = MP4GetNumberOfTracks(mp4File);
		  for (uint32_t ix = 0; ix < trackNum; ix++) {
			  MP4TrackId trackId = MP4FindTrackId(mp4File, ix);
			  
			  const char *type =
				  MP4GetTrackType(mp4File, trackId);
			  if (MP4HaveTrackIntegerProperty(mp4File, trackId,
				  "mdia.minf.stbl.stsd.*.esds.decConfigDescr.objectTypeId")) {
				  if (!strcmp(type, MP4_AUDIO_TRACK_TYPE)) { 
					  allMpeg4Streams &=
						  (MP4GetTrackEsdsObjectTypeId(mp4File, trackId) 
						  == MP4_MPEG4_AUDIO_TYPE);
					  
				  } else if (!strcmp(type, MP4_VIDEO_TRACK_TYPE)) { 
					  allMpeg4Streams &=
						  (MP4GetTrackEsdsObjectTypeId(mp4File, trackId) 
						  == MP4_MPEG4_VIDEO_TYPE);
				  }
			  }
		  }
	  }
	  
	  char *buffer;
	  char *value;
	  uint32_t newverbosity;
	  newverbosity = Verbosity & ~(MP4_DETAILS_ERROR);
	  MP4SetVerbosity(mp4File, newverbosity);
	  bool retval = MP4GetMetadataTool(mp4File, &value);
	  MP4SetVerbosity(mp4File, Verbosity);
	  if (retval && value != NULL) {
		  if (strncasecmp("mp4creator", value, strlen("mp4creator")) != 0) {
			  buffer = (char *)malloc(strlen(value) + 80);
			  sprintf(buffer, "%s mp4creator %s", value, MPEG4IP_VERSION);
			  MP4SetMetadataTool(mp4File, buffer);
			  free(buffer);
		  }
	  } else {
		  buffer = (char *)malloc(80);
		  sprintf(buffer, "mp4creator %s", MPEG4IP_VERSION);
		  MP4SetMetadataTool(mp4File, buffer);
	  }
	  MP4Close(mp4File);					
  } 
  else if (doEncrypt)
  { 
	  //加密选项接口
  } 
  else if (doExtract) {
	  if (!mp4FileExists) {
		  fprintf(stderr,
			  "%s: can't extract track in file that doesn't exist\n", 
			  ProgName);
		  exit(EXIT_CREATE_FILE);
	  }
	  
	  mp4File = MP4Read(mp4FileName, Verbosity);
	  if (!mp4File) {
		  // mp4 library should have printed a message
		  exit(EXIT_CREATE_FILE);
	  }
	  
	  char tempName[PATH_MAX];
	  if (outputFileName == NULL) {
		  snprintf(tempName, sizeof(tempName), 
			  "%s.t%u", mp4FileName, extractTrackId);
		  outputFileName = tempName;
	  }
	  
	  ExtractTrack(mp4File, extractTrackId, outputFileName);
	  
	  MP4Close(mp4File);
	  
  } else if (doDelete) {
	  if (!mp4FileExists) {
		  fprintf(stderr,
			  "%s: can't delete track in file that doesn't exist\n", 
			  ProgName);
		  exit(EXIT_CREATE_FILE);
	  }
	  
	  mp4File = MP4Modify(mp4FileName, Verbosity);
	  if (!mp4File) {
		  // mp4 library should have printed a message
		  exit(EXIT_CREATE_FILE);
	  }
	  
	  MP4DeleteTrack(mp4File, deleteTrackId);
	  
	  MP4Close(mp4File);
	  
	  doOptimize = true;	// to purge unreferenced track data
  }
  
  if (doIsma) {
	  MP4MakeIsmaCompliant(mp4FileName, Verbosity);
  }
  
  if (doOptimize) {
	  if (!MP4Optimize(mp4FileName, NULL, Verbosity)) {
		  // mp4 library should have printed a message
		  exit(EXIT_OPTIMIZE_FILE);
	  }
  }
  
  return(EXIT_SUCCESS);
}


MP4TrackId* CreateMediaTracks(MP4FileHandle mp4File, const char* inputFileName,
							  bool doEncrypt,bool doForcetimestamp)//Lirh20071113
{
	FILE* inFile = fopen(inputFileName, "rb");

//---------------------------------------------------------
	FILE* datFile = NULL;

	if(doForcetimestamp)
	{
		char datFileName[MAX_PATH];
		strcpy(datFileName, inputFileName);
		char* extension1 = strrchr(datFileName, '.');
		strcpy(extension1, ".dat");
		datFile = fopen(datFileName, "rb");	
	}
//---------------------------------------------------------//Lirh20071113

	//no file
	if (inFile == NULL) {
		fprintf(stderr, 
			"%s: can't open file %s: %s\n",
			ProgName, inputFileName, strerror(errno));
		return NULL;
	}
	//error expression
	struct stat s;
	if (fstat(fileno(inFile), &s) < 0) {
		fprintf(stderr, 
			"%s: can't stat file %s: %s\n",
			ProgName, inputFileName, strerror(errno));
		return NULL;
	}
	
	if (s.st_size == 0) {
		fprintf(stderr, 
			"%s: file %s is empty\n",
			ProgName, inputFileName);
		return NULL;
	}
	
	const char* extension = strrchr(inputFileName, '.');					
	if (extension == NULL) {
		fprintf(stderr, 
			"%s: no file type extension\n", ProgName);
		return NULL;
	}
	
	static MP4TrackId trackIds[2] = {
		MP4_INVALID_TRACK_ID, MP4_INVALID_TRACK_ID							
	};
	MP4TrackId* pTrackIds = trackIds;
	
	if (!strcasecmp(extension, ".mp3")) {
		trackIds[0] = Mp3Creator(mp4File, inFile, doEncrypt);
		
	} else if (strcasecmp(extension, ".avsm")== 0) {
		trackIds[0] = AVSMCreator(mp4File, inFile);										
	} else if (strcasecmp(extension, ".avs")== 0) {
		trackIds[0] = AVSCreator(mp4File, inFile, datFile);		//Lirh20071113								
	} else {
		fprintf(stderr, 
			"%s: unknown file type\n", ProgName);
		return NULL;
	}
	
	if (inFile) {
		fclose(inFile);
	}

//---------------------------------------------------------
	if (datFile) {
		fclose(datFile);
	}
//---------------------------------------------------------//Lirh20071113
	
	if (pTrackIds == NULL || pTrackIds[0] == MP4_INVALID_TRACK_ID) {
		return NULL;
	}
	
	return pTrackIds;
}


void CreateHintTrack(MP4FileHandle mp4File, MP4TrackId mediaTrackId,	
					 const char* payloadName, bool interleave,					
					 u_int16_t maxPayloadSize, bool doEncrypt)					
{
	
	bool rc = FALSE;
	
	if (MP4GetTrackNumberOfSamples(mp4File, mediaTrackId) == 0) {
		fprintf(stderr, 
			"%s: couldn't create hint track, no media samples\n", ProgName);
		MP4Close(mp4File);
		exit(EXIT_CREATE_HINT);
	}
	
	// vector out to specific hinters
	const char* trackType = MP4GetTrackType(mp4File, mediaTrackId);
	
	if (0)//if( 加密 || 媒体流是通过采用某种加密算法得到)
	{				
		
	
	}
   else if (!strcmp(trackType, MP4_AUDIO_TRACK_TYPE)) {							//audio track no encrypted
	  const char *media_data_name;
	  media_data_name = MP4GetTrackMediaDataName(mp4File, mediaTrackId);			
	  
	  if (strcasecmp(media_data_name, "mp4a") == 0) {
		  u_int8_t audioType = MP4GetTrackEsdsObjectTypeId(mp4File, mediaTrackId);
		  
		  switch (audioType) {
		  case MP4_MPEG4_AUDIO_TYPE:
		  case MP4_MPEG2_AAC_MAIN_AUDIO_TYPE:
		  case MP4_MPEG2_AAC_LC_AUDIO_TYPE:
		  case MP4_MPEG2_AAC_SSR_AUDIO_TYPE:
			  rc = MP4AV_RfcIsmaHinter(mp4File, mediaTrackId, 
				  interleave, maxPayloadSize);
			  break;
		  case MP4_MPEG1_AUDIO_TYPE:
		  case MP4_MPEG2_AUDIO_TYPE:
			  if (payloadName && 
				  (!strcasecmp(payloadName, "3119") 
				  || !strcasecmp(payloadName, "mpa-robust"))) {
				  rc = MP4AV_Rfc3119Hinter(mp4File, mediaTrackId, 
					  interleave, maxPayloadSize);
			  } else {
				  rc = MP4AV_Rfc2250Hinter(mp4File, mediaTrackId, 
					  false, maxPayloadSize);
			  }
			  break;
		  case MP4_PCM16_BIG_ENDIAN_AUDIO_TYPE:
		  case MP4_PCM16_LITTLE_ENDIAN_AUDIO_TYPE:
			  rc = L16Hinter(mp4File, mediaTrackId, maxPayloadSize);
			  break;
		  default:
			  fprintf(stderr, 
				  "%s: can't hint non-MPEG4/non-MP3 audio type\n", ProgName);
		  }
	  } else if (strcasecmp(media_data_name, "samr") == 0 ||
		  strcasecmp(media_data_name, "sawb") == 0) {
		  rc = MP4AV_Rfc3267Hinter(mp4File, mediaTrackId, maxPayloadSize);
	  }
  }
  else if (!strcmp(trackType, MP4_VIDEO_TRACK_TYPE)) {							//video track no encrypted
	  const char *media_data_name;
	  media_data_name = MP4GetTrackMediaDataName(mp4File, mediaTrackId);			
	  
	  if (strcasecmp(media_data_name, "avs1") == 0) {
		  // AVSM;
		  rc = MP4AV_AVSMHinter(mp4File, mediaTrackId, maxPayloadSize);					//AVSM_Hinter
	  }else if (strcasecmp(media_data_name, "avs2") == 0) {
		  // AVS1-P2;
		  rc = MP4AV_AVSHinter(mp4File, mediaTrackId, maxPayloadSize);					//AVS1-P2_Hinter
	  }
  } else {
	  fprintf(stderr, 
		  "%s: can't hint track type %s\n", ProgName, trackType);
  }
  
  if (!rc) {
	  fprintf(stderr, 
		  "%s: error hinting track %u\n", ProgName, mediaTrackId);
	  MP4Close(mp4File);
	  exit(EXIT_CREATE_HINT);
  }
}




void ExtractTrack (MP4FileHandle mp4File, MP4TrackId trackId, 
				   const char* outputFileName)
{
	int openFlags = O_WRONLY | O_TRUNC | OPEN_CREAT;
	u_int8_t amrType = AMR_TYPE_NONE;
	int outFd = open(outputFileName, openFlags, 0644);
	
	bool media_data_is_avs2 = false ;
	if (outFd == -1) {
		fprintf(stderr, "%s: can't open %s: %s\n",
			ProgName, outputFileName, strerror(errno));
		exit(EXIT_EXTRACT_TRACK);
	}
	
	
	const char* trackType =
		MP4GetTrackType(mp4File, trackId);
	const char *media_data_name = 
		MP4GetTrackMediaDataName(mp4File, trackId);
	
	if (!strcmp(trackType, MP4_VIDEO_TRACK_TYPE)) 
	{
        if (strcmp(media_data_name, "avs2") == 0) 
			media_data_is_avs2 = true ;
		else
		{
			fprintf(stderr, "%s:  media data name  %s  is unknown \n",
				ProgName, media_data_name);
			exit(EXIT_EXTRACT_TRACK);
		}
	} 
	else if (!strcmp(trackType, MP4_AUDIO_TRACK_TYPE))
	{
		if (strcmp(media_data_name, "mp4a") != 0)
		{
			fprintf(stderr, "%s:  media data name  %s  is unknown \n",
				ProgName, media_data_name);
			exit(EXIT_EXTRACT_TRACK);
		}
	}
	else
	{
        fprintf(stderr, "%s:  track type is unknown \n",ProgName);
		exit(EXIT_EXTRACT_TRACK);
	}
	
	MP4SampleId numSamples = 
		MP4GetTrackNumberOfSamples(mp4File, trackId);
	u_int8_t* pSample;
	u_int32_t sampleSize;

	MP4Timestamp pStartTime;

	// extraction loop
	for (MP4SampleId sampleId = 1 ; sampleId <= numSamples; sampleId++) 
	{
		int rc;
		pSample = NULL;
		sampleSize = 0;
		
		rc = MP4ReadSample(
			mp4File, 
			trackId, 
			sampleId, 
			&pSample, 
			&sampleSize,
			&pStartTime);

		if (rc == 0)
		{
				fprintf(stderr, "%s: read sample %u for %s failed\n",
					ProgName, sampleId, outputFileName);
				exit(EXIT_EXTRACT_TRACK);
		}

		if( media_data_is_avs2 )
		{
	        u_int8_t* dest_pSample=NULL;
        	u_int32_t dest_sampleSize;
	        dest_pSample = (uint8_t *)realloc(dest_pSample,sampleSize); 
//			printf("sampleId:%d  sampleSize:%u\n",sampleId,sampleSize);
            rc = ExtractAvsSlice( pSample , sampleSize, dest_pSample , &dest_sampleSize ) ;

			if( rc ) 
			{
//				printf("sampleId:%d  sampleSize:%u\n",sampleId,sampleSize);
				int t=dest_sampleSize;

				rc = write(outFd, dest_pSample, dest_sampleSize);
			
	    		if (rc == -1 || (u_int32_t)rc != dest_sampleSize)
				{
		    		fprintf(stderr, "%s: write to %s failed: %s\n",
	    				ProgName, outputFileName, strerror(errno));
		    		exit(EXIT_EXTRACT_TRACK);
				}
			}
            free( dest_pSample ) ;
		}
		else
		{
			rc = write(outFd, pSample, sampleSize);		
//			printf("sampleId:%d  pStartTime:%u\n",sampleId,pStartTime);
	    	if (rc == -1 || (u_int32_t)rc != sampleSize)
			{
				fprintf(stderr, "%s: write to %s failed: %s\n",
				ProgName, outputFileName, strerror(errno));
			    exit(EXIT_EXTRACT_TRACK);
			}
		}
		
		free(pSample);
	}
	
	// close ES file
	close(outFd);
}


bool  ExtractAvsSlice( u_int8_t *source , u_int32_t sourceSize , u_int8_t*dest , u_int32_t*destSize ) 
{
	u_int32_t len = 0 ;
	int index = 0 ;
	*destSize = 0 ;
	u_int8_t *pSample = source ;
	u_int32_t *NalLen;
	while( len < sourceSize )
	{
        pSample = &source[len] ;
		load_avs_start_code( dest , index ) ;
		index += 3 ;
		adjust_Len_Pos( pSample ) ;
        NalLen=(u_int32_t*)pSample ;
		int nallen = NalLen[0] ;
		if( nallen + len + 4 > sourceSize || nallen < 1)
		{ 
			return false;
		}
		memcpy( &dest[index] , &pSample[5] , nallen-1 ) ;
		*destSize += (nallen+2) ;
		index = *destSize ;
		len += (nallen+4) ;			
	}
	return true;
}

//插入起始码前缀
void load_avs_start_code( u_int8_t *dest , u_int32_t index )
{
	dest[index++] = 0x00 ;
	dest[index++] = 0x00 ;
	dest[index++] = 0x01 ;
}


void adjust_Len_Pos( u_int8_t *pSample ) 
{
	u_int8_t temp[4] ;
	for( int i = 0 ; i < 4 ; i++ )
		temp[3-i] = pSample[i] ;
	for( int j = 0 ; j < 4 ; j++ )
        pSample[j] = temp[j] ;
}


 

//---------------------------------------------------------//Lirh20080528
void ExtractTimestamp( const char* mp4FileName )
{
	uint32_t count,trackNum;
	bool media_data_is_avs2 = false ;
//	bool headerpart = true ;
	u_int8_t* pSample;
	u_int32_t sampleSize;
    u_int32_t CountTime = 0 ;
	char datFileName[MAX_PATH];
	FILE* datFile = NULL;
	
	MP4FileHandle mp4File; 
	MP4TrackId trackId;
	MP4SampleId numSamples;
	MP4Timestamp pStartTime;
	MP4Duration pDuration;
	MP4Timestamp iTime;

	mp4File = MP4Modify(mp4FileName, Verbosity);
	trackNum = MP4GetNumberOfTracks(mp4File);
	
	strcpy(datFileName, mp4FileName);
	char* extension1 = strrchr(datFileName, '.');
	strcpy(extension1, ".ex.dat");
	datFile = fopen(datFileName, "wb");	

	for(count=0; count<trackNum ;count++)
	{
		trackId = MP4FindTrackId(mp4File, count);
		const char* trackType =
			MP4GetTrackType(mp4File, trackId);
		const char *media_data_name = 
			MP4GetTrackMediaDataName(mp4File, trackId);
		
		if (!strcmp(trackType, MP4_VIDEO_TRACK_TYPE)) 
		{
			if (strcmp(media_data_name, "avs2") == 0) 
			{
				media_data_is_avs2 = true ;
				break;
			}
			else
			{
				fprintf(stderr, "%s:  media data name  %s  is unknown \n",
					ProgName, media_data_name);
				exit(EXIT_EXTRACT_TRACK);
			}
		} 
	}
	
	if(media_data_is_avs2 == false)
	{
		fprintf(stderr, "Can not find avs track in %s \n", mp4FileName );
		MP4Close(mp4File);
		exit(EXIT_EXTRACT_TRACK);
	}

	numSamples = MP4GetTrackNumberOfSamples(mp4File, trackId);

	// extraction loop
	for (MP4SampleId sampleId = 1 ; sampleId <= numSamples; sampleId++) 
	{
		int rc;
		pSample = NULL;
		sampleSize = 0;
		
		rc = MP4ReadSample(
			mp4File, 
			trackId, 
			sampleId, 
			&pSample, 
			&sampleSize,
			&pStartTime,
			&pDuration);

		if (rc == 0)
		{
				fprintf(stderr, "%s: read sample %u failed\n",
					ProgName, sampleId);
				exit(EXIT_EXTRACT_TRACK);
		}

//		if(pDuration > 0)
//			headerpart = false;

		if(pDuration > 0 || sampleId == 1)
		{
			CountTime += pDuration;
			iTime = CountTime * (double)1000 / Mp4TimeScale ;
//			printf("%d \n", pDuration);
			fwrite(&iTime, 4, 1, datFile);
		}

		free(pSample);
	}

	if (datFile) 
		fclose(datFile);
	MP4Close(mp4File);
}
//---------------------------------------------------------//Lirh20080528
