/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2004,  Audio Video Coding Standard, Part III
*
* This software module was originally developed by
*
* JungHoe Kim (kjh94@samsung.com), Samsung AIT
*
* edited by
*
* Lei Miao (win.miaolei@samsung.com), Samsung AIT
* Lei Miao, CBC Multi-channel extension, 2005-09-19
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

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>


static FILE *ptrBitstream;
static int buffer;		/* Bits buffered for output */
static int bits_to_go;
static unsigned char OutputBuffer[8191*50];
static int CurrentByte;


void output_byte(long byte,int len);
void start_outputing_bits();
void done_outputing_bits();
void FlushBuffer(void);
void FlushBufferWithLength(void);
int BitstreamOpen(FILE *fname);
void BitstreamClose(void);
int GetBitstreamSize(void);
int ByteAlign(void);


unsigned char OutputBuffer_ext[8191*8];
int nOutput_ext;


void FlushBufferInit();
void FlushBufferWithLength_ext(void);
void FlushFrame(int fill_bits);

// 2006-01-20 xusen
void ScanAATFFrame(void);

//for AASF
void FlushAASFHeaderBuffer();

//for AATF        
void Register_buffer(void);
int	EncodeBSHCGetAATFSize();
void Restore_buffer(void);
void FlushAATFHeaderBuffer();


static int buffer_last;
static int bits_to_go_last;
static int CurrentByte_last;

int BitstreamOpen(FILE *fname)
{
	ptrBitstream = fname;
	nOutput_ext = 0;
	if(ptrBitstream==NULL) return 1;
	return 0;
}
void BitstreamClose(void)
{
	fclose(ptrBitstream);
}
int GetBitstreamSize(void)
{
	return CurrentByte*8+(8-bits_to_go);
}
void FlushBuffer(void)
{
	fwrite(OutputBuffer,1,CurrentByte,ptrBitstream);
	CurrentByte = 0;
}
void FlushBufferWithLength(void)
{
	fwrite(OutputBuffer,1,CurrentByte,ptrBitstream);
}

void FlushBufferInit()
{
	nOutput_ext = 0;
}

void FlushBufferWithLength_ext(void)
{
	memcpy(OutputBuffer_ext+nOutput_ext,OutputBuffer,CurrentByte);
	nOutput_ext += CurrentByte;
}

void FlushFrame(int fill_bits)
{
	fwrite(OutputBuffer,1,CurrentByte-fill_bits/8,ptrBitstream);
	fwrite(OutputBuffer_ext,1,nOutput_ext,ptrBitstream);
	fwrite(OutputBuffer+CurrentByte-fill_bits/8, 1, fill_bits/8, ptrBitstream);
}

void start_outputing_bits()
{   
    buffer = 0;					/* Buffer is empty to start with*/
    bits_to_go= 8;				                   
	CurrentByte = 0;
}
int PutByte(unsigned char c)
{
	if(CurrentByte>8191)
	{
		fprintf(stderr,"\n\n\t\t\terr");
		return 1;
	}
	OutputBuffer[CurrentByte] = c;
	CurrentByte++;	
	return 0;
}
void output_bit(int bit)
{
	buffer <<= 1; if (bit) buffer |= 0x1;	/* Put bit in top of buffer.*/
    bits_to_go -= 1;
    if (bits_to_go==0) {			/* Output buffer if it is   */
		PutByte(buffer);			/* now full.                */
        bits_to_go = 8;
    }
}

void done_outputing_bits()
{   
	PutByte(buffer<<bits_to_go);
	buffer = 0;					/* Buffer is empty to start */
    bits_to_go= 8;				/* with.                    */
}
void output_byte(long byte,int len)
{
	int i;
	int mask;
	/* MSB first */
	mask = 1<<(len-1);
	for(i=0;i<len;i++)
	{
		if(byte & mask)
			output_bit(1);
		else
			output_bit(0);		
		mask >>= 1;
	}
}
int ByteAlign(void)
{
	if(bits_to_go!=8)
	{
		output_byte(0,bits_to_go);  /*byte align*/
		return bits_to_go;
	}
	else
		return 0;
}

/************************************************************************/
/*                              AASF                                    */
/************************************************************************/

void FlushAASFHeaderBuffer()
{
	OutputBuffer[4] = (CurrentByte>>16)&0xFF;
	OutputBuffer[5] = (CurrentByte>>8)&0xFF;
	OutputBuffer[6] = CurrentByte&0xFF;
	
	fwrite(OutputBuffer,1,CurrentByte,ptrBitstream);
}

void Register_buffer(void)
{
	buffer_last = buffer;					
	bits_to_go_last = bits_to_go;				                   
	CurrentByte_last = CurrentByte;
}

int	EncodeBSHCGetAATFSize()
{
	return (CurrentByte - CurrentByte_last)<<3;
}

void Restore_buffer(void)
{
	buffer = buffer_last;					
	bits_to_go = bits_to_go_last;				                   
	CurrentByte = CurrentByte_last;

}

void FlushAATFHeaderBuffer()
{
	// 2006-01-19 xusen deleted frame length in aatf header
//	OutputBuffer[3] |= (CurrentByte & 0x00001FE0) >> 5;
//	OutputBuffer[4] |= (CurrentByte & 0x0000001F) << 3;
	ScanAATFFrame();
	fwrite(OutputBuffer,1,CurrentByte,ptrBitstream);
}

void FlushAATFMCFrame()
{
	int frame_len;
	frame_len = CurrentByte+nOutput_ext;

	// 2006-01-19 xusen, not debugging in multichannel case yet
//	OutputBuffer[3] |= (CurrentByte & 0x00001FE0) >> 5;
//	OutputBuffer[4] |= (CurrentByte & 0x0000001F) << 3;
	memcpy(&OutputBuffer[CurrentByte], OutputBuffer_ext, sizeof(unsigned char) * nOutput_ext);
	ScanAATFFrame();

	fwrite(OutputBuffer,1,CurrentByte,ptrBitstream);
}

unsigned char mask1[8] = {0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
unsigned char mask0[8] = {0x7f, 0xbf, 0xdf, 0xef, 0xf7, 0xfb, 0xfd, 0xfe};

void ScanAATFFrame(void)
{
	int NewCurrentByte;					// Byte number after scan and interpolating.
	int TotalBits;
	int NewBits;						// number of the bit being scanned
	int AlignBits;
	unsigned char NewOutputBuffer[2048];
	int i;
	int zeroCounter = 0;
	int intCounter = 0;

	/* syncword which is not scanned */
	NewOutputBuffer[0] = 0x00;
	NewOutputBuffer[1] = 0x10;

	TotalBits = (CurrentByte<<3);	// number of bits need to scan.
	NewBits = 12;

	/* scan */
	for (i = 12; i < TotalBits; i ++) {
		if (NewBits%8 == 3 && zeroCounter >= 11) {
			NewOutputBuffer[NewBits/8] &= 0xEF;
			NewBits ++;
			intCounter++;
			zeroCounter++;
		}

		if (OutputBuffer[i/8]&mask1[i%8]) {
			NewOutputBuffer[NewBits/8] |= mask1[NewBits%8];
			NewBits++;
			zeroCounter = 0;
		}else{
			NewOutputBuffer[NewBits/8] &= mask0[NewBits%8];
			NewBits++;
			zeroCounter++;
		}
	}
	
	NewCurrentByte = (NewBits+7)/8;
	/* align */
	AlignBits = NewBits%8 == 0 ? 0 : 8-NewBits%8;
    if (AlignBits == 5 && zeroCounter >= 11) {
		NewOutputBuffer[NewCurrentByte-1] &= 0xe0;
		NewOutputBuffer[NewCurrentByte-1] |= 0x0f;
    }else if(AlignBits != 0)
		NewOutputBuffer[NewCurrentByte-1] |= (1<<AlignBits)-1;
	
	CurrentByte = NewCurrentByte;
	memcpy(OutputBuffer, NewOutputBuffer, sizeof(unsigned char)*CurrentByte);

	return;
}

