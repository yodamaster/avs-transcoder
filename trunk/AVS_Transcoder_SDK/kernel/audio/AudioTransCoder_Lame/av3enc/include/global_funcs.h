
extern int BitstreamOpen(char *fname);
extern void BitstreamClose(void);
extern void output_byte(long byte,int len);
extern void start_outputing_bits();
extern void done_outputing_bits();
extern void FlushBuffer(void);
extern void FlushBufferWithLength(void);
extern int GetBitstreamSize(void);
extern int ByteAlign(void);
//~bs
extern void BSHCInit(FILE *name);
extern void BSHCClose(void);
extern int EncodeBSHCQuant(int model,int upper,int cur);
extern int EncodeBSHCSi(int siVal);
extern int EncodeBSHCSf(int sfVal);
extern int EncodeBSHCBin(int ms);
extern void EncodeBSHCStart(void);
extern void EncodeBSHCEnd(void);
extern void EncodeBSHCPutBits(int val,int len);
extern int EncodeBSHCGetSize();
extern void EncodeBSHCHeader(int ch,int freq);
extern int EncodeBSHCByteAlign(void);
extern int BSHCModelSelect(int bshcModel,int bpl);
extern void EncodeBSHCFlush(void);
extern int WriteAASFHeader(int ch,int freq,int bitrate,int aasf_flag);		//for AASF




extern void output_byte(long byte,int len);
extern void start_outputing_bits();
extern void done_outputing_bits();
extern void FlushBuffer(void);
extern void FlushBufferWithLength(void);
extern int BitstreamOpen(FILE *fname);
extern void BitstreamClose(void);
extern int GetBitstreamSize(void);
extern int ByteAlign(void);


extern void FlushBufferInit();
extern void FlushBufferWithLength_ext(void);
extern void FlushFrame(int fill_bits);

extern void FlushAATFMCFrame()
extern void FlushBufferInit()
// 2006-01-20 xusen
extern void ScanAATFFrame(void);

//for AASF
extern void FlushAASFHeaderBuffer();

//for AATF        
extern void Register_buffer(void);
extern int	EncodeBSHCGetAATFSize();
extern void Restore_buffer(void);
extern void FlushAATFHeaderBuffer();