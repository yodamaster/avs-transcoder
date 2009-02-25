;*****************************************************************************
;* dct_sse.asm: avs encoder library
;*****************************************************************************

BITS 32

;=============================================================================
; Macros and other preprocessor constants
;=============================================================================
%macro cglobal 1
    %ifdef PREFIX
        global _%1
        %define %1 _%1
    %else
        global %1
    %endif
%endmacro

%macro SSE_SUM_SUBW 2
    paddw   %1, %2
    paddw   %2, %2
    psubw   %2, %1
%endmacro

%macro SSE_SUM_SUBW_ROUND 4
    movdqu  %4, %1
    paddsw  %2, %3
    paddsw  %1, %2
    psubsw  %2, %4
%endmacro

%macro SSE_LOAD_SUM_SUBW 4
    movdqu            %2, %3
    movdqu            %1, %4
    SSE_SUM_SUBW   %1, %2
%endmacro

%macro SSE_SUM_SUBD 2
    paddd   %1, %2
    paddd   %2, %2
    psubd   %2, %1
%endmacro

%macro SSE_LOAD_WORD2XMM 3
    movq              %3, %2
    movq2dq           %1, %3
    pshufd            %1, %1, 114
    pshufhw           %1, %1, 216
    pshuflw           %1, %1, 216
    psrad             %1, 16
%endmacro

%macro SSE_LOAD_SUM_SUBD 5
    SSE_LOAD_WORD2XMM %2, %3, %5
    SSE_LOAD_WORD2XMM %1, %4, %5
    SSE_SUM_SUBD      %1, %2
%endmacro


%macro SSE_STORE_XMM2MMX2MM 3
    pshufhw           %3, %3, 216
    pshuflw           %3, %3, 216
    pshufd            %3, %3, 216
    movdq2q           %2, %3
    movq              %1, %2
%endmacro

%macro SBUTTERFLYwd 3
    movq        %3, %1
    punpcklwd   %1, %2
    punpckhwd   %3, %2
%endmacro

%macro SBUTTERFLYdq 3
    movq        %3, %1
    punpckldq   %1, %2
    punpckhdq   %3, %2
%endmacro

%macro XMM_INTO_MMX 3
    movdq2q        %3, %1
    pshufd         %1, %1, 78
    movdq2q        %2, %1
%endmacro

%macro MMX_INTO_XMM 4
    movq2dq        %1, %3
    movq2dq        %2, %4
    pshufd         %1, %1, 78
    paddd          %1, %2
%endmacro


;-----------------------------------------------------------------------------
; input ABCD output ADTC
;-----------------------------------------------------------------------------
%macro TRANSPOSE_SSE 5
    SBUTTERFLYwd %1, %2, %5
    SBUTTERFLYwd %3, %4, %2
    SBUTTERFLYdq %1, %3, %4
    SBUTTERFLYdq %5, %2, %3
%endmacro

;=============================================================================
; Code
;=============================================================================

SECTION .text

;-----------------------------------------------------------------------------
;   void __cdecl avs_1st_vdct8_sse( __int16 data[8][8] );
;-----------------------------------------------------------------------------
ALIGN 16
avs_1st_vdct8_sse:
    mov         eax, [esp+8]           ; load data[8][8] pointer

    ;-------------------------------------------------------------------------
    ; vertical dct ( compute 4 columns at a time -> 2 loops )
    ;-------------------------------------------------------------------------
    SSE_LOAD_SUM_SUBW      xmm0,      xmm4,      [eax],           [eax+112]
    SSE_LOAD_SUM_SUBW      xmm1,      xmm5,      [eax + 16],      [eax+ 96]
    SSE_LOAD_SUM_SUBW      xmm2,      xmm6,      [eax + 32],      [eax+ 80]
    SSE_LOAD_SUM_SUBW      xmm3,      xmm7,      [eax + 48],      [eax+ 64]

    SSE_SUM_SUBW           xmm3,      xmm0       ;tmp0 tmp2
    SSE_SUM_SUBW           xmm2,      xmm1       ;tmp1 tmp3

    SSE_SUM_SUBW           xmm2,      xmm3
    psllw                 xmm2,         3       ;b0
    psllw                 xmm3,         3       ;b1

    movdqu                [eax],     xmm2      ;data[0][0:7] = xmm2
    movdqu                [eax+64],  xmm3      ;data[4][0:7] = xmm3

    movdqu                xmm2,      xmm0       ;tmp2
    movdqu                xmm3,      xmm1       ;tmp3

    psllw                 xmm2,         1
    paddw                 xmm3,      xmm2
    psllw                 xmm3,         2
    paddw                 xmm3,      xmm2       ;b[2]

    psllw                 xmm1,         1
    psubw                 xmm0,      xmm1
    psllw                 xmm0,         2
    psubw                 xmm0,      xmm1       ;b[3]

    movdqu                [eax+32],  xmm3      ;data[2][0:7] = b[2] = xmm3
    movdqu                [eax+96],  xmm0      ;data[6][0:7] = b[3] = xmm0

    movdqu                xmm0,      xmm4
    movdqu                xmm1,      xmm5
    movdqu                xmm2,      xmm6
    movdqu                xmm3,      xmm7

    SSE_SUM_SUBW           xmm7,      xmm4      ;tmp4+tmp7, tmp4-tmp7
    SSE_SUM_SUBW           xmm6,      xmm5      ;tmp5+tmp6, tmp5-tmp6

    psllw                 xmm4,         1
    psllw                 xmm5,         1
    psllw                 xmm6,         1
    psllw                 xmm7,         1

    paddw                 xmm4,      xmm0      ;tmp0
    paddw                 xmm6,      xmm1      ;tmp1
    psubw                 xmm5,      xmm2      ;tmp2
    paddw                 xmm7,      xmm3      ;tmp3

    movdqu                xmm0,      xmm4      ;tmp0
    paddw                 xmm0,      xmm6      ;tmp0+tmp1
    paddw                 xmm0,      xmm7      ;tmp0+tmp1+tmp3
    psllw                 xmm0,         1      ;(tmp0+tmp1+tmp3)<<1
    paddw                 xmm0,      xmm6      ;b4=(tmp0+tmp1+tmp3)<<1+tmp1

    movdqu                xmm1,      xmm4      ;tmp0
    psubw                 xmm1,      xmm6      ;tmp0-tmp1
    paddw                 xmm1,      xmm5      ;tmp0-tmp1+tmp2
    psllw                 xmm1,         1      ;(tmp0-tmp1+tmp2)<<1
    paddw                 xmm1,      xmm4      ;b5=(tmp0-tmp1+tmp2)<<1+tmp0

    pxor                  xmm2,      xmm2      ;0
    psubw                 xmm2,      xmm6      ;-tmp1
    psubw                 xmm2,      xmm5      ;-tmp1-tmp2
    paddw                 xmm2,      xmm7      ;-tmp1-tmp2+tmp3
    psllw                 xmm2,         1      ;(-tmp1-tmp2+tmp3)<<1
    paddw                 xmm2,      xmm7      ;b6=(-tmp1-tmp2+tmp3)<<1+tmp3

    movdqu                xmm3,      xmm4      ;tmp0
    psubw                 xmm3,      xmm5      ;tmp0-tmp2
    psubw                 xmm3,      xmm7      ;tmp0-tmp2-tmp3
    psllw                 xmm3,         1      ;(tmp0-tmp2-tmp3)<<1
    psubw                 xmm3,      xmm5      ;b7=(tmp0-tmp2-tmp3)<<1-tmp2
    

    
    movdqu                [eax+16],  xmm0      ;data[1][0:7] = xmm0
    movdqu                [eax+48],  xmm1      ;data[3][0:7] = xmm1
    movdqu                [eax+80],  xmm2      ;data[5][0:7] = xmm2
    movdqu                [eax+112], xmm3      ;data[7][0:7] = xmm3
    
    ret

;-----------------------------------------------------------------------------
;   void __cdecl avs_2nd_vdct8_sse( __int16 data[8][8] );
;-----------------------------------------------------------------------------
ALIGN 16
avs_2nd_vdct8_sse:
    mov         eax, [esp+8]           ; load data[8][8] pointer

    ;-------------------------------------------------------------------------
    ; vertical dct ( compute 4 columns at a time -> 2 loops )
    ;-------------------------------------------------------------------------
    %assign disp 0
    %rep 2
    add                    eax,       disp
    
    SSE_LOAD_SUM_SUBD      xmm0,      xmm4,      [eax ],          [eax + 112],  mm0
    SSE_LOAD_SUM_SUBD      xmm1,      xmm5,      [eax + 16],      [eax + 96],  mm0
    SSE_LOAD_SUM_SUBD      xmm2,      xmm6,      [eax + 32],      [eax + 80],  mm0
    SSE_LOAD_SUM_SUBD      xmm3,      xmm7,      [eax + 48],      [eax + 64],  mm0

    SSE_SUM_SUBD           xmm3,      xmm0       ;tmp0 tmp2
    SSE_SUM_SUBD           xmm2,      xmm1       ;tmp1 tmp3

    SSE_SUM_SUBD           xmm2,      xmm3
    pslld                 xmm2,         3       ;b0
    pslld                 xmm3,         3       ;b1

    XMM_INTO_MMX         xmm4,     mm2,       mm1       
    mov                   ecx,         16
    movd                  mm0,        ecx
    movq2dq              xmm4,        mm0
    pshufd               xmm4,       xmm4,     0
    paddd                xmm2,       xmm4
    psrad                xmm2,          5
    paddd                xmm3,       xmm4
    psrad                xmm3,          5
      
    
    SSE_STORE_XMM2MMX2MM  [eax],      mm0,       xmm2   ;data[0][0:7] = xmm2
    SSE_STORE_XMM2MMX2MM  [eax+64],   mm0,       xmm3   ;data[4][0:7] = xmm3

    movdqu                xmm2,      xmm0       ;tmp2
    movdqu                xmm3,      xmm1       ;tmp3

    pslld                 xmm2,         1
    paddd                 xmm3,      xmm2
    pslld                 xmm3,         2
    paddd                 xmm3,      xmm2       ;b[2]

    pslld                 xmm1,         1
    psubd                 xmm0,      xmm1
    pslld                 xmm0,         2
    psubd                 xmm0,      xmm1       ;b[3]

    paddd                xmm3,       xmm4
    psrad                xmm3,          5
    paddd                xmm0,       xmm4
    psrad                xmm0,          5    

    MMX_INTO_XMM         xmm4,     xmm1,       mm2,       mm1 
        
    SSE_STORE_XMM2MMX2MM  [eax+32],      mm0,    xmm3   ;data[2][0:7] = b[2] = xmm3
    SSE_STORE_XMM2MMX2MM  [eax+96],      mm0,    xmm0   ;data[6][0:7] = b[3] = xmm0


    movdqu                xmm0,      xmm4
    movdqu                xmm1,      xmm5
    movdqu                xmm2,      xmm6
    movdqu                xmm3,      xmm7

    SSE_SUM_SUBD           xmm7,      xmm4      ;tmp4+tmp7, tmp4-tmp7
    SSE_SUM_SUBD           xmm6,      xmm5      ;tmp5+tmp6, tmp5-tmp6

    pslld                 xmm4,         1
    pslld                 xmm5,         1
    pslld                 xmm6,         1
    pslld                 xmm7,         1

    paddd                 xmm4,      xmm0      ;tmp0
    paddd                 xmm6,      xmm1      ;tmp1
    psubd                 xmm5,      xmm2      ;tmp2
    paddd                 xmm7,      xmm3      ;tmp3

    movdqu                xmm0,      xmm4      ;tmp0
    paddd                 xmm0,      xmm6      ;tmp0+tmp1
    paddd                 xmm0,      xmm7      ;tmp0+tmp1+tmp3
    pslld                 xmm0,         1      ;(tmp0+tmp1+tmp3)<<1
    paddd                 xmm0,      xmm6      ;b4=(tmp0+tmp1+tmp3)<<1+tmp1

    movdqu                xmm1,      xmm4      ;tmp0
    psubd                 xmm1,      xmm6      ;tmp0-tmp1
    paddd                 xmm1,      xmm5      ;tmp0-tmp1+tmp2
    pslld                 xmm1,         1      ;(tmp0-tmp1+tmp2)<<1
    paddd                 xmm1,      xmm4      ;b5=(tmp0-tmp1+tmp2)<<1+tmp0

    pxor                  xmm2,      xmm2      ;0
    psubd                 xmm2,      xmm6      ;-tmp1
    psubd                 xmm2,      xmm5      ;-tmp1-tmp2
    paddd                 xmm2,      xmm7      ;-tmp1-tmp2+tmp3
    pslld                 xmm2,         1      ;(-tmp1-tmp2+tmp3)<<1
    paddd                 xmm2,      xmm7      ;b6=(-tmp1-tmp2+tmp3)<<1+tmp3

    movdqu                xmm3,      xmm4      ;tmp0
    psubd                 xmm3,      xmm5      ;tmp0-tmp2
    psubd                 xmm3,      xmm7      ;tmp0-tmp2-tmp3
    pslld                 xmm3,         1      ;(tmp0-tmp2-tmp3)<<1
    psubd                 xmm3,      xmm5      ;b7=(tmp0-tmp2-tmp3)<<1-tmp2
    
    mov                   ecx,         16
    movd                  mm0,        ecx
    movq2dq              xmm4,        mm0
    pshufd               xmm4,       xmm4,     0
    paddd                xmm0,       xmm4
    paddd                xmm1,       xmm4
    paddd                xmm2,       xmm4
    paddd                xmm3,       xmm4
    
    psrad                xmm0,          5
    psrad                xmm1,          5
    psrad                xmm2,          5
    psrad                xmm3,          5
 
    SSE_STORE_XMM2MMX2MM  [eax+16],      mm0,    xmm0    ;data[1][0:7] = xmm0
    SSE_STORE_XMM2MMX2MM  [eax+48],      mm0,    xmm1    ;data[3][0:7] = xmm1
    SSE_STORE_XMM2MMX2MM  [eax+80],      mm0,    xmm2    ;data[5][0:7] = xmm2
    SSE_STORE_XMM2MMX2MM  [eax+112],     mm0,    xmm3    ;data[7][0:7] = xmm3
    
    %assign disp disp+8
    %endrep  
     
    ret
 
;-----------------------------------------------------------------------------
;   void __cdecl avs_1st_vidct8_sse( __int16 data[8][8] );
;-----------------------------------------------------------------------------
ALIGN 16
avs_1st_vidct8_sse:
    mov         eax, [esp+8]           ; load data[8][8] pointer
    
    ;-------------------------------------------------------------------------
    ; vertical dct ( compute 4 columns at a time -> 2 loops )
    ;-------------------------------------------------------------------------
    movdqu               xmm4,          [ eax + 16 ]       ;tmp4
    movdqu               xmm5,          [ eax + 48 ]       ;tmp5
    movdqu               xmm6,          [ eax + 80 ]       ;tmp6
    movdqu               xmm7,          [ eax + 112]       ;tmp7
    
    movdqu               xmm0,           xmm4   ;tmp4
    movdqu               xmm1,           xmm5   ;tmp5
    movdqu               xmm2,           xmm6   ;tmp6
    movdqu               xmm3,           xmm7   ;tmp7
    
    SSE_SUM_SUBW         xmm7,           xmm4   ;tmp4+tmp7  tmp4-tmp7
    psllw                xmm4,           1
    paddw                xmm0,           xmm4   ;b0
    psllw                xmm7,           1
    paddw                xmm3,           xmm7   ;b3
    
    SSE_SUM_SUBW         xmm6,           xmm5   ;tmp5+tmp6  tmp5-tmp6
    psllw                xmm6,           1
    paddw                xmm1,           xmm6   ;b1
    psllw                xmm5,           1
    psubw                xmm5,           xmm2   ;b2
    
    movdqu               xmm2,           xmm5   ;b2
    movdqu               xmm4,           xmm0   ;b0
    movdqu               xmm6,           xmm1   ;b1
    movdqu               xmm7,           xmm3   ;b3
    
    paddw                xmm1,           xmm4   ;b0+b1
    paddw                xmm1,           xmm7   ;b0+b1+b3
    psllw                xmm1,           1      ;(b0+b1+b3)<<1
    paddw                xmm1,           xmm6   ;b4
    
    paddw                xmm2,           xmm4   ;b0+b2
    psubw                xmm2,           xmm6   ;b0-b1+b2
    psllw                xmm2,           1      ;(b0-b1+b2)<<1
    paddw                xmm2,           xmm4   ;b5
    
    psubw                xmm3,           xmm5   ;-b2+b3
    psubw                xmm3,           xmm6   ;-b1-b2+b3
    psllw                xmm3,           1      ;(-b1-b2+b3)<<1
    paddw                xmm3,           xmm7   ;b6
    
    psubw                xmm0,           xmm5   ;b0-b2
    psubw                xmm0,           xmm7   ;b0-b2-b3
    psllw                xmm0,           1      ;(b0-b2-b3)<<1
    psubw                xmm0,           xmm5   ;b7
    
    
    movdqu               xmm4,          [ eax + 32 ]       ;tmp2
    movdqu               xmm5,          [ eax + 96 ]       ;tmp3    
    
    movdqu               xmm6,           xmm4  ;tmp2
    psllw                xmm6,           1     ;2tmp2
    movdqu               xmm7,           xmm6  ;2tmp2
    paddw                xmm7,           xmm5  ;2tmp2+tmp3
    psllw                xmm7,           2     ;8tmp2+4tmp3
    paddw                xmm7,           xmm6  ;tmp2=10tmp2+4tmp3      
    
    movdqu               xmm6,           xmm5  ;tmp3
    psllw                xmm6,           1     ;2tmp3
    psubw                xmm4,           xmm6  ;tmp2-2tmp3
    psllw                xmm4,           2
    psubw                xmm4,           xmm6  ;tmp3=4tmp2-10tmp3    

    movdqu               xmm5,          [ eax      ]       ;tmp0
    movdqu               xmm6,          [ eax + 64 ]       ;tmp1 
                       
    SSE_SUM_SUBW         xmm6,           xmm5  ;tmp0+tmp1   tmp0-tmp1
    psllw                xmm5,           3     ;tmp1=(tmp0-tmp1)<<3
    psllw                xmm6,           3     ;tmp0=(tmp0+tmp1)<<3
 
    ;tmp0=xmm6  tmp1=xmm5  tmp2=xmm7  tmp3=xmm4   
    SSE_SUM_SUBW         xmm7,           xmm6  ;b0=tmp0+tmp2  b3=tmp0-tmp2
    SSE_SUM_SUBW         xmm4,           xmm5  ;b1=tmp1+tmp3  b2=tmp1-tmp3
    
    ;b0=xmm7  b1=xmm4  b2=xmm5  b3=xmm6
    ;b4=xmm1  b5=xmm2  b6=xmm3  b7=xmm0 
    SSE_SUM_SUBW         xmm1,           xmm7  ;b0+b4  b0-b4  
    
    XMM_INTO_MMX         xmm2,           mm1,      mm0  ;b5=[mm1,mm0] ;store b5 
      
    mov                  ecx,            40004h  ;0x0404
    movd                 mm2,            ecx
    movq2dq              xmm2,           mm2
    pshufd               xmm2,           xmm2,     0
    paddw                xmm1,           xmm2
    paddw                xmm7,           xmm2
    psraw                xmm1,           3              ;curr_blk[0][..]
    psraw                xmm7,           3              ;curr_blk[7][..]
     
    movdqu              [eax],           xmm1           ;data[0][0:7] = xmm1
    movdqu              [eax+112],       xmm7           ;data[7][0:7] = xmm7
    
    movdqu               xmm1,           xmm2           ;xmm1=[4,4,4,4,4,4,4,4]
       
    MMX_INTO_XMM         xmm2,    xmm7,    mm1,      mm0  ;b2=xmm5 ;restore b5  
    
    SSE_SUM_SUBW         xmm2,           xmm4  ;b1+b5  b1-b5
    paddw                xmm2,           xmm1
    paddw                xmm4,           xmm1
    psraw                xmm2,           3
    psraw                xmm4,           3 
    movdqu              [eax+16],        xmm2           ;data[1][0:7] = xmm2
    movdqu              [eax+96],        xmm4           ;data[6][0:7] = xmm4    
    
    SSE_SUM_SUBW         xmm3,           xmm5  ;b2+b6  b2-b6
    paddw                xmm3,           xmm1
    paddw                xmm5,           xmm1
    psraw                xmm3,           3
    psraw                xmm5,           3
    movdqu              [eax+32],        xmm3           ;data[2][0:7] = xmm3
    movdqu              [eax+80],        xmm5           ;data[5][0:7] = xmm5    
    
    SSE_SUM_SUBW         xmm0,           xmm6  ;b3+b7  b3-b7   
    paddw                xmm0,           xmm1
    paddw                xmm6,           xmm1
    psraw                xmm0,           3
    psraw                xmm6,           3
    movdqu              [eax+48],        xmm0           ;data[3][0:7] = xmm0
    movdqu              [eax+64],        xmm6           ;data[4][0:7] = xmm6       
    ret   

;-----------------------------------------------------------------------------
;   void __cdecl avs_2nd_vidct8_sse( __int16 data[8][8] );
;-----------------------------------------------------------------------------
ALIGN 16
avs_2nd_vidct8_sse:
    mov         eax, [esp+8]           ; load data[8][8] pointer
    
    ;-------------------------------------------------------------------------
    ; vertical dct ( compute 4 columns at a time -> 2 loops )
    ;-------------------------------------------------------------------------
    movdqu               xmm4,          [ eax + 16 ]       ;tmp4
    movdqu               xmm5,          [ eax + 48 ]       ;tmp5
    movdqu               xmm6,          [ eax + 80 ]       ;tmp6
    movdqu               xmm7,          [ eax + 112]       ;tmp7
    
    movdqu               xmm0,           xmm4   ;tmp4
    movdqu               xmm1,           xmm5   ;tmp5
    movdqu               xmm2,           xmm6   ;tmp6
    movdqu               xmm3,           xmm7   ;tmp7
    
    SSE_SUM_SUBW         xmm7,           xmm4   ;tmp4+tmp7  tmp4-tmp7
    psllw                xmm4,           1
    paddsw                xmm0,           xmm4   ;b0
    psllw                xmm7,           1
    paddsw                xmm3,           xmm7   ;b3
    
    SSE_SUM_SUBW         xmm6,           xmm5   ;tmp5+tmp6  tmp5-tmp6
    psllw                xmm6,           1
    paddsw                xmm1,           xmm6   ;b1
    psllw                xmm5,           1
    psubsw                xmm5,           xmm2   ;b2
    
    movdqu               xmm2,           xmm5   ;b2
    movdqu               xmm4,           xmm0   ;b0
    movdqu               xmm6,           xmm1   ;b1
    movdqu               xmm7,           xmm3   ;b3
    
    paddsw                xmm1,           xmm4   ;b0+b1
    paddsw                xmm1,           xmm7   ;b0+b1+b3
    psllw                xmm1,           1      ;(b0+b1+b3)<<1
    paddsw                xmm1,           xmm6   ;b4
    
    paddsw                xmm2,           xmm4   ;b0+b2
    psubsw                xmm2,           xmm6   ;b0-b1+b2
    psllw                xmm2,           1      ;(b0-b1+b2)<<1
    paddsw                xmm2,           xmm4   ;b5
    
    psubsw                xmm3,           xmm5   ;-b2+b3
    psubsw                xmm3,           xmm6   ;-b1-b2+b3
    psllw                xmm3,           1      ;(-b1-b2+b3)<<1
    paddsw                xmm3,           xmm7   ;b6
    
    psubsw                xmm0,           xmm5   ;b0-b2
    psubsw                xmm0,           xmm7   ;b0-b2-b3
    psllw                xmm0,           1      ;(b0-b2-b3)<<1
    psubsw                xmm0,           xmm5   ;b7
    
    
    movdqu               xmm4,          [ eax + 32 ]       ;tmp2
    movdqu               xmm5,          [ eax + 96 ]       ;tmp3    
    
    movdqu               xmm6,           xmm4  ;tmp2
    psllw                xmm6,           1     ;2tmp2
    movdqu               xmm7,           xmm6  ;2tmp2
    paddsw                xmm7,           xmm5  ;2tmp2+tmp3
    psllw                xmm7,           2     ;8tmp2+4tmp3
    paddsw                xmm7,           xmm6  ;tmp2=10tmp2+4tmp3      
    
    movdqu               xmm6,           xmm5  ;tmp3
    psllw                xmm6,           1     ;2tmp3
    psubsw                xmm4,           xmm6  ;tmp2-2tmp3
    psllw                xmm4,           2
    psubsw                xmm4,           xmm6  ;tmp3=4tmp2-10tmp3    

    movdqu               xmm5,          [ eax      ]       ;tmp0
    movdqu               xmm6,          [ eax + 64 ]       ;tmp1 
                       
    SSE_SUM_SUBW         xmm6,           xmm5  ;tmp0+tmp1   tmp0-tmp1
    psllw                xmm5,           3     ;tmp1=(tmp0-tmp1)<<3
    psllw                xmm6,           3     ;tmp0=(tmp0+tmp1)<<3
 
    ;tmp0=xmm6  tmp1=xmm5  tmp2=xmm7  tmp3=xmm4   
    SSE_SUM_SUBW         xmm7,           xmm6  ;b0=tmp0+tmp2  b3=tmp0-tmp2
    SSE_SUM_SUBW         xmm4,           xmm5  ;b1=tmp1+tmp3  b2=tmp1-tmp3
    
    ;b0=xmm7  b1=xmm4  b2=xmm5  b3=xmm6
    ;b4=xmm1  b5=xmm2  b6=xmm3  b7=xmm0 
    
    
    XMM_INTO_MMX         xmm2,           mm1,      mm0  ;b5=[mm1,mm0] ;store b5 
     
    mov                  ecx,            400040h  ;0x4040
    movd                 mm2,            ecx
    movq2dq              xmm2,           mm2
    pshufd               xmm2,           xmm2,     0
 
    ;--------------- 200803226 in case of 16 bits overflow------------   
    ;SSE_SUM_SUBW         xmm1,           xmm7  ;b0+b4  b0-b4   
    ;paddsw                xmm1,           xmm2
    ;paddsw                xmm7,           xmm2
     XMM_INTO_MMX         xmm3,           mm3,      mm2  ;b6=[mm3,mm2] ;store b6 
     SSE_SUM_SUBW_ROUND    xmm1,      xmm7,     xmm2,    xmm3
    ;-----------------------------------------------------------------  
       
    psraw                xmm1,           7              ;curr_blk[0][..]
    psraw                xmm7,           7              ;curr_blk[7][..]
     
    movdqu              [eax],           xmm1           ;data[0][0:7] = xmm1
    movdqu              [eax+112],       xmm7           ;data[7][0:7] = xmm7
    
    movdqu               xmm1,           xmm2           ;xmm1=[4,4,4,4,4,4,4,4]
       
    MMX_INTO_XMM         xmm2,    xmm7,    mm1,      mm0  ;b2=xmm5 ;restore b5  
    MMX_INTO_XMM         xmm3,    xmm7,    mm3,      mm2  ;b6=[mm3,mm2] ;restore b6 
    ;--------------- 200803226 in case of 16 bits overflow------------
    ;SSE_SUM_SUBW         xmm2,           xmm4     ;b1+b5  b1-b5    
    ;paddsw                xmm2,           xmm1
    ;paddsw                xmm4,           xmm1  
    SSE_SUM_SUBW_ROUND    xmm2,      xmm4,     xmm1,    xmm7
    ;-----------------------------------------------------------------  
    psraw                xmm2,           7
    psraw                xmm4,           7 
    movdqu              [eax+16],        xmm2           ;data[1][0:7] = xmm2
    movdqu              [eax+96],        xmm4           ;data[6][0:7] = xmm4    
    
    
    ;--------------- 200803226 in case of 16 bits overflow------------
    ;SSE_SUM_SUBW         xmm3,           xmm5  ;b2+b6  b2-b6
    ;paddsw                xmm3,           xmm1
    ;paddsw                xmm5,           xmm1
    SSE_SUM_SUBW_ROUND    xmm3,      xmm5,     xmm1,    xmm7
    ;----------------------------------------------------------------- 
    psraw                xmm3,           7
    psraw                xmm5,           7
    movdqu              [eax+32],        xmm3           ;data[2][0:7] = xmm3
    movdqu              [eax+80],        xmm5           ;data[5][0:7] = xmm5    
  
  
    ;--------------- 200803226 in case of 16 bits overflow------------  
    ;SSE_SUM_SUBW         xmm0,           xmm6  ;b3+b7  b3-b7   
    ;paddsw                xmm0,           xmm1
    ;paddsw                xmm6,           xmm1
    SSE_SUM_SUBW_ROUND    xmm0,      xmm6,     xmm1,    xmm7
    ;-----------------------------------------------------------------
    psraw                xmm0,           7
    psraw                xmm6,           7
    movdqu              [eax+48],        xmm0           ;data[3][0:7] = xmm0
    movdqu              [eax+64],        xmm6           ;data[4][0:7] = xmm6       
    ret   
    
;-----------------------------------------------------------------------------
;   void __cdecl transpose_sse( int16_t dest[8][8] );
;-----------------------------------------------------------------------------
ALIGN 16
transpose_sse:
    mov         eax, [esp+8]           ; load data[8][8] pointer
    
    movq  mm0, [eax    ]
    movq  mm1, [eax+ 16]
    movq  mm2, [eax+ 32]
    movq  mm3, [eax+ 48]
    TRANSPOSE_SSE  mm0, mm1, mm2, mm3, mm4
    movq  [eax    ], mm0
    movq  [eax+ 16], mm3
    movq  [eax+ 32], mm4
    movq  [eax+ 48], mm2

    movq  mm0, [eax+ 72]
    movq  mm1, [eax+ 88]
    movq  mm2, [eax+104]
    movq  mm3, [eax+120]
    TRANSPOSE_SSE  mm0, mm1, mm2, mm3, mm4
    movq  [eax+ 72], mm0
    movq  [eax+ 88], mm3
    movq  [eax+104], mm4
    movq  [eax+120], mm2

    movq  mm0, [eax+  8]
    movq  mm1, [eax+ 24]
    movq  mm2, [eax+ 40]
    movq  mm3, [eax+ 56]
    TRANSPOSE_SSE  mm0, mm1, mm2, mm3, mm4
    movq  mm1, [eax+ 64]
    movq  mm5, [eax+ 80]
    movq  mm6, [eax+ 96]
    movq  mm7, [eax+112]

    movq  [eax+ 64], mm0
    movq  [eax+ 80], mm3
    movq  [eax+ 96], mm4
    movq  [eax+112], mm2
    TRANSPOSE_SSE  mm1, mm5, mm6, mm7, mm4
    movq  [eax+  8], mm1
    movq  [eax+ 24], mm7
    movq  [eax+ 40], mm4
    movq  [eax+ 56], mm6
 

    ret
    
;-----------------------------------------------------------------------------
;   void __cdecl avs_dct_sse( __int16 data[8][8])
;-----------------------------------------------------------------------------
ALIGN 16
cglobal avs_dct_sse
avs_dct_sse:
    call transpose_sse
    call avs_1st_vdct8_sse
    call transpose_sse
    call avs_2nd_vdct8_sse
    emms
    ret
    
;-----------------------------------------------------------------------------
;   void __cdecl avs_idct_sse( __int16 data[8][8])
;-----------------------------------------------------------------------------
ALIGN 16
cglobal avs_idct_sse
avs_idct_sse:
    call transpose_sse
    call avs_1st_vidct8_sse
    call transpose_sse
    call avs_2nd_vidct8_sse
    emms
    ret
    
;-----------------------------------------------------------------------------
;   void __cdecl avs_idct_sse( __int16 data[8][8])
;-----------------------------------------------------------------------------
ALIGN 16
cglobal avs_idct_sse2
avs_idct_sse2:
    call transpose_sse
    call avs_1st_vidct8_sse
    call transpose_sse
    ;call avs_2nd_vidct8_sse
    emms
    ret