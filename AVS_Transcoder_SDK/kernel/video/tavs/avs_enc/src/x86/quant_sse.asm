;*****************************************************************************
;* quant_sse.asm: avs encoder library
;* by xzhao
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

%macro SSE_ABS_AND_SIGN_W 2
    pxor        %2,  %2
    pcmpgtw     %2,  %1
    pxor        %1,  %2
    psubw       %1,  %2
%endmacro

%macro SSE_COMBINE_D2W 3
    pshufd      %3,  %1,  144
    pshufd      %1,  %1,  216
    punpckhwd   %3,  %1
    pshufd      %3,  %3,  216  
     
    pshufd      %1,  %2,  144
    pshufd      %2,  %2,  216    
    punpckhwd   %1,  %2
    pshufd      %1,  %1,  141 
      
    por         %1,  %3
%endmacro


;=============================================================================
; Local Data (Read Only)
;=============================================================================

SECTION .rodata
;-----------------------------------------------------------------------------
; Scale Matrix & Quantization Table
;-----------------------------------------------------------------------------
ALIGN 16
avs_scale_matrix:  dw   32768,37958,36158,37958,37958,43969,41884,43969,36158,41884,39898,41884,37958,43969,41884,43969
avs_quant_table:   dw   32768,29775,27554,25268,23170,21247,19369,17770,
               dw   16302,15024,13777,12634,11626,10624,9742,8958,
               dw   8192,7512,6889,6305,5793,5303,4878,4467,
               dw   4091,3756,3444,3161,2894,2654,2435,2235,
               dw   2048,1878,1722,1579,1449,1329,1218,1117,
               dw   1024,939,861,790,724,664,609,558,
               dw   512,470,430,395,362,332,304,279,
               dw   256,235,215,197,181,166,152,140
avs_qp_const:      dw   5285,5285,5285,5285,10570,10570,10570,10570
rounding_scale     dw   4,4,4,4,4,4,4,4

avs_iq_table:      dw   32768,36061,38968,42495,46341,50535,55437,60424,
                 dw   32932,35734,38968,42495,46177,50535,55109,59933,    
               dw   65535,35734,38968,42577,46341,50617,55027,60097,
               dw   32809,35734,38968,42454,46382,50576,55109,60056,
               dw   65535,35734,38968,42495,46320,50515,55109,60076,
               dw   65535,35744,38968,42495,46341,50535,55099,60087,
               dw   65535,35734,38973,42500,46341,50535,55109,60097,
               dw   32771,35734,38965,42497,46341,50535,55109,60099

avs_iq_shift:      dw   15,15,15,15,15,15,15,15,
                 dw  14,14,14,14,14,14,14,14,
               dw  14,13,13,13,13,13,13,13,
               dw  12,12,12,12,12,12,12,12,
               dw  12,11,11,11,11,11,11,11,
               dw  11,10,10,10,10,10,10,10,
               dw  10,9,9,9,9,9,9,9,
               dw  8,8,8,8,8,8,8,8

;=============================================================================
; Code
;=============================================================================

SECTION .text
       
;-----------------------------------------------------------------------------
;   void __cdecl avs_quant_sse( int32 qp, int32 mode, int16_t data[8][8], );
;-----------------------------------------------------------------------------
ALIGN 16
cglobal avs_quant_sse
avs_quant_sse:
    ;mov         eax, [esp+4]           ; load qp
    ;mov         eax, [esp+8]           ; load mode
    ;mov         eax, [esp+12]          ; load data
    
    
    
    %assign disp 0
    %rep 8
    
    ; scaling "temp * ScaleM[yy&3][xx&3]"
    mov                   eax,      [esp+12]           ; load data[8][8] pointer
    movdqu                xmm0,     [eax+(disp<<1)]
    movq                   mm0,     [avs_scale_matrix+(disp&31)]
    movq2dq               xmm1,     mm0
    pshufd                xmm1,     xmm1,     68           ;[1,0,1,0]
    SSE_ABS_AND_SIGN_W    xmm0,     xmm2
    pmulhuw               xmm0,     xmm1
    
    ; rounding  " + (1<<18)"
    movdqu                xmm3,     [rounding_scale]
    paddw                 xmm0,     xmm3
    psrlw                 xmm0,     3
    movdqa                xmm1,     xmm0
    
    ; quantization "*Q_TAB[qp]"
    mov                   eax,      [esp+4]             ; load qp
    shl                   eax,      1
    mov                    cx,      [avs_quant_table+eax]
    mov                    dx,      cx
    shl                   ecx,      16
    mov                    cx,      dx
    movd                  mm0,      ecx
    movq2dq              xmm3,      mm0
    pshufd               xmm3,      xmm3,      0
    
    mov                   eax,      [esp+8]             ; load mode
    shl                   eax,      1
    pxor                 xmm4,      xmm4
    mov                    cx,      [avs_qp_const+eax]
    pinsrw               xmm4,      cx,       0
    pshufd               xmm4,      xmm4,     0
    
    
    pmullw               xmm0,      xmm3
    pmulhuw              xmm1,      xmm3
    movdqa               xmm3,      xmm0 
    punpcklwd            xmm0,      xmm1
    punpckhwd            xmm3,      xmm1    
    
    paddd                xmm0,      xmm4
    paddd                xmm3,      xmm4      
    
    psrld                xmm0,      15
    psrld                xmm3,      15
    
    packssdw             xmm0,      xmm3
    
    pxor                 xmm0,      xmm2
    psubw                xmm0,      xmm2
    
    mov                  eax,      [esp+12]           ; load data[8][8] pointer
    movdqu               [eax+(disp<<1)],  xmm0 
    
    %assign disp disp+8
    %endrep 
    
    emms
    ret
    
    
;-----------------------------------------------------------------------------
;   void __cdecl avs_dequant_sse( int32 qp, int16_t data[8][8], );
;-----------------------------------------------------------------------------
ALIGN 16
cglobal avs_dequant_sse
avs_dequant_sse:
    ;mov         eax, [esp+4]           ; load qp
    ;mov         eax, [esp+8]          ; load data
    
    mov                   eax,      [esp+4]             ; load qp
    shl                   eax,      1
    
    mov                    cx,      [avs_iq_table+eax]
    mov                    dx,      cx
    shl                   ecx,      16
    mov                    cx,      dx
    movd                  mm0,      ecx
    pxor                 xmm2,      xmm2
    movq2dq              xmm2,      mm0
    pshufd               xmm2,      xmm2,      0
    
    xor                   ecx,      ecx
    xor                   edx,      edx
    mov                    cx,      [avs_iq_shift+eax]
    pxor                 xmm3,      xmm3
    pxor                 xmm4,      xmm4
    mov                    dx,      cx
    sub                    cx,      2
    sub                    dx,      1
    mov                    ax,      1
    shl                    ax,      cl
    pinsrw               xmm3,      ax,       0
    pinsrw               xmm4,      dx,       0
    pshufd               xmm3,      xmm3,     0
    
    
    mov                   eax,      [esp+8]           ; load data[8][8] pointer
    %assign disp 0
    %rep 8
    
    ; val*QPI
    movdqu                xmm0,     [eax+(disp<<1)]
    SSE_ABS_AND_SIGN_W    xmm0,     xmm5
    movdqa                xmm1,     xmm0
    movdqa                xmm7,     xmm5
    movdqa                xmm6,     xmm5
    punpcklwd             xmm5,     xmm6
    punpckhwd             xmm7,     xmm6


    
    ; xmm2_w = QPI  xmm3_d = shift-2  xmm4_d = shift-1
    
    pmullw               xmm0,      xmm2
    pmulhuw              xmm1,      xmm2
    movdqa               xmm6,      xmm0 
    punpcklwd            xmm0,      xmm1
    punpckhwd            xmm6,      xmm1  

    pxor                 xmm0,      xmm5
    psubd                xmm0,      xmm5
    pxor                 xmm6,      xmm7
    psubd                xmm6,      xmm7
            
    paddd                xmm0,      xmm3
    paddd                xmm6,      xmm3
    
    psrad                xmm0,      xmm4
    psrad                xmm6,      xmm4
    
    packssdw             xmm0,      xmm6
    

    
    movdqu               [eax+(disp<<1)],  xmm0 
    
    %assign disp disp+8
    %endrep 
    
    emms
    ret