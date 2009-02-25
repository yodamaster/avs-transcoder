#include <stdio.h>
#include <math.h>
#include "global.h"



// Initiate rate control parameters
void c_avs_enc:: rc_init_seq()
{
  double L1,L2,L3,bpp;
  int_32_t qp;
  Xb=0;
  bit_rate   = input->bit_rate;
  frame_rate = img->framerate;
  bit_rate_per_frame = bit_rate / frame_rate;

  PreviousBit_Rate = bit_rate;

  if(input->basicunit > img->total_number_mb)
  {
    input->basicunit = img->total_number_mb;
    TotalNumberofBasicUnit = 0;
  }
  else if(input->basicunit < img->total_number_mb)
  {
    TotalNumberofBasicUnit = img->total_number_mb/input->basicunit;
  }

  MINVALUE=4.0;
  /*initialize the parameters of fluid flow traffic model*/

  BufferSize = bit_rate*2.56;
  CurrentBufferFullness = 0;
  GOPTargetBufferLevel=0;
  /*HRD consideration*/
  InitialDelayOffset = BufferSize*0.8;

  /*initialize the previous window size*/
  m_windowSize    = 0;
  MADm_windowSize = 0;
  img->NumberofCodedBFrame = 0;
  img->NumberofCodedPFrame = 0;
  img->NumberofGOP = 0;
  /*remaining bits in GOP */
  R = 0;
  /*control parameter */
  if(input->successive_Bframe>0)
  {
    BETAP=0.9;
  }
  else
  {
    GAMMAP=0.5;
    BETAP =0.5;
  }

  /*quadratic rate-distortion model*/
  PPreHeader = 0;
  Pm_X1 = bit_rate;
  Pm_X2 = 0.0;
  /* linear prediction model for P picture*/
  PMADPictureC1 = 1.0;
  PMADPictureC2 = 0.0;

  memset(PP_frm_Qstep,     0, sizeof(PP_frm_Qstep));
  memset(PP_frm_R,         0, sizeof(PP_frm_R));
  memset(PPictureMAD,      0, sizeof(PPictureMAD));

  //Define the largest variation of quantization parameters
  PDuantQp = 2;

  /*basic unit layer rate control*/
  PAveHeaderBits1=0;
  PAveHeaderBits3=0;
  if(TotalNumberofBasicUnit>=9)
    DDquant = 1;
  else
    DDquant = 2;
  /*adaptive field/frame coding*/
  img->FieldControl = 0;

  RC_MAX_QUANT = 63;// clipping  42 for LBC    35 for HD
  RC_MIN_QUANT = 0;// clipping  32 for LBC    22 for HD
  if (input->SeinitialQP==0)
  {
    /*compute the initial QP*/
    bpp = (float)bit_rate /(frame_rate*img->width*img->height);
    if (img->width == 176)
    {
      L1 = 0.1;
      L2 = 0.3;
      L3 = 0.6;
    }
    else if (img->width == 352)
    {
      L1 = 0.2;
      L2 = 0.6;
      L3 = 1.2;
    }
    else
    {
      L1 = 0.6;
      L2 = 1.4;
      L3 = 2.4;
    }

    if(bpp<= L1)
      qp = 35;
    else
    {
      if(bpp<=L2)
        qp = 25;
      else
      {
        if(bpp<=L3)
          qp = 20;
        else
          qp = 10;
      }
    }
    input->SeinitialQP = qp;
  }
  TargetBufferLevel = 0;
  PreviousFrameMAD  = 0;
  PreviousQp1 = PreviousQp2 = 0;
  PAverageQp = 0;
  AWp = AWb = 1;
}

// Initiate one GOP
void c_avs_enc::rc_init_GOP(int_32_t np, int_32_t nb)
{
  int_32_t AllocatedBits;
  int_32_t GOPDquant;

  /*initialize the lower bound and the upper bound for the target bits of each frame, HRD consideration*/
  LowerBound  = (long)(R+bit_rate/frame_rate);
  UpperBound1 = (long)(R+InitialDelayOffset);

  /*compute the total number of bits for the current GOP*/
  AllocatedBits = (int_32_t) floor((1 + np + nb) * bit_rate / frame_rate + 0.5);
  R  += AllocatedBits;
  Np  = np;
  Nb  = nb;

  GOPOverdue  = myfalse;

  /*field coding*/
  img->IFLAG=1;

  /*Compute InitialQp for each GOP*/
  TotalPFrame = np;
  img->NumberofGOP++;
  if(img->NumberofGOP==1)
  {
    MyInitialQp = input->SeinitialQP;
    PreviousQp2 = MyInitialQp-1;
    QPLastGOP   = MyInitialQp;
  }
  else
  {
    /*adaptive field/frame coding*/
    if((input->InterlaceCodingOption==2))
    {
      if (img->FieldFrame == 1)
      {
        img->TotalQpforPPicture += FrameQPBuffer;
        QPLastPFrame = FrameQPBuffer;
      }
      else
      {
        img->TotalQpforPPicture += FieldQPBuffer;
        QPLastPFrame = FieldQPBuffer;
      }

    }
    /*compute the average QP of P frames in the previous GOP*/
    PAverageQp=(int_32_t)(1.0*img->TotalQpforPPicture/img->NumberofPPicture+0.5);

    GOPDquant=(int_32_t)(0.5+1.0*(np+nb+1)/15);
    if(GOPDquant>2)
      GOPDquant=2;
    PAverageQp-=GOPDquant;
    if (PAverageQp > (QPLastPFrame - 2))
      PAverageQp--;
    PAverageQp = max(QPLastGOP-2,  PAverageQp);
    PAverageQp = min(QPLastGOP+2, PAverageQp);
    PAverageQp = min(RC_MAX_QUANT, PAverageQp);
    PAverageQp = max(RC_MIN_QUANT, PAverageQp);

    MyInitialQp=PAverageQp;
    QPLastGOP = MyInitialQp;
    Pm_Qp = PAverageQp;
    PAveFrameQP = PAverageQp;
    PreviousQp1 = PreviousQp2;
    PreviousQp2 = MyInitialQp-1;
  }

  img->TotalQpforPPicture = 0;
  img->NumberofPPicture   = 0;
}

void c_avs_enc::rc_init_pict(int_32_t fieldpic,int_32_t topfield,int_32_t targetcomputation)
{
  int_32_t i;

  /*compute the total number of basic units in a frame*/
  img->NumberofCodedMacroBlocks = 0;

  /*Normally, the bandwidth for the VBR case is estimated by
  a congestion control algorithm. A bandwidth curve can be predefined if we only want to
  test the proposed algorithm*/
  if(input->channel_type==1)
  {
    if(img->NumberofCodedPFrame==58)
      bit_rate *=1.5;
    else if(img->NumberofCodedPFrame==59)
      PreviousBit_Rate=bit_rate;
  }

  /*predefine a target buffer level for each frame*/
  if((fieldpic||topfield) && targetcomputation)
  {
    switch (img->type)
    {
    case INTER_IMG:
      /*Since the available bandwidth may vary at any time, the total number of bits is updated picture by picture*/
      if(PreviousBit_Rate != bit_rate)
        R +=(int_32_t) floor((bit_rate-PreviousBit_Rate)*(Np+Nb)/frame_rate+0.5);

      /* predefine the target buffer level for each picture. frame layer rate control*/
      if(img->BasicUnit==img->total_number_mb)
      {
        if(img->NumberofPPicture==1)
        {
          TargetBufferLevel = CurrentBufferFullness;
          DeltaP = (CurrentBufferFullness-GOPTargetBufferLevel)/(TotalPFrame-1);
          TargetBufferLevel -= DeltaP;
        }
        else if(img->NumberofPPicture>1)
          TargetBufferLevel -= DeltaP;
      }
      /*basic unit layer rate control*/
      else
      {
        if(img->NumberofCodedPFrame>0)
        {
          /*adaptive frame/filed coding*/
          if(((input->InterlaceCodingOption==PAFF_CODING))&&(img->FieldControl==1))
          {
            for(i=0; i<TotalNumberofBasicUnit; i++)
              FCBUPFMAD[i] = FCBUCFMAD[i];
          }
          else
          {
            for(i=0; i<TotalNumberofBasicUnit; i++)
              BUPFMAD[i] = BUCFMAD[i];
          }
        }

        if(img->NumberofGOP==1)
        {
          if(img->NumberofPPicture==1)
          {
            TargetBufferLevel=CurrentBufferFullness;
            DeltaP=(CurrentBufferFullness-GOPTargetBufferLevel)/(TotalPFrame-1);
            TargetBufferLevel -=DeltaP;
          }
          else if(img->NumberofPPicture>1)
            TargetBufferLevel -=DeltaP;
        }
        else if(img->NumberofGOP>1)
        {
          if(img->NumberofPPicture==0)
          {
            TargetBufferLevel=CurrentBufferFullness;
            DeltaP=(CurrentBufferFullness-GOPTargetBufferLevel)/TotalPFrame;
            TargetBufferLevel -=DeltaP;
          }
          else if(img->NumberofPPicture>0)
            TargetBufferLevel -=DeltaP;
        }
      }

      if(img->NumberofCodedPFrame==1)
        AWp = Wp;
      else if((img->NumberofCodedPFrame<8)&&(img->NumberofCodedPFrame>1))
        AWp = Wp*(img->NumberofCodedPFrame-1)/img->NumberofCodedPFrame+AWp/img->NumberofCodedPFrame;
      else if(img->NumberofCodedPFrame>1)
        AWp = Wp/8+7*AWp/8;

      //compute the average complexity of B frames
      if(input->successive_Bframe>0)
      {
        //compute the target buffer level
        TargetBufferLevel +=(AWp*(input->successive_Bframe+1)*bit_rate/(frame_rate*(AWp+AWb*input->successive_Bframe))-bit_rate/frame_rate);
      }
      break;
    case B_IMG:
      /* update the total number of bits if the bandwidth is changed*/
      if(PreviousBit_Rate!=bit_rate)
        R +=(int_32_t) floor((bit_rate-PreviousBit_Rate)*(Np+Nb)/frame_rate+0.5);
      if((img->NumberofCodedPFrame==1)&&(img->NumberofCodedBFrame==1))
      {
        AWp = Wp;
        AWb = Wb;
      }
      else if(img->NumberofCodedBFrame>1)
      {
        //compute the average weight
        if(img->NumberofCodedBFrame<8)
          AWb=Wb*(img->NumberofCodedBFrame-1)/img->NumberofCodedBFrame+AWb/img->NumberofCodedBFrame;
        else
          AWb=Wb/8+7*AWb/8;
      }
      break;
    }
    /*Compute the target bit for each frame*/
    if(img->type==INTER_IMG)
    {
      /*frame layer rate control*/
      if(img->BasicUnit==img->total_number_mb)
      {
        if(img->NumberofCodedPFrame>0)
        {
          T  = (long) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);
          T1 = (long) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1 =  max(0,T1);
          T  = (long)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
      }
      /*basic unit layer rate control*/
      else
      {
        if((img->NumberofGOP==1)&&(img->NumberofCodedPFrame>0))
        {
          T = (int_32_t) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);
          T1 = (int_32_t) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1=max(0,T1);
          T = (int_32_t)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
        else if(img->NumberofGOP>1)
        {
          T = (long) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);
          T1 = (long) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1 = max(0,T1);
          T = (long)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
      }

      /*reserve some bits for smoothing*/

      /*HRD consideration*/
      T = max(T, LowerBound);
      T = min(T, UpperBound2);

      if((topfield)||(fieldpic&&((input->InterlaceCodingOption==2))))
        T_field=T;
    }
  }

  if(fieldpic||topfield)
  {
    /*frame layer rate control*/
    img->NumberofHeaderBits  = 0;
    img->NumberofTextureBits = 0;

    /*basic unit layer rate control*/
    if(img->BasicUnit < img->total_number_mb)
    {
      TotalFrameQP=0;
      img->NumberofBasicUnitHeaderBits  = 0;
      img->NumberofBasicUnitTextureBits = 0;
      img->TotalMADBasicUnit = 0;
      if(img->FieldControl == 0)
        NumberofBasicUnit=TotalNumberofBasicUnit;
      else
        NumberofBasicUnit=TotalNumberofBasicUnit/2;
    }
  }

  if((img->type==INTER_IMG)&&(img->BasicUnit<img->total_number_mb)&&(img->FieldControl==1))
  {
    /*top filed at basic unit layer rate control*/
    if(topfield)
    {
      bits_topfield=0;
      T=(long)(T_field*0.6);
    }
    /*bottom filed at basic unit layer rate control*/
    else
    {
      T=T_field-bits_topfield;
      img->NumberofBasicUnitHeaderBits=0;
      img->NumberofBasicUnitTextureBits=0;
      img->TotalMADBasicUnit=0;
      NumberofBasicUnit=TotalNumberofBasicUnit/2;
    }
  }
}

//calculate MAD for the current macroblock
double c_avs_enc::calc_MAD()
{
  int_32_t k,l;
  int_32_t s = 0;
  double MAD;

  for (k = 0; k < 16; k++)
  {
    for (l = 0; l < 16; l++)
    {
      s += abs(diffy[k][l]);
    }
  }

  MAD = s*1.0/256;
  return MAD;
}

// update one picture after frame/field encoding
void c_avs_enc::rc_update_pict(int_32_t nbits)
{
  R -= nbits; /* remaining bits in GOP */
  //printf("remaining bits:%4d\n", R);
  CurrentBufferFullness += (nbits - bit_rate_per_frame);

  /*update the lower bound and the upper bound for the target bits of each frame, HRD consideration*/
  LowerBound  += (long)(bit_rate_per_frame-nbits);
  UpperBound1 += (long)(bit_rate_per_frame-nbits);
  UpperBound2  = (long)(OMEGA*UpperBound1);

  return;
}

// update after frame encoding
void c_avs_enc::rc_update_pict_frame(int_32_t nbits)
{
  /*update the complexity weight of I, P, B frame*/
  int_32_t Avem_Qc;
  int_32_t X;
  if (img->type == INTRA_IMG)
  {
    return;
  }
  /*frame layer rate control*/
  if(img->BasicUnit==img->total_number_mb)
  {
    X = (int_32_t) floor(nbits*m_Qc+ 0.5);
  }
  else
  {
    if(img->type==INTER_IMG)
    {
      if(((img->IFLAG==0)&&(img->FieldControl==1))||(img->FieldControl==0))
      {
        Avem_Qc=TotalFrameQP/TotalNumberofBasicUnit;
        X=(int_32_t)floor(nbits*Avem_Qc+0.5);
      }
    }
    else if(img->type==B_IMG)
    {
      X = (int_32_t) floor(nbits*m_Qc+ 0.5);
    }
  }
  switch (img->type)
  {
  case INTER_IMG:
    if(((img->IFLAG==0)&&(img->FieldControl==1))||(img->FieldControl==0))
    {
      Np--;
      Wp=X;
      img->NumberofCodedPFrame++;
      img->NumberofPPicture++;
    }
    else if((img->IFLAG!=0)&&(img->FieldControl==1))
    {
      img->IFLAG=0;
    }
    break;
  case B_IMG:
    Xb = X;
    Nb--;
    Wb = Xb/THETA;
    img->NumberofCodedBFrame++;
    break;
  }
}

// coded bits for top field
void c_avs_enc::setbitscount(int_32_t nbits)
{
  bits_topfield = nbits;
}

//compute a  quantization parameter for each frame
int_32_t c_avs_enc::updateQuantizationParameter(int_32_t topfield)
{
  double dtmp;
  int_32_t m_Bits;
  int_32_t BFrameNumber;
  int_32_t StepSize;
  int_32_t PAverageQP;
  int_32_t SumofBasicUnit;
  int_32_t i;

  /*frame layer rate control*/
  if(img->BasicUnit==img->total_number_mb)
  {
    /*fixed quantization parameter is used to coded I frame, the first P frame and the first B frame
    the quantization parameter is adjusted according the available channel bandwidth and
    the type of vide*/
    /*top field*/
    if((topfield)||(img->FieldControl==0))
    {
      if(img->type == INTRA_IMG)
      {
        m_Qc=MyInitialQp;
        return m_Qc;
      }
      else if(img->type==B_IMG)
      {
        if(input->successive_Bframe==1)
        {
          BFrameNumber=(Bframe_ctr+1)%input->successive_Bframe;
          if(BFrameNumber==0)
            BFrameNumber=input->successive_Bframe;
          /*adaptive field/frame coding*/
          else if(BFrameNumber==1)
          {
            if((input->InterlaceCodingOption==PAFF_CODING))
            {
              if(img->FieldControl==0)
              {
                /*previous choice is frame coding*/
                if(img->FieldFrame==1)
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FrameQPBuffer;
                }
                /*previous choice is field coding*/
                else
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FieldQPBuffer;
                }
              }
            }
          }
          if(PreviousQp1==PreviousQp2)
            m_Qc=PreviousQp1+2;
          else
            m_Qc=(PreviousQp1+PreviousQp2)/2+1;
          m_Qc = min(m_Qc, RC_MAX_QUANT); // clipping
          m_Qc = max(RC_MIN_QUANT, m_Qc);//clipping
        }
        else
        {
          BFrameNumber = (Bframe_ctr+1)%input->successive_Bframe;
          if(BFrameNumber==0) // the last b frame
            BFrameNumber = input->successive_Bframe;
          /*adaptive field/frame coding*/
          else if(BFrameNumber==1) // the first b frame
          {
            if((input->InterlaceCodingOption==PAFF_CODING))
            {
              if(img->FieldControl==0)
              {
                /*previous choice is frame coding*/
                if(img->FieldFrame==1)
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FrameQPBuffer;
                }
                /*previous choice is field coding*/
                else
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FieldQPBuffer;
                }
              }
            }
          }

          if((PreviousQp2-PreviousQp1)<=(-2*input->successive_Bframe-3))
            StepSize = -3;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe-2))
            StepSize = -2;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe-1))
            StepSize = -1;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe))
            StepSize = 0;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe+1))
            StepSize = 1;
          else
            StepSize = 2;

          m_Qc  = PreviousQp1+StepSize;
          m_Qc += min(2*(BFrameNumber-1),max(-2*(BFrameNumber-1), (BFrameNumber-1)*(PreviousQp2-PreviousQp1)/(input->successive_Bframe-1)));
          m_Qc  = clamp(m_Qc, RC_MIN_QUANT, RC_MAX_QUANT);
        }
        return m_Qc;
      }
      else if((img->type==INTER_IMG)&&(img->NumberofPPicture==0))
      {
        m_Qc = MyInitialQp;

        if(img->FieldControl==0)
        {
          if(input->InterlaceCodingOption==0)
          {
            img->TotalQpforPPicture += m_Qc;
            PreviousQp1=PreviousQp2;
            PreviousQp2=m_Qc;
            Pm_Qp=m_Qc;
          }
          /*adaptive field/frame coding*/
          else
            FrameQPBuffer=m_Qc;
        }

        return m_Qc;
      }
      else
      {
        /*adaptive field/frame coding*/
        if(((input->InterlaceCodingOption==PAFF_CODING))&&(img->FieldControl==0))
        {
          /*previous choice is frame coding*/
          if(img->FieldFrame==1)
          {
            img->TotalQpforPPicture += FrameQPBuffer;
            Pm_Qp = FrameQPBuffer;
          }
          /*previous choice is field coding*/
          else
          {
            img->TotalQpforPPicture +=FieldQPBuffer;
            Pm_Qp=FieldQPBuffer;
          }
        }

        m_X1 = Pm_X1;
        m_X2 = Pm_X2;
        m_Hp = PPreHeader;
        m_Qp = Pm_Qp;
        DuantQp = PDuantQp;
        MADPictureC1 = PMADPictureC1;
        MADPictureC2 = PMADPictureC2;
        PreviousPictureMAD = PPictureMAD[0];
        /* predict the MAD of current picture*/
        CurrentFrameMAD = MADPictureC1*PreviousPictureMAD+MADPictureC2;
        /*compute the number of bits for the texture*/
        if(T<0)
        {
          m_Qc = m_Qp+DuantQp;
          m_Qc = min(m_Qc, RC_MAX_QUANT); // clipping
        }
        else
        {
          m_Bits = T-m_Hp;
          m_Bits = max(m_Bits, (int_32_t)(bit_rate/(MINVALUE*frame_rate)));
          dtmp   = CurrentFrameMAD * m_X1 * CurrentFrameMAD * m_X1 + 4 * m_X2 * CurrentFrameMAD * m_Bits;
          if ((m_X2 == 0.0) || (dtmp < 0) || ((sqrt (dtmp) - m_X1 * CurrentFrameMAD) <= 0.0)) // fall back 1st order mode
            m_Qstep = (float) (m_X1 * CurrentFrameMAD / (double) m_Bits);
          else // 2nd order mode
            m_Qstep = (float) ((2 * m_X2 * CurrentFrameMAD) / (sqrt (dtmp) - m_X1 * CurrentFrameMAD));

          m_Qc = Qstep2QP(m_Qstep);
          m_Qc = min(m_Qp+DuantQp, m_Qc);  // control variation
          m_Qc = min(m_Qc, RC_MAX_QUANT);  // clipping
          m_Qc = max(m_Qp-DuantQp, m_Qc);  // control variation
          m_Qc = max(RC_MIN_QUANT, m_Qc);
        }

        if(img->FieldControl==0)
        {
          /*frame coding*/
          if(input->InterlaceCodingOption == FRAME_CODING)
          {
            img->TotalQpforPPicture += m_Qc;
            PreviousQp1 = PreviousQp2;
            PreviousQp2 = m_Qc;
            Pm_Qp       = m_Qc;
          }
          else
          {
            FrameQPBuffer = m_Qc;
          }
        }
        return m_Qc;
      }
    }
    /*bottom field*/
    else
    {
      if((img->type==INTER_IMG)&&(img->IFLAG==0))
      {
        /*field coding*/
        if(input->InterlaceCodingOption==1)
        {
          img->TotalQpforPPicture +=m_Qc;
          PreviousQp1=PreviousQp2+1;
          PreviousQp2=m_Qc;
          Pm_Qp=m_Qc;
        }
        else
          FieldQPBuffer=m_Qc;
      }
      return m_Qc;
    }
  }
  /*basic unit layer rate control*/
  else
  {
    /*top filed of I frame*/
    if(img->type==INTRA_IMG)
    {
      m_Qc=MyInitialQp;
      return m_Qc;
    }
    /*bottom field of I frame*/
    else if((img->type==INTER_IMG)&&(img->IFLAG==1)&&(img->FieldControl==1))
    {
      m_Qc=MyInitialQp;
      return m_Qc;
    }
    else if(img->type==B_IMG)
    {
      /*top filed of B frame*/
      if((topfield)||(img->FieldControl==0))
      {
        if(input->successive_Bframe==1)
        {
          BFrameNumber=(Bframe_ctr+1)%input->successive_Bframe;
          if(BFrameNumber==0)
            BFrameNumber=input->successive_Bframe;
          /*adaptive field/frame coding*/
          else if(BFrameNumber==1)
          {
            if((input->InterlaceCodingOption==2))
            {
              if(img->FieldControl==0)
              {
                /*previous choice is frame coding*/
                if(img->FieldFrame==1)
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FrameQPBuffer;
                }
                /*previous choice is field coding*/
                else
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FieldQPBuffer;
                }
              }
            }
          }
          if(PreviousQp1==PreviousQp2)
            m_Qc=PreviousQp1+2;
          else
            m_Qc=(PreviousQp1+PreviousQp2)/2+1;
          m_Qc = min(m_Qc, RC_MAX_QUANT); // clipping
          m_Qc = max(RC_MIN_QUANT, m_Qc);//clipping
        }
        else
        {
          BFrameNumber=(Bframe_ctr+1)%input->successive_Bframe;
          if(BFrameNumber==0)
            BFrameNumber=input->successive_Bframe;
          /*adaptive field/frame coding*/
          else if(BFrameNumber==1)
          {
            if((input->InterlaceCodingOption==PAFF_CODING))
            {
              if(img->FieldControl==0)
              {
                /*previous choice is frame coding*/
                if(img->FieldFrame==1)
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FrameQPBuffer;
                }
                /*previous choice is field coding*/
                else
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FieldQPBuffer;
                }
              }
            }
          }

          if((PreviousQp2-PreviousQp1)<=(-2*input->successive_Bframe-3))
            StepSize=-3;
          else  if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe-2))
            StepSize=-2;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe-1))
            StepSize=-1;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe))
            StepSize=0;//0
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe+1))
            StepSize=1;//1
          else
            StepSize=2;//2
          m_Qc=PreviousQp1+StepSize;
          m_Qc += min(2*(BFrameNumber-1),max(-2*(BFrameNumber-1), (BFrameNumber-1)*(PreviousQp2-PreviousQp1)/(input->successive_Bframe-1)));
          m_Qc = min(m_Qc, RC_MAX_QUANT); // clipping
          m_Qc = max(RC_MIN_QUANT, m_Qc);//clipping
        }
        return m_Qc;
      }
      /*bottom field of B frame*/
      else
        return m_Qc;
    }
    else if(img->type==INTER_IMG)
    {
      if((img->NumberofGOP==1)&&(img->NumberofPPicture==0))
      {
        if((img->FieldControl==0)||((img->FieldControl==1)&&(img->IFLAG==0)))
        {
          /*top field of the first P frame*/
          m_Qc=MyInitialQp;
          img->NumberofBasicUnitHeaderBits=0;
          img->NumberofBasicUnitTextureBits=0;
          NumberofBasicUnit--;
          /*bottom field of the first P frame*/
          if((!topfield)&&(NumberofBasicUnit==0))
          {
            /*frame coding or field coding*/
            if((input->InterlaceCodingOption==FRAME_CODING)||(input->InterlaceCodingOption==FIELD_CODING))
            {
              img->TotalQpforPPicture +=m_Qc;
              PreviousQp1=PreviousQp2;
              PreviousQp2=m_Qc;
              PAveFrameQP=m_Qc;
              PAveHeaderBits3=PAveHeaderBits2;
            }
            /*adaptive frame/field coding*/
            else if((input->InterlaceCodingOption==2))
            {
              if(img->FieldControl==0)
              {
                FrameQPBuffer=m_Qc;
                FrameAveHeaderBits=PAveHeaderBits2;
              }
              else
              {
                FieldQPBuffer=m_Qc;
                FieldAveHeaderBits=PAveHeaderBits2;
              }
            }
          }
          Pm_Qp=m_Qc;
          TotalFrameQP +=m_Qc;
          return m_Qc;
        }
      }
      else
      {
        m_X1=Pm_X1;
        m_X2=Pm_X2;
        m_Hp=PPreHeader;
        m_Qp=Pm_Qp;
        DuantQp=PDuantQp;
        MADPictureC1=PMADPictureC1;
        MADPictureC2=PMADPictureC2;

        if(img->FieldControl==0)
          SumofBasicUnit=TotalNumberofBasicUnit;
        else
          SumofBasicUnit=TotalNumberofBasicUnit/2;

        /*the average QP of the previous frame is used to coded the first basic unit of the current frame or field*/
        if(NumberofBasicUnit==SumofBasicUnit)
        {
          /*adaptive field/frame coding*/
          if(((input->InterlaceCodingOption==2))&&(img->FieldControl==0))
          {
            /*previous choice is frame coding*/
            if(img->FieldFrame==1)
            {
              if(img->NumberofPPicture>0)
                img->TotalQpforPPicture +=FrameQPBuffer;
              PAveFrameQP=FrameQPBuffer;
              PAveHeaderBits3=FrameAveHeaderBits;
            }
            /*previous choice is field coding*/
            else
            {
              if(img->NumberofPPicture>0)
                img->TotalQpforPPicture +=FieldQPBuffer;
              PAveFrameQP=FieldQPBuffer;
              PAveHeaderBits3=FieldAveHeaderBits;
            }
          }

          if(T<=0)
          {
            m_Qc=PAveFrameQP+2;
            if(m_Qc>RC_MAX_QUANT)
              m_Qc=RC_MAX_QUANT;
            if(topfield||(img->FieldControl==0))
              GOPOverdue=mytrue;
          }
          else
          {
            m_Qc=PAveFrameQP;
          }
          TotalFrameQP +=m_Qc;
          NumberofBasicUnit--;
          Pm_Qp=PAveFrameQP;
          return m_Qc;
        }
        else
        {
          /*compute the number of remaining bits*/
          TotalBasicUnitBits=img->NumberofBasicUnitHeaderBits+img->NumberofBasicUnitTextureBits;
          T -=TotalBasicUnitBits;
          img->NumberofBasicUnitHeaderBits=0;
          img->NumberofBasicUnitTextureBits=0;
          if(T<0)
          {
            if(GOPOverdue==mytrue)
              m_Qc=m_Qp+2;
            else
              m_Qc=m_Qp+DDquant;//2
            m_Qc = min(m_Qc, RC_MAX_QUANT);  // clipping
            if(input->basicunit>=img->img_width_in_mb)
              m_Qc = min(m_Qc, PAveFrameQP+6);
            else
              m_Qc = min(m_Qc, PAveFrameQP+3);

            TotalFrameQP +=m_Qc;
            NumberofBasicUnit--;
            if(NumberofBasicUnit==0)
            {
              if((!topfield)||(img->FieldControl==0))
              {
                /*frame coding or field coding*/
                if((input->InterlaceCodingOption==0)||(input->InterlaceCodingOption==1))
                {
                  PAverageQP=(int_32_t)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                  if (img->NumberofPPicture == (input->intra_period - 2))
                    QPLastPFrame = PAverageQP;

                  img->TotalQpforPPicture +=PAverageQP;
                  if(GOPOverdue==mytrue)
                  {
                    PreviousQp1=PreviousQp2+1;
                    PreviousQp2=PAverageQP;
                  }
                  else
                  {
                    if((img->NumberofPPicture==0)&&(img->NumberofGOP>1))
                    {
                      PreviousQp1=PreviousQp2;
                      PreviousQp2=PAverageQP;
                    }
                    else if(img->NumberofPPicture>0)
                    {
                      PreviousQp1=PreviousQp2+1;
                      PreviousQp2=PAverageQP;
                    }
                  }
                  PAveFrameQP=PAverageQP;
                  PAveHeaderBits3=PAveHeaderBits2;
                }
                /*adaptive field/frame coding*/
                else if((input->InterlaceCodingOption==2))
                {
                  if(img->FieldControl==0)
                  {
                    PAverageQP=(int_32_t)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                    FrameQPBuffer=PAverageQP;
                    FrameAveHeaderBits=PAveHeaderBits2;
                  }
                  else
                  {
                    PAverageQP=(int_32_t)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                    FieldQPBuffer=PAverageQP;
                    FieldAveHeaderBits=PAveHeaderBits2;
                  }
                }
              }
            }
            if(GOPOverdue==mytrue)
              Pm_Qp=PAveFrameQP;
            else
              Pm_Qp=m_Qc;
            return m_Qc;
          }
          else
          {
            /*predict the MAD of current picture*/
            if(((input->InterlaceCodingOption==2))&&(img->FieldControl==1))
            {
              CurrentFrameMAD=MADPictureC1*FCBUPFMAD[TotalNumberofBasicUnit-NumberofBasicUnit]+MADPictureC2;
              TotalBUMAD=0;
              for(i=TotalNumberofBasicUnit-1; i>=(TotalNumberofBasicUnit-NumberofBasicUnit);i--)
              {
                CurrentBUMAD=MADPictureC1*FCBUPFMAD[i]+MADPictureC2;
                TotalBUMAD +=CurrentBUMAD*CurrentBUMAD;
              }
            }
            else
            {
              CurrentFrameMAD=MADPictureC1*BUPFMAD[TotalNumberofBasicUnit-NumberofBasicUnit]+MADPictureC2;
              TotalBUMAD=0;
              for(i=TotalNumberofBasicUnit-1; i>=(TotalNumberofBasicUnit-NumberofBasicUnit);i--)
              {
                CurrentBUMAD=MADPictureC1*BUPFMAD[i]+MADPictureC2;
                TotalBUMAD +=CurrentBUMAD*CurrentBUMAD;
              }
            }

            /*compute the total number of bits for the current basic unit*/
            m_Bits =(int_32_t)(T*CurrentFrameMAD*CurrentFrameMAD/TotalBUMAD);
            /*compute the number of texture bits*/
            m_Bits -=PAveHeaderBits2;

            m_Bits=max(m_Bits,(int_32_t)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));

            dtmp = CurrentFrameMAD * m_X1 * CurrentFrameMAD * m_X1 + 4 * m_X2 * CurrentFrameMAD * m_Bits;
            if ((m_X2 == 0.0) || (dtmp < 0) || ((sqrt (dtmp) - m_X1 * CurrentFrameMAD) <= 0.0))  // fall back 1st order mode
              m_Qstep = (float)(m_X1 * CurrentFrameMAD / (double) m_Bits);
            else // 2nd order mode
              m_Qstep = (float) ((2 * m_X2 * CurrentFrameMAD) / (sqrt (dtmp) - m_X1 * CurrentFrameMAD));

            m_Qc=Qstep2QP(m_Qstep);
            m_Qc = min(m_Qp+DDquant,  m_Qc); // control variation

            if(input->basicunit>=img->img_width_in_mb)
              m_Qc = min(PAveFrameQP+6, m_Qc);
            else
              m_Qc = min(PAveFrameQP+3, m_Qc);

            m_Qc = min(m_Qc, RC_MAX_QUANT);  // clipping
            m_Qc = max(m_Qp-DDquant, m_Qc);  // control variation
            if(input->basicunit>=img->img_width_in_mb)
              m_Qc = max(PAveFrameQP-6, m_Qc);
            else
              m_Qc = max(PAveFrameQP-3, m_Qc);

            m_Qc = max(RC_MIN_QUANT, m_Qc);
            TotalFrameQP +=m_Qc;
            Pm_Qp=m_Qc;
            NumberofBasicUnit--;
            if((NumberofBasicUnit==0)&&(img->type==INTER_IMG))
            {
              if((!topfield)||(img->FieldControl==0))
              {
                /*frame coding or field coding*/
                if((input->InterlaceCodingOption==0)||(input->InterlaceCodingOption==1))
                {
                  PAverageQP=(int_32_t)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                  if (img->NumberofPPicture == (input->intra_period - 2))
                    QPLastPFrame = PAverageQP;

                  img->TotalQpforPPicture +=PAverageQP;
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=PAverageQP;
                  PAveFrameQP=PAverageQP;
                  PAveHeaderBits3=PAveHeaderBits2;
                }
                else if((input->InterlaceCodingOption==2))
                {
                  if(img->FieldControl==0)
                  {
                    PAverageQP=(int_32_t)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                    FrameQPBuffer=PAverageQP;
                    FrameAveHeaderBits=PAveHeaderBits2;
                  }
                  else
                  {
                    PAverageQP=(int_32_t)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                    FieldQPBuffer=PAverageQP;
                    FieldAveHeaderBits=PAveHeaderBits2;
                  }
                }
              }
            }
            return m_Qc;
          }
        }
      }
    }
  }
}

//update the parameters of quadratic R-D model
void c_avs_enc::updateRCModel ()
{
  int_32_t n_windowSize;
  int_32_t i;
  double error[20], std = 0.0, threshold;
  int_32_t m_Nc;
  myboolean MADModelFlag = myfalse;
  if(img->type==INTER_IMG)
  {
    if(img->BasicUnit==img->total_number_mb)
    {
      CurrentFrameMAD = ComputeFrameMAD();
      m_Nc = img->NumberofCodedPFrame;
    }
    else
    {
      /*compute the MAD of the current basic unit*/
      CurrentFrameMAD = img->TotalMADBasicUnit/img->BasicUnit;
      img->TotalMADBasicUnit = 0;
      /* compute the average number of header bits*/
      CodedBasicUnit = TotalNumberofBasicUnit-NumberofBasicUnit;
      if(CodedBasicUnit>0)
      {
        PAveHeaderBits1 = (int_32_t)(1.0*(PAveHeaderBits1*(CodedBasicUnit-1)+img->NumberofBasicUnitHeaderBits)/CodedBasicUnit+0.5);
        if(PAveHeaderBits3 == 0)
          PAveHeaderBits2 = PAveHeaderBits1;
        else
          PAveHeaderBits2 = (int_32_t)(1.0*(PAveHeaderBits1*CodedBasicUnit+PAveHeaderBits3*NumberofBasicUnit)/TotalNumberofBasicUnit+0.5);
      }
      /*update the record of MADs for reference*/
      if(((input->InterlaceCodingOption==PAFF_CODING)) && (img->FieldControl==1))
        FCBUCFMAD[TotalNumberofBasicUnit-1-NumberofBasicUnit] = CurrentFrameMAD;
      else
        BUCFMAD[TotalNumberofBasicUnit-1-NumberofBasicUnit] = CurrentFrameMAD;
      if(NumberofBasicUnit!=0)
        m_Nc = img->NumberofCodedPFrame*TotalNumberofBasicUnit+CodedBasicUnit;
      else
        m_Nc = (img->NumberofCodedPFrame-1)*TotalNumberofBasicUnit+CodedBasicUnit;
    }

    if(m_Nc>1)
      MADModelFlag = mytrue;
    PPreHeader = img->NumberofHeaderBits;
    for (i = 19; i > 0; i--)
    {// update the history
      PP_frm_Qstep[i] = PP_frm_Qstep[i - 1];
      P_frm_Qstep[i]  = PP_frm_Qstep[i];
      PP_frm_R[i]     = PP_frm_R[i - 1];
      P_frm_R[i]      = PP_frm_R[i];
    }
    PP_frm_Qstep[0] = QP2Qstep(m_Qc);
    if(img->BasicUnit==img->total_number_mb)
      PP_frm_R[0] = img->NumberofTextureBits*1.0/CurrentFrameMAD;
    else
      PP_frm_R[0] = img->NumberofBasicUnitTextureBits*1.0/CurrentFrameMAD;

    P_frm_Qstep[0] = PP_frm_Qstep[0];
    P_frm_R[0] = PP_frm_R[0];
    m_X1 = Pm_X1;
    m_X2 = Pm_X2;
    /*compute the size of window*/
    n_windowSize = (CurrentFrameMAD>PreviousFrameMAD)?(int_32_t)(PreviousFrameMAD/CurrentFrameMAD*20):(int_32_t)(CurrentFrameMAD/PreviousFrameMAD*20);
    n_windowSize = max(n_windowSize, 1);
    n_windowSize = min(n_windowSize,m_Nc);
    n_windowSize = min(n_windowSize,m_windowSize+1);
    n_windowSize = min(n_windowSize,20);
    /*update the previous window size*/
    m_windowSize = n_windowSize;

    memset(m_rgRejected, myfalse, 20*sizeof(myboolean));
    // initial RD model estimator
    RCModelEstimator (n_windowSize);
    for (i = 0; i < n_windowSize; i++)
    {
      error[i] = m_X1 / P_frm_Qstep[i] + m_X2 / (P_frm_Qstep[i] * P_frm_Qstep[i]) - P_frm_R[i];
      std += error[i] * error[i];
    }
    threshold = (n_windowSize == 2) ? 0 : sqrt (std / n_windowSize);
    for (i = 0; i < n_windowSize; i++)
    {
      if (fabs(error[i]) > threshold)
        m_rgRejected[i] = mytrue;
    }
    // always include the last data point
    m_rgRejected[0] = myfalse;
    // second RD model estimator
    RCModelEstimator (n_windowSize);
    if(MADModelFlag)
      updateMADModel();
    else if(img->type==INTER_IMG)
      PPictureMAD[0] = CurrentFrameMAD;
  }
}


void c_avs_enc::RCModelEstimator (int_32_t n_windowSize)
{
  int_32_t n_real_win_size = n_windowSize;
  int_32_t i;
  double oneSampleQ;
  double a00 = 0.0, a01 = 0.0, a10 = 0.0, a11 = 0.0, b0 = 0.0, b1 = 0.0;
  double MatrixValue;
  myboolean estimateX2 = myfalse;
  // find the number of samples which are not rejected
  for (i = 0; i < n_windowSize; i++)
  {
    if (m_rgRejected[i])
      n_real_win_size--;
  }
  // default RD model estimation results
  m_X1 = m_X2 = 0.0;
  for (i = 0; i < n_windowSize; i++)
  {
    if (!m_rgRejected[i])
      oneSampleQ = P_frm_Qstep[i];
  }
  // if all non-rejected Q are the same, take 1st order model
  for (i = 0; i < n_windowSize; i++)
  {
    if (!m_rgRejected[i])
    {
      if (P_frm_Qstep[i] != oneSampleQ)
      {
        estimateX2 = mytrue;
      }
      m_X1 += (P_frm_Qstep[i] * P_frm_R[i]) / n_real_win_size;
    }
  }

  // take 2nd order model to estimate X1 and X2
  if ((n_real_win_size >= 1) && estimateX2)
  {
    for (i = 0; i < n_windowSize; i++)
    {
      if (!m_rgRejected[i])
      {
        a00 = a00 + 1.0;
        a01 += 1.0 / P_frm_Qstep[i];
        a10 = a01;
        a11 += 1.0 / (P_frm_Qstep[i] * P_frm_Qstep[i]);
        b0 += P_frm_Qstep[i] * P_frm_R[i];
        b1 += P_frm_R[i];
      }
    }
    // solve the equation of AX = B
    MatrixValue=a00*a11-a01*a10;
    if(fabs(MatrixValue)>0.000001)
    {
      m_X1=(b0*a11-b1*a01)/MatrixValue;
      m_X2=(b1*a00-b0*a10)/MatrixValue;
    }
    else
    {
      m_X1=b0/a00;
      m_X2=0.0;
    }
  }
  if(img->type==INTER_IMG)
  {
    Pm_X1 = m_X1;
    Pm_X2 = m_X2;
  }
}

double c_avs_enc::ComputeFrameMAD()
{
  double TotalMAD;
  int_32_t i;
  TotalMAD = 0.0;
  for(i=0;i<img->total_number_mb;i++)
    TotalMAD += img->MADofMB[i];
  TotalMAD /= img->total_number_mb;
  return TotalMAD;
}


//update the parameters of linear prediction model
void c_avs_enc::updateMADModel ()
{
  int_32_t n_windowSize;
  int_32_t i;
  double error[20], std = 0.0, threshold;
  int_32_t m_Nc;

  if(img->NumberofCodedPFrame>0)
  {
    if(img->type==INTER_IMG)
    {
      /*frame layer rate control*/
      if(img->BasicUnit==img->total_number_mb)
        m_Nc=img->NumberofCodedPFrame;
      /*basic unit layer rate control*/
      else
        m_Nc=img->NumberofCodedPFrame*TotalNumberofBasicUnit+CodedBasicUnit;
      for (i = 19; i > 0; i--) {// update the history
        PPictureMAD[i] = PPictureMAD[i - 1];
        PictureMAD[i]=PPictureMAD[i];
        ReferenceMAD[i]= ReferenceMAD[i-1];
      }
      PPictureMAD[0] = CurrentFrameMAD;
      PictureMAD[0]=PPictureMAD[0];
      if(img->BasicUnit==img->total_number_mb)
        ReferenceMAD[0]=PictureMAD[1];
      else
      {
        if(((input->InterlaceCodingOption==2))&&(img->FieldControl==1))
          ReferenceMAD[0]=FCBUPFMAD[TotalNumberofBasicUnit-1-NumberofBasicUnit];
        else
          ReferenceMAD[0]=BUPFMAD[TotalNumberofBasicUnit-1-NumberofBasicUnit];
      }
      MADPictureC1=PMADPictureC1;
      MADPictureC2=PMADPictureC2;
    }
    /*compute the size of window*/
    n_windowSize = (CurrentFrameMAD>PreviousFrameMAD)?(int_32_t)(PreviousFrameMAD/CurrentFrameMAD*20)\
      :(int_32_t)(CurrentFrameMAD/PreviousFrameMAD*20);
    n_windowSize = min(n_windowSize,(m_Nc-1));
    n_windowSize = max(n_windowSize, 1);
    n_windowSize = min(n_windowSize,MADm_windowSize+1);
    n_windowSize = min(20,n_windowSize);
    /*update the previous window size*/
    MADm_windowSize=n_windowSize;

    memset(PictureRejected, myfalse, 20*sizeof(myboolean));
    //update the MAD for the previous frame
    if(img->type==INTER_IMG)
      PreviousFrameMAD=CurrentFrameMAD;

    // initial MAD model estimator
    MADModelEstimator (n_windowSize);
    for (i = 0; i < n_windowSize; i++) {
      error[i] = MADPictureC1*ReferenceMAD[i]+MADPictureC2-PictureMAD[i];
      std += error[i] * error[i];
    }
    threshold = (n_windowSize == 2) ? 0 : sqrt (std / n_windowSize);
    for (i = 0; i < n_windowSize; i++) {
      if (fabs(error[i]) > threshold)
        PictureRejected[i] = mytrue;
    }
    // always include the last data point
    PictureRejected[0] = myfalse;

    // second MAD model estimator
    MADModelEstimator (n_windowSize);
  }
}

void c_avs_enc::MADModelEstimator (int_32_t n_windowSize)
{
  int_32_t n_realSize = n_windowSize;
  int_32_t i;
  double oneSampleQ;
  double a00 = 0.0, a01 = 0.0, a10 = 0.0, a11 = 0.0, b0 = 0.0, b1 = 0.0;
  double MatrixValue;
  myboolean estimateX2 = myfalse;

  for (i = 0; i < n_windowSize; i++) {// find the number of samples which are not rejected
    if (PictureRejected[i])
      n_realSize--;
  }

  // default MAD model estimation results

  MADPictureC1 = MADPictureC2 = 0.0;

  for (i = 0; i < n_windowSize; i++)  {
    if (!PictureRejected[i])
      oneSampleQ = PictureMAD[i];
  }
  for (i = 0; i < n_windowSize; i++)  {// if all non-rejected MAD are the same, take 1st order model
    if ((PictureMAD[i] != oneSampleQ) && !PictureRejected[i])
      estimateX2 = mytrue;
    if (!PictureRejected[i])
      MADPictureC1 += PictureMAD[i] / (ReferenceMAD[i]*n_realSize);
  }

  // take 2nd order model to estimate X1 and X2
  if ((n_realSize >= 1) && estimateX2)
  {
    for (i = 0; i < n_windowSize; i++)
    {
      if (!PictureRejected[i]) {
        a00 = a00 + 1.0;
        a01 += ReferenceMAD[i];
        a10 = a01;
        a11 += ReferenceMAD[i]*ReferenceMAD[i];
        b0 += PictureMAD[i];
        b1 += PictureMAD[i]*ReferenceMAD[i];
      }
    }
    // solve the equation of AX = B
    MatrixValue=a00*a11-a01*a10;
    if(fabs(MatrixValue)>0.000001)
    {
      MADPictureC2=(b0*a11-b1*a01)/MatrixValue;
      MADPictureC1=(b1*a00-b0*a10)/MatrixValue;
    }
    else
    {
      MADPictureC1=b0/a01;
      MADPictureC2=0.0;
    }
  }
  if(img->type==INTER_IMG)
  {
    PMADPictureC1=MADPictureC1;
    PMADPictureC2=MADPictureC2;
  }
}


double c_avs_enc::QP2Qstep( int_32_t QP )
{
  int_32_t i;
  double Qstep;
  TLS static const double QP2QSTEP[8] = { 1.0, 1.0905, 1.189, 1.297, 1.414, 1.542, 1.682, 1.834};

  Qstep = QP2QSTEP[QP % 8];
  for( i=0; i<(QP/8); i++)
    Qstep *= 2;

  return Qstep;
}

int_32_t c_avs_enc::Qstep2QP( double Qstep )
{
  int_32_t q_per = 0, q_rem = 0;

  //  assert( Qstep >= QP2Qstep(0) && Qstep <= QP2Qstep(51) );
  if( Qstep < QP2Qstep(0))
    return 0;
  else if (Qstep > QP2Qstep(63) )
    return 63;

  while( Qstep > QP2Qstep(7) )
  {
    Qstep /= 2;
    q_per += 1;
  }

  if (Qstep <= (1.0+1.0905)/2)
  {
    Qstep = 1.0;
    q_rem = 0;
  }
  else if (Qstep <= (1.0905+1.1895)/2)
  {
    Qstep = 1.0905;
    q_rem = 1;
  }
  else if (Qstep <= (1.189+1.297)/2)
  {
    Qstep = 1.189;
    q_rem = 2;
  }
  else if (Qstep <= (1.297+1.414)/2)
  {
    Qstep = 1.297;
    q_rem = 3;
  }
  else if (Qstep <= (1.414+1.542)/2)
  {
    Qstep = 1.414;
    q_rem = 4;
  }
  else if (Qstep <= (1.542+1.682)/2)
  {
    Qstep = 1.414;
    q_rem = 5;
  }
  else if (Qstep <= (1.682+1.834)/2)
  {
    Qstep = 1.682;
    q_rem = 6;
  }
  else
  {
    Qstep = 1.834;
    q_rem = 7;
  }

  return (q_per * 8 + q_rem);
}
