
#include <mp4av_common.h>
#include <mp4av_avs.h>




extern "C" MP4TrackId MP4AV_AVS_HintTrackCreate (MP4FileHandle mp4File,
												 MP4TrackId mediaTrackId)
{
	MP4TrackId hintTrackId =
		MP4AddHintTrack(mp4File, mediaTrackId);	//添加hint相关的原子mdia.minf.hmhd;mdia.minf.stbl.stsd.rtp
	//mdia.minf.stbl.stsd.rtp .tims.timeScale;tref.hint;udta.hnti.sdp
	//udta.hinf
	
	if (hintTrackId == MP4_INVALID_TRACK_ID) {					
		return MP4_INVALID_TRACK_ID;
	}
	
	u_int8_t payloadNumber = MP4_SET_DYNAMIC_PAYLOAD;
	
	// don't include mpeg4-esid
	MP4SetHintTrackRtpPayload(mp4File, hintTrackId, 
		"AVS1-P2", &payloadNumber, 0,						//改变名称“H264”
		NULL, true, false);
				
	// get the mpeg4 video configuration 
	u_int8_t **pSeq, **pPict ;
	u_int32_t *pSeqSize, *pPictSize;								//参数改变
	char *base64;
	uint32_t profile_level;
	char *sprop = NULL;
	uint32_t ix = 0;
	
	MP4GetTrackAVSSeqHeaders(mp4File,
		mediaTrackId,
		&pSeq,
		&pSeqSize);										//**获得序列头相关参数pSeq，pSeqSize
	
	if (pSeqSize && pSeqSize[0] != 0) {								//添加seq头
		// we have valid sequence and picture headers
		uint8_t *p = pSeq[0];
		if (*p == 0 && p[1] == 0 && p[2] == 1 ) {
			p += 3;
		}
		profile_level = p[1] << 8 | p[2];
		while (pSeqSize[ix] != 0) {
			base64 = MP4BinaryToBase64(pSeq[ix], pSeqSize[ix]);
			if (sprop == NULL) {
				sprop = strdup(base64);
			} else {
				sprop = (char *)realloc(sprop, strlen(sprop) + strlen(base64) + 1 + 1);
				strcat(sprop, ",");
				strcat(sprop, base64);
			}
			free(base64);
			free(pSeq[ix]);
			ix++;
		}
		free(pSeq);
		free(pSeqSize);
		
		
		// create the appropriate SDP attribute *
		char* sdpBuf = (char*)malloc(strlen(sprop) + 128);				//sdp
		
		u_int16_t svideoWidth = MP4GetTrackVideoWidth(mp4File, mediaTrackId);		
		u_int16_t svideoHeight = MP4GetTrackVideoHeight(mp4File, mediaTrackId);	//********添加********
		
		
		sprintf(sdpBuf,
			"a=cliprect:0,0,%d,%d\015\012"
			"a=fmtp:%u profile-level-id=%04x; sprop-parameter-sets=%s; packetization-mode=1\015\012",
			svideoWidth,
			svideoHeight,
			payloadNumber,
			profile_level,
			sprop); 
		
		/* add this to the track's sdp */
		MP4AppendHintTrackSdp(mp4File, hintTrackId, sdpBuf);
		
		free(sprop);
		free(sdpBuf);
	}
	return hintTrackId;
}


static uint32_t avs_get_nal_size (uint8_t *pData,
								  uint32_t sizeLength)
{
	if (sizeLength == 1) {
		return *pData;
	} else if (sizeLength == 2) {
		return (pData[0] << 8) | pData[1];
	} else if (sizeLength == 3) {
		return (pData[0] << 16) |(pData[1] << 8) | pData[2];
	}
	return (pData[0] << 24) |(pData[1] << 16) |(pData[2] << 8) | pData[3];
}



static uint8_t avs_get_sample_nal_type (uint8_t *pSampleBuffer, 
										uint32_t sampleSize,
										uint32_t sizeLength)
{
	uint8_t nal_type = pSampleBuffer[sizeLength] & 0x1f;
	return nal_type;
}




extern "C" void MP4AV_AVS_HintAddSample (MP4FileHandle mp4File,
										 MP4TrackId hintTrackId,
										 MP4SampleId sampleId,
										 uint8_t *pSampleBuffer,
										 uint32_t sampleSize,				//sampleSize：整个sample多少个byte
										 uint32_t sizeLength,
										 MP4Duration duration,
										 MP4Duration renderingOffset,
										 bool isSyncSample,
										 uint16_t maxPayloadSize)
{
	uint8_t nal_type = avs_get_sample_nal_type(pSampleBuffer, 
		sampleSize, 
		sizeLength);
	bool pic_is_idr = false;
	
	if (nal_type==AVS_NAL_TYPE_I_PIC_HEADER) 
	{
		pic_is_idr=true;
	}
	
	// for now, we don't know if we can drop frames, so don't indiate
	// that any are "b" frames
	bool isBFrame = false;
	uint32_t nal_size;									//nal_size：一个nalu多少个byte，不包括前面打包的4byte，包括nalu头
	uint32_t offset = 0;
	uint32_t remaining = sampleSize;	
	/*#ifdef DEBUG_H264_HINT			
	printf("hint for sample %d %u\n", sampleId, remaining);
#endif*/
	MP4AddRtpVideoHint(mp4File, hintTrackId, isBFrame, renderingOffset);
	
	if (sampleSize - sizeLength < maxPayloadSize) {				//sample小于MTU	sizelength相当于是头信息不包括在sample的净荷内				
		uint32_t first_nal = avs_get_nal_size(pSampleBuffer, sizeLength);
		if (first_nal + sizeLength == sampleSize) {					//一个sample只有一个nalu，且小于MTU时的情况
			// we have a single nal, less than the maxPayloadSize,	//只要打一个rtp包就可以了
			// so, we have Single Nal unit mode
			MP4AddRtpPacket(mp4File, hintTrackId, true);
			MP4AddRtpSampleData(mp4File, hintTrackId, sampleId,
				sizeLength, sampleSize - sizeLength);
			MP4WriteRtpHint(mp4File, hintTrackId, duration, 
				pic_is_idr);//nal_type改变规则
			return;
		}
	}
	
	// TBD should scan for resync markers (if enabled in ES config)
	// and packetize on those boundaries
	while (remaining) {				//sample数据大于MTU时，只要这个sample还有数据就执行********sample层**循环*******************
		
		nal_size = avs_get_nal_size(pSampleBuffer + offset,
			sizeLength);								//sizelength=4
		//******* skip the sizeLength								
		/*#ifdef DEBUG_H264_HINT//不执行
		printf("offset %u nal size %u remain %u\n", offset, nal_size, remaining);
#endif*/
		
		offset += sizeLength;
		remaining -= sizeLength;
		
		if (nal_size > maxPayloadSize) {				//如果nalu_size大于最大允许值MTU 需要分割 FU
#ifdef DEBUG_H264_HINT								//FU头8bit含义见说明
			printf("fragmentation units\n");
#endif
			uint8_t head = pSampleBuffer[offset];//************************************************************************
			offset++;
			nal_size--;
			remaining--;
			
			uint8_t fu_header[2];
			// uint8_t *final_fu_header;
			fu_header[0] = (head&0xe0)|0x1c;			//FU指示
			fu_header[1] = head;						//FU头
			fu_header[1] |= 0x80;
			fu_header[1] &= 0x9f;				//100+type5bit
			
			
			while (nal_size > 0) {			//******************只要nalu还有值************************nalu层**循环*******************
				
				
				uint32_t write_size;
				if (nal_size + 2 <= maxPayloadSize) {//如果nalu净荷加nalu头和FU头小于等于MTU则结束分割单元
					fu_header[1] |= 0x40;//分割单元结束
					write_size = nal_size;
				} else {
					write_size = maxPayloadSize - 2;//write_size就是这个分割单元包括的payload=MTU-2
				}
#ifdef DEBUG_H264_HINT
				printf("frag off %u write %u left in nal %u remain %u\n",
					offset, write_size, nal_size, remaining);
#endif
				remaining -= write_size;//remaining为分完一个MTU后这个sample剩余的数据
				
				MP4AddRtpPacket(mp4File, hintTrackId, remaining == 0);
				MP4AddRtpImmediateData(mp4File, hintTrackId, 
					fu_header, 2);
				fu_header[1] &= 0x7f;//表示结束第一个分割单元
				
				MP4AddRtpSampleData(mp4File, hintTrackId, sampleId, 
					offset, write_size);
				offset += write_size;
				nal_size -= write_size;//nal_size为这个nalu分割后剩余的数据，如果还有数据就在此循环，组成新的分割单元
			}
		}								//************FU****end********
		
		else {												//如果nalu_size 小于 最大允许值MTU，则要考虑是否打复合包stap
			// we have a smaller than MTU nal.  check the next sample			//这句话什么意思,为什么是next sample，应该是nalu吧
			// see if the next one fits;				
			uint32_t next_size_offset;
			bool have_stap = false;
			next_size_offset = offset + nal_size;
			if (next_size_offset < sampleSize) {//if (next_size_offset < remaining) {//是否有下一个nalu？有，则执行下面 remaining改成samplesize
				// we have a remaining NAL
				uint32_t next_nal_size = 
					avs_get_nal_size(pSampleBuffer + next_size_offset, sizeLength);
				if (next_nal_size + nal_size + 4 + 1 <= maxPayloadSize) {					//判断下一个两个nalu大小加起来是否小于MTU
					have_stap = true;															//小于MTU的话,我们就用STAP复合包
				}																		//4 两个NALU尺寸 1 一个STAP-A净载头 (NALU尺寸为1byte)
			} 
			if (have_stap == false) {									//如果两个nalu大于MTU则不用STAP,直接把这个nalu打成一个rtp包，下一个再重新判断
				MP4AddRtpPacket(mp4File, hintTrackId, next_size_offset >= remaining);
				MP4AddRtpSampleData(mp4File, hintTrackId, sampleId, 
					offset, nal_size);
				offset += nal_size;
				remaining -= nal_size;
			} 
			
			else {													//****************如果两个NALU小于MTU则使用STAP-A复合包实现******
				uint32_t bytes_in_stap = 1 + 2 + nal_size;					//1byte STAP-A净载头;2byte是STAP的nalu尺寸:表明随后nalu(包括nalu头)大小
				
				
				uint8_t max_nri = pSampleBuffer[offset] & 0x60;	//0 11 00000	//**重要问题:一个sample包含几个nalu**
				while (next_size_offset <= sampleSize && bytes_in_stap <= maxPayloadSize)//while (next_size_offset < sampleSize && bytes_in_stap <= maxPayloadSize) //remaining改成samplesize//while (next_size_offset < remaining && bytes_in_stap < maxPayloadSize) 
				{													
					uint8_t nri;										//在一个stap内比较出nalu最大的NRI优先级，在净载头中NRI单元就用最大优先级
					nri = pSampleBuffer[next_size_offset + sizeLength] & 0x60;		//0 11 00000
					if (nri > max_nri) max_nri = nri;
					
					uint32_t next_nal_size = 
						avs_get_nal_size(pSampleBuffer + next_size_offset, sizeLength);
					bytes_in_stap += 2 + next_nal_size;
					next_size_offset += sizeLength + next_nal_size;
				}
				
				//下面验证是否是最后一个nalu
				bool last;
				if (next_size_offset > sampleSize)// && bytes_in_stap <= maxPayloadSizey)//if (next_size_offset <= remaining && bytes_in_stap <= maxPayloadSize)		
					// stap is last frame						//已经到了这个sample最后一个nalu，M标志位置1，表示sample结束
				{
					last = true;
				} else last = false;						
				MP4AddRtpPacket(mp4File, hintTrackId, last);//last=true时z说明这个stap是这个sample最后一个stap，把M标志位置1
				uint8_t data[3];
				data[0] = max_nri | 24;						//stap-A 净载头
				data[1] = nal_size >> 8;
				data[2] = nal_size & 0xff;					//NALU尺寸2 byte
				MP4AddRtpImmediateData(mp4File, hintTrackId, 
					data, 3);
				MP4AddRtpSampleData(mp4File, hintTrackId, sampleId, 
					offset, nal_size);
				offset += nal_size;
				remaining -= nal_size;
				bytes_in_stap = 1 + 2 + nal_size;							//1byte STAP-A净载头;2byte是STAP的nalu尺寸:表明随后nalu(包括nalu头)大小
				nal_size = avs_get_nal_size(pSampleBuffer + offset, sizeLength);
				while (bytes_in_stap + nal_size + 2 <= maxPayloadSize &&remaining) {			//***当两个NALU以后还小于MTU则下面的nalu再打入stap**循环*
					offset += sizeLength;
					remaining -= sizeLength;
					data[0] = nal_size >> 8;
					data[1] = nal_size & 0xff;
					MP4AddRtpImmediateData(mp4File, hintTrackId, data, 2);
					MP4AddRtpSampleData(mp4File, hintTrackId, sampleId, offset, nal_size);
					offset += nal_size;
					remaining -= nal_size;
					bytes_in_stap += nal_size + 2;							//再次将nalu打入stap时就不用再加1byte净载头了，只需要加2byte nalu尺寸
					if (remaining) {
						nal_size = avs_get_nal_size(pSampleBuffer + offset, sizeLength);
					}
				} // end while stap
			} // end have stap	  
		} // end check size
   }
   MP4WriteRtpHint(mp4File, hintTrackId, duration,				//整个sample读完后写？？！
	   pic_is_idr);
   
}


extern "C" bool MP4AV_AVSHinter(
								MP4FileHandle mp4File, 
								MP4TrackId mediaTrackId, 
								u_int16_t maxPayloadSize)
{
	u_int32_t numSamples = MP4GetTrackNumberOfSamples(mp4File, mediaTrackId);
	u_int32_t maxSampleSize = MP4GetTrackMaxSampleSize(mp4File, mediaTrackId);
	
	uint32_t sizeLength;
	
	if (numSamples == 0 || maxSampleSize == 0) {
		return false;
	}
	
	/*if (MP4GetTrackAVSLengthSize(mp4File, mediaTrackId, &sizeLength) == false) {
    return false;
}*/
	sizeLength=4;						//why?
	
	MP4TrackId hintTrackId = 
		MP4AV_AVS_HintTrackCreate(mp4File, mediaTrackId);				//****AVSspecial****
	
	if (hintTrackId == MP4_INVALID_TRACK_ID) {
		return false;
	}
	
	u_int8_t* pSampleBuffer = (u_int8_t*)malloc(maxSampleSize);
	if (pSampleBuffer == NULL) {
		MP4DeleteTrack(mp4File, hintTrackId);
		return false;
	}
	for (MP4SampleId sampleId = 1; sampleId <= numSamples; sampleId++) {
		u_int32_t sampleSize = maxSampleSize;
		MP4Timestamp startTime;
		MP4Duration duration;
		MP4Duration renderingOffset;
		bool isSyncSample;//stss指定同步帧
		
		bool rc = MP4ReadSample(
			mp4File, mediaTrackId, sampleId, 
			&pSampleBuffer, &sampleSize, 
			&startTime, &duration, 
			&renderingOffset, &isSyncSample);
		
		if (!rc) {
			MP4DeleteTrack(mp4File, hintTrackId);
			CHECK_AND_FREE(pSampleBuffer);
			return false;
		}
		
		MP4AV_AVS_HintAddSample(mp4File,								//****AVSspecial****
			hintTrackId,
			sampleId,
			pSampleBuffer,
			sampleSize,
			sizeLength,
			duration,
			renderingOffset,
			isSyncSample,
			maxPayloadSize);
		
	}
	CHECK_AND_FREE(pSampleBuffer);
	
	return true;
}

