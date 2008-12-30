/*
 * The contents of this file are subject to the Mozilla Public
 * License Version 1.1 (the "License"); you may not use this file
 * except in compliance with the License. You may obtain a copy of
 * the License at http://www.mozilla.org/MPL/
 * 
 * Software distributed under the License is distributed on an "AS
 * IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * rights and limitations under the License.
 * 
 * The Original Code is MPEG4IP.
 * 
 * The Initial Developer of the Original Code is Cisco Systems Inc.
 * Portions created by Cisco Systems Inc. are
 * Copyright (C) Cisco Systems Inc. 2000-2002.  All Rights Reserved.
 * 
 * Contributor(s): 
 *		Dave Mackie		dmackie@cisco.com
 */

/* 
 * Notes:
 *  - file formatted with tabstops == 4 spaces 
 */

#include <mp4av_common.h>

bool MP4AV_RfcIsmaConcatenator(
	MP4FileHandle mp4File, 
	MP4TrackId mediaTrackId, 
	MP4TrackId hintTrackId,
	u_int8_t samplesThisHint, 
	MP4SampleId* pSampleIds, 
	MP4Duration hintDuration,
	u_int16_t maxPayloadSize)
{
	// handle degenerate case
	if (samplesThisHint == 0) {
		return true;
	}

	u_int8_t auPayloadHdrSize;

	// LATER would be more efficient if this were a parameter
	u_int8_t mpeg4AudioType =
		MP4GetTrackAudioMpeg4Type(mp4File, mediaTrackId);

	if (mpeg4AudioType == MP4_MPEG4_CELP_AUDIO_TYPE) {
		auPayloadHdrSize = 1;
	} else {
		auPayloadHdrSize = 2;
	}

	// construct the new hint
	MP4AddRtpHint(mp4File, hintTrackId);
	MP4AddRtpPacket(mp4File, hintTrackId, true);

	u_int8_t payloadHeader[2];

	u_int16_t numHdrBits = samplesThisHint * auPayloadHdrSize * 8;
	payloadHeader[0] = numHdrBits >> 8;
	payloadHeader[1] = numHdrBits & 0xFF;

	MP4AddRtpImmediateData(mp4File, hintTrackId,
		(u_int8_t*)&payloadHeader, sizeof(payloadHeader));

	u_int8_t i;

	// first the headers
	for (i = 0; i < samplesThisHint; i++) {
		MP4SampleId sampleId = pSampleIds[i];

		u_int32_t sampleSize = 
			MP4GetSampleSize(mp4File, mediaTrackId, sampleId);

		if (auPayloadHdrSize == 1) {
			// AU payload header is 6 bits of size
			// follow by 2 bits of the difference between sampleId's - 1
			payloadHeader[0] = sampleSize << 2;

		} else { // auPayloadHdrSize == 2
			// AU payload header is 13 bits of size
			// follow by 3 bits of the difference between sampleId's - 1
			payloadHeader[0] = sampleSize >> 5;
			payloadHeader[1] = (sampleSize & 0x1F) << 3;
		}

		if (i > 0) {
			payloadHeader[auPayloadHdrSize - 1] 
				|= ((sampleId - pSampleIds[i-1]) - 1); 
		}
#if 0
		printf("sample %u size %u %02x %02x prev sample %d\n", 
		       sampleId, sampleSize, payloadHeader[0],
		       payloadHeader[1], pSampleIds[i-1]);
#endif

		MP4AddRtpImmediateData(mp4File, hintTrackId,
			(u_int8_t*)&payloadHeader, auPayloadHdrSize);
	}

	// then the samples
	for (i = 0; i < samplesThisHint; i++) {
		MP4SampleId sampleId = pSampleIds[i];

		u_int32_t sampleSize = 
			MP4GetSampleSize(mp4File, mediaTrackId, sampleId);

		MP4AddRtpSampleData(mp4File, hintTrackId, sampleId, 0, sampleSize);
	}

	// write the hint
	MP4WriteRtpHint(mp4File, hintTrackId, hintDuration);

	return true;
}

bool MP4AV_RfcIsmaFragmenter(
	MP4FileHandle mp4File, 
	MP4TrackId mediaTrackId, 
	MP4TrackId hintTrackId,
	MP4SampleId sampleId, 
	u_int32_t sampleSize, 
	MP4Duration sampleDuration,
	u_int16_t maxPayloadSize)
{
	MP4AddRtpHint(mp4File, hintTrackId);

	MP4AddRtpPacket(mp4File, hintTrackId, false);

	// Note: CELP is never fragmented
	// so we assume the two byte AAC-hbr payload header
	u_int8_t payloadHeader[4];
	payloadHeader[0] = 0;
	payloadHeader[1] = 16;
	payloadHeader[2] = sampleSize >> 5;
	payloadHeader[3] = (sampleSize & 0x1F) << 3;

	MP4AddRtpImmediateData(mp4File, hintTrackId,
		(u_int8_t*)&payloadHeader, sizeof(payloadHeader));

	u_int16_t sampleOffset = 0;
	u_int16_t fragLength = maxPayloadSize - 4;

	do {
		MP4AddRtpSampleData(mp4File, hintTrackId,
			sampleId, sampleOffset, fragLength);

		sampleOffset += fragLength;

		if (sampleSize - sampleOffset > maxPayloadSize) {
			fragLength = maxPayloadSize; 
			MP4AddRtpPacket(mp4File, hintTrackId, false);
		} else {
			fragLength = sampleSize - sampleOffset; 
			if (fragLength) {
				MP4AddRtpPacket(mp4File, hintTrackId, true);
			}
		}
	} while (sampleOffset < sampleSize);

	MP4WriteRtpHint(mp4File, hintTrackId, sampleDuration);

	return true;
}

extern "C" bool MP4AV_RfcIsmaHinter(
	MP4FileHandle mp4File, 
	MP4TrackId mediaTrackId, 
	bool interleave,
	u_int16_t maxPayloadSize)
{
	// gather information, and check for validity
	return true;
}

