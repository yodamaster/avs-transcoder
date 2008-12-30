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
 * Copyright (C) Cisco Systems Inc. 2004.  All Rights Reserved.
 * 
 * Contributor(s): 
 *		Bill May wmay@cisco.com
 */

#include "mp4common.h"

/*
 * SizeTableProperty is a special version of the MP4TableProperty - 
 * the BytesProperty will need to set the value before it can read
 * from the file
 */
class SizeTableProperty : public MP4TableProperty 
{
 public:
  SizeTableProperty(char *name, MP4IntegerProperty *pCountProperty) :
    MP4TableProperty(name, pCountProperty) {};
 protected:
  void ReadEntry(MP4File *pFile, u_int32_t index) {
    // Each table has a size, followed by the length field
    // first, read the length
    m_pProperties[0]->Read(pFile, index);
    MP4IntegerProperty *pIntProp = (MP4IntegerProperty *)m_pProperties[0];
    // set the size in the bytes property
    MP4BytesProperty *pBytesProp = (MP4BytesProperty *)m_pProperties[1];
    pBytesProp->SetValueSize(pIntProp->GetValue(index), index);
    // And read the bytes
    m_pProperties[1]->Read(pFile, index);
  };
};

MP4AvssAtom::MP4AvssAtom() 
	: MP4Atom("avss")
{
	MP4BitfieldProperty *pCount;
	MP4TableProperty *pTable;
	
	AddProperty( new MP4Integer8Property("configurationVersion")); /* 0 */
	
	AddProperty( new MP4Integer8Property("AvsProfileIndication")); /* 1 */
	
	AddProperty( new MP4Integer8Property("AvsLevelIndication")); /* 2 */
	
	AddProperty( new MP4BitfieldProperty("reserved", 6)); /* 3 *///111111
	AddProperty( new MP4BitfieldProperty("lengthSizeMinusOne", 2)); /* 4 */
	AddProperty( new MP4BitfieldProperty("reserved1", 4)); /* 5 *///1111
	pCount = new MP4BitfieldProperty("numOfSequenceHeader", 4);
	AddProperty(pCount); /* 6 */
	
	pTable = new SizeTableProperty("sequenceEntries", pCount);
	AddProperty(pTable); /* 7 */
	pTable->AddProperty(new MP4Integer16Property("SequenceHeaderLength"));
	pTable->AddProperty(new MP4BytesProperty("SequenceHeaderNALUnit"));
	
}

void MP4AvssAtom::Generate()
{
	MP4Atom::Generate();

	((MP4Integer8Property*)m_pProperties[0])->SetValue(1);

	m_pProperties[3]->SetReadOnly(false);
	((MP4BitfieldProperty*)m_pProperties[3])->SetValue(0x3f);
	m_pProperties[3]->SetReadOnly(true);

	m_pProperties[5]->SetReadOnly(false);
	((MP4BitfieldProperty*)m_pProperties[5])->SetValue(0xf);
	m_pProperties[5]->SetReadOnly(true);
}

