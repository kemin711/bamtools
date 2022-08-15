// ***************************************************************************
// BamReader.cpp (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 29 July 2013 (DB)
// ---------------------------------------------------------------------------
// Provides read access to BAM files.
// ***************************************************************************

#include "api/BamReader.h"
#include "api/internal/bam/BamReader_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace std;

BamReader::BamReader() 
: d(new BamTools::Internal::BamReaderPrivate(this)) { }

BamReader::~BamReader(void) {
   delete d;
   d = 0;
}

bool BamReader::Close(void) {
   return d->Close();
}

const string BamReader::GetFilename(void) const {
   return d->Filename();
}
bool BamReader::IsOpen(void) const {
   return d->IsOpen();
}

bool BamReader::Jump(int refID, int position) {
   return d->SetRegion( BamRegion(refID, position) );
}
bool BamReader::Open(const std::string& filename) {
   return d->Open(filename);
}
bool BamReader::Rewind(void) {
   return d->Rewind();
}
bool BamReader::SetRegion(const BamRegion& region) { 
   return d->SetRegion(region);
}
bool BamReader::SetRegion(const int& leftRefID, const int& leftPos,
               const int& rightRefID, const int& rightPos) 
{
   return d->SetRegion(BamRegion(leftRefID, leftPos, rightRefID, rightPos));
}

bool BamReader::GetNextAlignment(BamAlignment& alignment) {
   // uses Internal::BamReaderPrivate to do the work
   return d->GetNextAlignment(alignment);
}

BamAlignment* BamReader::next() {
   BamAlignment* ptr = new BamAlignment();
   if (d->GetNextAlignment(*ptr)) {
      return ptr;
   }
   delete ptr;
   return nullptr;
}

bool BamReader::GetNextAlignmentCore(BamAlignment& alignment) {
   return d->GetNextAlignmentCore(alignment);
}
const SamHeader& BamReader::GetConstSamHeader(void) const {
   return d->GetConstSamHeader();
}
SamHeader BamReader::GetHeader(void) const {
   return d->GetSamHeader();
}
std::string BamReader::GetHeaderText(void) const {
   return d->GetHeaderText();
}
int BamReader::GetReferenceCount(void) const {
   return d->GetReferenceCount();
}
const RefVector& BamReader::GetReferenceData(void) const {
   return d->GetReferenceData();
}

int BamReader::GetReferenceID(const std::string& refName) const {
   return d->GetReferenceID(refName);
}

bool BamReader::CreateIndex(const BamIndex::IndexType& type) {
   return d->CreateIndex(type);
}
bool BamReader::HasIndex(void) const {
   return d->HasIndex();
}
bool BamReader::LocateIndex(const BamIndex::IndexType& preferredType)
{
   return d->LocateIndex(preferredType);
}
bool BamReader::OpenIndex(const std::string& indexFilename) {
   return d->OpenIndex(indexFilename);
}
void BamReader::SetIndex(BamIndex* index) {
   d->SetIndex(index);
}

//std::string BamReader::GetErrorString(void) const {
//   return d->GetErrorString();
//}
        
vector<pair<string,int>> BamReader::getReferenceMetaData() const {
   RefVector tmp = GetReferenceData();
   vector<pair<string, int> > res(tmp.size());
   for (size_t i=0; i<tmp.size(); ++i) {
      res[i]=make_pair(tmp[i].getRefname(), tmp[i].getReflength());
   }
   return res;
}

const string& BamReader::getReferenceName(int refid) const {
   const RefVector& rd = GetReferenceData();
   if (refid < 0 || refid >= static_cast<int>(rd.size())) {
      throw logic_error("reference id " + to_string(refid) + " is invalid");
   }
   return rd[refid].getRefname();
}
