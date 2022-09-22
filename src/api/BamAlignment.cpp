// ***************************************************************************
// BamAlignment.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 4 December 2012 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#include "api/BamAlignment.h"
#include "api/BamConstants.h"
#include <functional>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <sstream>

//#define DEBUG

using namespace BamTools;
using namespace std;

///////// static member and functions /////////////////////
int BamAlignment::TRIMLEN_MAX = 6;
int BamAlignment::GAP_CUT = 3;
vector<pair<string,int>> BamAlignment::rsname;
map<string,int> BamAlignment::refname2id;

void BamAlignment::setRefvector(vector<pair<string,int>>&& refvec) {
   rsname=refvec;
   for (size_t i=0; i<refvec.size(); ++i) {
      refname2id[refvec[i].first] = i;
   }
}

vector<pair<char,int>> BamAlignment::parseCigar(const string& cigarstr) {
   string::size_type i=0, b;
   int len;
   char CO;
   vector<pair<char,int>> res;
   while (i < cigarstr.size()) {
      b=i;
      while (i < cigarstr.size()-1 && isdigit(cigarstr[i])) { ++i; }
      len = stoi(cigarstr.substr(b,i-b));
      CO = cigarstr[i];
      if ( CO == 'M' || CO == 'I' || CO == 'D' || CO == 'S' || CO == 'H' || CO == 'N') {
         res.push_back(make_pair(CO, len));
      }
      else {
         throw runtime_error(string(__FILE__) + to_string(__LINE__) + ":ERROR Illigal cigar op: " + CO);
      }
      ++i;
   }
   return res;
}

string BamAlignment::cigarToString(const vector<pair<char,int>>& cg) {
   ostringstream osr;
   for (auto& p : cg) {
      osr << p.second << p.first;
   }
   return osr.str();
}


/// Regular member function part /////

/*! \fn BamAlignment::BamAlignment(const BamAlignment& other)
    \brief copy constructor
*/
BamAlignment::BamAlignment(const BamAlignment& other)
    : Name(other.Name)
    , QueryBases(other.QueryBases)
    , AlignedBases(other.AlignedBases)
    , Qualities(other.Qualities)
    , TagData(other.TagData)
    , RefID(other.RefID)
    , Position(other.Position)
    , Bin(other.Bin)
    , MapQuality(other.MapQuality)
    , AlignmentFlag(other.AlignmentFlag)
    , CigarData(other.CigarData)
    , MateRefID(other.MateRefID)
    , MatePosition(other.MatePosition)
    , InsertSize(other.InsertSize)
    , SupportData(other.SupportData)
{ }

BamAlignment::BamAlignment(BamAlignment&& other)
    : Name(std::move(other.Name)), 
      QueryBases(std::move(other.QueryBases)),
      AlignedBases(std::move(other.AlignedBases)),
      Qualities(std::move(other.Qualities)),
      TagData(std::move(other.TagData)), 
      RefID(other.RefID), Position(other.Position), Bin(other.Bin), 
      MapQuality(other.MapQuality), AlignmentFlag(other.AlignmentFlag), 
      CigarData(std::move(other.CigarData)), 
      MateRefID(other.MateRefID), MatePosition(other.MatePosition), 
      InsertSize(other.InsertSize),  
      SupportData(std::move(other.SupportData))
{ }

BamAlignment& BamAlignment::operator=(const BamAlignment& other) {
   if (this != &other) {
      Name=other.Name;
      QueryBases=other.QueryBases;
      AlignedBases=other.AlignedBases;
      Qualities=other.Qualities;
      TagData=other.TagData;
      RefID=other.RefID;
      Position=other.Position;
      Bin=other.Bin;
      MapQuality=other.MapQuality;
      AlignmentFlag=other.AlignmentFlag;
      CigarData=other.CigarData;
      MateRefID=other.MateRefID;
      MatePosition=other.MatePosition;
      InsertSize=other.InsertSize;
      SupportData=other.SupportData;
   }
   return *this;
}

BamAlignment& BamAlignment::operator=(BamAlignment&& other) {
   if (this != &other) {
      Name=std::move(other.Name);
      QueryBases=std::move(other.QueryBases);
      AlignedBases=std::move(other.AlignedBases);
      Qualities=std::move(other.Qualities);
      TagData=std::move(other.TagData);
      RefID=other.RefID;
      Position=other.Position;
      Bin=other.Bin;
      MapQuality=other.MapQuality;
      AlignmentFlag=other.AlignmentFlag;
      CigarData=std::move(other.CigarData);
      MateRefID=other.MateRefID;
      MatePosition=other.MatePosition;
      InsertSize=other.InsertSize;
      SupportData=std::move(other.SupportData);
   }
   return *this;
}

//std::ostream& BamTools::operator<<(std::ostream &ous, const BamTools::BamAlignment &ba) {
// error: ‘std::ostream& BamTools::operator<<(std::ostream&, const BamTools::BamAlignment&)’ should have been declared inside ‘BamTools’
//  std::ostream& BamTools::operator<<(std::ostream &ous, const BamTools::BamAlignment &ba) {
//
/**
 * For human to read
 */
namespace BamTools {
std::ostream& operator<<(std::ostream &ous, const BamAlignment &ba) {
   const string sep="\t";
   ous << ba.getQueryName() << sep << ba.getQueryLength() << sep
      << ba.getPosition() << sep
      << "primary: " << ba.IsPrimaryAlignment() << sep
      //<< "reverseStrand: " << ba.IsReverseStrand() << sep
      << "strand: ";
   if (ba.IsReverseStrand()) ous << '-'; 
   else ous << '+';
   ous << sep;
      //<< "mateReverseStrand: " << ba.IsMateReverseStrand() << sep
   if (ba.IsPaired()) {
      if (ba.IsMateReverseStrand()) ous << '-';
      else ous << '+';
      ous << sep;
   }
   ous << "duplicate: " << ba.IsDuplicate() << sep
      << "mate: ";
   if (ba.IsFirstMate()) ous << 1;
   else ous << 2;
   ous << sep
      //<< "secondMate: " << ba.IsSecondMate() << sep
      //<< "mapped: " << ba.IsMapped() << sep 
      //<< "mateMapped: " << ba.IsMateMapped() << sep
      << "paired: " << ba.IsPaired() << sep
      << "properPair: " << ba.IsProperPair() << sep
      << "passedQC: " << !ba.IsFailedQC() << sep
      << "refid: " << ba.getReferenceId() << sep;
   if (ba.IsPaired()) {
      ous << "mateRefid: " << ba.getMateReferenceId() << sep
         << ba.getMatePosition() << sep;
   }
   ous << "insertSize: " << ba.getInsertSize() << sep
      << "mapQuality: " << ba.getMapQuality() << sep
      << ba.getAlignedQueryBases() << sep;
   for (size_t i=0; i<ba.CigarData.size(); ++i) {
      ous << ba.CigarData[i];
   }
   ous << sep;
   ous << ba.getQueryBases() << sep;
   vector<int> qs = ba.getQualityScore();
   ous << sep;
   copy(qs.begin(), qs.end(), ostream_iterator<int>(ous, "|"));
   ous << sep;
   /* problem regardless of which int type to use
    * int32_t or int, cause strange character to 
    * be printed to the terminal
    */
   vector<std::string> tagNames = ba.GetTagNames();
   // automatic probe is causing some problems
   // removing it.
   /*
   for (auto& t : tagNames) {
      char tagtype;
      ba.GetTagType(t, tagtype);
      if (tagtype == 'i' || tagtype == 'c' || tagtype == 's') { // BamConstants.h define the tag type
         int32_t intval; // must use this type, int is wrong
         ba.GetTag(t, intval);
         ous << t << ":" << tagtype << ":" << intval << "; ";
      }
      else if (tagtype == 'I' || tagtype == 'C' || tagtype == 'S') { // BamConstants.h define the tag type
         uint32_t intval; // must use this type, int is wrong
         ba.GetTag(t, intval);
         ous << t << ":" << tagtype << ":" << intval << "; ";
      }
      else if (tagtype == 'Z') {
         string tagval;
         ba.GetTag(t, tagval);
         ous << t << ": " << tagval << "; ";
      }
      else {
         ous << "tag type: " << tagtype << " for tag: " << t
            << " not considered yet\n";
      }
   }
   */
   // raw output
   for (const char c : ba.TagData) {
      if (c == '\0') ous << "|";
      else if (isprint(c)) 
         ous << c;
      else 
         ous << '~';
   }
   // print integer types
   for (auto& t : tagNames) {
      char tagtype;
      ba.GetTagType(t, tagtype);
      if (tagtype == 'i' || tagtype == 'c' || tagtype == 's') { // int 32, 16, 8 respectively BamConstants.h define the tag type
         //int32_t intval; // must use this type, int is wrong
         //ba.GetTag(t, intval);
         auto [ival, hastag] = ba.getTag<int32_t>(t);
         if (!hastag) {
            throw runtime_error("tag: " + t + " getTag() failed");
         }
         ous << t << ":" << tagtype << ":" << ival << "; ";
      }
      else if (tagtype == 'I' || tagtype == 'C' || tagtype == 'S') { // BamConstants.h define the tag type
        // uint32_t intval; // must use this type, int is wrong
        // ba.GetTag(t, intval);
         auto [uival, hastag] = ba.getTag<uint32_t>(t);
         if (!hastag) {
            throw runtime_error("Failed getTag call on " + t);
         }
         ous << t << ":" << tagtype << ":" << uival << "; ";
      }
      else if (t == "XM" || t == "XW" || t == "XD" || t == "YD" || t == "ZD") {
         vector<int32_t> tmpv = ba.getArrayTag<int32_t>(t);
         assert(tmpv.size() > 1);
         ous << t << ":" << tmpv[0] << ',' << tmpv[1];
         for (size_t x=2; x < tmpv.size(); x += 2) {
            ous << '|' << tmpv[x] << ',' << tmpv[x+1];
         }
         ous << " ";
      }
   }
   ous << endl;

   return ous;
}}

bool BamAlignment::operator<(const BamAlignment& other) const {
   if (getPosition() < other.getPosition()) return true;
   if (getPosition() > other.getPosition()) return false;
   if (getEndPosition() < other.getEndPosition()) return true;
   if (getEndPosition() > other.getEndPosition()) return false;
   return getMate() < other.getMate();
}

bool BamAlignment::operator>(const BamAlignment& other) const {
   if (getPosition() > other.getPosition()) return true;
   if (getPosition() < other.getPosition()) return false;
   if (getEndPosition() > other.getEndPosition()) return true;
   if (getEndPosition() < other.getEndPosition()) return false;
   return getMate() > other.getMate();
}

void BamAlignment::changePosition(int32_t alnstart) { 
   if (alnstart == getPosition()) return;
   if (CigarData.front().getType() == 'S') {
      CigarData.front().setLength(alnstart);
      CigarData[1].setLength(CigarData[1].getLength() - alnstart + Position);
   }
   else {
      CigarData[0].setLength(CigarData[0].getLength() - alnstart + Position);
   }
   Position = alnstart; 
   AlignedBases.clear();
}

bool BamAlignment::operator==(const BamAlignment& other) const {
   return getQueryName() == other.getQueryName()
      && getMate() == other.getMate();
}

vector<pair<char,int> > BamAlignment::getCigarOperation() const {
   vector<pair<char,int> > tmp(CigarData.size());
   transform(CigarData.begin(), CigarData.end(), tmp.begin(), mem_fn(&CigarOp::topair));
   return tmp;
}

string BamAlignment::getCigarString() const {
   string tmp;
   for (auto& cd : CigarData) {
      tmp += (to_string(cd.Length) + cd.Type);
   }
   return tmp;
}

bool BamAlignment::hasDCigar() const {
   for (auto& co : CigarData) {
      if (co.getType() == 'D') return true;
   }
   return false;
}

bool BamAlignment::lackDCigar() const {
   for (auto& co : CigarData) {
      if (co.getType() == 'D') return false;
   }
   return true;
}

bool BamAlignment::lackICigar() const {
   for (auto& co : CigarData) {
      if (co.getType() == 'I') return false;
   }
   return true;
}

bool BamAlignment::hasICigar() const {
   for (auto& co : CigarData) {
      if (co.getType() == 'I') return true;
   }
   return false;
}

bool BamAlignment::valid() const {
   //cerr << QueryBases.size() << " " << Qualities.size()
   //   << " " << getReferenceWidth() << " " << getReferenceWidth() <<  " ";
   if (QueryBases.size() != Qualities.size()) {
      cerr << __LINE__ << ": query sequence and quality not the same length\n";
      return false;
   }
   //return validCigar() && refwidthAgreeWithMD() && QueryBases.size() == Qualities.size();
   return validCigar() && refwidthAgreeWithMD();
}

bool BamAlignment::validCigar() const {
   if (isUnmapped() || lackCigar())
      return true;
   int cigarQL=0;
   int cigarRL=0;
   for (auto& c : CigarData) {
      if (c.getType() == 'S' || c.getType() == 'M' || c.getType() == 'I') {
         cigarQL += c.getLength();
      }
      if (c.getType() == 'M' || c.getType() == 'D') {
         cigarRL += c.getLength();
      }
   }
   //cerr << cigarQL << " " << cigarRL << endl;
   if (cigarQL != getLength()) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR Query Length inconsistent queryLength from cigar=" 
         << cigarQL << " queryLength=" << getLength() << endl;
      cerr << *this << endl;
   }
   if (cigarRL != getReferenceWidth()) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR reference length contradict cigar computed: "
         << cigarRL << endl;
      cerr << *this << endl;
   }
   return cigarQL == getLength() && cigarRL == getReferenceWidth();
}

bool BamAlignment::sameCigar(const vector<pair<char,int>>& cigar) const {
   if (cigar.size() != CigarData.size()) 
      return false;
   for (size_t i=0; i<cigar.size(); ++i) {
      if (cigar[i].first != CigarData[i].getType() 
            || cigar[i].second != (int)CigarData[i].getLength()) 
         return false;
   }
   return true;
}

bool BamAlignment::hasEndIndel() const {
   if (CigarData.size() < 3 || CigarData.front().getType() == 'S' 
         || CigarData.back().getType() == 'S') 
   {
      return false;
   }
   int i = 0;
   while (i < (int)CigarData.size() && CigarData[i].getType() != 'M') ++i;
   if (i < (int)CigarData.size()-1 && CigarData[i].getLength() < 22U && (CigarData[i+1].getType() == 'I' || CigarData[i+1].getType() == 'D'))
      return true;
   i = CigarData.size()-1;
   while (i > 0 && CigarData[i].getType() != 'M') --i;
   if (i > 0 && CigarData[i].getLength() < 22U && (CigarData[i-1].getType() == 'I' || CigarData[i-1].getType() == 'D'))
      return true;
   return false;
}

bool BamAlignment::hasAmbiguousBase() const {
   for (char b : getQuerySequence()) {
      b = toupper(b);
      if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
         return true;
      }
   }
   return false;
}

int BamAlignment::getMatchedReferenceLength() const {
   int sum = 0;
   for (auto& c : CigarData) {
      if (c.getType() == 'M') sum += c.getLength();
   }
   return sum;
}

string BamAlignment::getFirstSoftclip() const {
   if (!startWithSoftclip()) return "";
   return getQuerySequence().substr(0, getCigar().front().Length);
}

string BamAlignment::getFirstSoftquality() const {
   if (!startWithSoftclip()) return "";
   return getQuality().substr(0, getCigar().front().getLength());
}

int BamAlignment::getFirstSoftclipLength() const {
   if (!startWithSoftclip()) return 0;
   return getCigar().front().Length;
}

string BamAlignment::getLastSoftclip() const {
   if (!endWithSoftclip()) return "";
   return getQueryBases().substr(getQueryLength() - getCigar().back().Length);
}
int BamAlignment::getLastSoftclipLength() const {
   if (!endWithSoftclip()) return 0;
   return getCigar().back().Length;
}
string BamAlignment::getLastSoftquality() const {
   if (!endWithSoftclip()) return "";
   return getQuality().substr(getQueryLength() - getCigar().back().getLength());
}

int BamAlignment::getSoftclipLength() const {
   if (CigarData.empty()) return 0;
   int res = 0;
   if (getCigar().front().getType() == 'S')
      res += getCigar().front().getLength();
   if (getCigar().back().getType() == 'S')
      res += getCigar().back().getLength();
   return res;
}

int BamAlignment::getMaxSoftclipLength() const {
   if (CigarData.empty()) return 0;
   int res = 0;
   if (getCigar().front().getType() == 'S')
      res = getCigar().front().getLength();
   if (getCigar().back().getType() == 'S' &&
         static_cast<int>(getCigar().back().getLength()) > res)
      res = getCigar().back().getLength();
   return res;
}

// 221M4I2M1D38M
void BamAlignment::setCigar(const string& cstr) {
   if (!CigarData.empty()) CigarData.clear();
   unsigned int i=0;
   unsigned int b;
   while (i < cstr.size()) {
      b=i;
      while (isdigit(cstr[i])) ++i;
      CigarData.push_back(CigarOp(cstr[i], stoi(cstr.substr(b,i-b))));
      ++i;
   }
   SupportData.NumCigarOperations=CigarData.size();
}

void BamAlignment::setCigarOperation(const std::vector<pair<char,int> > &cd) {
   CigarData.resize(cd.size());
   for (size_t i=0; i < CigarData.size(); ++i) {
      CigarData[i].fromPair(cd[i]);
   }
   SupportData.NumCigarOperations=cd.size();
}

// only do one correction!
// deprecated, use fix1M to fix 1 or 2 M flanked by indel
// in all 4 possible combinations
void BamAlignment::fixStaggerGap() {
   if (CigarData.size() < 5) return;
   for (size_t i=0; i+4 < CigarData.size(); ++i) {
      if ((CigarData[i].Type  == 'M' && 
          CigarData[i+1].Type == 'D' && CigarData[i+1].Length == CigarData[i+3].Length &&
          CigarData[i+2].Type == 'M' && CigarData[i+2].Length == 1 &&
          CigarData[i+3].Type == 'I' &&
          CigarData[i+4].Type == 'M' ) ||
         (CigarData[i].Type   == 'M' && 
          CigarData[i+1].Type == 'I' && CigarData[i+1].Length == CigarData[i+3].Length &&
          CigarData[i+2].Type == 'M' && CigarData[i+2].Length == 1 &&
          CigarData[i+3].Type == 'D' &&
          CigarData[i+4].Type == 'M' ))
      {
         int gaplen = CigarData[i+1].Length;
         CigarData[i].Length += (1 + gaplen + CigarData[i+4].Length);
         CigarData.erase(CigarData.begin()+i+1, CigarData.begin()+i+5);
         // here is a devil in duplicating data: we have to remember to change 
         // the duplicated information, design flaw
         SupportData.NumCigarOperations = CigarData.size();
         uint16_t edit=0;
         try {
            pair<uint16_t,bool> res = getTag<uint16_t>("NM");
            if (res.second) edit = res.first;
         }
         catch (const BamTypeException& ler) {
            pair<int,bool> res = getTag<int>("NM");
            if (res.second) edit = res.first;
         }
         if (gaplen > edit) {
            throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR gap lenth greater than edit");
         }
         edit -= gaplen;
         EditTag("NM", "C", edit);
         return;
      }
   }
}

// make this also fix stagger M left and right different
bool BamAlignment::fix1M() {
   if (CigarData.size() < 5) return false;
   bool changed=false;
   for (int i=2; i+2 < getCigarSize(); ++i) {
      if (CigarData[i].getType() == 'M' && CigarData[i].getLength() < 3) 
      { // candidate for doing fixing, imply i-2>-1 && i+2 < CigarData.size()
         //cerr << " i=" << i << " id fix1M candidate\n";
         if (CigarData[i-1].getType() == CigarData[i+1].getType()) { // same ins or del
            if (CigarData[i-1].getType() == 'I' || CigarData[i-1].getType() == 'D') {
               //int nmvalue = getNMValue();
               uint32_t nmvalue = getNMValue();
               if (getCigarType(i-2) == 'M') { // merge 1M with previous M
                  // the last segment could be S  MD1MS is fine, last done no need M
                  nmvalue += getCigarLength(i);
                  CigarData[i-2].expand(getCigarLength(i));
                  CigarData[i-1].expand(getCigarLength(i+1));
                  CigarData.erase(CigarData.begin()+i, CigarData.begin()+i+2);
                  changed=true;
               }
               else if (CigarData[i+2].getType() == 'M') {
                  nmvalue += getCigarLength(i);
                  CigarData[i+2].expand(getCigarLength(i));
                  CigarData[i+1].expand(CigarData[i-1].getLength());
                  CigarData.erase(CigarData.begin()+i-1, CigarData.begin()+i+1);
                  changed=true;
               }
               else {
                  cerr << __FILE__ << ":" << __LINE__ << ":WARN Cannot fix 1M somthing is wrong!\n";
                  throw runtime_error("failed to fix1M");
               }
               //EditTag("NM", "i", nmvalue);
               if (nmvalue < UINT8_MAX) 
                  EditTag("NM", "C", nmvalue);
               else {
                  throw runtime_error("nmvalue is greater than UINT8_MAX");
               }
               SupportData.NumCigarOperations = CigarData.size();
            }
            //cerr << "fixed 1M\n";
         }
         else if (((getCigarType(i-1) == 'I' && getCigarType(i+1) == 'D') ||
                     (getCigarType(i-1) == 'D' && getCigarType(i+1) == 'I')) &&
               getCigarType(i-2) == 'M' && getCigarType(i+2) == 'M')
         { // M[DI]M[ID]M
            //cerr << "old cigar: ";
            //for (auto& c : CigarData) {
            //   cerr << c.getLength() << c.getType();
            //}
            //cerr << endl;
            uint32_t nmvalue = getNMValue();
            if (getCigarLength(i-1) < getCigarLength(i+1)) {
               CigarData[i-2].expand(getCigarLength(i-1) + getCigarLength(i));
               CigarData[i+1].setLength(getCigarLength(i+1) - getCigarLength(i-1));
               nmvalue += (getCigarLength(i) - getCigarLength(i-1));
               CigarData.erase(CigarData.begin()+i-1, CigarData.begin()+i+1);
               // erase [DI]M 2 segments to the left of i
               SupportData.NumCigarOperations -= 2;
            }
            else if (getCigarLength(i-1) > getCigarLength(i+1)) {
               CigarData[i+2].expand(getCigarLength(i+1) + getCigarLength(i));
               CigarData[i-1].setLength(getCigarLength(i-1) - getCigarLength(i+1));
               nmvalue += (getCigarLength(i) - getCigarLength(i+1));
               CigarData.erase(CigarData.begin()+i, CigarData.begin()+i+2);
               // erase M[ID] 2 segments to the right of i
               SupportData.NumCigarOperations -= 2;
            }
            else { // lengths of segment i-1 and i+1 are identical
               CigarData[i-2].expand(getCigarLength(i) + getCigarLength(i+1) + getCigarLength(i+2));
               nmvalue += getCigarLength(i) - getCigarLength(i+1);
               SupportData.NumCigarOperations -= 4;
               if ((int)CigarData.size() == i+3) {
                  CigarData.resize(getCigarSize() - 4);
               }
               else if ((int)CigarData.size() > i+3) {
                  CigarData.erase(CigarData.begin()+i-1, CigarData.begin()+i+3);
               }
               else {
                  throw logic_error("cigar segment too few");
               }
            }
            if (nmvalue < UINT8_MAX) {
               EditTag("NM", "i", nmvalue); // approximate, exact computation too expensive
            }
            else {
               throw runtime_error("nmvalue greater than UMI8_MAX");
            }
            changed=true;
            //cerr << "fixed stagger: ";
            //for (auto& c : CigarData) {
            //   cerr << c.getLength() << c.getType();
            //}
            //cerr << endl;
         }
      }
   }
   if (changed && !validCigar()) { // only check if changed
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR CigarData and len mismatch\n" << *this << endl;
      throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR CigarData and query sequence length does not match");
   }
   return changed;
}

/* must be inline for performance
void BamAlignment::fixCigarError() {
   if (getCigarType(0) == 'I') {
      CigarData[0].setType('S');
   }
}
*/

void BamAlignment::moveDeletion(int oldloc, int newloc) {
   assert(oldloc != newloc);
   pair<int,bool> rtn = isDeletionAtRefloc(newloc);
   if (rtn.second) {
      //cerr << __FILE__ << ":" << __LINE__ << ":WARN there is already a deletion at "
      //   << newloc << " skip deletion reloacation from " << oldloc << endl
      //   << *ba << endl;
      return;
   }
   //cerr << __LINE__ << ": no deletion at " << newloc << " c=" << rtn.first << endl;
   rtn = isInsertionAtRefloc(newloc, 0);
   if (rtn.second) {
      //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG cannot move deletion to an insertion place\n";
      //throw runtime_error(string(__func__) + ":DEBUG cannot move deletion to location with insertion");
      return;
   }
   int c=0, r=0, q=0; 
   while (c < getCigarSize() && (getCigarType(c) == 'S' || getCigarType(c) == 'H')) { 
      q += getCigarLength(c);
      ++c; 
   }
   while (c < getCigarSize() && r < oldloc && q < getLength()) {
      if (getCigarType(c) == 'M') {
         r += getCigarLength(c);
         q += getCigarLength(c);
      }    
      else if (getCigarType(c) == 'I') {
         q += getCigarLength(c);
      }    
      else if (getCigarType(c) == 'D') {
         r += getCigarLength(c);
      }    
      else {
         cerr << __FILE__ << ":" << __LINE__ << ":ERROR trying to find Deletion at " << oldloc << endl;
         throw runtime_error(string("wrong CigarData state: ") + string(1, getCigarType(c)));
      }    
      ++c; 
   }
   if (getCigarType(c) != 'D' || r != oldloc) {
      cerr << *this << endl;
      throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR cannot find Deletion at " + to_string(oldloc));
   }
   if (c < 1 || c+1 >= getCigarSize() || getCigarType(c+1) != 'M' || getCigarType(c-1) != 'M') {
      throw runtime_error("D must be flanked by M on both sides");
   }
   if (newloc > oldloc) { // for newloc > oldloc
      // r is not first base of right M segment
      if (r + (int)getCigarLength(c+1) - 1 < newloc + (int)getCigarLength(c)) {
         //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG new del ouside alignment\n"
         //   << *this << endl << " r=" << r << " next M seglen=" << getCigarLength(c+1)
         //   << " newloc=" << newloc << " dellen=" << getCigarLength(c) 
         //   << " after adding deletion outside of alignment\n";
         return;
      }
      // actually newloc can be even in the delete segment
      // now do the actual move
      CigarData[c-1].expand(newloc - oldloc);
      CigarData[c+1].shrink(newloc - oldloc);
   }
   else { // newloc < oldloc
      if (newloc <= r- (int)CigarData[c-1].getLength()) {
         return;
      }
      CigarData[c-1].shrink(oldloc - newloc);
      CigarData[c+1].expand(oldloc - newloc);
   }
   if (!validCigar()) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR CigarData and len mismatch\n" << *this << endl;
      throw logic_error("CigarData and query sequence length does not match");
   }
}

// just silently do nothing if not feasible
// need to limit only repeat regions
void BamAlignment::moveInsertion(int oldloc, int newloc) {
   assert(oldloc != newloc);
   if (oldloc < 0 || newloc < 0) {
      if (oldloc == -1) {
        if (newloc < (int)getCigarLength(1)-1 && getCigarType(0) == 'I') {
            // special case 2I145M => 1M2I144M
            // add new M at front
            CigarData[1].shrink(newloc+1);
            CigarData.insert(CigarData.begin(), CigarOp('M', newloc+1));
            SupportData.NumCigarOperations += 1;
            return;
        }
        return; // do nothing
      }
      cerr << *this << endl;
      throw logic_error(string(__FILE__) + ":" + string(__func__) + ":ERROR negative oldloc="
           + to_string(oldloc) + " newloc=" + to_string(newloc));
   }
   pair<int,bool> rtn = isInsertionAtRefloc(newloc);
   if (rtn.second) {
      //cerr << __FILE__ << ":" << __LINE__ << ":WARN " << newloc
      //   << " already has insertion skip moving of insertion\n";
      return;
   }
   rtn = isDeletionAtRefloc(newloc);
   if (rtn.second) {
      //cerr << __FILE__ << ":" << __LINE__ << ":WARN " << newloc
      //   << " cannot move insertion into a deletion\n";
      return;
   }
   rtn = isDeletionAtRefloc(newloc+1);
   if (rtn.second) {
      //cerr << __FILE__ << ":" << __LINE__ << ":WARN " << newloc+1
      //   << " cannot move insertion to the sart of deletion\n";
      return;
   }
   int c=0, r=0, q=0;
   while (c < getCigarSize() && (getCigarType(c) == 'S' || getCigarType(c) == 'H')) {
      q += getCigarLength(c); ++c;
   }
   while (c < getCigarSize() && r <= oldloc && q < getLength()) {
      if (getCigarType(c) == 'M') {
         r += getCigarLength(c);
         q += getCigarLength(c);
      }
      else if (getCigarType(c) == 'I') {
         if (r-1 == oldloc) { // never reach this sate
            // TODO: remove after some testing
            cerr << "?Unreachable code found insertion at " << oldloc << endl;
            --r; // put r at last base of previous M
            break;
         }
         q += getCigarLength(c);
      }
      else if (getCigarType(c) == 'D') {
         r += getCigarLength(c);
      }
      else {
         cerr << "trying to find Insertion at " << oldloc << endl;
         throw runtime_error(string("wrong CigarData state: ") + string(1, getCigarType(c)));
      }
      ++c;
   }
   if (getCigarType(c) != 'I' || r-1 != oldloc) {
      cerr << endl << *this << endl << " oldloc=" << oldloc << " newloc=" << newloc
         << " r=" << r << " c=" << c << endl;
      throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) +
            ":DEBUG cannot find Insertion at " + to_string(oldloc));
   }
   else --r;
   if (c >= getCigarSize()-1 || c < 1 || getCigarType(c-1) != 'M' || getCigarType(c+1) != 'M') {
      cerr << __LINE__ << ":DEBUG c=" << c << " r=" << r << " q=" << q << endl;
      throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) + ":DEBUG I must be flanked by M on both sides");
   }
   if (newloc > oldloc) { // for newloc > oldloc
      if (newloc < r+1 || newloc >= r + (int)getCigarLength(c+1)) {
         //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG newloc not in next M segment. oldloc=" << oldloc
         //   << " newloc=" << newloc << " r=" << r << " q=" << q << " c=" << c << endl
         //   << *this << endl;;
         //throw logic_error(string(__func__) + ":ERROR newloc must be located in the right M segment");
         return;
      }
      try {
         // now do the actual move
         CigarData[c-1].expand(newloc - oldloc);
         CigarData[c+1].shrink(newloc - oldloc);
      }
      catch (exception& er) {
         cerr << __FILE__ << ":" << __LINE__ << ":DEBUG failed to move insertion from "
            << oldloc << " to " << newloc << endl;
         cerr << er.what() << endl << *this << endl;
         throw;
         //exit(1);
      }
   }
   else { // newloc < oldloc
      if (newloc <= r - (int)getCigarLength(c-1) + 2) {
         //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG newloc outof range\n"
         //   << " newloc=" << newloc << " oldloc=" << oldloc
         //   << " q=" << q << " r=" << r << " c=" << c << endl
         //   << *this << endl;
         //throw logic_error(string(__func__) + ":ERROR newloc before start of left M segment");
         return;
      }
      CigarData[c-1].shrink(oldloc - newloc);
      CigarData[c+1].expand(oldloc - newloc);
   }
   if (!validCigar()) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR CigarData and len mismatch\n" << *this << endl;
      throw logic_error("CigarData and query sequence length does not match");
   }
}


pair<int,int> BamAlignment::getMismatchCount() const {
   int32_t numdiff = -1;
   try { // bwa use uint8_t for NM tag, so try this first to
      // reduce the chance of trying signed integer type
      pair<uint32_t,bool> res = getTag<uint32_t>("NM");
      if (res.second) numdiff = res.first;
   }
   catch (const BamTypeException& err) {
      //cerr << __FILE__ << ":" << __LINE__ << ":ERROR " 
      //   << err.what() << endl;
      //uint32_t num=0;
      pair<int32_t,bool> res = getTag<int32_t>("NM");
      if (res.second) numdiff = res.first;
   }
   if (numdiff == -1) {
      cerr << *this << endl;
      throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) 
            +  ":ERROR No NM tag in bam alignment " + getName());
   }
   int alnlen=0;
   int indel=0;
   for (auto& cd : CigarData) {
      if (cd.Type == 'M') {
         alnlen += cd.Length;
      }
      else if (cd.Type == 'D' || cd.Type == 'I') {
         indel += cd.Length;
      }
   }
   if (indel > numdiff) {
      throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR inde greater than numdiff");
   }
   return make_pair(numdiff-indel, alnlen);
}

// if no NM tag then assume 0% identity, could be unaligned read
float BamAlignment::getNGIdentity() const {
   int32_t numMis = -1;
   try {
      //pair<uint32_t,bool> res = getTag<uint32_t>("NM");
      pair<uint8_t,bool> res = getTag<uint8_t>("NM");
      if (res.second) numMis = res.first;
   }
   catch (const BamTypeException& err) {
      pair<int32_t,bool> res = getTag<int32_t>("NM");
      if (res.second) numMis = res.first;
   }
   if (numMis == -1) {
      cerr << *this << endl;
      //throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) +  ":ERROR No NM tag in bam file");
      cerr  << __FILE__ << ":" << __LINE__  << ":ERROR No NM tag in bam alignment\n";
      return 0;
   }
   if (numMis < -1) {
      cerr << *this << endl;
      cerr << __FILE__ << ":" << __LINE__ << ": numMis=" << numMis << endl; 
      throw logic_error("check getTag function");
   }
   int alnlen=0;
   int indel=0;
   for (auto& cd : CigarData) {
      if (cd.Type == 'M') {
         alnlen += cd.Length;
      }
      else if (cd.Type == 'D' || cd.Type == 'I') {
         indel += cd.Length;
      }
   }
   if (indel > numMis) {
      cerr << *this << endl;
      throw logic_error("indel " + to_string(indel) + " greater than total edit distance: " + to_string(numMis));
   }
   return (1-(numMis-indel)/(float)alnlen);
}

float BamAlignment::getIdentity() const {
   int32_t numdiff = -1;
   try {
      pair<uint8_t,bool> res = getTag<uint8_t>("NM");
      if (res.second) numdiff = res.first;
   }
   catch (const BamTypeException& err) {
      //cerr << __LINE__ << ":ERROR " << err.what() << endl;
      //uint32_t num=0;
      pair<int32_t,bool> res = getTag<int32_t>("NM");
      if (res.second) numdiff = res.first;
   }
   if (numdiff == -1) {
      // cerr << *this << endl << " has no NM tag assuming no aligned\n";
      return 0;
      //throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) +  ":ERROR No NM tag in bam file");
   }
   int alnlen=0;
   for (auto& cd : CigarData) {
      if (cd.Type == 'M' || cd.Type == 'D' || cd.Type == 'I') {
         alnlen += cd.Length;
      }
      else if (cd.Type == 'S' || cd.Type == 'H') {
         // ignored
      }
      else {
         cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
            << " Cigarop: " << cd.Type << " not added to alignment length\n";
      }
   }
   return (1-(numdiff)/(float)alnlen);
}

void BamAlignment::setQuality(const vector<int> &qual) {
   if (!Qualities.empty()) Qualities.clear();
   //cout << "Quality values:\n";
   for (size_t i=0; i<qual.size(); ++i) {
      //cout << qual[i] << " ";
      Qualities.append(1, char(qual[i]+33));
   }
   //cout << endl;
   // for debug
   //cout << "quality string after setQuality() call\n"
   //   << Qualities << endl;
}

void BamAlignment::fillQuality(int score) {
   if (Qualities.size() < QueryBases.size()) {
      size_t B=Qualities.size();
      Qualities.resize(QueryBases.size());
      for (size_t i=B; i < QueryBases.size(); ++i) {
         Qualities[i]= char(score+33);
      }
   }
}

/*! \fn bool BamAlignment::BuildCharData(void)
    \brief Populates alignment string fields (read name, bases, qualities, tag data).

    An alignment retrieved using BamReader::GetNextAlignmentCore() lacks this data.
    Using that method makes parsing much quicker when only positional data is required.

    However, if you later want to access the character data fields from such an alignment,
    use this method to populate those fields. Provides ability to do 'lazy evaluation' of
    alignment parsing.

    \return \c true if character data populated successfully (or was already available to begin with)
*/
bool BamAlignment::BuildCharData(void) {
    // skip if char data already parsed
    if ( !SupportData.HasCoreOnly )
        return true; // already parsed, no repeat work
    // check system endianness
    bool IsBigEndian = BamTools::SystemIsBigEndian();
    // calculate character lengths/offsets
    const unsigned int dataLength     = SupportData.BlockLength - Constants::BAM_CORE_SIZE;
    const unsigned int seqDataOffset  = SupportData.QueryNameLength + (SupportData.NumCigarOperations*4);
    const unsigned int qualDataOffset = seqDataOffset + (SupportData.QuerySequenceLength+1)/2;
    const unsigned int tagDataOffset  = qualDataOffset + SupportData.QuerySequenceLength;
    const unsigned int tagDataLength  = dataLength - tagDataOffset;

    // check offsets to see what char data exists
    const bool hasSeqData  = ( seqDataOffset  < qualDataOffset );
    const bool hasQualData = ( qualDataOffset < tagDataOffset );
    const bool hasTagData  = ( tagDataOffset  < dataLength );

    // store alignment name (relies on null char in name as terminator)
    Name.assign(SupportData.AllCharData.data());

    // save query sequence
    QueryBases.clear();
    if ( hasSeqData ) {
        const char* seqData = SupportData.AllCharData.data() + seqDataOffset;
        QueryBases.reserve(SupportData.QuerySequenceLength);
        for ( size_t i = 0; i < SupportData.QuerySequenceLength; ++i ) {
            const char singleBase = Constants::BAM_DNA_LOOKUP[ ( (seqData[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
            QueryBases.append(1, singleBase);
        }
    }
    // save qualities, quality stored in the data structure have 33 added
    Qualities.clear();
    if ( hasQualData ) {
        const char* qualData = SupportData.AllCharData.data() + qualDataOffset;
        // if marked as unstored (sequence of 0xFF) - don't do conversion, just fill with 0xFFs
        if ( qualData[0] == (char)0xFF )
            Qualities.resize(SupportData.QuerySequenceLength, (char)0xFF);
        // otherwise convert from numeric QV to 'FASTQ-style' ASCII character
        else {
            Qualities.reserve(SupportData.QuerySequenceLength);
            for ( size_t i = 0; i < SupportData.QuerySequenceLength; ++i )
                //Qualities.append(1, qualData[i]); // no change to ASCII value
                Qualities.append(1, qualData[i]+33);
        }
    }

    // clear previous AlignedBases
    AlignedBases.clear();

    // if QueryBases has data, build AlignedBases using CIGAR data
    // otherwise, AlignedBases will remain empty (this case IS allowed)
    if ( !QueryBases.empty() && QueryBases != "*" ) {
        // resize AlignedBases
        AlignedBases.reserve(SupportData.QuerySequenceLength);

        // iterate over CigarOps
        int k = 0;
        vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
        vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
        for ( ; cigarIter != cigarEnd; ++cigarIter ) {
            const CigarOp& op = (*cigarIter);

            switch ( op.Type ) {

                // for 'M', 'I', '=', 'X' - write bases
                case (Constants::BAM_CIGAR_MATCH_CHAR)    :
                case (Constants::BAM_CIGAR_INS_CHAR)      :
                case (Constants::BAM_CIGAR_SEQMATCH_CHAR) :
                case (Constants::BAM_CIGAR_MISMATCH_CHAR) :
                    AlignedBases.append(QueryBases.substr(k, op.Length));
                    // fall through
                // for 'S' - soft clip, do not write bases
                // but increment placeholder 'k'
                case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) :
                    k += op.Length;
                    break;
                // for 'D' - write gap character
                case (Constants::BAM_CIGAR_DEL_CHAR) :
                    AlignedBases.append(op.Length, Constants::BAM_DNA_DEL);
                    break;
                // for 'P' - write padding character
                case (Constants::BAM_CIGAR_PAD_CHAR) :
                    AlignedBases.append( op.Length, Constants::BAM_DNA_PAD );
                    break;
                // for 'N' - write N's, skip bases in original query sequence
                case (Constants::BAM_CIGAR_REFSKIP_CHAR) :
                    AlignedBases.append( op.Length, Constants::BAM_DNA_N );
                    break;
                // for 'H' - hard clip, do nothing to AlignedBases, move to next op
                case (Constants::BAM_CIGAR_HARDCLIP_CHAR) :
                    break;
                // invalid CIGAR op-code
                default:
                    cerr << __FILE__ << ":" << __LINE__ << ":ERROR invalid CIGAR operation type: "
                       << op.Type << endl;
                    return false;
            }
        }
    }

    // save tag data
    TagData.clear();
    if ( hasTagData ) {

        char* tagData = (((char*)SupportData.AllCharData.data()) + tagDataOffset);

        if ( IsBigEndian ) {
            size_t i = 0;
            while ( i < tagDataLength ) {
               //if (tagData[i] == 'N' && tagData[i+1] == 'M') {
               //   cerr << __func__ << ": NM type: " << tagData[i+2] << endl;
               //}
                i += Constants::BAM_TAG_TAGSIZE;  // skip tag chars (e.g. "RG", "NM", etc.)
                const char type = tagData[i];     // get tag type at position i
                ++i;                              // move i past tag type

                switch (type) {

                    case(Constants::BAM_TAG_TYPE_ASCII) :
                    case(Constants::BAM_TAG_TYPE_INT8)  :
                    case(Constants::BAM_TAG_TYPE_UINT8) :
                        // no endian swapping necessary for single-byte data
                        ++i;
                        break;

                    case(Constants::BAM_TAG_TYPE_INT16)  :
                    case(Constants::BAM_TAG_TYPE_UINT16) :
                        BamTools::SwapEndian_16p(&tagData[i]);
                        i += sizeof(uint16_t);
                        break;

                    case(Constants::BAM_TAG_TYPE_FLOAT)  :
                    case(Constants::BAM_TAG_TYPE_INT32)  :
                    case(Constants::BAM_TAG_TYPE_UINT32) :
                        BamTools::SwapEndian_32p(&tagData[i]);
                        i += sizeof(uint32_t);
                        break;

                    case(Constants::BAM_TAG_TYPE_HEX) :
                    case(Constants::BAM_TAG_TYPE_STRING) :
                        // no endian swapping necessary for hex-string/string data
                        while ( tagData[i] )
                            ++i;
                        // increment one more for null terminator
                        ++i;
                        break;

                    case(Constants::BAM_TAG_TYPE_ARRAY) :

                    {
                        // read array type
                        const char arrayType = tagData[i];
                        ++i;

                        // swap endian-ness of number of elements in place, then retrieve for loop
                        BamTools::SwapEndian_32p(&tagData[i]);
                        uint32_t numElements;
                        memcpy(&numElements, &tagData[i], sizeof(uint32_t));
                        i += sizeof(uint32_t);

                        // swap endian-ness of array elements
                        for ( size_t j = 0; j < numElements; ++j ) {
                            switch (arrayType) {
                                case (Constants::BAM_TAG_TYPE_INT8)  :
                                case (Constants::BAM_TAG_TYPE_UINT8) :
                                    // no endian-swapping necessary
                                    ++i;
                                    break;
                                case (Constants::BAM_TAG_TYPE_INT16)  :
                                case (Constants::BAM_TAG_TYPE_UINT16) :
                                    BamTools::SwapEndian_16p(&tagData[i]);
                                    i += sizeof(uint16_t);
                                    break;
                                case (Constants::BAM_TAG_TYPE_FLOAT)  :
                                case (Constants::BAM_TAG_TYPE_INT32)  :
                                case (Constants::BAM_TAG_TYPE_UINT32) :
                                    BamTools::SwapEndian_32p(&tagData[i]);
                                    i += sizeof(uint32_t);
                                    break;
                                default:
                                    cerr << __FILE__ << ":" << __LINE__ << ":WARN invalid binary array type : "
                                       << arrayType << endl;
                                    return false;
                            }
                        }

                        break;
                    }

                    // invalid tag type-code
                    default :
                        cerr << __FILE__ << ":" << __LINE__ << "ERROR: invalid tag type: " << type << endl;
                        return false;
                }
            }
        }

        // store tagData in alignment
        TagData.resize(tagDataLength);
        memcpy((char*)(TagData.data()), tagData, tagDataLength);
    }

    // clear core-only flag & return success
    SupportData.HasCoreOnly = false;
    return true;
}

bool BamAlignment::FindTag(const std::string& tag, char*& pTagData,
      const unsigned int& tagDataLength, unsigned int& numBytesParsed) const
{
    while (numBytesParsed < tagDataLength) {
        const char* pTagType        = pTagData;
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;
        // check the current tag, return true on match
        if (strncmp(pTagType, tag.c_str(), 2) == 0)
            return true;
        // get the storage class and find the next tag
        if (*pTagStorageType == '\0') return false;
        //if (!SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed)) 
        if (!SkipToNextTag(pTagData, numBytesParsed)) 
           return false;
        if (*pTagData == '\0') return false;
    }
    return false;
}

size_t BamAlignment::getTagWidth(const char* p) const {
   try {
      p += 2; // move to data type char
      size_t w;
      if (*p == Constants::BAM_TAG_TYPE_ARRAY) {
         w = getArrayTagLength(p+1) + 3;
      }
      else {
         w = getBasicTagLength(p) + 3;
      }
      return w;
   }
   catch (const BamTypeException& ler) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR failed inside "
         << __func__ << " " << ler.what() << endl;
      throw BamTypeException(string(ler.what()) + " failed getTagWidth()");
   }
}

void BamAlignment::showTagData(ostream& ous) const {
   for (const char c : TagData) {
      if (c == '\0') ous << ".";
      else if (isprint(c)) 
         ous << c;
      else 
         ous << '~';
   }
   ous << endl;
}

char* BamAlignment::findTag(const std::string& tag) {
   //cerr << "finding tag " << tag << endl;
   char* p = TagData.data();
   while (*p != '\0') {
      //cerr << "at " << p << endl;
      if (tag[0] == *p && tag[1] == *(p+1)) {
         return p;
      }
      p += getTagWidth(p);
   }
   return nullptr;
}

// assume tag data is constructed
const char* BamAlignment::findTag(const std::string& tag) const
{
   //cerr << __LINE__ << ": finding tag " << tag << endl;
   const char* p = TagData.data();
   while (*p != '\0') {
      //cerr << "at " << p << endl;
      if (tag[0] == *p && tag[1] == *(p+1))
         return p;
      p += getTagWidth(p);
   }
   return nullptr;
}


/*! \fn bool BamAlignment::GetArrayTagType(const std::string& tag, char& type) const
    \brief Retrieves the BAM tag type-code for the array elements associated with requested tag name.

    \param[in]  tag  2-character tag name
    \param[out] type retrieved (1-character) type-code

    \return \c true if found. False if not found, or if tag is not an array type.
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::GetArrayTagType(const std::string& tag, char& type) const {

    // skip if alignment is core-only
    if ( SupportData.HasCoreOnly ) {
        // TODO: set error string?
        return false;
    }

    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag not found, return failure
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ){
        // TODO: set error string?
        return false;
    }

    // check that tag type code is array
    type = *(pTagData - 1);
    if ( type != Constants::BAM_TAG_TYPE_ARRAY ) {
        // TODO: set error string
        return false;
    }

    // fetch element type
    const char elementType = *pTagData;
    switch ( elementType ) {

        // allowable types
        case (Constants::BAM_TAG_TYPE_INT8)   :
        case (Constants::BAM_TAG_TYPE_UINT8)  :
        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
            type = elementType;
            break;

        default:
            //TODO: set error string
            return false;
    }

    // if we get here, return success
    return true;
}


/*! \fn int BamAlignment::GetEndPosition(bool usePadded = false, bool closedInterval = false) const
    \brief Calculates alignment end position, based on its starting position and CIGAR data.

    \warning The position returned now represents a zero-based, HALF-OPEN interval.
    In previous versions of BamTools (0.x & 1.x) all intervals were treated
    as zero-based, CLOSED.

    \param[in] usePadded      Allow inserted bases to affect the reported position. Default is
                              false, so that reported position stays synced with reference
                              coordinates.
    \param[in] closedInterval Setting this to true will return a 0-based end coordinate. Default is
                              false, so that his value represents a standard, half-open interval.

    \return alignment end position
*/
int BamAlignment::GetEndPosition(bool usePadded, bool closedInterval) const {
    // initialize alignment end to starting position
    int alignEnd = Position;
    // iterate over cigar operations
    //vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    //vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    //for ( ; cigarIter != cigarEnd; ++cigarIter) {
    for (auto& op : CigarData) {
        //const CigarOp& op = (*cigarIter);
        //switch ( op.Type ) {
        switch ( op.Type ) {
            // increase end position on CIGAR chars [DMXN=]
            case Constants::BAM_CIGAR_DEL_CHAR      :
            case Constants::BAM_CIGAR_MATCH_CHAR    :
            case Constants::BAM_CIGAR_MISMATCH_CHAR :
            case Constants::BAM_CIGAR_REFSKIP_CHAR  :
            case Constants::BAM_CIGAR_SEQMATCH_CHAR :
                alignEnd += op.Length;
                break;
            // increase end position on insertion, only if @usePadded is true
            case Constants::BAM_CIGAR_INS_CHAR :
                if ( usePadded )
                    alignEnd += op.Length;
                break;
            // all other CIGAR chars do not affect end position
            default :
                break;
        }
    }
    // adjust for closedInterval, if requested
    if ( closedInterval ) alignEnd -= 1;
    // return result
    return alignEnd;
}

/*! \fn std::string BamAlignment::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
//std::string BamAlignment::GetErrorString(void) const {
//    return ErrorString;
//}

/*! \fn bool BamAlignment::GetSoftClips(std::vector<int>& clipSizes, std::vector<int>& readPositions, std::vector<int>& genomePositions, bool usePadded = false) const
    \brief Identifies if an alignment has a soft clip. If so, identifies the
           sizes of the soft clips, as well as their positions in the read and reference.

    \param[out] clipSizes       vector of the sizes of each soft clip in the alignment
    \param[out] readPositions   vector of the 0-based read locations of each soft clip in the alignment.
                                These positions are basically indexes within the read, not genomic positions.
    \param[out] genomePositions vector of the 0-based genome locations of each soft clip in the alignment
    \param[in]  usePadded       inserted bases affect reported position. Default is false, so that
                                reported position stays 'sync-ed' with reference coordinates.

    \return \c true if any soft clips were found in the alignment
*/
bool BamAlignment::GetSoftClips(vector<int>& clipSizes,
                                vector<int>& readPositions,
                                vector<int>& genomePositions,
                                bool usePadded) const
{
    // initialize positions & flags
    int refPosition  = Position;
    int readPosition = 0;
    bool softClipFound = false;
    bool firstCigarOp  = true;

    // iterate over cigar operations
    vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const CigarOp& op = (*cigarIter);

        switch ( op.Type ) {

            // increase both read & genome positions on CIGAR chars [DMXN=]
            case Constants::BAM_CIGAR_DEL_CHAR      :
            case Constants::BAM_CIGAR_MATCH_CHAR    :
            case Constants::BAM_CIGAR_MISMATCH_CHAR :
            case Constants::BAM_CIGAR_REFSKIP_CHAR  :
            case Constants::BAM_CIGAR_SEQMATCH_CHAR :
                refPosition  += op.Length;
                readPosition += op.Length;
                break;

            // increase read position on insertion, genome position only if @usePadded is true
            case Constants::BAM_CIGAR_INS_CHAR :
                readPosition += op.Length;
                if ( usePadded )
                    refPosition += op.Length;
                break;

            case Constants::BAM_CIGAR_SOFTCLIP_CHAR :

                softClipFound = true;

                //////////////////////////////////////////////////////////////////////////////
                // if we are dealing with the *first* CIGAR operation
                // for this alignment, we increment the read position so that
                // the read and genome position of the clip are referring to the same base.
                // For example, in the alignment below, the ref position would be 4, yet
                //              the read position would be 0. Thus, to "sync" the two,
                //              we need to increment the read position by the length of the
                //              soft clip.
                // Read:  ATCGTTTCGTCCCTGC
                // Ref:   GGGATTTCGTCCCTGC
                // Cigar: SSSSMMMMMMMMMMMM
                //
                // NOTE: This only needs to be done if the soft clip is the _first_ CIGAR op.
                //////////////////////////////////////////////////////////////////////////////
                if ( firstCigarOp )
                    readPosition += op.Length;

                // track the soft clip's size, read position, and genome position
                clipSizes.push_back(op.Length);
                readPositions.push_back(readPosition);
                genomePositions.push_back(refPosition);

            // any other CIGAR operations have no effect
            default :
                break;
        }

        // clear our "first pass" flag
        firstCigarOp = false;
    }

    // return whether any soft clips found
    return softClipFound;
}

string BamAlignment::getStringTag(const std::string& tag) const {
    if (TagData.empty() || SupportData.HasCoreOnly) {
        return string();
    }
    const char* pTag = findTag(tag);
    if (pTag == nullptr) {
       return string();
    }
    else {
       pTag += (Constants::BAM_TAG_TAGSIZE + Constants::BAM_TAG_TYPESIZE);
       return string(pTag);
    }
    /*
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // return failure if tag not found
    if (!FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) +
              ":ERROR Cannot find tag: " + tag);
    }
    // otherwise copy data into destination
    const unsigned int dataLength = strlen(pTagData);
    return TagData.substr(numBytesParsed, dataLength);
    */
}


/* \fn std::vector<std::string> BamAlignment::GetTagNames(void) const
    \brief Retrieves the BAM tag names.

    When paired with GetTagType() and GetTag(), this method allows you
    to iterate over an alignment's tag data without knowing the names (or types)
    beforehand.

    \return \c vector containing all tag names found (empty if none available)
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
std::vector<std::string> BamAlignment::GetTagNames(void) const {
    std::vector<std::string> result;
    if ( SupportData.HasCoreOnly || TagData.empty() )
        return result;
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    while ( numBytesParsed < tagDataLength ) {
        // get current tag name & type
        const char* pTagName = pTagData;
        const char* pTagType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed +=3;
        // store tag name
        result.push_back( std::string(pTagName, 2)  );
        // find the next tag
        if ( *pTagType == '\0' ) break; // end of input
        //if (!SkipToNextTag(*pTagType, pTagData, numBytesParsed)) 
        if (!SkipToNextTag(pTagData, numBytesParsed)) 
           break;
        if (*pTagData == '\0') break;
    }
    return result;
}

/*! \fn bool BamAlignment::GetTagType(const std::string& tag, char& type) const
    \brief Retrieves the BAM tag type-code associated with requested tag name.

    \param[in]  tag  2-character tag name
    \param[out] type retrieved (1-character) type-code

    \return \c true if found
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
bool BamAlignment::GetTagType(const std::string& tag, char& type) const {
    // skip if alignment is core-only
    if ( SupportData.HasCoreOnly ) {
        // TODO: set error string?
        return false;
    }
    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }
    if (!isValidTagName(tag)) {
       cerr << *this << endl;
       cerr << __FILE__ << ":" << __LINE__ << ":ERROR tag " << tag 
          << " is not valid\n";
       //showTagData(cerr);
       throw BamNotagException("invalide tag name");
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // if tag not found, return failure
    if (!FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        // TODO: set error string?
        return false;
    }
    // otherwise, retrieve & validate tag type code
    type = *(pTagData - 1);
    switch (type) {
        case (Constants::BAM_TAG_TYPE_ASCII)  :
        case (Constants::BAM_TAG_TYPE_INT8)   :
        case (Constants::BAM_TAG_TYPE_UINT8)  :
        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            return true;
        // unknown tag type
        default:
            cerr << __FILE__ << ":" << __LINE__ << ":ERROR invalid tag type: " << type << endl;
            return false;
    }
}

bool BamAlignment::HasTag(const std::string& tag) const {
    // return false if no tag data present
    if ( SupportData.HasCoreOnly || TagData.empty() )
        return false;
    // localize the tag data for lookup
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if result of tag lookup
    return FindTag(tag, pTagData, tagDataLength, numBytesParsed);
}

bool BamAlignment::hasTag(const std::string& tag) const {
    // return false if no tag data present
    if ( SupportData.HasCoreOnly || TagData.empty() )
        return false;
    // localize the tag data for lookup
    const char* p = findTag(tag);
    return p != nullptr;
}

bool BamAlignment::IsDuplicate(void) const {
    return (AlignmentFlag & Constants::BAM_ALIGNMENT_DUPLICATE) == Constants::BAM_ALIGNMENT_DUPLICATE;
}

/*! \fn bool BamAlignment::IsFailedQC(void) const
    \return \c true if this read failed quality control
*/
bool BamAlignment::IsFailedQC(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_QC_FAILED) != 0 );
}

bool BamAlignment::IsFirstMate(void) const {
   //bool x = ((AlignmentFlag & Constants::BAM_ALIGNMENT_READ_1) != 0);
   //if (x) {
   //   cerr << " x is true\n";
   //}
   //else cerr << " x is false\n";
    return (AlignmentFlag & Constants::BAM_ALIGNMENT_READ_1) == Constants::BAM_ALIGNMENT_READ_1;
}

/*! \fn bool BamAlignment::IsMapped(void) const
    \return \c true if alignment is mapped
*/
bool BamAlignment::IsMapped(void) const {
    return !((AlignmentFlag & Constants::BAM_ALIGNMENT_UNMAPPED) == Constants::BAM_ALIGNMENT_UNMAPPED);
}

/*! \fn bool BamAlignment::IsMateMapped(void) const
    \return \c true if alignment's mate is mapped
*/
bool BamAlignment::IsMateMapped(void) const {
    return !((AlignmentFlag & Constants::BAM_ALIGNMENT_MATE_UNMAPPED) == Constants::BAM_ALIGNMENT_MATE_UNMAPPED);
}

/*! \fn bool BamAlignment::IsMateReverseStrand(void) const
    \return \c true if alignment's mate mapped to reverse strand
*/
bool BamAlignment::IsMateReverseStrand(void) const {
    return (AlignmentFlag & Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND) == Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND;
}

double BamAlignment::getFractionStrand() const {
   //int32_t overlap = -1;
   pair<int32_t,bool> res = getTag<int32_t>("XO");
   if (res.second) {
      if (res.first == getReferenceWidth()) {
         return 0;
      }
      else {
         return 1-static_cast<double>(res.first)/getReferenceWidth();
      }
   }
   else {
      if (isReverseStrand()) return -1;
      else return 1;
   }
}

/*! \fn bool BamAlignment::IsPaired(void) const
    \return \c true if alignment part of paired-end read
*/
bool BamAlignment::IsPaired(void) const {
    return (AlignmentFlag & Constants::BAM_ALIGNMENT_PAIRED) == Constants::BAM_ALIGNMENT_PAIRED;
}

/*! \fn bool BamAlignment::IsPrimaryAlignment(void) const
    \return \c true if reported position is primary alignment
*/
bool BamAlignment::IsPrimaryAlignment(void) const  {
    return !((AlignmentFlag & Constants::BAM_ALIGNMENT_SECONDARY) == Constants::BAM_ALIGNMENT_SECONDARY);
}

/*! \fn bool BamAlignment::IsProperPair(void) const
    \return \c true if alignment is part of read that satisfied paired-end resolution
*/
bool BamAlignment::IsProperPair(void) const {
    //return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_PROPER_PAIR) != 0 );
    return (AlignmentFlag & PROPER_PAIR) == PROPER_PAIR;
}

bool BamAlignment::IsReverseStrand(void) const {
    return (AlignmentFlag & Constants::BAM_ALIGNMENT_REVERSE_STRAND) == Constants::BAM_ALIGNMENT_REVERSE_STRAND;
}

/*! \fn bool BamAlignment::IsSecondMate(void) const
    \return \c true if alignment is second mate on read
*/
bool BamAlignment::IsSecondMate(void) const {
    return (AlignmentFlag & Constants::BAM_ALIGNMENT_READ_2) == Constants::BAM_ALIGNMENT_READ_2;
}

/*! \fn bool BamAlignment::IsValidSize(const std::string& tag, const std::string& type) const
    \internal

    Checks that tag name & type strings are expected sizes.

    \param tag[in]  BAM tag name
    \param type[in] BAM tag type-code
    \return \c true if both input strings are valid sizes
*/
bool BamAlignment::IsValidSize(const std::string& tag, const std::string& type) const {
    return (tag.size()  == Constants::BAM_TAG_TAGSIZE) &&
           (type.size() == Constants::BAM_TAG_TYPESIZE);
}

bool BamAlignment::isValidArrayTag(const string& tag) const {
   if (tag.size() != Constants::BAM_TAG_TAGSIZE) 
      return false;
   const char* p = findTag(tag);
   if (p == nullptr) {
      return false;
   }
   return isValidArrayTag(p);
}

void BamAlignment::removeTag(const std::string& tag) {
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();
    // skip if no tags available
    if (TagData.empty()) {
       cerr << __FILE__ << ":" << __LINE__ << ":WARN no tag data while removing tag: "
         << tag << endl;
       return;
    }
    char* pTag = findTag(tag); // pTag pointer to first char of TAG
    if (pTag == nullptr) {
       //throw BamNotagException("there is no tag " + tag + " to remove");
       //cerr << __FILE__ << ":" << __LINE__ << ":WARN there is no tag " + tag + " to remove\n";
       return;
    }
    size_t tglen = getTagWidth(pTag); // num char of whole tag entry
    size_t newtgdLength = TagData.size() - tglen;
    if (*(pTag + tglen) != '\0') { // tag is not the last one
       memmove(pTag, pTag+tglen, newtgdLength - (pTag - TagData.data()));
    }
    TagData.resize(newtgdLength);
}

/*! \fn void BamAlignment::SetIsDuplicate(bool ok)
    \brief Sets value of "PCR duplicate" flag to \a ok.
*/
void BamAlignment::SetIsDuplicate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_DUPLICATE;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_DUPLICATE;
}

/*! \fn void BamAlignment::SetIsFailedQC(bool ok)
    \brief Sets "failed quality control" flag to \a ok.
*/
void BamAlignment::SetIsFailedQC(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_QC_FAILED;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_QC_FAILED;
}

/*! \fn void BamAlignment::SetIsFirstMate(bool ok)
    \brief Sets "alignment is first mate" flag to \a ok.
*/
void BamAlignment::SetIsFirstMate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_READ_1;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_READ_1;
}

/*! \fn void BamAlignment::SetIsMapped(bool ok)
    \brief Sets "alignment is mapped" flag to \a ok.
*/
void BamAlignment::SetIsMapped(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_UNMAPPED;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_UNMAPPED;
}

/*! \fn void BamAlignment::SetIsMateMapped(bool ok)
    \brief Sets "alignment's mate is mapped" flag to \a ok.
*/
void BamAlignment::SetIsMateMapped(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_MATE_UNMAPPED;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_MATE_UNMAPPED;
}

/*! \fn void BamAlignment::SetIsMateReverseStrand(bool ok)
    \brief Sets "alignment's mate mapped to reverse strand" flag to \a ok.
*/
void BamAlignment::SetIsMateReverseStrand(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND;
}

/*! \fn void BamAlignment::SetIsPaired(bool ok)
    \brief Sets "alignment part of paired-end read" flag to \a ok.
*/
void BamAlignment::SetIsPaired(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_PAIRED;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_PAIRED;
}

/*! \fn void BamAlignment::SetIsPrimaryAlignment(bool ok)
    \brief Sets "position is primary alignment" flag to \a ok.
*/
void BamAlignment::SetIsPrimaryAlignment(bool ok) {
    if (ok) AlignmentFlag &= ~Constants::BAM_ALIGNMENT_SECONDARY;
    else    AlignmentFlag |=  Constants::BAM_ALIGNMENT_SECONDARY;
}

/*! \fn void BamAlignment::SetIsProperPair(bool ok)
    \brief Sets "alignment is part of read that satisfied paired-end resolution" flag to \a ok.
*/
void BamAlignment::SetIsProperPair(bool ok) {
    //if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_PROPER_PAIR;
    //else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_PROPER_PAIR;
    if (ok) AlignmentFlag |=  PROPER_PAIR;
    else    AlignmentFlag &= ~PROPER_PAIR;
}

/*! \fn void BamAlignment::SetIsReverseStrand(bool ok)
    \brief Sets "alignment mapped to reverse strand" flag to \a ok.
*/
void BamAlignment::SetIsReverseStrand(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_REVERSE_STRAND;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_REVERSE_STRAND;
}

/*! \fn void BamAlignment::SetIsSecondMate(bool ok)
    \brief Sets "alignment is second mate on read" flag to \a ok.
*/
void BamAlignment::SetIsSecondMate(bool ok) {
    if (ok) AlignmentFlag |=  Constants::BAM_ALIGNMENT_READ_2;
    else    AlignmentFlag &= ~Constants::BAM_ALIGNMENT_READ_2;
}

short BamAlignment::getAtomicTagLength(const char tagt) const {
   switch (tagt) {
      case Constants::BAM_TAG_TYPE_ASCII:
      case Constants::BAM_TAG_TYPE_INT8:
      case Constants::BAM_TAG_TYPE_UINT8:
         return 1; break;
      case Constants::BAM_TAG_TYPE_INT16:
      case Constants::BAM_TAG_TYPE_UINT16:
         return sizeof(uint16_t); break;
      case Constants::BAM_TAG_TYPE_INT32:
      case Constants::BAM_TAG_TYPE_UINT32:
      case Constants::BAM_TAG_TYPE_FLOAT:
         return sizeof(uint32_t); break;
      default:
         cerr << __FILE__ << ":" << __LINE__ << ":ERROR atomic tag error\n" 
            << *this << endl;
         throw BamTypeException(string(__func__) + ":ERROR " + string(tagt, 1) + " is not elemental Bam TAG type");
   }
}

// p point to the type char field 3rd char
size_t BamAlignment::getBasicTagLength(const char* p) const {
   switch (*p) {
      case Constants::BAM_TAG_TYPE_ASCII:
      case Constants::BAM_TAG_TYPE_INT8:
      case Constants::BAM_TAG_TYPE_UINT8:
         return 1; 
      case Constants::BAM_TAG_TYPE_INT16:
      case Constants::BAM_TAG_TYPE_UINT16:
         return sizeof(uint16_t);
      case Constants::BAM_TAG_TYPE_INT32:
      case Constants::BAM_TAG_TYPE_UINT32:
      case Constants::BAM_TAG_TYPE_FLOAT:
         return sizeof(uint32_t); 
      case (Constants::BAM_TAG_TYPE_STRING) :
      case (Constants::BAM_TAG_TYPE_HEX)    :
      {
         const char* x=p+1;
         while(*x != '\0') ++x; // x now at null-terminator
         return x - p; // length of data includding \0
      }
      default:
         throw logic_error(string(*p, 1) + " is not basic Bam TAG type");
   }
}

bool BamAlignment::isValidTag(const char* p) const {
   for (auto i=0; i < Constants::BAM_TAG_TAGSIZE; ++i) {
      if (!isalpha(*p)) return false;
      ++p;
   }
   if (*p == Constants::BAM_TAG_TYPE_ARRAY) {
      return Constants::isAtomicBamTagType(*(p+1));
   }
   else {
      return Constants::isBasicBamTagType(*p);
   }
}

bool BamAlignment::isValidArrayTag(const char* p) const {
   try {
      size_t w = getArrayTagLength(p+3) + 3;
      p += w;
      if (*p == '\0' || isValidTag(p)) {
         return true;
      }
      return false;
   }
   catch (const BamTypeException& ler) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR failed inside isValidArrayTag()\n"
         << ler.what() << endl;
      throw BamTypeException(string(ler.what()) + " inside isValidArrayTag()");
   }
}

// p ponter at array element type char (4th char)
size_t BamAlignment::getArrayTagLength(const char* p) const {
   try {
      short elemlen = getAtomicTagLength(*p);
      const char* x=p+1; // [ 4 byte int32_t ] for number of elements
      int32_t numelement;
      memcpy(&numelement, x, sizeof(uint32_t)); // already endian-swapped, if needed
      return numelement*elemlen + Constants::BAM_TAG_ARRAYBASE_SIZE - 3;
   }
   catch (const BamTypeException& ler) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR failed to getArrayTagLength()\n"
         << ler.what() << endl << p << endl;
      showTagData(cerr);
      cerr << endl;
      throw BamTypeException(string(ler.what()) + " Failed getArrayTagLength()");
   }
   return 0;
}

// pTagData at first char of data, for array is at element type char.
//bool BamAlignment::SkipToNextTag(const char storageType, char*& pTagData, unsigned int& numBytesParsed) const
bool BamAlignment::SkipToNextTag(char*& pTagData, unsigned int& numBytesParsed) const
{
   size_t dlen;
   try {
      if (*(pTagData-1) != Constants::BAM_TAG_TYPE_ARRAY) {
         dlen = getBasicTagLength(pTagData-1);
      }
      else { 
         dlen = getArrayTagLength(pTagData);
      }
      numBytesParsed += dlen;
      pTagData += dlen; // if string type will be at \0
   }
   catch (logic_error& err) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR " << err.what() 
         << " cannot skip to next tag " << endl;
      return false;
   }
   return true;
}

// BamAlignment stores the ASCII values, Phred+33
vector<int> BamAlignment::getQualityScore() const {
   vector<int> qual(Qualities.size());
   for (string::size_type i=0; i<Qualities.size(); ++i) {
      qual[i] = Qualities[i] - 33;
   }
   return qual;
}

// there is a potential for overflow for long sequences
int BamAlignment::getAverageQualityScore() const {
   vector<int> q = getQualityScore();
   return accumulate(q.begin(), q.end(), 0)/float(q.size());
}

std::pair<int,int> BamAlignment::getInterval() const { 
   pair<int,int> tmp = std::pair<int,int>(getPosition(), GetEndPosition(false, true)); 
   if (tmp.first > tmp.second) {
      swap(tmp.first, tmp.second);
   }
   return tmp;
}

bool BamAlignment::sameInterval(const BamAlignment& ba) const {
   pair<int,int> itv1 = getInterval();
   pair<int,int> itv2 = ba.getInterval();
   return itv1.first == itv2.first && itv1.second == itv2.second;
}

std::pair<int,int> BamAlignment::getSoftInterval() const {
   auto res = getInterval();
   int tmp = getFirstSoftclipLength();
   if (tmp > 0) {
      //assert(res.first > tmp);
      res.first -= tmp;
   }
   tmp = getLastSoftclipLength();
   if (tmp > 0) {
      res.second += tmp;
   }
   return res;
}

std::pair<int,int> BamAlignment::getPairedRange() const {
   if (!mateOnSameReference()) {
      //return getRange();
      return getInterval();
   }
   int b, e;
   if (IsReverseStrand()) { // read1 - 
      if (isMateReverseStrand()) { // read 2 - <--R1--  <--R2--
         if (getPosition() < getMatePosition()) {
            b=getEndPosition();
            e = b + getInsertSize() - 1;
         }
         else if (getPosition() == getMatePosition()) {
            // <--R1--      or   <--R1--------
            // <--R2-----        <--R2--
            b=getEndPosition();
            e=b+getInsertSize();
         }
         else {  // <--R2-- <--R1--
            e = getEndPosition();
            b= e-abs(getInsertSize()) + 1; // insert Size is negative for the last one
         }
      }
      else { // read 2 +   -/+
         if (getEndPosition() <= getMatePosition()) {
            // <--R1--     --R2-->
            //       |--I--|
            b=getEndPosition();
            e=getMatePosition();
            //cerr << "improper mapped read pair <--R1-- <--R2--\n"
            //   << b << "-" << e << endl;
            //cerr << *this << endl;
            //exit(1);
         }
         else { // read 2 + --R2--> <--R1--
            b=getMatePosition();
            e=getEndPosition();
         }
      }
   }
   else { // + strand --R1-->
      if (isMateReverseStrand()) { // Read 2 -
        if (getPosition() <= getMatePosition() + getQueryLength()) { // --R1--> <--R2--
            // this is an estimate of Read 2 length with read 1 length
            // this is the defect of Sam/Bam schema
            b=getPosition();
            e=b+getInsertSize()-1;
        }
        else { // <--R2--   --R1--> 
           //           b   e   cannot get e exactly without R2 length
           e=getPosition();
           b=e-abs(getInsertSize()) + 1; // how bwa compute it
        }
      }
      else { // --R1--> --R2-->
         if (getPosition() < getMatePosition()) {
            b=getPosition();
            e=getMatePosition();
         }
         else {
            e=getPosition();
            b=getMatePosition();
         }
      }
   }
   return make_pair(b,e);
}

// compute the reference length from the MC tag
int BamAlignment::getMateRefwidth() const {
   pair<string,bool> mcval = getTag<string>("MC");
   if (!mcval.second) {
      cerr << *this << endl << __FILE__ << ":" << __LINE__ 
         << ": WARN BamAlignment has no MC tag\n";
      throw runtime_error("BamAlignment has no MC tag");
   }
   vector<pair<char,int>> mateCigar = parseCigar(mcval.first);
   int w =0;
   for (auto& x : mateCigar) {
      if (x.first == 'M' || x.first == 'D') {
         w += x.second;
      }
   }
   return w;
}

std::pair<int,int> BamAlignment::getPairedInterval() const {
   if (!mateOnSameReference() || getInsertSize() == 0) {
      return getInterval();
   }
   /*
   pair<int,int> tmp = getPairedRange();
   if (tmp.first > tmp.second) {
      return make_pair(tmp.second, tmp.first);
   }
   return tmp;
   */
   int b = getPosition();
   int b2 = getMatePosition();
   //if (getQueryName() == "S1464448") {
   //   cerr << "b=" << b << " b2=" << b2 << " templen=" << getInsertSize() << endl;
   //}
   if (isForwardStrand()) { // --R-->
      assert(getInsertSize() > 0);
      if (isMateForwardStrand()) { // both this and mate are forward direction
         // now figure out which one is on the left (smaller)
         //if (b < b2) { // --R--> --M-->
         //  // this is one the left, mate on right
         //  return make_pair(b, b2 + getMateRefwidth() - 1);
         //}
         //else { // --M-->  --R-->
         //   return make_pair(b2, getEndPosition());
         //}
         return make_pair(min(b,b2), max(getEndPosition(), b2 + getMateRefwidth() - 1));
      }
      else { // mate on Reverse strand
         //if (b <= b2) { // --R--> <--M--, properly mapped case
         //   // special case ==R==>
         //   //              <=M====
         //   return make_pair(b, b+getInsertSize()-1);
         //}
         //else { // <--M-- --R--> improper head-to-head case
         //   //return make_pair(b2, getEndPosition());
         //   return make_pair(b2, b2+getInsertSize()-1);
         //}
         int minB=min(b,b2);
         int maxE = minB + getInsertSize() - 1;
         if (getEndPosition() > maxE) {
            //cerr << *this << endl;
            //cerr << __FILE__ << ":" << __LINE__ << ":WARN: Warning MC tag and InsertSize inconsistent consider fix bug somewhere\n";
            maxE = getEndPosition();
         }
         if (getMatePosition() + getMateRefwidth() - 1 > maxE)
            maxE = getMatePosition() + getMateRefwidth() - 1;
         return make_pair(minB, maxE);
      }
   }
   else { // <--R-- this on reverse direction
      if (isMateForwardStrand()) { // <--R-- --M-->
         assert(getInsertSize() < 0);
         /*
         if (b < b2) { // <--R-- --M--> this is left
            // <==R==
            //   =M===>
            //return make_pair(b, b2 + getMateRefwidth() - 1);
            return make_pair(b, b-getInsertSize()-1);
         }
         else { // --M--> <--R-- Proper pair
            //return make_pair(b2, getEndPosition());
            return make_pair(b2, b2-getInsertSize()-1);
         }
         */
         int minB=min(b,b2);
         if (getEndPosition() + getInsertSize() +1 < minB)
            minB = getEndPosition() + getInsertSize() +1;
         int maxE = minB - getInsertSize() - 1;
         if (getEndPosition() > maxE) {
            //cerr << *this << endl;
            //cerr << __FILE__ << ":" << __LINE__ << ":WARN: Warning MC tag and InsertSize inconsistent consider fix bug somewhere\n";
            maxE = getEndPosition();
         }
         if (minB + getMateRefwidth() - 1 > maxE)
            maxE = minB + getMateRefwidth() - 1;
         return make_pair(minB, maxE);
      }
      else { //both this and mate on reverse strand
         //if (b < b2) { // <--R-- <--M--
         //   assert(getInsertSize() > 0);
         //   return make_pair(b, b2+getInsertSize() -1);
         //}
         //else { // <--M-- <--R--  this on the right hand side
         //   return make_pair(b2, getEndPosition());
         //}
         return make_pair(min(b,b2), max(getEndPosition(), b2 + getMateRefwidth() - 1));
      }
   }
}

int BamAlignment::getPairedEndPosition() const {
   if (!mateOnSameReference()) {
      return getEndPosition();
   }
   int b;
   if (IsReverseStrand()) {
      b=getMatePosition();
      return GetEndPosition(false,true);
   }
   b=getPosition();
   return b+getInsertSize()-1;
}
      
// [b, e] on the query sequence
BamAlignment BamAlignment::subsequence(int b, int e) const {
   int len = e - b +1;
   BamAlignment tmp(*this);
   tmp.setQueryLength(len);
   tmp.QueryBases = tmp.QueryBases.substr(b, len);
   tmp.Qualities = tmp.Qualities.substr(b,len);
   tmp.AlignedBases.clear();
   // Position
   // AlignmentFlag needs to be modified
   // InsertSize, CigarData
   int gi=Position, i=0;
   char cigarState='M';
   unsigned int cigarIdx=0, ci=0;

   while (i < b) {
      if (cigarIdx < CigarData[ci].Length) { // in one cigar segment
         if (cigarState == 'M' && cigarIdx < CigarData[ci].Length) {
            ++i; ++gi; ++cigarIdx;
         } // i at b
         else if (cigarState == 'D') {
            ++gi; ++cigarIdx;
         }
         else if (cigarState == 'I' || cigarState == 'S') {
            ++i; ++cigarIdx;
         }
         else {
            cerr << "wrong cigarop: " << cigarState 
               << __FILE__ << ":" << __LINE__ << endl;
            exit(1);
         }
      }
      else { // next cigar segment
         ++ci;
         if (ci >= CigarData.size()) {
            cerr << __FILE__ << ":" << __LINE__ 
               << " walked off the cigar string\n";
            exit(1);
         }
         char newState = CigarData[ci].Type;
         if ((cigarState == 'I' && newState == 'D')
               || (cigarState == 'D' && newState == 'I')) {
            cerr << "I/D transition in cigarop not permitted\n";
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << endl;
            exit(1);
         }
         cigarIdx = 0;
         cigarState = newState;
      }
   }
   int cigarIdx_b = cigarIdx;
   tmp.Position = gi; // new Position on genomic DNA
   vector<pair<char,int> > newcigarOp;
   // first cigar segment from usually a match_segment
   // if insert state then softclip
   if (cigarState == 'I' || cigarState == 'D') {
      cerr << __func__ << " Cannot stop inside an indel state!\n";
      exit(1);
   }
   while (i < e) {
      if (cigarIdx < CigarData[ci].Length) { // in one cigar segment
         if (cigarState == 'M' && cigarIdx < CigarData[ci].Length) {
            ++i; ++gi; ++cigarIdx;
         } // i at b
         else if (cigarState == 'D') {
            ++gi; ++cigarIdx;
         }
         else if (cigarState == 'I' || cigarState == 'S') {
            ++i; ++cigarIdx;
         }
         else {
            cerr << "wrong cigarop: " << cigarState 
               << __FILE__ << ":" << __LINE__ << endl;
            exit(1);
         }
      }
      else { // next cigar segment
         ++ci;
         if (ci >= CigarData.size()) {
            cerr << __FILE__ << ":" << __LINE__ 
               << " walked off the cigar string\n";
            exit(1);
         }
         newcigarOp.push_back(make_pair(cigarState, cigarIdx - cigarIdx_b));
         cigarIdx_b=0; // in most cases
         char newState = CigarData[ci].Type;
         if ((cigarState == 'I' && newState == 'D')
               || (cigarState == 'D' && newState == 'I')) {
            cerr << "I/D transition in cigarop not permitted\n";
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << endl;
            exit(1);
         }
         cigarIdx = 0;
         cigarState = newState;
      }
   }
   // last cigar segment
   newcigarOp.push_back(make_pair(cigarState, cigarIdx - cigarIdx_b));
   tmp.setCigarOperation(newcigarOp);
   return tmp;
} 

// if x is after position then we need to advance both
// i and j to the starting point
// must not start with S state, should be in the middle
// passed the Beginning state
/* tried to update but failed
void BamAlignment::advanceIndex(int &i, int &j, int &x, unsigned int &cigarIdx, unsigned int &ci) const
{
   char oldCigarState;
   while (i < x) {
      if (cigarIdx < CigarData[ci].Length) { // index in one cigar segment [0, length-1]
         int diff = x-i;
         cigarIdx += diff;
         if (CigarData[ci].Type == 'M') {
            i = x; j += diff;
         } 
         else if (CigarData[ci].Type == 'D') {
            i = x;
         }
         else if (CigarData[ci].Type == 'I') {
            j += diff;
         }
         else {
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__)
                  + ":ERROR wrong cigarop: " + string(1,CigarData[ci].Type));
         }
      }
      else { // next cigar segment
         oldCigarState = CigarData[ci].Type;
         ++ci; // advance cigar segment
         if (ci >= CigarData.size()) {
            throw out_of_range(string( __FILE__) + ":" + to_string(__LINE__)
               + ":ERROR walked off the cigar");
         }
         if ((CigarData[ci].Type == 'I' && oldCigarState == 'D')
               || (CigarData[ci].Type == 'D' && oldCigarState == 'I')) {
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__)
                  + ":ERROR Cigar I/D transition not permitted");
         }
         cigarIdx = 0;
      }
   }
   // ended at a new cigar segment
   if (cigarIdx == CigarData[ci].Length) {
      oldCigarState = CigarData[ci].Type;
      ++ci;
      cigarIdx = 0;
      if ((CigarData[ci].Type == 'I' && oldCigarState == 'D')
            || (CigarData[ci].Type == 'D' && oldCigarState == 'I')) {
         throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__)
                  + ":ERROR Cigar I/D transition not permitted");
      }
   }
   if (CigarData[ci].Type == 'I') { // skip I state to M
      j += CigarData[ci].Length;
      //while (cigarIdx < CigarData[ci].Length) {
      //   ++cigarIdx;
      //   ++j;
      //}
      ++ci;
      if (CigarData[ci].Type != 'M') {
         throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__)
               + ":ERROR M not follow I state!");
      }
      cigarIdx = 0;
   }
   // first cigar segment from usually a match_segment
   // if insert state then softclip
   if (CigarData[ci].Type == 'D') { // it is possible that the consensus is D
      // but one of the read has a base, we will skip the D and 
      // make the subsequence shorter
      //cerr << *this << endl;
      cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
         << " WARN: stop inside a D state! We will advance x to the next M\n";
      //cerr << "Subsequence is shortened because consensus favors deletion\n";
      //throw BamAlignmentException(string(__FILE__) + to_string(__LINE__)
      //      + string(__func__) + " begin index start inside D");
      //while (cigarIdx < CigarData[ci].Length) {
      //   ++cigarIdx; ++i; ++x;
      //}
      i += CigarData[i].Length;
      x += CigarData[i].Length;
      ++ci;
      cigarIdx = 0;
      if (CigarData[ci].Type == 'I') {
         throw BamAlignmentException(string(__FILE__) + to_string(__LINE__) 
               + string(__func__) + " ERROR: D/I transistion");
      }
   }
}
*/

// old implementation still works
void BamAlignment::advanceIndex(int &i, int &j, int &b, 
      unsigned int &cigarIdx, unsigned int &ci, char &cigarState) const
{
   // move i to position b if i is not at b
   // j to the start of the query index
   while (i < b) {
      //cout << "ci: " << ci << " cigarIdx: " << cigarIdx << " i: " << i 
      //   << " j: " << j << " cigarState: " << cigarState << endl;
      if (cigarIdx < CigarData[ci].Length) { // in one cigar segment
         //cout << cigarIdx << " < " << CigarData[ci].Length << endl;
         if (cigarState == 'M') {
            ++i; ++j; ++cigarIdx;
         } // i at b
         else if (cigarState == 'D') {
            ++i; ++cigarIdx;
         }
         else if (cigarState == 'I') {
            ++j; ++cigarIdx;
         }
         else {
            cerr << "wrong cigarop: " << cigarState 
               << __FILE__ << ":" << __LINE__ << endl;
            exit(1);
         }
      }
      else { // next cigar segment
         //cout << cigarIdx << " == " << CigarData[ci].Length << endl;
         ++ci;
         if (ci >= CigarData.size()) {
            cerr << __FILE__ << ":" << __LINE__ 
               << " walked off the cigar string\n";
            exit(1);
         }
         char newState = CigarData[ci].Type;
         if ((cigarState == 'I' && newState == 'D')
               || (cigarState == 'D' && newState == 'I')) {
            cerr << "I/D transition in cigarop not permitted\n";
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << endl;
            exit(1);
         }
         cigarIdx = 0;
         cigarState = newState;
      }
   }
   // ended at a new cigar segment
   if (cigarIdx == CigarData[ci].Length) { // need one more transition
      ++ci;
      cigarIdx = 0;
      char newState = CigarData[ci].Type;
      if ((cigarState == 'I' && newState == 'D')
            || (cigarState == 'D' && newState == 'I')) {
         cerr << "I/D transition in cigarop not permitted\n";
         cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
            << endl;
         exit(1);
      }
      cigarState = newState;
   }
   if (cigarState == 'I') { // skip I state to M
      //cout << "landed in an query insert state\n"
      //   << "i: " << i << " j: " << j << endl;
      while (cigarIdx < CigarData[ci].Length) {
         ++cigarIdx;
         ++j;
      }
      ++ci;
      char newState = CigarData[ci].Type;
      if (newState != 'M') {
         cerr << "M must follow I state!\n";
         exit(1);
      }
      cigarState=newState;
      cigarIdx = 0;
      //cout << "i: " << i << " j: " << j << endl;
   }
   // first cigar segment from usually a match_segment
   // if insert state then softclip
   if (cigarState == 'D') { // it is possible that the consensus is D
      // but one of the read has a base, we will skip the D and 
      // make the subsequence shorter
      //cerr << *this << endl;
      //cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
      //   << " WARN: stop inside a D state!\n";
      //cerr << "Subsequence is shortened because consensus favors deletion\n";
      //throw BamAlignmentException(string(__FILE__) + to_string(__LINE__)
      //      + string(__func__) + " begin index start inside D");
      while (cigarIdx < CigarData[ci].Length) {
         ++cigarIdx; ++i; ++b;
      }
      ++ci;
      cigarIdx = 0;
      char newState = CigarData[ci].Type;
      if (newState == 'I') {
         cerr << "D/I transition in cigarop not permitted\n";
         cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
            << endl;
         throw BamAlignmentException(string(__FILE__) + to_string(__LINE__) 
               + string(__func__) + " ERROR: D/I transistion");
         //exit(1);
      }
      cigarState = newState;
   }
}

void BamAlignment::nextCigar(int& i, int& j, unsigned int& ci) const {
   if  (getCigarType(ci) == 'M') {
      i += getCigarLength(ci);
      j += getCigarLength(ci);
   }
   else if (getCigarType(ci) == 'D') {
      i += getCigarLength(ci);
   }
   else if (getCigarType(ci) == 'I') {
      j += getCigarLength(ci);
   }
   else {
      throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) 
            + ":ERROR unexpected CIGAR Type: " + string(1, getCigarType(ci)));
   }
   ++ci;
}

// default starR=0
// for M stop at the start of the segment
// for I stop at the end of the Previous M segment
// for D stop at the start of D segment
pair<int,bool> BamAlignment::isInsertionAtRefloc(int desiredR, int startR) const {
   int q=0, c=0; // q is only used for checking condition
   int r=startR;
   while (c < (int)CigarData.size() && (CigarData[c].getType() == 'S' || CigarData[c].getType() == 'H')) {
      q += CigarData[c].getLength();
      ++c;
   }
   while (c < (int)CigarData.size() && r < desiredR && q < getLength() && CigarData[c].getType() != 'S') {
      if (CigarData[c].getType() == 'M') {
         r += CigarData[c].getLength();
         q += CigarData[c].getLength();
      }
      else if (CigarData[c].getType() == 'I') {
         if (r-1 == desiredR) { // insert attach to the last base of previous M segment
            return make_pair(c,true);
         }
         q += CigarData[c].getLength();
      }
      else if (CigarData[c].getType() == 'D') {
         r += CigarData[c].getLength();
      }
      else {
         cerr << "trying to find Insertion at " << desiredR << endl;
         throw runtime_error(string("wrong CigarData state: ") + string(1, CigarData[c].getType()));
      }
      ++c;
   }
   if (r-1 == desiredR && CigarData[c].getType() == 'I') {
      return make_pair(c, true);
   }
   //if (c >= CigarData.size()) c=CigarData.size()-1;
   return make_pair(c, false);
}

pair<int,bool> BamAlignment::isDeletionAtRefloc(int desiredR, int startR) const {
   int q=0, c=0; // don't need q yet, may need it later.
   int r=startR;
   while (c < (int)CigarData.size() && (CigarData[c].getType() == 'S' || CigarData[c].getType() == 'H')) {
      q += CigarData[c].getLength();
      ++c;
   }
   while (c < (int)CigarData.size() && r <= desiredR && q < getLength() && CigarData[c].getType() != 'S') {
      if (CigarData[c].getType() == 'M') {
         r += CigarData[c].getLength();
         q += CigarData[c].getLength();
      }
      else if (CigarData[c].getType() == 'I') {
         q += CigarData[c].getLength();
      }
      else if (CigarData[c].getType() == 'D') {
         if (r == desiredR) { // r is the first base of the deleted Refseq
            return make_pair(c, true);
         }
         r += CigarData[c].getLength();
         if (desiredR < r && desiredR > r-(int)getCigarLength(c)) {
            //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG " << desiredR
            //   << " is not the first base of the deletion\n";
            return make_pair(c,true);
         }
      }
      else {
         cerr << "trying to find Insertion at " << desiredR << endl;
         throw runtime_error(string("wrong CigarData state: ") + string(1, CigarData[c].getType()));
      }
      ++c;
   }
   //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG r=" << r << " q=" << q
   //   << " c=" << c << " no deletion at " << desiredR << endl;
   //if (c >= CigarData.size()) c=CigarData.size()-1;
   return make_pair(c, false);
}

// if i is in any position inside D segment, then
// j will be the first base index of the next M segment.
// if i is in any insert segment, the j will be
// the first base-index fn the insert
// insertion can only be addressed by the Base before or after
// on the reference.
int BamAlignment::indexRef2Query(int ri) const {
   int i=getPosition();
   assert(ri>= i && ri <= getEndPosition());
   int j=0; //i index on reference, j index on query
   unsigned int ci=0; // ci is index in the cigaroperation: CigarData
   if (getCigarType(0) == 'S') { // S => next state
      j += getCigarLength(0);
      ci=1;
   }
   while (i < ri && ci < getCigarOperationCount()) {
      if ((unsigned int)ri < i+getCigarLength(ci)) { // fall within this segment, done
         if  (getCigarType(ci) == 'M') {
            // this is expected segment
            return j + ri -i; // if i falls on the last base of M
            // and next segment is I, you known 
         }
         else if (getCigarType(ci) == 'D') {
            return j;
         }
         else if (getCigarType(ci) == 'I') { // should never reach this state
            //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG should not land inside Insertion\n";
            throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":DEBUG unreachable code section inside Insertion");
         }
         else {
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) 
                  + ":ERROR unexpected CIGAR Type: " + string(1, getCigarType(ci)));
         }
         //return QueryBases[j];
      }
      else { // advance to the next segment
         nextCigar(i,j,ci);
      }
   }
   throw logic_error(string(__FILE__) + ":" + to_string(__LINE__)
         + ":ERROR coding error, cannot fine char at " + to_string(ri));
}

// [b,e] are reference coordinate
/* new implementation not working yet
BamAlignment BamAlignment::subsequenceByRef(int b, int e) const {
   assert(b>=Position);
   int i=Position, j=0; //i on reference, j index on query
   unsigned int cigarIdx=0, ci=0; // cigarIdx is index within each cigar segment
   int subqseqBegin = 0; // begin index in query bases

   // skip initial softclip if start after the soft clip
   vector<pair<char,int> > newcigarOp;
   if (CigarData[ci].Type == 'S') { // S => next state
      if (b == Position) {
         newcigarOp.push_back(
               make_pair(CigarData[ci].Type, CigarData[ci].Length));
      }
      else {
         subqseqBegin = CigarData[ci].Length;
      }
      j += CigarData[ci].Length;
      cigarIdx = 0;
      ++ci;
   }
   if (i<b) {
      advanceIndex(i, j, b, cigarIdx, ci);
      advanceIndex(i, j, b, cigarIdx, ci);
      subqseqBegin=j;
   }
   int cigarIdx_b = cigarIdx;
   while (i < e) {
      if (cigarIdx < CigarData[ci].Length) { // in one cigar segment
         if (CigarData[ci].Type == 'M') {
            ++i; ++j; ++cigarIdx;
         } // i at b
         else if (CigarData[ci].Type == 'D') {
            ++i; ++cigarIdx;
         }
         else if (CigarData[ci].Type == 'I') {
            ++j; ++cigarIdx;
         }
         else {
            cerr << "wrong cigarop: " << CigarData[ci].Type 
               << __FILE__ << ":" << __LINE__ << endl;
            exit(1);
         }
      }
      else { // next cigar segment, end of last segment
         //cout << "Next cigar segment\n";
         char oldCigarState = CigarData[ci].Type;
         ++ci;
         if (ci >= CigarData.size()) {
            cerr << __FILE__ << ":" << __LINE__ 
               << " walked off the cigar string: ci=" << ci 
               << endl;
            exit(1);
         }
         newcigarOp.push_back(make_pair(oldCigarState, cigarIdx - cigarIdx_b));
         cigarIdx_b=0; // in most cases
         if ((CigarData[ci].Type == 'I' && oldCigarState == 'D')
               || (CigarData[ci].Type == 'D' && oldCigarState == 'I')) {
            cerr << "I/D transition in cigarop not permitted\n";
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << endl;
            exit(1);
         }
         cigarIdx = 0;
         //cigarState = newState;
      }
   }
   //cout << "cigarIdx_b " << cigarIdx_b << endl;
   // last cigar segment
   newcigarOp.push_back(make_pair(CigarData[ci].Type, cigarIdx - cigarIdx_b + 1));
   // if got softclip, we need to add it
   if (cigarIdx == CigarData[ci].Length && ci+1 < CigarData.size()
         && CigarData[ci+1].Type == 'S') 
   {
      newcigarOp.push_back(make_pair(CigarData[ci+1].Type, CigarData[ci+1].Length));
      // query index needs to be push further
      j += CigarData[ci+1].Length;
   }

   BamAlignment tmp(*this);
   tmp.Position = b; // new Position on genomic DNA
   tmp.QueryBases = tmp.QueryBases.substr(subqseqBegin, j-subqseqBegin+1);
   tmp.Qualities = tmp.Qualities.substr(subqseqBegin, j-subqseqBegin+1);
   tmp.setQueryLength(tmp.QueryBases.length());
   tmp.AlignedBases.clear();
   tmp.setCigarOperation(newcigarOp);
   return tmp;
} 

// TODO: need further testing and verification
std::string BamAlignment::substringByRef(int b, int e) const {
   assert(b>=Position);
   // InsertSize, CigarData
   int i=Position, j=0; //i on reference, j index on query
   //char cigarState = 'M';
   unsigned int cigarIdx=0, ci=0; // cigarIdx is index within each cigar segment
   // ci is the cigar segment index
   int subqseqBegin = 0; // begin index in query bases
   // skip initial softclip if start after the soft clip
   if (CigarData[ci].Type == 'S') { // S => next state
      subqseqBegin = CigarData[ci].Length;
      j += CigarData[ci].Length;
      ++ci;
   }
   if (i<b) {
      advanceIndex(i, j, b, cigarIdx, ci);
      subqseqBegin=j;
   }
   // bring i to e on reference, j follows on query
   while (i < e) {
      if (cigarIdx < CigarData[ci].Length) { // in this cigar segment
         if (CigarData[ci].Type == 'M') {
            ++i; ++j; ++cigarIdx;
         } // i at b
         else if (CigarData[ci].Type == 'D') {
            ++i; ++cigarIdx;
         }
         else if (CigarData[ci].Type == 'I') {
            ++j; ++cigarIdx;
         }
         else {
            cerr << "wrong cigarop: " << CigarData[ci].Type 
               << __FILE__ << ":" << __LINE__ << endl;
            throw runtime_error("while obtaining subseq unknown cigar state");
         }
      }
      else { // next cigar segment, end of last segment
         //cout << "Next cigar segment\n";
         char oldCigarState = CigarData[ci].Type;
         ++ci;
         if (ci >= CigarData.size()) {
            cerr << __FILE__ << ":" << __LINE__ 
               << " walked off the cigar string: ci=" << ci 
               << endl;
            exit(1);
         }
         if ((CigarData[ci].Type == 'I' && oldCigarState == 'D')
               || (CigarData[ci].Type == 'D' && oldCigarState == 'I')) {
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << ":WARN I/D or D/I transition in cigarop need more coding.\n";
            throw runtime_error("Cigar I|D or D|I transition");
         }
         cigarIdx = 0;
      }
   }
   if (cigarIdx == CigarData[ci].Length && ci+1 < CigarData.size()
         && CigarData[ci+1].Type == 'S') 
   {
      // query index needs to be push further
      j += CigarData[ci+1].Length;
   }
   //cout << "query sub: " << subqseqBegin << "-" << j+1 << endl;
   return QueryBases.substr(subqseqBegin, j-subqseqBegin+1);
}
*/

BamAlignment BamAlignment::subsequenceByRef(int b, int e) const {
   assert(b>=Position);
   // AlignmentFlag needs to be modified
   // InsertSize, CigarData
   int i=Position, j=0; //i on reference, j index on query
   char cigarState = 'M';
   unsigned int cigarIdx=0, ci=0; // cigarIdx is index within each cigar segment
   int subqseqBegin = 0; // begin index in query bases
   //cout << "Taking subsequence of BamAlignment by ref index:\n"
   //   << " [" << b << "," << e << "]: " 
   //   << e-b+1 << " i: " << i << "\n";
   //cout << "parent sequence\n" << *this << endl;

   // skip initial softclip if start after the soft clip
   vector<pair<char,int> > newcigarOp;
   if (CigarData[ci].Type == 'S') { // S => next state
      if (b == Position) {
         newcigarOp.push_back(
               make_pair(CigarData[ci].Type, CigarData[ci].Length));
      }
      else {
         subqseqBegin = CigarData[ci].Length;
      }
      j += CigarData[ci].Length;
      cigarIdx = 0;
      ++ci;
   }
   if (i<b) {
      //cout << " i(" << i << ") less than b(" << b << ") advaicing i,j to\n";
      advanceIndex(i, j, b, cigarIdx, ci, cigarState);
      subqseqBegin=j;
      //cout << "i,j" << i << "," << j << " cigarIdx: " << cigarIdx 
      //   << " b " << b << endl;
   }
   //else {
   //   cout << "Not need to advance: i,j" << i << "," << j << " cigarIdx: " << cigarIdx 
   //      << " b " << b << endl;
   //}
   int cigarIdx_b = cigarIdx;
   // bring i to e
   //cout << "Now i: " << i << " should be at b: " << b << "\n";
   while (i < e) {
      //cout << "ci: " << ci << " cigarIdx: " << cigarIdx << " i: " << i 
      //      << " j: " << j << " cigarState: " << cigarState << endl;
      if (cigarIdx < CigarData[ci].Length) { // in one cigar segment
         if (cigarState == 'M') {
            ++i; ++j; ++cigarIdx;
         } // i at b
         else if (cigarState == 'D') {
            ++i; ++cigarIdx;
         }
         else if (cigarState == 'I') {
            ++j; ++cigarIdx;
         }
         else {
            cerr << "wrong cigarop: " << cigarState 
               << __FILE__ << ":" << __LINE__ << endl;
            exit(1);
         }
      }
      else { // next cigar segment, end of last segment
         //cout << "Next cigar segment\n";
         ++ci;
         if (ci >= CigarData.size()) {
            cerr << __FILE__ << ":" << __LINE__ 
               << " walked off the cigar string: ci=" << ci 
               << endl;
            exit(1);
         }
         newcigarOp.push_back(make_pair(cigarState, cigarIdx - cigarIdx_b));
         cigarIdx_b=0; // in most cases
         char newState = CigarData[ci].Type;
         if ((cigarState == 'I' && newState == 'D')
               || (cigarState == 'D' && newState == 'I')) {
            cerr << "I/D transition in cigarop not permitted\n";
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << endl;
            exit(1);
         }
         cigarIdx = 0;
         cigarState = newState;
      }
   }
   //cout << "cigarIdx_b " << cigarIdx_b << endl;
   // last cigar segment
   newcigarOp.push_back(make_pair(cigarState, cigarIdx - cigarIdx_b + 1));
   // if got softclip, we need to add it
   if (cigarIdx == CigarData[ci].Length && ci+1 < CigarData.size()
         && CigarData[ci+1].Type == 'S') 
   {
      newcigarOp.push_back(make_pair(CigarData[ci+1].Type, CigarData[ci+1].Length));
      // query index needs to be push further
      j += CigarData[ci+1].Length;
   }
   //cout << "At end: i,j" << i << "," << j << " cigarIdx: " << cigarIdx 
   //   << " b " << b << " subqseqBegin: " << subqseqBegin << endl;

   BamAlignment tmp(*this);
   tmp.Position = b; // new Position on genomic DNA
   tmp.QueryBases = tmp.QueryBases.substr(subqseqBegin, j-subqseqBegin+1);
   tmp.Qualities = tmp.Qualities.substr(subqseqBegin, j-subqseqBegin+1);
   tmp.setQueryLength(tmp.QueryBases.length());
   tmp.AlignedBases.clear();
   tmp.setCigarOperation(newcigarOp);
   //cout << "parent query bases:\n"
   //   << getQueryBases() << "\nnew query bases\n"
   //   << tmp.getQueryBases() << endl
   //   << "new cigarstring: " << tmp.getCigarString() << endl;
   return tmp;
} 

std::string BamAlignment::substringByRef(int b, int e) const {
   assert(b>=Position);
   // InsertSize, CigarData
   int i=Position, j=0; //i on reference, j index on query
   char cigarState = 'M';
   unsigned int cigarIdx=0, ci=0; // cigarIdx is index within each cigar segment
   // ci is the cigar segment index
   int subqseqBegin = 0; // begin index in query bases
   // skip initial softclip if start after the soft clip
   if (CigarData[ci].Type == 'S') { // S => next state
      subqseqBegin = CigarData[ci].Length;
      j += CigarData[ci].Length;
      ++ci;
   }
   if (i<b) {
      advanceIndex(i, j, b, cigarIdx, ci, cigarState);
      subqseqBegin=j;
   }
   // bring i to e on reference, j follows on query
   while (i < e) {
      if (cigarIdx < CigarData[ci].Length) { // in this cigar segment
         if (cigarState == 'M') {
            ++i; ++j; ++cigarIdx;
         } // i at b
         else if (cigarState == 'D') {
            ++i; ++cigarIdx;
         }
         else if (cigarState == 'I') {
            ++j; ++cigarIdx;
         }
         else {
            cerr << "wrong cigarop: " << cigarState 
               << __FILE__ << ":" << __LINE__ << endl;
            throw runtime_error("while obtaining subseq unknown cigar state");
         }
      }
      else { // next cigar segment, end of last segment
         //cout << "Next cigar segment\n";
         ++ci;
         if (ci >= CigarData.size()) {
            cerr << __FILE__ << ":" << __LINE__ 
               << " walked off the cigar string: ci=" << ci 
               << endl;
            exit(1);
         }
         char newState = CigarData[ci].Type;
         if ((cigarState == 'I' && newState == 'D')
               || (cigarState == 'D' && newState == 'I')) {
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
               << ":WARN I/D or D/I transition in cigarop need more coding.\n";
            throw runtime_error("Cigar I|D or D|I transition");
            //exit(1);
         }
         cigarIdx = 0;
         cigarState = newState;
      }
   }
   if (cigarIdx == CigarData[ci].Length && ci+1 < CigarData.size()
         && CigarData[ci+1].Type == 'S') 
   {
      // query index needs to be push further
      j += CigarData[ci+1].Length;
   }
   //cout << "query sub: " << subqseqBegin << "-" << j+1 << endl;
   return QueryBases.substr(subqseqBegin, j-subqseqBegin+1);
}

char BamAlignment::charAtByRef(int ri) const {
   int i=getPosition();
   assert(ri>=i && ri <= getEndPosition());
   int j=0; //i index on reference, j index on query
   if (i == ri) { // special case where first base requested
      return QueryBases[j];
   }
   unsigned int ci=0; // ci is index in the cigaroperation: CigarData
   if (getCigarType(0) == 'S') { // S => next state
      j += getCigarLength(0);
      ci=1;
   }
   while (i <= ri && ci < getCigarOperationCount()) {
      if (i < ri && ci < getCigarOperationCount() && getCigarType(ci) == 'I') {
         nextCigar(i,j,ci);
      }
      else if ((unsigned int)ri < i+getCigarLength(ci)) { // fall within this segment, done
         if  (getCigarType(ci) == 'M') {
            // this is expected segment
            return QueryBases[j + ri - i];
         }
         else if (getCigarType(ci) == 'D') {
            return '-';
            //throw logic_error("char at " + to_string(ri) + " is inside deletion");
         }
         else if (getCigarType(ci) == 'I') { // unnecessary check
            // 177M2I60M p=109398, ri=190575 falls inside INSERT
            // need to advance to next M segment.
            nextCigar(i,j,ci);
            // the nextCigar() method will skip this one
            // TODO: replace with assert()
            //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG should not land inside Insertion\n";
            //cerr << *this << endl;
            //cerr << "ri=" << ri << " i=" << i << " j=" << j << " ci=" << ci << endl;
            //throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":DEBUG should not land inside Insertion");
         }
         else {
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) 
                  + ":ERROR unexpected CIGAR Type: " + string(1, getCigarType(ci)));
         }
      }
      else { // advance to the next segment
         nextCigar(i,j,ci);
      }
   }
   cerr << "i=" << i << " j=" << j << " ci=" << ci << endl;
   cerr << *this << endl;
   throw logic_error(string(__FILE__) + ":" + to_string(__LINE__)
         + ":ERROR coding error, cannot find char at " + to_string(ri));
}

//check whether this alignment has a query deletion
//of length len at index i
//TODO: convert exceptions to assert after debug
bool BamAlignment::isDeletionAt(int ri, int len) const {
   if (lackDCigar()) return false;
   int i = getPosition();
   assert(ri>=i && ri <= getEndPosition());
   int j=0; //i index on reference, j index on query
   unsigned int ci=0; // ci is index in the cigaroperation: CigarData
   if (getCigarType(0) == 'S') { // S => next state
      j += getCigarLength(0);
      ci=1;
   }
   while (i <= ri && ci < getCigarOperationCount()) {
      //cout << "i=" << i << " j=" << j << " ci=" << ci << endl;
      if (i < ri && ci < getCigarOperationCount() && getCigarType(ci) == 'I') {
         nextCigar(i,j,ci);
      }
      else if ((unsigned int)ri < i+getCigarLength(ci)) { // fall within this segment, done
         //cout << "ri=" << ri << " less than i+getCigarLength(ci)=" << i+getCigarLength(ci) << endl;
         //cout << "cigar type: " << getCigarType(ci) << endl;
         if  (getCigarType(ci) == 'M') {
            //cout << __FILE__ << ":" << __LINE__ << ":DEBUG Looking for deletion but fall inside Match\n";
            return false;
         }
         else if (getCigarType(ci) == 'D') {
            if (i == ri) { // must fall on the first base of the DEL segment
               if (getCigarLength(ci) == (unsigned int)len) {
                  return true;
               }
               else { // is deletion of different length
                  //cout << __FILE__ << ":" << __LINE__ << ":DEBUG Looking for deletion of length="
                  //   << len << " but saw DEL length=" << getCigarLength(ci) << endl;
                  return false;
               }
            }
            else { // not on first base of Deletion
               //cout << *this << endl;
               //cout << "ri=" << ri << " len=" << len << " i=" << i 
               //   << " j=" << j << " ci=" << ci << endl;
               //throw logic_error(string(__FILE__) + ":" + to_string(__LINE__)
               //      + ":DEBUG index " + to_string(ri) + " does not fall on the first base of deleted region");
               return false;
            }
         }
         else if (getCigarType(ci) == 'I') {
            //cout << *this << endl;
            //cout << "ri=" << ri << " len=" << len << " i=" << i 
            //   << " j=" << j << " ci=" << ci << endl;
            //throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":DEBUG should not land inside Insertion");
            return false;
         }
         else {
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) 
                  + ":ERROR unexpected CIGAR Type: " + string(1, getCigarType(ci)));
         }
      }
      else { // advance to the next segment
         //cout << "advance cigar\n";
         nextCigar(i,j,ci);
      }
   }
   throw logic_error(string(__FILE__) + ":" + to_string(__LINE__)
         + ":ERROR coding error, cannot fine char at " + to_string(ri));
}

// the insertion sequence is usually the lastBase of M + insert_sequence
bool BamAlignment::isInsertionAt(int ri, const string& seq) const {
   if (lackICigar()) return false;
   int i = getPosition();
   assert(ri>=i && ri <= getEndPosition());
   int j=0; //i index on reference, j index on query
   unsigned int ci=0; // ci is index in the cigaroperation: CigarData
   if (getCigarType(0) == 'S') { // S => next state
      j += getCigarLength(0);
      ci=1;
   }
   while (i <= ri && ci < getCigarOperationCount()) {
      if (i < ri && ci < getCigarOperationCount() && getCigarType(ci) == 'I') {
         nextCigar(i,j,ci);
      }
      else if ((unsigned int)ri < i+getCigarLength(ci)) { // fall within this segment, done
         if (getCigarType(ci) == 'M') {
            // this is expected segment, ri should be at last base
            // of M exactly
            if (i+getCigarLength(ci)-1 == (unsigned int)ri) { // last base of M
               assert(ci+1 < getCigarOperationCount());
               if (getCigarLength(ci+1)+1 == seq.size()) { // same insert size
                  if (QueryBases.substr(j+getCigarLength(ci)-1, seq.size()) == seq) {
                     return true;
                  }
                  else {
                     //cout << __FILE__ << ":" << __LINE__ << ":DEBUG saw insertion of same length but different sequence: "
                     //   << QueryBases.substr(j, seq.size()) << "\n";
                     return false;
                  }
               }
               else {
                  //cerr << *this << endl
                  //   << "ri=" << ri << " i=" << i << " j=" << j
                  //   << " ci=" << ci << endl;
                  //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG saw insertion of different length: "
                  //   << QueryBases.substr(j, seq.size()) << "\n";
                  return false;
               }
            }
            else {
               //cerr << *this << endl
               //   << "ri=" << ri << " i=" << i << " j=" << j
               //   << " ci=" << ci << endl;
               return false;
               //throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + to_string(ri) + " does not fall at the base before the insert segment");
            }
         }
         else if (getCigarType(ci) == 'D') { // inside D
            //cerr << *this << endl
            //   << "ri=" << ri << " i=" << i << " j=" << j
            //   << " ci=" << ci << endl;
            //throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Insert location " + to_string(ri) + " is inside deletion");
            return false;
         }
         else if (getCigarType(ci) == 'I') { // insertion location should always be
            // the last base of M, not inside INS
            //cerr << *this << endl
            //   << "ri=" << ri << " i=" << i << " j=" << j
            //   << " ci=" << ci << endl;
            //cerr << __FILE__ << ":" << __LINE__ << ":DEBUG should not land inside Insertion\n";
            //throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":DEBUG should not land inside Insertion");
            return false;
         }
         else {
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__) 
                  + ":ERROR unexpected CIGAR Type: " + string(1, getCigarType(ci)));
         }
         //return QueryBases[j];
      }
      else { // advance to the next segment
         nextCigar(i,j,ci);
      }
   }
   throw logic_error(string(__FILE__) + ":" + to_string(__LINE__)
         + ":ERROR coding error, cannot fine char at " + to_string(ri));
}

// due to redundancy, extra work is needed
// bad design
void BamAlignment::chopFirstSoftclip() {
   if (CigarData.front().Type == 'S') {
      int tmplen = CigarData.front().getLength();
      QueryBases=QueryBases.substr(tmplen);
      SupportData.QuerySequenceLength -= tmplen;
      Qualities=Qualities.substr(tmplen);
      SupportData.NumCigarOperations = CigarData.size();
      CigarData.erase(CigarData.begin());
   }
}

void BamAlignment::chopLastSoftclip() {
   //assert(CigarData.back().getType() == 'S');
   if (CigarData.back().getType() == 'S') {
      int tmplen = CigarData.back().getLength();
      SupportData.QuerySequenceLength -= tmplen;
      QueryBases.resize(SupportData.QuerySequenceLength);
      Qualities.resize(SupportData.QuerySequenceLength);
      SupportData.NumCigarOperations = CigarData.size();
      CigarData.resize(CigarData.size()-1);
   }
}
// same operaton as chopFirstSoftclip except will make sure
// the soft clip is dangling off the reference.
void BamAlignment::chopDangleFrontSoft() {
   assert(CigarData.front().getType() == 'S' && getPosition() == 0);
   SupportData.QuerySequenceLength -= CigarData.front().getLength();
   QueryBases = QueryBases.substr(CigarData.front().getLength());
   Qualities = Qualities.substr(CigarData.front().getLength());
   --SupportData.NumCigarOperations;
   SupportData.QuerySequenceLength = QueryBases.size(); 
   // update cigar last
   CigarData.erase(CigarData.begin());
}

void BamAlignment::chopDangleBackSoft() {
   //pair<string,int> nl = getRefnameFromId(getReferenceId());
   int reflen = getReferenceLength();
   assert(CigarData.back().getType() == 'S' && 
         static_cast<int>(getEndPosition() + CigarData.back().getLength()) >= reflen);
   SupportData.QuerySequenceLength -= CigarData.back().getLength();
   QueryBases.resize(QueryBases.size() - CigarData.back().getLength());
   Qualities.resize(QueryBases.size());
   --SupportData.NumCigarOperations;
   SupportData.QuerySequenceLength = QueryBases.size(); 
   // update cigar last
   CigarData.resize(CigarData.size()-1);
}

bool BamAlignment::validMD() const {
   string tmp = getStringTag("MD");
   auto i = tmp.find('^');
   if (i == 0) return false;
   while (i != string::npos) {
      if (i >= tmp.size()-1) return false;
      if (isalpha(tmp[i-1]) || !isdigit(tmp[i-1]))
         return false;
      i = tmp.find('^', i+1);
   }
   return true;
}

// the read should not be longer than 255 nt
// for longer reads we need to use integer as type
// the first character being ^ means deletion
pair<vector<int>, vector<string>> BamAlignment::getMDArray() const {
   vector<int> mdseg;
   vector<string> mdref;
   auto[md, hasmd] = getTag<string>("MD");
   if (hasmd) {
      // even position number, odd position letter last 3 digits
      // for bases 00 A, 01 C, 10 G, 11 T, 100 N, first bit del
      string::size_type i=0;
      while (i < md.size()) {
         auto b=i;
         while (i < md.size() && isdigit(md[i])) ++i;
         if (i == md.size()) {
            mdseg.push_back(std::stoi(md.substr(b)));
            break;
         }
         else {
            mdseg.push_back(std::stoi(md.substr(b, i-b)));
         }
         b = i;
         while (i < md.size() && !isdigit(md[i])) ++i;
         if (i == md.size()) {
            throw runtime_error("improper md tag: " + md);
         }
         mdref.push_back(md.substr(b, i-b));
      }
   }
   else {
      cerr << __FILE__ << ":" << __LINE__ << ":WARN no MD tag will return empty object\n";
   }
   for (const string& s : mdref) {
      if (s.front() == '^') {
         if (s.find('^', 1) != string::npos) {
            throw runtime_error(string(__func__) + ": Invalid MD tag: " + md);
         }
      }
      else {
         if (s.find('^') != string::npos) {
            throw runtime_error(string(__func__) + ": Invalid MD tag: " + md);
         }
      }
   }
   return make_pair(mdseg, mdref);
}

int BamAlignment::getMDWidth() const {
   int len=0;
   string mdval = getStringTag("MD");
   if (mdval.empty()) {
      if (isUnmapped()) {
         //cerr << "empty tag on unmapped is fine\n";
         return 0;
      }
      else {
         throw logic_error("MD tag is empty");
      }
   }
   size_t i=0; 
   while (i < mdval.size()) {
      size_t j=i+1;
      while (j < mdval.size() && isdigit(mdval[j])) ++j;
      try {
         len += stoi(mdval.substr(i, j-i));
      }
      catch (const exception& err) {
         cerr << " mdval=" << mdval << " i=" << i << " j-i=" << j-i
            << " " << mdval.substr(i, j-i) << endl;
         throw logic_error(string(__func__) + ": invalid MD " + mdval);
      }
      i = j;  // first base in seg
      if (i >= mdval.size()) break;
      if (mdval[i] == '^') ++i;
      j = i+1;
      while (j < mdval.size() && isalpha(mdval[j])) ++j;
      len += (j-i);
      i=j;
   }
   return len;
}

void BamAlignment::updateMDTag(const pair<vector<int>, vector<string>>& mdvec) {
   ostringstream oust;
   size_t i=0;
   while (i < mdvec.first.size()) {
      oust << mdvec.first[i];
      if (i < mdvec.second.size()) {
         oust << mdvec.second[i];
      }
      ++i;
   }
   EditTag("MD", "Z", oust.str());
}

// for debug
//void BamAlignment::recalMD(const string& refsq, mutex& mtx) {
void BamAlignment::recalMD(const string& refsq) {
   vector<int> matchPart;
   vector<string> mismatchPart;
   int matchcnt=0; int mismatchcnt=0;
   int ri = getPosition();
   int qi = 0;
   unsigned int c = 0;
   bool sawins=false;
   while (c < CigarData.size() && qi < static_cast<int>(QueryBases.size()) && ri <= getEndPosition()) {
      if (getCigarType(c) == 'S' || getCigarType(c) == 'H') {
         qi += getCigarLength(c);
      }
      else if (getCigarType(c) == 'I') { // will not record query insertion
         qi += getCigarLength(c);
         mismatchcnt += getCigarLength(c);
         sawins = true;
      }
      else if (getCigarType(c) == 'D') {
         if (matchPart.empty()) {
            throw logic_error("Cigar cannot start with D");
         }
         else if (matchPart.size() == mismatchPart.size())
            matchPart.push_back(0);
         mismatchcnt += getCigarLength(c);
         string tmpseq=refsq.substr(ri, getCigarLength(c));
         for (char& c : tmpseq) c = toupper(c);
         mismatchPart.push_back("^" + tmpseq);
         ri += getCigarLength(c);
      }
      else if (getCigarType(c) == 'M') {
         int Iend = qi + getCigarLength(c);
         while (qi < Iend) {
            char rch = toupper(refsq[ri]);
            if (rch != QueryBases[qi]) { // start differ string
               // need to add 0 to match part
               string diffseq(1, rch);
               ++qi; ++ri;
               while (qi < Iend && (rch=toupper(refsq[ri])) != QueryBases[qi]) { 
                  diffseq.append(1, rch);
                  ++qi; ++ri;
               }
               if (matchPart.size() == mismatchPart.size()) 
                  matchPart.push_back(0);
               mismatchcnt += diffseq.size();
               mismatchPart.push_back(std::move(diffseq));
            }
            else { // if (refseq[ri] == QueryBases[qi]) {
               int numSame=1;
               ++qi; ++ri;
               while (qi < Iend && toupper(refsq[ri]) == QueryBases[qi]) {
                  ++numSame;
                  ++qi; ++ri;
               }
               matchcnt += numSame;
               if (sawins) {
                  matchPart.back() += numSame;
                  sawins=false;
               }
               else {
                  matchPart.push_back(numSame);
               }
            }
         }
      }
      else {
         throw logic_error("un expected cigar type: " + string(1, getCigarType(c)));
      }
      ++c;
   }
   float iden = static_cast<float>(matchcnt)/(matchcnt+mismatchcnt);
   if (iden < 0.75) {
      cerr << __FILE__ << ":" << ": identity=" << iden << " too low after MD recalculation\n";
      throw logic_error("MD compuation need further testing");
   }
   if (matchPart.size() == mismatchPart.size()) {
      // need padding
      matchPart.push_back(0);
   }
   else if (matchPart.size() - mismatchPart.size() != 1) {
      //lock_guard<mutex> lg(mtx);
      cerr << matchPart.size() << " match segment\n";
      for (auto& m : matchPart)
         cerr << m << " ";
      cerr << endl;
      cerr << mismatchPart.size() << " mismatch segment\n";
      for (auto& m : mismatchPart) 
         cerr << m << " ";
      cerr << endl;
      // this is bad
      throw logic_error(string(__func__) + ": match should have one more element than mismatch");
   }
   updateMDTag(make_pair(matchPart, mismatchPart));
}


// len is the length to trim from query
void BamAlignment::chopFront(size_t len, int numMismatch) {
   // if InsertSize is zero do nothing
   if (numMismatch > 0) {
      uint8_t NMval = getNMValue();
      if (static_cast<int>(NMval) >= numMismatch) {
         NMval -= static_cast<uint8_t>(numMismatch);
         EditTag("NM", "C", NMval);
      }
      else {
         throw runtime_error("NM value update error");
      }
   }
   unsigned int alnchop=len;
   unsigned int querychop=len;
   unsigned int refadvance = len;
   // with and without softclip different part
   if (!startWithSoftclip()) { // no softclip
      if (CigarData.front().getType() == 'I') {
         alnchop += CigarData.front().getLength();
         querychop += CigarData.front().getLength();
         CigarData.erase(CigarData.begin());
      }
      if (CigarData.front().getType() != 'M') { // exception
         cerr << *this << endl;
         throw logic_error(getQueryName() + " expecting the front one either as I (removed) or M " + getCigarString());
      }
      //if (CigarData.front().getType() == 'M') { // such as 147M, 114M2I19M
      if (CigarData.front().getLength() > len) {
         CigarData.front().shrink(len);
      }
      else if (CigarData.front().getLength() == len) {
         // 2M6D145M needs special treatment
         if (CigarData[1].getType() == 'D') {
            // remove two cigardata
            alnchop += CigarData[1].getLength();
            refadvance += CigarData[1].getLength();
            CigarData.erase(CigarData.begin(), CigarData.begin()+2);
         }
         else if (CigarData[1].getType() == 'I') {
            alnchop += CigarData[1].getLength();
            querychop += CigarData[1].getLength();
            CigarData.erase(CigarData.begin(), CigarData.begin()+2);
         }
         else {
            throw logic_error(to_string(__LINE__) + ":DEBUG unexpected cigar pattern " + getCigarString());
         }
      }
      else {  // front cigar len < len
         // 3M1D135M MD: 3^G1G133 len=5 3match, 1D, 1G
         unsigned int c=0;
         size_t x=0;
         while (c < CigarData.size()) {
            if (CigarData[c].getType() == 'M') {
               if (x + CigarData[c].getLength() < len) {
                  x += CigarData[c].getLength();
                  ++c; 
               }
               else if (x + CigarData[c].getLength() == len) {
                  throw logic_error(to_string(__LINE__) + ":DEBUG write more code cigar insertion state " + getCigarString() + " c=" + to_string(c) + " x=" + to_string(x));
               }
               else { // larger than 
                  CigarData[c].shrink(len-x);
                  break; // finished
               }
            }
            else if (CigarData[c].getType() == 'D') {
               refadvance += CigarData[c].getLength();
               alnchop += CigarData[c].getLength();
               ++c;
            }
            else if (CigarData[c].getType() == 'I') {
               throw logic_error("write more code cigar insertion state " + getCigarString()
                     + " c=" + to_string(c) + " x=" + to_string(x));
            }
            else {
               throw logic_error("unexpected cigar op " + getCigarString()  + " c=" + to_string(c) + " x=" + to_string(x));
            }
         }
         CigarData.erase(CigarData.begin(), CigarData.begin()+c);
      }
      SupportData.QuerySequenceLength -= querychop;
      QueryBases = QueryBases.substr(querychop);         
      Qualities = Qualities.substr(querychop);          
   }
   else {  // front is softclip
      if (CigarData[1].getType() != 'M') { 
         cerr << __FILE__ << ":" << __LINE__ << " write more code for Cigar modification\n";
         throw logic_error(string(__FILE__) + to_string(__LINE__) + ":ERROR need more work on first not M operation in Cigar: " + getCigarString());
      }
      if (CigarData[1].getLength() > len) { // enough length to remove
         CigarData[1].shrink(len);
         CigarData.front().expand(len);
      }
      else if (CigarData[1].getLength() == len) {
         if (CigarData.size() <= 3) {
            throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":DEBUG not enough cigar segments at front: " + getCigarString() + " write more code");
         }
         if (CigarData[2].getType() == 'D') { // 21S2M15D91M5I28M
            alnchop += CigarData[2].getLength();
            refadvance += CigarData[2].getLength();
            CigarData[0].expand(len); // query missing 15D, not affecting softclip len
            CigarData.erase(CigarData.begin()+1, CigarData.begin()+3);
         }
         else if (CigarData[2].getType() == 'I') {
            alnchop += CigarData[2].getLength();
            querychop += CigarData[2].getLength();
            CigarData[0].expand(querychop);
            CigarData.erase(CigarData.begin()+1, CigarData.begin()+3);
         }
         else {
            throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":DEBUG unexpected cigar pattern at front: " + getCigarString() + " write more code");
         }
      }
      else { // < len
         unsigned int c=1;
         size_t x=0;
         while (c < CigarData.size()) {
            if (CigarData[c].getType() == 'M') {
               if (x + CigarData[c].getLength() < len) {
                  x += CigarData[c].getLength();
                  ++c; 
               }
               else if (x + CigarData[c].getLength() == len) {
                  throw logic_error(to_string(__LINE__) + ":DEBUG write more code cigar insertion state " + getCigarString() + " c=" + to_string(c) + " x=" + to_string(x));
               }
               else { // larger than 
                  CigarData[c].shrink(len-x);
                  break; // finished
               }
            }
            else if (CigarData[c].getType() == 'D') {
               refadvance += CigarData[c].getLength();
               alnchop += CigarData[c].getLength();
               ++c;
            }
            else if (CigarData[c].getType() == 'I') {
               throw logic_error("write more code cigar insertion state " + getCigarString()
                     + " c=" + to_string(c) + " x=" + to_string(x));
            }
            else {
               throw logic_error("unexpected cigar op " + getCigarString()  + " c=" + to_string(c) + " x=" + to_string(x));
            }
         }
         CigarData.erase(CigarData.begin()+1, CigarData.begin()+c);
         CigarData.front().expand(len);
         //throw logic_error("len " + to_string(len) + " is longer than cigar match at front with softclip " + getCigarString() + " write more code");
      }
   }
   // common operation in different cases
   if (AlignedBases.size() > alnchop) {
      string::size_type x=0;
      size_t i=0;
      while (x < alnchop) {
         if (AlignedBases[i] == '-') {
            ++i;
         }
         else {
            ++i; ++x;
         }
      }
      AlignedBases=AlignedBases.substr(i);
   }
   Position += refadvance;
   if (getInsertSize() > 0) {
      InsertSize -= refadvance;        
   }
   else if (getInsertSize() < 0) {
      InsertSize += refadvance;
   }
}

void BamAlignment::chopBack(size_t len, int numMismatch) {
   // if InsertSize is zero do nothing
   if (numMismatch > 0) {
      int8_t NMval = getNMValue();
      if (static_cast<int>(NMval) >= numMismatch) {
         NMval -= static_cast<int8_t>(numMismatch);
         EditTag("NM", "C", NMval);
      }
      else {
         throw runtime_error("NM value update error");
      }
   }
   unsigned int chopQuery=len;
   unsigned int chopAlign=len;
   unsigned int chopRef = len; 
   // with and without softclip different part
   if (!endWithSoftclip()) { // no softclip at 3' end
      if (CigarData.back().getType() != 'M') {
         throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR expected M as last cigarop " + getCigarString());
      }
      if (CigarData.back().getLength() > len) {
         CigarData.back().shrink(len);
      }
      else if (CigarData.back().getLength() == len) {
         if (CigarData.size() <= 2) {
            throw runtime_error("cigar " + getCigarString() + " too small");
         }
         unsigned int c=CigarData.size()-1;
         if (CigarData.back().getType() == 'M' && CigarData[c-1].getType() == 'D') {
            //123M5D3M
            chopAlign += CigarData[c-1].getLength();
            chopRef = chopAlign;
            CigarData.resize(CigarData.size()-2);
         }
         else if (CigarData.back().getType() == 'M' && CigarData[c-1].getType() == 'I') {
            //72M6D68M1I10M MD: 72^GTGTGA1T52C13T1G3T3 3' end of Query CCTTTTTTNNTNTTTNTTT
            chopAlign += CigarData[c-1].getLength();
            chopQuery += CigarData[c-1].getLength();
            CigarData.resize(CigarData.size()-2);
         }
         else {
            cerr << *this << endl;
            throw runtime_error(string(__FILE__) + ":" + to_string(__LINE__)
                  + " write more code for " + getQueryName() + " " + getCigarString());
         }
      }
      else { //if (CigarData.back().getLength() < len) 
         //unsigned int c=CigarData.size()-1;
         int c=CigarData.size()-1;
         size_t x=0;
         while (c >= 0) {
            if (CigarData[c].getType() == 'M') {
               if (x + CigarData[c].getLength() < len) {
                  x += CigarData[c].getLength();
                  --c; 
               }
               else if (x + CigarData[c].getLength() == len) {
                  throw logic_error(to_string(__LINE__) + ":DEBUG write more code cigar insertion state " + getCigarString() + " c=" + to_string(c) + " x=" + to_string(x));
               }
               else { // larger than 
                  CigarData[c].shrink(len-x);
                  break; // finished
               }
            }
            else if (CigarData[c].getType() == 'D') {
               chopAlign += CigarData[c].getLength();
               chopRef += CigarData[c].getLength();
               --c;
            }
            else if (CigarData[c].getType() == 'I') {
               throw logic_error("write more code cigar insertion state " + getCigarString()
                     + " c=" + to_string(c) + " x=" + to_string(x));
            }
            else {
               throw logic_error(to_string(__LINE__) + ":DEBUG unexpected cigar op " + getCigarString()  + " c=" + to_string(c) + " x=" + to_string(x));
            }
         }
         CigarData.erase(CigarData.begin()+c+1, CigarData.end());
      }
      SupportData.QuerySequenceLength -= chopQuery;
      QueryBases.resize(getLength());
      Qualities.resize(getLength());
   }
   else {  // end with softclip
      //unsigned int c=CigarData.size()-2; // last M index
      int c=CigarData.size()-2; // last M index
      if (CigarData[c].getType() != 'M') { //such as 121M50S
         cerr << __FILE__ << ":" << __LINE__ << " write more code for "
            << getCigarString() << " Cigar modification\n" << *this << endl;
         throw runtime_error("need more work on last before S not M operation in Cigar=" + getCigarString());
      }
      //if (CigarData[c].getType() == 'M') { //such as 121M50S
      if (CigarData[c].getLength() > len) {
         CigarData[c].shrink(len);
         CigarData.back().expand(len);
      }
      else if (CigarData[c].getLength() == len) {
         if (CigarData.size() <= 3) {
            throw logic_error("cigar " + getCigarString() + " not enough segment for len "
                  + to_string(len) + " trimming from back");
         }
         if (CigarData[c-1].getType() == 'D') {
            // such as 112M5D2M50S => 112M(2+50)S Note 5D does not
            // contribute to softclip
            chopAlign += CigarData[c-1].getLength();
            chopRef = chopAlign;
            CigarData.back().expand(len);
            CigarData.erase(CigarData.end()-3, CigarData.end()-1);
         }
         else if (CigarData[c-1].getType() == 'I') {
            chopAlign += CigarData[c-1].getLength();
            chopQuery += CigarData[c-1].getLength();
            CigarData.back().expand(chopQuery);
            CigarData.erase(CigarData.end()-3, CigarData.end()-1);
         }
         else {
            throw runtime_error("write more code for cigar " + getCigarString());
         }
      }
      else { //if (len > CigarData[lastIdx].getLength()) {
         if (CigarData.size() < 3) {
            throw logic_error("cigar must have 3 or more segments " + getCigarString());
         }
         size_t x=0;
         while (c >= 0) {
            if (CigarData[c].getType() == 'M') {
               if (x + CigarData[c].getLength() < len) {
                  x += CigarData[c].getLength();
                  --c; 
               }
               else if (x + CigarData[c].getLength() == len) {
                  throw logic_error(to_string(__LINE__) + ":DEBUG write more code cigar insertion state " + getCigarString() + " c=" + to_string(c) + " x=" + to_string(x));
               }
               else { // larger than 
                  CigarData[c].shrink(len-x);
                  break; // finished
               }
            }
            else if (CigarData[c].getType() == 'D') {
               chopAlign += CigarData[c].getLength();
               chopRef += CigarData[c].getLength();
               --c;
            }
            else if (CigarData[c].getType() == 'I') {
               throw logic_error("write more code cigar insertion state " + getCigarString()
                     + " c=" + to_string(c) + " x=" + to_string(x));
            }
            else {
               throw logic_error(to_string(__LINE__) + ":DEBUG unexpected cigar op " + getCigarString()  + " c=" + to_string(c) + " x=" + to_string(x));
            }
         }
         CigarData.back().expand(len);
         CigarData.erase(CigarData.begin()+c+1, CigarData.end()-1);
      }
   }
   if (AlignedBases.size() > chopAlign) {
      string::size_type x=0;
      size_t i=AlignedBases.size()-1;
      while (x < chopAlign) {
         if (AlignedBases[i] == '-') {
            --i;
         }
         else {
            --i; ++x;
         }
      }
      AlignedBases=AlignedBases.substr(0, i+1);
   }
   if (getInsertSize() > 0) {
      InsertSize -= chopRef;        
   }
   else if (getInsertSize() < 0) {
      InsertSize += chopRef;
   }
}

void BamAlignment::chopBefore(int idx) {
   if (idx <= getPosition()) {
      cerr << __FILE__ << ":" << __LINE__ << ": " << getName() << " idx=" 
         << idx << " not greater than position " << getPosition() << endl;
      throw logic_error("idx not greater than position in ChopBefore()");
   }
   /*
   bool inbug=false;
   if (getQueryName() == "S32897") {
      //cerr << *this << "before chopBefore idx=" << idx << endl;
      cerr << __LINE__ << "before chopBefore idx=" << idx << endl;
      inbug=true;
   }
   */
   int ri = getPosition();
   int qi = 0;
   vector<CigarOp> newcigar;
   int inscnt = 0;
   auto it = CigarData.begin();
   while (it != CigarData.end()) {
      if (it->getType() == 'S' || it->getType() == 'H') { // will ignore 
         qi += it->getLength();
      }
      else if (it->getType() == 'M') {
         if (ri+static_cast<int>(it->getLength()) > idx) { // idx inside this cigar segment
            // ri  idx   ri+Cgrlen idx will be the new position
            // |    |    |
            // -----=====
            newcigar.push_back(CigarOp('M', ri + it->getLength() - idx)); // idx will be retained in new object
            qi += (idx - ri); // advance to postion corresponds to idx
            ++it;
            break;
         }
         else {
            // ri ri+Cgrlen idx
            // |         |  | 
            // ==========  
            ri += it->getLength();
            qi += it->getLength();
         }
      }
      else if (it->getType() == 'I') { // insertion is between postion 
         inscnt += it->getLength();
         qi += it->getLength();
      }
      else if (it->getType() == 'D') { 
         if (ri + static_cast<int>(it->getLength()) <= idx) {
            ri += it->getLength();
         }
         else { // idx fall inside deletion of query
            //cerr << endl << *this;
            //cerr << __LINE__ << ": c=" << it - CigarData.begin() << " ri=" << ri << " qi=" << qi
            //   << " idx=" << idx << " is inside deletion will be ignored" << endl;
            // alignment cannot end with D state
            // will not add D to newcigar, new idx will be moved to first base of next M
            idx = ri + it->getLength();
            ++it;
            break;
            //throw logic_error("BamAlignment::chopBefore() idx is inside deletion");
            //break;
         }
      } 
      else { // cigar type N or other not considered
         throw logic_error("Cigar type: " + string(1, it->getType()) + " not considered in BamAlignment::chopBefore()");
      }
      ++it;
   }
   while (it != CigarData.end()) {
      newcigar.push_back(*it);
      ++it;
   }
   try { // must do this first before editing other fields
      //if (inbug)
      //   cerr << "before chopMDBefore() idx=" << idx << endl;
      int mismatchCnt = chopMDBefore(idx) + inscnt; // update MD tag
      reduceNMTag(mismatchCnt);
   }
   catch (const logic_error& ler) {
      cerr << __LINE__ << ": idx=" << idx << " pos=" << getPosition() << endl;
      cerr << ler.what() << endl;
      throw;
   }
   // insert size, position will not be affected
   SupportData.QuerySequenceLength -= qi;
   QueryBases = QueryBases.substr(qi);
   Qualities = Qualities.substr(qi);
   CigarData = std::move(newcigar); // need to update MC tag of mate
   // only need to count D and mismatch from MD tag after idx
   // must chop MD first then update position!
   setPosition(idx); // update position if chopping from head
   clearAlignedBases();
   if (!validCigar()) {
      cerr << *this;
      throw logic_error("invalid cigar after chopBefore operation()");
   }
   if (!refwidthAgreeWithMD()) {
      cerr << *this;
      throw logic_error("bad MD after chopBefore()");
   }
}

void BamAlignment::chopAfter(int idx) {
   //assert(idx < getEndPosition());
   if (idx >= getEndPosition()) {
      throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ": " + to_string(idx) + " after end " + to_string(getEndPosition()) + " invalid operation");
   }
   size_t c = 0;
   int ri = getPosition();
   int qi = 0;
   vector<CigarOp> newcigar;
   //if (getName() == "S16634") {
   //   cerr << __LINE__ << ":  at start of chopAfter() ri=" << ri
   //      << " idx=" << idx << " current end " << getEndPosition() << endl;
   //}
   while (c < CigarData.size()) {
      if (CigarData[c].getType() == 'S' || CigarData[c].getType() == 'H') {
         newcigar.push_back(CigarData[c]);
         qi += CigarData[c].getLength();
      }
      else if (CigarData[c].getType() == 'M') {
         if (ri+static_cast<int>(CigarData[c].getLength()) <= idx) {
            // ri ri+cigarlen idx
            // |            |  |
            // =============   
            ri += CigarData[c].getLength();
            qi += CigarData[c].getLength();
            newcigar.push_back(CigarData[c]);
         }
         else { // idx inside this Cigar segment
            // ri   idx      ri+CigarLen > idx
            // |     |        |
            // =======--------
            newcigar.push_back(CigarOp('M', idx-ri+1)); // idx will be retained in new object
            //ri = idx; // move ri to this location, ri not used no need to do it
            qi += (idx - ri);
            ++c;
            break;
         }
      }
      else if (CigarData[c].getType() == 'I') {
         qi += CigarData[c].getLength();
         newcigar.push_back(CigarData[c]);
      }
      else if (CigarData[c].getType() == 'D') { 
         if (ri + static_cast<int>(getCigarLength(c)) <= idx) {
            ri += getCigarLength(c);
            newcigar.push_back(CigarData[c]);
         }
         else { // idx fall inside deletion of query, this is legal
            //   idx
            //   ri      idx  ri+cigarlen 
            //   |       |    |   
            //===-------------=======
            //      QDEL      qi
            // alignment cannot end with D state
            // will not add D to newcigar, D will be discarded
            ++c;
            //cerr << endl << *this;
            //cerr << __LINE__ << ": c=" << c << " ri=" << ri << " qi=" << qi
            //   << " idx=" << idx << " inside deletion will discard" << endl;
            // new idx will be ri-1, ri is now the first base in the DEL
            idx = ri-1;
            --qi;
            //throw logic_error("BamAlignment::chopAfter() idx fall inside D of query sequence");
            break;
         }
      } 
      else { // cigar type N or other not considered
         throw logic_error("Cigar type: " + string(1, CigarData[c].getType()) + " not considered in BamAlignment::chopAfter()");
      }
      ++c;
   }
   int inscnt=0; // number of query inserted total bases
   while (c < CigarData.size()) {
      if (CigarData[c].getType() == 'I')
         inscnt += CigarData[c].getLength();
      ++c;
   }
   try {
      // only need to count D and mismatch from MD tag after idx
      int mismatchCnt = chopMDAfter(idx) + inscnt; // update MD tag
      reduceNMTag(mismatchCnt);
   }
   catch (const logic_error& ler) {
      //cerr << endl << *this << __LINE__ << ": idx=" << idx << endl << ler.what() << endl;
      cerr << __LINE__ << ": idx=" << idx << endl << ler.what() << endl;
      throw;
   }
   // insert size, position will not be affected
   SupportData.QuerySequenceLength = qi+1;
   QueryBases.resize(qi+1);
   Qualities.resize(qi+1);
   // check last setment is D or I
   if (newcigar.back().getType() == 'D') {
      cerr << "warning last segment is D is not permissted will remove\n";
      newcigar.resize(newcigar.size()-1);
   }
   if (newcigar.back().getType() == 'I') {
      cerr << "Warning last segment is I will be converted to S\n";
      newcigar.back().setType('S');
   }
   CigarData = std::move(newcigar);
   clearAlignedBases();
   //if (getName() == "S858") {
   //   cerr << *this << "  at end of chopAfter()\n";
   //}
   if (!refwidthAgreeWithMD()) {
      throw logic_error("bad MD at end of chopAfter()");
   }
}

int BamAlignment::chopMDBefore(int idx) {
   if (idx <= getPosition()) {
      cerr << __FILE__ << ":" << __LINE__ << ": idx=" << idx << " not greater than position "
         << getPosition() << endl;
      throw logic_error(string(__func__) + ": idx and position are the same cannot chop");
   }
   // mdvec second has one fewer element
   pair<vector<int>, vector<string>> mdvec = getMDArray();
   Matchdiff mdchopper(std::move(mdvec.first), std::move(mdvec.second));
   idx -= getPosition(); // convert to 0-index from chromosome index
   //cerr << __LINE__ << ": 0-based idx=" << idx << endl;
   int diff_inhead = mdchopper.removeBefore(idx);
   string newmdtag = mdchopper.toString();
   EditTag("MD", "Z", newmdtag);
   return diff_inhead;
}

// count the number of mismatch after idx
// return the number of mismatched after idx
int BamAlignment::chopMDAfter(int idx) {
   // mdvec second has one fewer element
   pair<vector<int>, vector<string>> mdvec = getMDArray();
   Matchdiff mdchopper(std::move(mdvec.first), std::move(mdvec.second));
   idx -= getPosition(); // convert to 0-index from chromosome index
   int diff_intail = mdchopper.removeAfter(idx);
   string newmdtag = mdchopper.toString();
   EditTag("MD", "Z", newmdtag);
   return diff_intail;
}

void BamAlignment::reduceNMTag(int diff) {
   pair<uint8_t, bool> nmv;
   try {
      nmv = getTag<uint8_t>("NM");
      if (!nmv.second) {
         cerr << endl << *this << endl;
         throw logic_error(to_string(__LINE__) + ": No NM tag");
      }
   }
   catch (const BamTypeException& ter) {
      //cerr << ter.what() << endl;
      auto [nmval, hasNM] = getTag<int32_t>("NM"); // number of mismatch
      if (!hasNM) {
         cerr << *this << __FILE__ << ":" << __LINE__ << "Failed to get NM tag\n";
         throw logic_error(string(__func__) + ":ERROR Bam dose not have NM tab");
      }
      if (nmval < static_cast<int32_t>(UINT8_MAX)) {
         nmv.first = static_cast<uint8_t>(nmval);
      }
      else {
         throw runtime_error("NM value " + to_string(nmval) + " is greater than UINT8_MAX");
      }
   }
   //if (!nmv.second) {
   //   cerr << *this << __FILE__ << ":" << __LINE__ << "Failed to get NM tag\n";
   //   throw logic_error(string(__func__) + ":ERROR Bam dose not have NM tab");
   //}
   if (diff > static_cast<int>(UINT8_MAX)) {
      throw runtime_error("mismatch value too large to be store as uint8_t");
   }
   if (diff > static_cast<int>(nmv.first)) {
      throw logic_error(string(__func__) + ":ERROR diff mismatch more than entire NM tag value");
   }
   //nmval -= diff;
   nmv.first -= static_cast<uint8_t>(diff);
   if (!EditTag("NM", "C", nmv.first)) {
      throw logic_error(string(__func__) + ":ERROR Failed to edit NM tag");
   }
}

// len is ref count only
int BamAlignment::countFrontMismatch(int len) const {
   // mdvec second has one fewer element
   pair<vector<int>, vector<string>> mdvec = getMDArray();
   size_t d = 0; // index into mdvec.first
   int nummismatch = 0;
   while (len > 0 && d < mdvec.first.size()) {
      if (d >= mdvec.second.size()) {
         throw logic_error("d out of bound in BamAlignment::countFrontMismatch()");
      }
      if (mdvec.second[d].front() == '^') {
         if (len == static_cast<int>(mdvec.second[d].size())) {
            throw logic_error("after " + to_string(len) + " from left is a ref deletion: "
                  + mdvec.second[d]);
         }
         if (static_cast<int>(mdvec.second[d].size()) <= len) {
            len -= mdvec.second[d].size();
            nummismatch += (mdvec.second[d].size()-1);
         }
         else { // left match longer than len done, no mismatch
            break;
         }
      }
      else {
         len -= mdvec.second[d].size();
         if (len <= 0) break; // no mismatch
         if (static_cast<int>(mdvec.second[d].size()) <= len) {
            nummismatch += mdvec.second[d].size();
            len -= mdvec.second[d].size();
         }
         else { // the mismatch is longer then the remaining len
            len = 0;
         }
      }
      ++d;
   }
   return nummismatch;
}

int BamAlignment::countBackMismatch(int len) const {
   // mdvec second has one fewer element
   pair<vector<int>, vector<string>> mdvec = getMDArray();
   //int i=0;
   int d = mdvec.first.size() - 1; // index into mdvec.first
   int nummismatch = 0;
   while (len > 0 && d > -1) {
      if (d-1 < 0) {
         throw logic_error("d out of bound in BamAlignment::countBackMismatch()");
      }
      if (mdvec.second[d-1].front() == '^') {
         if (len == static_cast<int>(mdvec.second[d].size())) {
            throw logic_error("after " + to_string(len) + " is a ref deletion: "
                  + mdvec.second[d-1]);
         }
         if (static_cast<int>(mdvec.second[d].size()) <= len) {
            len -= mdvec.second[d].size();
            nummismatch += (mdvec.second[d].size()-1);
         }
         else { // left match longer than len done, no mismatch
            break;
         }
      }
      else {
         len -= mdvec.second[d].size();
         if (len < 0) break; // no mismatch
         if (static_cast<int>(mdvec.second[d].size()) <= len) {
            nummismatch += mdvec.second[d].size();
            len -= mdvec.second[d].size();
         }
         else len = 0;
      }
      ++d;
   }
   return nummismatch;
}

// remove the first len bases
// assume has MD tag for performance if no MD then
// need to write more code
// Caller need to update mateposition from the other mate
bool BamAlignment::trimFront() {
   pair<vector<int>, vector<string>> mdvec = getMDArray();
   int trimlen=0, mismatch=0;
   size_t i = 0;
   while (i < mdvec.first.size() && mdvec.first[i] < 4) {
      if (i >= mdvec.second.size()) {
         throw runtime_error("out of range in trimFront()");
      }
      if (mdvec.second[i].front() != '^') {
         trimlen += (mdvec.first[i] + 1);
         ++mismatch;
      }
      else {
         trimlen += mdvec.first[i];
      }
      ++i;
   }
   // trim raw data if trimlen > 0
   if (trimlen > 0) {
#ifdef DEBUG
      cerr << "before trimFront()\n" << *this << endl;
#endif
      chopFront(trimlen, mismatch);
      updateMDTag(mdvec);
#ifdef DEBUG
      cerr << "after trimFront()\n" << *this << endl;
#endif
      return true;
   }
   return false;
}

// MD tag such as 130A1C0
// 130, 1, 0 | A, C
bool BamAlignment::trimBack() {
   pair<vector<int>, vector<string>> mdvec = getMDArray();
   int trimlen=0;
   int mismatch=0;
   int i = mdvec.first.size()-1;
   while (i > -1 && mdvec.first[i] < 4) {
      if (i-1 < 0) {
         throw runtime_error("i index out of range in trimBack()");
      }
      if (mdvec.second[i-1].front() == '^') {
         trimlen += mdvec.first[i];
      }
      else {
         trimlen += (mdvec.first[i] + 1);
         ++mismatch;
      }
      --i;
   }
   // trim raw data if trimlen > 0
   if (trimlen > 0) {
#ifdef DEBUG
      cerr << "before trimBack():\n" << *this << endl;
#endif
      chopBack(trimlen, mismatch);
      updateMDTag(mdvec);
#ifdef DEBUG
      cerr << "after trimBack():\n" << *this << endl;
#endif
      return true;
   }
   return false;
}

pair<bool,bool> BamAlignment::trim() {
   pair<bool,bool> trimmed;
   trimmed.first = trimFront();
   trimmed.second = trimBack();
   return trimmed;
}

// in case A tailing extra base could be added now
void BamAlignment::patchEnd() {
   pair<vector<int>, vector<string>> mdvec = getMDArray();

   static const int TRIMLEN_MAX = 6;
   static const int GAP_CUT = 3;
   //map<int, char> correction; // only one, pos, char of ref
   size_t trimlen=0, m=0, qi=0, c=0;
   bool patchFront = true;
   //cerr << CigarData.size() << "Cigar size\n";
   if (CigarData.front().getType() == 'M') {
   }
   else if (CigarData.front().getType() == 'S') {
      qi = CigarData.front().getLength();
      c=1;
   }
   else {
      // need to do more reserch not possible to do patch
      patchFront=false;
   }
   if (patchFront) {
      while (trimlen < TRIMLEN_MAX && mdvec.first[m] < GAP_CUT && mdvec.second[m].front() != '^'
            && CigarData[c].getLength() > trimlen + mdvec.first[m]) 
      {
         qi += mdvec.first[m];
         QueryBases[qi] = mdvec.second[m].front();
         trimlen += (mdvec.first[m] + 1);
         ++m; ++qi;
      }
      if (m > 0) {
         mdvec.first.erase(mdvec.first.begin(), mdvec.first.begin()+m);
         mdvec.first.front() += trimlen;
         mdvec.second.erase(mdvec.second.begin(), mdvec.second.begin()+m);
         updateMDTag(mdvec);
         uint8_t nmval = getNMValue();
         assert(m < nmval);
         nmval -= m;
         EditTag("NM", "C", nmval);
      }
   }

   if (mdvec.first.size() <= 1) { return; }
   trimlen=0;
   m=mdvec.first.size()-1;
   qi=getLength()-1;
   c=CigarData.size()-1;
   if (CigarData[c].getType() == 'S') {
      qi -= CigarData[c].getLength();
      --c;
   }
   else if (CigarData[c].getType() == 'M') { }
   else return;

   while (trimlen < TRIMLEN_MAX && m > 0 && mdvec.first[m] < GAP_CUT && mdvec.second[m-1].front() != '^'
            && CigarData[c].getLength() > trimlen + mdvec.first[m]) 
   {
      qi -= mdvec.first[m];
      QueryBases[qi] = mdvec.second[m-1].front();
      trimlen += (mdvec.first[m]+1);
      --m; --qi;
//#ifdef DEBUG
//      cerr << "m=" << m << " qi=" << qi << " trimlen=" << trimlen << endl;
//#endif
   }
   if (m < mdvec.first.size()-1) {
      uint8_t nmval = getNMValue();
      nmval -= (mdvec.first.size()-1 - m);
      EditTag("NM", "C", nmval);
      mdvec.first.erase(mdvec.first.begin()+m+1, mdvec.first.end());
      mdvec.first.back() += trimlen;
      mdvec.second.erase(mdvec.second.begin()+m, mdvec.second.end());
      updateMDTag(mdvec);
   }
}

// regenerate AlignedBases
void BamAlignment::updateNMTag(const string& refseq) {
   int b = getPosition();
   int e = GetEndPosition(); // one passed the end [b,e)
   string subseq = refseq.substr(b, e-b);
   uint8_t edit=0;
   size_t ci=0, ri=0, qi=0, cigarIdx;
   AlignedBases.clear();
   //cout << "reference sequence:\n"
   //   << subseq << "\nquery bases\n"
   //   << getQueryBases() << endl;

   while (ci < CigarData.size()) {
      //cout << "ri: " << ri << " qi: " << qi << endl;
      if (CigarData[ci].Type == 'S' || CigarData[ci].Type == 'H') {
         qi += CigarData[ci].Length;
      }
      else if (CigarData[ci].Type == 'M') {
         AlignedBases += QueryBases.substr(qi, CigarData[ci].Length);
         cigarIdx=0;
         while (cigarIdx < CigarData[ci].Length) {
            if (toupper(subseq[ri]) != QueryBases[qi]) {
               //cout << "diff ref: " << subseq[ri] << " | " << QueryBases[qi] << endl;
               ++edit;
            }
            // assume bam file use upper cases for all bases!
            ++cigarIdx; ++ri; ++qi;
         }
      }
      else if (CigarData[ci].Type == 'I') {
         AlignedBases += QueryBases.substr(qi, CigarData[ci].Length);
         //cout << "query insert: " << QueryBases.substr(qi, CigarData[ci].Length) << endl;
         edit += CigarData[ci].Length;
         qi += CigarData[ci].Length;
      }
      else if (CigarData[ci].Type == 'D') {
         AlignedBases += string(CigarData[ci].Length, '-');
         //cout << "query deletion: " << CigarData[ci].Length << endl;
         edit += CigarData[ci].Length;
         ri += CigarData[ci].Length;
      }
      else {
         throw runtime_error("CigarOP " + string(1, CigarData[ci].Type) 
               + " not considered inside " + string(__func__));
      }
      ++ci;
   }
   //cout << "recalculated edit distance: " << edit << endl;
   if (hasTag("NM")) {
      EditTag("NM", "C", edit);
   }
   else {
      AddTag("NM", "C", edit);
   }
}

bool BamAlignment::hasSoftclip() const {
   if (CigarData.empty()) return false;
   return CigarData.front().Type == 'S' || CigarData.back().Type == 'S';
}

int BamAlignment::getASValue() const {
   if (!hasTag("AS")) return -1;
   int val=-1;
   try {
      pair<int,bool> res = getTag<int>("AS");
      if (res.second) val = res.first;
   }
   catch (const BamTypeException& ler) {
      pair<uint32_t,bool> res = getTag<uint32_t>("AS");
      if (res.second) val = res.first;
   }
   return val;
}

uint8_t BamAlignment::getNMValue() const {
   //if (!hasTag("NM")) return -1;
   uint8_t val = 0;
   //int val = -1;
   try {
      pair<uint8_t,bool> res = getTag<uint8_t>("NM");
      if (res.second) val = res.first;
   }
   catch (const BamTypeException& ler) {
      pair<int,bool> res = getTag<int32_t>("NM");
      //if (res.second) val = res.first;
      if (res.first > static_cast<int>(UINT8_MAX)) {
         throw logic_error("NM value " + to_string(res.first) + " cannot be stored as uint8_t");
      }
      val = static_cast<uint8_t>(res.first);
   }
   catch (const exception& err) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR failed to get NM tag value\n";
      throw;
   }
   //if (val < -1) {
   //   cerr << *this << endl;
   //   cerr << __FILE__ << ":" << __LINE__ << ": val=" << val << endl;
   //   throw logic_error("check getTag");
   //}
   return val;
}

int BamAlignment::getTemplateLength() const {
   int tlen = getInsertSize();
   if (tlen == 0 && HasTag("XO")) {
      tlen = getReferenceWidth();
   }
   return abs(tlen);
}

string BamAlignment::getMatchedQuerySequence() const {
   int b = 0;
   if (CigarData.front().getType() == 'S' || CigarData.front().getType() == 'H') {
      b = CigarData.front().getLength();
   }
   int e = getLength();
   if (CigarData.back().getType() == 'S' || CigarData.back().getType() == 'H') {
      e -= CigarData.back().getLength();
   }
   return QueryBases.substr(b, e-b);
}

pair<int,int> BamAlignment::getMatchBound() const {
   int b = 0;
   if (CigarData.front().getType() == 'S' || CigarData.front().getType() == 'H') {
      b = CigarData.front().getLength();
   }
   int e = getLength();
   if (CigarData.back().getType() == 'S' || CigarData.back().getType() == 'H') {
      e -= CigarData.back().getLength();
   }
   return make_pair(b, e);
}

int BamAlignment::getMatchedQueryLength() const {
   if (CigarData.empty() || !hasTag("NM")) return 0;
   int len=getLength();
   if (CigarData.front().getType() == 'S' || CigarData.front().getType() == 'H') {
      len -= CigarData.front().getLength();
   }
   if (CigarData.back().getType() == 'S' || CigarData.back().getType() == 'H') {
      len -= CigarData.back().getLength();
   }
   return len;
}

int BamAlignment::numberBaseAligned() const {
   int len=getLength();
   for (auto& cd : CigarData) {
      if (cd.getType() == 'H' || cd.getType() == 'S' || cd.getType() == 'D' || cd.getType() == 'I') 
         len -= cd.getLength();
   }
   return len;
}

void BamAlignment::makeUnmapped() {
   setUnmapped(); setMateUnmapped();
   setImproperPair();
   setReferenceId(-1); setMateReferenceId(-1);
   setPosition(-1); setMatePosition(-1);
   CigarData.clear();
   clearAlignedBases();
   // remove tags
   removeTag("NM"); removeTag("MD"); removeTag("MC");
}



