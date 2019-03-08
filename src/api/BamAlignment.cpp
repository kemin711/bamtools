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

using namespace BamTools;
using namespace std;

/*! \class BamTools::BamAlignment
    \brief The main BAM alignment data structure.

    Provides methods to query/modify BAM alignment data fields.
*/
/*! \var BamAlignment::Name
    \brief read name
*/
/*! \var BamAlignment::Length
    \brief length of query sequence
*/
/*! \var BamAlignment::QueryBases
    \brief 'original' sequence (as reported from sequencing machine)

    \note Setting this field to "*" indicates that the sequence is not to be stored on output.
    In this case, the contents of the Qualities field should be invalidated as well (cleared or marked as "*").
*/
/*! \var BamAlignment::Qualities
    \brief FASTQ qualities (ASCII characters, not numeric values)

    \note Setting this field to "*" indicates to BamWriter that the quality scores are not to be stored,
    but instead will be output as a sequence of '0xFF'. Otherwise, QueryBases must not be a "*" and
    the length of this field should equal the length of QueryBases.
*/
/*! \var BamAlignment::TagData
    \brief tag data (use the provided methods to query/modify)
*/
/*! \var BamAlignment::RefID
    \brief ID number for reference sequence
*/
/*! \var BamAlignment::Position
    \brief position (0-based) where alignment starts
*/
/*! \var BamAlignment::Bin
    \brief BAM (standard) index bin number for this alignment
*/
/*! \var BamAlignment::MapQuality
    \brief mapping quality score
*/
/*! \var BamAlignment::CigarData
    \brief CIGAR operations for this alignment
*/
/*! \var BamAlignment::MateRefID
    \brief ID number for reference sequence where alignment's mate was aligned
*/
/*! \var BamAlignment::MatePosition
    \brief position (0-based) where alignment's mate starts
*/
/*! \var BamAlignment::Filename
    \brief name of BAM file which this alignment comes from
*/

/*! \fn BamAlignment::BamAlignment(void)
    \brief constructor
*/
BamAlignment::BamAlignment(void)
    : Length(0)
    , RefID(-1)
    , Position(-1)
    , Bin(0)
    , MapQuality(0)
    , AlignmentFlag(0)
    , MateRefID(-1)
    , MatePosition(-1)
    , InsertSize(0)
{ }

/*! \fn BamAlignment::BamAlignment(const BamAlignment& other)
    \brief copy constructor
*/
BamAlignment::BamAlignment(const BamAlignment& other)
    : Name(other.Name)
    , Length(other.Length)
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
    , Filename(other.Filename)
    , SupportData(other.SupportData)
{ }

BamAlignment::BamAlignment(BamAlignment&& other)
    : Name(std::move(other.Name)), 
      Length(other.Length), 
      QueryBases(std::move(other.QueryBases)),
      AlignedBases(std::move(other.AlignedBases)),
      Qualities(std::move(other.Qualities)),
      TagData(std::move(other.TagData)), 
      RefID(other.RefID), Position(other.Position), Bin(other.Bin), 
      MapQuality(other.MapQuality), AlignmentFlag(other.AlignmentFlag), 
      CigarData(std::move(other.CigarData)), 
      MateRefID(other.MateRefID), MatePosition(other.MatePosition), 
      InsertSize(other.InsertSize),  Filename(other.Filename), 
      SupportData(std::move(other.SupportData))
{ }

BamAlignment& BamAlignment::operator=(const BamAlignment& other) {
   if (this != &other) {
      Name=other.Name;
      Length=other.Length;
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
      Filename=other.Filename;
      SupportData=other.SupportData;
      // ErrorString is not copied
   }
   return *this;
}

BamAlignment& BamAlignment::operator=(BamAlignment&& other) {
   if (this != &other) {
      Name=std::move(other.Name);
      Length=other.Length;
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
      Filename=std::move(other.Filename);
      SupportData=std::move(other.SupportData);
      // ErrorString is not copied
   }
   return *this;
}

/*! \fn BamAlignment::~BamAlignment(void)
    \brief destructor
*/
BamAlignment::~BamAlignment(void) { }

//std::ostream& BamTools::operator<<(std::ostream &ous, const BamTools::BamAlignment &ba) {
// error: ‘std::ostream& BamTools::operator<<(std::ostream&, const BamTools::BamAlignment&)’ should have been declared inside ‘BamTools’
//  std::ostream& BamTools::operator<<(std::ostream &ous, const BamTools::BamAlignment &ba) {
//
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
   vector<std::string> tagNames = ba.GetTagNames();
   // automatic probe is causing some problems
   // removing it.
   for (auto& t : tagNames) {
      char tagtype;
      ba.GetTagType(t, tagtype);
      if (tagtype == 'i') {
         int32_t intval; // must use this type, int is wrong
         ba.GetTag(t, intval);
         ous << t << ": " << intval << "; ";
      }
      else {
         string tagval;
         ba.GetTag(t, tagval);
         ous << t << ": " << tagval << "; ";
      }
   }
   */
   // the following is fine
   string val;
   if (ba.HasTag("BC")) {
      ba.GetTag("BC", val);
      ous << "BC: " << val << sep;
   }
   int ival=0;
   //int32_t ival=0;
   if (ba.HasTag("NM")) {
      if (ba.GetTag("NM", ival)) {
         ous << "NM: " << ival << sep;
      }
      else {
         ous << "NM: FAIL" << sep;
      }
   }
   else {
      cerr << "there is not NM tag!\n";
   }
   if (ba.HasTag("AS")) {
      if (ba.GetTag("AS", ival)) {
         ous << "AS: " << ival << sep;
      }
      else {
         ous << "AS: FAIL" << sep;
      }
   }
   if (ba.HasTag("XS")) {
      if (ba.GetTag("XS", ival)) {
         ous << "XS: " << ival << sep;
      }
      else {
         ous << "XS: FAIL" << sep;
      }
   }

   return ous;
}
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

int BamAlignment::getSoftclipLength() const {
   int res = 0;
   if (getCigar().front().getType() == 'S')
      res += getCigar().front().getLength();
   if (getCigar().back().getType() == 'S')
      res += getCigar().back().getLength();
   return res;
}

void BamAlignment::setCigarOperation(const std::vector<pair<char,int> > &cd) {
   CigarData.resize(cd.size());
   for (size_t i=0; i < CigarData.size(); ++i) {
      CigarData[i].fromPair(cd[i]);
   }
}

// only do one correction!
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
         int edit;
         GetTag("NM", edit);
         edit -= gaplen;
         EditTag("NM", "i", edit);
         return;
      }
   }
}

pair<int,int> BamAlignment::getMismatchCount() const {
   int num_mismatch = 0;
   if (HasTag("NM")) {
      GetTag("NM", num_mismatch);
   }
   else {
      throw runtime_error("No NM tag in bam file");
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
   return make_pair(num_mismatch-indel, alnlen);
}

float BamAlignment::getNGIdentity() const {
   int num_mismatch = 0;
   if (HasTag("NM")) {
      GetTag("NM", num_mismatch);
   }
   else {
      throw runtime_error("No NM tag in bam file");
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
   return (1-(num_mismatch-indel)/(float)alnlen);
}

float BamAlignment::getIdentity() const {
   int num_mismatch = 0;
   if (HasTag("NM")) {
      GetTag("NM", num_mismatch);
   }
   else {
      throw runtime_error("No NM tag in bam file");
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
   return (1-(num_mismatch)/(float)alnlen);
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
                    const string message = string("invalid CIGAR operation type: ") + op.Type;
                    SetErrorString("BamAlignment::BuildCharData", message);
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
                                    const string message = string("invalid binary array type: ") + arrayType;
                                    SetErrorString("BamAlignment::BuildCharData", message);
                                    return false;
                            }
                        }

                        break;
                    }

                    // invalid tag type-code
                    default :
                        const string message = string("invalid tag type: ") + type;
                        SetErrorString("BamAlignment::BuildCharData", message);
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
    while ( numBytesParsed < tagDataLength ) {
        const char* pTagType        = pTagData;
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;
        // check the current tag, return true on match
        if ( strncmp(pTagType, tag.c_str(), 2) == 0 )
            return true;
        // get the storage class and find the next tag
        if ( *pTagStorageType == '\0' ) return false;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return false;
        if ( *pTagData == '\0' ) return false;
    }
    return false;
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
    vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const CigarOp& op = (*cigarIter);

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
    if ( closedInterval )
        alignEnd -= 1;

    // return result
    return alignEnd;
}

/*! \fn std::string BamAlignment::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
std::string BamAlignment::GetErrorString(void) const {
    return ErrorString;
}

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

/*! \fn std::vector<std::string> BamAlignment::GetTagNames(void) const
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
        if ( *pTagType == '\0' ) break;
        if ( !SkipToNextTag(*pTagType, pTagData, numBytesParsed) ) break;
        if ( *pTagData == '\0' ) break;
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

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag not found, return failure
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ){
        // TODO: set error string?
        return false;
    }

    // otherwise, retrieve & validate tag type code
    type = *(pTagData - 1);
    switch ( type ) {
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
            const string message = string("invalid tag type: ") + type;
            SetErrorString("BamAlignment::GetTagType", message);
            return false;
    }
}

/*! \fn bool BamAlignment::HasTag(const std::string& tag) const
    \brief Returns true if alignment has a record for requested tag.

    \param[in] tag 2-character tag name
    \return \c true if alignment has a record for tag
*/
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

bool BamAlignment::IsDuplicate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_DUPLICATE) != 0 );
}

/*! \fn bool BamAlignment::IsFailedQC(void) const
    \return \c true if this read failed quality control
*/
bool BamAlignment::IsFailedQC(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_QC_FAILED) != 0 );
}

bool BamAlignment::IsFirstMate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_READ_1) != 0 );
}

/*! \fn bool BamAlignment::IsMapped(void) const
    \return \c true if alignment is mapped
*/
bool BamAlignment::IsMapped(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_UNMAPPED) == 0 );
}

/*! \fn bool BamAlignment::IsMateMapped(void) const
    \return \c true if alignment's mate is mapped
*/
bool BamAlignment::IsMateMapped(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_MATE_UNMAPPED) == 0 );
}

/*! \fn bool BamAlignment::IsMateReverseStrand(void) const
    \return \c true if alignment's mate mapped to reverse strand
*/
bool BamAlignment::IsMateReverseStrand(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_MATE_REVERSE_STRAND) != 0 );
}

double BamAlignment::getFractionStrand() const {
   if (HasTag("XO")) {
      int32_t overlap;
      GetTag("XO", overlap);
      if (overlap == getReferenceWidth()) {
         return 0;
      }
      else {
         return (double)overlap/getReferenceWidth();
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
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_PAIRED) != 0 );
}

/*! \fn bool BamAlignment::IsPrimaryAlignment(void) const
    \return \c true if reported position is primary alignment
*/
bool BamAlignment::IsPrimaryAlignment(void) const  {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_SECONDARY) == 0 );
}

/*! \fn bool BamAlignment::IsProperPair(void) const
    \return \c true if alignment is part of read that satisfied paired-end resolution
*/
bool BamAlignment::IsProperPair(void) const {
    //return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_PROPER_PAIR) != 0 );
    return ( (AlignmentFlag & PROPER_PAIR) != 0 );
}

bool BamAlignment::IsReverseStrand(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_REVERSE_STRAND) != 0 );
}

/*! \fn bool BamAlignment::IsSecondMate(void) const
    \return \c true if alignment is second mate on read
*/
bool BamAlignment::IsSecondMate(void) const {
    return ( (AlignmentFlag & Constants::BAM_ALIGNMENT_READ_2) != 0 );
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

/*! \fn void BamAlignment::RemoveTag(const std::string& tag)
    \brief Removes field from BAM tags.

    \param[in] tag 2-character name of field to remove
*/
void BamAlignment::RemoveTag(const std::string& tag) {
  
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // skip if no tags available
    if ( TagData.empty() )
        return;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;

    // skip if tag not found
    if  (!FindTag(tag, pTagData, originalTagDataLength, numBytesParsed))
        return;

    // otherwise, remove it
    RaiiBuffer newTagData(originalTagDataLength);

    // copy original tag data up til desired tag
    pTagData       -= 3;
    numBytesParsed -= 3;
    const unsigned int beginningTagDataLength = numBytesParsed;
    newTagDataLength += beginningTagDataLength;
    memcpy(newTagData.Buffer, pOriginalTagData, numBytesParsed);

    // attemp to skip to next tag
    const char* pTagStorageType = pTagData + 2;
    pTagData       += 3;
    numBytesParsed += 3;
    if ( SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) {

        // squeeze remaining tag data
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagDataLength = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData.Buffer + beginningTagDataLength, pTagData, endTagDataLength );

        // save modified tag data in alignment
        TagData.assign(newTagData.Buffer, beginningTagDataLength + endTagDataLength);
    }
}

/*! \fn void BamAlignment::SetErrorString(const std::string& where, const std::string& what) const
    \internal

    Sets a formatted error string for this alignment.

    \param[in] where class/method where error occurred
    \param[in] what  description of error
*/
void BamAlignment::SetErrorString(const std::string& where, const std::string& what) const {
    static const string SEPARATOR = ": ";
    ErrorString = where + SEPARATOR + what;
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

/*! \fn bool BamAlignment::SkipToNextTag(const char storageType, char*& pTagData, unsigned int& numBytesParsed) const
    \internal

    Moves to next available tag in tag data string

    \param[in]     storageType    BAM tag type-code that determines how far to move cursor
    \param[in,out] pTagData       pointer to current position (cursor) in tag string
    \param[in,out] numBytesParsed report of how many bytes were parsed (cumulatively)

    \return \c if storageType was a recognized BAM tag type

    \post \a pTagData       will point to the byte where the next tag data begins.
          \a numBytesParsed will correspond to the cursor's position in the full TagData string.
*/
bool BamAlignment::SkipToNextTag(const char storageType,
                                 char*& pTagData,
                                 unsigned int& numBytesParsed) const
{
    switch (storageType) {

        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            ++numBytesParsed;
            ++pTagData;
            break;

        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            numBytesParsed += sizeof(uint16_t);
            pTagData       += sizeof(uint16_t);
            break;

        case (Constants::BAM_TAG_TYPE_FLOAT)  :
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
            numBytesParsed += sizeof(uint32_t);
            pTagData       += sizeof(uint32_t);
            break;

        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
            while( *pTagData ) {
                ++numBytesParsed;
                ++pTagData;
            }
            // increment for null-terminator
            ++numBytesParsed;
            ++pTagData;
            break;

        case (Constants::BAM_TAG_TYPE_ARRAY) :

        {
            // read array type
            const char arrayType = *pTagData;
            ++numBytesParsed;
            ++pTagData;

            // read number of elements
            int32_t numElements;
            memcpy(&numElements, pTagData, sizeof(uint32_t)); // already endian-swapped, if needed
            numBytesParsed += sizeof(uint32_t);
            pTagData       += sizeof(uint32_t);

            // calculate number of bytes to skip
            int bytesToSkip = 0;
            switch (arrayType) {
                case (Constants::BAM_TAG_TYPE_INT8)  :
                case (Constants::BAM_TAG_TYPE_UINT8) :
                    bytesToSkip = numElements;
                    break;
                case (Constants::BAM_TAG_TYPE_INT16)  :
                case (Constants::BAM_TAG_TYPE_UINT16) :
                    bytesToSkip = numElements*sizeof(uint16_t);
                    break;
                case (Constants::BAM_TAG_TYPE_FLOAT)  :
                case (Constants::BAM_TAG_TYPE_INT32)  :
                case (Constants::BAM_TAG_TYPE_UINT32) :
                    bytesToSkip = numElements*sizeof(uint32_t);
                    break;
                default:
                    const string message = string("invalid binary array type: ") + arrayType;
                    SetErrorString("BamAlignment::SkipToNextTag", message);
                    return false;
            }

            // skip binary array contents
            numBytesParsed += bytesToSkip;
            pTagData       += bytesToSkip;
            break;
        }

        default:
            const string message = string("invalid tag type: ") + storageType;
            SetErrorString("BamAlignment::SkipToNextTag", message);
            return false;
    }

    // if we get here, tag skipped OK - return success
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

std::pair<int,int> BamAlignment::getPairedRange() const {
   if (!mateOnSameReference()) {
      return getRange();
   }
   int b, e;
   if (IsReverseStrand()) {
      b=getMatePosition();
      e=GetEndPosition(false,true);
   }
   else { // + strand
      b=getPosition();
      e=b+getInsertSize()-1;
   }
   return make_pair(b,e);
}

int BamAlignment::getPairedEndPosition() const {
   if (!mateOnSameReference()) {
      return getEndPosition();
   }
   if (IsReverseStrand()) {
      return GetEndPosition(false,true);
   }
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

// if b is after position then we need to advance both
// i and j to the starting point
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

// [b,e] are reference coordinate
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

// due to redundancy, extra work is needed
// bad design
void BamAlignment::chopFirstSoftclip() {
   // remove the first cigar operation
   //assert(CigarData.front().Type == 'S');
   if (CigarData.front().Type == 'S') {
      int tmplen = CigarData.front().getLength();
      QueryBases=QueryBases.substr(tmplen);
      Length -= tmplen;
      Qualities=Qualities.substr(tmplen);
      SupportData.QuerySequenceLength = Length;
      SupportData.NumCigarOperations = CigarData.size()-1;
      CigarData.erase(CigarData.begin());
   }
}

void BamAlignment::chopLastSoftclip() {
   //assert(CigarData.back().getType() == 'S');
   if (CigarData.back().getType() == 'S') {
      int tmplen = CigarData.back().getLength();
      Length -= tmplen;
      QueryBases.resize(Length);
      Qualities.resize(Length);
      SupportData.QuerySequenceLength = Length;
      SupportData.NumCigarOperations = CigarData.size()-1;
      CigarData.resize(CigarData.size()-1);
   }
}

void BamAlignment::updateNMTag(const string& refseq) {
   int b = getPosition();
   int e = GetEndPosition(); // one passed the end [b,e)
   string subseq = refseq.substr(b, e-b);
   int edit=0;
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
   if (HasTag("NM")) {
      int tmp;
      if (!GetTag("NM", tmp)) {
         cerr << "Failed to get NM tag!\n";
         exit(1);
      }
      //cout << "Old edit distance: " << tmp << endl;
      if (abs(tmp - edit) > 50) {
         cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
            << " old " << tmp << " and  new " << edit << " edit distance too big check logic\n"
            << *this << endl;
         exit(1);
      }
      EditTag("NM", "i", edit);
   }
   else {
      AddTag("NM", "i", edit);
   }
}

bool BamAlignment::hasSoftclip() const {
   if (CigarData.empty()) return false;
   return CigarData.front().Type == 'S' || CigarData.back().Type == 'S';
}

int BamAlignment::getASValue() const {
   if (!HasTag("AS")) return -1;
   int val=-1;
   if (!GetTag("AS", val)) {
      cerr << "Failed to get AS tag vaule returning -1\n";
      return -1;
   }
   return val;
}

int BamAlignment::getNMValue() const {
   if (!HasTag("NM")) return -1;
   int val=-1;
   if (!GetTag("NM", val)) {
      cerr << "Failed to get NM value returning -1\n";
      return -1;
   }
   return val;
}

int BamAlignment::getTemplateLength() const {
   int tlen = getInsertSize();
   if (tlen == 0 && HasTag("XO")) {
      tlen = getReferenceWidth();
   }
   return abs(tlen);
}

