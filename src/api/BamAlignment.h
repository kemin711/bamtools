// ***************************************************************************
// BamAlignment.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 July 2013 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#ifndef BAMALIGNMENT_H
#define BAMALIGNMENT_H

#include "api/api_global.h"
#include "api/BamAux.h"
#include "api/BamConstants.h"
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <utility>

namespace BamTools {

//! \cond
// forward declaration of BamAlignment's "friends"
namespace Internal {
    class BamReaderPrivate;
    class BamWriterPrivate;
} // namespace Internal
//! \endcond

// BamAlignment data structure
class API_EXPORT BamAlignment {
   // API_EXPORT are constructed used for Window DDL

    // constructors & destructor
    public:
        /** 
         *   Default constructor
         */
        BamAlignment(void);
        /** 
         *  Copy constructor
         */
        BamAlignment(const BamAlignment& other);
        /**
         * Move constructor
         */
        BamAlignment(BamAlignment&& other);
        ~BamAlignment(void);
        /**
         * Test version for looking at the critical information about the
         * alignment.  Should output a human readable tab-dlimited forms with
         * most of the essential information about an alignment.
         */
        friend std::ostream& operator<<(std::ostream& ous, const BamAlignment& ba);
        /**
         * Assignment operator
         */
        BamAlignment& operator=(const BamAlignment& other);
        /**
         * move assignment operator
         */
        BamAlignment& operator=(BamAlignment&& other);

    // queries against alignment flags
    // // BAM alignment flags
    // const int BAM_ALIGNMENT_PAIRED              = 0x0001;
    // const int BAM_ALIGNMENT_PROPER_PAIR         = 0x0002;
    // const int BAM_ALIGNMENT_UNMAPPED            = 0x0004;
    // const int BAM_ALIGNMENT_MATE_UNMAPPED       = 0x0008;
    // const int BAM_ALIGNMENT_REVERSE_STRAND      = 0x0010;
    // const int BAM_ALIGNMENT_MATE_REVERSE_STRAND = 0x0020;
    // const int BAM_ALIGNMENT_READ_1              = 0x0040;
    // const int BAM_ALIGNMENT_READ_2              = 0x0080;
    // const int BAM_ALIGNMENT_SECONDARY           = 0x0100;
    // const int BAM_ALIGNMENT_QC_FAILED           = 0x0200;
    // const int BAM_ALIGNMENT_DUPLICATE           = 0x0400;
    //
    // The naming in this is different from the Bam documentation
    //                     Description
    // 1    0x1   template having multiple segments in sequencing (paired)
    // 2    0x2   each segment properly aligned according to the aligner
    // 4    0x4   segment unmapped
    // 8    0x8   next segment in the template unmapped
    // 16   0x10  SEQ being reverse complemented
    // 32   0x20  SEQ of the next segment in the template being reverse complemented
    // 64   0x40  the first segment in the template (first mate or read)
    // 128  0x80  the last segment in the template (second mate or read)
    // 256  0x100 secondary alignment
    // 512  0x200 not passing filters, such as platform/vendor quality controls
    // 1024 0x400 PCR or optical duplicate
    // 2048 0x800 supplementary alignment
    // the following tests the flag field against different bits
    public:        
        bool IsDuplicate(void) const;         // returns true if this read is a PCR duplicate
        bool IsFailedQC(void) const;          // returns true if this read failed quality control
        /**
         *  @returns true if alignment is first mate on paired-end read
         */
        bool IsFirstMate(void) const;         
        bool isFirstMate(void) const { return AlignmentFlag & 0x40; }
        /** 
         * @returns true if alignment is second mate on paired-end read
         */
        bool IsSecondMate(void) const;        
        bool isSecondMate(void) const { return AlignmentFlag & 0x80; }
        int getMate() const { if (isFirstMate()) return 1; 
           else if (isSecondMate()) return 2;
           else return 0;
        }
        bool IsMapped(void) const;            // returns true if alignment is mapped
        bool IsMateMapped(void) const;        // returns true if alignment's mate is mapped
        bool IsMateReverseStrand(void) const; // returns true if alignment's mate mapped to reverse strand
        bool IsPaired(void) const;            // returns true if alignment part of paired-end read
        bool IsPrimaryAlignment(void) const;  // returns true if reported position is primary alignment
        /**
         * @return true if is secondary alignment
         */
        bool isSecondaryAlignment() const { return AlignmentFlag & 0x100; }
        bool isSupplementaryAlignment() const { return AlignmentFlag & 0x800; }
        bool IsProperPair(void) const;        // returns true if alignment is part of read that satisfied paired-end resolution
        bool IsReverseStrand(void) const;     // returns true if alignment mapped to reverse strand

    // manipulate alignment flags
    public:        
        void SetIsDuplicate(bool ok);         // sets value of "PCR duplicate" flag
        void SetIsFailedQC(bool ok);          // sets value of "failed quality control" flag
        void SetIsFirstMate(bool ok);         // sets value of "alignment is first mate" flag
        void SetIsMapped(bool ok);            // sets value of "alignment is mapped" flag
        void SetIsMateMapped(bool ok);        // sets value of "alignment's mate is mapped" flag
        void SetIsMateReverseStrand(bool ok); // sets value of "alignment's mate mapped to reverse strand" flag
        void SetIsPaired(bool ok);            // sets value of "alignment part of paired-end read" flag
        void SetIsPrimaryAlignment(bool ok);  // sets value of "position is primary alignment" flag
        void SetIsProperPair(bool ok);        // sets value of "alignment is part of read that satisfied paired-end resolution" flag
        void SetIsReverseStrand(bool ok);     // sets value of "alignment mapped to reverse strand" flag
        void SetIsSecondMate(bool ok);        // sets value of "alignment is second mate on read" flag

        // convenient constants
       static const int PAIRED              = 0x0001;
       static const int PROPER_PAIR         = 0x0002;
       static const int UNMAPPED            = 0x0004;
       static const int MATE_UNMAPPED       = 0x0008;
       static const int REVERSE_STRAND      = 0x0010;
       static const int MATE_REVERSE_STRAND = 0x0020;
       static const int READ_1              = 0x0040;
       static const int READ_2              = 0x0080;
       static const int SECONDARY           = 0x0100;
       static const int QC_FAILED           = 0x0200;
       static const int DUPLICATE           = 0x0400;
       static const int SUPPLEMENTARY       = 0x0800; 

    // tag data access methods
    public:

      /** 
       * \brief Adds a field to the BAM tags.

       * Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

       * @param[in] tag   2-character tag name
       * @param[in] type  1-character tag type
       * @param[in] value data to store
       * @return \c true if the \b new tag was added successfully
       * @see \samSpecURL for more details on reserved tag names, supported tag types, etc.
      */
        template<typename T> bool AddTag(const std::string& tag, const std::string& type, const T& value);
        template<typename T> bool AddTag(const std::string& tag, const std::vector<T>& values);

        // edit (or append) tag
        /** 
         *  \brief Edits a BAM tag field.
         *
         *  If \a tag does not exist, a new entry is created.
         *
         *  @param tag[in]   2-character tag name
         *  @param type[in]  1-character tag type (must be "Z" or "H")
         *     Z for string. H for byte array in Hex format.
         *  @param value[in] new data value
         *
         *  @return \c true if the tag was modified/created successfully
         *
         *  @see BamAlignment::RemoveTag()
         *  @see \samSpecURL for more details on reserved tag names, supported tag types, etc.
        */
        template<typename T> bool EditTag(const std::string& tag, const std::string& type, const T& value);
        template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values);

        // retrieves tag data
        /** 
            Retrieves the value associated with a BAM tag.

            @param tag[in]          2-character tag name
            @param destination[out] retrieved value
            @return true if found
        */
        template<typename T> bool GetTag(const std::string& tag, T& destination) const;
        template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const;

        // retrieves all current tag names
        std::vector<std::string> GetTagNames(void) const;

        // retrieves the SAM/BAM type-code for requested tag name
        bool GetTagType(const std::string& tag, char& type) const;

        // retrieves the SAM/BAM type-code for the data elements in an array tag
        bool GetArrayTagType(const std::string& tag, char& type) const;

        // returns true if alignment has a record for this tag name
        bool HasTag(const std::string& tag) const;

        // removes a tag
        void RemoveTag(const std::string& tag);

    // additional methods
    public:
        // populates alignment string fields
        bool BuildCharData(void);
        /** 
         *  Calculates alignment end position, based on its starting position and CIGAR data.
         * 
         *  @warning The position returned now represents a zero-based, HALF-OPEN interval.
         *  In previous versions of BamTools (0.x & 1.x) all intervals were treated
         *  as zero-based, CLOSED.
         *
         *  @param[in] usePadded      Allow inserted bases to affect the reported position. Default is
         *                               false, so that reported position stays synced with reference
         *                               coordinates.
         *  @param[in] closedInterval Setting this to true will return a 0-based end coordinate. Default is
         *                               false, so that his value represents a standard, half-open interval.
         * 
         *  @return alignment end position
         */
        int GetEndPosition(bool usePadded = false, bool closedInterval = false) const;
        /**
         * return the [start, end] range of the mapping 
         * of reads on the reference.
         */
        std::pair<int,int> getRange() const { return std::pair<int,int>(getPosition(), GetEndPosition(false, true)); }

        // returns a description of the last error that occurred
        std::string GetErrorString(void) const;

        // retrieves the size, read locations and reference locations of soft-clip operations
        bool GetSoftClips(std::vector<int>& clipSizes,
                          std::vector<int>& readPositions,
                          std::vector<int>& genomePositions,
                          bool usePadded = false) const;

    // getter methods, adding these getter methods for security
    // using public member is a security issue and may not pass FDA
    // quality
        /**
         * get the name of the query squence
         */
        const std::string& getQueryName() const { return Name; }
        /**
         * getter method for the length of the query (read)
         */
        int32_t getQueryLength() const { return Length; }
        /**
         * 'original' sequence (contained in BAM file)
         */
        const std::string& getQueryBases() const { return QueryBases; }
        /**
         * 'aligned' sequence (QueryBases plus deletion, padding, clipping chars)
         */
        const std::string& getAlignedQueryBases() const { return AlignedBases; }
        /** 
         * @return the FASTQ qualities (ASCII characters, not numeric values)
         * Values are ASCII 33-93
         */
        std::string getQuality() const { return Qualities; }
        /**
         * @return the fastq quality as integer value from 0 to 63
         * This is the Phred score after ASCII - 33
         */
        vector<int> getQualityScore() const;
        /**
         * @return ID number for reference sequence
         * use this id and RefVector to get reference name.
         * The RefVector is only available in BamReader's
         * header section.
         */
        int32_t getReferenceId() const { return RefID; }
        /**
         * get the fist mapping position in 0-based index on the reference
         * sequence
         * Not sure what happens with soft clips at the beginning.
         */
        int32_t getPosition() const { return Position; }
        int16_t getMapQuality() const { return MapQuality; }
        int32_t getMateReferenceId() const { return MateRefID; }
        int32_t getMatePosition() const { return MatePosition; }
        int32_t getInsertSize() const { return InsertSize; }
        /**
         * @return a const reference to the CIGAR operations for this alignment
         */
        const std::vector<CigarOp>& getCigar() const { return CigarData; } 
        /**
         * Provide a more user-friendly interface
         * for working with other applications.
         */
        vector<pair<char,int> > getCigarOperation() const;
        /**
         * Some alignment's cigar entry is *
         * this is the same as no cigar. This function test this situation
         * @return true if no cigar
         */
        bool lackCigar() const { return CigarData.empty(); }
        /**
         * To fix certain aligner's tendency to put two gap
         * when a small region of the sequence has more mismatches
         */
        void fixStaggerGap();

        /// setter methods
        void setQueryName(const std::string &qname) { Name = qname; }
        void setQuerySequenceLength(int32_t qlen) {  Length = qlen; }
        void setQueryBases(const std::string &qseq) { QueryBases = qseq; }
        void setAlignedBases(const std::string &alnedseq) { AlignedBases = alnedseq; }
        /** set quality from string data */
        void setQuality(const std::string &qual) { Qualities = qual; }
        /** integer version
         * @param qual is a vector of Phred scores from 0-63
         */
        void setQuality(const std::vector<int> &qual);
        void setRefID(int32_t refid) { RefID = refid; }
        void setPosition(int32_t alnstart) { Position = alnstart; }
        /** alias for setPosition */
        void setStart(int32_t alnstart) { Position = alnstart; }
        void setBin(uint16_t indexbin) { Bin = indexbin; }
        void setMapQuality(uint16_t mqual) { MapQuality = mqual; }
        void setCigarData(const std::vector<CigarOp> &cd) { CigarData = cd; } 
        /** set cigarop from a vector a pair */
        void setCigarOperation(const std::vector<pair<char,int> > &cd); 
        /**
         * @param materefid Mate reference id, set to -1 if mate unmapped
         */
        void setMateRefID(int32_t materefid)  {  MateRefID = materefid; } 
        void setMatePosition(int32_t matepos) { MatePosition = matepos; } 
        void setInsertSize(int32_t insize) { InsertSize = insize; }  

    // public data fields, these fileds should all become private in the future
    public:
        /** read or query name */
        std::string Name;    
        int32_t     Length;             // length of query sequence
        std::string QueryBases;         // 'original' sequence (contained in BAM file)
        /** not sure this is useful */
        std::string AlignedBases;       // 'aligned' sequence (QueryBases plus deletion, padding, clipping chars)
        std::string Qualities;          // FASTQ qualities (ASCII characters, not numeric values)
        std::string TagData;            // tag data (use provided methods to query/modify)
        int32_t     RefID;              // ID number for reference sequence
        int32_t     Position;           // position (0-based) where alignment starts
        uint16_t    Bin;                // BAM (standard) index bin number for this alignment
        uint16_t    MapQuality;         // mapping quality score
        /**
         * This is sam/bam file field #2 containing 12 bit information
         * alignment bit-flag. 
         * use the provided methods to query/modify.
         */
        uint32_t    AlignmentFlag;      // alignment bit-flag (use provided methods to query/modify)
       /**  
        * CIGAR operations for this alignment.
        * CigarOp has Type,Length public field
        */
        std::vector<CigarOp> CigarData; 
        int32_t     MateRefID;          // ID number for reference sequence where alignment's mate was aligned
        int32_t     MatePosition;       // position (0-based) where alignment's mate starts
        int32_t     InsertSize;         // mate-pair insert size
        // alignment should not store its file name
        // information repetation, remove in future version
        // TODO: remove in next release
        std::string Filename;           // name of BAM file which this alignment comes from

    //! \internal
    // internal utility methods
    private:
        bool FindTag(const std::string& tag,
                     char*& pTagData,
                     const unsigned int& tagDataLength,
                     unsigned int& numBytesParsed) const;
        bool IsValidSize(const std::string& tag, const std::string& type) const;
        void SetErrorString(const std::string& where, const std::string& what) const;
        bool SkipToNextTag(const char storageType,
                           char*& pTagData,
                           unsigned int& numBytesParsed) const;

    // internal data
    private:

        // nested class TODO: simplify in future versions
        struct BamAlignmentSupportData {
            // data members
            std::string AllCharData;
            uint32_t    BlockLength;  // not sure what this is
            uint32_t    NumCigarOperations; // should be calculated on the fly
            uint32_t    QueryNameLength;  // duplicate data discard in the future
            uint32_t    QuerySequenceLength; // is this duplicate of QueryLength?
            bool        HasCoreOnly;
            
            // constructor
            BamAlignmentSupportData(void)
                : BlockLength(0)
                , NumCigarOperations(0)
                , QueryNameLength(0)
                , QuerySequenceLength(0)
                , HasCoreOnly(false)
            { }
            BamAlignmentSupportData(const BamAlignmentSupportData& o) 
               : AllCharData(o.AllCharData), BlockLength(o.BlockLength),
                 NumCigarOperations(o.NumCigarOperations),
                 QueryNameLength(o.QueryNameLength),
                 QuerySequenceLength(o.QuerySequenceLength),
                 HasCoreOnly(o.HasCoreOnly) { }
            BamAlignmentSupportData(BamAlignmentSupportData&& o) 
               : AllCharData(std::move(o.AllCharData)), BlockLength(o.BlockLength),
                 NumCigarOperations(o.NumCigarOperations),
                 QueryNameLength(o.QueryNameLength),
                 QuerySequenceLength(o.QuerySequenceLength),
                 HasCoreOnly(o.HasCoreOnly) { }
            BamAlignmentSupportData& operator=(const BamAlignmentSupportData& o) {
               if (this != &o) {
                 AllCharData=o.AllCharData;
                 BlockLength=o.BlockLength;
                 NumCigarOperations=o.NumCigarOperations;
                 QueryNameLength=o.QueryNameLength;
                 QuerySequenceLength=o.QuerySequenceLength;
                 HasCoreOnly=o.HasCoreOnly; 
               }
               return *this;
            }
            BamAlignmentSupportData& operator=(BamAlignmentSupportData&& o) {
               if (this != &o) {
                 AllCharData=std::move(o.AllCharData);
                 BlockLength=o.BlockLength;
                 NumCigarOperations=o.NumCigarOperations;
                 QueryNameLength=o.QueryNameLength;
                 QuerySequenceLength=o.QuerySequenceLength;
                 HasCoreOnly=o.HasCoreOnly; 
               }
               return *this;
            }
        };
        BamAlignmentSupportData SupportData;
        friend class Internal::BamReaderPrivate;
        friend class Internal::BamWriterPrivate;

        mutable std::string ErrorString; // mutable to allow updates even in logically const methods
    //! \endinternal
};

// ---------------------------------------------------------
// BamAlignment tag access methods

template<typename T>
inline bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const T& value) {

    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // check tag/type size
    if ( !IsValidSize(tag, type) ) {
        // TODO: set error string?
        return false;
    }

    // check that storage type code is OK for T
    if ( !TagTypeHelper<T>::CanConvertTo(type.at(0)) ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // otherwise, convert value to string
    union { T value; char valueBuffer[sizeof(T)]; } un;
    un.value = value;

    // copy original tag data to temp buffer
    const std::string newTag = tag + type;
    const size_t newTagDataLength = tagDataLength + newTag.size() + sizeof(T); // leave room for new T
    RaiiBuffer originalTagData(newTagDataLength);
    memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term

    // append newTag
    strcat(originalTagData.Buffer + tagDataLength, newTag.data());
    memcpy(originalTagData.Buffer + tagDataLength + newTag.size(), un.valueBuffer, sizeof(T));

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer;
    TagData.assign(newTagData, newTagDataLength);
    return true;
}

template<>
inline bool BamAlignment::AddTag<std::string>(const std::string& tag,
                                              const std::string& type,
                                              const std::string& value)
{
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // check tag/type size
    if ( !IsValidSize(tag, type) ) {
        // TODO: set error string?
        return false;
    }

    // check that storage type code is OK for string
    if ( !TagTypeHelper<std::string>::CanConvertTo(type.at(0)) ) {
        // TODO: set error string?
        return false;
    }

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // otherwise, copy tag data to temp buffer
    const std::string newTag = tag + type + value;
    const size_t newTagDataLength = tagDataLength + newTag.size() + 1; // leave room for null-term
    RaiiBuffer originalTagData(newTagDataLength);
    memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term

    // append newTag (removes original null-term, then appends newTag + null-term)
    strcat(originalTagData.Buffer + tagDataLength, newTag.data());

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer;
    TagData.assign(newTagData, newTagDataLength);
    return true;
}

/*! \fn template<typename T> bool AddTag(const std::string& tag, const std::vector<T>& values)
    \brief Adds a numeric array field to the BAM tags.

    Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.

    \param[in] tag    2-character tag name
    \param[in] values vector of data values to store
    \return \c true if the \b new tag was added successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
template<typename T>
inline bool BamAlignment::AddTag(const std::string& tag, const std::vector<T>& values) {

    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // check for valid tag name length
    if ( tag.size() != Constants::BAM_TAG_TAGSIZE )
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // build new tag's base information
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    memcpy( newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE );
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY;
    newTagBase[3] = TagTypeHelper<T>::TypeCode();

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    memcpy(newTagBase + 4, &numElements, sizeof(int32_t));

    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const size_t newTagDataLength = tagDataLength +
                                    Constants::BAM_TAG_ARRAYBASE_SIZE +
                                    numElements*sizeof(T);
    RaiiBuffer originalTagData(newTagDataLength);
    memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term

    // write newTagBase (removes old null term)
    strcat(originalTagData.Buffer + tagDataLength, (const char*)newTagBase);

    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for ( int i = 0 ; i < numElements; ++i ) {
        const T& value = values.at(i);
        memcpy(originalTagData.Buffer + elementsBeginOffset + i*sizeof(T), &value, sizeof(T));
    }

    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer;
    TagData.assign(newTagData, newTagDataLength);
    return true;
}

template<typename T>
inline bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const T& value) {

    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // remove existing tag if present, then append tag with new value
    if ( HasTag(tag) )
        RemoveTag(tag);
    return AddTag(tag, type, value);
}

/*! \fn template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values)
    \brief Edits a BAM tag field containing a numeric array.

    If \a tag does not exist, a new entry is created.

    \param tag[in]   2-character tag name
    \param value[in] vector of data values

    \return \c true if the tag was modified/created successfully
    \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
*/
template<typename T>
inline bool BamAlignment::EditTag(const std::string& tag, const std::vector<T>& values) {

    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();

    // remove existing tag if present, then append tag with new values
    if ( HasTag(tag) )
        RemoveTag(tag);
    return AddTag(tag, values);
}


template<typename T>
inline bool BamAlignment::GetTag(const std::string& tag, T& destination) const {

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
    // return failure if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }
    // fetch data type
    const char type = *(pTagData - 1);
    if ( !TagTypeHelper<T>::CanConvertFrom(type) ) {
        // TODO: set error string ?
        return false;
    }

    // determine data length
    int destinationLength = 0;
    switch ( type ) {
        // 1 byte data
        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            destinationLength = 1;
            break;
        // 2 byte data
        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            destinationLength = 2;
            break;
        // 4 byte data
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
            destinationLength = 4;
            break;
        // var-length types not supported for numeric destination
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            SetErrorString("BamAlignment::GetTag",
                           "cannot store variable length tag data into a numeric destination");
            return false;
        // unrecognized tag type
        default:
            const std::string message = std::string("invalid tag type: ") + type;
            SetErrorString("BamAlignment::GetTag", message);
            return false;
    } // using a string specialization version for string data
    // store data in destination
    destination = 0;
    memcpy(&destination, pTagData, destinationLength);
    // return success
    return true;
}

template<>
inline bool BamAlignment::GetTag<std::string>(const std::string& tag,
                                              std::string& destination) const
{
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
    // return failure if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // otherwise copy data into destination
    const unsigned int dataLength = strlen(pTagData);
    destination.clear();
    destination.resize(dataLength);
    memcpy( (char*)destination.data(), pTagData, dataLength );

    // return success
    return true;
}

/*! \fn template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const
    \brief Retrieves the numeric array associated with a BAM tag.

    \param tag[in]          2-character tag name
    \param destination[out] retrieved values
    \return \c true if found
*/
template<typename T>
inline bool BamAlignment::GetTag(const std::string& tag, std::vector<T>& destination) const {

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
    // return false if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }
    // check that tag is array type
    const char tagType = *(pTagData - 1);
    if ( tagType != Constants::BAM_TAG_TYPE_ARRAY ) {
        SetErrorString("BamAlignment::GetTag", "cannot store a non-array tag in array destination");
        return false;
    }
    // fetch element type
    const char elementType = *pTagData;
    if ( !TagTypeHelper<T>::CanConvertFrom(elementType) ) {
        // TODO: set error string ?
        return false;
    }
    ++pTagData;

    // calculate length of each element in tag's array
    int elementLength = 0;
    switch ( elementType ) {
        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            elementLength = sizeof(uint8_t);
            break;
        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            elementLength = sizeof(uint16_t);
            break;
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
            elementLength = sizeof(uint32_t);
            break;
        // var-length types not supported for numeric destination
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            SetErrorString("BamAlignment::GetTag",
                           "invalid array data, variable-length elements are not allowed");
            return false;
        // unknown tag type
        default:
            const std::string message = std::string("invalid array element type: ") + elementType;
            SetErrorString("BamAlignment::GetTag", message);
            return false;
    }
    // not using elementLength
    // TODO: remove above code block
    std::cerr << "BamTagData arrary element data width (Byte): " << elementLength << std::endl;

    // get number of elements
    int32_t numElements;
    memcpy(&numElements, pTagData, sizeof(int32_t));
    pTagData += 4;
    destination.clear();
    destination.reserve(numElements);

    // read in elements
    T value;
    for ( int i = 0 ; i < numElements; ++i ) {
        memcpy(&value, pTagData, sizeof(T));
        pTagData += sizeof(T);
        destination.push_back(value);
    }

    // return success
    return true;
}

typedef std::vector<BamAlignment> BamAlignmentVector;

} // namespace BamTools

#endif // BAMALIGNMENT_H
