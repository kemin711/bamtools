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
#include <typeinfo>
#include <map>
//#include <mutex>

// to turn off debug mode uncomment the following line
//#define NDEBUG 
//above will disable assert()
namespace BamTools {

//! \cond
// forward declaration of BamAlignment's "friends"
namespace Internal {
    class BamReaderPrivate;
    class BamWriterPrivate;
} // namespace Internal
//! \endcond

class BamAlignmentException : public exception {
   public:
      BamAlignmentException(const string& msg) {
         message=msg;
      }

      virtual const char* what() noexcept {
         return message.c_str();
      }

   private:
      string message;

};

/**
 *  BamAlignment data structure to hold a single
 *  Alignment data with summary information to it partner 
 *  (mate or read2) alignment data.
 */
class API_EXPORT BamAlignment {
   // API_EXPORT are constructed used for Window DDL
    // constructors & destructor
    public:
        /** 
         *   Default constructor
         */
        BamAlignment()
            : Name(), QueryBases(), AlignedBases(), Qualities(), 
            TagData(), RefID(-1), Position(-1), Bin(0), MapQuality(0),
            AlignmentFlag(0), CigarData(), MateRefID(-1), MatePosition(-1), 
            InsertSize(0), SupportData()           
         { }

        /**
         * Convenient constructor for testing
         * No Cigar string thus no actual alignment information.
         * Also Filename is not given.
         */
        BamAlignment(const std::string& qname, int32_t refid, int32_t refpos, uint32_t alnflag, 
              int32_t mrefid, int32_t mrefpos, const std::string& queryseq, const std::string& qstring)
           : 
              Name(qname), 
             QueryBases(queryseq), AlignedBases(), Qualities(qstring),
             TagData(), RefID(refid), Position(refpos), Bin(0), MapQuality(0),
             AlignmentFlag(alnflag), CigarData(),
            MateRefID(mrefid), MatePosition(mrefpos), InsertSize(0),
            SupportData()
        {}
        /**
         * Convenient constructor with Cigar input for testing
         */
        BamAlignment(const std::string& qname, int32_t refid, int32_t refpos, uint32_t alnflag, 
              int32_t mrefid, int32_t mrefpos, const std::string& queryseq, const std::string& qstring,
              const string& cigarstr)
           : Name(qname), 
             QueryBases(queryseq), AlignedBases(), Qualities(qstring),
             TagData(), RefID(refid), Position(refpos), Bin(0), MapQuality(0),
             AlignmentFlag(alnflag), CigarData(),
            MateRefID(mrefid), MatePosition(mrefpos), InsertSize(0),
            SupportData()
        { setCigar(cigarstr); }

        /** 
         *  Copy constructor
         */
        BamAlignment(const BamAlignment& other);
        /**
         * Move constructor
         */
        BamAlignment(BamAlignment&& other);
        /**
         * This class does not have dynamicly allocated data.
         * Destructor plannng to derive from this class:
         * MultiAln, ChimearAln
         */
        virtual ~BamAlignment() { };
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
        /**
         * First compare by [begin,end], then first mate < second mate
         * Mainly used for sorting out alignments on the reference.
         */
        bool operator<(const BamAlignment& other) const;
        bool operator>(const BamAlignment& other) const;
        /**
         * Two objects are considered identical if
         * same name, same mate. Mainly using this to sort out
         * same query mapped to multiple locations.
         */
        bool operator==(const BamAlignment& other) const;
        /**
         * interface function should never be used by this class
         * Should be implemented by derived classes.
         */
        virtual bool operator()(BamAlignment* ba) { 
           if (ba != nullptr) {
              throw runtime_error("base class cannot add true pointer");
           }
            return true;
        }
        /**
         * Should be derived by subclasses
         */
        virtual bool isFull() const {
           return true;
        }
        /**
         * Should be derived by subclasses
         */
        virtual void findBestLocation() {
        }

        /**
         * This object has not been populated with real data.
         */
        bool empty() const {
           return getLength() == 0;
        }
        bool isEmpty() const { return getLength() == 0; }
        /**
         * Not sure the condition is sufficient, need some testing.
         * The caller needs to make sure the query name is the same or 
         * not. This function only cares about the reference location.
         * @return true if two aligns have the same refid, strand, Postion, endPosition
         */
        bool sameHit(const BamAlignment& other) const {
           return getReferenceId() == other.getReferenceId()
              && getStrand() == other.getStrand()
              && getPosition() == other.getPosition()
              && getEndPosition() == other.getEndPosition();
        }
        uint32_t getAlignmentFlag() const {
            return AlignmentFlag;
        }
         /** 
          * This cannot be relied on.
          * @return true if this read is a PCR duplicate
         */
        bool IsDuplicate(void) const;
        bool IsFailedQC(void) const;          // returns true if this read failed quality control
        /**
         *  @return true if alignment is first mate on paired-end read
         *  First mate is also first read: Two words meaning the same thing.
         *  First read usually have better quality than second read.
         */
        bool IsFirstMate(void) const;         
        /**
         * @return true if is the First Read in a paired end read.
         * Alias for IsFirstMate.
         * Use this version for carmel casing.
         * @see isFirstRead
         */
        bool isFirstMate(void) const { 
           return (AlignmentFlag & READ_1) == READ_1; 
        }
        /**
         * Alias for isFirstMate()
         */
        bool isFirstRead(void) const { 
           return (AlignmentFlag & READ_1) == READ_1; 
        }
        /** 
         * @returns true if alignment is second mate on paired-end read
         */
        bool IsSecondMate(void) const;        
        /**
         * @return true if it is the second read (mate)
         */
        bool isSecondMate(void) const { 
           return (AlignmentFlag & READ_2) == READ_2; 
        }
        /**  in C++ true is 1 false is 0
         * Alias for isSecondMate()
         */
        bool isSecondRead(void) const { 
           return (AlignmentFlag & READ_2) == READ_2; 
        }
        /**
         * @return 1 for first mate, 2 for second mate, 
         *     and 0 for unknown mate or not pair-end read.
         */
        int getMate() const { 
           if (isFirstMate()) return 1; 
           else if (isSecondMate()) return 2;
           else return 0;
        }
        /**
         * @return 0 if is unpaired or one of the mate is unmapped
         *    otherwise return 1 for read (mate) 1, 2 for read (mate) 2
         */
        int getMappingStatus() const {
           if (!isPaired() || isMateUnmapped() || isUnmapped()
                 || !mateOnSameReference()) return 0;
           return getMate();
        }
        /**
         * @return true if query is mapped to references
         * @see isMapped is an alias conforming to convension.
         */
        bool IsMapped(void) const;            
        /** 
         * @return true if alignment is mapped
         */
        bool isMapped(void) const {
           return !((AlignmentFlag & UNMAPPED) == UNMAPPED);            
        }
        /** 
         * @return true if alignment is not mapped
         */
        bool isUnmapped() const {
           return (AlignmentFlag & UNMAPPED) == UNMAPPED;
        }
        /** 
         * Is a flag check
         * @return true if alignment's mate is mapped
         */
        bool IsMateMapped(void) const;        
        bool isMateMapped(void) const {
            return !((AlignmentFlag & MATE_UNMAPPED) == MATE_UNMAPPED);
        }
        /**
         * @return true if Mate is unmapped
         */
        bool isMateUnmapped() const {
           return (AlignmentFlag & MATE_UNMAPPED) == MATE_UNMAPPED;
        }
        /** 
         *  If the read is mapped to the reference then there is
         *  only two choices +/-.  - is the reverse strand.
         *
         *  @return true if alignment mapped to reverse strand of refseq.
         */
        bool IsReverseStrand(void) const;     
        /**
         * @return true if alignment mapped to reverse strand of refseq.
         */
        bool isReverseStrand(void) const {
           return (AlignmentFlag & REVERSE_STRAND) == REVERSE_STRAND; 
        }
        bool isForwardStrand() const {
           return !((AlignmentFlag & REVERSE_STRAND) == REVERSE_STRAND); 
        }
        void setReverseStrand() {
           AlignmentFlag |= REVERSE_STRAND;
        }
        void setForwardStrand() {
           AlignmentFlag &= ~(REVERSE_STRAND);
        }
        /**
         * There is only two state, cannot be zero.
         * @return -1 reverse strand, +1 for forward strand, and
         *      0 for both strand or strand unknown.
         */
        int getStrand() const {
           if (isReverseStrand()) return -1;
           return 1;
           //if (isForwardStrand()) return 1;
           //return 0; // never reach this line
        }
        char getStrandChar() const {
           if (isReverseStrand()) return '-';
           return '+';
           //if (isForwardStrand()) return '+';
           //return '?'; // nover reach this line
        }
        /** 
         * @return true if alignment's mate mapped to reverse strand
         */
        bool IsMateReverseStrand(void) const; 
        bool isMateReverseStrand() const {
            return (AlignmentFlag & MATE_REVERSE_STRAND) == MATE_REVERSE_STRAND;
        }
        bool isMateForwardStrand() const {
            return !((AlignmentFlag & MATE_REVERSE_STRAND) == MATE_REVERSE_STRAND);
        }
        bool isMateOppositeStrand() const {
           return  (isForwardStrand() && isMateReverseStrand())
              || (isReverseStrand() && isMateForwardStrand()); 
        }
        bool isMateSameStrand() const {
           return  (isForwardStrand() && isMateForwardStrand())
              || (isReverseStrand() && isMateReverseStrand()); 
        }
        /**
         * Get a value to represent the strand.
         * For simple alignment of single sequence it is either
         * +1 or -1.  For Merged sequence which is reprented as
         * having the XO tag this value could be [0, 1)
         * a fractional number between 0 and 1.
         * If the overlap is 100% of the read length, then it is zero.
         * The higher the value approaching 1 the lower the overlap.
         * We use a single return value to record
         * this multiple situations.
         */
        double getFractionStrand() const;
        /** @returns true if alignment part of paired-end read
         */
        bool IsPaired(void) const; 
        /**
         * If is not paired then READ1 or READ2 is meaningless
         * @return this alignment is paired or single.
         */
        bool isPaired() const { 
            return (AlignmentFlag & PAIRED) == PAIRED;
        }
        /** 
         * @return true if reported position is primary alignment
         */
        bool IsPrimaryAlignment(void) const;  
        /**
         * Test SECONDARY flag is not set.
         */
         bool isPrimaryAlignment(void) const  {
            return !((AlignmentFlag & SECONDARY) == SECONDARY);
         }
        /**
         * @return true if is secondary alignment
         *    This indicates that the query was mapped
         *    multiple times in the genome. The choice of
         *    primary/secondary is usually arbitrary!
         */
        bool isSecondaryAlignment() const { 
           return (AlignmentFlag & SECONDARY) == SECONDARY; }
        /**
         * Alignment is part of a chimera. The choice
         * of representative/supplementary is arbitrary.
         */
        bool isSupplementaryAlignment() const { 
           return (AlignmentFlag & SUPPLEMENTARY) == SUPPLEMENTARY; 
        }
        /**
         * Alias for isSupplementaryAlignment
         */
        bool isSupplementary() const { 
           return (AlignmentFlag & SUPPLEMENTARY) == SUPPLEMENTARY; 
        }
        /** 
         * @return true if alignment is part of read that satisfied paired-end resolution
         */
        bool IsProperPair(void) const;        
        /**
         * Test against PROPER_PAIR flag. This is usually set by an aligner.
         * For a negation test the result will include:
         *    1. Singleton
         *    2. pairs located to multi-references
         *    3. Inert size way out of range compared to the
         *       mean std (3 or 4 std away from mean?) This is
         *       all determined by the aligner.
         */
        bool isProperPair() const {
            return (AlignmentFlag & PROPER_PAIR) == PROPER_PAIR;
        }
        bool isNotProperPair() const {
            return !((AlignmentFlag & PROPER_PAIR) == PROPER_PAIR);
        }
        bool isImproperPair() const {
            return !((AlignmentFlag & PROPER_PAIR) == PROPER_PAIR);
        }

    // manipulate alignment flags
    public:        
        void setAlignmentFlag(uint32_t flag) {
           AlignmentFlag = flag;
        }
        /** sets value of "PCR duplicate" flag
         */
        void SetIsDuplicate(bool ok);         
        void SetIsFailedQC(bool ok);          // sets value of "failed quality control" flag
        void SetIsFirstMate(bool ok);         // sets value of "alignment is first mate" flag
        void SetIsMapped(bool ok);            // sets value of "alignment is mapped" flag
        void SetIsMateMapped(bool ok);        // sets value of "alignment's mate is mapped" flag
        /** 
         * sets value of "alignment mapped to reverse strand" flag
         */
        void SetIsReverseStrand(bool ok);     
        /**  
         * sets value of "alignment's mate mapped to reverse strand" flag
         */
        void SetIsMateReverseStrand(bool ok); 
        /** sets value of "alignment part of paired-end read" flag
         */
        void SetIsPaired(bool ok);            
        void SetIsPrimaryAlignment(bool ok);  // sets value of "position is primary alignment" flag
        /** sets value of "alignment is part of read that satisfied paired-end 
         *  resolution" flag
         */
        void SetIsProperPair(bool ok);        
        void SetIsSecondMate(bool ok);        // sets value of "alignment is second mate on read" flag

        // convenient constants octal number, first digit is 0
         static const uint32_t PAIRED              = 0x0001;
         static const uint32_t PROPER_PAIR         = 0x0002;
         static const uint32_t UNMAPPED            = 0x0004;
         static const uint32_t MATE_UNMAPPED       = 0x0008;
         static const uint32_t REVERSE_STRAND      = 0x0010;
         static const uint32_t MATE_REVERSE_STRAND = 0x0020;
         static const uint32_t READ_1              = 0x0040;
         static const uint32_t READ_2              = 0x0080;
         static const uint32_t SECONDARY           = 0x0100;
         static const uint32_t QC_FAILED           = 0x0200;
         static const uint32_t DUPLICATE           = 0x0400;
         static const uint32_t SUPPLEMENTARY       = 0x0800; 

    // tag data access methods
    public:
        /** 
         * Adds a field to the BAM tags.
         * Does NOT modify an existing tag - use \link BamAlignment::EditTag() 
         *    \endlink instead.
         * @param[in] tag   2-character tag name
         * @param[in] type  1-character tag type such as Z, H, or i
         * @param[in] value data to store
         * @return \c true if the \b new tag was added successfully
         * @see \samSpecURL for more details on reserved tag names, 
         *     supported tag types, etc.
         */
        template<typename T> bool AddTag(const std::string& tag, const std::string& type, const T& value);
        /**
         * Adds a numeric array field to the BAM tags.
         *
         *  Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.
         *
         *  \param[in] tag    2-character tag name
         *  \param[in] values vector of data values to store
         *  \return \c true if the \b new tag was added successfully
         *  \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
         */
        template<typename T> bool AddTag(const std::string& tag, const std::vector<T>& values);
        /** 
         *  edit (or append) tag
         *  \brief Edits a BAM tag field.
         *
         *  If \a tag does not exist, a new entry is created.
         *
         *  @param tag[in]   2-character tag name
         *  @param type[in]  1-character tag type (must be "Z", "i", or "H")
         *     Z for string. H for byte array in Hex format.
         *  @param value[in] new data value
         *
         *  @return \c true if the tag was modified/created successfully
         *
         *  @see BamAlignment::RemoveTag()
         *  @see \samSpecURL for more details on reserved tag names, supported tag types, etc.
        */
        template<typename T> bool EditTag(const std::string& tag, const std::string& type, const T& value);
        /**
         *  Edits a BAM tag field containing a numeric array.
         *
         *  If \a tag does not exist, a new entry is created.
         *
         *  \param tag[in]   2-character tag name
         *  \param value[in] vector of data values
         *
         *  \return \c true if the tag was modified/created successfully
         *  \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
         */
        template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values);

        // retrieves tag data
        /** 
         *  Retrieves the value associated with a BAM tag.
         *  There is no excetion, when the function return false
         *  it means no such tag.
         *
         *  @param tag[in]          2-character tag name
         *  @param destination[out] retrieved value, will clear previous values.
         *  @return true if found.
         *
        */
        template<typename T> bool GetTag(const std::string& tag, T& destination) const;
	/**
	 * @param destination will clear previous values and fill with tag value.
	 * @return true if found
	 */
        template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const;
        /**
         * More convenient version return value directly. If not found
         * for string then return empty string. For integer types
         * will return the smallest integer.
         * TODO: Need testing
         * More convenient method than the reference version
         * @return the tag string value if exists otherwise return empty string
         */
        string getStringTag(const std::string& tag) const;
        /** 
         * retrieves all current tag names
         *
         * When paired with GetTagType() and GetTag(), this method allows you
         * to iterate over an alignment's tag data without knowing the names (or types)
         * beforehand.

         * @return \c vector containing all tag names found (empty if none available)
         * \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
         */
        std::vector<std::string> GetTagNames(void) const;
        /** 
         * retrieves the SAM/BAM type-code for requested tag name
         * @param[in]  tag  2-character tag name
         * @param[out] type retrieved (1-character) type-code
         * @return \c true if found
         * \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
         */
        bool GetTagType(const std::string& tag, char& type) const;
        // retrieves the SAM/BAM type-code for the data elements in an array tag
        bool GetArrayTagType(const std::string& tag, char& type) const;
        /** 
         * @returns true if alignment has a record for this tag name
         */
        bool HasTag(const std::string& tag) const;
        /** 
         * removes a tag
         */
        void RemoveTag(const std::string& tag);

    // additional methods
    public:
        // populates alignment string fields
        bool BuildCharData(void);
        /** 
         *  Calculates alignment end position, based on its starting 
         *  position and CIGAR data.
         * 
         *  @warning The position returned now represents a zero-based,
         *      HALF-OPEN interval.  In previous versions of BamTools 
         *      (0.x & 1.x) all intervals were treated as zero-based, CLOSED.
         *  @param[in] usePadded Allow inserted bases to affect the reported
         *      position. Default is false, so that reported position
         *      stays synced with reference coordinates.
         *  @param[in] closedInterval Setting this to true will return a
         *      0-based end coordinate. Default is false, so that his value
         *      represents a standard, half-open interval.
         *  @return alignment end position on the reference.
         */
        int GetEndPosition(bool usePadded = false, bool closedInterval = false) const;
        /**
         * @return end of the mapping index (0-based) on reference.
         *    This is the last mapped base index on reference. 
         *    Closed end of the interval.
         */
        int getEndPosition() const { return GetEndPosition(false, true); }
        /**
         * @return the end position of the read pair.
         *   if calling from this which is mapped to the forward
         *   strand, mate will be mapped to the reverse strand
         *   then the end of the end of mate mapping.
         *  if this read is the reverse strand, then the end is 
         *  this mapping's end.
         *  If the mate is mapped to a different chromosome
         *  then just the end position of this mapping for this read.
         */
        int getPairedEndPosition() const;
        /**
         * @return true if this alignment contains a closed end
         *   range [b,e] on the genomic reference coordinate.
         */
        bool contains(int b, int e) const {
            return getPosition() <= b && getEndPosition() >= e;
        }

        /**
         * @return the [start, end] range of the mapping 
         * of reads on the reference in 0-based index.
         * start < end
         */
        std::pair<int,int> getRange() const { 
           return std::pair<int,int>(getPosition(), GetEndPosition(false, true)); 
        }
        /**
         * the first is always smaller than the second
         * @return the [begin,end] of the closed range for the mapping of this sequence
         *    on the genomic DNA. 0-based index. 
         */
        std::pair<int,int> getInterval() const;
        bool sameInterval(const BamAlignment& ba) const;
        /**
         * @return the interval including soft-clipped region if they
         *   exist. A more generic version of getInterval() but 
         *   more costly. So only in special situation you need
         *   to use this version.
         */
        std::pair<int,int> getSoftInterval() const;
        /**
         * The distance coverted by the alignment on the reference.
         * Should be the same as the sum of cigar M+D length.
         */
        int getReferenceWidth() const {
           return GetEndPosition(false,false) - getPosition();
        }
        /**
         * For properly mapped pairs: +/- or -/+ configuration
         * <pre>
         * Get the range if it is paired on the same reference
         *    b1      e1     b2      e2
         * 1) ==this==> ---- <==mate==  return [b1, b2] we don't know e2
         *    b2      e2     b1      e1
         * 2) ==mate==> ---- <==this==   return [b2, e1] we know the whole range
         * Otherwise just the range of this alignment.
         * </pre>
         * b could be larger than e
         */
        std::pair<int,int> getPairedRange() const;
        /**
         * Same as getPairedRange
         * Should always be [small, large]
         */
        std::pair<int,int> getPairedInterval() const;

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
         * @return the name of the query.
         */
        const std::string& getName() const { return Name; }
        /**
         * Alias for shorter method getLength()
         * getter method for the length of the query sequence (read)
         * @return the query sequence length
         */
        int32_t getQueryLength() const { 
           return getLength(); 
        }
        /**
         * @return the length of the matched part of the
         *    the reference.  This is the sum of cigar
         *    string M 
         *  @see getReferenceWidth
         */
        int getMatchedReferenceLength() const;
        /**
         * @return true if More match than softclip
         */
        bool matchDominate() const {
           return getMatchedReferenceLength() > getSoftclipLength();
        }
        /**
         * 'original' sequence (contained in BAM file)
         */
        const std::string& getQueryBases() const { return QueryBases; }
        /**
         * @return the query sequence same as getQueryBases.
         * If there is Hard-clip, the only the alignmed part will be 
         * returned. So should tell BWA never use hardclip.
         */
        const std::string& getQuerySequence() const { return QueryBases; }
        std::string& getQueryBases() { return QueryBases; }
        std::string& accessSequence() {
           return QueryBases;
        }
        bool sameQuerySequence(const BamAlignment& ba) const {
           return getQueryBases() == ba.getQueryBases();
        }
        /**
         * aligned sequence (QueryBases plus deletion, padding, clipping chars)
         * Not good for sequence analysis, need a separate method to 
         * extract the aligned without indel symbols.
         * @return the formated aligned query sequence.
         *   if not present, could return empty string.
         * see: getMatchedQuerySequence()
         */
        const std::string& getAlignedQueryBases() const { 
           //if (AlignedBases.empty()) {
           //   cerr << "empty aligned bases!\n";
           //}
           return AlignedBases; 
        }
        /**
         * Will not depend on whether the AlignedBases has been formated or not.
         * @return the portion of query sequences that matched the reference
         *    without any deletion characters. QueryBases excluding the 
         *    softclip sequences if they exist.
         */
        string getMatchedQuerySequence() const;
        /**
         * alignment length - (H or S)
         * did not exclude D, could measure entire aligned part
         */
        int getMatchedQueryLength() const;
        /**
         * Query bases in alignment excluding S, H, I, and D segments.
         */
        int numberBaseAligned() const;
        void clearAlignedBases() {
           AlignedBases.clear();
        }
        /**
         * Update AlignedBases.
         * This function needs to be called if the alignment
         * part of the object is changed.
         * For simplicity, the alignment can simply be cleared.
         * Most operations will not need this part of the computed
         * result. The result can be computed from the Cigar string.
         * The designer should not have included it as a member of the 
         * object in the first place.
         */
        void setAlignedBases(const std::string &alnedseq) { AlignedBases = alnedseq; }
        /** 
         * @return the FASTQ qualities (ASCII characters, not numeric values)
         * Values are ASCII 33-93
         */
        const std::string& getQuality() const { return Qualities; }
        std::string& getQuality() { return Qualities; }
        /**
         * @return the fastq quality as integer value from 0 to 63
         * This is the Phred score after ASCII - 33
         */
        vector<int> getQualityScore() const;
        int getAverageQualityScore() const;
        /**
         * @return ID number for reference sequence
         * use this id and RefVector to get reference name.
         * The RefVector is only available in BamReader's
         * header section.
         * RefVector is vector<RefData>, RefData is { RefName, RefLength }
         * unnecessary typedef, more confusing than help.
         * RefID -> index of the RefVector
         */
        int32_t getReferenceId() const { return RefID; }
        /**
         * get the fist mapping position in 0-based index on the reference
         * sequence.
         * Soft clips at the beginning does not count, the 
         * position is that of the reference.
         */
        int32_t getPosition() const { return Position; }
        /**
         * @return mapping quality.
         */
        int16_t getMapQuality() const { return MapQuality; }
        /**
         * @return the reference if of the mate
         */
        int32_t getMateReferenceId() const { return MateRefID; }
        /**
         * @return mate mapped position
         */
        int32_t getMatePosition() const { return MatePosition; }
        /**
         * @return the strand as a single char + or -
         */
        char getMateStrandChar() const {
           if (IsMateReverseStrand()) return '-';
           else return '+';
        }

        /**
         * @return true if mate is mapped to the same reference.
         *   If no mate or mate is not mapped will return false.
         */
        bool mateOnSameReference() const { return MateRefID == RefID; }
        bool mateOnDifferentReference() const { return MateRefID != RefID; }
        /**
         * If read is unpaired, the insert size is also zero.
         * @return length of the template positive for
         *    plus strand mate and negative for minus strand mate.
         *    If the two reads map to different references then
         *    the insert size is set to zero.
         *
         *  BamWriter use the InsertSize directly for output.
         */
        int32_t getInsertSize() const { 
           return InsertSize; 
        }
        /**
         * The returned value is >= 0. not sure better
         * to use unsigned int or not.
         *
         * @return the length of the InsertSize if it is not zero
         *   otherwise return the length of the projected length
         *   of the read on the reference.
         *
         *  This function used getInsertSize()
         */
        int getTemplateLength() const;
        /**
         * @return a const reference to the CIGAR operations for this alignment
         * CigarData is a typedef of vector<CigarOp>
         * DigarOp: struc { char Type; uint32_t Length; }
         */
        const std::vector<CigarOp>& getCigar() const { return CigarData; } 
        std::vector<CigarOp>& getCigar() { return CigarData; } 
        /**
         * Provide a more user-friendly interface
         * for working with other applications.
         */
        vector<pair<char,int> > getCigarOperation() const;
        /**
         * @return a string version of Cigar.
         */
        string getCigarString() const;
        /**
         * @return true if the CigarData of this object is the same as 
         *    the argument cigar
         */
        bool sameCigar(const vector<pair<char,int>>& cigar) const;
        bool sameCigar(const BamAlignment& ba) const {
           return CigarData == ba.CigarData;
        }
        /**
         * Some alignment's cigar entry is *
         * this is the same as no cigar. This function test this situation
         * @return true if no cigar
         */
        bool lackCigar() const { return CigarData.empty(); }
        bool hasCigar() const { return !CigarData.empty(); }
        /**
         * There is no Deletion segment in the Cigar
         */
        bool lackDCigar() const;
        bool hasDCigar() const;
        /**
         * There is not Insertion segment in the Cigar
         */
        bool lackICigar() const;
        bool hasICigar() const;
        bool validCigar() const;
        /**
         * @return the  Cigar operation type as a single char
         *    at index i.
         */
        char getCigarType(unsigned int i) const {
           return CigarData[i].getType();
        }
        /**
         * @return the length of the Cigar setment at index i
         */
        unsigned int getCigarLength(unsigned int i) const {
           return CigarData[i].getLength();
        }
        unsigned int getCigarOperationCount() const {
           return CigarData.size();
        }
        /**
         * @return Number of cigar operations
         */
        int getCigarSize() const {
           return CigarData.size();
        }
        /**
         * Has indel near < 22 nt from the end
         */
        bool hasEndIndel() const;
        bool hasAmbiguousBase() const;
        /**
         * To fix certain aligner's tendency to put two gap
         * when a small region of the sequence has more mismatches
         * @deprecated use fix1M.
         */
        void fixStaggerGap();
        /**
         * Fix cigar with pattern like 57M1I1M1I69M 1M flank I or D
         */
        bool fix1M();
        /**
         * will convert 12I135M into 12S135M
         */
        void fixCigarError() {
            if (getCigarType(0) == 'I') {
               CigarData[0].setType('S');
            }
            else if (getCigarType(0) == 'S' && getCigarType(1) == 'I') {
               CigarData[0].expand(getCigarLength(1));
               CigarData.erase(CigarData.begin()+1);
            }
            if (CigarData.back().getType() == 'I') {
               CigarData.back().setType('S');
            }
            else if (CigarData.size() > 2 && CigarData.back().getType() == 'S'
                  && getCigarType(CigarData.size()-2) == 'I') {
               CigarData.back().expand(getCigarLength(CigarData.size()-2));
               CigarData.erase(CigarData.begin()+CigarData.size()-2);
            }
        }

        /**
         * @return true if start with softclip
         */
        bool startWithSoftclip() const {
            return !CigarData.empty() && getCigar().front().Type == 'S';
        }
        /**
         * @return true if ends with softclip
         */
        bool endWithSoftclip() const {
            return !CigarData.empty() && getCigar().back().Type == 'S';
        }
        /**
         * @param startR is the starting value for reference index.
         *    default is 0.  Can be any value. But the alternative
         *    meanful value is the absolute value of the genome
         *    reference such as the postion value for the alignment.
         * @return pair[cigar_index, bool]. The second value will be true
         *    if reference deisredR index has an insertion (attached to the
         *    last base of the previous M segment). The first value is the
         *    cigar segment index (0-based). If is not insertion, then, the
         *    cigar index is the one whose start index of the segment is
         *    greater than desiedR.
         */
        pair<int,bool> isInsertionAtRefloc(int desiredR, int startR=0) const;
        /**
         * @return pair[cigar_idx, bool] when desiredR is a deletion then
         *   the second value is true. The desiredR is the index of the 
         *   first base of the D segment.
         */
        pair<int,bool> isDeletionAtRefloc(int desiredR, int startR=0) const;

        /**
         * Alignment has soft clip on either start
         * or end. Not checking middle, which may not make sense.
         * @return true of has soft clip otherwise, including no cigar,
         *    then return false.
         */
        bool hasSoftclip() const;
        bool noSoftclip() const {
            if (CigarData.empty() || CigarData.size() == 1) return true;
            return CigarData.front().Type != 'S' && CigarData.back().Type != 'S';
        }
        /**
         * Calculate identity over the aligned part excluding
         * insert/delete and soft/hard clips.
         * @return ungapped identity.
         */
        float getNGIdentity() const;
        /**
         * This function requires the NM tag being present
         * NM is the mismath+del total base count.
         * @return number of mismatch, alngnment length excluding gaps
         */
        pair<int,int> getMismatchCount() const;
        /**
         * @return NM/alnlen as a fraction number that
         *   represents the local alignment identity.
         *   Soft clipped regions are not counted.
         */
        float getIdentity() const;
        /**
         * Same as horizontal coverage of the query.
         * @return the portion of aligned query sequence.
         */
        float getFractionAligned() const {
            return static_cast<float>(getMatchedQueryLength())/getLength();
        }
        /**
         * @return horizontal coverage of the query between 0 and 1
         */
        float getQCoverage() const {
            return static_cast<float>(getMatchedQueryLength())/getLength();
        }
        pair<int,int> getMatchBound() const;
        int getQueryMatchBegin() const {
           if (CigarData.front().getType() == 'S' ||
               CigarData.front().getType() == 'H')
              return CigarData.front().getLength();
           return 0;
        }
        int getQueryMatchEnd() const {
           if (CigarData.back().getType() == 'S' ||
               CigarData.back().getType() == 'H')
              return getLength() - CigarData.front().getLength() - 1;
           return getLength()-1;
        }
        /**
         * @return the query sequence for the first soft clip.
         *    If there is no soft clip then an empty string is returned.
         */
        string getFirstSoftclip() const;
        string getFirstSoftquality() const;
        /**
         * @return the length of the first softclip.
         *   if not start with softclip, then return 0.
         */
        int getFirstSoftclipLength() const;
        /**
         * @return the last soft clip in query sequence.
         */
        string getLastSoftclip() const;
        string getLastSoftquality() const;
        /**
         * @return the length of the last softclip
         *   if no softclip return 0.
         */
        int getLastSoftclipLength() const;
        /**
         * @return sum of softclip length if both are present.
         */
        int getSoftclipLength() const;
        /**
         * @return the longer of the softclip length.
         *   If no softlcip, the return 0
         */
        int getMaxSoftclipLength() const;
        /**
         * Remove the first soft clip so that the alignment
         * appears to be better. The query sequence will also 
         * be changed.
         */
        void chopFirstSoftclip();
        /**
         * Remove the last soft clip
         */
        void chopLastSoftclip();
        /**
         * Same as chopFirstSoftclip() except for checking 
         * that first soft is off the start of the reference.
         */
        void chopDangleFrontSoft();
        void chopDangleBackSoft();

        ///// setter methods ////
        /**
         * change the name of the query
         */
        void setQueryName(const std::string &qname) { Name = qname; }
        void setQuerySequenceLength(int32_t qlen) {  
           SupportData.QuerySequenceLength = qlen; 
        }
        void setQueryLength(int32_t qlen) {  
           SupportData.QuerySequenceLength = qlen; }
        /**
         * Once query bases is changed the length will also
         * change. There is no need to call setQueryLength()
         * after calling this method.
         */
        void setQueryBases(const std::string &qseq) { 
           QueryBases = qseq; 
           setQueryLength(QueryBases.size());
        }
        void setQueryBases(std::string &&qseq) { 
           QueryBases = std::move(qseq); 
           setQueryLength(QueryBases.size());
        }
        /** 
         * set quality from string data. Quality and Query sequence
         * modification should always be done in sync.
         * */
        void setQuality(const std::string &qual) { Qualities = qual; }
        void setQuality(std::string &&qual) { Qualities = std::move(qual); }
        /** 
         * integer version.
         * Since the internal representation is with char, you cannot
         * user rreference.
         * @param qual is a vector of Phred scores from 0-63
         */
        void setQuality(const std::vector<int> &qual);
        void setRefID(int32_t refid) { RefID = refid; }
        void setReferenceId(int32_t refid) { RefID = refid; }
        /**
         * Update the refgenome mapping start position
	      * If Cigar string needs to be updated then also do that.
         * Will also clear AlignedBases member.
         */
        void changePosition(int32_t alnstart); 
        /**
         * Simply update the Position member without changing anything
         * else.
         */
        void setPosition(int32_t alnstart) {
            Position = alnstart; 
        }
        /** alias for setPosition */
        void setStart(int32_t alnstart) { Position = alnstart; }
        void setBin(uint16_t indexbin) { Bin = indexbin; }
        void setMapQuality(uint16_t mqual) { MapQuality = mqual; }
        /**
         * update CigarData with a new value.
         * TODO: SupportData.NumCigarOperatons is redundant with
         *       CigarData.size() need to remove
         */
        void setCigarData(const std::vector<CigarOp> &cd) { 
           CigarData = cd; 
           SupportData.NumCigarOperations=cd.size();
        } 
        /** 
         * set CigarData from a vector a pair {char, int}
         *  for easy communication with external world.
         *  @param cd CigarData in more universal std data type
         * */
        void setCigarOperation(const std::vector<pair<char,int> > &cd); 
        void setCigar(const string& cstr);
        void setCigar(const std::vector<pair<char,int> > &cd) {
           setCigarOperation(cd);
        }
        void setCigar(const std::vector<CigarOp> &cd) {
           setCigarData(cd);
        }
        /**
         * @param materefid Mate reference id, set to -1 if mate unmapped
         */
        void setMateRefID(int32_t materefid)  {  MateRefID = materefid; } 
        /**
         * update mate position
         */
        void setMatePosition(int32_t matepos) { MatePosition = matepos; } 
        /**
         * Sets the insert size which is the 
         * length of the template.
         */
        void setInsertSize(int32_t insize) { 
           InsertSize = insize; 
        }  

        //// end of setter methods ////

        // mutation functions
        /**
         * helper function to iterate over the cigar object
         * @param i is the index on the reference. The index
         *    can be based on any value: 0-based on the first
         *    base of the subreference sequence, or the absolute
         *    index for the full chromosome. This function will
         *    make relative move.
         * @param j index on the query. Usually starting from 0
         *   for the first base of the query.
         * @param ci is the index for the cigar segment. Start
         *    from zero for the first cigar element.
         */
        void nextCigar(int& i, int& j, unsigned int& ci) const;
        int indexRef2Query(int ri) const;
        /**
         * Pick the subsequence based on the 0-based index
         * of query sequence. If the [b,e] fall in a 
         * Softclip region then the returned alignment may not
         * make sense. This function does not check for this
         * special case.
         */
        BamAlignment subsequence(int b, int e) const;
        /**
         * Use the reference (genomic) coordinate to pick subsequence
         * of the query alignment.
         * @return a new BamAlignment from b to e on the reference
         *   [b,e] is an closed end interval.
         */
        BamAlignment subsequenceByRef(int b, int e) const;
        /**
         * @return the substring of the query sequence according
         * to closed range [b,e]
         */
        std::string substringByRef(int b, int e) const;
        /**
         * @return the char at index b according to the reference
         *    sequence index coordinate system.
         */
        char charAtByRef(int ri) const;
        /**
         * @param len length of the deletion from query.
         * @return true of index is a deletion of query sequence for
         *    length of len.
         */
        bool isDeletionAt(int ri, int len) const;
        /**
         * The position is followed by a insertion of particular
         *   sequence.
         * @param seq is the sequence of last Base of M
         *   plus the entire bases of I.
         */
        bool isInsertionAt(int ri, const string& seq) const;
        /** 
         * even position number, odd position letter last 3 digits
         * for bases 00 A, 01 C, 10 G, 11 T, 100 N, first bit del
         * MD:Z:20^A127
         * MD:Z:108^TTCTAAGGCCAGCTCCTGCACC39 =>108 identical deletion of TT...ACC match of 39
         * MD tag ^AC to represent deletion of AC in query
         * But insertion in query is not recorded.
         * Difference in base is represented by Base
         * identical residues are represented by number.
        */ 
        pair<vector<int>, vector<string>> getMDArray();
        void updateMDTag(const pair<vector<int>, vector<string>>& mdvec);
        /**
         * Helper method used by trimFront()
         */
        void chopFront(size_t len, int numMismatch);
        /**
         * Helper used by trimBack()
         */
        void chopBack(size_t len, int numMismatch);
        /**
         * Remove fuzzy end from the front of the alignment
         * @return true if trimming happened.
         */
        bool trimFront();
        /**
         * Remove fuzzy end from the back of the alignment
         * @return true if trimming happened.
         */
        bool trimBack();
        /**
         * @return the trimming status for front and back.
         *   The first if for front, second for back.
         */
        pair<bool,bool> trim();
        /**
         * Move deletion from oldloc to newloc
         * the sequence in between must be repeat sequence.
         */
        void moveDeletion(int oldloc, int newloc);
        void moveInsertion(int oldloc, int newloc);
        /**
         * Only change the bases near the end < 7 nt
         * if the Query base is different from the reference base
         * on both ends. MD tag will also needs to be updated.
         * This might be a bad idea and cause trouble.
         */
        void patchEnd();
        /**
         * Recalculate the value for tag NM
         * This is needed after taking a subsequence or 
         * NM tag is missing from the object.
         * @param refseq the reference sequence whole string
         *  not just the subsequence for this object.
         */
         void updateNMTag(const string& refseq);
         /**
         * friendly wrapper for better code.
         * @return -1 if no tag
         */
         int getASValue() const;
         /**
         * @return the NM tag value or -1 if not found NM tag
         */
         int getNMValue() const;

         //// static methods ////
         static void setPolishMax(int len) {
            TRIMLEN_MAX = len;
         }
         static void setPolishGap(int gap) {
            GAP_CUT = gap;
         }
         /**
          * Helper function to be used by other functions to convert
          * string version of cigar to vector version.
          * @param cigarstr cigar string.
          * @return vector of cigar operations [code,length]
          */
         static vector<pair<char,int>> parseCigar(const string& cigarstr);
         /**
          * convert cigar data to string version.
          */
         static string cigarToString(const vector<pair<char,int>>& cg);
         /**
         * The header part of BamReader
         * for indexing readname,readlength.
         * Bam file refid, matid directly look up into this table
         * and you can get REFNAME (first) and seqlen (second)
         * You can swithc during BamAlignment processing to a different
         * file header.
         * @param refvec is the refernece vector to be obtained from 
         *    BamReader::getReferenceMetaData()
         */
         static void setRefvector(vector<pair<string,int>>&& refvec);
         /**
          * Before using this function you must load the reference look up table
          * with: setRefvector(BamReader::getReferenceMetaData())
          * You need a specific object to use getReferenceMetaData().
          * @name is the referece string name such as 1 or chr1 depends on 
          * the version of reference genome used.
          * @return the reference id given name for looking into 
          *   the id used by Bamfile.
          */
         static int referenceIdFromName(const string& name) {
            return refname2id[name];
         }
         /**
          * Given a reference id, what's the string name.
          * This is used to extract the genomic sequence
          * from external sources.
          */
         static pair<string,int> getRefnameFromId(int refid) {
            return rsname[refid];
         }

    private:
         //// static members to be shared by derived class ////
         /**
          * [refname, reflength] indexed on refid
          */
         static vector<pair<string,int>> rsname;
         static map<string,int> refname2id;
         /** 
          * Helper function apply to generic situation where
          * both Ref or Query may be advanced.
          *
          * move i to position b if i is not at b
          * j to the start of the query index
          * @param i index in reference seq 0-based index
          * @param j index in query sequence 0-based index
          * @param x desired index to advance i to on refseq. If x happen
          *     to reside in a Deletion, then x will be automatically
          *     advanced to the first base of the next M. Not sure
          *     this behavior is good or bad. Need to debug.
          * @param cigarIdx the index inside one Cigar segment
          *     range from [0, segment.length-1]
          * @param ci cigar segment index [0, NumberOfCigar-1]
          */
         //void advanceIndex(int &i, int &j, int &x, unsigned int &cigarIdx, unsigned int &ci) const;
         // above new implementation did not pass test yet
         void advanceIndex(int &i, int &j, int &x, unsigned int &cigarIdx, unsigned int &ci, char &cigarState) const;

    // public data fields, TODO: these fileds should all become private in the future
    public:
        /** 
         * read or query name 
         * Use getName() to read this one
         * */
        std::string Name;    
        int32_t     Length() const {
           return SupportData.QuerySequenceLength;
        }             
        int32_t     length() const {
           return SupportData.QuerySequenceLength;
        }             
        /**
         * @return the length of the query sequence.
         */
        int32_t getLength() const {
           assert(QueryBases.size() == SupportData.QuerySequenceLength);
           return SupportData.QuerySequenceLength;
        }
        /**
         * length setter function.
         * This should only be use in conjunction of query base 
         * insertion and deletion.
         */
        void setLength(int32_t len) {
           SupportData.QuerySequenceLength=len;
           if ((int)QueryBases.size() != len) {
              cerr << __FILE__ << ":" << __LINE__ << ":WARN QueryBase length being changed\n";
              QueryBases.resize(len);
           }
        }
        /** 'original' sequence (contained in BAM file)
         */
        std::string QueryBases;         
        /** 
         * not sure this is useful 
         * 'aligned' sequence (QueryBases plus deletion, padding, clipping
         * chars) This field will be completely empty after reading from
         * BamReader/BamMultiReader when QueryBases is empty.
         */
        std::string AlignedBases;       
        /** 
         * FASTQ qualities (ASCII characters, not numeric values)
         * String representation from char 33 to 93
         */
        std::string Qualities;          
        /** 
         * tag data (use provided methods to query/modify)
         * String encoding of everything. Should consider using
         * a proper data structure. 
         * TODO: rewrite this part.
         */
        std::string TagData;            
        /** 
         * ID number for reference sequence
         *  -1 for unmapped reads
         */
        int32_t     RefID;              
        /**
         * (0-based) where alignment starts on reference.
         */
        int32_t     Position;
        /** 
         * BAM (standard) index bin number for this alignment
         */
        uint16_t    Bin;                
        uint16_t    MapQuality;         // mapping quality score
        /**
         * This is sam/bam file field #2 containing 12 bit information
         * alignment bit-flag. 
         * use the provided methods to query/modify.
         */
        uint32_t    AlignmentFlag;
       /**  
        * CIGAR operations for this alignment.
        * CigarOp has Type,Length public field
        * Defined in BamAux.h
        */
        std::vector<CigarOp> CigarData; 
        /** ID number for reference sequence where alignment's mate was aligned */
        int32_t MateRefID;          
        /** position (0-based) where alignment's mate starts.
         *  SEELF position [a, b], MATE[c, d]
         *  If mate on the same reference, [a, b] and [c, d]
         *  could be identical.
         *  */
        int32_t MatePosition;       
        /**
         * Field 9: TLEN in bam/sam format.  signed observed Template
         *          LENgth.
         *     Observed template length equals the number of bases
         *     from the leftmost mapped base to the rightmost mapped base.
         *     It is set as zero for singl-segment template or when the
         *     information is unavailable.
         *
         * Note: even when I set this member to > 0, when the 
         * read is marked as unpaired, the BamWriter will still
         * output the filed as ZERO. Have not figured why.
         *
         * This value may be misleading in computing real fragment
         * length of the read. Need another function.
         * For both reads in the negative strand the insert size is wrong
         * <pre>
         * BWA problem:
         * <--R1-- <--R2--
         * |--I----| insert size is the number of bases in the range
         *       |--I----| should the the correct value
         * </pre>
         * Need some correction for the -/- configuration
         * Then you can get the range properly
         */
        int32_t InsertSize;        

    //! \internal
    // internal utility methods
    private:
        /** 
         *  Searches for requested tag in BAM tag data.
         *  @param  tag            requested 2-character tag name
         *  @param  pTagData       pointer to current position in BamAlignment::TagData
         *  @param  tagDataLength  length of BamAlignment::TagData
         *  @param  numBytesParsed number of bytes parsed so far before TAG
         *  @return true if found
         *  If tag is found, pTagData will point to the byte where the tag data begins.
         *        numBytesParsed will correspond to the position in the full TagData string.
        */
        bool FindTag(const std::string& tag, char*& pTagData, 
              const unsigned int& tagDataLength, unsigned int& numBytesParsed) const;
        bool IsValidSize(const std::string& tag, const std::string& type) const;
        bool SkipToNextTag(const char storageType,
                           char*& pTagData,
                           unsigned int& numBytesParsed) const;
         /**
          * Max length from either end to to polishing
          * default 6
          */
         static int TRIMLEN_MAX; // = 6;
         /**
          * Max length between mismatch to do polishing
          * default 3
          */
         static int GAP_CUT; // = 3;


    // internal data member
    private:
        // nested class TODO: simplify in future versions
        struct BamAlignmentSupportData {
            /** 
             * cigarop, tags are stored here
             * Method to populate this field:
             *    BamReaderPrivate::LoadNextAlignment(BamAlignment& alignment) 
             */
            std::string AllCharData;
            uint32_t    BlockLength;  // not sure what this is
            /**
             * TODO: make sure This is the same as Cigar.size()?
             */
            uint32_t    NumCigarOperations; // should be calculated on the fly
            uint32_t    QueryNameLength;  // duplicate data discard in the future
            /**
             * TODO: this is rudundant with QueryBases.size()
             * consider remove this field in the future.
             */
            uint32_t    QuerySequenceLength; // is this duplicate of QueryLength?
            bool        HasCoreOnly;
            
            // constructor
            BamAlignmentSupportData()
                : BlockLength(0) , NumCigarOperations(0)
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
    //! \endinternal
};

// ---------------------------------------------------------
// BamAlignment tag access methods

// should use a data structure other than playing with fire here
// The optimization is useless.  TODO: replace with a proper 
// data structure
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

template<typename T>
inline bool BamAlignment::AddTag(const std::string& tag, const std::vector<T>& values) {
    // if char data not populated, do that first
    if (SupportData.HasCoreOnly) BuildCharData();
    // check for valid tag name length
    if (tag.size() != Constants::BAM_TAG_TAGSIZE) {
        cerr << __FILE__ << ":" << __LINE__ << ":WARN tag " + tag + " TAGSIZE wrong\n";
        return false;
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if (FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        cerr << __FILE__ << ":" << __LINE__ << ":WARN tag " + tag + " already exists\n";
        // TODO: set error string?
        return false;
    }
    // build new tag's base information, BAM_TAG_ARRAYBASE_SIZE=8
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    std::memcpy(newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE); // BAM_TAG_TAGSIZE=2
    newTagBase[2] = Constants::BAM_TAG_TYPE_ARRAY; // 'B'
    newTagBase[3] = TagTypeHelper<T>::TypeCode(); // uint32_t => 'I'

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    std::memcpy(newTagBase + 4, &numElements, sizeof(int32_t));
    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const size_t newTagDataLength = 
       tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE + numElements*sizeof(T);
    RaiiBuffer originalTagData(newTagDataLength);
    std::memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength+1); // '+1' for TagData's null-term
    // write newTagBase (removes old null term)
    std::strcat(originalTagData.Buffer + tagDataLength, (const char*)newTagBase);
    // add vector elements to tag
    int elementsBeginOffset = tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for (int i = 0 ; i < numElements; ++i) {
        const T& value = values.at(i);
        //cerr << "copying value: " << value << endl;
        memcpy(originalTagData.Buffer + elementsBeginOffset + i*sizeof(T), &value, sizeof(T));
    }
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer;
    TagData.assign(newTagData, newTagDataLength);
    return true;
}

template<typename T>
inline bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const T& value) 
{
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();
    // remove existing tag if present, then append tag with new value
    if ( HasTag(tag) )
        RemoveTag(tag);
    return AddTag(tag, type, value);
}

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
/**
 * Obtain the tag value tag and put into destination
 * @param tag name of the tag such as NM
 * @param destination value of which will be updated from this 
 *    object.
 */
template<typename T>
inline bool BamAlignment::GetTag(const std::string& tag, T& destination) const {
    // skip if alignment is core-only
    if (SupportData.HasCoreOnly) {
        // TODO: set error string?
        return false;
    }
    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data(); // string's low level representation
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // return failure if tag not found
    if (!FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        // TODO: set error string?
        return false;
    }
    char type = *(pTagData - 1);
    if (type == 'C') // C is uint32_t
    { 
      destination=(uint32_t)(*pTagData);
      return true;
    }
    size_t destinationLength = 0;
    switch ( type ) { // TODO: replace case with object design
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
            //SetErrorString("BamAlignment::GetTag",
            //              "cannot store variable length tag data into a numeric destination");
            cerr << __FILE__ << ":" << __LINE__ 
               << ":ERROR cannot store variable length tag data into a numeric destination\n";
            return false;
        // unrecognized tag type
        default:
            //const std::string message = std::string("invalid tag type: ") + type;
            //SetErrorString("BamAlignment::GetTag", message);
            cerr << __FILE__ << ":" << __LINE__ << ":ERROR invalid tag type: " << type << endl;
            return false;
    } // using a string specialization version for string data
    // store data in destination
    destination = 0;
    memcpy(&destination, pTagData, destinationLength);
    //cerr << "Just after memcpy destination value: " << destination << endl;
    return true;
}

/**
 * Specialization for string
 */
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
    if (SupportData.HasCoreOnly ) {
        cerr << __FILE__ << ":" << __LINE__ << ": bam has only core data failed to get tag: "
            << tag << "\n";
        // TODO: set error string?
        return false;
    }
    // skip if no tags present
    if (TagData.empty()) {
        cerr << __FILE__ << ":" << __LINE__ << ": empty tag data, failed to get tag value\n";
        return false;
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // return false if tag not found
    if (!FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        //cerr << __FILE__ << ":" << __LINE__ << ": cannot find tag: " << tag 
        //   << ", failed to get tag value\n";
        return false;
    }
    // check that tag is array type
    const char tagType = *(pTagData - 1);
    if (tagType != Constants::BAM_TAG_TYPE_ARRAY) {
        cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
           << ":ERROR cannot store a non-array tag in array destination\n";
        return false;
    }
    // fetch element type
    const char elementType = *pTagData;
    if (!TagTypeHelper<T>::CanConvertFrom(elementType) ) {
        cerr << __FILE__ << ":" << __LINE__ << ": cannot convert from "
            << elementType << endl;
        return false;
    }
    ++pTagData;

    // calculate length of each element in tag's array
    //int elementLength = 0;
    switch (elementType ) {
        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            //elementLength = sizeof(uint8_t);
            break;
        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            //elementLength = sizeof(uint16_t);
            break;
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
            //elementLength = sizeof(uint32_t);
            break;
        // var-length types not supported for numeric destination
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
               << ":ERROR invalid array data, variable-length elements are not allowed\n";
            return false;
        // unknown tag type
        default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
               << ":ERROR invalid array element type: " << elementType << endl;
            return false;
    }
    // not using elementLength
    // TODO: remove above code block
    //std::cerr << "BamTagData arrary element data width (Byte): " << elementLength << std::endl;
    // get number of elements
    int32_t numElements;
    memcpy(&numElements, pTagData, sizeof(int32_t));
    pTagData += 4;
    destination.clear();
    destination.reserve(numElements);

    // read in elements
    T value;
    for (int i = 0 ; i < numElements; ++i ) {
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
