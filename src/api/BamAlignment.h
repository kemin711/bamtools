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
//#include <mutex>

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
         * Convenient constructor for testing
         */
        BamAlignment(const std::string& qname, int32_t refid, int32_t refpos, uint32_t alnflag, 
              int32_t mrefid, int32_t mrefpos, const std::string& queryseq, const std::string& qstring)
           : Name(qname), Length(queryseq.size()),
             QueryBases(queryseq), AlignedBases(), Qualities(qstring),
             TagData(), RefID(refid), Position(refpos), Bin(0), MapQuality(0),
             AlignmentFlag(alnflag), CigarData(),
            MateRefID(mrefid), MatePosition(mrefpos), InsertSize(0),
            Filename()
        {}
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
    // // BAM alignment flags Designe problem, replaced with in class static constant
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
    //             This is determined by the order in the input file (READ1)
    // 128  0x80  the last segment in the template (second mate or read)
    //             Appears after READ1 in the input order.
    // 256  0x100 secondary alignment
    // 512  0x200 not passing filters, such as platform/vendor quality controls
    // 1024 0x400 PCR or optical duplicate
    // 2048 0x800 supplementary alignment
    // the following tests the flag field against different bits
    public:        
         /** 
          * This cannot be relied on.
          * @return true if this read is a PCR duplicate
         */
        bool IsDuplicate(void) const;
        bool IsFailedQC(void) const;          // returns true if this read failed quality control
        /**
         *  @return true if alignment is first mate on paired-end read
         *  First mate is also first read: Two words meaning the same thing.
         */
        bool IsFirstMate(void) const;         
        void setAlignmentFlag(uint32_t flag) {
           AlignmentFlag = flag;
        }
        /**
         * @return true if is the First Read in a paired end read.
         * Alias for IsFirstMate.
         * Use this version for carmel casing.
         * @see isFirstRead
         */
        bool isFirstMate(void) const { return (AlignmentFlag & READ_1) != 0; }
        bool isFirstRead(void) const { return (AlignmentFlag & READ_1) != 0; }
        /** 
         * @returns true if alignment is second mate on paired-end read
         */
        bool IsSecondMate(void) const;        
        /**
         * @return true if it is the second read (mate)
         */
        bool isSecondMate(void) const { return AlignmentFlag & READ_2; }
        // in C++ true is 1 false is 0
        bool isSecondRead(void) const { return AlignmentFlag & READ_2; }
        /**
         * @return 1 for first mate, 2 for second mate, 
         *     and 0 for unknown mate.
         */
        int getMate() const { if (isFirstMate()) return 1; 
           else if (isSecondMate()) return 2;
           else return 0;
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
           return !(AlignmentFlag & UNMAPPED);            
        }
        /** 
         * @return true if alignment is not mapped
         */
        bool isUnmapped() const {
           return AlignmentFlag & UNMAPPED;
        }
        /** 
         * Is a flag check
         * @return true if alignment's mate is mapped
         */
        bool IsMateMapped(void) const;        
        bool isMateMapped(void) const {
            return (AlignmentFlag & MATE_UNMAPPED) == 0;
        }
        /**
         * @return true if Mate is unmapped
         */
        bool isMateUnmapped() const {
           return (AlignmentFlag & MATE_UNMAPPED) != 0;
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
           return (AlignmentFlag & REVERSE_STRAND) != 0; 
        }
        bool isForwardStrand() const {
           return (AlignmentFlag & REVERSE_STRAND) == 0; 
        }
        /**
         * @return -1 reverse strand, +1 for forward strand, and
         *      0 for both strand or strand unknown.
         */
        int getStrand() const {
           if (isReverseStrand()) return -1;
           if (isForwardStrand()) return 1;
           return 0;
        }
        char getStrandChar() const {
           if (isReverseStrand()) return '-';
           if (isForwardStrand()) return '+';
           return '?';
        }
        /** 
         * @return true if alignment's mate mapped to reverse strand
         */
        bool IsMateReverseStrand(void) const; 
        bool isMateReverseStrand() const {
            return (AlignmentFlag & MATE_REVERSE_STRAND) != 0;
        }
        bool isMateForwardStrand() const {
            return (AlignmentFlag & MATE_REVERSE_STRAND) == 0;
        }
        /**
         * Get a value to represent the strand.
         * For simple alignment of single sequence it is either
         * +1 or -1.  For Merged sequence which is reprented as
         * having the XO tag this value could be [0, 1)
         * a fractional number between 0 and 1.
         * If the overlap is 100% of the read length, then it is zero.
         */
        double getFractionStrand() const;
        /** @returns true if alignment part of paired-end read
         */
        bool IsPaired(void) const; 
        bool isPaired() const { 
            return (AlignmentFlag & PAIRED) != 0;
        }
        /** 
         * @return true if reported position is primary alignment
         */
        bool IsPrimaryAlignment(void) const;  
        /**
         * Test SECONDARY flag is not set.
         */
         bool isPrimaryAlignment(void) const  {
            return !(AlignmentFlag & SECONDARY);
         }
        /**
         * @return true if is secondary alignment
         *    This indicates that the query was mapped
         *    multiple times in the genome. The choice of
         *    primary/secondary is usually arbitrary!
         */
        bool isSecondaryAlignment() const { 
           return AlignmentFlag & SECONDARY; }
        /**
         * Alignment is part of a chimera. The choice
         * of representative/supplementary is arbitrary.
         */
        bool isSupplementaryAlignment() const { 
           return AlignmentFlag & SUPPLEMENTARY; 
        }
        /**
         * Alias for isSupplementaryAlignment
         */
        bool isSupplementary() const { 
           return AlignmentFlag & SUPPLEMENTARY; 
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
            return (AlignmentFlag & PROPER_PAIR) != 0;
        }
        bool isNotProperPair() const {
            return (AlignmentFlag & PROPER_PAIR) == 0;
        }
        bool isImproperPair() const {
            return (AlignmentFlag & PROPER_PAIR) == 0;
        }

    // manipulate alignment flags
    public:        
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

        // convenient constants octal number
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
        template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values);

        // retrieves tag data
        /** 
         *  Retrieves the value associated with a BAM tag.
         *
         *  @param tag[in]          2-character tag name
         *  @param destination[out] retrieved value
         *  @return true if found.
         *
         *  Documents for TAGs
         *
         *  Type  | Regexp matching VALUE  |  Description
         *  ------------------------------------------------------
         *  A [!-~]                          Printable character
         *  i [-+]?[0-9]+                    Signed integer
         *  f [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? Single-precision 
         *                                   floating number
         *  Z [ !-~]*                        Printable string, including space
         *  H ([0-9A-F][0-9A-F])*            Byte array in the Hex format
         *  B [cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+ Integer or 
         *                                   numeric array
         *  ==============================================================
         *
         *  Standard tags
         *  ----------------------------------------------
         *  AM i The smallest template-independent mapping quality of segments in the rest
         *  AS i Alignment score generated by aligner
         *  BC Z Barcode sequence
         *  BQ Z Offset to base alignment quality (BAQ)
         *  CC Z Reference name of the next hit
         *  CM i Edit distance between the color sequence and the color reference (see also NM) 
         *  CO Z Free-text comments
         *  CP i Leftmost coordinate of the next hit
         *  CT Z Complete read annotation tag, used for consensus annotation dummy features.
         *  E2 Z The 2nd most likely base calls
         *  FI i The index of segment in the template
         *  FS Z Segment suffix
         *  FZ B,S Flow signal intensities
         *  GC ?  Reserved for backwards compatibility reasons
         *  GQ ?  Reserved for backwards compatibility reasons
         *  GS ?  Reserved for backwards compatibility reasons
         *  H0 i Number of perfect hits
         *  H1 i Number of 1-difference hits (see also NM)
         *  H2 i Number of 2-difference hits
         *  HI i Query hit index
         *  IH i Number of stored alignments in SAM that contains the query in the current record
         *  LB Z Library
         *  MC Z CIGAR string for mate/next segment
         *  MD Z String for mismatching positions
         *  MF ? Reserved for backwards compatibility reasons
         *  MQ i Mapping quality of the mate/next segment
         *  NH i Number of reported alignments that contains the query in the current record
         *  NM i Edit distance to the reference. This is use by BWA. MissMatch+Indel.
         *  OC Z Original CIGAR
         *  OP i Original mapping position
         *  OQ Z Original base quality
         *  PG Z Program
         *  PQ i Phred likelihood of the template
         *  PT Z Read annotations for parts of the padded read sequence
         *  PU Z Platform unit
         *  QT Z Barcode ( BC or RT) phred-scaled base qualities
         *  Q2 Z Phred quality of the mate/next segment sequence in the R2 tag
         *  R2 Z Sequence of the mate/next segment in the template
         *  RG Z Read group
         *  RT Z Barcode sequence (deprecated; use BC instead)
         *  SA Z Other canonical alignments in a chimeric alignment
         *  SM i Template-independent mapping quality
         *  SQ ?  Reserved for backwards compatibility reasons
         *  S2 ?  Reserved for backwards compatibility reasons
         *  TC i The number of segments in the template
         *  U2 Z Phred probility of the 2nd call being wrong conditional on the best being wrong
         *  UQ i Phred likelihood of the segment, conditional on the mapping being correct
         *  X?  ?  Reserved for end users
         *  Y?  ?  Reserved for end users
         *  Z?  ?  Reserved for end users
         * ==========================================
         *
         * Additional tags
         * -----------------------------------------------
         *  1.1 Additional Template and Mapping data
         *  ------------------------
         * AM:i:int The smallest template-independent mapping quality of segments in the rest.
         * AS:i:score Alignment score generated by aligner.
         * BQ:Z:qualities Offset to base alignment quality (BAQ), of the same length as the read 
         *    sequence. At the i-th read base, BAQi = Qi − (BQi−64) where Qi is the i-th base quality.
         * CC:Z:rname Reference name of the next hit; ‘=’ for the same chromosome.
         * CP:i:pos Leftmost coordinate of the next hit.
         * E2:Z:qualities The 2nd most likely base calls.  Same encoding and same length as QUAL.
         * FI:i:int The index of segment in the template.
         * FS:Z:str Segment suffix.
         * H0:i:count Number of perfect hits.
         * H1:i:count Number of 1-difference hits (see also NM).
         * H2:i:count Number of 2-difference hits.
         * HI:i:i Query hit index, indicating the alignment record is the i-th one stored in SAM.
         * IH:i:count Number of stored alignments in SAM that contains the query in the current record.
         * MC:Z:cigar CIGAR string for mate/next segment.
         * MD:Z:[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)* String for mismatching positions.
         *    The MD field aims to achieve SNP/indel calling without looking at
         *    the reference.  For example, a string ‘10A5^AC6’ means  from
         *    the  leftmost reference  base in the  alignment,  there  are 10
         *    matches followed by an A on the reference which is different from
         *    the aligned read base; the next 5 reference bases are matches
         *    followed by a 2bp deletion from the reference; the deleted
         *    sequence is AC; the last 6 bases are matches.  The MD field ought
         *    to match the CIGAR string.
         * MQ:i: Mapping quality of the mate/next segment.
         * NH:i: Number of reported alignments that contains the query in the 
         *    current record.
         * NM:i: Edit distance to the reference, including ambiguous bases but 
         *    excluding clipping.
         * PQ:i: Phred likelihood of the template, conditional on both the mapping 
         *    being correct.
         * Q2:Z: Phred quality of the mate/next segment sequence in the R2 tag.
         *    Same encoding as QUAL.
         * R2:Z: Sequence of the mate/next segment in the template.
         * SA:Z: ( rname , pos , strand , CIGAR , mapQ , NM ;)+ Other
         *     canonical alignments in a chimeric alignment, for- matted as a
         *     semicolon-delimited list.  Each element in the list represents a
         *     part of the chimeric align- ment.  Conventionally, at a
         *     supplementary line, the first element points to the primary line.
         * SM:i: Template-independent mapping quality.
         * TC:i: The number of segments in the template.
         * U2:Z: Phred probility of the 2nd call being wrong conditional on the best being wrong. 
         *    The same encoding as QUAL.
         * UQ:i: Phred likelihood of the segment, conditional on the mapping being correct.
         * --------------------------------
         * 1.2    Metadata
         * ------------------------------
         * RG:Z: readgroup The read group to which the read belongs.  If @RG headers are 
         *    present,  then readgroup must match the RG-ID field of one of
         *    the headers.
         * LB:Z: library The library from which the read has been sequenced.
         *    If @RG headers are present, then library must match the RG-LB
         *    field of one of the headers.  
         * PG:Z: Program.  Value matches the header PG-ID tag if @PG is present.
         * PU:Z:platformunit The platform unit in which the read was
         *    sequenced.  If @RG headers are present, then platformunit must
         *    match the RG-PU field of one of the headers.
         * CO:Z:text Free-text comments.
         *  ----------------------------
         *  1.3    Barcodes
         * ----------------------------
         * BC:Z:sequence Barcode sequence, with any quality scores stored in the
         *    QT tag.
         * QT:Z:qualities Phred quality of the barcode sequence in the BC (or RT) tag.  
         *    Same encoding as QUAL .
         * RT:Z:sequence Deprecated alternative to BC tag originally used at Sanger.
         * ----------------------------
         * 1.4    Original data
         * ----------------------------
         * OC:Z:cigar Original CIGAR, usually before realignment.
         * OP:i:pos Original mapping position, usually before realignment.
         * OQ:Z:qualities Original base quality, usually before recalibration. 
         *    Same encoding as QUAL
         * ----------------------------
         * 1.5    Annotation and Padding
         * ----------------------------
         * CT:Z:strand;type(;key(=value))* Complete read annotation tag, used for consensus 
         *    annotation dummy features.  The CT tag is intended primarily for
         *    annotation dummy reads, and consists of a strand , type and zero
         *    or more key=value pairs, each separated with semicolons.  The
         *    strand field has four values as in GFF3, and supplements FLAG
         *    bit 0x10 to allow unstranded (‘.’),  and stranded but unknown
         *    strand (‘?’) annotation.  For these and annotation on the
         *    forward strand (strand set to ‘+’), do not set FLAG bit 0x10.
         *    For annotation on the reverse strand, set the strand to ‘-’ and
         *    set FLAG bit 0x10.  The type and  any keys and  their optional
         *    values are  all  percent  encoded  according  to RFC3986  to
         *    escape meta-characters ‘=’,‘%’,‘;’,‘|’ or non-printable
         *    characters not matched by the isprint() macro (with the C
         *    locale).  For example a percent sign becomes ‘%2C’.
         * PT:Z:start;end;strand;type(;key(=value))*(\|start;end;strand;type(;key(=value))*)*
         *   Read annotations for parts of the padded read sequence.  The PT
         *   tag value has the format of a series of tags separated by ‘|’,
         *   each annotating a sub-region of the read.  Each tag consists
         *   of start, end, strand , type and zero or more key = value pairs,
         *   each separated with semicolons.  Start and end are 1-based
         *   positions between one and the sum of the M/I/D/P/S/=/X CIGAR
         *   operators, i.e.  SEQ length plus any pads.  Note any editing of
         *   the CIGAR string may require updating the ‘PT’ tag coordinates,
         *   or even invalidate them.  As in GFF3, strand is one of ‘+’ for
         *   forward strand tags, ‘-’ for reverse strand, ‘.’  for unstranded
         *   or ‘?’  for stranded but unknown strand.  The type and any keys
         *   and their optional values are all percent encoded as in the CT
         *   tag.  1.6 Technology-specific data FZ:B,S: intensities Flow
         *   signal intensities  on  the original strand of  the  read, stored
         *   as (uint16 t) round(value*100.0).
         *   -------------------------
         *       1.6.1  Color space
         *   -------------------------
         * CM:i:distance Edit distance between the color sequence and the
         *    color reference (see also NM).
         * CS:Z:sequence Color read sequence on the original strand of the
         *    read.  The primer base must be included.
         * CQ:Z:qualities Color  read  quality  on  the  original  strand
         *    of  the  read.   Same  encoding  as QUAL ;  same length as CS .
         * ----------------------------------
         * 2    Locally-defined tags
         * ----------------------------------
         * You can freely add new tags.  Note that tags starting with ‘X’, ‘Y’,
         * or ‘Z’ and tags containing lowercase letters in either position are
         * reserved for local use and will not be formally defined in any
         * future version of this specification.  If a new tag may be of
         * general interest, it may be useful to have it added to this 
         * specification.  Additions can be proposed by opening a new issue 
         * at https://github.com/samtools/hts-specs/issues
         * and/or by sending email to samtools-devel@lists.sourceforge.net.
        */
        template<typename T> bool GetTag(const std::string& tag, T& destination) const;
        template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const;
        // retrieves all current tag names
        std::vector<std::string> GetTagNames(void) const;
        // retrieves the SAM/BAM type-code for requested tag name
        bool GetTagType(const std::string& tag, char& type) const;
        // retrieves the SAM/BAM type-code for the data elements in an array tag
        bool GetArrayTagType(const std::string& tag, char& type) const;
        /** 
         * returns true if alignment has a record for this tag name
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
        /**
         * The distance coverted by the alignment on the reference.
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
        const std::string& getName() const { return Name; }
        /**
         * getter method for the length of the query (read)
         */
        int32_t getQueryLength() const { return Length; }
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
         * returned.
         */
        const std::string& getQuerySequence() const { return QueryBases; }
        /**
         * 'aligned' sequence (QueryBases plus deletion, padding, clipping chars)
         */
        const std::string& getAlignedQueryBases() const { 
           //if (AlignedBases.empty()) {
           //   cerr << "empty aligned bases!\n";
           //}
           return AlignedBases; 
        }
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
        std::string getQuality() const { return Qualities; }
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
         * RefVector <= vector<RefData>, RefData { RefName, RefLength }
         * unnecessary typedef, more confusing than help.
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
         * @return mape position
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
         */
        const std::vector<CigarOp>& getCigar() const { return CigarData; } 
        /**
         * Provide a more user-friendly interface
         * for working with other applications.
         */
        vector<pair<char,int> > getCigarOperation() const;
        string getCigarString() const;
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
         * Alignment has soft clip on either start
         * or end. Not checking middle, which may not make sense.
         * @return true of has soft clip otherwise, including no cigar,
         *    then return false.
         */
        bool hasSoftclip() const;
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
         * @return the query sequence for the first soft clip.
         *    If there is no soft clip then an empty string is returned.
         */
        string getFirstSoftclip() const;
        /**
         * @return the length of the first softclip
         */
        int getFirstSoftclipLength() const;
        /**
         * @return the last soft clip in query sequence.
         */
        string getLastSoftclip() const;
        int getLastSoftclipLength() const;
        /**
         * @return sum of softclip length if both are present.
         */
        int getSoftclipLength() const;
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

        /// setter methods
        /**
         * change the name of the query
         */
        void setQueryName(const std::string &qname) { Name = qname; }
        void setQuerySequenceLength(int32_t qlen) {  Length = qlen; }
        void setQueryLength(int32_t qlen) {  Length = qlen; }
        /**
         * Once query bases is changed the length will also
         * change. There is no need to call setQueryLength()
         * after calling this method.
         */
        void setQueryBases(const std::string &qseq) { 
           QueryBases = qseq; 
           setQueryLength(QueryBases.size());
        }
        /** set quality from string data */
        void setQuality(const std::string &qual) { Qualities = qual; }
        /** integer version
         * @param qual is a vector of Phred scores from 0-63
         */
        void setQuality(const std::vector<int> &qual);
        void setRefID(int32_t refid) { RefID = refid; }
        void setReferenceId(int32_t refid) { RefID = refid; }
        void setPosition(int32_t alnstart) { Position = alnstart; }
        /** alias for setPosition */
        void setStart(int32_t alnstart) { Position = alnstart; }
        void setBin(uint16_t indexbin) { Bin = indexbin; }
        void setMapQuality(uint16_t mqual) { MapQuality = mqual; }
        /**
         * update CigarData with a new value.
         */
        void setCigarData(const std::vector<CigarOp> &cd) { CigarData = cd; } 
        /** 
         * set CigarData from a vector a pair {char, int}
         *  for easy communication with external world.
         *  @param cd CigarData in more universal std data type
         * */
        void setCigarOperation(const std::vector<pair<char,int> > &cd); 
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

        // mutation functions
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
        /** even position number, odd position letter last 3 digits
         * for bases 00 A, 01 C, 10 G, 11 T, 100 N, first bit del
         * MD:Z:20^A127
         * MD:Z:108^TTCTAAGGCCAGCTCCTGCACC39
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
         * Only change the bases near the end < 7 nt
         * if the Query base is different from the reference base
         * on both ends. MD tag will also needs to be updated.
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
         * @return -1 if not found NM tag
         */
         int getNMValue() const;
         static void setPolishMax(int len) {
            TRIMLEN_MAX = len;
         }
         static void setPolishGap(int gap) {
            GAP_CUT = gap;
         }

    private:
      void advanceIndex(int &i, int &j, int &b, unsigned int &cigarIdx, unsigned int &ci, char &cigarState) const;

    // public data fields, these fileds should all become private in the future
    public:
        /** read or query name 
         * Use getName() to read this one
         * */
        std::string Name;    
        /** 
         * length of query sequence
         * Design flaw: This field is redundant with
         *   SupporData::QuerySequenceLength. Both fields
         *   needs to be updated!
         *   This filed is also redundant with the QueryBases.size()
         *
         *  TODO: consider redesign 
         */
        int32_t     Length;             
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
        /** ID number for reference sequence
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
        // alignment should not store its file name
        // information repetation, remove in future version
        // TODO: remove in next release
        // this is used in multiple file input operations
        std::string Filename;           // name of BAM file which this alignment comes from
        //static mutex gmtx;

    //! \internal
    // internal utility methods
    private:
        /** 
         *  Searches for requested tag in BAM tag data.
         *  @param  tag            requested 2-character tag name
         *  @param  pTagData       pointer to current position in BamAlignment::TagData
         *  @param  tagDataLength  length of BamAlignment::TagData
         *  @param  numBytesParsed number of bytes parsed so far
         *  @return true if found
         *  If tag is found, pTagData will point to the byte where the tag data begins.
         *        numBytesParsed will correspond to the position in the full TagData string.
        */
        bool FindTag(const std::string& tag, char*& pTagData, 
              const unsigned int& tagDataLength, unsigned int& numBytesParsed) const;
        bool IsValidSize(const std::string& tag, const std::string& type) const;
        void SetErrorString(const std::string& where, const std::string& what) const;
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


    // internal data
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
    const size_t newTagDataLength = 
       tagDataLength + Constants::BAM_TAG_ARRAYBASE_SIZE + numElements*sizeof(T);
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
    char* pTagData = (char*)TagData.data(); // string's low leve representation
    //cout << __FILE__ << ":" << __LINE__ << ": pTagData: "
    //   << pTagData << endl;
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // return failure if tag not found
    if (!FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }
    //cout << __FILE__ << ":" << __LINE__ << ": after FindTag operation. pTagData-1: "
    //   << pTagData-1 << endl;
    // fetch data type
    //const char type = *(pTagData - 1);
    char type = *(pTagData - 1);
    if (type == 'C' && (tag == "NM" || tag == "AS" || tag == "XS")) { // patch a bug
       //type = Constants::BAM_TAG_TYPE_INT32;
      //cout << __FILE__ << ":" << __LINE__ << ": pTagData: "
      //   << pTagData << endl;
      //cerr << "NM int val at pTagData: " << (int)(*pTagData) << " at 2 bytes later: "
      //   << (int)(*(pTagData+2)) << endl;
       // this is a short term patc for the bug
      destination=(int)(*pTagData);
      return true;
    }
    else if ( !TagTypeHelper<T>::CanConvertFrom(type) ) {
        // TODO: set error string ?
       cerr << __FILE__ << ":" << __LINE__ << ":"
          << "Failed to convert to " << typeid(destination).name() 
          << " from " << type << endl;
        return false;
    }
    // determine data length
    //int destinationLength = 0;
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
    //cerr << "Just after memcpy destination value: " << destination << endl;
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
