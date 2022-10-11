// ***************************************************************************
// BamReader_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 18 November 2012 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#ifndef BAMREADER_P_H
#define BAMREADER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/BamAlignment.h"
#include "api/BamIndex.h"
#include "api/BamReader.h"
#include "api/SamHeader.h"
#include "api/internal/bam/BamHeader_p.h"
#include "api/internal/bam/BamRandomAccessController_p.h"
#include "api/internal/io/BgzfStream_p.h"
#include <string>

namespace BamTools {
class BamReader;
namespace Internal {

using namespace BamTools;

/**
 * Should have a hiarchy of different readers
 * Base: Sequential readers
 * Random readers for different situations.
 */
class BamReaderPrivate {
    // ctor & dtor
    public:
       /**
        * Hold a pointer to the BamReader parent
        * So that BamReader could have multiple 
        * objects such as random controller in addtion to this one?
        */
        BamReaderPrivate(BamTools::BamReader* parent);
        ~BamReaderPrivate(void);

    // BamReader interface
    public:

        // file operations
        bool Close(void);
        const std::string Filename(void) const;
        bool IsOpen(void) const;
        bool Open(const std::string& filename);
        bool Rewind(void);
        bool SetRegion(const BamRegion& region);

        /**
         * Take alignment data from file stream and 
         * fill up the Alignment object.
         * @return true if obtained the next valid alignment
         *    false otherwise.
         */
        bool GetNextAlignment(BamAlignment& alignment);
        /**
         * Collect only the core part of the alignment.
         */
        bool GetNextAlignmentCore(BamAlignment& alignment);

        // access auxiliary data
        std::string GetHeaderText(void) const;
        const SamHeader& GetConstSamHeader(void) const;
        SamHeader GetSamHeader(void) const;
        int GetReferenceCount(void) const;
        const RefVector& GetReferenceData(void) const;
        int GetReferenceID(const std::string& refName) const;

        // index operations
        bool CreateIndex(const BamIndex::IndexType& type);
        bool HasIndex(void) const;
        bool LocateIndex(const BamIndex::IndexType& preferredType);
        bool OpenIndex(const std::string& indexFilename);
        void SetIndex(BamIndex* index);

    // internal methods, but available as a BamReaderPrivate 'interface'
    //
    // these methods should only be used by BamTools::Internal classes
    // (currently only used by the BamIndex subclasses)
    public:
        // retrieves header text from BAM file
        void LoadHeaderData(void);
        /**
         * retrieves BAM alignment under file pointer
         * (does no overlap checking or character data parsing)
         * Parse the Cigar data into CigarData field.
         *
         * @return true if success false if failure.
         */
        bool LoadNextAlignment(BamAlignment& alignment);
        // builds reference data structure from BAM file
        bool LoadReferenceData(void);
        // seek reader to file position
        bool Seek(const int64_t& position);
        // return reader's file position
        int64_t Tell(void) const;

    // data members
    public:
        // general BAM file data
        int64_t     m_alignmentsBeginOffset;
        std::string m_filename;
        /**
         * vector of [refname, reflen] index by refid
         */
        RefVector   m_references;
        /** 
         * system data
         * TODO: Should auto detect.
         */
        bool m_isBigEndian;
        /**
         * parent BamReader
         * TODO: poor design should remove.
         */
        BamReader* m_parent;

        // BamReaderPrivate components
        /**
         * Holds the BamFile header
         */
        BamHeader m_header;
        /**
         * Random access controller
         */
        BamRandomAccessController m_randomAccessController;
        /**
         * Bgzfile stream
         */
        BgzfStream m_stream;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMREADER_P_H
