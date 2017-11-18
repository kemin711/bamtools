// ***************************************************************************
// BamWriter_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#ifndef BAMWRITER_P_H
#define BAMWRITER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#include "api/BamAux.h"
#include "api/internal/io/BgzfStream_p.h"
#include <string>
#include <vector>

namespace BamTools {

class BamAlignment;

namespace Internal {

class BamWriterPrivate {

    // ctor & dtor
    public:
        BamWriterPrivate(void);
        ~BamWriterPrivate(void);

    // interface methods
    public:
        void Close(void);
        std::string GetErrorString(void) const;
        bool IsOpen(void) const;
        bool Open(const std::string& filename,
                  const std::string& samHeaderText,
                  const BamTools::RefVector& referenceSequences);
        /**
         * Used by BamWriter
         * Implementation of the BamWriter::SaveAlignment() interface
         * saves the alignment to the alignment archive
         *
         * @param al BamAlignment object
         *
         * The Open method opens up an output stream
         *
         * For a full alignment object, this method
         * calls WriteAlignment(al)
         *
         * @see WriteAlignment
         */
        bool SaveAlignment(const BamAlignment& al);
        void SetWriteCompressed(bool ok);

    // 'internal' methods
    public:
        /**
         * Don't know what does this do?
         */
        uint32_t CalculateMinimumBin(const int begin, int end) const;
        /**
         * Convert the vector version of the CigarOperaton into
         * a string representation.
         */
        void CreatePackedCigar(const std::vector<BamTools::CigarOp>& cigarOperations, std::string& packedCigar);
        /**
         * Not sure what this does, save space?
         */
        void EncodeQuerySequence(const std::string& query, std::string& encodedQuery);
        /**
         * Method used by SaveAlignment()
         * This is a lengthy function.
         * Use GgzfStream::write() method.
         * I am trying to figure out how it writing the insert size
         * as zero even when I set it to a non-zero!
         */
        void WriteAlignment(const BamAlignment& al);
        void WriteCoreAlignment(const BamAlignment& al);
        void WriteMagicNumber(void);
        void WriteReferences(const BamTools::RefVector& referenceSequences);
        void WriteSamHeaderText(const std::string& samHeaderText);

    // data members
    private:
        /**
         * Need to read the details of this class
         */
        BgzfStream m_stream;
        /**
         * Machine type
         */
        bool m_isBigEndian;
        std::string m_errorString;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMWRITER_P_H
