// ***************************************************************************
// BamWriter.h (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#ifndef BAMWRITER_H
#define BAMWRITER_H

#include "api/api_global.h"
#include "api/BamAux.h"
#include "api/internal/bam/BamWriter_p.h"
#include <string>
#include "api/SamHeader.h"

using namespace BamTools::Internal;

namespace BamTools {

class BamAlignment;
struct SamHeader;

//! \cond
namespace Internal {
    class BamWriterPrivate;
} // namespace Internal
//! \endcond

class API_EXPORT BamWriter {
    public:
        enum CompressionMode { Compressed = 0, Uncompressed };
         /**
          * Default constructor.
         */
         BamWriter(void) : d(new BamWriterPrivate) { }
        /**
         *  destructor. 
         *  Close d. Deallocate memory and set d to nullptr
         */
        ~BamWriter(void);

    // public interface
    public:
        //  closes the current BAM file
        void Close(void);
        // returns a human-readable description of the last error that occurred
        std::string GetErrorString(void) const;
        // returns true if BAM file is open for writing
        bool IsOpen(void) const;
        /**
         * Open the Bam file from a string representation of the sam header
         * for writing.
         *  Will overwrite the BAM file if it already exists.
         *  \param[in] filename           name of output BAM file
         *  \param[in] samHeaderText      header data, as SAM-formatted string
         *  \param[in] referenceSequences list of reference entries
         *  @return \c true if opened successfully
         *  \sa Close(), IsOpen(), BamReader::GetHeaderText(), BamReader::GetReferenceData()
         */
        bool Open(const std::string& filename, 
                  const std::string& samHeaderText,
                  const RefVector& referenceSequences)
        {
            return d->Open(filename, samHeaderText, referenceSequences);
        }
        /** 
         * Typical usage:
         *  Open(outBamFile, br.getHeader(), br.getReferenceData());
         *
         *  Opens a BAM file for writing from a SamHeaer object.
         *  Will overwrite the BAM file if it already exists.
         *  \param[in] filename           name of output BAM file
         *  \param[in] samHeader          header data, wrapped in SamHeader object
         *  \param[in] referenceSequences list of reference entries
         *  \return \c true if opened successfully
         *  \sa Close(), IsOpen(), BamReader::GetHeader(), BamReader::GetReferenceData()
         */
        bool Open(const std::string& filename,
                  const SamHeader& samHeader,
                  const RefVector& referenceSequences)
        {
            return d->Open(filename, samHeader.ToString(), referenceSequences);
        }
        /** 
         * Saves an alignment to the BAM file.
         * 
         * @param[in] alignment BamAlignment record to save
         * @see  BamReader::GetNextAlignment(), BamReader::GetNextAlignmentCore()
         * @return true of success false if fail.
         *
         * Note: This method delegates to BamWriterPrivate::SaveAlignment
         */
        bool SaveAlignment(const BamAlignment& alignment) {
            return d->SaveAlignment(alignment);
        }
        bool saveAlignment(const BamAlignment& alignment) {
            return d->SaveAlignment(alignment);
        }
        // sets the output compression mode
        void SetCompressionMode(const BamWriter::CompressionMode& compressionMode);

    private:
        Internal::BamWriterPrivate* d;
};

} // namespace BamTools

#endif // BAMWRITER_H
