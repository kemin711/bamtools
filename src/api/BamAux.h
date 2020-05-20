// ***************************************************************************
// BamAux.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides data structures & utility methods that are used throughout the API.
// ***************************************************************************
// This is a header-only file, no implementation

#ifndef BAMAUX_H
#define BAMAUX_H

#include "api/api_global.h"
#include <cstring>
#include <fstream> 
#include <iostream>
#include <string>
#include <vector>

using namespace std;

/*! \file BamAux.h

    Provides data structures & utility methods that are used throughout the API.
*/

/*! \namespace BamTools
    \brief Contains all BamTools classes & methods.

    The BamTools API contained in this namespace contains classes and methods
    for reading, writing, and manipulating BAM alignment files.
*/
namespace BamTools {

// ----------------------------------------------------------------
// CigarOp

/*! \struct BamTools::CigarOp
    \brief Represents a CIGAR alignment operation.

    \sa \samSpecURL for more details on using CIGAR operations.
*/
struct API_EXPORT CigarOp {
    char     Type;   //!< CIGAR operation type (MIDNSHPX=)
    uint32_t Length; //!< CIGAR operation length (number of bases)
    
    /** 
     * default constructor 
     */
    CigarOp() : Type('\0'), Length(0) { }
    /**
     * Constructor from known value
     */
    CigarOp(const char type, const uint32_t& length)
        : Type(type) , Length(length) { }
    friend ostream& operator<<(ostream& ous, const CigarOp& co) {
       //ous << co.Type << co.Length; return ous; 
       ous << co.Length << co.Type; return ous; 
    }
    /**
     * Convert this to a more universal type
     * If this class does not provide any special operation,
     * then there is no use to use a special type. Just 
     * more noise.
     */
    pair<char, int> topair() const { return pair<char,int>(Type, Length); }
    void fromPair(const pair<char,int> &p) { Type=p.first; Length = p.second; }
    /**
     * @return the length as 32-bits integer
     */
    uint32_t getLength() const { return Length; }
    /**
     * @return the cigar operation one of MIDNSHPX=
     */
    char getType() const { return Type; }
    /**
     * change the length for the cigar segment to l.
     * @param l new length.
     */
    void setLength(uint32_t l) { Length = l; }
    void shrink(uint32_t l) { 
       if (Length <= l) throw logic_error("Length become zoer or negative");
       Length -= l; 
    }
    void expand(uint32_t l) { Length += l; }
    void setType(char t) { Type = t; }
};

// ----------------------------------------------------------------
// RefData

/*! \struct BamTools::RefData
    \brief Represents a reference sequence entry
*/
struct API_EXPORT RefData {
    std::string RefName;    //!< name of reference sequence
    int32_t     RefLength;  //!< length of reference sequence
    
    //! constructor
    RefData(const std::string& name = "",
            const int32_t& length = 0)
        : RefName(name)
        , RefLength(length)
    { }
    pair<string, int32_t> asPair() const {
       return make_pair(RefName, RefLength);
    }
    const string& getRefname() const {
       return RefName;
    }
    int32_t getReflength() const {
       return RefLength;
    }
};

/**
 * Convenience typedef for vector of RefData entries.
 * Alias for vector of RefData.
 */
typedef std::vector<RefData> RefVector;

// ----------------------------------------------------------------
// BamRegion

/** 
 *  Represents a sequential genomic region
 *
 *  Allowed to span multiple (sequential) references.
 *
 *  Warning: BamRegion now represents a zero-based, HALF-OPEN interval.
 *  In previous versions of BamTools (0.x & 1.x) all intervals were treated
 *  as zero-based, CLOSED.
 *
 *  Region can span multiple references if the references are
 *  sorted in some order.
*/
struct API_EXPORT BamRegion {
    int LeftRefID;      //!< reference ID for region's left boundary
    int LeftPosition;   //!< position for region's left boundary
    int RightRefID;     //!< reference ID for region's right boundary
    int RightPosition;  //!< position for region's right boundary
    
    //! constructor
    BamRegion(const int& leftID   = -1, const int& leftPos  = -1,
              const int& rightID  = -1, const int& rightPos = -1)
        : LeftRefID(leftID), LeftPosition(leftPos)
        , RightRefID(rightID), RightPosition(rightPos)
    { }
    /**
     * Constructor for only one reference
     * @param refid reference id, chr1 is 1, chrm is 0, chrX is 23
     * @param leftp left position
     * @param rightp right position
     */
    BamRegion(int refid, int leftp, int rightp) 
       : LeftRefID(refid), LeftPosition(leftp), RightRefID(refid), RightPosition(rightp)
    { }
    
    /** 
      * copy constructor
    */
    BamRegion(const BamRegion& other)
        : LeftRefID(other.LeftRefID), LeftPosition(other.LeftPosition)
        , RightRefID(other.RightRefID), RightPosition(other.RightPosition)
    { }
    BamRegion& operator=(const BamRegion& o) {
	if (this != &o) {
	LeftRefID=o.LeftRefID;
	LeftPosition=o.LeftPosition;
	RightRefID=o.RightRefID;
	RightPosition=o.RightPosition;
	}
	return *this;
   }
    
    //! Clears region boundaries
    void clear(void) {
        LeftRefID  = -1; LeftPosition  = -1;
        RightRefID = -1; RightPosition = -1;
    }

    //! Returns true if region has a left boundary
    bool isLeftBoundSpecified(void) const {
        return ( LeftRefID >= 0 && LeftPosition >= 0 );
    }

    //! Returns true if region boundaries are not defined
    bool isNull(void) const {
        return ( !isLeftBoundSpecified() && !isRightBoundSpecified() );
    }

    //! Returns true if region has a right boundary
    bool isRightBoundSpecified(void) const {
        return ( RightRefID >= 0 && RightPosition >= 1 );
    }
    bool isSingleReference() const {
       return LeftRefID == RightRefID && LeftRefID != -1;
    }
};

struct CustomHeaderTag {
  std::string TagName;
  std::string TagValue;
};


// ----------------------------------------------------------------
// General utility methods

/*! \fn bool FileExists(const std::string& filename)
    \brief returns true if the file exists
*/
API_EXPORT inline bool FileExists(const std::string& filename) {
    std::ifstream f(filename.c_str(), std::ifstream::in);
    return !f.fail();
}

/*! \fn void SwapEndian_16(int16_t& x)
    \brief swaps endianness of signed 16-bit integer, in place
*/
API_EXPORT inline void SwapEndian_16(int16_t& x) {
    x = ((x >> 8) | (x << 8));
}

/*! \fn void SwapEndian_16(uint16_t& x)
    \brief swaps endianness of unsigned 16-bit integer, in place
*/
API_EXPORT inline void SwapEndian_16(uint16_t& x) {
    x = ((x >> 8) | (x << 8));
}

/*! \fn void SwapEndian_32(int32_t& x)
    \brief swaps endianness of signed 32-bit integer, in place
*/
API_EXPORT inline void SwapEndian_32(int32_t& x) {
    x = ( (x >> 24) | 
         ((x << 8) & 0x00FF0000) | 
         ((x >> 8) & 0x0000FF00) | 
          (x << 24)
        );
}

/*! \fn void SwapEndian_32(uint32_t& x)
    \brief swaps endianness of unsigned 32-bit integer, in place
*/
API_EXPORT inline void SwapEndian_32(uint32_t& x) {
    x = ( (x >> 24) | 
         ((x << 8) & 0x00FF0000) | 
         ((x >> 8) & 0x0000FF00) | 
          (x << 24)
        );
}

/*! \fn void SwapEndian_64(int64_t& x)
    \brief swaps endianness of signed 64-bit integer, in place
*/
API_EXPORT inline void SwapEndian_64(int64_t& x) {
    x = ( (x >> 56) | 
         ((x << 40) & 0x00FF000000000000ll) |
         ((x << 24) & 0x0000FF0000000000ll) |
         ((x << 8)  & 0x000000FF00000000ll) |
         ((x >> 8)  & 0x00000000FF000000ll) |
         ((x >> 24) & 0x0000000000FF0000ll) |
         ((x >> 40) & 0x000000000000FF00ll) |
          (x << 56)
        );
}

/*! \fn void SwapEndian_64(uint64_t& x)
    \brief swaps endianness of unsigned 64-bit integer, in place
*/
API_EXPORT inline void SwapEndian_64(uint64_t& x) {
    x = ( (x >> 56) | 
         ((x << 40) & 0x00FF000000000000ll) |
         ((x << 24) & 0x0000FF0000000000ll) |
         ((x << 8)  & 0x000000FF00000000ll) |
         ((x >> 8)  & 0x00000000FF000000ll) |
         ((x >> 24) & 0x0000000000FF0000ll) |
         ((x >> 40) & 0x000000000000FF00ll) |
          (x << 56)
        );
}

/*! \fn void SwapEndian_16p(char* data)
    \brief swaps endianness of the next 2 bytes in a buffer, in place
*/
API_EXPORT inline void SwapEndian_16p(char* data) {
    uint16_t& value = (uint16_t&)*data; 
    SwapEndian_16(value);
}

/*! \fn void SwapEndian_32p(char* data)
    \brief swaps endianness of the next 4 bytes in a buffer, in place
*/
API_EXPORT inline void SwapEndian_32p(char* data) {
    uint32_t& value = (uint32_t&)*data; 
    SwapEndian_32(value);
}

/*! \fn void SwapEndian_64p(char* data)
    \brief swaps endianness of the next 8 bytes in a buffer, in place
*/
API_EXPORT inline void SwapEndian_64p(char* data) {
    uint64_t& value = (uint64_t&)*data; 
    SwapEndian_64(value);
}

/*! \fn bool SystemIsBigEndian(void)
    \brief checks host architecture's byte order
    \return \c true if system uses big-endian ordering
*/
API_EXPORT inline bool SystemIsBigEndian(void) {
   const uint16_t one = 0x0001;
   return ((*(char*) &one) == 0 );
}

/*! \fn void PackUnsignedInt(char* buffer, unsigned int value)
    \brief stores unsigned integer value in a byte buffer

    \param[out] buffer destination buffer
    \param[in]  value  value to 'pack' in buffer
*/
API_EXPORT inline void PackUnsignedInt(char* buffer, unsigned int value) {
    buffer[0] = (char)value;
    buffer[1] = (char)(value >> 8);
    buffer[2] = (char)(value >> 16);
    buffer[3] = (char)(value >> 24);
}

/*! \fn void PackUnsignedShort(char* buffer, unsigned short value)
    \brief stores unsigned short integer value in a byte buffer

    \param[out] buffer destination buffer
    \param[in]  value  value to 'pack' in buffer
*/
API_EXPORT inline void PackUnsignedShort(char* buffer, unsigned short value) {
    buffer[0] = (char)value;
    buffer[1] = (char)(value >> 8);
}

/*! \fn double UnpackDouble(const char* buffer)
    \brief reads a double value from byte buffer

    \param[in] buffer source byte buffer
    \return the (double) value read from the buffer
*/
API_EXPORT inline double UnpackDouble(const char* buffer) {
    union { double value; unsigned char valueBuffer[sizeof(double)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    un.valueBuffer[4] = buffer[4];
    un.valueBuffer[5] = buffer[5];
    un.valueBuffer[6] = buffer[6];
    un.valueBuffer[7] = buffer[7];
    return un.value;
}

/*! \fn double UnpackDouble(char* buffer)
    \brief reads a double value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (double) value read from the buffer
*/
API_EXPORT inline double UnpackDouble(char* buffer) {
    return UnpackDouble( (const char*)buffer );
}

/*! \fn double UnpackFloat(const char* buffer)
    \brief reads a float value from byte buffer

    \param[in] buffer source byte buffer
    \return the (float) value read from the buffer
*/
API_EXPORT inline float UnpackFloat(const char* buffer) {
    union { float value; unsigned char valueBuffer[sizeof(float)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn double UnpackFloat(char* buffer)
    \brief reads a float value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (float) value read from the buffer
*/
API_EXPORT inline float UnpackFloat(char* buffer) {
    return UnpackFloat( (const char*)buffer );
}

/*! \fn signed int UnpackSignedInt(const char* buffer)
    \brief reads a signed integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (signed int) value read from the buffer
*/
API_EXPORT inline signed int UnpackSignedInt(const char* buffer) {
    union { signed int value; unsigned char valueBuffer[sizeof(signed int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn signed int UnpackSignedInt(char* buffer)
    \brief reads a signed integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (signed int) value read from the buffer
*/
API_EXPORT inline signed int UnpackSignedInt(char* buffer) {
    return UnpackSignedInt( (const char*) buffer );
}

/*! \fn signed short UnpackSignedShort(const char* buffer)
    \brief reads a signed short integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (signed short) value read from the buffer
*/
API_EXPORT inline signed short UnpackSignedShort(const char* buffer) {
    union { signed short value; unsigned char valueBuffer[sizeof(signed short)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    return un.value;
}

/*! \fn signed short UnpackSignedShort(char* buffer)
    \brief reads a signed short integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (signed short) value read from the buffer
*/
API_EXPORT inline signed short UnpackSignedShort(char* buffer) {
    return UnpackSignedShort( (const char*)buffer );
}

/*! \fn unsigned int UnpackUnsignedInt(const char* buffer)
    \brief reads an unsigned integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (unsigned int) value read from the buffer
*/
API_EXPORT inline unsigned int UnpackUnsignedInt(const char* buffer) {
    union { unsigned int value; unsigned char valueBuffer[sizeof(unsigned int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn unsigned int UnpackUnsignedInt(char* buffer)
    \brief reads an unsigned integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (unsigned int) value read from the buffer
*/
API_EXPORT inline unsigned int UnpackUnsignedInt(char* buffer) {
    return UnpackUnsignedInt( (const char*)buffer );
}

/*! \fn unsigned short UnpackUnsignedShort(const char* buffer)
    \brief reads an unsigned short integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (unsigned short) value read from the buffer
*/
API_EXPORT inline unsigned short UnpackUnsignedShort(const char* buffer) {
    union { unsigned short value; unsigned char valueBuffer[sizeof(unsigned short)]; } un;
    un.value = 0;
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    un.valueBuffer[0] = buffer[1];
    un.valueBuffer[1] = buffer[0];
#else
    #error "Unsupported hardware"
#endif
    return un.value;
}

/*! \fn unsigned short UnpackUnsignedShort(char* buffer)
    \brief reads an unsigned short integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (unsigned short) value read from the buffer
*/
API_EXPORT inline unsigned short UnpackUnsignedShort(char* buffer) {
    return UnpackUnsignedShort( (const char*)buffer );
}

// ----------------------------------------------------------------
// 'internal' helper structs

/**
 * struct RaiiBuffer
 *   internal
 */
struct RaiiBuffer {
    // data members
    char* Buffer;
    const size_t NumBytes;

    // ctor & dtor
    RaiiBuffer(const size_t n)
        : Buffer( new char[n]() )
        , NumBytes(n)
    { }
    ~RaiiBuffer(void) {
        delete[] Buffer;
    }
    // add'l methods
    void Clear(void) {
        memset(Buffer, 0, NumBytes);
    }
};

} // namespace BamTools

#endif // BAMAUX_H
