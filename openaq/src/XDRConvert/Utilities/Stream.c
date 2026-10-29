
/******************************************************************************
PURPOSE: Stream.c - Define streams for reading and writing to ASCII and
         portable binary (XDR/IEEE-754/MSB) files, pipes and sockets.
         Provides an efficient, portable (multi-platform, cross-language
         compatible), convenient and compatible alternative to fopen, popen
         and socket and associated fread/fwrite/xdr_vector calls.

NOTES:   Portability and Efficiency (and Correctness):

         In general, programs using fread/fwrite are not portable since
         different platforms (and different Fortran compilers) have different
         binary formats. For example, some platforms represent an integer as
         32-bits while others represent it as 64-bits. Also, some platforms
         store the most-significant-byte first (MSB aka "big-endian") while
         others do the reverse (LSB aka "little-endian"). And Fortran compilers
         are notorious for inserting byte-count words around blocks of data
         written and each compiler may use a different scheme.

         This is the reason the XDR (eXternal Data Representation) library was
         invented (by Sun). Using it instead of fread/fwrite/READ/WRITE allows
         data to be shared across different platforms. (Byte-count separators
         are not used by XDR or this library.)

         However, using XDR on platforms whose native binary format happens to
         match XDR (i.e., MSB/IEEE-754) results in an unnecessary performance
         penalty (a factor of about 4 to 10). Also, while the XDR library is
         ubiquitous, implementations on some platforms are defective. Some
         routines in this library that do use XDR detect and correct for these
         defects.

         Using this library instead of fread/fwrite/xdr_vector provides the
         best of both worlds: portability and efficiency (and correctness
         work-arounds).

         Convenience:

         Furthermore, this library provides a higher-level (i.e., simpler, more
         convenient) interface for I/O compared to fread, xdr_vector, socket,
         etc. (See the man pages for these routines for evidence.)

         Cross-platform and Cross-language:

         This library uses only 64-bit integer and real data as routine
         arguments (see Integer and Real from BasicNumerics.h). This is because
         64-bits is the only word size that is implemented across all targeted
         platforms (including Pentium PCs to most workstations to Cray
         Supercomputers) and all targeted languages (C, C++, Fortran-90,
         Eiffel, Java).

         Large Files and Large Data Range:

         Also, 64-bit integers are needed to represent offsets into files
         larger than 2GB. And some integer data such as hydrologic units
         catalog (HUC) values have 14 decimal digits which require a 64-bit
         representation. Using Integer and Real from BasicNumerics.h avoids
         the subtle truncation and alignment defects and subsequent debugging
         associated with porting to new platforms and interfacing with other
         languages.

         Compatability, Efficiency and Functionality:

         Real (floating-point) data is streamed as MSB (most significant byte
         first - "big-endian") IEEE-754 32-bit and 64-bit format just as with
         XDR (xdr_vector/xdr_float/xdr_double). In fact, this library uses
         these XDR routines on non-native IEEE platforms (e.g., Pentium, Cray).
         This is necessary to achieve the required byte-swapping and/or convert
         non-finite values (+/-inf, nan) to a representable form (+/-HUGE_VAL,
         zero).

         On IEEE-754 platforms, real data is streamed using fread/fwrite since
         the result is the same, but it is an order-of-magnitude faster than
         using xdr_vector.

         Therefore, since real data streamed is formatted identically to XDR
         programs using this library can read data written by programs using
         XDR and vice-versa.

         Integer data is streamed as MSB format just as with XDR. However, XDR
         only supports 32-bit integer data - that is, xdr_char/xdr_short/
         xdr_int/xdr_long are all 32-bits. (What a waste of a spec!). In
         contrast, this library supports the streaming of integer data values
         from 64-bits to 8, 16, 32 and 64-bits. It also handles clamping of
         the 64-bit integer values to the targeted range (rather than modulo
         truncation). For example, when writing 64-bit integer data as 32-bits,
         the values outside the range [-2147483648, 2147483647] (2^31) are
         clamped to that range so, the value 4294967297 is written as
         2147483647 (instead of the modulo result 1).

         So while it is compatible with XDR (programs using this library can
         read 32-bit integer data written by programs using XDR and
         vice-versa), this library provides additional functionality.

         Interleaving Calls:

         This library provides routines to access the underlying file pointer
         and descriptor so, while this violates encapsulation, it does provide
         a mechanism to allow the stream to be used by other routines based on
         fread/fwrite (e.g., external libraries that the application uses).
         That is, calls to fread/fwrite/fseek/fcntl, etc. can be used on the
         stream. However, the underlying file pointer/descriptor should not be
         closed since this will result in a 'dangling reference' with ensuing
         'undefined behavior'.

         Convenient Correct Buffering:

         Streaming large arrays of data using lower-level routines such as
         fread and xdr_vector requires checking that the size does not exceed
         the representation of the parameter type (e.g., size_t or unsigned
         int). In such cases, multiple calls must made to avoid the integer
         overflows than can occur when such types are not large enough to
         represent the actual size (e.g., 32-bits are not enough to represent
         more than a 2GB array). This library detects and correctly implements
         the multiple calls in such cases so client code is simpler.

         Similarly, UNIX sockets have a relatively small buffer (e.g., 256kb)
         associated with them which limits the amount of data that can be
         streamed in a single call using read/write on the socket descriptor.
         That is, read will return a short byte count. This library uses
         fread/fwrite which handles this buffering correctly (fread detects
         short counts then uses multiple calls and does not return until it has
         read all that was requested) so clients don't have to be burdened with
         discovering, debugging and implementing such work-arounds.

         Reliable Socket Streaming

         The sockets created with this library use TCP/IP (AF_INET,SOCK_STREAM)
         for reliable streaming. They also allow reusability of ports i.e.,
         the same port can be reused repeatedly so applications can open and
         close connections on-the-fly without running out of available ports.

         Relative Efficiency:

         Approximate relative efficiency cost factor compared to fread/fwrite:

         Routine           Native IEEE  Non-IEEE    Non-IEEE       xdr_vector
                           big-endian   big-endian  little-endian  compared to
                           (SGI)        (Cray C90)  (DEC)          fread/fwrite
         ----------------------------------------------------------------------
         readBytes            1            1           1           -
         read8BitIntegers     6.4          700?        4.1         -
         read16BitIntegers    4.2          600?        3.3         -
         read32BitIntegers    3            500?        2.4         4.4
         read64BitIntegers    1            1           1.7         -
         read32BitReals       2.8          14?         7.3         4.4
         read64BitReals       1            12?         6.3         10

         writeBytes           1            1           1           -
         write8BitIntegers    7            700?        1.3         -
         write16BitIntegers   4.1          600?        1.2         -
         write32BitIntegers   2.9          500?        1           4.8
         write64BitIntegers   1            1           1.1         -
         write32BitReals      2.9          30?         1.4         4.2
         write64BitReals      1            30?         1.2         8.3

         Performance conclusions:

           1. If you need to stream byte data use read/writeBytes.
              Only use read/write8BitIntegers if you need to subsequently
              process the data as integers (e.g., for doing arithmetic).

           2. Use read/write64BitIntegers/Reals to maximize i/o performance
              of numeric data since doing so will avoid allocation of temporary
              buffers and copying when clamping the data during reads and
              expanding the data during writes needed to translate to/from a
              smaller (8, 16, 32-bit) word size. Note that only a 64-bit word
              size, for both intergers and reals, is supported on all platforms
              (e.g., Crays can't represent 16-bit or 32-bit word sizes).

           3. All non-IEEE platforms (e.g., Crays, DEC alpha, Intel Pentium)
              will have suboptimal i/o performance (compared to fread/fwrite)
              of real data since XDR must be used to translate the format
              to/from the IEEE-754 standard (handling Nans/Infinity, and
              clamping).

           4. All little-endian (LSB: least-significant byte first) platforms
              (DEC alpha, Intel Pentium) will have suboptimal i/o performance
              (compared to fread/fwrite) of integer data since byte-swapping
              must be done (in place - i.e., without allocation of temporary
              buffers) to translate the format to/from the standard big-endian
              (MSB: most significant byte first) format.

           5. There is an unexpectedly large penalty on Crays for non-64-bit
              integer data i/o. Perhaps the temporary buffer allocation and/or
              copy/clamp pass throw the job off the fast path? The measurements
              were made in interactive mode. A more thorough investigation is
              warranted.

           6. There is an unexpectedly small writing penalty on the DECs which
              is a mystery since xdr or byte-swapping must be performed which
              should add significant overhead but does not appear to. A more
              thorough investigation is warranted.

HISTORY: 1995/12 Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>        /* For FILE, fopen(), fdopen(), vfprintf(), *getc()*/
#include <ctype.h>        /* For isspace().                                  */
#include <string.h>       /* For strncpy(), memset().                        */
#include <stdarg.h>       /* For va_list, va_start, va_end.                  */
#include <unistd.h>       /* For close(), sleep().                           */
#include <netdb.h>        /* For gethostbyname().                            */
#include <sys/types.h>    /* For struct stat.                                */
#include <sys/stat.h>     /* For struct stat.                                */
#include <sys/socket.h>   /* For socket(), connect(), getsockopt().          */
#include <netinet/in.h>   /* For struct sockaddr, struct sockaddr_in, htons()*/
#include <rpc/rpc.h>      /* For XDR.                                        */
#include <rpc/xdr.h>      /* For XDR.                                        */

#ifdef _CRAY
#include <math.h>         /* HACK: CRAY BUG: For HUGE_VAL, HUGE_VALF.        */
#include <float.h>        /* For FLT_MIN, FLT_MAX to detect/avoid overflow.  */
#define BROKEN_CRAY_XDR 1
#else
#define BROKEN_CRAY_XDR 0
#endif

#if defined(__linux) && defined(__ia64__)
extern int fseeko64(FILE*, unsigned long long, int);
extern unsigned long long ftello64(FILE*);
#endif

#include <Assertions.h>    /* For macros PRE(), POST(), AND*().              */
#include <BasicNumerics.h> /* For Integer, Real.                             */
#include <Failure.h>       /* For failureMessage().                          */
#include <Memory.h>        /* For NEW_ZERO(), FREE_ZERO().                   */
#include <Stream.h>        /* For public interface.                          */

/*================================= MACROS ==================================*/

/* Determine the existence of suitable 64-bit fseek/ftell routines: */

#if _CRAY || __alpha
/*
 * Only ANSI fseek() and ftell() routines exist (taking long arguments) but
 * longs are 64-bits:
 */
#define FSEEK fseek
#define FTELL ftell
#define FSEEK_FTELL_ARE_64_BITS 1
#elif __sgi || __sun || __linux__ || _AIX
/*
 * Non-ANSI fseeko64() and ftello64() routines exist taking 64-bit long long
 * arguments:
 */
#define FSEEK fseeko64
#define FTELL ftello64
#define FSEEK_FTELL_ARE_64_BITS 1
#elif __APPLE__
/*
 * Non-ANSI fseeko() and ftello() routines exist taking 64-bit quad_t
 * arguments:
 */
#define FSEEK fseeko
#define FTELL ftello
#define FSEEK_FTELL_ARE_64_BITS 1
#else
/*
 * Only ANSI fseek() and ftell() routines exist (taking long arguments) and
 * longs are 32-bits:
 */
#define FSEEK fseek
#define FTELL ftell
#define FSEEK_FTELL_ARE_64_BITS 0
#endif

/*
 * HACK for AIX BUG: even though size_t is 64-bits and ftell/fseek/fread
 * accept 64-bit arguments and work with values greater than 2^31,
 * the broken fwrite implementation only uses the lower 31-bits!
 * Also, APPLE BUG: fwrite is broken: an 8GB write results in a 3.7GB file!
 */
#if defined( _AIX ) || defined( __APPLE__ )
#define BROKEN_FWRITE 1
#else
#define BROKEN_FWRITE 0
#endif


/* Define read/write macros that allow optimized (non-xdr) implementations: */

#define READ_ITEM(data_,type_,item_) \
  (USES_NATIVE_XDR(type_) ? \
    fread( (item_), sizeof (type_), 1, (data_)->file ) == 1 \
  : xdr_##type_( &(data_)->xdr, (item_) ) )

#define WRITE_ITEM(data_,type_,item_) \
  (USES_NATIVE_XDR(type_) ? \
    fwrite( (item_), sizeof (type_), 1, (data_)->file ) == 1 \
  : xdr_##type_( &(data_)->xdr, (item_) ) )

#define READ_ITEMS(data_,type_,count_,items_) \
  (USES_NATIVE_XDR(type_) ? \
    fread( (items_), sizeof (type_), (count_), (data_)->file ) == (count_) \
  : xdr_vector( &(data_)->xdr, (char*) (items_), (count_), sizeof (type_), \
                (xdrproc_t) xdr_##type_ ) )

#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b)?(a):(b))

/*================================== TYPES ==================================*/

enum { FILE_STREAM, PIPE_STREAM, SOCKET_STREAM }; /* Stream->type. */

struct StreamData {
  char*   name;        /* Pathed file or command or host name.              */
  char*   mode;        /* I/O mode: "r", "w", "a", "r+", etc.               */
  Integer ok;          /* Did last command succeed?                         */
  Integer type;        /* FILE_STREAM, PIPE_STREAM, SOCKET_STREAM           */
  Integer port;        /* Port number, if socket.                           */
  enum xdr_op xdrMode; /* XDR_DECODE when reading, XDR_ENCODE when writing. */
  FILE*   file;        /* Stream file structure.                            */
  XDR     xdr;         /* XDR buffer structure for encoding/decoding.       */
};

/*
 * Type of routine to copy/clamp integers to subranges.
 * See clampTo8BitInteger(), clampTo16BitInteger(), clampTo32BitInteger()
 * which are called implicitly by writeClampedCopy().
 */

typedef void (*Clamper)( const Integer*, Integer, unsigned char* );

/*========================== FORWARD DECLARATIONS ===========================*/

/* Commands: */

static void closeStream( Stream* self ); /* Destructor called by ->free(). */

static void flushStream(      Stream* self );
static void seekFromStart(    Stream* self, Integer byteOffset );
static void seekFromEnd(      Stream* self, Integer byteOffset );
static void seekFromCurrent(  Stream* self, Integer byteOffset );

static void readString( Stream* self, char* s, Integer n );
static void readWord(   Stream* self, char* s, Integer n );
static void readByte(   Stream* self, void* x );

static void read8BitInteger(    Stream* self, Integer* x );
static void read16BitInteger(   Stream* self, Integer* x );
static void read32BitInteger(   Stream* self, Integer* x );
static void read64BitInteger(   Stream* self, Integer* x );

static void read32BitReal(      Stream* self, Real*    x );
static void read64BitReal(      Stream* self, Real*    x );

static void readBytes(          Stream* self, void*    a, Integer n );
static void readUpToNBytes(     Stream* self, void*    a, Integer n,
                                Integer* actualCount );

static void read8BitIntegers(   Stream* self, Integer* a, Integer n );
static void read16BitIntegers(  Stream* self, Integer* a, Integer n );
static void read32BitIntegers(  Stream* self, Integer* a, Integer n );
static void read64BitIntegers(  Stream* self, Integer* a, Integer n );

static void read32BitReals(     Stream* self, Real*    a, Integer n );
static void read64BitReals(     Stream* self, Real*    a, Integer n );

static void writeString(        Stream* self, const char* s, ... );
static void writeByte(          Stream* self, unsigned char x );

static void write8BitInteger(   Stream* self, Integer x );
static void write16BitInteger(  Stream* self, Integer x );
static void write32BitInteger(  Stream* self, Integer x );
static void write64BitInteger(  Stream* self, Integer x );

static void write32BitReal(     Stream* self, Real    x );
static void write64BitReal(     Stream* self, Real    x );

static void writeBytes(         Stream* self, const void*    a, Integer n );

static void write8BitIntegers(  Stream* self, const Integer* a, Integer n );
static void write16BitIntegers( Stream* self, const Integer* a, Integer n );
static void write32BitIntegers( Stream* self, const Integer* a, Integer n );
static void write64BitIntegers( Stream* self, const Integer* a, Integer n );

static void write32BitReals(    Stream* self, const Real*    a, Integer n );
static void write64BitReals(    Stream* self, const Real*    a, Integer n );

/* Queries: */

static Integer     invariant(  const Stream* self );
static Integer     ok(         const Stream* self );
static Integer     isReadable( const Stream* self );
static Integer     isWritable( const Stream* self );
static Integer     isSeekable( const Stream* self );
static Integer     isBlocking( const Stream* self );
static Integer     isAtEnd(    const Stream* self );

static Integer     offset(     const Stream* self );
static Integer     size(       const Stream* self );
static const char* name(       const Stream* self );
static FILE*       file(       const Stream* self );
static Integer     descriptor( const Stream* self );


/* Helpers: */

static void deallocate( Stream* self );

static Stream* newSocketStream( Integer port, const char* host );

static const char* typeName( const Stream* self );

CHECKING( static Integer isReadMode( const Stream* self ); )

CHECKING( static Integer isWriteMode( const Stream* self ); )

static void ensureReadMode( Stream* self );

static void ensureWriteMode( Stream* self );

static void seekStream( Stream* self, Integer byteOffset, Integer whence );

static void streamArrayBuffered( Stream* self,
                                 Integer count,
                                 Integer sizeOfItem,
                                 xdrproc_t elproc,
                                 /* const */ void* array );

static void writeClampedCopy( Stream* self, const Integer* a, Integer n,
                              Integer bytesPerItem, Clamper clamper,
                              const char* kindOfValues );

static void checkAndReport( const Stream* self, const char* readOrWrite,
                            Integer count, const char* dataType );

static Stream* allocate( Integer nameLength, Integer modeLength );


/* Generic integer copy/clamp routines: */

static void clampTo8BitInteger( const Integer* src, Integer count,
                                unsigned char* dst );

static void clampTo16BitInteger( const Integer* src, Integer count,
                                 unsigned char* dst );

static void clampTo32BitInteger( const Integer* src, Integer count,
                                 unsigned char* dst );

static void clampTo64BitInteger( const Integer* src, Integer count,
                                 unsigned char* dst );


/* Generic socket routines: */

static FILE* createSocket( Integer port, const char* host );

static Integer createServerSocket( Integer port );

static Integer createClientSocket( Integer port, const char* host );

static Integer enablePortReusability( Integer theSocket );

static Integer setSocketBufferSize( Integer theSocket );

static Integer establishSocketProperty( Integer theSocket,
                                        Integer property,
                                        Integer resetValue,
                                        const char* propertyName );

static FILE* fileOfSocket( Integer theSocket );


/* Generic string routines: */

static Integer matchesWord( const char* target, const char* candidates );


/* Generic byte-swapping routines: */

static Integer swapped8Bytes( const Integer x );

static void swap8Bytes( Integer* dst, const Integer* src, Integer count );


#ifdef _CRAY
/* HACK: CRAY BUG: */
static void convertXDRRead32BitNonFinite( float* x );
static void convertXDRRead64BitNonFinite( double* x );
#endif


/*============================ PUBLIC FUNCTIONS =============================*/


/******************************************************************************
PURPOSE: newFileStream - Open a file for reading and/or writing.
INPUTS:  const char* fileName  The name of the file to open or one of:
                               "-stdin", "-stdout", "-stderr".
         const char* mode      Valid and appropriate fopen() mode:
                               Read: "r", "rb"
                               Truncate or create: "w", "wb".
                               Append: "a", "ab".
                               Read or write: "r+", "r+b", "rb+".
                               Truncate/create for update: "w+", "w+b", "wb+".
                               Append or update at EOF: "a+", "a+b", "ab+".
                               Note that the file stream will be opened for
                               binary mode regardless of whether or not the
                               given mode includes 'b'.
RETURNS: Stream* Initialized stream structure, or 0 if failed.
NOTES:   If unsuccessful then failureMessage() is called and 0 is returned.

         Example:

           Stream* stream = newFileStream( "/home/plessel/.login", "r" );

           if ( stream ) {

             while ( stream->ok( stream ) && ! stream->isAtEnd( stream ) ) {

               char string[ 80 ];

               stream->readString( stream, string, 80 );

               if ( stream->ok( stream ) ) {
                 fprintf( stderr, "%s\n", string );
               }
             }

             FREE_OBJECT( stream );
           }

******************************************************************************/

Stream* newFileStream( const char* fileName, const char* mode ) {
  PRE06( fileName, *fileName, mode,
         matchesWord(mode," r rb w wb a ab r+ r+b rb+ w+ w+b wb+ a+ a+b ab+ "),
         IMPLIES( strcmp( fileName, "-stdin" ) == 0,
                  matchesWord( mode, " r rb " ) ),
         IMPLIES( matchesWord( fileName, " -stdout -stderr " ), *mode != 'r'));

  char binaryMode[ 4 ]; /* To ensure a binary mode. */
  Stream* result = 0;

  /* Generate a mode that includes 'b' for binary i/o: */

  if ( strchr( mode, 'b' ) ) {
    CHECK( IN3( strlen( mode ), 2, 3 ) );
    strncpy( binaryMode, mode, 4 );
  } else {
    CHECK( IN3( strlen( mode ), 1, 2 ) );
    strncpy( binaryMode, mode, 3 );
    strncat( binaryMode, "b", 2 );
  }

  CHECK2( IN3( strlen( binaryMode ), 2, 3 ), strchr( binaryMode, 'b' ) );
  result = allocate( strlen( fileName ), strlen( binaryMode ) );

  if ( result ) {

    StreamData* const data = result->data;

    if ( strcmp( fileName, "-stdin"  ) == 0 ) {
      data->file = stdin;
    } else if ( strcmp( fileName, "-stdout" ) == 0 ) {
      data->file = stdout;
    } else if ( strcmp( fileName, "-stderr" ) == 0 ) {
      data->file = stderr;
    } else {
      data->file = fopen( fileName, binaryMode );
    }

    if ( data->file ) {
      strcpy( data->name, fileName );
      strcpy( data->mode, binaryMode );
      data->ok         = 1;
      data->type       = FILE_STREAM;
      data->xdrMode    = *mode == 'r' ? XDR_DECODE : XDR_ENCODE;
      xdrstdio_create( &data->xdr, data->file, data->xdrMode );
    } else {
      deallocate( result );
      result = 0;
      failureMessage( "Can't open file '%s' for %s.", fileName,
                      *mode == 'r' ? "reading" : "writing" );
    }
  }

/*
 * Must bypass invariant in post-condition too since result may be 0.
 * However, in cases where result is not 0, the invariant will be evaluated
 * during calls to ok() anyway.
 */

  POST0( IMPLIES( result,
                  AND6( ok( result ),
                        IMPLIES_ELSE( strcmp( fileName, "-stdin" ) == 0,
                                      AND2( isReadable( result ),
                                            isBlocking( result ) ),
                                      IMPLIES( isReadable( result ),
                                               ! isBlocking( result ) ) ),
                        IMPLIES( *mode == 'r', isReadable( result ) ),
                        IMPLIES( OR3( IN3( mode[ 0 ], 'w', 'a' ),
                                      mode[ 1 ] == '+',
                                      AND2( mode[ 1 ] != '\0',
                                            mode[ 2 ] == '+' ) ),
                                 isWritable( result ) ),
                        OR2( isReadable( result ), isWritable( result ) ),
                        IMPLIES( ! matchesWord( fileName,
                                     " -stdin -stdout -stderr /dev/null " ),
                                 isSeekable( result ) ) ) ) );
  return result;
}


/******************************************************************************
PURPOSE: newPipeStream - Open a pipe for reading or writing.
INPUTS:  const char* command   Command to launch the process to open a pipe to.
         const char* mode      Valid and appropriate popen() mode: "r" or "w".
RETURNS: Stream*               Initialized stream structure, or 0 if failed.
NOTES:   If unsuccessful then failureMessage() is called and 0 is returned.
         Blocks on FREE_OBJECT() until connected-to process exits.
         Example:

           Stream* stream = newPipeStream( "ls", "r" );

           if ( stream ) {

             while ( stream->ok( stream ) && ! stream->isAtEnd( stream ) ) {

               char string[ 10 ];

               stream->readString( stream, string, 10 );

               if ( stream->ok( stream ) ) {
                 fprintf( stderr, "%s\n", string );
               }
             }

             FREE_OBJECT( stream );
           }

******************************************************************************/

Stream* newPipeStream( const char* command, const char* mode ) {
  PRE04( command, *command, mode, matchesWord( mode, " r w " ) );

  Stream* result = allocate( strlen( command ), strlen( mode ) );

  if ( result ) {

    StreamData* const data = result->data;

    /*
     * Must flush stdout & stderr to sync subsequent output in cases where
     * pipe command outputs to stdout or stderr. E.g., "cat", "echo", etc.
     */

    fflush( stdout );
    fflush( stderr );
    data->file = popen( command, mode );

    if ( data->file ) {
      strcpy( data->name, command );
      strcpy( data->mode, mode );
      data->ok         = 1;
      data->type       = PIPE_STREAM;
      data->xdrMode    = *mode == 'r' ? XDR_DECODE : XDR_ENCODE;
      xdrstdio_create( &data->xdr, data->file, data->xdrMode );
    } else {
      deallocate( result );
      result = 0;
      failureMessage( "Can't open pipe with command '%s' for %s.", command,
                      *mode == 'r' ? "reading" : "writing" );
    }
  }

/*
 * Must bypass invariant in post-condition too since result may be 0.
 * However, in cases where result is not 0, the invariant will be evaluated
 * during calls to ok() anyway.
 */

  POST0( IMPLIES( result, AND4( ok( result ),
                                IMPLIES( isReadable( result ),
                                         isBlocking( result ) ),
                                ! isSeekable( result ),
                                IMPLIES_ELSE( *mode == 'r',
                                              isReadable( result ),
                                              isWritable( result ) ) ) ) );
  return result;
}


/******************************************************************************
PURPOSE: newServerSocketStream - Create a server socket for reading & writing.
INPUTS:  Integer     port  The port number to use. (E.g., 4321).
RETURNS: Stream*           Initialized file structure, or 0 if failed.
NOTES:   If unsuccessful then failureMessage() is called and 0 is returned.
         This routine blocks awaiting a (single) connection to another
         process - the "client" using the same port number. Launch the server
         first then launch the client (with a delay, sleep( 3 )). This allows
         the server time to establish the socket on the port before the client
         attempts to connect to it (otherwise the client's attempt to connect
         will fail). Multi-client processing requires separate socket calls
         (with unique ports) for each client.

         Example:

           server:

             Stream* theSocket = newServerSocketStream( 4321 );

             if ( theSocket ) {
               readFromClient( theSocket );
             }

           client:

             sleep( 3 );
             theSocket = newClientSocketStream( 4321, "localhost" );

             if ( theSocket ) {
               writeToServer( theSocket );
             }

******************************************************************************/

Stream* newServerSocketStream( Integer port ) {
  PRE0( IN_RANGE( port, 1, 65535 ) );

  Stream* const result = newSocketStream( port, 0 );

  POST0( IMPLIES( result, AND5( ok( result ),
                                isReadable( result ),
                                isWritable( result ),
                                isBlocking( result ),
                                ! isSeekable( result ) ) ) );
  return result;
}


/******************************************************************************
PURPOSE: newClientSocketStream - Create a client socket for reading & writing.
INPUTS:  Integer     port  The port number to use. (E.g., 4321).
         const char* host  Name of host to connect to, e.g., "localhost",
                           "barney", "fred.xyz.epa.gov".
RETURNS: Stream*           Initialized file structure, or 0 if failed.
NOTES:   If unsuccessful then failureMessage() is called and 0 is returned.
         This routine blocks awaiting a (single) connection to another
         process - the "server" using the same port number. Launch the server
         first then launch the client (with a delay, sleep( 3 )). This allows
         the server time to establish the socket on the port before the client
         attempts to connect to it (otherwise the client's attempt to connect
         will fail). Multi-client processing requires separate socket calls
         (with unique ports) for each client.

         Example:

           server:

             Stream* theSocket = newServerSocketStream( 4321 );

             if ( theSocket ) {
               readFromClient( theSocket );
             }

           client:

             sleep( 3 );
             theSocket = newClientSocketStream( 4321, "localhost" );

             if ( theSocket ) {
               writeToServer( theSocket );
             }

******************************************************************************/

Stream* newClientSocketStream( Integer port, const char* host ) {
  PRE03( IN_RANGE( port, 1, 65535 ), host, *host );

  Stream* const result = newSocketStream( port, host );

  POST0( IMPLIES( result, AND5( ok( result ),
                                isReadable( result ),
                                isWritable( result ),
                                isBlocking( result ),
                                ! isSeekable( result ) ) ) );
  return result;
}


/******************************************************************************
PURPOSE: closeStream - Destruct (close and deallocate) a stream.
INPUTS:  Stream* self  The stream to close.
OUTPUTS: Stream* self  Deallocated and zeroed stream structure.
NOTES:   Use FREE_OBJECT() macro instead - it calls this routine indirectly,
         but also zeros-out the argument to avoid dangling references to
         deallocated memory.
         Also flushes the file, though for failure detection, flush() should be
         explicitly called first and ok() checked before calling FREE_OBJECT().
         Blocks on pipe streams until the other process exits.
         If called system routines (e.g., pclose(), fclose()) return "failure
         codes" then a message is printed to stderr in lieu of calling
         failureMessage(), since destructors are not allowed to fail.
******************************************************************************/

static void closeStream( Stream* self ) {

  PRE( self );

  const char* const type = typeName( self );
  StreamData* const data = self->data;

  data->ok = 1;
  fflush( stdout );
  fflush( stderr );
  fflush( data->file );

  if ( data->type == PIPE_STREAM ) {
    const Integer status = pclose( data->file ); /* To dbx check exit status.*/
    data->ok = status != -1;
  } else if ( ! matchesWord( data->name, " -stdin -stdout -stderr " ) ) {
    data->ok = fclose( data->file ) == 0;
  }

  data->file = 0; /* Invalidates invariant. */

  if ( ! data->ok ) {
    fprintf( stderr, "\n\n\aWarning: Abnormal close of %s '%s' detected.\n"
             "Some data may not have been completely written.\n\n",
             type, data->name );
  }

  xdr_destroy( &data->xdr ); /* Free XDR buffer. */
  deallocate( self );
  self = 0;

  POST0( self == 0 );
}


/******************************************************************************
PURPOSE: flushStream - Flush the output buffer associated with a stream.
INPUTS:  Stream* self   The stream to flush.
OUTPUTS: Stream* self   Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called.
         This should be called (on stream opened for writing) before
         close() so that ok() can be used to determine if the last
         write succeeded. Also flushes stdout and stderr first.
******************************************************************************/

static void flushStream( Stream* self ) {

  PRE( isWritable( self ) );

  StreamData* const data = self->data;

  fflush( stdout );
  fflush( stderr );
  data->ok = fflush( data->file ) == 0;

  if ( ! data->ok ) {
    failureMessage( "Can't flush file '%s'.", data->name );
  }

  POST( isWritable( self ) );
}


/******************************************************************************
PURPOSE: seekFromStart - Seek to a given byte offset from the start of a file.
INPUTS:  Stream* self        The stream to seek on.
         Integer byteOffset  The offset to seek by.
OUTPUTS: Stream* self        Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called.
******************************************************************************/

static void seekFromStart( Stream* self, Integer byteOffset ) {

  PRE2( isSeekable( self ), byteOffset >= 0 );
  CHECKING( Integer OLD( offset ) = offset( self ); )

  seekStream( self, byteOffset, SEEK_SET );

  POST2( isSeekable( self ),
         IMPLIES_ELSE( ok( self ),
                       offset( self ) == byteOffset,
                       offset( self ) == OLD( offset ) ) );
}


/******************************************************************************
PURPOSE: seekFromEnd - Seek to a given byte offset from the end of a file.
INPUTS:  Stream* self        The stream to seek on.
         Integer byteOffset  The offset to seek by.
OUTPUTS: Stream* self        Updated stream structure.
NOTES:   If unsuccessful (e.g., attempt to seek beyond the start of a file)
         then failureMessage() is called.
         If offset is positive and a write occurs later then intermediate
         data between the previous end-of-file and the new end-of-file will
         be zero'd by the write (man fseek).
******************************************************************************/

static void seekFromEnd( Stream* self, Integer byteOffset ) {

  PRE( isSeekable( self ) );
  CHECKING( Integer OLD( offset ) = offset( self ); )

  seekStream( self, byteOffset, SEEK_END );

  POST2( isSeekable( self ),
         IMPLIES_ELSE( ok( self ),
                       offset( self ) == size( self ) + byteOffset,
                       offset( self ) == OLD( offset ) ) );
}


/******************************************************************************
PURPOSE: seekFromCurrent - Seek to a given byte offset from the current offset
         within a file.
INPUTS:  Stream* self        The stream to seek on.
         Integer byteOffset  The offset to seek by.
OUTPUTS: Stream* self        Updated stream structure.
NOTES:   If unsuccessful (e.g., attempt to seek beyond the start of a file)
         then failureMessage() is called.
******************************************************************************/

static void seekFromCurrent( Stream* self, Integer byteOffset ) {

  PRE( isSeekable( self ) );
  CHECKING( Integer OLD( offset ) = offset( self ); )

  seekStream( self, byteOffset, SEEK_CUR );

  POST2( isSeekable( self ),
         IMPLIES_ELSE( ok( self ),
                       offset( self ) == (OLD( offset ) + byteOffset),
                       offset( self ) == OLD( offset ) ) );
}


/******************************************************************************
PURPOSE: readString - Read a one-line string from a stream.
INPUTS:  Stream* self  The stream to read from.
         Integer n     The size of s[].
OUTPUTS: char    s[n]  The string read.
         Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called and *s is set to '\0'.
         Reading stops when n - 1 characters are read or a newline character
         is read (and stored in s) - just like fgets().
******************************************************************************/

static void readString( Stream* self, char s[], Integer n ) {

  PRE3( isReadable( self ), s, n > 0 );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  data->ok = 0;
  *s = '\0';

  ensureReadMode( self );

  if ( n <= INT_MAX ) { /* One read ok since within fgets() limit: */
    data->ok = fgets( s, n, data->file ) != 0;
    checkAndReport( self, "read", n, "size string" );
  } else {
    /* Must loop over multiple buffered reads: */
    const Integer maxCharsReadable = INT_MAX / sizeof (char);
    Integer charsRead = 0; /* Number of chars already read so far. */
    Integer done = 0;

    CHECK( IN_RANGE( maxCharsReadable, 1, INT_MAX ) );

    do {
      const Integer maxCharsRemaining = n - charsRead;

      const Integer maxCharsToRead    = maxCharsRemaining < maxCharsReadable ?
                                         maxCharsRemaining : maxCharsReadable;

      char* const partialString = s + charsRead;

      CHECK( maxCharsRemaining >= 0 );
      CHECK2( maxCharsToRead > 0, maxCharsToRead * sizeof (char) <= INT_MAX );

      data->ok = fgets( partialString, maxCharsToRead, data->file ) != 0;

      if ( ! data->ok ) {
        done = 1;

        if ( charsRead == 0 ) /* No chars were ever read. */ {
          failureMessage( "Can't read up to %"INTEGER_FORMAT" characters of"
                          " a partial string from '%s'.",
                          maxCharsToRead - 1, data->name );
        } else {
          data->ok = 1;  /* False alarm. Return short s w/o '\n' */
        }
      } else {
        const Integer partialStringLength = strlen( partialString );

        CHECK( partialStringLength >= 0 ); /* Might read a string w/ '\0'? */
        charsRead += partialStringLength;  /* In which case we're done.    */
        done = OR2( charsRead == n, partialStringLength == 0 );
      }
    } while ( ! done );
  }

  if ( ! data->ok ) {
    *s = '\0';
  }

  POST4( isReadMode( self ),
         strlen( s ) < n,
         IMPLIES( OR2( n == 1, ! ok( self ) ), *s == '\0' ),
         IMPLIES( AND2( isSeekable( self ), ok( self ) ),
                  offset( self ) >= OLD( offset ) + strlen( s ) ) );
}


/******************************************************************************
PURPOSE: readWord - Read a whitespace-delimited word from a stream.
INPUTS:  Stream* self  The stream to read from.
         Integer n     The capacity
          of s, the string to read into.
OUTPUTS: char    s[n]  The string read.
         Stream* self  Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called and s[ 0 ] is '\0'.
         Skips leading whitespace characters - like scanf( "%5s", s ).
         Reading stops when n - 1 characters are stored or a whitespace
         character is read (and putback) or EOF is reached.
******************************************************************************/

static void readWord( Stream* self, char s[], Integer n ) {

  PRE3( isReadable( self ), s, n > 0 );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  ensureReadMode( self );

  if ( n == 1 ) {
    *s = '\0';
    data->ok = 1;
  } else {
    char format[ 80 ];
    /* HACK DEC BUG: snprintf() is missing! */
    CHECKING( const int formatLength = )
      sprintf( format, "%c%"INTEGER_FORMAT"s", '%', n - 1 );

    CHECK( IN_RANGE( formatLength, 3, sizeof format - 2 ) ); /* Untruncated. */
    data->ok = fscanf( data->file, format, s ) == 1;

    if ( ! data->ok ) {
      *s = '\0';
    }

    checkAndReport( self, "read", 1, "word" );
  }

  POST4( isReadMode( self ),
         strlen( s ) < n,
         IMPLIES( OR2( n == 1, ! ok( self ) ), *s == '\0' ),
         IMPLIES( AND2( isSeekable( self ), ok( self ) ),
                  offset( self ) >= OLD( offset ) + strlen( s ) ) );
}


/******************************************************************************
PURPOSE: readByte - Read a byte from a stream.
INPUTS:  Stream* self  The stream to read from.
OUTPUTS: void* x       The byte read.
         Stream* self  Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called and x is '\0'.
******************************************************************************/

static void readByte( Stream* self, void* x ) {

  PRE2( isReadable( self ), x );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  char* c = (char*) x;

  ensureReadMode( self );
  *c = '\0'; /* Indirectly initialize output in case fread() fails. */
  data->ok = fread( x, 1, 1, data->file ) == 1;
  checkAndReport( self, "read", 1, "byte" );

  POST3( isReadMode( self ),
         IMPLIES( ! ok( self ), *((char*) x) == '\0' ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 1,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 1) ) ) );
}


/******************************************************************************
PURPOSE: read8BitInteger - Read a signed 8-bit integer.
INPUTS:  Stream* self  The stream to read from.
OUTPUTS: Integer* x    The 8-bit integer read.
         Stream* self  Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called and x is zero.
******************************************************************************/

static void read8BitInteger( Stream* self, Integer* x ) {

  PRE2( isReadable( self ), x );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  signed char c = 0;

  ensureReadMode( self );
  data->ok = fread( &c, 1, 1, data->file ) == 1;
  *x = c;
  checkAndReport( self, "read", 1, "8-bit integer" );

  POST4( isReadMode( self ),
         IN_RANGE( *x, -128, 127 ),
         IMPLIES( ! ok( self ), *x == 0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 1,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 1 ) ) ) );
}


/******************************************************************************
PURPOSE: read16BitInteger - Read a big-endian signed 16-bit integer.
INPUTS:  Stream* self  The stream to read from.
OUTPUTS: Integer* x    The value read.
         Stream* self  Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called and x is zero.
******************************************************************************/

static void read16BitInteger( Stream* self, Integer* x ) {

  PRE2( isReadable( self ), x );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  unsigned char bytes[ 2 ];

  *x = 0;
  ensureReadMode( self );
  data->ok = fread( bytes, 2, 1, data->file ) == 1;

  if ( data->ok ) {
    const unsigned char lowByte  = bytes[ 1 ]; /* Prevent sign-bit extension.*/
    const   signed char highByte = bytes[ 0 ]; /* Cause   sign-bit extension.*/

    *x = ( highByte << 8 ) + lowByte;

    DEBUG( fprintf( stderr, "highByte = %d (%x), lowByte = %d (%x), "
                   "64-bit value = %"INTEGER_FORMAT"\n",
                   highByte, highByte, lowByte, lowByte, *x ); )
  } else {
    checkAndReport( self, "read", 1, "16-bit integer" );
  }

  POST4( isReadMode( self ),
         IN_RANGE( *x, -32768, 32767 ),
         IMPLIES( ! ok( self ), *x == 0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 2,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 2 ) ) ) );
}


/******************************************************************************
PURPOSE: read32BitInteger - Read a big-endian signed 32-bit integer.
INPUTS:  Stream*  self  The stream to read from.
OUTPUTS: Integer* x     The value read.
         Stream*  self  Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called and x is zero.
******************************************************************************/

static void read32BitInteger( Stream* self, Integer* x ) {

  PRE2( isReadable( self ), x );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data  = self->data;
  unsigned char bytes[ 4 ];

  *x = 0;
  ensureReadMode( self );
  data->ok = fread( bytes, 4, 1, data->file ) == 1;

  if ( data->ok ) {
    const unsigned char byte1 = bytes[ 3 ]; /* Prevent sign-bit extension. */
    const unsigned char byte2 = bytes[ 2 ]; /* Prevent sign-bit extension. */
    const unsigned char byte3 = bytes[ 1 ]; /* Prevent sign-bit extension. */
    const   signed char byte4 = bytes[ 0 ]; /* Cause   sign-bit extension. */

    *x = ( byte4 << 24 ) + ( byte3 << 16 ) + ( byte2 << 8 ) + byte1;

    DEBUG( fprintf( stderr, "byte4 = %x, byte3 = %x, byte2 = %x, byte1 = %x, "
                   "64-bit value = %"INTEGER_FORMAT"\n",
                    byte4, byte3, byte2, byte1, *x ); )
  } else {
    checkAndReport( self, "read", 1, "32-bit integer" );
  }

  POST4( isReadMode( self ),
         IN_RANGE( *x, -2147483647 - 1, 2147483647 ),
         IMPLIES( ! ok( self ), *x == 0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 4,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 4 ) ) ) );
}


/******************************************************************************
PURPOSE: read64BitInteger - Read a big-endian signed 64-bit integer.
INPUTS:  Stream*  self  The stream to read from.
OUTPUTS: Integer* x     The value read.
         Stream*  self  Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called and x is zero.
******************************************************************************/

static void read64BitInteger( Stream* self, Integer* x ) {

  PRE2( isReadable( self ), x );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  ensureReadMode( self );
  data->ok = fread( x, 8, 1, data->file ) == 1;

  if ( ! data->ok ) {
    *x = 0;
    checkAndReport( self, "read", 1, "64-bit integer" );
  } else if ( IS_LITTLE_ENDIAN ) {
    *x = swapped8Bytes( *x );
  }

  DEBUG( fprintf( stderr, "value = %"INTEGER_FORMAT"\n", *x ); )

  POST4( isReadMode( self ),
        IN_RANGE(*x,INTEGER_CONSTANT(-9223372036854775807)-INTEGER_CONSTANT(1),
                     INTEGER_CONSTANT( 9223372036854775807) ),
         IMPLIES( ! ok( self ), *x == 0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 8,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 8 ) ) ) );
}


/******************************************************************************
PURPOSE: read32BitReal - Read a (big-endian) IEEE-754 32-bit floating-point
         value and convert it to a (64-bit) Real.
INPUTS:  Stream* self  The stream to read from.
OUTPUTS: Real* x       The value read and converted to a Real.
         Stream* self  Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called and x is zero.
         On IEEE platforms, x may be Nan, -Inf or +Inf otherwise
         non-finite values may be returned as 0, -HUGE_VAL, or +HUGE_VAL.
         See man pages on XDR and IEEE floating-point.
         Uses fread() on native IEEE platforms and XDR (xdr_float) elsewhere.
******************************************************************************/

static void read32BitReal( Stream* self, Real* x ) {

  PRE2( isReadable( self ), x );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  float f = 0.0;

  *x = 0.0;
  ensureReadMode( self );
  data->ok = READ_ITEM( data, float, &f );

#ifdef _CRAY
  if ( data->ok && BROKEN_CRAY_XDR ) {
    convertXDRRead32BitNonFinite( &f );
  }
#endif

  if ( data->ok ) {
    *x = f;
  } else {
    checkAndReport( self, "read", 1, "32-bit real" );
  }

  POST3( isReadMode( self ),
         IMPLIES( ! ok( self ), *x == 0.0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 4,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 4 ) ) ) );
}


/******************************************************************************
PURPOSE: read64BitReal - Read a (big-endian) IEEE-754 64-bit floating-point
         value.
INPUTS:  Stream* self  The stream to read from.
OUTPUTS: Real* x       The value read.
         Stream* self  Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called and x is zero.
         On IEEE platforms, x may be Nan, -Inf or +Inf otherwise
         non-finite values may be returned as 0, -HUGE_VAL, or +HUGE_VAL.
         See man pages on XDR and IEEE floating-point.
         Uses fread() on native IEEE platforms and XDR (xdr_double) elsewhere.
******************************************************************************/

static void read64BitReal( Stream* self, Real* x ) {

  PRE2( isReadable( self ), x );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  double d = 0.0;

  *x = 0.0;
  ensureReadMode( self );
  data->ok = READ_ITEM( data, double, &d );

#ifdef _CRAY
  if ( data->ok && BROKEN_CRAY_XDR ) {
    convertXDRRead64BitNonFinite( &d );
  }
#endif

  if ( data->ok ) {
    *x = d;
  } else {
    checkAndReport( self, "read", 1, "64-bit real" );
  }

  POST3( isReadMode( self ),
         IMPLIES( ! ok( self ), *x == 0.0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 8,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 8 ) ) ) );
}


/******************************************************************************
PURPOSE: readBytes - Read a block of bytes.
INPUTS:  Stream* self  The stream to read from.
         Integer n     The number of bytes to read.
OUTPUTS: void*   a     The block of bytes read.
         Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called and byte a[0] is zero.
******************************************************************************/

static void readBytes( Stream* self, void* a, Integer n ) {

  PRE3( isReadable( self ), a, n > 0 );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  Integer actualCount = 0;

  readUpToNBytes( self, a, n, &actualCount );
  data->ok = actualCount == n;

  if ( ! data->ok ) {
    *((char*) a) = '\0';
  }

  checkAndReport( self, "read", n, "bytes" );

  POST3( isReadMode( self ),
         IMPLIES( ! ok( self ), *((char*) a) == '\0' ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + n ) ) ) );
}


/******************************************************************************
PURPOSE: readUpToNBytes - Read a block of up to n bytes.
INPUTS:  Stream*  self         The stream to read from.
         Integer  n            The maximum number of bytes to read.
OUTPUTS: void*    a            The block of bytes read.
         Integer* actualCount  The actual number of bytes read.
         Stream*  self         Updated file structure.
RETURNS: None.
NOTES:   This is used when the number of bytes expected is unknown, e.g., as
         with pipes, and may be less than a given maximum without indicating a
         failure. If used with a socket, blocking will occur until exactly 'n'
         bytes are read. This routine does not fail - i.e., call
         failureMessage().
         An actualCount value of 0 may indicate an end-of-file condition for a
         file or pipe. If actualCount is zero then the first byte a[0] is zero.
******************************************************************************/

static void readUpToNBytes( Stream* self, void* a, Integer n,
                            Integer* actualCount ) {

  PRE4( isReadable( self ), a, n > 0, actualCount );
  CHECKING( const Integer previousOffset = isSeekable(self) ? offset(self): 0;)
  CHECKING( const Integer OLD(offset) = previousOffset; )

  StreamData* const data = self->data;
  Integer bytesRead = 0; /* Number of bytes actually read. */
  unsigned char* const c = (unsigned char*) a; /* For init and arithmetic. */

  *c = '\0';
  *actualCount = 0;
  ensureReadMode( self );

  if ( isSizet( n ) ) { /* One read will work since within fread() limit: */
    bytesRead = fread( a, 1, n, data->file );
  } else {
    /* Multiple reads may be necessary: */
    /* Use largest buffer size allowed by fread(): */

    const Integer bufferSize = sizeof (size_t) < 8 ? SIZET_MAX : INTEGER_MAX;

    Integer done = 0;

    CHECK2( bufferSize > 0, isSizet( bufferSize ) );

    do {
      const Integer bytesRemaining = n - bytesRead;

      const Integer bytesToRead    =   bytesRemaining < bufferSize ?
                                       bytesRemaining
                                     : bufferSize;

      const Integer charsAlreadyRead = bytesRead / sizeof (char);

      const Integer newBytesRead = fread( c + charsAlreadyRead, 1,
                                          bytesToRead, data->file );

      CHECK( isSizet( bytesToRead ) ); /* Too late... */
      bytesRead += newBytesRead;
      done = OR2( newBytesRead == 0, bytesRead == n );
    } while ( ! done );
  }

  *actualCount = bytesRead;
  data->ok = 1; /* This routine never fails. */

  POST5( isReadMode( self ),
         ok( self ),
         IN_RANGE( *actualCount, 0, n ),
         IMPLIES( *actualCount == 0, *((char*) a) == '\0' ),
         IMPLIES( isSeekable( self ),
                  offset( self ) == OLD( offset ) + *actualCount ) );
}


/******************************************************************************
PURPOSE: read8BitIntegers - Read an array of 8-bit signed integers and
         expand them into 64-bit integers.
INPUTS:  Stream*  self  The stream to read from.
         Integer  n     The number of values to read.
OUTPUTS: Integer* a     The array of integers read.
         Stream*  self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called and a[ 0 ] is zero.
******************************************************************************/

static void read8BitIntegers( Stream* self, Integer* a, Integer n ) {

  PRE3( isReadable( self ), a, n > 0 );
  CHECKING( const Integer previousOffset = isSeekable(self) ? offset(self): 0;)
  CHECKING( const Integer OLD(offset) = previousOffset; )

  readBytes( self, a, n );

  /* If ok, loop backwards over each byte and copy/expand/sign-extend it: */

  if ( self->data->ok ) {
    const signed char* const c = (const signed char*) a;
    Integer i;

    for ( i = n - 1; i >= 0; --i ) {
      a[ i ] = c[ i ];
    }
  } else {
    *a = 0; /* Zero all 8 bytes. */
  }

  POST4( isReadMode( self ),
         IN_RANGE( a[ 0 ], -128, 127 ),
         IMPLIES_ELSE( ok( self ), IN_RANGE( a[ n - 1 ], -128, 127 ),
                       a[ 0 ] == 0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n,
                                IN_RANGE( offset( self ),
                                          OLD( offset ),
                                          OLD( offset ) + n ) ) ) );
}


/******************************************************************************
PURPOSE: read16BitIntegers - Read an array of MSB 16-bit signed integers and
         expand them into 64-bit integers.
INPUTS:  Stream*  self  The stream to read from.
         Integer  n     The number of values to read [1, INTEGER_MAX / 8].
OUTPUTS: Integer* a     The array of integers read.
         Stream*  self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called and a[ 0 ] is zero.
******************************************************************************/

static void read16BitIntegers( Stream* self, Integer* a, Integer n ) {

  PRE3( isReadable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer previousOffset = isSeekable(self) ? offset(self): 0;)
  CHECKING( const Integer OLD(offset) = previousOffset; )

  const Integer bytesPerValue = 2; /* Two bytes per 16 bit value. */

  readBytes( self, a, n * bytesPerValue );

  /* If ok, loop backwards over each 2-bytes and copy/expand/sign-extend it: */

  if ( self->data->ok ) {
    const unsigned char* c   = (const unsigned char*) a;
    const unsigned char* src = c + n * bytesPerValue;
    Integer*             dst = a + n;

    do {
      const unsigned char lowByte  = *--src;
      const   signed char highByte = *--src;
      *--dst = ( highByte << 8 ) + lowByte;
      CHECK( IN_RANGE( *dst, -32768, 32767 ) );
    } while ( dst != a );

    CHECK2( dst == a, src == c );
  } else {
    *a = 0; /* Zero all 8 bytes. */
  }

  POST4( isReadMode( self ),
         IN_RANGE( a[ 0 ], -32768, 32767 ),
         IMPLIES_ELSE( ok( self ), IN_RANGE( a[ n - 1 ], -32768, 32767 ),
                       a[ 0 ] == 0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 2,
                                IN_RANGE( offset( self ),
                                          OLD( offset ),
                                          OLD( offset ) + n * 2 ) ) ) );
}


/******************************************************************************
PURPOSE: read32BitIntegers - Read an array of MSB 32-bit signed integers and
         expand them into 64-bit integers.
INPUTS:  Stream*  self  The stream to read from.
         Integer  n     The number of values to read [1, INTEGER_MAX / 8].
OUTPUTS: Integer* a     The array of integers read.
         Stream*  self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called and a[ 0 ] is zero.
******************************************************************************/

static void read32BitIntegers( Stream* self, Integer* a, Integer n ) {

  PRE3( isReadable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer previousOffset = isSeekable(self) ? offset(self): 0;)
  CHECKING( const Integer OLD(offset) = previousOffset; )

  const Integer bytesPerValue = 4; /* Four bytes per 32 bit value. */

  readBytes( self, a, n * bytesPerValue );

  /* If ok, loop backwards over each 4-bytes and copy/expand/sign-extend it: */

  if ( self->data->ok ) {
    const unsigned char* c   = (const unsigned char*) a;
    const unsigned char* src = c + n * bytesPerValue;
    Integer*             dst = a + n;

    do {
      const unsigned char byte1 = *--src;
      const unsigned char byte2 = *--src;
      const unsigned char byte3 = *--src;
      const   signed char byte4 = *--src;
      *--dst = ( byte4 << 24 ) + ( byte3 << 16 ) + ( byte2 << 8 ) + byte1;
      CHECK( IN_RANGE( *dst, -2147483647 - 1, 2147483647 ) );
    } while ( dst != a );

    CHECK2( dst == a, src == c );
  } else {
    *a = 0; /* Zero all 8 bytes. */
  }

  POST4( isReadMode( self ),
         IN_RANGE( a[ 0 ], -2147483647 - 1, 2147483647 ),
         IMPLIES_ELSE( ok( self ),
                       IN_RANGE( a[ n - 1 ], -2147483647 - 1, 2147483647 ),
                       a[ 0 ] == 0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 4,
                                IN_RANGE( offset( self ),
                                          OLD( offset ),
                                          OLD( offset ) + n * 4 ) ) ) );
}


/******************************************************************************
PURPOSE: read64BitIntegers - Read an array of MSB 64-bit signed integers.
INPUTS:  Stream*  self  The stream to read from.
         Integer  n     The number of values to read [1, INTEGER_MAX / 8].
OUTPUTS: Integer* a     The array of integers read.
         Stream*  self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called and a[ 0 ] is zero.
******************************************************************************/

static void read64BitIntegers( Stream* self, Integer* a, Integer n ) {

  PRE3( isReadable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer previousOffset = isSeekable(self) ? offset(self): 0;)
  CHECKING( const Integer OLD(offset) = previousOffset; )

  readBytes( self, a, n * 8 ); /* Eight bytes per 64-bit value. */

  if ( ! self->data->ok ) {
    *a = 0; /* Zero all 8 bytes. */
  } else if ( IS_LITTLE_ENDIAN ) {
    swap8Bytes( a, a, n );
  }

  POST4( isReadMode( self ),
         IN_RANGE( a[ 0 ],
                   INTEGER_CONSTANT(-9223372036854775807)-INTEGER_CONSTANT(1),
                   INTEGER_CONSTANT( 9223372036854775807) ),
         IMPLIES_ELSE( ok( self ),
                       IN_RANGE( a[ n - 1 ],
                    INTEGER_CONSTANT(-9223372036854775807)-INTEGER_CONSTANT(1),
                                 INTEGER_CONSTANT( 9223372036854775807) ),
                       a[ 0 ] == 0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 8,
                                IN_RANGE( offset( self ),
                                          OLD( offset ),
                                          OLD( offset ) + n * 8 ) ) ) );
}


/******************************************************************************
PURPOSE: read32BitReals - Read an array of (big-endian) IEEE-754 32-bit
         floating-point values and expand them into 64-bit Reals.
INPUTS:  Stream* self  The stream to read from.
         Integer n     The number of values to read [1, INTEGER_MAX / 8].
OUTPUTS: Real*   a     The array of reals read.
         Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called and a[ 0 ] is zero.
         On IEEE platforms, a[ i ] may be Nan, -Inf or +Inf otherwise
         non-finite values may be returned as 0, -HUGE_VAL, or +HUGE_VAL.
         See man pages on XDR and IEEE floating-point.
         Uses fread() on native IEEE platforms and XDR (xdr_vector/xdr_float)
         elsewhere.
******************************************************************************/

static void read32BitReals( Stream* self, Real* a, Integer n ) {

  PRE3( isReadable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  ensureReadMode( self );

  if ( IMPLIES_ELSE( USES_NATIVE_XDR(float), isSizet(n), isUnsignedInt(n) ) ) {
    data->ok = READ_ITEMS( data, float, n, a );
  } else {
    streamArrayBuffered( self, n, sizeof (float),
                         USES_NATIVE_XDR( float ) ? 0 : (xdrproc_t) xdr_float,
                         a );
  }

  if ( ! data->ok ) {
    *a = 0.0; /* Zero all 8 bytes. */
  } else if ( sizeof (float) != sizeof (Real) || BROKEN_CRAY_XDR ) {

    /* Expand floats to Reals (unless they are equivalent, e.g., on Cray): */

    const float* fa    = (const float*) a;
    const float* src   = fa + n;
    Real*        dst   = a  + n;

    do {
      float value = *--src;
#ifdef _CRAY
      convertXDRRead32BitNonFinite( &value );
#endif
      *--dst = value;
    } while ( dst != a );

    CHECK2( dst == a, src == fa );
  }

  checkAndReport( self, "read", n, "32-bit reals" );

  POST3( isReadMode( self ), IMPLIES( ! ok( self ), a[ 0 ] == 0.0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 4,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + n * 4 ) ) ) );
}


/******************************************************************************
PURPOSE: read64BitReals - Read an array of (big-endian) IEEE-754 64-bit
         floating-point values.
INPUTS:  Stream* self  The stream to read from.
         Integer n     The number of values to read [1, INTEGER_MAX / 8].
OUTPUTS: Real*   a     The array of reals read.
         Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called and a[ 0 ] is zero.
         On IEEE platforms, a[ i ] may be Nan, -Inf or +Inf otherwise
         non-finite values may be returned as 0, -HUGE_VAL, or +HUGE_VAL.
         See man pages on XDR and IEEE floating-point.
         Uses fread() on native IEEE platforms and XDR (xdr_vector/xdr_double)
         elsewhere.
******************************************************************************/

static void read64BitReals( Stream* self, Real* a, Integer n ) {

  PRE3( isReadable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  ensureReadMode( self );


  if ( IMPLIES_ELSE( USES_NATIVE_XDR(double), isSizet(n), isUnsignedInt(n) )) {
    data->ok = READ_ITEMS( data, double, n, a );
  } else {
    streamArrayBuffered( self, n, sizeof (double),
                         USES_NATIVE_XDR(double) ? 0 : (xdrproc_t) xdr_double,
                         a );
  }

#ifdef _CRAY
  if ( data->ok && BROKEN_CRAY_XDR ) {
    Integer i;

    for ( i = 0; i < n; ++i ) {
      convertXDRRead64BitNonFinite( a + i );
    }
  }
#endif

  if ( ! data->ok ) {
    *a = 0.0;
  }

  checkAndReport( self, "read", n, "64-bit reals" );

  POST3( isReadMode( self ), IMPLIES( ! ok( self ), a[ 0 ] == 0.0 ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 8,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + n * 8 ) ) ) );
}


/******************************************************************************
PURPOSE: writeString - Write a string to a stream.
INPUTS:  Stream*     self   The stream to write to.
         const char* format The string to write. May be a printf()-like format
                            string like "%s %c\n" followed by appropriate args.
OUTPUTS: Stream*     self   Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Implemented in terms of stdarg and vfprintf (man vfprintf) so no
         additional checking is done to verify the appropriateness of the
         format specifications with respect to the actual arguments passed.
         Use appropriate portable macro definitions (from BasicNumerics.h)
         for formats for Integer and Real types, for example:
         Integer i = 123456789012345; Real r = 0.123456789012345E+100;
         stream->writeString( stream, "i = %"INTEGER_FORMAT"\n", i );
         stream->writeString( stream, "r = %"REAL_E_FORMAT"\n", r );
******************************************************************************/

static void writeString( Stream* self, const char* format, ... ) {

  PRE2( isWritable( self ), format );

  StreamData* const data = self->data;
  va_list args; /* For stdarg magic. */

  ensureWriteMode( self );
  va_start( args, format );                          /* Begin stdarg magic.  */
  data->ok = vfprintf( data->file, format, args ) != -1; /* Pass args along. */
  va_end( args );                                    /* End of stdarg magic. */
  data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  checkAndReport( self, "write", 1, "string" );

  POST( isWriteMode( self ) );
}


/******************************************************************************
PURPOSE: writeByte - Write a byte to a stream.
INPUTS:  Stream*       self  The stream to write to.
         unsigned char x     The byte to write.
OUTPUTS: Stream*       self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
******************************************************************************/

static void writeByte( Stream* self, unsigned char x ) {

  PRE( isWritable( self ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  ensureWriteMode( self );
  data->ok = fwrite( &x, 1, 1, data->file ) == 1;
  data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  checkAndReport( self, "write", 1, "byte" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 1,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 1 ) ) ) );
}


/******************************************************************************
PURPOSE: write8BitInteger - Write an 8-bit signed integer to a stream.
INPUTS:  Stream* self  The stream to write to.
         Integer x     The integer to write.
OUTPUTS: Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes signed value clamped to range [-128, 127].
******************************************************************************/

static void write8BitInteger( Stream* self, Integer x ) {

  PRE( isWritable( self ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  const signed char c = CLAMPED_TO_RANGE( x, -128, 127 );

  ensureWriteMode( self );
  data->ok = fwrite( &c, 1, 1, data->file ) == 1;
  data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  checkAndReport( self, "write", 1, "8-bit integer" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 1,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 1 ) ) ) );
}


/******************************************************************************
PURPOSE: write16BitInteger - Write a signed MSB 16-bit integer to a stream.
INPUTS:  Stream* self  The stream to write to.
         Integer x     The integer to write.
OUTPUTS: Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes signed value clamped to [-32768, 32767].
******************************************************************************/

static void write16BitInteger( Stream* self, Integer x ) {

  PRE( isWritable( self ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data  = self->data;
  const Integer     value = CLAMPED_TO_RANGE( x, -32768, 32767 );
  unsigned char bytes[ 2 ];

  bytes[ 1 ] =   value        & 0xff; /* Low-byte.  */
  bytes[ 0 ] = ( value >> 8 ) & 0xff; /* High-byte. */

  DEBUG( fprintf( stderr, "value = %"INTEGER_FORMAT"\n"
                  "writing high-byte = %d (%x), low-byte = %d (%x)\n",
                  value, bytes[ 0 ], bytes[ 0 ], bytes[ 1 ], bytes[ 1 ] ); )

  ensureWriteMode( self );
  data->ok = fwrite( bytes, 2, 1, data->file ) == 1;
  data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  checkAndReport( self, "write", 1, "16-bit integer" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 2,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 2 ) ) ) );
}


/******************************************************************************
PURPOSE: write32BitInteger - Write a signed MSB 32-bit integer to a stream.
INPUTS:  Stream* self  The stream to write to.
         Integer x     The integer to write.
OUTPUTS: Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes signed value clamped to [-2147483648, 2147483647].
******************************************************************************/

static void write32BitInteger( Stream* self, Integer x ) {

  PRE( isWritable( self ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data  = self->data;
  const Integer     value = CLAMPED_TO_RANGE( x, -2147483647 - 1, 2147483647 );
  unsigned char bytes[ 4 ];

  bytes[ 3 ] =   value         & 0xff;
  bytes[ 2 ] = ( value >>  8 ) & 0xff;
  bytes[ 1 ] = ( value >> 16 ) & 0xff;
  bytes[ 0 ] = ( value >> 24 ) & 0xff;

  ensureWriteMode( self );
  data->ok = fwrite( bytes, 4, 1, data->file ) == 1;
  data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  checkAndReport( self, "write", 1, "32-bit integer" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 4,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 4 ) ) ) );
}


/******************************************************************************
PURPOSE: write64BitInteger - Write a signed MSB 64-bit integer to a stream.
INPUTS:  Stream* self  The stream to write to.
         Integer x     The integer to write.
OUTPUTS: Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes signed value in the 64-bit representable range
         [-9223372036854775808, -9223372036854775807].
******************************************************************************/

static void write64BitInteger( Stream* self, Integer x ) {

  PRE( isWritable( self ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
  const Integer value    = IS_LITTLE_ENDIAN ? swapped8Bytes( x ) : x;

  ensureWriteMode( self );
  data->ok = fwrite( &value, 8, 1, data->file ) == 1;
  data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  checkAndReport( self, "write", 1, "64-bit integer" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 8,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 8 ) ) ) );
}


/******************************************************************************
PURPOSE: write32BitReal - Write an IEEE-754/XDR 32-bit real to a stream.
INPUTS:  Stream* self  The stream to write to.
         Real    x     The real to write.
OUTPUTS: Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes value resulting from assignment to a float then written as a
         32-bit XDR float. This may yield +/- infinity or HUGE_VAL for values
         that are beyond the range of a 32-bit representation.
         See man pages on XDR and IEEE floating-point.
         Uses fwrite() on native IEEE platforms and XDR (xdr_float) elsewhere.
******************************************************************************/

static void write32BitReal( Stream* self, Real x ) {

  PRE( isWritable( self ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;
#if _CRAY || __osf__
  /* Clamp to avoid floating-point overflow exception. */
  /* const */ float floatX = CLAMPED_TO_RANGE( x, -FLT_MAX, FLT_MAX );
#else
  /* const */ float floatX = x; /* Otherwise allow inf, etc. */
#endif
  /* Note: xdr_float() takes a non-const argument. */

  ensureWriteMode( self );
  data->ok = WRITE_ITEM( data, float, &floatX );
  data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  checkAndReport( self, "write", 1, "32-bit real" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 4,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 4 ) ) ) );
}


/******************************************************************************
PURPOSE: write64BitReal - Write an IEEE-754/XDR 64-bit real to a stream.
INPUTS:  Stream* self  The stream to write to.
         Real    x     The real to write.
OUTPUTS: Stream* self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes value as a 64-bit XDR double.
         This may yield +/- infinity or HUGE_VAL for values that are beyond
         the range of a 64-bit representation.
         See man pages on XDR and IEEE floating-point.
         Uses fwrite() on native IEEE platforms and XDR (xdr_double) elsewhere.
******************************************************************************/

static void write64BitReal( Stream* self, Real x ) {

  PRE( isWritable( self ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  CHECK( IMPLIES( USES_NATIVE_XDR( double ),
                  AND2( sizeof (Real) == 8, sizeof (double) == 8 ) ) );

  ensureWriteMode( self );
  data->ok = WRITE_ITEM( data, double, &x );
  data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  checkAndReport( self, "write", 1, "64-bit real" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + 8,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + 8 ) ) ) );
}


/******************************************************************************
PURPOSE: writeBytes - Write a block of bytes to a stream.
INPUTS:  Stream*     self  The stream to write to.
         const void* a     The block of bytes to write.
         Integer     n     The number of bytes to write.
OUTPUTS: Stream*     self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
******************************************************************************/

static void writeBytes( Stream* self, const void* a, Integer n ) {

  PRE3( isWritable( self ), a, n > 0 );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  ensureWriteMode( self );

  if ( isSizet( n ) && ! BROKEN_FWRITE ) {
    data->ok = fwrite( a, 1, n, data->file ) == n;
    data->ok = AND2( data->ok, fflush( data->file ) == 0 );
  } else {
    streamArrayBuffered( self, n, 1, 0, (void*) a );
  }

  checkAndReport( self, "write", n, "bytes" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + n ) ) ) );
}


/******************************************************************************
PURPOSE: write8BitIntegers - Write an array of signed 8-bit integers.
INPUTS:  Stream*        self  The stream to write to.
         const Integer* a     The integers to write.
         Integer        n     The number of integers to write.
OUTPUTS: Stream*        self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes signed values clamped to the range [-128, 127].
         Allocates/copies/frees a temporary buffer copy of the data and
         performs multiple writes as needed.
******************************************************************************/

static void write8BitIntegers( Stream* self, const Integer* a, Integer n ) {

  PRE3( isWritable( self ), a, n > 0 );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  writeClampedCopy( self, a, n, 1, clampTo8BitInteger, "8-bit integers" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n,
                                IN_RANGE( offset( self ),
                                          OLD( offset ),
                                          OLD( offset ) + n ) ) ) );
}


/******************************************************************************
PURPOSE: write16BitIntegers - Write an array of signed MSB 16-bit integers.
INPUTS:  Stream*        self  The stream to write to.
         const Integer* a     The integers to write.
         Integer        n     The number of values to write in the range:
                              [1, INTEGER_MAX / 8].
OUTPUTS: Stream*        self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes signed values clamped to the range [-32768, 32767].
         Allocates/copies/frees a temporary buffer copy of the data and
         performs multiple writes as needed.
******************************************************************************/

static void write16BitIntegers( Stream* self, const Integer* a, Integer n ) {

  PRE3( isWritable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  writeClampedCopy( self, a, n, 2, clampTo16BitInteger, "16-bit integers" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 2,
                                IN_RANGE( offset( self ),
                                          OLD( offset ),
                                          OLD( offset ) + n * 2 ) ) ) );
}


/******************************************************************************
PURPOSE: write32BitIntegers - Write an array of signed MSB 32-bit integers.
INPUTS:  Stream*        self  The stream to write to.
         const Integer* a     The integers to write.
         Integer        n     The number of values to write in the range:
                              [1, INTEGER_MAX / 8].
OUTPUTS: Stream*        self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes signed values clamped to the range [-2147483648, 2147483647].
         Allocates/copies/frees a temporary buffer copy of the data and
         performs multiple writes as needed.
******************************************************************************/

static void write32BitIntegers( Stream* self, const Integer* a, Integer n ) {

  PRE3( isWritable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  writeClampedCopy( self, a, n, 4, clampTo32BitInteger, "32-bit integers" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 4,
                                IN_RANGE( offset( self ),
                                          OLD( offset ),
                                          OLD( offset ) + n * 4 ) ) ) );
}


/******************************************************************************
PURPOSE: write64BitIntegers - Write an array of MSB 64-bit integers.
INPUTS:  Stream*        self  The stream to write to.
         const Integer* a     The integers to write.
         Integer        n     The number of values to write in the range:
                              [1, INTEGER_MAX / 8].
OUTPUTS: Stream*        self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         On little-endian platforms or if n > MAX_SIZET, this routine
         allocates/copies/frees a temporary buffer copy of the data and
         performs multiple writes as needed.
******************************************************************************/

static void write64BitIntegers( Stream* self, const Integer* a, Integer n ) {

  PRE3( isWritable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  if ( IS_LITTLE_ENDIAN ) {
    writeClampedCopy( self, a, n, 8, clampTo64BitInteger, "64-bit integers" );
  } else {
    StreamData* const data = self->data;

    ensureWriteMode( self );

    if ( isSizet( n ) && ! BROKEN_FWRITE ) {
      data->ok = fwrite( a, 8, n, data->file ) == n;
      data->ok = AND2( data->ok, fflush( data->file ) == 0 );
    } else {
      streamArrayBuffered( self, n, 8, 0, (void*) a );
    }

    checkAndReport( self, "write", n, "64-bit integers" );
  }

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 8,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + n * 8 ) ) ) );
}


/******************************************************************************
PURPOSE: write32BitReals - Write an array of MSB 32-bit IEEE-754/XDR reals.
INPUTS:  Stream*     self  The stream to write to.
         const Real* a     The reals to write.
         Integer     n     The number of values to write in the range:
                           [1, INTEGER_MAX / 8].
OUTPUTS: Stream*     self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes value truncated to 32-bit precision.
         Uses fwrite() on native IEEE platforms and XDR xdr_vector/xdr_float
         elsewhere.
         Allocates/copies/frees a temporary buffer copy of the data
         on platforms where sizeof (float) != 8, for example, on non-Crays.
******************************************************************************/

static void write32BitReals( Stream* self, const Real* a, Integer n ) {

  PRE3( isWritable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  ensureWriteMode( self );
  data->ok = 0;

  if ( sizeof (float) == sizeof (Real) ) { /* Avoid (explicit) copy on Cray: */

    if ( isSizet( n ) ) { /* A single write will work: */
      data->ok = xdr_vector( &data->xdr, (char*) a, n, sizeof (Real),
                             (xdrproc_t) xdr_float);
      data->ok = AND2( data->ok, fflush( data->file ) == 0 );
    } else { /* Must use multiple writes: */
      streamArrayBuffered( self, n, 8, (xdrproc_t) xdr_float, (void*) a );
    }

    checkAndReport( self, "write", n, "32-bit reals" );
  } else {
    /* Must copy data and write the copy: */
    const Integer maximumBufferSize = 1024 * 1024;
    const Integer bufferSize        = MIN( n, maximumBufferSize );
    float* copy                     = 0;

    CHECK( IN_RANGE( maximumBufferSize, 1024, INT_MAX ) );
    CHECK( IN_RANGE( bufferSize, 1, INTEGER_MAX / sizeof (float) ) );
    copy = NEW_ZERO( float, bufferSize );

    if ( copy ) {
      Integer itemsWritten = 0; /* Number of items written so far. */

      do {
        const Integer itemsRemaining  = n - itemsWritten;
        const Integer itemsToWriteNow = MIN( itemsRemaining, bufferSize );
        const Real* src = a + itemsWritten;
        float*      dst = copy;
        Integer count   = itemsToWriteNow;

        CHECK2( itemsWritten < n, IN_RANGE( itemsToWriteNow, 1, bufferSize ) );

        while ( count-- ) {
#if _CRAY || __osf__
          /* Must clamp to avoid floating-point overflow exception. */
          const float value = CLAMPED_TO_RANGE( *src, -FLT_MAX, FLT_MAX );
#else
          const float value = *src; /* Otherwise allow inf, etc. */
#endif
          ++src;
          *dst++ = value; /* Copy to 32-bit representation. */
        }

        if ( USES_NATIVE_XDR( float ) ) /* Use fwrite() if native XDR: */ {
          CHECK2( sizeof (float) == 4, isSizet( itemsToWriteNow ) );
          data->ok = fwrite( copy, sizeof (float), itemsToWriteNow, data->file)
                     == itemsToWriteNow;
        } else {
          /* Must use XDR: */
          CHECK ( isUnsignedInt( itemsToWriteNow ) );
          data->ok = xdr_vector( &data->xdr, (char*) copy, itemsToWriteNow,
                                 sizeof (float), (xdrproc_t) xdr_float );
        }

        data->ok = AND2( data->ok, fflush( data->file ) == 0 );

        if ( data->ok ) {
          itemsWritten += itemsToWriteNow;
        }

      } while ( AND2( data->ok, itemsWritten != n ) );

      FREE_ZERO( copy );
      checkAndReport( self, "write", n, "32-bit reals" );
    }
  }

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 4,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + n * 4 ) ) ) );
}


/******************************************************************************
PURPOSE: write64BitReals - Write an array of MSB 64-bit IEEE-754/XDR reals.
INPUTS:  Stream*     self  The stream to write to.
         const Real* a     The reals to write.
         Integer     n     The number of values to write in the range:
                           [1, INTEGER_MAX / 8].
OUTPUTS: Stream*     self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Uses fwrite() on native IEEE platforms and XDR xdr_vector/xdr_double
         elsewhere.
******************************************************************************/

static void write64BitReals( Stream* self, const Real* a, Integer n ) {

  PRE3( isWritable( self ), a,
        IN_RANGE( n, 1, INTEGER_CONSTANT( 1152921504606846975 ) ) );
  CHECKING( const Integer OLD(offset) = isSeekable(self) ? offset(self) : 0; )

  StreamData* const data = self->data;

  ensureWriteMode( self );

  if ( USES_NATIVE_XDR( double ) ) {

    CHECK2( sizeof (Real) == 8, sizeof (double) == 8 );

    /* If small enough, write data in one call otherwise buffer it: */

    if ( isSizet( n ) && ! BROKEN_FWRITE ) {
      data->ok = fwrite( a, 8, n, data->file ) == n;
      data->ok = AND2( data->ok, fflush( data->file ) == 0 );
    } else {
      streamArrayBuffered( self, n, sizeof (Real), 0, (void*) a );
    }

  } else { /* Must use XDR: */

    /* If small enough, write data in one call otherwise buffer it: */

    if ( isUnsignedInt( n ) ) {
      data->ok = xdr_vector( &data->xdr, (char*) a, n, sizeof (Real),
                             (xdrproc_t) xdr_double );
      data->ok = AND2( data->ok, fflush( data->file ) == 0 );
    } else {
      streamArrayBuffered( self, n, sizeof (Real), (xdrproc_t) xdr_double,
                           (void*) a );
    }
  }

  checkAndReport( self, "write", n, "64-bit reals" );

  POST2( isWriteMode( self ),
         IMPLIES( isSeekable( self ),
                  IMPLIES_ELSE( ok( self ),
                                offset( self ) == OLD( offset ) + n * 8,
                                IN_RANGE( offset( self ),
                                          0,
                                          OLD( offset ) + n * 8 ) ) ) );
}


/******************************************************************************
PURPOSE: invariant - Checks if a Stream structure is valid.
INPUTS:  const Stream* self  The stream to check.
RETURNS: Integer 1 if non-null and valid, else 0.
NOTES:   If this routine ever returns a value other than 1, then this indicates
         that the code contains a defect.
******************************************************************************/

static Integer invariant( const Stream* self ) {

  const StreamData* const data = self ? self->data : 0;
  const Integer result =
    AND18( self,
           AND5(  self->free               == closeStream,
                  self->flush              == flushStream,
                  self->seekFromStart      == seekFromStart,
                  self->seekFromEnd        == seekFromEnd,
                  self->seekFromCurrent    == seekFromCurrent ),
           AND17( self->readString         == readString,
                  self->readWord           == readWord,
                  self->readByte           == readByte,
                  self->read8BitInteger    == read8BitInteger,
                  self->read16BitInteger   == read16BitInteger,
                  self->read32BitInteger   == read32BitInteger,
                  self->read64BitInteger   == read64BitInteger,
                  self->read32BitReal      == read32BitReal,
                  self->read64BitReal      == read64BitReal,
                  self->readBytes          == readBytes,
                  self->readUpToNBytes     == readUpToNBytes,
                  self->read8BitIntegers   == read8BitIntegers,
                  self->read16BitIntegers  == read16BitIntegers,
                  self->read32BitIntegers  == read32BitIntegers,
                  self->read64BitIntegers  == read64BitIntegers,
                  self->read32BitReals     == read32BitReals,
                  self->read64BitReals     == read64BitReals ),
           AND15( self->writeString        == writeString,
                  self->writeByte          == writeByte,
                  self->write8BitInteger   == write8BitInteger,
                  self->write16BitInteger  == write16BitInteger,
                  self->write32BitInteger  == write32BitInteger,
                  self->write64BitInteger  == write64BitInteger,
                  self->write32BitReal     == write32BitReal,
                  self->write64BitReal     == write64BitReal,
                  self->writeBytes         == writeBytes,
                  self->write8BitIntegers  == write8BitIntegers,
                  self->write16BitIntegers == write16BitIntegers,
                  self->write32BitIntegers == write32BitIntegers,
                  self->write64BitIntegers == write64BitIntegers,
                  self->write32BitReals    == write32BitReals,
                  self->write64BitReals    == write64BitReals ),
           AND12( self->invariant          == invariant,
                  self->ok                 == ok,
                  self->isReadable         == isReadable,
                  self->isWritable         == isWritable,
                  self->isSeekable         == isSeekable,
                  self->isBlocking         == isBlocking,
                  self->isAtEnd            == isAtEnd,
                  self->offset             == offset,
                  self->size               == size,
                  self->name               == name,
                  self->file               == file,
                  self->descriptor         == descriptor ),
           self->data,
           data == self->data,
           data->name,
           data->name[ 0 ],
           data->mode,
           data->mode[ 0 ],
           matchesWord( data->mode,
                        " r rb w wb a ab r+ r+b rb+ w+ w+b wb+ a+ a+b ab+ " ),
           IN4( data->type, FILE_STREAM, PIPE_STREAM, SOCKET_STREAM ),
           IMPLIES_ELSE( data->type == SOCKET_STREAM,
                         AND2( data->port > 0, strcmp(data->mode, "r+b") == 0),
                         data->port == 0 ),
           IMPLIES( strcmp( data->name, "-stdin" ) == 0,
                    matchesWord( data->mode, " r rb " ) ),
           IMPLIES( matchesWord( data->name, " -stdout -stderr " ),
                    data->mode[ 0 ] != 'r' ),
           IN3( data->xdrMode, XDR_ENCODE, XDR_DECODE ),
           data->file );

  CHECK( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: ok - Determine if the last command succeeded.
INPUTS:  const Stream* self  The stream to examine.
RETURNS: Integer 1 if the last command was successful, else 0.
NOTES:   This status is set to 1 at the beginning of every command.
******************************************************************************/

static Integer ok( const Stream* self ) {

  PRE( self );

  const Integer result = self->data->ok;

  POST( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: isReadable - Checks if a stream's mode (established at creation)
         permits reading.
INPUTS:  const Stream* self  The stream to check.
RETURNS: Integer 1 if true, else 0.
******************************************************************************/

static Integer isReadable( const Stream* self ) {

  PRE( self );

  const Integer result = self->data->mode[ 0 ] == 'r';

  POST( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: isWritable - Checks if a stream's mode (established at creation)
         permits writing.
INPUTS:  const Stream* self  The stream to check.
RETURNS: Integer 1 if true, else 0.
******************************************************************************/

static Integer isWritable( const Stream* self ) {

  PRE( self );

  const char* const mode = self->data->mode;

  const Integer result = OR3( IN3( mode[ 0 ], 'w', 'a' ),
                              mode[ 1 ] == '+',
                              AND2( mode[ 1 ] != '\0', mode[ 2 ] == '+' ) );

  POST( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: isSeekable - Checks if a stream's type (established at creation)
         permits seeking.
INPUTS:  const Stream* self  The stream to check.
RETURNS: Integer 1 if true, else 0.
******************************************************************************/

static Integer isSeekable( const Stream* self ) {

  PRE( self );

  const StreamData* const data = self->data;
  const Integer result = AND2( data->type == FILE_STREAM,
                               ! matchesWord( data->name,
                                     " -stdin -stdout -stderr /dev/null " ) );

  POST( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: isBlocking - Checks if a stream may block if read from (established
         at creation based upon type).
INPUTS:  Stream* self  The stream to check.
RETURNS: Integer 1 if true, else 0.
******************************************************************************/

static Integer isBlocking( const Stream* self ) {

  PRE( isReadable( self ) );

  const StreamData* const data = self->data;
  const Integer result = OR2( data->type != FILE_STREAM,
                              strcmp( data->name, "-stdin" ) == 0 );

  POST( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: isAtEnd - Checks if a stream's offset is at (or beyond) the end.
INPUTS:  const Stream* self  The stream to check.
RETURNS: Integer 1 if EOF, else 0.
NOTES:   May block if used on pipes or sockets (with empty buffers).
         feof() is not used since it does not work with fseek().
         ftell() (compared with a cached file size) is five times slower
         than the fgetc(); ungetc() approach used here! Consequently,
         const-cast-away is used therefore this routine is not thread-safe.
******************************************************************************/

static Integer isAtEnd( const Stream* self ) {

  PRE( isReadable( self ) );

  FILE* const file = (FILE*) self->data->file; /* HACK: const-cast-away. */
  const int ch = fgetc( file );
  const Integer result = ch == EOF;

  ungetc( ch, file );

  POST( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: offset - Get the current byte offset from the beginning of the file.
INPUTS:  const Stream* self   The stream to check.
RETURNS: Integer the current byte offset from the beginning of the file.
NOTES:   Const-cast-away is used since ftell() is non-const (why?) therefore
         this routine is possibly not thread-safe.
******************************************************************************/

static Integer offset( const Stream* self ) {

  PRE( isSeekable( self ) );

  /* HACK: const-cast-away. Why isn't ftell() const? */

  Integer result = FTELL( (FILE*) self->data->file );

  if ( result < 0 ) {
    result = 0;
  }

  POST( result >= 0 );

  return result;
}


/******************************************************************************
PURPOSE: size - Get the size (in bytes) of an existing file.
INPUTS:  const Stream* self  The stream to check.
RETURNS: Integer             Size, in bytes, of the file or 0 if indeterminate.
******************************************************************************/

static Integer size( const Stream* self ) {

  PRE( isSeekable( self ) );

  Integer result = 0;
  struct stat buffer;

  if ( fstat( descriptor( self ), &buffer ) == 0 ) {
    result = buffer.st_size;
  }

  POST( result >= 0 );

  return result;
}


/******************************************************************************
PURPOSE: name - Get the name of a stream.
INPUTS:  const Stream* self  The stream to check.
RETURNS: const char*         Name of the stream.
NOTES:   HACK! This violates encapsulation but is needed to interface with
         external software. The returned pointer should be treated as read-only
         - i.e., not modified or deallocated.
******************************************************************************/

static const char* name( const Stream* self ) {

  PRE( self );

  const char* const result = self->data->name;

  POST2( result, *result );

  return result;
}


/******************************************************************************
PURPOSE: file - Get the underlying FILE pointer.
INPUTS:  const Stream* self  The stream to violate.
RETURNS: FILE*               Underlying FILE pointer.
NOTES:   HACK! This violates encapsulation but is needed to interface with
         external software. This file pointer should not be closed.
******************************************************************************/

static FILE* file( const Stream* self ) {

  PRE( self );

  FILE* const result = self->data->file;

  POST( result );

  return result;
}


/******************************************************************************
PURPOSE: descriptor - Get the underlying file descriptor.
INPUTS:  const Stream* self  The stream to violate.
RETURNS: Integer             Underlying file descriptor number.
NOTES:   HACK! This violates encapsulation but is needed to interface with
         external software. This file descriptor should not be closed.
******************************************************************************/

static Integer descriptor( const Stream* self ) {

  PRE( self );

  const Integer result = fileno( self->data->file );

  POST( result >= 0 );

  return result;
}


/*============================ PRIVATE FUNCTIONS ============================*/


/******************************************************************************
PURPOSE: deallocate - Deallocates (and zeros-out) a stream structure.
INPUTS:  Stream* self The partially or fully allocated stream structure.
OUTPUTS: Stream* self The deallocated and zero'd out stream structure.
NOTE:    Must coordinate the implementation of this routine with allocate().
******************************************************************************/

static void deallocate( Stream* self ) {
  PRE0( self );

  const Integer charsInSegment = sizeof (Stream)     / sizeof (char) +
                                 sizeof (StreamData) / sizeof (char) +
                                 strlen( self->data->name ) + 1 +
                                 strlen( self->data->mode ) + 1;

  memset( self, 0, charsInSegment * sizeof (char) );
  FREE_ZERO( self );

  POST0( self == 0 );
}





/******************************************************************************
PURPOSE: newSocketStream - Create a socket for reading and writing.
INPUTS:  Integer     port  The port number to use. (E.g., 4321).
         const char* host  The server process passes 0 for host and
                           the client process passes "localhost" (or other
                           host name, e.g., "fred").
RETURNS: Stream*           Initialized file structure, or 0 if failed.
NOTES:   If unsuccessful then failureMessage() is called and 0 is returned.
         This routine blocks awaiting a (single) connection to another
         process (using the same port number). Launch the server first then
         launch the client - perhaps even with a 2 second delay (sleep(3)) -
         to allow the server time to establish the socket on the port before
         the client attempts to connect to it (otherwise the client's attempt
         to connect will fail).
         Multi-client processing requires separate socket calls (with unique
         ports) for each client.

         Example:

           server:

             Stream* theSocket = newSocketStream( 4321, 0 );

             if ( theSocket ) {
               readFromClient( theSocket );
             }

           client:

             sleep( 3 );
             theSocket = newSocketStream( 4321, "localhost" );

             if ( theSocket ) {
               writeToServer( theSocket );
             }

******************************************************************************/

static Stream* newSocketStream( Integer port, const char* host ) {
  PRE02( IN_RANGE( port, 1, 65535 ), IMPLIES( host, *host ) );

  const char* const nameFormat   = "-socket(%s:%"INTEGER_FORMAT")";
  const Integer digitsInInteger  = 20; /* Decimal digits per 64-bit integer. */
  const Integer nameFormatLength = strlen( nameFormat ) + digitsInInteger;
  const Integer hostLength       = host ? strlen( host ) : 0;
  const Integer nameLength       = nameFormatLength + hostLength;
  const char* const mode         = "r+b"; /* Sockets support both read/write.*/

  Stream* result = allocate( nameLength, strlen( mode ) );

  if ( result ) {

    StreamData* const data = result->data;
    const Integer callbackEnabled = failureCallingEnabled();

    failureDisableCalling(); /* Disable possible non-returning client handler*/
    data->file = createSocket( port, host );

    if ( callbackEnabled ) {
      failureEnableCalling(); /* Re-enable for failure. */
    }

    if ( data->file ) {
      sprintf( data->name, nameFormat, host ? host : "", port );
      strcpy( data->mode, mode );
      data->ok      = 1;
      data->type    = SOCKET_STREAM;
      data->port    = port;
      data->xdrMode = XDR_DECODE;
      xdrstdio_create( &data->xdr, data->file, data->xdrMode );
    } else {
      deallocate( result );
      result = 0;
      failureMessage( "Can't create and open socket on"
                      " port %"INTEGER_FORMAT"%s%s.",
                      port, host ? " to host " : "", host ? host : "" );
    }
  }

/*
 * Must bypass invariant in post-condition too since result may be 0.
 * However, in cases where result is not 0, the invariant will be evaluated
 * during calls to ok() anyway.
 */

  POST0( IMPLIES( result, AND5( ok( result ),
                                isReadable( result ),
                                isWritable( result ),
                                isBlocking( result ),
                                ! isSeekable( result ) ) ) );
  return result;
}


/******************************************************************************
PURPOSE: typeName - Get the name of the type of stream.
INPUTS:  const Stream* self  The stream to check.
RETURNS: const char*         Name of stream type.
******************************************************************************/

static const char* typeName( const Stream* self ) {

  PRE( self );

  const Integer type = self->data->type;
  const char* const result =   type == FILE_STREAM   ? "file"
                             : type == PIPE_STREAM   ? "pipe"
                             : type == SOCKET_STREAM ? "socket"
                             : 0;

  POST3( result, *result, matchesWord( result, " file pipe socket " ) );

  return result;
}


/******************************************************************************
PURPOSE: isReadMode - Determine if the stream is in read mode.
INPUTS:  const Stream* self  The stream to check.
RETURNS: Integer 1 if xdrMode is XDR_DECODE.
******************************************************************************/

CHECKING(
static Integer isReadMode( const Stream* self ) {

  PRE( isReadable( self ) );

  const Integer result = self->data->xdrMode == XDR_DECODE;

  POST( IS_BOOL( result ) );

  return result;
}
)


/******************************************************************************
PURPOSE: isWriteMode - Determine if the stream is in write mode.
INPUTS:  const Stream* self  The stream to check.
RETURNS: Integer 1 if xdrMode is XDR_ENCODE.
******************************************************************************/

CHECKING(
static Integer isWriteMode( const Stream* self ) {

  PRE( isWritable( self ) );

  const Integer result = self->data->xdrMode == XDR_ENCODE;

  POST( IS_BOOL( result ) );

  return result;
}
)


/******************************************************************************
PURPOSE: ensureReadMode - Make sure the stream is in read mode.
INPUTS:  Stream* self  The stream to put in read mode (if not already).
OUTPUTS: Stream* self  The stream put in read mode.
******************************************************************************/

static void ensureReadMode( Stream* self ) {

  PRE( isReadable( self ) );

  StreamData* const data = self->data;

  if ( data->xdrMode != XDR_DECODE ) {
    fflush( data->file );      /* Flush the output buffer. */
    xdr_destroy( &data->xdr ); /* Free XDR buffer. */
    data->xdrMode = XDR_DECODE;
    xdrstdio_create( &data->xdr, data->file, data->xdrMode );
  }

  POST( isReadMode( self ) );
}


/******************************************************************************
PURPOSE: ensureWriteMode - Make sure the stream is in write mode.
INPUTS:  Stream* self  The stream to put in write mode (if not already).
OUTPUTS: Stream* self  The stream put in write mode.
******************************************************************************/

static void ensureWriteMode( Stream* self ) {

  PRE( isWritable( self ) );

  StreamData* const data = self->data;

  if ( data->xdrMode != XDR_ENCODE ) {
    fflush( data->file );      /* Flush the input buffer. */
    xdr_destroy( &data->xdr ); /* Free XDR buffer. */
    data->xdrMode = XDR_ENCODE;
    xdrstdio_create( &data->xdr, data->file, data->xdrMode );
  }

  POST( isWriteMode( self ) );
}


/******************************************************************************
PURPOSE: seekStream - Seek to a byte offset in a file.
INPUTS:  Stream* self        The stream to seek on.
         Integer byteOffset  The byte offset to seek by.
         Integer whence      SEEK_SET, SEEK_CUR, SEEK_END.
OUTPUTS: Stream* self        Updated stream structure.
NOTES:   If unsuccessful then failureMessage() is called.
******************************************************************************/

static void seekStream( Stream* self, Integer byteOffset, Integer whence ) {

  PRE3( isSeekable( self ), IN4( whence, SEEK_SET, SEEK_CUR, SEEK_END ),
        IMPLIES( whence == SEEK_SET, byteOffset >= 0 ) );

  StreamData* const data = self->data;

  data->ok = OR2( FSEEK_FTELL_ARE_64_BITS, isSignedLong( byteOffset ) );

  if ( data->ok ) {
    const Integer oldOffset = offset( self ); /*Remember current valid offset*/
    data->ok = FSEEK( data->file, byteOffset, whence ) == 0;

    if ( ! data->ok ) {
      FSEEK( data->file, oldOffset, SEEK_SET ); /* Try to reset position. */
    }
  }

  if ( ! data->ok ) {
    failureMessage( "Can't seek to byte %"INTEGER_FORMAT"%s in file '%s'.",
                    byteOffset,
                    whence == SEEK_CUR ? " from current location" :
                    whence == SEEK_END ? " from end" : "",
                    data->name );
  }

  POST( isSeekable( self ) );
}


/******************************************************************************
PURPOSE: streamArrayBuffered - Read or write an array of data to a stream using
         multiple reads/writes to avoid overflow in cases where more than
         UINT_MAX or SIZET_MAX bytes are to be read/written.
INPUTS:  Stream*         self        Stream to read/write from/to.
         Integer         count       Number of items to read/write.
         Integer         sizeOfItem  Size in bytes of a basic type item.
                                     E.g., sizeof (float)
         xdrproc_t       elproc      Optional: XDR element read/write procedure
                                     e.g., xdr_float. If 0 then fwrite() is
                                     used instead of xdr_vector().
         const void*     array       Array of items to write (if writing mode).
OUTPUTS: void*           array       Array of items read (if reading mode).
         Stream*         self        Updated offset in file written to.
NOTES:   Mode is determined by self->data->xdrMode.
         If unsuccessful then self->data->ok is 0.
******************************************************************************/

static void streamArrayBuffered( Stream* self,
                                 Integer count,
                                 Integer sizeOfItem,
                                 xdrproc_t elproc,
                                 /* const */ void* array ) {

  PRE4( self, IN5( sizeOfItem, 1, 2, 4, 8 ),
        IN_RANGE( count, 1, INTEGER_MAX / sizeOfItem ), array );

  /* Use the largest buffer size allowed by the underlying implementation -
   * either xdr_vector() (u_int) or fwrite() (size_t) provided it is also
   * representable without loss as an Integer, i.e., signed 64-bits: 2^63.
   */

  StreamData* const data     = self->data;
  const Integer readMode     = self->data->xdrMode == XDR_DECODE;
  const Integer sizeOfBuffer =
    elproc ? ( sizeof (unsigned int) < 8 ? UINT_MAX  : INTEGER_MAX )
    :        ( sizeof (size_t)       < 8 ? SIZET_MAX
               : ( ( ! readMode && BROKEN_FWRITE ) ? INT_MAX : INTEGER_MAX ) );

  const Integer maxItemsStreamable = sizeOfBuffer / sizeOfItem;
  const Integer stride             = sizeOfItem   / sizeof (char);
  char* const   charArray          = (char*) array; /* For pointer arith. */
  Integer       itemsStreamed      = 0; /* # of items streamed so far. */

  CHECK3( sizeOfBuffer >= sizeOfItem,
          IN_RANGE( stride, 1, sizeOfItem ),
          IN_RANGE( maxItemsStreamable, 1, sizeOfBuffer ) );

  do {
    const Integer itemsRemaining   = count - itemsStreamed;
    const Integer itemsToStreamNow = MIN( itemsRemaining, maxItemsStreamable );
    const Integer offset           = itemsStreamed * stride; /* To next item.*/
    char* const   address          = charArray + offset; /* Of next access. */

    CHECK2( IN_RANGE( itemsStreamed, 0, count - 1 ), itemsToStreamNow > 0 );

    if ( elproc ) {
      CHECK( isUnsignedInt( itemsToStreamNow ) );
      data->ok = xdr_vector( &data->xdr, address, itemsToStreamNow,
                             sizeOfItem, elproc );
    } else {
      const Integer itemsActuallyStreamed =
        readMode ? fread(  address, sizeOfItem, itemsToStreamNow, data->file )
                 : fwrite( address, sizeOfItem, itemsToStreamNow, data->file );

      CHECK( isSizet( itemsToStreamNow ) ); /* Better late than never. */
      data->ok = itemsActuallyStreamed == itemsToStreamNow;
    }

    if ( ! readMode ) {
      data->ok = AND2( data->ok, fflush( data->file ) == 0 );
    }

    if ( data->ok ) {
      itemsStreamed += itemsToStreamNow;
    }

  } while ( AND2( data->ok, itemsStreamed != count ) );

  CHECK2( IN_RANGE( itemsStreamed, 0, count ),
          IMPLIES( data->ok, itemsStreamed == count ) );
}


/******************************************************************************
PURPOSE: writeClampedCopy - Write an array of Integers using a temporary small
         buffer copy and clamping them to 8, 16, 32 or 64-bit ranges.
INPUTS:  Stream*        self          The stream to write to.
         const Integer* a             The integers to write.
         Integer        n             The number of integers to write.
         Integer        bytesPerItem  The size (in bytes) of each clamped item.
         Clamper        clamper       A routine that copies clamped values
                                      into a temporary buffer.
         const char*    kindOfValues  E.g., "8-bit integers".
OUTPUTS: Stream*        self  Updated file structure.
NOTES:   If unsuccessful then failureMessage() is called.
         Writes value truncated to a subrange as implemented by the given
         called clamper routine.
         Allocates/copies/frees a temporary buffer copy of the data and
         performs multiple writes as needed.
         For 64-bit data, clamping is not actually needed, but this routine is
         still needed to handle byte-swapping on little-endian platforms.
******************************************************************************/

static void writeClampedCopy( Stream* self, const Integer* a, Integer n,
                              Integer bytesPerItem, Clamper clamper,
                              const char* kindOfValues ) {

  PRE7( isWritable( self ), a, IN5( bytesPerItem, 1, 2, 4, 8 ),
        IN_RANGE( n, 1, INTEGER_MAX / bytesPerItem ),
        clamper, kindOfValues, *kindOfValues );
  StreamData* const data = self->data;
  const Integer  maximumBufferSize = 1024 * 1024;
  const Integer  bufferSize        = MIN( n, maximumBufferSize );
  unsigned char* copy              = 0;

  ensureWriteMode( self );
  data->ok = 0;
  CHECK( IN_RANGE( maximumBufferSize, 1024, INT_MAX ) );
  CHECK( IN_RANGE( bufferSize, 1, INTEGER_MAX / bytesPerItem ) );
  copy = NEW_ZERO( unsigned char, bytesPerItem * bufferSize );

  if ( copy ) {

    Integer itemsWritten = 0; /* Number of items written so far. */

    do {
      const Integer itemsRemaining  = n - itemsWritten;
      const Integer itemsToWriteNow = MIN( itemsRemaining, bufferSize );

      CHECK2( IN_RANGE( itemsWritten, 0, n - 1 ),
              IN_RANGE( itemsToWriteNow, 1, bufferSize ) );

      clamper( a + itemsWritten, itemsToWriteNow, copy );
      CHECK( isSizet( itemsToWriteNow ) );

      data->ok = fwrite( copy, bytesPerItem, itemsToWriteNow, data->file )
                 == itemsToWriteNow;
      data->ok = AND2( data->ok, fflush( data->file ) == 0 );

      if ( data->ok ) {
        itemsWritten += itemsToWriteNow;
      }

    } while ( AND2( data->ok, itemsWritten != n ) );

    CHECK2( IN_RANGE( itemsWritten, 0, n ),
            IMPLIES( data->ok, itemsWritten == n ) );

    FREE_ZERO( copy );
    checkAndReport( self, "write", n, kindOfValues );
  }

  POST( isWriteMode( self ) );
}


/******************************************************************************
PURPOSE: checkAndReport - Check if last command failed and if so call failure.
INPUTS:  const Stream* self         The file being read/written.
         const char*   readOrWrite  The kind of operation being performed:
                                    "read" or "write".
         Integer       count        The number of items read/written.
         const char*   dataType     The kind of data being read/written.
                                    e.g., "16-bit integers".
******************************************************************************/

static void checkAndReport( const Stream* self, const char* readOrWrite,
                            Integer count, const char* dataType ) {

  PRE6( self, readOrWrite, matchesWord( readOrWrite, " read write " ),
        count > 0, dataType, *dataType );

  const StreamData* const data = self->data;

  if ( ! data->ok ) {
    failureMessage( "Can't %s %"INTEGER_FORMAT" %s %s %s '%s'.",
                    readOrWrite, count, dataType,
                    *readOrWrite == 'r' ? "from" : "to",
                    typeName( self ), data->name );
  }
}


/******************************************************************************
PURPOSE: allocate - Allocates and zeros a Stream structure.
INPUTS:  Integer nameLength  Number of characters for name (excluding NULL).
         Integer modeLength  Number of characters for mode (excluding NULL).
RETURNS: Stream* result      The allocated, zeroed stream if successful else
                             0 and failure is called (by allocator).
NOTE:    Must coordinate the implementation of this routine with deallocate().
******************************************************************************/

static Stream* allocate( Integer nameLength, Integer modeLength ) {

  PRE02( nameLength > 0, modeLength > 0 );

  /*
   * Allocate just one segment of memory and establish sub-pointers within it:
   *
   *  result =             [ sizeof (Stream)     ]
   *  result->data =       [ sizeof (StreamData) ]
   *  result->data->name = [ nameLength + 1      ]
   *  result->data->mode = [ modeLength + 1      ]
   */

  const Integer charsInSegment = sizeof (Stream)     / sizeof (char) +
                                 sizeof (StreamData) / sizeof (char) +
                                 nameLength + 1 + modeLength + 1;

  char* const segment = NEW_ZERO( char, charsInSegment );
  Stream*     result  = (Stream*) segment;

  if ( result ) {

    result->data       = (StreamData*) ( result       + 1 );
    result->data->name = (char*)       ( result->data + 1 );
    result->data->mode = (char*)       ( result->data->name + nameLength + 1 );

    result->free               = closeStream;
    result->flush              = flushStream;
    result->seekFromStart      = seekFromStart;
    result->seekFromEnd        = seekFromEnd;
    result->seekFromCurrent    = seekFromCurrent;
    result->readString         = readString;
    result->readWord           = readWord;
    result->readByte           = readByte;
    result->read8BitInteger    = read8BitInteger;
    result->read16BitInteger   = read16BitInteger;
    result->read32BitInteger   = read32BitInteger;
    result->read64BitInteger   = read64BitInteger;
    result->read32BitReal      = read32BitReal;
    result->read64BitReal      = read64BitReal;
    result->readBytes          = readBytes;
    result->readUpToNBytes     = readUpToNBytes;
    result->read8BitIntegers   = read8BitIntegers;
    result->read16BitIntegers  = read16BitIntegers;
    result->read32BitIntegers  = read32BitIntegers;
    result->read64BitIntegers  = read64BitIntegers;
    result->read32BitReals     = read32BitReals;
    result->read64BitReals     = read64BitReals;
    result->writeString        = writeString;
    result->writeByte          = writeByte;
    result->write8BitInteger   = write8BitInteger;
    result->write16BitInteger  = write16BitInteger;
    result->write32BitInteger  = write32BitInteger;
    result->write64BitInteger  = write64BitInteger;
    result->write32BitReal     = write32BitReal;
    result->write64BitReal     = write64BitReal;
    result->writeBytes         = writeBytes;
    result->write8BitIntegers  = write8BitIntegers;
    result->write16BitIntegers = write16BitIntegers;
    result->write32BitIntegers = write32BitIntegers;
    result->write64BitIntegers = write64BitIntegers;
    result->write32BitReals    = write32BitReals;
    result->write64BitReals    = write64BitReals;
    result->invariant          = invariant;
    result->ok                 = ok;
    result->isReadable         = isReadable;
    result->isWritable         = isWritable;
    result->isSeekable         = isSeekable;
    result->isBlocking         = isBlocking;
    result->isAtEnd            = isAtEnd;
    result->offset             = offset;
    result->size               = size;
    result->name               = name;
    result->file               = file;
    result->descriptor         = descriptor;
  }

  POST0( IMPLIES( result, AND6( result->data,
                                result->data->name,
                                result->data->mode,
                                result->data->name[ 0 ] == '\0',
                                result->data->mode[ 0 ] == '\0',
                                result->descriptor == descriptor ) ) );

  return result;
}


/******************************************************************************
PURPOSE: clampTo8BitInteger - Copy each value of src to dst clamped to
         signed 8-bit integers in the range [-128, 127].
INPUTS:  const Integer* src   The integers to copy.
         Integer        count The number of integers (src[count], dst[count]).
OUTPUTS: unsigned char* dst   The clamped copy of src integers.
******************************************************************************/

static void clampTo8BitInteger( const Integer* src, Integer count,
                                unsigned char* dst ) {

  PRE04( src, count > 0, dst, (const void*) src != (const void*) dst );

  signed char* const copy = (signed char*) dst;
  Integer index;

  /* Copy clamped to signed 8-bit representation: */

  for ( index = 0; index < count; ++index ) {
    const Integer value        = src[ index ];
    const Integer clampedValue = CLAMPED_TO_RANGE( value, -128, 127 );
    copy[ index ]              = clampedValue;
    CHECK( IN_RANGE( (int) copy[ index ], -128, 127 ) );
  }
}


/******************************************************************************
PURPOSE: clampTo16BitInteger - Copy each value of src to dst clamped to
         signed 16-bit integers in the range [-32768, 32767].
INPUTS:  const Integer* src   The integers to copy.
         Integer        count The number of integers (src[count], dst[count]).
OUTPUTS: unsigned char* dst   The clamped copy of src integers.
******************************************************************************/

static void clampTo16BitInteger( const Integer* src, Integer count,
                                 unsigned char* dst ) {

  PRE04( src, IN_RANGE( count, 1, INTEGER_MAX / 8 ),
         dst, (const void*) src != (const void*) dst );

  unsigned char* copy = dst;
  Integer index;

  /* Copy clamped to signed 16-bit representation: */

  for ( index = 0; index < count; ++index ) {
    const Integer value          = src[ index ];
    const Integer clampedValue   = CLAMPED_TO_RANGE( value, -32768, 32767 );
    const   signed char highByte = ( clampedValue >> 8 ) & 0xff; /* w/  sign */
    const unsigned char lowByte  =   clampedValue        & 0xff; /* w/o sign */
    *copy++ = highByte;
    *copy++ = lowByte;
    DEBUG(fprintf(stderr, "clampedValue = %"INTEGER_FORMAT"\n", clampedValue);)
    DEBUG(fprintf(stderr, "highByte = %d, lowByte = %d\n", highByte, lowByte);)
  }
}


/******************************************************************************
PURPOSE: clampTo32BitInteger - Copy each value of src to dst clamped to
         signed 32-bit integers in the range [-2147483648, 2147483647].
INPUTS:  const Integer* src   The integers to copy.
         Integer        count The number of integers (src[count], dst[count]).
OUTPUTS: unsigned char* dst   The clamped copy of src integers.
******************************************************************************/

static void clampTo32BitInteger( const Integer* src, Integer count,
                                 unsigned char* dst ) {

  PRE04( src, IN_RANGE( count, 1, INTEGER_MAX / 8 ),
         dst, (const void*) src != (const void*) dst );

  unsigned char* copy = dst;
  Integer index;

  /* Copy clamped to signed 32-bit representation: */

  for ( index = 0; index < count; ++index ) {
    const Integer value        = src[ index ];
    const Integer clampedValue =
      CLAMPED_TO_RANGE( value, -2147483647 - 1, 2147483647 );

    CHECK( IN_RANGE( clampedValue, -2147483647 - 1, 2147483647 ) );
    *copy++ = ( clampedValue >> 24 ) & 0xff; /* Highest byte. */
    *copy++ = ( clampedValue >> 16 ) & 0xff;
    *copy++ = ( clampedValue >>  8 ) & 0xff;
    *copy++ =   clampedValue         & 0xff; /* Lowest byte. */
  }
}


/******************************************************************************
PURPOSE: clampTo64BitInteger - Copy byte-swapped values of src to dst.
INPUTS:  const Integer* src   The integers to copy.
         Integer        count The number of integers (src[count], dst[count]).
OUTPUTS: unsigned char* dst   The byte-swapped copy of src integers.
NOTE:    Unlike other clamp routines, this one's purpose is actually only to
         perform byte-swapping (since clamping is unnecessary). Therefore it
         need only be called on little-endian platforms.
******************************************************************************/

static void clampTo64BitInteger( const Integer* src, Integer count,
                                 unsigned char* dst ) {

  PRE05( src, IN_RANGE( count, 1, INTEGER_MAX / 8 ),
         dst, (const void*) src != (const void*) dst, IS_LITTLE_ENDIAN );

  swap8Bytes( (Integer*) dst, src, count );
}




/******************************************************************************
PURPOSE: createSocket - Creates an opened server or client socket.
INPUTS:  Integer     port        Number of the port to open the socket on.
         const char* host        Name of host if client, else 0 if server.
RETURNS: FILE*       result      Pointer to opened socket or 0 if unsuccessful.
NOTES:   If unsuccessful, failureMessage() is called.
         Returned pointer can be used with fread() and fwrite() and when no
         longer needed, it should be closed with fclose().
******************************************************************************/

static FILE* createSocket( Integer port, const char* host ) {

  PRE02( IN_RANGE( port, 1, 65535 ), IMPLIES( host, *host ) );

  FILE* result       = 0;
  int   socketNumber = -1;

  if ( host == 0 ) {
    socketNumber = createServerSocket( port );
  } else {
    socketNumber = createClientSocket( port, host );
  }

  if ( socketNumber != -1 ) {
    const Integer calling = failureCallingEnabled();
    char message[ 256 ]; /* To hold error messages when reporting is delayed.*/

    message[ 0 ] = '\0';
    failureDisableCalling();
    result = fileOfSocket( socketNumber );

    if ( result == 0 ) {
      sprintf( message, "Can't obtain a FILE pointer to opened socket %d.",
               socketNumber );
    }

    if ( calling ) {
      failureEnableCalling();
    }

    if ( result == 0 ) {
      close( socketNumber ); /* Unusable so close it. */
      socketNumber = -1;
      failureMessage( message );
    }
  }

  return result;
}


/******************************************************************************
PURPOSE: createServerSocket - Create a connected server socket.
INPUTS:  Integer port  The port number to use. E.g., 4321.
RETURNS: Integer       A file descriptor of an open, connected socket,
                       or -1 if failed.
                       Retured file descriptor should be closed using close()
                       when no longer needed.
NOTES:   If unsuccessful then failureMessage() is called.
         This routine blocks awaiting a (single) connection to another
         process (using the same port number). The socket created is uses
         TCP/IP (AF_INET, SOCK_STREAM) for reliable streaming.
******************************************************************************/

static Integer createServerSocket( Integer port ) {

  PRE0( IN_RANGE( port, 1, 65535 ) );

  Integer result = -1;
  int theSocket  = socket( AF_INET, SOCK_STREAM, 0 );
  char message[ 256 ]; /* To hold error messages when reporting is delayed. */

  message[ 0 ] = '\0';

  if ( theSocket == -1 ) {
    sprintf( message, "Can't create server socket." );
  } else {
    const Integer calling = failureCallingEnabled();
    Integer reusable = 0, resized = 0;

    failureDisableCalling();
    reusable = enablePortReusability( theSocket );
    resized = AND2( reusable, setSocketBufferSize( theSocket ) );

    if ( calling ) {
      failureEnableCalling();
    }

    if ( ! reusable ) {
      sprintf( message, "Can't enable reusability for server socket (%d).",
               theSocket );
    } else if ( ! resized ) {
      sprintf( message, "Can't resize buffers of socket (%d).", theSocket );
    } else {
      struct sockaddr_in server;

      CHECK( isUnsignedShort( port ) );
      memset( &server, 0, sizeof server );
      server.sin_family      = AF_INET;
      server.sin_addr.s_addr = INADDR_ANY;
      server.sin_port        = htons( port );

      if ( bind( theSocket, (struct sockaddr*) &server, sizeof server) == -1) {
        sprintf( message,
                 "Can't bind server socket %d to port %"INTEGER_FORMAT".",
                 theSocket, port );
      } else {
        const int backlog = 5; /* HACK? Using just 1 sometimes fails... */

        if ( listen( theSocket, backlog ) == -1 ) {
          sprintf( message, "Can't listen on server socket %d bound to "
                            "port %"INTEGER_FORMAT".", theSocket, port );
        } else {
          socklen_t sizeOfServer = sizeof server;

          result = accept(theSocket,(struct sockaddr*) &server,&sizeOfServer);

          if ( result == -1 ) {
            sprintf( message,
                     "Can't accept connection on server socket %d bound to "
                     "port %"INTEGER_FORMAT".", theSocket, port );
          }
        }
      }
    }

    /* Close first opened socket: */

    CHECK( theSocket != -1 );
    close( theSocket ); /* Ignore failed close. */
    theSocket = -1;
  }

  if ( result == -1 ) {
    failureMessage( message );
  }

  POST0( IN_RANGE( result, -1, INT_MAX ) );

  return result;
}


/******************************************************************************
PURPOSE: createClientSocket - Create a connected client socket.
INPUTS:  Integer port          The port number to use. E.g., 4321.
         const char* hostName  The name of the host to connect to.
                               E.g., "fred" or "localhost".
RETURNS: Integer  File descriptor of an open connected socket, or -1 if failed.
                  Retured file descriptor should be closed using close()
                  when no longer needed.
NOTES:   This routine blocks awaiting a (single) connection to another
         process (using the same port number). The socket created is uses
         TCP/IP (AF_INET, SOCK_STREAM) for reliable streaming.
         If the first attempt to connect to the server fails then this process
         sleeps a few seconds and retries a number of times before giving up.
         This is necessary since due to race conditions the client may attempt
         to connect to the server before the server is ready to listen.
         If unsuccessful then failureMessage() is called.
******************************************************************************/

static Integer createClientSocket( Integer port, const char* hostName ) {

  PRE03( IN_RANGE( port, 1, 65535 ), hostName, *hostName );

  Integer result = -1;
  const struct hostent* const host = gethostbyname( (char*) hostName );
  /* Why isn't gethostbyname() const-correct on all platforms? */

  if ( host == 0 ) {
    failureMessage( "Can't get entry for host '%s'.", hostName );
  } else {
    enum { RETRY_AFTER_SECONDS = 2, RETRIES = 30 }; /* Tunable. */
    Integer failed = 0, retries = RETRIES;

    do {
      result = socket( AF_INET, SOCK_STREAM, 0 );

      if ( result == -1 ) {
        failureMessage( "Can't create client socket." );
        retries = 0;
      } else {
        struct sockaddr_in server;
        CHECK2( isUnsignedShort( port ), IN_RANGE( result, 0, INT_MAX ) );
        memset( &server, 0, sizeof server );
        server.sin_family      = AF_INET;
        server.sin_addr.s_addr = INADDR_ANY;
        server.sin_port        = htons( port );
        memcpy( &server.sin_addr, host->h_addr, host->h_length );
        failed =
          connect( result, (struct sockaddr*) &server, sizeof server ) == -1;

        if ( failed ) {
          close( result ); /* Can't use it so close it (ignore failed close)*/
          sleep( RETRY_AFTER_SECONDS );
          --retries;
        } else {
          retries = 0;
        }
      }
    } while ( retries );

    if ( failed ) {
      const int theUnusableSocket = result;
      result = -1;
      failureMessage( "Can't connect client socket %d"
                      " to port %"INTEGER_FORMAT" on host '%s'.",
                      theUnusableSocket, port, hostName );
    }
  }

  POST0( IN_RANGE( result, -1, INT_MAX ) );

  return result;
}


/******************************************************************************
PURPOSE: enablePortReusability - Specify that a socket should allow re-binding
         on the same port number without delay so that subsequent program runs
         using the same port number will not fail to bind.
INPUTS:  Integer theSocket The socket to enable.
RETURNS: Integer 1 if successful, else 0.
NOTES:   If unsuccessful, failureMessage() is called.
         SGI BUG: This routine is required because SGI sockets, by default,
         disallow re-use of a socket/port soon after (within a minute or so)
         the port was used by a (now defunct) process. Other platforms,
         e.g., Sun, DEC, Cray, etc. do not have this annoying default
         behavior. Nonetheless, this routine establishes a consistent more
         usable behavior on all platforms. Perhaps there is a potential risk
         of receiving old/delayed packets from a previous now defunct process?
         In any case, applications this library was originally written for
         typically establish many short-lived connections and without enabling
         port reusability would soon run out of available ports (at most only
         65,535 ports are available).
******************************************************************************/

static Integer enablePortReusability( Integer theSocket ) {

  PRE0( IN_RANGE( theSocket, 0, INT_MAX ) );

  const Integer result = establishSocketProperty( theSocket, SO_REUSEADDR, 1,
                                                  "reuse mode" ) != 0;

  POST0( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: setSocketBufferSize - Set the socket buffer size to 1MB.
INPUTS:  Integer theSocket The socket to size.
RETURNS: Integer 1 if successful, else 0.
NOTES:   If unsuccessful, failureMessage() is called.
******************************************************************************/

static Integer setSocketBufferSize( Integer theSocket ) {

  PRE0( IN_RANGE( theSocket, 0, INT_MAX ) );

  const Integer size = 262144; /* 256KB is a portable maximum. */

  Integer result =
    establishSocketProperty( theSocket, SO_SNDBUF, size, "send buffer size" );

  if ( result > 0 ) { /* If set to something, set read buffer to same: */
    result =
      establishSocketProperty( theSocket, SO_RCVBUF, result,
                               "receive buffer size" ) > 0;
  }

  POST0( IS_BOOL( result ) );

  return result;
}


/******************************************************************************
PURPOSE: establishSocketProperty - Set a property of a socket.
INPUTS:  Integer theSocket  The socket to establish property of.
         Integer property   The name of a socket property (e.g., SO_REUSEADDR,
                            SO_RCVBUF, etc.)
         Integer value      The value to set the property to.
         const char* propertyName The name of the property - e.g., "reuse mode"
                                  or "buffer size" - used for failureMessage()
                                  messages.
RETURNS: Integer  The established property value (may still be different from
                  value) or 0 upon failure.
NOTES:   If the named property can't be reset then failureMessage() is called.
******************************************************************************/

static Integer establishSocketProperty( Integer theSocket,
                                        Integer property,
                                        Integer value,
                                        const char* propertyName ) {

  PRE05( IN_RANGE( theSocket, 0, INT_MAX ),
         isSignedInt( property ),
         IN_RANGE( value, 1, INT_MAX ),
         propertyName, *propertyName );

  Integer result    = 0;
  int propertyValue = value;
  socklen_t sizeOfArg = sizeof propertyValue;

#ifdef DEBUGGING
  propertyValue = 0;

  if ( getsockopt( theSocket, SOL_SOCKET, property,
                   (char*) &propertyValue, &sizeOfArg ) == 0 ) {
    fprintf( stderr, "\nDefault %s on socket %"INTEGER_FORMAT" is %d.\n",
             propertyName, theSocket, propertyValue );
  }

  propertyValue = value;
#endif

  if ( setsockopt( theSocket, SOL_SOCKET, property,
                   (char*) &propertyValue, sizeOfArg ) != 0 ) {
    failureMessage( "Can't set %s on socket %"INTEGER_FORMAT".",
                    propertyName, theSocket );
  } else if ( getsockopt( theSocket, SOL_SOCKET, property,
                          (char*) &propertyValue, &sizeOfArg ) != 0 ) {
    failureMessage( "Can't verify %s on socket %"INTEGER_FORMAT".",
                    propertyName, theSocket );
  } else if ( propertyValue <= 0 ) {
    failureMessage( "Can't establish %s on socket %"INTEGER_FORMAT".",
                    propertyName, theSocket );
  } else {
    result = propertyValue;

#ifdef DEBUGGING
    if ( result != value ) {
      fprintf( stderr, "\nWarning: Could only set %s on socket %"
               INTEGER_FORMAT" to %"INTEGER_FORMAT".\n",
               propertyName, theSocket, result );
    }
#endif
  }

  DEBUG( fprintf(stderr, "established result = %"INTEGER_FORMAT"\n", result);)

  POST0( IN_RANGE( result, 0, INT_MAX ) );

  return result;
}


/******************************************************************************
PURPOSE: fileOfSocket - Obtain a FILE structure for an opened socket.
INPUTS:  Integer theSocket  Socket descriptor to obtain a FILE for.
RETURNS: FILE*  FILE structure if successful, else 0.
NOTES:   If successful, the returned file structure should be closed with
         fclose() when no longer needed. If unsuccessful, failureMessage()
         is called.
******************************************************************************/

static FILE* fileOfSocket( Integer theSocket ) {

  PRE0( IN_RANGE( theSocket, 0, INT_MAX ) );

  FILE* const result = fdopen( theSocket, "r+" ); /* Sockets are read/write. */

  if ( result == 0 ) {
    failureMessage( "Can't open a file buffer for the"
                    " socket %"INTEGER_FORMAT".", theSocket );
  }

  return result;
}




/******************************************************************************
PURPOSE: matchesWord - Determine if a word matches a member of a
         space-delimited set of words.
INPUTS:  const char* word   The word to search for. E.g., "bar".
         const char* words  The space-delimited set of words to search.
                            E.g., " foo bar baz ".
RETURNS: Integer 1 if word matches one of the words, else 0.
******************************************************************************/

 static Integer matchesWord( const char* word, const char* words ) {

  PRE05(word, *word, words, *words == ' ', words[ strlen(words) - 1 ] == ' ');

  Integer result          = 0;
  const size_t wordLength = strlen( word );
  const char* pointer     = strchr( word, ' ' ) == 0 ? words : 0;

  while ( pointer ) {
    pointer = strstr( pointer, word );

    if ( pointer ) {
      result = AND2( *( pointer - 1 ) == ' ', *( pointer + wordLength) == ' ');

      if ( result ) {
        pointer = 0; /* If whole word found, we're done. */
      } else {
        ++pointer; /* Else keep looking. */
      }
    }
  }

  POST0( IS_BOOL( result ) );

  return result;
}




/******************************************************************************
PURPOSE: swapped8Bytes - Returns the 8-byte Integer with the bytes swapped
         so that byte1 is byte8, byte2 is byte7, etc.
INPUTS:  const Integer x  The 8-byte Integer to swap bytes of.
RETURNS: Integer the byte-swapped value of x.
******************************************************************************/

static Integer swapped8Bytes( Integer x ) {

  const Integer result =
    ( x & INTEGER_CONSTANT( 0xff00000000000000 ) ) >> 56 |
    ( x & INTEGER_CONSTANT( 0x00ff000000000000 ) ) >> 40 |
    ( x & INTEGER_CONSTANT( 0x0000ff0000000000 ) ) >> 24 |
    ( x & INTEGER_CONSTANT( 0x000000ff00000000 ) ) >>  8 |
    ( x & INTEGER_CONSTANT( 0x00000000ff000000 ) ) <<  8 |
    ( x & INTEGER_CONSTANT( 0x0000000000ff0000 ) ) << 24 |
    ( x & INTEGER_CONSTANT( 0x000000000000ff00 ) ) << 40 |
    ( x & INTEGER_CONSTANT( 0x00000000000000ff ) ) << 56;

  CHECK( sizeof (Integer) == 8 );

  POST0( IMPLIES( OR2( x == 0, x == -1 ), result == x ) );

  return result;
}


/******************************************************************************
PURPOSE: swap8Bytes - Swaps the 8-bytes of each Integer.
INPUTS:  const Integer* src   The array of 8-byte Integers to swap.
         Integer        count The number of items in the arrays.
OUTPUTS: Integer*       dst   The array of 8-byte Integers swapped.
******************************************************************************/

static void swap8Bytes( Integer* dst, const Integer* src, Integer count ) {

  PRE03( dst, src, count > 0 );

  for ( ; count--; ++src, ++dst ) {
    *dst = swapped8Bytes( *src );
  }
}



#ifdef _CRAY
/* HACK CRAY BUG: */

/******************************************************************************
PURPOSE: convertXDRRead32BitNonFinite - Detect and correct non-finite values.
INPUTS:  float* x  The possibly non-finite value read by xdr_float.
OUTPUTS: float* x  The finite approximation to the original value.
NOTES:   HACK: CRAY BUG: As of 2000/02/24, Cray C-90, T3D and T3E UNICOS 10
         platforms have broken xdr_vector/xdr_float/xdr_double routines that
         fail to convert valid XDR 32-bit -inf, +inf, nan bit patterns to a
         reasonable bit pattern (e.g., -HUGE_VALF, HUGE_VALF, 0.0). Instead,
         the broken xdr routines (in /usr/lib/libc.a, man xdr_float) convert it
         from the bit pattern on disk to a different bit pattern in memory but
         such 'half-converted' values print as '************' and when they are
         compared to 0.0 they cause a floating-point exception!

         The purpose of this routine is to finish the job - by converting the
         in-memory 'half-converted' bit pattern to a reasonable value,
         specifically: -inf becomes -HUGE_VALF, +inf becomes HUGE_VALF and
         nan becomes 0.0. Except that on Cray C90 platforms, xdr_float converts
         both -inf and +inf to the same in-memory bit-pattern, as a result this
         routine will map both -inf and +inf to -HUGE_VALF.

         It is hoped that in the future Cray will fix the xdr_float routine
         so this  routine can be deleted.
******************************************************************************/

static void convertXDRRead32BitNonFinite( float* x ) {

  DEBUG( fprintf( stderr, "original x = %e\n", *x ); )

#if _CRAYMPP
  {
    const short* const asx = (const short*) x;
    const short         sx = *asx;

    CHECK( sizeof (short) == sizeof (float) );
    DEBUG( fprintf( stderr, "sx = %d\n", sx ); )

    if ( sx == -8388608 ) {
      *x = -HUGE_VALF;
    } else if ( sx == 2143289343 ) {
      *x =  0.0;
    } else if ( sx == 2139095040 ) {
      *x =  HUGE_VALF;
    }
  }
#elif _CRAY && ! _CRAYMPP
  {
    const long long* const allx = (const long long*) x;
    const long long         llx = *allx;

    CHECK( sizeof (long long) == sizeof (float) );
    DEBUG( fprintf( stderr, "llx = %lld\n", llx ); )

    if ( llx == -2305702271725338624L ) {
      *x = -HUGE_VALF;
    } else if ( llx ==  6917529286654119486L ) {
      *x = 0.0;
    } else if ( llx ==  6917669765129437184L ) {
      *x = HUGE_VALF;
    }
  }
#endif

  DEBUG( fprintf( stderr, "converted x = %e\n", *x ); )
}


/******************************************************************************
PURPOSE: convertXDRRead32BitNonFinite - Detect and correct non-finite values.
INPUTS:  float* x  The possibly non-finite value read by xdr_float.
OUTPUTS: float* x  The finite approximation to the original value.
NOTES:   HACK: CRAY BUG: As of 2000/02/24, Cray C-90, T3D and T3E UNICOS 10
         platforms have broken xdr_vector/xdr_float/xdr_double routines that
         fail to convert valid XDR 64-bit -inf, +inf, nan bit patterns to a
         reasonable bit pattern (e.g., -HUGE_VAL, HUGE_VAL, 0.0). Instead, the
         broken xdr routines (in /usr/lib/libc.a, man xdr_double) convert it
         from the bit pattern on disk to a different bit pattern in memory but
         such 'half-converted' values print as '************' and when they are
         compared to 0.0 they cause a floating-point exception!

         The purpose of this routine is to finish the job - by converting the
         in-memory 'half-converted' bit pattern to a reasonable value,
         specifically: -inf becomes -HUGE_VAL, +inf becomes HUGE_VAL and
         nan becomes 0.0.

         It is hoped that in the future Cray will fix the xdr_double routine
         so this  routine can be deleted.
******************************************************************************/

static void convertXDRRead64BitNonFinite( double* x ) {

  PRE0( sizeof (long long) == sizeof (double) );

  const long long* const allx = (const long long*) x;
  const long long         llx = *allx;

  DEBUG( fprintf( stderr, "original x = %le\n", *x ); )
  DEBUG( fprintf( stderr, "llx = %lld\n", llx ); )

#if _CRAYMPP
  if ( llx ==   -4503599627370496L ) {
    *x = -HUGE_VAL;
  } else if ( llx == 9221120237041090559L ) {
    *x = 0.0;
  } else if ( llx == 9218868437227405312L ) {
    *x = HUGE_VAL;
  }
#elif _CRAY && ! _CRAYMPP
  if ( llx == -2305702271725338624L ) {
    *x = -HUGE_VAL;
  } else if ( llx == 6917529286654119486L ) {
    *x = 0.0;
  } else if ( llx == 6917669765129437184L ) {
    *x = HUGE_VAL;
  }
#endif

  DEBUG( fprintf( stderr, "converted x = %le\n", *x ); )
}


/* HACK: CRAY BUG */
#endif


