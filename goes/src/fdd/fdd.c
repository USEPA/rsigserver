/******************************************************************************
PURPOSE: fdd.c - Implements a fast/flexible dd and od.
NOTES:   man dd. man od. fdd help.
To compile:
if ( `uname` == "IRIX64" ) cc -64 -mips4 -xansi -fullwarn -g -o fdd fdd.c
if ( `uname` == "SunOS"  ) cc -xarch=v9 -g -o fdd fdd.c
if ( `uname` == "AIX"    ) cc -q64 -qlanglvl=ansi -qfloat=nomaf -U__STR__ -g -o fdd fdd.c
if ( `uname` == "Linux"  ) cc -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -g -o fdd fdd.c
if ( `uname` == "Darwin" ) cc -mcpu=G5 -mtune=G5 -mpowerpc64 -g -o fdd fdd.c
HISTORY: 2005/04, plessel@computer.org, Created.
STATUS:  unreviewed, untested.
******************************************************************************/

/*=============================== INCLUDES ==================================*/

#include <assert.h>    /* For assert(). */
#include <stdlib.h>    /* For malloc(), free(), strtoul(), strtod(). */
#include <stdarg.h>    /* For va_list, va_start(), va_end(). */
#include <errno.h>     /* For errno. */
#include <stdio.h>     /* For FILE,stdin,fseek(),fread(),vfprintf(),perror()*/
#include <string.h>    /* For memset(), strchr(). */
#include <ctype.h>     /* For isspace(). */
#include <limits.h>    /* For LONG_MAX, ULONG_MAX. */
#include <unistd.h>    /* For unlink(). */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

/*================================== MACROS =================================*/

/* assert_static() is for assertions that can be checked by the compiler: */

#define assert_static( a ) extern int unused_assert_static_[ (a) ? 1 : -1 ]

/* This code won't work on CRAY, for example: */

assert_static( sizeof (char)      == 1 );
assert_static( sizeof (short)     == 2 );
assert_static( sizeof (int)       == 4 );
assert_static( sizeof (long long) == 8 );
assert_static( sizeof (float)     == 4 );
assert_static( sizeof (double)    == 8 );

/* Macros used to simplify assertions and other boolean expressions: */

#ifdef MIN
#undef MIN
#endif
#define MIN( a, b ) ( (a) < (b) ? (a) : (b) )

#define OR2( a, b ) ( (a) || (b) )
#define AND2( a, b ) ( (a) && (b) )
#define AND3( a, b, c ) ( (a) && (b) && (c) )
#define IN3( x, a, b ) ( (x) == (a) || (x) == (b) )
#define IS_BOOL( x ) ( (x) == 0 || (x) == 1 )
#define IMPLIES( p, c ) ( !(p) || (c) )
#define IMPLIES_ELSE( p, c1, c2 ) ( ( (p) && (c1) ) || ( ( !(p) ) && (c2) ) )
#define IN_RANGE( x, low, high ) ( (low) <= (x) && (x) <= (high) )
#define ZERO_OBJECT( x ) memset( (x), 0, sizeof *(x) )
#define IS_SEEKABLE( file ) AND2( (file) != stdin, (file) != stdout )


/* Define appropriate 64-bit integer format specification for printf/scanf: */

#if defined(__osf__)
#define INTEGER_FORMAT "ld"
#elif defined(__OPENNT)
/*
 * strtoll() is missing,
 * but strtoq() exists in libc.a and works for 64-bit integers (typedef quad_t)
 * however, it is not declared in a header so declare it & alias it to strtoll.
 * scanf( "%lld", &a_long_long ) is broken (as is "%qd") for 64-bit integers.
 * Therefore, instead of scanf/sscanf/fscanf() of 64-bit integers,
 * use (aliased) strtoll(), reading a string first if necessary.
 * (Note: printf( "%lld", a_long_long ) works.)
 */
extern long long strtoq( const char*, char**, int );
#define strtoll strtoq
#elif defined(__APPLE__)
#define INTEGER_FORMAT "qd"
#else
#define INTEGER_FORMAT "lld"
#endif


#if defined( _AIX )

/* AIX BUG: fwrite() modulos size_t arguments to 2^31 so we must spoon-feed: */

static size_t buffered_fwrite( const void* array, size_t item_size,
                               size_t count, FILE* stream ) {
  size_t result = 0;
  assert( array ); assert( item_size ); assert( count ); assert( stream );

  /* If there is no overflow then call fwrite multiple times w/ partial data:*/

  if ( item_size <= ULONG_MAX / count && count <= ULONG_MAX / item_size ) {
    const char* data = array;
    const size_t maximum_size = 2147483647; /* 2^31 = 2GB. */
    size_t remainder = item_size * count; /* Won't overflow. */
    int ok = 0;
    assert( remainder >= item_size ); assert( remainder >= count );

    do {
      const size_t bytes = MIN( remainder, maximum_size );
      ok = fwrite( data, bytes, 1, stream ) == 1;
      remainder -= bytes;
      data      += bytes;
    } while ( remainder && ok );

    result = ok * count;
  }

  return result;
}

#define fwrite buffered_fwrite
#endif /* _AIX */

/*================================== TYPES ==================================*/

/* Data conversion modes (all binary, ASCII to/from binary integers/reals): */

enum {
  BINARY,
  ASCII_INTEGER1, ASCII_INTEGER2, ASCII_INTEGER4, ASCII_INTEGER8,
  ASCII_REAL4, ASCII_REAL8,
  INTEGER1_ASCII, INTEGER2_ASCII, INTEGER4_ASCII, INTEGER8_ASCII,
  REAL4_ASCII, REAL8_ASCII,
  MODES
};

#define IS_VALID_MODE( mode ) IN_RANGE( (mode), 0, MODES - 1 )

#define IS_READ_ASCII_MODE( mode ) \
IN_RANGE( (mode), ASCII_INTEGER1, ASCII_REAL8 )

#define IS_WRITE_ASCII_MODE( mode ) \
IN_RANGE( (mode), INTEGER1_ASCII, REAL8_ASCII )

#define MODE_WORD_SIZE(mode) \
  ( IN3( (mode), ASCII_INTEGER1, INTEGER1_ASCII ) ? 1 : \
    IN3( (mode), ASCII_INTEGER2, INTEGER2_ASCII ) ? 2 : \
    IN3( (mode), ASCII_INTEGER4, INTEGER4_ASCII ) ? 4 : \
    IN3( (mode), ASCII_INTEGER8, INTEGER8_ASCII ) ? 8 : \
    IN3( (mode), ASCII_REAL4, REAL4_ASCII ) ? 4 : \
    IN3( (mode), ASCII_REAL8, REAL8_ASCII ) ? 8 : 1 )

#define MAXIMUM_INTEGER_VALUE(mode) \
  ( IN3( (mode), ASCII_INTEGER1, INTEGER1_ASCII ) ? 127LL : \
    IN3( (mode), ASCII_INTEGER2, INTEGER2_ASCII ) ? 32767LL : \
    IN3( (mode), ASCII_INTEGER4, INTEGER4_ASCII ) ? 2147483647LL : \
    IN3( (mode), ASCII_INTEGER8, INTEGER8_ASCII ) ? 9223372036854775807LL:0LL)

#define MINIMUM_INTEGER_VALUE(mode) \
  ( IN3( (mode), ASCII_INTEGER1, INTEGER1_ASCII ) ? -128LL : \
    IN3( (mode), ASCII_INTEGER2, INTEGER2_ASCII ) ? -32768LL : \
    IN3( (mode), ASCII_INTEGER4, INTEGER4_ASCII ) ? -2147483648LL : \
    IN3( (mode), ASCII_INTEGER8, INTEGER8_ASCII ) ?-9223372036854775807LL-1LL\
    : 0LL )

#define MAXIMUM_REAL_VALUE(mode) \
  ( IN3( (mode), ASCII_REAL4, REAL4_ASCII ) ? 3.40282347E+38: \
    IN3( (mode), ASCII_REAL8, REAL8_ASCII ) ? 1.7976931348623157E+308 : 0.0 )

#define MINIMUM_REAL_VALUE(mode) \
  ( IN3( (mode), ASCII_REAL4, REAL4_ASCII ) ? -3.40282347E+38 : \
    IN3( (mode), ASCII_REAL8, REAL8_ASCII ) ? -1.7976931348623157E+308 : 0.0 )

typedef struct Parameters Parameters;

/* Routine that swaps a subset of buffer bytes: */

typedef void (*Swapper)( Parameters*, size_t );

/* ADT/"class" for this program, most routines are "member functions": */

struct Parameters {
  FILE*   inputFile;    /* Input file to read.  */
  FILE*   outputFile;   /* Output file to write. */
  void*   buffer;       /* Buffer allocated for I/O. */
  size_t  bufferSize;   /* Bytes to allocate for I/O buffer. */
  size_t  inputOffset;  /* Byte offset to seek/skip on input file. */
  size_t  outputOffset; /* Byte offset to seek/skip on output file. */
  size_t  count;        /* Bytes to read/write or 0 to read until EOF. */
  Swapper swapper;      /* Routine that swaps buffer bytes. */
  int     mode;         /* BINARY ... ASCII_REAL8. */
  int     ok;           /* Did last command succeed? */
  int truncateOutputFile;     /* Remove it before opening it? */
  const char* outputFileName; /* Needed for re-opening. */
};

/* Table-driven helper routines: */

/* Helpers for parsing command-line arguments: */

typedef void (*ArgumentParser)( const char*, Parameters* );

static void ifParser(    const char* option, Parameters* self );
static void ofParser(    const char* option, Parameters* self );
static void iseekParser( const char* option, Parameters* self );
static void oseekParser( const char* option, Parameters* self );
static void countParser( const char* option, Parameters* self );
static void cbsParser(   const char* option, Parameters* self );
static void convParser(  const char* option, Parameters* self );
static void asciiParser( const char* option, Parameters* self );

typedef struct {
  const char* const option; /* Second part of argument. */
  ArgumentParser    parser; /* Routine to parse option. */
} DispatchEntry;

static const DispatchEntry parsers[] = {
  { "if=",    ifParser },
  { "of=",    ofParser },
  { "iseek=", iseekParser },
  { "oseek=", oseekParser },
  { "count=", countParser },
  { "cbs=",   cbsParser },
  { "conv=",  convParser }
};

/* Helpers for storing ASCII read values into buffer: */

typedef void (*Storer)( void*, size_t, const void* );

/* They're nearly the same, so define a parameterized macro "template": */

#define STORER( name, buffer_type, value_type ) \
static void name##Storer( void* buffer, size_t index, const void* value ) { \
  assert( buffer ); assert( value ); \
  { \
  buffer_type* const tbuffer = buffer; \
  const value_type* const pvalue = value; \
  const value_type tvalue = *pvalue; \
  tbuffer[ index ] = tvalue; \
  } \
}

/* Instantiate macro to define type-specific routines: */

STORER( integer1, signed char, long long )
STORER( integer2, short,       long long )
STORER( integer4, int,         long long )
STORER( integer8, long long,   long long )
STORER( real4,    float,       double )
STORER( real8,    double,      double )

/* Store these routines in a table for later "dynamic dispatch": */

static const Storer storers[ MODES ] = {
  0,
  integer1Storer, integer2Storer, integer4Storer, integer8Storer,
  real4Storer, real8Storer,
  0, 0, 0, 0, 0, 0
};

/* Helpers for printing ASCII-format values in the buffer: */

typedef void (*Printer)( const void*, size_t );

#define PRINTER( name, buffer_type, value_type, format ) \
static void name##Printer( const void* buffer, size_t index ) { \
  assert( buffer ); \
  { \
  const buffer_type* const tbuffer = buffer; \
  const value_type tvalue = tbuffer[ index ]; \
  printf( format, tvalue ); \
  } \
}

PRINTER( integer1, signed char, long long, "%"INTEGER_FORMAT"\n" )
PRINTER( integer2, short,       long long, "%"INTEGER_FORMAT"\n" )
PRINTER( integer4, int,         long long, "%"INTEGER_FORMAT"\n" )
PRINTER( integer8, long long,   long long, "%"INTEGER_FORMAT"\n" )
PRINTER( real4,    float,       double,    "%0.16le\n" )
PRINTER( real8,    double,      double,    "%0.16le\n" )

static const Printer printers[ MODES ] = {
  0, 0, 0, 0, 0, 0, 0,
  integer1Printer, integer2Printer, integer4Printer, integer8Printer,
  real4Printer, real8Printer
};

/* Size in bytes of i/o buffer evenly divisible by largest data word size: */

enum { LARGEST_WORD_SIZE = 16 };
assert_static( LARGEST_WORD_SIZE >= sizeof (double) );
assert_static( LARGEST_WORD_SIZE %  sizeof (double) == 0 );

static const size_t minimumBufferSize = 1024 * 1024; /* 1MB. */
static const size_t maximumBufferSize =
  ULONG_MAX / LARGEST_WORD_SIZE - ULONG_MAX % LARGEST_WORD_SIZE;

static const char* programName = 0;
static int failures = 0; /* Number of program failures. */

/*=========================== FORWARD DECLARATIONS ==========================*/

/* "Member functions" of "class" Parameters: */

static int invariant( const Parameters* self );
static void deallocate( Parameters* self );
static void allocate( Parameters* self );
static void processArguments( int argc, char* argv[], Parameters* self );
static void processArgument( const char* argument, Parameters* self );
static void processFiles( Parameters* self );
static void processSubset( Parameters* self );
static void processAll( Parameters* self );
static void checkWordSizes( Parameters* self, const char* option );
static void readASCII( Parameters* self );
static void readASCIIInteger( Parameters* self, size_t bufferIndex,
                              size_t* bufferBytes );
static void readASCIIReal( Parameters* self, size_t bufferIndex,
                           size_t* bufferBytes );
static void writeASCII( Parameters* self, size_t bytes );
static void swapper8( Parameters* self, size_t count );
static void swapper4( Parameters* self, size_t count );
static void swapper2( Parameters* self, size_t count );
static void seekFiles( Parameters* self );

/* Helpers: */

static void usage( const char* programName );
static int seekFile( FILE* file, size_t offset, void* buffer, size_t size );
static size_t toSizet( const char* string, size_t lower, int* ok );
static long long toInteger( const char* string, long long lower,
                            long long upper, int* ok );
static double toReal(const char* string, double lower, double upper, int* ok);
static const char* skipWhitespace( const char* string );
static void failure( const char* message, ... );

#define SWAPPER_WORD_SIZE( swapper ) \
  ( (swapper) == swapper8 ? 8 \
  : (swapper) == swapper4 ? 4 \
  : (swapper) == swapper2 ? 2 : 1 )

/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: main - Process command-line arguments and files.
INPUTS:  int argc      Number of command-line argument strings.
         char* argv[]  Command-line argument strings.
RETURNS: int 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  int ok = 0;
  Parameters self;
  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  programName = argv[ 0 ];
  processArguments( argc, argv, &self );

  if ( self.ok ) {
    processFiles( &self );
    ok = self.ok;
  }

  deallocate( &self );
  return ! ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: usage - Print program documentation.
INPUTS:  const char* programName  Name of this executable.
******************************************************************************/

static void usage( const char* programName ) {
  assert( programName ); assert( *programName );
  fprintf( stderr, "\a\n\n%s - ", programName );
  fprintf( stderr, "Fast/flexible data dump (like UNIX dd/od commands).\n" );
  fprintf( stderr, "\nusage: %s [option=value] ...\n\n", programName );
  fprintf( stderr, "  Option              Description                     " );
  fprintf( stderr, "[default]\n" );
  fprintf( stderr, "  ----------------------------------------------------" );
  fprintf( stderr, "-------------\n" );
  fprintf( stderr, "  help                Print these instructions.\n" );
  fprintf( stderr, "  if=file             Name of input file.  " );
  fprintf( stderr, "           [stdin]\n" );
  fprintf( stderr, "  of=file             Name of output file. " );
  fprintf( stderr, "           [stdout]\n" );
  fprintf( stderr, "  iseek=bytes         # of input  bytes to skip. " );
  fprintf( stderr, "     [0]\n");
  fprintf( stderr, "  oseek=bytes         # of output bytes to skip. " );
  fprintf( stderr, "     [0]\n");
  fprintf( stderr, "  count=bytes         # of bytes/words to read/write." );
  fprintf( stderr, " [all]\n");
  fprintf( stderr, "                      (In words only if conv=ascii-*)\n");
  fprintf( stderr, "  cbs=bytes           Size of i/o buffer. " );
  fprintf( stderr, "            [1048576]\n");
  fprintf( stderr, "  conv=swab           Byte swap 2-byte words." );
  fprintf( stderr, "         [no swap]\n");
  fprintf( stderr, "  conv=swab2          Byte swap 2-byte words." );
  fprintf( stderr, "         [no swap]\n");
  fprintf( stderr, "  conv=swab4          Byte swap 4-byte words." );
  fprintf( stderr, "         [no swap]\n");
  fprintf( stderr, "  conv=swab8          Byte swap 8-byte words." );
  fprintf( stderr, "         [no swap]\n");
  fprintf( stderr, "  conv=notrunc        Don't remove output file first." );
  fprintf( stderr, " [truncate]\n" );
  fprintf( stderr, "  conv=ascii-integer1 ASCII to 1-byte integer." );
  fprintf( stderr, "        [binary]\n" );
  fprintf( stderr, "  conv=ascii-integer2 ASCII to 2-byte integer." );
  fprintf( stderr, "        [binary]\n" );
  fprintf( stderr, "  conv=ascii-integer4 ASCII to 4-byte integer." );
  fprintf( stderr, "        [binary]\n" );
  fprintf( stderr, "  conv=ascii-integer8 ASCII to 8-byte integer." );
  fprintf( stderr, "        [binary]\n" );
  fprintf( stderr, "  conv=ascii-real4    ASCII to 4-byte real." );
  fprintf( stderr, "           [binary]\n");
  fprintf( stderr, "  conv=ascii-real8    ASCII to 8-byte real." );
  fprintf( stderr, "           [binary]\n");
  fprintf( stderr, "  conv=integer1-ascii 1-byte integer to ASCII." );
  fprintf( stderr, "        [binary]\n" );
  fprintf( stderr, "  conv=integer2-ascii 2-byte integer to ASCII." );
  fprintf( stderr, "        [binary]\n" );
  fprintf( stderr, "  conv=integer4-ascii 4-byte integer to ASCII." );
  fprintf( stderr, "        [binary]\n" );
  fprintf( stderr, "  conv=integer8-ascii 8-byte integer to ASCII." );
  fprintf( stderr, "        [binary]\n" );
  fprintf( stderr, "  conv=real4-ascii    4-byte real to ASCII." );
  fprintf( stderr, "           [binary]\n");
  fprintf( stderr, "  conv=real8-ascii    8-byte real to ASCII." );
  fprintf( stderr, "           [binary]\n");
  fprintf( stderr, "\nExamples:\n\n" );
  fprintf( stderr, "  %s if=data.xdr iseek=123456789", programName );
  fprintf( stderr, " count=4000000000 cbs=104857600 conv=swab4 | readlittle");
  fprintf( stderr, "\n\n" );
  fprintf( stderr, "  Skips 123456789 bytes and then reads 1 billion 4-byte");
  fprintf( stderr, " words\n  using a 100MB buffer and byte swaps each" );
  fprintf( stderr, " 4-byte word\n" );
  fprintf( stderr, "  (presumably for the little-endian host) and writes\n" );
  fprintf( stderr, "  the values to stdout which is piped to a program" );
  fprintf( stderr, " 'readlittle'.\n\n" );
  fprintf( stderr, "  streamer | %s iseek=123456789", programName );
  fprintf( stderr, " count=8000000000 cbs=104857600 conv=swab8 | readlittle");
  fprintf( stderr, "\n\n  Like above but reads from stdin (pipe) and" );
  fprintf( stderr, " processes 8-byte words.\n\n" );
  fprintf( stderr, "  cat data.xdr | %s conv=swab8 conv=integer8-ascii",
           programName );
  fprintf( stderr, " | head\n\n" );
  fprintf( stderr, "  Examine binary files containing 64-bit integers.\n\n" );
  fprintf( stderr, "  echo '104 88 52 1' | %s conv=ascii-integer4",
           programName );
  fprintf( stderr, " > header.bin\n\n" );
  fprintf( stderr, "  Convert ASCII integers to 32-bit binary integers.\n" );
  fprintf( stderr, "\nSupport: plessel@computer.org\n" );
  fprintf( stderr, "\n\n" );
}



/******************************************************************************
PURPOSE: invariant - Is object valid?
INPUTS:  const Parameters* self  Object to validate.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int invariant( const Parameters* self ) {
  int result = self != 0;
  result = AND2( result, self->inputFile );
  result = AND2( result, self->outputFile );
  result = AND2( result, self->bufferSize );
  result = AND2( result, self->bufferSize % LARGEST_WORD_SIZE == 0 );
  result = AND2( result, self->buffer );
  result = AND2( result, IS_VALID_MODE( self->mode ) );
  result = AND2( result, OR2( IS_READ_ASCII_MODE( self->mode ),
                              self->count % MODE_WORD_SIZE(self->mode) == 0));
  result = AND2( result, IS_BOOL( self->ok ) );
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: deallocate - Deallocate object resources and zero object.
INPUTS:  Parameters* self  Object with resources to deallocate.
OUTPUTS: Parameters* self  Object with deallocated and zeroed resources.
******************************************************************************/

static void deallocate( Parameters* self ) {
  assert( self );

  if ( self->buffer ) {
    memset( self->buffer, 0, self->bufferSize );
    free( self->buffer );
  }

  if ( self->inputFile ) {
    fclose( self->inputFile );
  }

  if ( self->outputFile ) {
    fclose( self->outputFile );
  }

  ZERO_OBJECT( self );
  assert( self );
}



/******************************************************************************
PURPOSE: allocate - Allocate object buffer, etc.
INPUTS:  Parameters* self  Object requiring allocated resources.
OUTPUTS: Parameters* self  Allocated initialized object or else self->ok is 0.
******************************************************************************/

static void allocate( Parameters* self ) {
  assert( self ); assert( self->ok ); assert( self->bufferSize );
  assert( self->buffer == 0 ); assert( self->outputFile );
  assert( IMPLIES( AND2( self->truncateOutputFile > 0,
                         IS_SEEKABLE( self->outputFile ) ),
                   self->outputFileName ) );

  self->buffer = malloc( self->bufferSize );

  if ( ! self->buffer ) {
    self->ok = 0;
    failure( "Could not allocate memory for %lu-byte buffer.",
             self->bufferSize );
  } else {
    memset( self->buffer, 0, self->bufferSize );

    /* Truncate output file if specified and possible: */

    if ( AND2( self->truncateOutputFile > 0, IS_SEEKABLE(self->outputFile))) {
      fclose( self->outputFile );
      unlink( self->outputFileName );
      self->outputFile = fopen( self->outputFileName, "wb" );
      self->ok = self->outputFile != 0;

      if ( ! self->ok ) {
        failure( "Could not open/truncate file '%s' for writing.",
                 self->outputFileName );
      }
    }
  }

  assert( IMPLIES( self->ok, invariant( self ) ) );
}



/******************************************************************************
PURPOSE: processArguments - Process command-line arguments.
INPUTS:  int argc      Number of command-line argument strings.
         char* argv[]  Command-line argument strings.
OUTPUTS: Parameters* self  Initialized object.
******************************************************************************/

static void processArguments( int argc, char* argv[], Parameters* self ) {
  size_t argument = 1;
  assert( argc > 0 ); assert( argv ); assert( *argv ); assert( *argv[ 0 ] );
  assert( IMPLIES( argc > 1, AND2( argv[ 1 ], *argv[ 1 ] ) ) );
  assert( self );
  assert( minimumBufferSize >= LARGEST_WORD_SIZE );
  assert( minimumBufferSize < maximumBufferSize );
  assert( minimumBufferSize % LARGEST_WORD_SIZE == 0 );
  assert( maximumBufferSize % LARGEST_WORD_SIZE == 0 );

  ZERO_OBJECT( self );
  self->inputFile  = stdin;
  self->outputFile = stdout;
  self->bufferSize = minimumBufferSize;
  self->ok = 1;

  for ( argument = 1; AND2( self->ok, argument < argc ); ++argument ) {
    processArgument( argv[ argument ], self );
  }

  if ( self->ok ) {
    allocate( self );
  } else {
    usage( argv[ 0 ] );
  }

  assert( IMPLIES( self->ok, invariant( self ) ) );
}



/******************************************************************************
PURPOSE: processFiles - Process input/output files.
INPUTS:  Parameters* self  Object containing parameters for processing.
******************************************************************************/

static void processFiles( Parameters* self ) {
  assert( invariant( self ) ); assert( self->ok );

  seekFiles( self ); /* Seek/skip to byte offsets in input/output files. */

  if ( self->ok ) {

    if ( IS_READ_ASCII_MODE( self->mode ) ) {
      readASCII( self );
    } else if ( self->count ) { /* Read a specified subset of bytes: */
      processSubset( self );
    } else { /* Read until the end of the input file: */
      processAll( self );
    }
  }

  if ( AND2( ! self->ok, failures == 0 ) ) {
    failure( "Failed to read/write all bytes." );
  }

  assert( invariant( self ) );
}



/******************************************************************************
PURPOSE: processSubset - Process a specified subset of input file data.
INPUTS:  Parameters* self  Object containing parameters for processing.
******************************************************************************/

static void processSubset( Parameters* self ) {
  assert( invariant( self ) ); assert( self->ok ); assert( self->count );

  { /* Read/write exactly count bytes: */
    size_t remainder = self->count;

    do {
      const size_t readNow = MIN( remainder, self->bufferSize );
      self->ok = fread( self->buffer, readNow, 1, self->inputFile ) == 1;

      if ( self->ok ) {

        if ( self->swapper ) { /* Apply swapper to buffer if specified: */
          assert( readNow % SWAPPER_WORD_SIZE( self->swapper ) == 0 );
          self->swapper( self, readNow );
        }

        if ( self->mode == BINARY ) {
          self->ok = fwrite( self->buffer, readNow, 1, self->outputFile) == 1;
        } else {
          writeASCII( self, readNow );
        }

        fflush( 0 );
      }

      remainder -= readNow;
    } while ( AND2( self->ok, remainder ) );
  }

  assert( invariant( self ) );
}



/******************************************************************************
PURPOSE: processAll - Process all of input file data.
INPUTS:  Parameters* self  Object containing parameters for processing.
******************************************************************************/

static void processAll( Parameters* self ) {
  assert( invariant( self ) ); assert( self->ok ); assert( self->count == 0 );

  { /* Read/write all bytes: */
    size_t bytesProcessed = 0;
    int done = 0;

    do {
      const size_t count =
        fread( self->buffer, 1, self->bufferSize, self->inputFile );

      if ( count ) {

        if ( self->swapper ) { /* Apply swapper to buffer if specified: */

          if ( count % SWAPPER_WORD_SIZE( self->swapper ) == 0 ) {
            self->swapper( self, count );
          } else {
            self->ok = 0; /* Didn't read enough bytes to make whole words. */
          }
        }

        if ( self->ok ) {

          if ( self->mode == BINARY ) {
            self->ok = fwrite( self->buffer, count, 1, self->outputFile) == 1;
          } else {
            writeASCII( self, count );
          }
        }

        fflush( 0 );
        bytesProcessed += self->ok * count;
      } else {
        done = 1;
      }

    } while ( AND2( self->ok, done == 0 ) );

    self->ok = AND2( self->ok, bytesProcessed );
  }

  assert( invariant( self ) );
}



/******************************************************************************
PURPOSE: processArgument - Process command-line argument.
INPUTS:  const char* argument  Command-line argument to process.
OUTPUTS: Parameters* self      Partially initialized object.
******************************************************************************/

static void processArgument( const char* argument, Parameters* self ) {
  const size_t parserCount = sizeof parsers / sizeof *parsers;
  size_t parser = 0;
  assert( argument ); assert( self ); assert( self->ok );

  /* Find parser by matching first part of argument: */

  do {
    const DispatchEntry* const entry = parsers + parser;
    const char* const option = entry->option;
    const size_t length = strlen( option );

    if ( strncmp( argument, option, length ) == 0 ) {
      const char* const parameter = argument + length;
      assert( entry->parser );
      entry->parser( parameter, self );
      parser = parserCount; /* Found & called parser so finish looping. */
    }

    ++parser;
  } while ( parser < parserCount );

  /* If not found and not already reported as invalid or unsuccessful: */

  if ( parser == parserCount ) {
    self->ok = 0;

    if ( strcmp( argument, "help" ) ) {
      failure( "Invalid argument '%s'.", argument );
    }
  }
}



/******************************************************************************
PURPOSE: ifParser - Process command-line argument "if=".
INPUTS:  const char* option  Portion of command-line argument after "=".
         Parameters* self    Partially initialized object.
OUTPUTS: Parameters* self    Partially initialized object.
******************************************************************************/

static void ifParser( const char* option, Parameters* self ) {
  assert( option ); assert( self ); assert( self->ok );
  self->ok = self->inputFile == stdin;

  if ( ! self->ok ) {
    failure( "Invalid redundant if= argument '%s'.", option );
  } else {
    self->inputFile = fopen( option, "rb" );
    self->ok = self->inputFile != 0;

    if ( ! self->ok ) {
      failure( "Could not open file '%s' for reading.", option );
    }
  }

  assert( IMPLIES( self->ok, self->inputFile ) );
}


/******************************************************************************
PURPOSE: ofParser - Process command-line argument "of=".
INPUTS:  const char* option  Portion of command-line argument after "=".
         Parameters* self    Partially initialized object.
OUTPUTS: Parameters* self    Partially initialized object.
******************************************************************************/

static void ofParser( const char* option, Parameters* self ) {
  assert( option ); assert( self ); assert( self->ok );
  self->ok = self->outputFile == stdout;

  if ( ! self->ok ) {
    failure( "Invalid redundant of= argument '%s'.", option );
  } else {
    struct stat unused;
    self->outputFileName = option;

    if ( stat( self->outputFileName, &unused ) == 0 ) {
      self->outputFile = fopen( option, "rb+" ); /* Seekable to beginning. */

      if ( self->outputFile ) {
        rewind( self->outputFile );
        self->truncateOutputFile += 1;
        /* Truncate unless conv=notrunc occurs before/after this argument. */
      }
    } else {
      self->outputFile = fopen( option, "wb" );
    }

    self->ok = self->outputFile != 0;

    if ( ! self->ok ) {
      failure( "Could not open file '%s' for writing.", option );
    }
  }

  assert( IMPLIES( self->ok, self->outputFile ) );
}



/******************************************************************************
PURPOSE: iseekParser - Process command-line argument "iseek=".
INPUTS:  const char* option  Portion of command-line argument after "=".
         Parameters* self    Partially initialized object.
OUTPUTS: Parameters* self    Partially initialized object.
******************************************************************************/

static void iseekParser( const char* option, Parameters* self ) {
  assert( option ); assert( self ); assert( self->ok );
  self->ok = self->inputOffset == 0;

  if ( ! self->ok ) {
    failure( "Invalid redundant iseek= argument '%s'.", option );
  } else {
    self->inputOffset = toSizet( option, 0, &self->ok );
  }
}



/******************************************************************************
PURPOSE: oseekParser - Process command-line argument "oseek=".
INPUTS:  const char* option  Portion of command-line argument after "=".
         Parameters* self    Partially initialized object.
OUTPUTS: Parameters* self    Partially initialized object.
******************************************************************************/

static void oseekParser( const char* option, Parameters* self ) {
  assert( option ); assert( self ); assert( self->ok );
  self->ok = self->outputOffset == 0;

  if ( ! self->ok ) {
    failure( "Invalid redundant oseek= argument '%s'.", option );
  } else {
    self->outputOffset = toSizet( option, 0, &self->ok );
  }
}



/******************************************************************************
PURPOSE: countParser - Process command-line argument "count=".
INPUTS:  const char* option  Portion of command-line argument after "=".
         Parameters* self    Partially initialized object.
OUTPUTS: Parameters* self    Partially initialized object.
******************************************************************************/

static void countParser( const char* option, Parameters* self ) {
  assert( option ); assert( self ); assert( self->ok );
  self->ok = self->count == 0;

  if ( ! self->ok ) {
    failure( "Invalid redundant count= argument '%s'.", option );
  } else {
    self->count = toSizet( option, 1, &self->ok );
    checkWordSizes( self, option );
  }
}



/******************************************************************************
PURPOSE: cbsParser - Process command-line argument "cbs=".
INPUTS:  const char* option  Portion of command-line argument after "=".
         Parameters* self    Partially initialized object.
OUTPUTS: Parameters* self    Partially initialized object.
******************************************************************************/

static void cbsParser( const char* option, Parameters* self ) {
  assert( option ); assert( self ); assert( self->ok );
  self->ok = self->bufferSize == minimumBufferSize;

  if ( ! self->ok ) {
    failure( "Invalid redundant cbs= argument '%s'.", option );
  } else {
    self->bufferSize = toSizet( option, 1, &self->ok );

    if ( self->bufferSize < minimumBufferSize ) {
      self->bufferSize = minimumBufferSize;
    } else {
      self->bufferSize -= self->bufferSize % LARGEST_WORD_SIZE;
    }

    if ( self->bufferSize > maximumBufferSize ) {
      self->bufferSize = maximumBufferSize;
      assert( maximumBufferSize % LARGEST_WORD_SIZE == 0 );
    }
  }

  assert( IMPLIES( self->ok, self->bufferSize ) );
  assert( IMPLIES( self->ok, self->bufferSize % LARGEST_WORD_SIZE == 0 ) );
}



/******************************************************************************
PURPOSE: convParser - Process command-line argument "conv=".
INPUTS:  const char* option  Portion of command-line argument after "=".
         Parameters* self    Partially initialized object.
OUTPUTS: Parameters* self    Partially initialized object.
******************************************************************************/

static void convParser( const char* option, Parameters* self ) {
  assert( option ); assert( self ); assert( self->ok );

  if ( AND2( ! strcmp( option, "notrunc" ), self->truncateOutputFile >= 0 )) {
    self->truncateOutputFile -= 1;
  } else if ( AND2( self->swapper == 0, ! strncmp( option, "swab", 4 ) ) ) {

    if ( OR2( strcmp( option, "swab" ) == 0, strcmp( option, "swab2") == 0)) {
      self->swapper = swapper2;
    } else if ( strcmp( option, "swab4" ) == 0 ) {
      self->swapper = swapper4;
    } else if ( strcmp( option, "swab8" ) == 0 ) {
      self->swapper = swapper8;
    }

    if ( self->swapper == 0 ) {
      failure( "Invalid value for conv=swab argument '%s'.", option );
      self->ok = 0;
    }

    assert( IMPLIES_ELSE( self->ok, self->swapper, ! self->swapper ) );
  } else if ( AND2( self->mode == BINARY, strstr( option, "ascii" ) != 0 ) ) {
    asciiParser( option, self );
  } else {
    self->ok = 0;
    failure( "Invalid value for conv= argument '%s'.", option );
  }

  checkWordSizes( self, option );
}



/******************************************************************************
PURPOSE: asciiParser - Process command-line argument "conv=*ascii*".
INPUTS:  const char* option  Portion of command-line argument after "=".
         Parameters* self    Partially initialized object.
OUTPUTS: Parameters* self    Partially initialized object.
******************************************************************************/

static void asciiParser( const char* option, Parameters* self ) {
  typedef struct {
    const char* const name;
    const int         value;
  } Pair;
  const Pair table[] = {
    { "ascii-integer1", ASCII_INTEGER1 },
    { "ascii-integer2", ASCII_INTEGER2 },
    { "ascii-integer4", ASCII_INTEGER4 },
    { "ascii-integer8", ASCII_INTEGER8 },
    { "ascii-real4",    ASCII_REAL4    },
    { "ascii-real8",    ASCII_REAL8    },
    { "integer1-ascii", INTEGER1_ASCII },
    { "integer2-ascii", INTEGER2_ASCII },
    { "integer4-ascii", INTEGER4_ASCII },
    { "integer8-ascii", INTEGER8_ASCII },
    { "real4-ascii",    REAL4_ASCII    },
    { "real8-ascii",    REAL8_ASCII    }
  };
  const size_t count = sizeof table / sizeof *table;
  size_t index = 0;
  assert( option ); assert( self ); assert( self->ok );
  assert( strstr( option, "ascii" ) ); assert( self->mode == BINARY );

  do {
    const Pair* const pair = table + index;
    assert( pair ); assert( pair->name );

    if ( ! strcmp( option, pair->name ) ) {
      self->mode = pair->value;
      index = count;
    }

    ++index;
  } while ( index < count );

  if ( self->mode == BINARY ) {
    self->ok = 0;
    failure( "Invalid value for conv= argument '%s'.", option );
  }

  assert( IMPLIES_ELSE(self->ok, self->mode != BINARY, self->mode == BINARY));
}



/******************************************************************************
PURPOSE: checkWordSizes - Verify that word sizes implied by count and conv
         options are compatible.
INPUTS:  Parameters* self   Object containing parameters to check.
         const char* option Argument option being processed.
******************************************************************************/

static void checkWordSizes( Parameters* self, const char* option ) {
  assert( self ); assert( IS_BOOL( self->ok ) );
  assert( IS_VALID_MODE( self->mode ) );

  if ( self->ok ) {
    const int modeWordSize = MODE_WORD_SIZE( self->mode );
    assert( option ); assert( *option );

    if ( IS_READ_ASCII_MODE( self->mode ) ) {
      const int swapWordSize = SWAPPER_WORD_SIZE( self->swapper );

      if ( AND2( self->swapper, swapWordSize != modeWordSize ) ) {
        self->ok = 0;
        failure( "Invalid value for argument option '%s' - "
                 "mismatched implied word sizes (%d) vs (%d).",
                 option, swapWordSize, modeWordSize );
      }
    } else if ( self->count % modeWordSize != 0 ) {
      self->ok = 0;
      failure( "Invalid value for argument option '%s' - "
               "indivisible implied word sizes (%lu) vs (%d).",
               option, self->count, modeWordSize );
    }
  }
}



/******************************************************************************
PURPOSE: readASCII - Read ASCII words and convert/write to binary.
INPUTS:  Parameters* self  Object containing parameters for processing.
NOTES:   ASCII I/O is several orders-of-magnitude slower than BINARY I/O.
******************************************************************************/

static void readASCII( Parameters* self ) {
  size_t bufferBytes = 0;
  size_t bufferIndex = 0;
  size_t count = 0;
  assert( invariant( self ) ); assert( self->ok );
  assert( IS_READ_ASCII_MODE( self->mode ) );

  do {
    int writeBuffer = 0;

    if ( IN_RANGE( self->mode, ASCII_INTEGER1, ASCII_INTEGER8 ) ) {
      readASCIIInteger( self, bufferIndex, &bufferBytes );
    } else {
      readASCIIReal( self, bufferIndex, &bufferBytes );
    }

    bufferIndex += self->ok;
    count       += self->ok;

    /* If the buffer is full or read input EOF then convert/write buffer: */

    writeBuffer += bufferBytes == self->bufferSize;
    writeBuffer += AND2( self->count, count == self->count );
    writeBuffer += AND2( feof( self->inputFile ), bufferBytes );

    if ( writeBuffer) {

      if ( self->swapper ) { /* Apply swapper to buffer if specified: */
        assert( bufferBytes % SWAPPER_WORD_SIZE( self->swapper ) == 0 );
        self->swapper( self, bufferBytes );
      }

      if ( fwrite( self->buffer, bufferBytes, 1, self->outputFile ) != 1 ) {
        self->ok = 0;
      }

      fflush( 0 );
      bufferBytes = 0;
      bufferIndex = 0;
    }

  } while ( AND2( self->ok, IMPLIES( self->count, count < self->count ) ) );

  self->ok = OR2( self->ok,
                  AND2( feof( self->inputFile ),
                        IMPLIES( self->count, count == self->count ) ) );

  assert( invariant( self ) );
}



/******************************************************************************
PURPOSE: readASCIIInteger - Read an ASCII integer.
INPUTS:  Parameters* self         Object containing parameters for processing.
         size_t      bufferIndex  Current index in buffer to store value.
         size_t*     bufferBytes  Current number of bytes stored in buffer.
OUTPUTS: size_t*     bufferBytes  Increased number of bytes stored in buffer.
NOTES:   ASCII I/O is several orders-of-magnitude slower than BINARY I/O.
******************************************************************************/

static void readASCIIInteger( Parameters* self, size_t bufferIndex,
                              size_t* bufferBytes ) {
  char string[ 40 ] = "";
  assert( invariant( self ) ); assert( self->ok );
  assert( IN_RANGE( self->mode, ASCII_INTEGER1, ASCII_INTEGER8 ) );
  assert( MODE_WORD_SIZE( self->mode ) );
  assert( IN_RANGE( self->mode, 0, sizeof storers / sizeof *storers - 1 ) );
  assert( bufferIndex <= self->bufferSize / MODE_WORD_SIZE( self->mode ) );
  assert( bufferBytes );
  assert( *bufferBytes <= self->bufferSize - MODE_WORD_SIZE( self->mode ) );
  memset( string, 0, sizeof string );

  self->ok = AND2( fscanf( self->inputFile, "%22s", string ) == 1,
                   strlen( string ) <= 21 );

  if ( self->ok ) {
    const long long lower = MINIMUM_INTEGER_VALUE( self->mode );
    const long long upper = MAXIMUM_INTEGER_VALUE( self->mode );
    const long long value = toInteger( string, lower, upper, &self->ok );

    if ( self->ok ) {
      const size_t bytes = MODE_WORD_SIZE( self->mode );
      Storer storer = storers[ self->mode ];
      assert( storer );
      storer( self->buffer, bufferIndex, &value );
      *bufferBytes += bytes;
    }
  }

  assert( IN_RANGE( *bufferBytes,
                    self->ok * MODE_WORD_SIZE( self->mode ),
                    self->bufferSize ) );

  assert( invariant( self ) );
}



/******************************************************************************
PURPOSE: readASCIIReal - Read an ASCII real.
INPUTS:  Parameters* self         Object containing parameters for processing.
         size_t      bufferIndex  Current index in buffer to store value.
         size_t*     bufferBytes  Current number of bytes stored in buffer.
OUTPUTS: size_t*     bufferBytes  Increased number of bytes stored in buffer.
NOTES:   ASCII I/O is several orders-of-magnitude slower than BINARY I/O.
******************************************************************************/

static void readASCIIReal( Parameters* self, size_t bufferIndex,
                           size_t* bufferBytes ) {
  char string[ 40 ] = "";
  assert( invariant( self ) ); assert( self->ok );
  assert( IN_RANGE( self->mode, ASCII_REAL4, ASCII_REAL8 ) );
  assert( MODE_WORD_SIZE( self->mode ) );
  assert( IN_RANGE( self->mode, 0, sizeof storers / sizeof *storers - 1 ) );
  assert( bufferIndex <= self->bufferSize / MODE_WORD_SIZE( self->mode ) );
  assert( bufferBytes );
  assert( *bufferBytes <= self->bufferSize - MODE_WORD_SIZE( self->mode ) );
  memset( string, 0, sizeof string );

  self->ok = AND2( fscanf( self->inputFile, "%26s", string ) == 1,
                   strlen( string ) <= 25 );

  if ( self->ok ) {
    const double lower = MINIMUM_REAL_VALUE( self->mode );
    const double upper = MAXIMUM_REAL_VALUE( self->mode );
    const double value = toReal( string, lower, upper, &self->ok );

    if ( self->ok ) {
      const size_t bytes = MODE_WORD_SIZE( self->mode );
      Storer storer = storers[ self->mode ];
      assert( storer );
      storer( self->buffer, bufferIndex, &value );
      *bufferBytes += bytes;
    }
  }

  assert( invariant( self ) );
  assert( IN_RANGE( *bufferBytes,
                    self->ok * MODE_WORD_SIZE( self->mode ),
                    self->bufferSize ) );
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII words from converted binary buffer.
INPUTS:  Parameters* self  Object containing parameters for processing.
         size_t      bytes Number of bytes filled in buffer.
NOTES:   ASCII I/O is several orders-of-magnitude slower than BINARY I/O.
******************************************************************************/

static void writeASCII( Parameters* self, size_t bytes ) {
  size_t count = 0, index = 0;
  Printer printer = 0;
  assert( invariant( self ) ); assert( bytes ); assert( self->ok );
  assert( IS_WRITE_ASCII_MODE( self->mode ) );
  assert( MODE_WORD_SIZE( self->mode ) );
  assert( IN_RANGE( self->mode, 0, sizeof printers / sizeof *printers - 1 ) );

  count = bytes / MODE_WORD_SIZE( self->mode );
  self->ok = count != 0;
  printer = printers[ self->mode ];
  assert( printer );

  for ( index = 0; index < count; ++index ) {
    printer( self->buffer, index );
  }

  assert( invariant( self ) );
}



/******************************************************************************
PURPOSE: swapper8 - Swap byte order of each 8-byte word in the buffer.
INPUTS:  Parameters* self  Object with buffer to convert.
         size_t      bytes Subset of initialized bytes in buffer.
OUTPUTS: Parameters* self  Object with converted buffer.
******************************************************************************/

static void swapper8( Parameters* self, size_t bytes ) {
  assert_static( sizeof (unsigned long long) == 8 );
  assert( invariant( self ) ); assert( self->swapper == swapper8 );
  assert( IN_RANGE( bytes, 8, self->bufferSize ) );
  assert( bytes % sizeof (unsigned long long) == 0 );
  {
    unsigned long long* word = self->buffer;
    size_t count = bytes / sizeof *word;

    do {
      const unsigned long long value = *word;
      const unsigned long long swapped =
        ( value & 0xff00000000000000ULL ) >> 56 |
        ( value & 0x00ff000000000000ULL ) >> 40 |
        ( value & 0x0000ff0000000000ULL ) >> 24 |
        ( value & 0x000000ff00000000ULL ) >>  8 |
        ( value & 0x00000000ff000000ULL ) <<  8 |
        ( value & 0x0000000000ff0000ULL ) << 24 |
        ( value & 0x000000000000ff00ULL ) << 40 |
        ( value & 0x00000000000000ffULL ) << 56;
      *word = swapped;
      ++word;
      --count;
    } while ( count );
  }
}



/******************************************************************************
PURPOSE: swapper4 - Swap byte order of each 4-byte word in the buffer.
INPUTS:  Parameters* self  Object with buffer to convert.
         size_t      bytes Subset of initialized bytes in buffer.
OUTPUTS: Parameters* self  Object with converted buffer.
******************************************************************************/

static void swapper4( Parameters* self, size_t bytes ) {
  assert_static( sizeof (unsigned int) == 4 );
  assert( invariant( self ) ); assert( self->swapper == swapper4 );
  assert( IN_RANGE( bytes, 4, self->bufferSize ) );
  assert( bytes % sizeof (unsigned int) == 0 );
  {
    unsigned int* word = self->buffer;
    size_t count = bytes / sizeof *word;

    do {
      const unsigned int value = *word;
      const unsigned int swapped =
        ( value & 0xff000000 ) >> 24 |
        ( value & 0x00ff0000 ) >>  8 |
        ( value & 0x0000ff00 ) <<  8 |
        ( value & 0x000000ff ) << 24;
      *word = swapped;
      ++word;
      --count;
    } while ( count );
  }
}



/******************************************************************************
PURPOSE: swapper2 - Swap byte order of each 2-byte word in the buffer.
INPUTS:  Parameters* self  Object with buffer to convert.
         size_t      bytes Subset of initialized bytes in buffer.
OUTPUTS: Parameters* self  Object with converted buffer.
******************************************************************************/

static void swapper2( Parameters* self, size_t bytes ) {
  assert_static( sizeof (unsigned short) == 2 );
  assert( invariant( self ) ); assert( self->swapper == swapper2 );
  assert( IN_RANGE( bytes, 2, self->bufferSize ) );
  assert( bytes % sizeof (unsigned short) == 0 );
  {
    unsigned short* word = self->buffer;
    size_t count = bytes / sizeof *word;

    do {
      const unsigned short value = *word;
      const unsigned short swapped =
        ( value & 0xff00 ) >> 8 | ( value & 0x00ff ) << 8;
      *word = swapped;
      ++word;
      --count;
    } while ( count );
  }
}



/******************************************************************************
PURPOSE: seekFiles - Seek/skip to specified byte offset in input/output files.
INPUTS:  Parameters* self  Object containing parameters for processing.
******************************************************************************/

static void seekFiles( Parameters* self ) {
  assert( invariant( self ) ); assert( self->ok );

  if ( self->inputOffset ) {
    self->ok = seekFile( self->inputFile, self->inputOffset,
                         self->buffer, self->bufferSize );
  }

  if ( self->ok ) {

    if ( ! IMPLIES( self->outputOffset, IS_SEEKABLE( self->outputFile ) ) ) {
      failure( "Can't seek on output file." );
      self->ok = 0;
    } else {
      self->ok = seekFile( self->outputFile, self->outputOffset,
                           self->buffer, self->bufferSize );
    }
  }
}



/* Helpers: */



/******************************************************************************
PURPOSE: seekFile - Seek/skip to a given byte offset in a rewound file.
INPUTS:  FILE* file     File to seek. Assumed to be at offset 0 already.
         size_t offset  Byte offset from beginning to seek to / skip.
         void* buffer   Buffer used if file is non-seekable (stdin).
         size_t size    Size in bytes of buffer.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int seekFile( FILE* file, size_t offset, void* buffer, size_t size ) {
  int result = offset == 0;
  assert( file );
  assert( IMPLIES( ! IS_SEEKABLE( file ), AND2( buffer, size ) ) );

  if ( offset ) {
    size_t remainder = offset;

    if ( IS_SEEKABLE( file ) ) {

      do {
        const size_t seekNow = MIN( remainder, LONG_MAX );
        assert( seekNow <= LONG_MAX );
        result = fseek( file, seekNow, SEEK_CUR ) == 0;
        remainder -= seekNow;
      } while ( AND2( result, remainder ) );
    } else {

      do {
        const size_t skipNow = MIN( remainder, size );
        result = fread( buffer, skipNow, 1, file ) == 1;
        remainder -= skipNow;
      } while ( AND2( result, remainder ) );
    }

    if ( ! result ) {
      failure( "Failed to seek/skip to byte offset %lu.", offset );
    }
  }

  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: toSizet - Integer value of string if within range [lower, ULONG_MAX].
INPUTS:  const char* string  - The string to convert.
         size_t      lower   - Minimum acceptable value.
OUTPUTS: int*        ok      - Does string represent an integer within range?
RETURNS: size_t value of string within range, else 0.
NOTES:   Unlike atol(), strtoul() and scanf(), this routine rejects strings
         that would overflow or that contain non-digit characters
         or lack digit characters, or contain multiple whitespace-separated
         tokens. A failure message is printed if string is invalid.
******************************************************************************/

static size_t toSizet( const char* string, size_t lower, int* ok ) {
  size_t result = 0;
  assert( string ); assert( ok );
  *ok = 0;
  errno = 0; /*strtoul sets errno upon failure but won't clear it on success!*/
  string = skipWhitespace( string );

  if ( *string ) {
    char* terminator = 0;

    if ( isdigit( *string ) ) {
      const size_t convertionResult = strtoul( string, &terminator, 10 );

      if ( AND2( terminator, terminator != string ) ) {
       const char* const remainder = skipWhitespace( terminator );

        *ok = AND3( errno != ERANGE, /* No overflow. */
                    *remainder == '\0', /* No remaining characters. */
                    IN_RANGE( convertionResult, lower, ULONG_MAX ) );

        if ( *ok ) {
          result = convertionResult;
        }
      }
    }
  }

  if ( ! *ok ) {
    failure( "Invalid/out-of-range non-negative integer '%s'", string );
  }

  assert( IMPLIES_ELSE( *ok, IN_RANGE(result, lower, ULONG_MAX), result == 0));
  return result;
}



/******************************************************************************
PURPOSE: toInteger - Integer value of string if within range [lower, upper].
INPUTS:  const char* string  - The string to convert.
         long long lower  - The lower limit of valid range.
         long long upper  - The upper limit of valid range.
OUTPUTS: int* ok          - Does string represent an integer in [lower, upper]?
RETURNS: long long value of string within range [lower, upper], else 0.
NOTES:   Unlike atoll(), strtoll() and scanf(), this routine rejects strings
         that would overflow or that contain non-digit characters (except
         an optional leading sign) or lack digit characters, or contain
         multiple whitespace-separated tokens.
******************************************************************************/

static long long toInteger( const char* string, long long lower,
                            long long upper, int* ok ) {
  long long result = 0;
  assert( string ); assert( lower <= upper ); assert( ok );
  *ok = 0;
  errno = 0; /*strtoll sets errno upon failure but won't clear it on success!*/

  if ( *string ) {
    char* terminator = 0;
    const long long convertionResult = strtoll( string, &terminator, 10 );

    if ( AND2( terminator, terminator != string ) ) {
      const char* const remainder = skipWhitespace( terminator );

      *ok = AND3( errno != ERANGE, /* No overflow. */
                  *remainder == '\0', /* No remaining characters. */
                  IN_RANGE( convertionResult, lower, upper ) );

      if ( *ok ) {
        result = convertionResult;
      }
    }
  }

  if ( ! *ok ) {
    failure( "Invalid/out-of-range integer '%s'.", string );
  }

  assert( IMPLIES_ELSE( *ok, IN_RANGE( result, lower, upper ), result == 0 ));
  return result;
}



/******************************************************************************
PURPOSE: toReal - Real value of string if within range [lower, upper].
INPUTS:  const char* string  - The string to convert.
         double lower  - The lower limit of valid range.
         double upper  - The upper limit of valid range.
OUTPUTS: int* ok       - Does string represent a real in [lower, upper]?
RETURNS: double value of string within range [lower, upper], else 0.0.
******************************************************************************/

static double toReal( const char* string, double lower, double upper,
                      int* ok ) {
  double result = 0.0;
  assert( string ); assert( lower <= upper ); assert( ok );
  *ok = 0;

  if ( *string ) {
    char* terminator = 0;
    const double convertionResult = strtod( string, &terminator );

    if ( AND2( terminator, terminator != string ) ) {
      const char* const remainder = skipWhitespace( terminator );

      *ok = AND2( *remainder == '\0', /* No remaining characters. */
                  IN_RANGE( convertionResult, lower, upper ) );

      if ( *ok ) {
        result = convertionResult;
      }
    }
  }

  if ( ! *ok ) {
    failure( "Invalid/out-of-range real '%s'.", string );
  }

  assert( IMPLIES_ELSE( *ok, IN_RANGE( result, lower, upper ), result == 0.0));
  return result;
}



/******************************************************************************
PURPOSE: skipWhitespace - Skip leading spaces in string.
INPUTS:  const char* string  - The string to examine.
RETURNS: const char* pointer into string after leading spaces, if any.
******************************************************************************/

static const char* skipWhitespace( const char* string ) {
  assert( string );

  while ( isspace( *string ) ) {
    ++string;
  }

  assert( ! isspace( *string ) );
  return string;
}



/******************************************************************************
PURPOSE: failure - Print an annotated failure message to stderr and update the
         number of failures and reset errno.
INPUTS:  const char* message  Printf-like format string.
         ... Optional arguments implied by format string.
******************************************************************************/

static void failure( const char* message, ... ) {
  va_list args; /* For stdarg magic. */
  assert( message ); assert( *message );
  assert( programName ); assert( *programName );
  fprintf( stderr, "\a\n\n%s: ", programName );
  va_start( args, message );         /* Begin stdarg magic. */
  vfprintf( stderr, message, args ); /* Forward to vfprintf(). */
  va_end( args );                    /* End of stdarg magic. */

  if ( errno != 0 ) {
    perror( " " );
    errno = 0; /* Clear errno since not all library routines that set it do.*/
  }

  fprintf( stderr, "\n\n" );
  ++failures;
}





