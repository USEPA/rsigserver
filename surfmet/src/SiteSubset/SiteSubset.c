
/******************************************************************************
PURPOSE: SiteSubset.c - Read a subset of a site file (e.g., Airnow, AQS) and
         write it to stdout as XDR (IEEE-754) format binary or ASCII tab-
         delimited spreadsheet.

NOTES:   To compile (debug or optimized versions):

         IRIX64:
           cc -64 -mips4 -xansi -fullwarn -g -o SiteSubset SiteSubset.c

           cc -64 -mips4 -xansi -fullwarn -DNO_ASSERTIONS -O \
              -o SiteSubset SiteSubset.c

         SunOS:
           cc -xarch=v9 -g -o SiteSubset SiteSubset.c

           cc -xarch=v9 -DNO_ASSERTIONS -O -o SiteSubset SiteSubset.c

         Linux.i686:
           gcc -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 \
               -g -o SiteSubset SiteSubset.c

           gcc -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 \
               -DNO_ASSERTIONS -O -o SiteSubset SiteSubset.c

         AIX:
           cc -q64 -qlanglvl=ansi -qfloat=nomaf -U__STR__ \
              -g -o SiteSubset SiteSubset.c

           cc -q64 -qlanglvl=ansi -qfloat=nomaf -U__STR__ \
              -DNO_ASSERTIONS -O -o SiteSubset SiteSubset.c

         Darwin:
           cc -mcpu_type=G5 -mpowerpc64 -g -o SiteSubset SiteSubset.c

           cc -mcpu_type=G5 -mpowerpc64 -DNO_ASSERTIONS -O \
              -o SiteSubset SiteSubset.c

HISTORY: 2005/10 plessel.todd@epa.gov, Created.
STATUS: unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifndef NO_ASSERTIONS
#include <assert.h> /* For function void __assert(). */
#endif
#include <stdlib.h>    /* For malloc(), free(), qsort(). */
#include <stdio.h>     /* For FILE*, fwrite(). */
#include <string.h>    /* For strlen(). */
#include <ctype.h>     /* For isdigit(), isspace(). */
#include <limits.h>    /* For INT_MAX. */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

/*================================== MACROS =================================*/

/* Macro to append 'L' or 'LL' to form appropriate 64-bit integer constant: */

#if _CRAY || __alpha
#define INTEGER_CONSTANT(x) x##L
#define ULONGLONG_CONSTANT(x) x##UL
#else
#define INTEGER_CONSTANT(x) x##LL
#define ULONGLONG_CONSTANT(x) x##ULL
#endif

/*
 * Select the appropriate strtoll() implementation and 64-bit integer format
 * specification for printf/scanf:
 */

#if defined(__osf__) || defined(_CRAY)
/*
 * strtoll() is missing and "%lld" is invalid,
 * but strtol() is 64-bits as is "%ld" so just use them instead.
 */
#define strtoll strtol
#define INTEGER_FORMAT "ld"
#elif defined(__OPENNT)
/*
 * strtoll() is missing,
 * but strtoq() exists in libc.a and works for 64-bit integers (typedef quad_t)
 * however, it is not declared in a header so declare it & alias it to strtoll.
 * scanf( "%lld", &a_long_long ) is BROKEN (as is "%qd") for 64-bit integers.
 * Therefore, instead of scanf/sscanf/fscanf() of 64-bit integers,
 * use (aliased) strtoll(), reading a string first if necessary.
 * (Note: printf( "%lld", a_long_long ) works.)
 */
extern long long strtoq( const char* s, char** unused1, int unused2 );
#define strtoll strtoq
#elif defined(__APPLE__)
/*
 * "%lld" does not work for scanf of 64-bit integers, but "%qd" does so use it.
 * strtoll() is missing but strtoq() works so alias it.
 */
#define INTEGER_FORMAT "qd"
#define strtoll strtoq
#else
#define INTEGER_FORMAT "lld"
#endif

/* Use atoI() instead of scanf() for 64-bit Integer type: */

#define atoI(s) strtoll((s),0,10)



/* Determine if byte-swapping is needed upon output: */

#ifndef IS_LITTLE_ENDIAN
#define IS_LITTLE_ENDIAN \
  ( \
    ( defined(__BYTE_ORDER__) && \
      defined(__ORDER_LITTLE_ENDIAN__) && \
      __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ ) || \
    defined(__x86_64__) || \
    defined(__i386__)   || \
    defined(__i486__)   || \
    defined(__i586__)   || \
    defined(__i686__)   || \
    defined(__alpha)    || \
    defined(__ia64__)   || \
    defined(__ARMEL__)  || \
    defined(__MIPSEL)   || \
    defined(_WIN32)     || \
    defined(_WIN64)     || \
    defined(_M_IX86)    || \
    defined(_M_AMD64)      \
  )
#endif


/* Memory macros: */

#define NEW( type, count ) ((type*) allocate( (count) * sizeof (type) ))
#define FREE( p ) ( (p) ? free( p ) : (void) 0 ), (p) = 0
#define ZERO_OBJECT( x ) memset( (x), 0, sizeof *(x) )



/* Assertions macros for DBC: */

/* assert_static() is for assertions that can be checked by the compiler: */

#define assert_static( a ) extern int unused_assert_static_[ (a) ? 1 : -1 ]

/*
 * Define a custom assert-like macro that includes a prefix tag that identifies
 * the assertion type (PRE, POST, CHECK).
 * This macro calls the "Standard" C Library function void __assert().
 */

#ifndef NO_ASSERTIONS
#define ASSERT2_(prefix,expression) \
((void) ((expression) ? 0 \
:(__assert(#prefix": "#expression,__FILE__,__LINE__), 0)))
#define PRE_ASSERT_(expression) ASSERT2_(PRE,expression)
#define POST_ASSERT_(expression) ASSERT2_(POST,expression)
#define CHECK_ASSERT_(expression) ASSERT2_(CHECK,expression)
#else
#define ASSERT2_(unused1,unused2)
#define PRE_ASSERT_(unused)
#define POST_ASSERT_(unused)
#define CHECK_ASSERT_(unused)
#endif

#ifdef NO_ASSERTIONS

/* Disable all checking. */

#define CHECK(c)
#define CHECK2(c1,c2)
#define CHECK3(c1,c2,c3)
#define CHECK4(c1,c2,c3,c4)
#define CHECK5(c1,c2,c3,c4,c5)
#define CHECK6(c1,c2,c3,c4,c5,c6)
#define CHECK7(c1,c2,c3,c4,c5,c6,c7)
#define CHECK8(c1,c2,c3,c4,c5,c6,c7,c8)
#define CHECK9(c1,c2,c3,c4,c5,c6,c7,c8,c9)

/* An external variable declaration is needed in C so pre can be first line. */

#define NOP_DECLARATION_ extern int errno

#define PRE(c1) NOP_DECLARATION_
#define PRE2(c1,c2) NOP_DECLARATION_
#define PRE3(c1,c2,c3) NOP_DECLARATION_
#define PRE4(c1,c2,c3,c4) NOP_DECLARATION_
#define PRE5(c1,c2,c3,c4,c5) NOP_DECLARATION_
#define PRE6(c1,c2,c3,c4,c5,c6) NOP_DECLARATION_
#define PRE7(c1,c2,c3,c4,c5,c6,c7) NOP_DECLARATION_
#define PRE8(c1,c2,c3,c4,c5,c6,c7,c8) NOP_DECLARATION_
#define PRE9(c1,c2,c3,c4,c5,c6,c7,c8,c9) NOP_DECLARATION_
#define POST CHECK
#define POST2 CHECK2
#define POST3 CHECK3
#define POST4 CHECK4
#define POST5 CHECK5
#define POST6 CHECK6
#define POST7 CHECK7
#define POST8 CHECK8
#define POST9 CHECK9

#define CHECKING(s)

#define OLD(variable)
#define REMEMBER(type,variable) NOP_DECLARATION_
#define REMEMBER_F(type,function_name) NOP_DECLARATION_

#else /* ! NO_ASSERTIONS. */

/* Enable all checking (the default). */

#define CHECK(c) CHECK_ASSERT_(c)

#define CHECK2(c1,c2) CHECK_ASSERT_(c1),CHECK_ASSERT_(c2)

#define CHECK3(c1,c2,c3) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3)

#define CHECK4(c1,c2,c3,c4) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4)

#define CHECK5(c1,c2,c3,c4,c5) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5)

#define CHECK6(c1,c2,c3,c4,c5,c6) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6)

#define CHECK7(c1,c2,c3,c4,c5,c6,c7) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7)

#define CHECK8(c1,c2,c3,c4,c5,c6,c7,c8) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8)

#define CHECK9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
CHECK_ASSERT_(c1),CHECK_ASSERT_(c2),CHECK_ASSERT_(c3),CHECK_ASSERT_(c4),\
CHECK_ASSERT_(c5),CHECK_ASSERT_(c6),CHECK_ASSERT_(c7),CHECK_ASSERT_(c8),\
CHECK_ASSERT_(c9)


#define PRE(c) const int unused_assert_ = (PRE_ASSERT_(c),&unused_assert_!=0)

#define PRE2(c1,c2) \
const int unused_assert_ = (PRE_ASSERT_(c1),PRE_ASSERT_(c2),&unused_assert_!=0)

#define PRE3(c1,c2,c3) \
const int unused_assert_ = \
(PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),&unused_assert_!=0)

#define PRE4(c1,c2,c3,c4) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
&unused_assert_!=0)

#define PRE5(c1,c2,c3,c4,c5) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),&unused_assert_!=0)

#define PRE6(c1,c2,c3,c4,c5,c6) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),&unused_assert_!=0)

#define PRE7(c1,c2,c3,c4,c5,c6,c7) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),&unused_assert_!=0)

#define PRE8(c1,c2,c3,c4,c5,c6,c7,c8) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
&unused_assert_!=0)

#define PRE9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
const int unused_assert_ = (\
PRE_ASSERT_(c1),PRE_ASSERT_(c2),PRE_ASSERT_(c3),PRE_ASSERT_(c4),\
PRE_ASSERT_(c5),PRE_ASSERT_(c6),PRE_ASSERT_(c7),PRE_ASSERT_(c8),\
PRE_ASSERT_(c9),&unused_assert_!=0)


#define POST(c) POST_ASSERT_(c)

#define POST2(c1,c2) POST_ASSERT_(c1),POST_ASSERT_(c2)

#define POST3(c1,c2,c3) POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3)

#define POST4(c1,c2,c3,c4) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4)

#define POST5(c1,c2,c3,c4,c5) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5)

#define POST6(c1,c2,c3,c4,c5,c6) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6)

#define POST7(c1,c2,c3,c4,c5,c6,c7) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7)

#define POST8(c1,c2,c3,c4,c5,c6,c7,c8) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8)

#define POST9(c1,c2,c3,c4,c5,c6,c7,c8,c9) \
POST_ASSERT_(c1),POST_ASSERT_(c2),POST_ASSERT_(c3),POST_ASSERT_(c4),\
POST_ASSERT_(c5),POST_ASSERT_(c6),POST_ASSERT_(c7),POST_ASSERT_(c8),\
POST_ASSERT_(c9)

#define CHECKING(s) s

#define OLD(variable) variable##_old_
#define REMEMBER(type,variable) type OLD(variable) = variable
#define REMEMBER_F(type,function_name) type OLD(function_name)=function_name()


#endif /* NO_ASSERTIONS */


/* Handle optional debugging statements separately. */

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(s)
#endif



/* Logic macros useful for assertions and avoiding precedence defects: */

#define OR2(a,b) ((a)||(b))
#define AND2(a,b) ((a)&&(b))
#define AND3(a,b,c) ((a)&&(b)&&(c))
#define AND4(a,b,c,d) ((a)&&(b)&&(c)&&(d))
#define AND5(a,b,c,d,e) ((a)&&(b)&&(c)&&(d)&&(e))
#define AND6(a,b,c,d,e,f) ((a)&&(b)&&(c)&&(d)&&(e)&&(f))
#define AND7(a,b,c,d,e,f,g) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g))
#define AND8(a,b,c,d,e,f,g,h) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h))
#define AND9(a,b,c,d,e,f,g,h,i) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i))
#define IMPLIES(p,c) (!(p)||(c))
#define IMPLIES_ELSE(p,c1,c2) (((p)&&(c1))||((!(p))&&(c2)))
#define IN_RANGE(x,min,max) ((min)<=(x)&&(x)<=(max))
#define IN3(x,a,b) ((x)==(a)||(x)==(b))
#define IN4(x,a,b,c) ((x)==(a)||(x)==(b)||(x)==(c))
#define IS_BOOL(x) ((x)==0||(x)==1)

/* Convert character to decimal digit: */

#define DIGIT( c ) ( c - '0' )


#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b)?(a):(b))

/*================================== TYPES ==================================*/

typedef long long Integer; /* Exactly 64-bits on all supported platforms. */

/* Output formats: */

enum { OUTPUT_HEADER, OUTPUT_XDR, OUTPUT_ASCII, OUTPUT_FORMATS };
#define IS_VALID_OUTPUT_FORMAT( mode ) IN_RANGE( mode, 0, OUTPUT_FORMATS - 1 )
static const char* const outputFormats[ OUTPUT_FORMATS ] = {
  "-header", "-xdr", "-ascii"
};

/* Dimensions of subset domain[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ]: */

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };

enum { NAME_LENGTH = 32, NOTE_LENGTH = 79 };
typedef char Name[ NAME_LENGTH + 1 ]; /* Variable name. */
typedef char Note[ NOTE_LENGTH + 1 ]; /* Station name/description. */

/* YYYY-MM-DDTHH:MM:SS-ZZZZ */

enum { UTC_TIMESTAMP_LENGTH = 24 };
typedef char UTCTimestamp[ UTC_TIMESTAMP_LENGTH + 1 ];


typedef struct {
  Integer id;
  double longitude;
  double latitude;
  Note note;
} Station;

typedef struct {
  Station station;
  Integer timestamp;
  double value;
  double value2;
} Line;

/* User-supplied command-line arguments: */

typedef struct {
  const char** fileNames;       /* Site ASCII spreadsheet files. */
  const char*  description;     /* User-supplied description. */
  int          days;            /* Number of data files. */
  int          outputFormat;    /* OUTPUT_XDR, OUTPUT_ASCII, OUTPUT_HEADER. */
  Integer      firstTimestamp;  /* YYYYDDDHH00 of subset. */
  size_t       timesteps;       /* Number of hours in subset. */
  double       domain[ 2 ][ 2 ];/*domain[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]*/
} Arguments;

/* Data read/computed: */

typedef struct {
  Integer     firstTimestamp; /* YYYYDDDHHMM of subset. */
  Integer     lastTimestamp;  /* YYYYDDDHHMM of subset. */
  size_t      timesteps;      /* firstTimestamp...lastTimestamp, inclusive. */
  size_t      lineCount;      /* Number of data lines in subset. */
  size_t      stationCount;   /* Number of stations in subset. */
  size_t      fileDataLength; /* Number of characters in fileData. */
  Name        variableName;   /* Name of variable in data file. */
  Name        units;          /* Units of variable in data file. */
  char*       fileData;       /* Contents of data file as a string, edited. */
  Line*       lines;          /* lines[ lineCount ]. */
  Station*    stations;       /* stations[ stationCount ]. */
  double*     data;           /* Subset data[ timesteps ][ stations ]. */
  int         ok;             /* Did last command succeed? */
} Data;

static const double missingValue = -9999.0; /* When no timestamp or station. */

static Integer failureCountDown = 0; /* Simulate allocation failure when > 1.*/

/*
 * 30 days hath September, April, June and November, all the rest have 31,
 * except February which has either 28 or 29 (on a leap year).
 */

static const int daysPerMonth[ 2 ][ 12 ] =
{
  { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
  { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
};

/* This code won't work on CRAY, for example: */

assert_static( sizeof (long long) == 8 );
assert_static( sizeof (double)    == 8 );

/*========================== FORWARD DECLARATIONS ===========================*/

/* Commands: */

static void printUsage( const char* programName );

static void checkForTest( int* argc, char* argv[] );

static int parseArguments( int argc, char* argv[], Arguments* arguments );

static int parseDomain( int argc, char* argv[],
                        int* arg, Arguments* arguments );

static void deallocateData( Data* data );

static void computeTimeRange( const Arguments* arguments, Data* data );

static int parseVariableNameAndUnits( Data* data );

static char* dataLineLength( char* fileData, int* lineLength );

static char* readFiles( const char** fileNames, int count, size_t* length );

static void subsetFileData( const Arguments* arguments, Data* data );

static int lineInSubset( char* dataLine,
                         Integer firstTimestamp, Integer lastTimestamp,
                         const double domain[ 2 ][ 2 ],
                         Integer* timestamp,
                         double* longitude, double* latitude );

static void uniqueStations( Data* data );

static void extractDataValues( Data* data );

static void writeHeader( const Arguments* arguments, const Data* data );

static void writeXDR( const Arguments* arguments, Data* data );

static void writeASCII( const Arguments* arguments, const Data* data );

/* Queries: */

CHECKING( static int isValidArguments( const Arguments* arguments ); )

CHECKING( static int isValidData( const Data* data ); )

static int isValidDomain( const double domain[ 2 ][ 2 ] );

static int stationComparer( const void* a, const void* b );

static int lineComparer( const void* a, const void* b );

static double findValue( const Data* data,
                        Integer yyyydddhhmm,
                        Integer stationId,
                        double* value2 );



/* General utility routines: */

/* Memory routines: */

static void* allocate( size_t bytes );

/* File routines: */

static size_t fileSize( const char* name );
static char* readFile( const char* name, size_t* length );


#if defined( _AIX )

/* AIX BUG: fwrite() modulos size_t arguments to 2^31 so we must spoon-feed: */

static size_t buffered_fwrite( const void* array, size_t item_size,
                               size_t count, FILE* stream ) {
  size_t result = 0;
  PRE4( array, item_size, count, stream );

  /* If there is no overflow then call fwrite multiple times w/ partial data:*/

  if ( item_size <= ULONG_MAX / count && count <= ULONG_MAX / item_size ) {
    const char* data = array;
    const size_t maximum_size = 2147483647; /* 2^31 = 2GB. */
    size_t remainder = item_size * count; /* Won't overflow. */
    int ok = 0;
    CHECK2( remainder >= item_size, remainder >= count );

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

#if IS_LITTLE_ENDIAN
/* static void swapBytes4( void* array, size_t count ); */
static void swapBytes8( void* array, size_t count );
#endif

/* Date/time routines: */

static Integer parseTimestamp( const char* string, Integer* value );
static void toUTCTimestamp( Integer yyyydddhhmm, UTCTimestamp string );
CHECKING( static int isValidTimestamp( Integer yyyydddhhmm ); )
static int isValidDate( Integer yyyymmdd );
static int isLeapYear( int yyyy );
static int convertDate( int yyyymmdd );
static void monthDay( int yyyyddd, int* mm, int* dd );
static void advanceTimestamp( const Integer hours, Integer* yyyydddhhmm );
static void incrementTimestamp( Integer* yyyydddhhmm );
static size_t commasToSpaces( char* string );
static int indexOfString( const char* string, const char* const strings[],
                          int count );
static int linesMatch( const char* line1, const char* line2 );
static const char* skipLine( const char* line, size_t* skipped );
static Integer isValidArgs( Integer argc, const char* argv[] );
static char* skipWords( char* string, int count );

static int asDigit( char c ) {
  int result = c;

  if ( ! isdigit( result ) ) {
    result = '0';
  }

  return result;
}

/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Read a subset of an Site file and write it to stdout in
         XDR or ASCII format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  int ok = 0;

  if ( ! isValidArgs( argc, (const char**) argv ) ) {
    fprintf( stderr, "\a\n\nInvalid command-line arguments.\n" );
  } else {
    Arguments arguments;

    if ( parseArguments( argc, argv, &arguments ) ) {
      Data data;
      ZERO_OBJECT( &data );
      data.fileData =
        readFiles( arguments.fileNames, arguments.days, &data.fileDataLength );

      if ( data.fileData ) {
        data.ok = 1;
        computeTimeRange( &arguments, &data );
        subsetFileData( &arguments, &data );

        if ( data.ok ) {
          uniqueStations( &data );

          if ( data.ok ) {

            /* Sort file data lines by timestamp and then stationId: */

            qsort( data.lines, data.lineCount, sizeof data.lines[ 0 ],
                   lineComparer );

            extractDataValues( &data );

            if ( data.ok ) {

              switch ( arguments.outputFormat ) {
              case OUTPUT_HEADER:
                writeHeader( &arguments, &data );
                break;
              case OUTPUT_XDR:
                writeXDR( &arguments, &data );
                break;
              default:
                CHECK( arguments.outputFormat == OUTPUT_ASCII );
                writeASCII( &arguments, &data );
                break;
              }

              ok = data.ok;
            }
          }
        }
      }

      deallocateData( &data );
    }
  }

  POST( IS_BOOL( ok ) );

  return ! ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* programName  Name of program.
******************************************************************************/

static void printUsage( const char* programName ) {
  PRE( programName );
  fprintf( stderr, "\n\n\n%s - Read Site (e.g., Airnow, AQS, etc.) files\n",
           programName );
  fprintf( stderr, "and write the specified subset of stations/data\n" );
  fprintf( stderr, "to stdout in XDR or ASCII format.\n" );
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-data <file_name> [<file_name> ... <file_name>] \\\n" );
  fprintf( stderr, "-header | -xdr | -ascii \\\n" );
  fprintf( stderr, "-desc \"description text\" \\\n" );
  fprintf( stderr, "-timestamp <yyyymmddhh> -hours <count> \\\n" );
  fprintf( stderr, "[ -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> ] \\\n" );
  fprintf( stderr, "Note: timestamp is in UTC (GMT)\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example #1:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-data ../../../../data/20050826.PM25.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050827.PM25.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050828.PM25.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050829.PM25.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050830.PM25.txt \\\n" );
  fprintf( stderr, "-xdr -desc \"New Orleans, LA\" \\\n" );
  fprintf( stderr, "-timestamp 2005082600 -hours 120 \\\n" );
  fprintf( stderr, "-domain -92 28 -89 31 > subset.xdr\n" );
  fprintf( stderr, "\nSubset of days 2005, August 26-30 for\n" );
  fprintf( stderr, "\narea near New Orleans, LA.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays\n" );
  fprintf( stderr, "lonlats[stations]\n");
  fprintf( stderr, "data[timesteps][stations]\n");
  fprintf( stderr, "For example:\n" );
  fprintf( stderr, "SITE 2.0\n" );
  fprintf( stderr, "New Orleans, LA.\n" );
  fprintf( stderr, "2005-08-26T00:00:00-0000\n" );
  fprintf( stderr, "# data dimensions: timesteps stations:\n" );
  fprintf( stderr, "120 7\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "PM25\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "ug/m3\n" );
  fprintf( stderr, "# char notes[stations][80] and\n" );
  fprintf( stderr, "# MSB 32-bit integers ids[stations] and\n" );
  fprintf( stderr, "# IEEE-754 32-bit reals " );
  fprintf( stderr, "sites[stations][2=<longitude,latitude>] and\n" );
  fprintf( stderr, "# IEEE-754 32-bit reals data[timesteps][stations]:\n");
  fprintf( stderr, "<binary data arrays here>\n\n\n");
  fprintf( stderr, "Example #2:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-data ../../../../data/20050826.PM25.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050827.PM25.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050828.PM25.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050829.PM25.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050830.PM25.txt \\\n" );
  fprintf( stderr, "-header -desc \"New Orleans, LA\" \\\n" );
  fprintf( stderr, "-timestamp 2005082600 -hours 120 \\\n" );
  fprintf( stderr, "-domain -92 28 -89 31 > subset.txt\n" );
  fprintf( stderr, "\nSame as above but only outputs ASCII header.\n\n\n" );
  fprintf( stderr, "Example #3:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-data ../../../../data/20050826.Ozone.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050827.Ozone.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050828.Ozone.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050829.Ozone.txt \\\n" );
  fprintf( stderr, "      ../../../../data/20050830.Ozone.txt \\\n" );
  fprintf( stderr, "-ascii -desc US \\\n" );
  fprintf( stderr, "-timestamp 2005082600 -hours 120 \\\n" );
  fprintf( stderr, "> subset.txt\n" );
  fprintf( stderr, "\nLike above but outputs ozone in a spreadsheet\n");
  fprintf( stderr, "importable format (tab-separated values).\n");
  fprintf( stderr, "\n\n\n");
}



/******************************************************************************
PURPOSE: isValidArgs() - Are each of the strings non-zero and non-zero-length?
INPUTS:  Integer argc              Number of strings to check.
         const char* argv[ argc ]  Strings to check.
RETURNS: Integer 1 if non-zero, else 0.
******************************************************************************/

static Integer isValidArgs( Integer argc, const char* argv[] ) {
  Integer result = 0;

  if ( AND2( IN_RANGE( argc, 1, INT_MAX - 1 ), argv != 0 ) ) {
    Integer index = 0;

    do {

      if ( OR2( argv[ index ] == 0, argv[ index ][ 0 ] == '\0' ) ) {
        index = argc;
      }

      ++index;
    } while ( index < argc );

    result = index == argc;
  }

  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: checkForTest - Check for and set-up for testing.
INPUTS:  int* argc     Number of command-line arguments.
         char* argv[]  Command-line argument strings.
OUTPUTS: int* argc     Possibly reduced count.
NOTES:   Used to set-up simulated memory allocation failures to attain
         coverage of failure-handling code.
******************************************************************************/

static void checkForTest( int* argc, char* argv[] ) {

  PRE( isValidArgs( *argc, (const char**) argv ) );

  if ( AND3( *argc >= 3, strcmp( argv[ *argc - 2 ], "-test" ) == 0,
             atoI( argv[ *argc - 1 ] ) > 0 ) ) {
    const Integer count = atoI( argv[ *argc - 1 ] );
    *argc -= 2; /* Remove last two arguments: -test 3 */
    failureCountDown = count; /* Simulate failure on count NEW(). */
  }

  POST( *argc > 0 );
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
INPUTS:  int argc              Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int parseArguments( int argc, char* argv[], Arguments* arguments ) {

  PRE2( isValidArgs( argc, (const char**) argv ), arguments );

  int result = 0;

  checkForTest( &argc, argv );
  ZERO_OBJECT( arguments );

  if ( argc >= 10 ) {
    int arg = 1;

    if ( AND4( argv[ 1 ], strcmp( argv[ 1 ], "-data" ) == 0,
               argv[ 2 ][ 0 ] != '-',
               strcmp( argv[ 2 ], argv[ 0 ] ) != 0 ) ) {
      arg = 2;
      arguments->fileNames = (const char**) argv + arg;
      arguments->days = 0;

      while ( AND4( arg < argc, argv[ arg ], *argv[ arg ] != '-',
                    strcmp( argv[ arg ], argv[ 0 ] ) ) ) {
        arguments->days += 1;
        ++arg;
      }

      if ( argc - arg >= 7 ) {
        arguments->outputFormat =
          indexOfString( argv[ arg ], outputFormats,
                         sizeof outputFormats / sizeof *outputFormats );

        if ( ! IS_VALID_OUTPUT_FORMAT( arguments->outputFormat ) ) {
          fprintf( stderr, "\a\n\nInvalid output format '%s'.\n", argv[ 3 ] );
        } else {
          ++arg;

          if ( AND4( argv[ arg ], strcmp( argv[ arg ], "-desc" ) == 0,
                     argv[ arg + 1 ][ 0 ] != '-',
                     ! strchr( argv[ arg + 1 ], '\n' ) ) ) {
            ++arg;
            arguments->description = argv[ arg ];
            ++arg;

            if ( AND3( argv[ arg ], strcmp( argv[ arg ], "-timestamp" ) == 0,
                       argv[ arg + 1 ][ 0 ] != '-' ) ) {
              ++arg;

              if ( parseTimestamp( argv[ arg ],
                                   &arguments->firstTimestamp ) ) {
                ++arg;

                if ( AND3( argv[ arg ], strcmp( argv[ arg ], "-hours" ) == 0,
                           argv[ arg + 1 ][ 0 ] != '-' ) ) {
                  ++arg;
                  arguments->timesteps = atoI( argv[ arg ] );

                  if ( arguments->timesteps < 1 ) {
                    fprintf( stderr, "\a\n\nInvalid hours specified '%s'.\n",
                             argv[ arg ] );
                  } else {
                    ++arg;
                    result = parseDomain( argc, argv, &arg, arguments );
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    ZERO_OBJECT( arguments );
    printUsage( argv[ 0 ] );
  }

  POST2( IS_BOOL( result ), IMPLIES( result, isValidArguments( arguments ) ));
  return result;
}



/******************************************************************************
PURPOSE: parseDomain - Parse command-line arguments for domain, if specified,
         else initialize domain to entire Earth.
INPUTS:  int argc              Number of command-line arguments.
         char* argv[]          Command-line argument strings.
         int* arg              Current argument to parse.
OUTPUTS: int* arg              Updated current argument to parse.
         Arguments* arguments  Initialized arguments.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int parseDomain( int argc, char* argv[], int* arg,
                        Arguments* arguments ) {

  PRE5( isValidArgs( argc, (const char**) argv ),
        argc >= 10, arguments, arg,
        IN_RANGE( *arg, 10, argc ) );

  int result = 0;

  if ( AND2( *arg + 4 < argc, strcmp( argv[ *arg ], "-domain" ) == 0 ) ) {
    sscanf( argv[ *arg + 1 ], "%lf", &arguments->domain[ LONGITUDE ][ MINIMUM]);
    sscanf( argv[ *arg + 2 ], "%lf", &arguments->domain[ LATITUDE  ][ MINIMUM]);
    sscanf( argv[ *arg + 3 ], "%lf", &arguments->domain[ LONGITUDE ][ MAXIMUM]);
    sscanf( argv[ *arg + 4 ], "%lf", &arguments->domain[ LATITUDE  ][ MAXIMUM]);

    if ( ! isValidDomain( (const double(*)[2]) arguments->domain ) ) {
      fprintf( stderr, "\a\n\nInvalid domain specified.\n" );
    } else {
      *arg += 5;
      result = 1;
    }
  } else {
    arguments->domain[ LONGITUDE ][ MINIMUM ] = -180.0;
    arguments->domain[ LONGITUDE ][ MAXIMUM ] =  180.0;
    arguments->domain[ LATITUDE  ][ MINIMUM ] =  -90.0;
    arguments->domain[ LATITUDE  ][ MAXIMUM ] =   90.0;
    result = 1;
  }

  POST3( IS_BOOL( result ),
         IMPLIES( result, isValidArguments( arguments ) ),
         IN_RANGE( *arg, 10, argc ) );
  return result;
}



/******************************************************************************
PURPOSE: deallocateData - Deallocate and zero contents of data.
INPUTS:  Data* data Data structure to examine and deallocate/zero members.
OUTPUTS: Data* data Data structure with deallocated/zeroed members.
******************************************************************************/

static void deallocateData( Data* data ) {
  PRE( data );

  FREE( data->fileData );
  FREE( data->lines );
  FREE( data->stations );
  FREE( data->data );

  ZERO_OBJECT( data );
  POST( data );
}



/******************************************************************************
PURPOSE: computeTimeRange - Compute last timestamp of subset.
INPUTS:  Arguments* arguments  The firstTimestamp and timesteps.
OUTPUTS: Data* data            The firstTimestamp, lastTimestamp, timesteps.
******************************************************************************/

static void computeTimeRange( const Arguments* arguments, Data* data ) {

  PRE5( isValidArguments( arguments ),
        data, data->ok, data->firstTimestamp == 0, data->lastTimestamp == 0 );

  int timestep = 0;
  data->timesteps      = arguments->timesteps;
  data->firstTimestamp = arguments->firstTimestamp;
  data->lastTimestamp  = data->firstTimestamp;

  for ( timestep = 1; timestep < data->timesteps; ++timestep ) {
    incrementTimestamp( &data->lastTimestamp );
  }

  POST5( data->ok, data->timesteps > 0,
         isValidTimestamp( data->firstTimestamp ),
         isValidTimestamp( data->lastTimestamp ),
         data->lastTimestamp >= data->firstTimestamp );
}



/******************************************************************************
PURPOSE: readFiles - Read Site ASCII data files into buffer.
INPUTS:  const char** names   Names of data files to read.
         int          count   Number of data files.
OUTPUTS: size_t*      length  Number of characters in result.
RETURNS: char* concatenated string buffer of files if successful, else 0.
NOTES:   Result must be FREE'd by caller. LEAK.
******************************************************************************/

static char* readFiles( const char** names, int count, size_t* length ) {
  PRE5( names, names[ 0 ], count > 0, names[ count - 1 ], length );
  char* result = readFile( names[ 0 ], length );
  int ok = 1;

  if ( result ) {
    int file = 1;

    for ( file = 1; AND2( ok, file < count ); ++file ) {
      size_t size = 0;
      char* contents = readFile( names[ file ], &size );
      ok = 0;

      if ( contents ) {

        if ( ! linesMatch( result, contents ) ) { /* Same header line? */
          fprintf( stderr, "\a\nInvalid/mismatched input data file '%s'.\n\n",
                   names[ file ] );
        } else {
          size_t skipped = 0;
          const char* data = skipLine( contents, &skipped ); /* Skip header.*/
          char* buffer = 0;
          *length += size - skipped;
          buffer = NEW( char, *length + 1 );

          if ( buffer ) {
            strcpy( buffer, result );
            strcat( buffer, data );
            FREE( contents );
            FREE( result );
            result = buffer;
            buffer = 0;
            data = 0;
            ok = 1;
          }
        }

        if ( ! ok ) {
          FREE( contents );
        }
      }
    }
  }

  if ( ! ok ) {
    FREE( result );
    *length = 0;
  }

  POST( IMPLIES_ELSE( result, AND2( *result, *length ), *length == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: parseVariableNameAndUnits - Parse the name and units of the variable
         from the file header.
INPUTS:  Data* data Data->fileData header line.
OUTPUTS: Data* data Data->variableName initialized name.
                    Data->units initialized name.
RETURNS: int 1 if successful, else 0.
NOTES: Header line looks like this:
SITE LATITUDE LONGITUDE YEAR JUL_DAY GMT_HR PM25_1HR ug/m3 SITE_NAME
or, for wind:
SITE LATITUDE LONGITUDE YEAR JUL_DAY GMT_HR wind_u wind_v m/s SITE_NAME
******************************************************************************/

static int parseVariableNameAndUnits( Data* data ) {
  PRE5( data->fileData, data->variableName, data->units,
        *data->variableName == '\0', *data->units == '\0' );
  int result = 0;
  char* newline = strchr( data->fileData, '\n' );

  if ( newline ) {
    char* word = 0;
    *newline = '\0'; /* Terminate header line for speed. */
    word = skipWords( data->fileData, 6 );

    if ( word ) { /* Copy the variable name: */
      int i = 0;

      for ( i = 0;
            AND3( i < NAME_LENGTH, word[ i ], ! isspace( word[ i ] ) );
            ++i ) {
        data->variableName[ i ] = word[ i ];
      }

      data->variableName[ NAME_LENGTH ] = '\0';



      if ( data->variableName[ 0 ] ) {

        if ( strstr( data->variableName, "wind_u" ) == data->variableName ) {
          data->variableName[ 4 ] = '\0'; /* "wind". units = "m\s". */
          data->units[ 0 ] = 'm';
          data->units[ 1 ] = '/';
          data->units[ 2 ] = 's';
          data->units[ 3 ] = '\0';
          result = 1;
        } else { /* Parse the variable units: */
          word = skipWords( word, 1 );

          if ( word ) {

            for ( i = 0;
                  AND3( i < NAME_LENGTH, word[ i ], ! isspace( word[ i ] ) );
                  ++i ) {
              data->units[ i ] = word[ i ];
            }

            data->units[ NAME_LENGTH ] = '\0';
            result = data->units[ 0 ] != '\0';
          }
        }
      }
    }

    *newline = '\n'; /* Restore newline at the end of the header line: */
  }

  if ( ! result ) {
    fprintf( stderr, "\a\nInvalid input data file header line.\n\n" );
  }

  POST2( IS_BOOL( result ),
         IMPLIES_ELSE( result,
                       AND2( IN_RANGE( strlen( data->variableName ),
                                       1, NAME_LENGTH ),
                             IN_RANGE( strlen( data->units ),
                                       1, NAME_LENGTH ) ),
                       AND2( *data->variableName == '\0',
                             *data->units == '\0' ) ) );
  return result;
}


/******************************************************************************
PURPOSE: dataLineLength - Get the first data line and its length.
INPUTS:  char* fileData  File data string.
OUTPUTS: int& lineLength Length of first data line.
RETURNS: char* first data line if successful, else 0.
******************************************************************************/

static char* dataLineLength( char* fileData, int* lineLength ) {
  PRE2( fileData, lineLength );
  char* result = 0;
  char* newline = strchr( fileData, '\n' );
  *lineLength = 0;

  if ( newline ) {
    result = newline + 1;
    newline = strchr( result, '\n' );

    if ( newline ) {
      *lineLength = (int) ( newline - result ) + 1;
    } else {
      result = 0;
    }
  }

  POST( IMPLIES( result, lineLength > 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: subsetFileData - Extract the lines of data within the timetamp/domain.
INPUTS:  Arguments* arguments  domain of subset.
         Data*      data       fileData, firstTimestamp, lastTimestamp.
OUTPUTS: Data*      data       fileData        Deallocated.
                    data       fileDataLength  Zero.
                    data       variableName    Initialized name.
                    data       units           Initialized units.
                    lines      Allocated array of data lines.
                    lineCount  Number of items in lines.
                    ok         1 if successful, else 0.
NOTES: See reformat_airnow script, /data/FAQSD/extract, etc.
******************************************************************************/

static void subsetFileData( const Arguments* arguments, Data* data ) {

  PRE9( isValidArguments( arguments ), data, data->fileData,
        isValidTimestamp( data->firstTimestamp ),
        isValidTimestamp( data->lastTimestamp ),
        data->firstTimestamp <= data->lastTimestamp,
        isValidDomain( arguments->domain ),
        data->lines == 0,
        data->lineCount == 0 );

  size_t lines = commasToSpaces( data->fileData );
  size_t subsetLineCount = 0;
  data->ok = 0;

  if ( AND2( lines > 1, parseVariableNameAndUnits( data ) ) ) {

    /*
     * For all subsequent (data) lines (that must be the same length),
     *   change time fields into a single (long long) timestamp.
     * So:
     *   01234567890123456789012345678901234567890123456789012
     *   350130021   31.7961 -106.5839 2004 244 17.5     -1.0
     * becomes:
     *   350130021   31.7961 -106.5839 20042441700       -1.0
     */

    int lineLength = 0;
    char* const firstDataLine = dataLineLength( data->fileData, &lineLength );

    if ( firstDataLine ) {
      const char* const end = data->fileData + data->fileDataLength;
      char* dataLine = firstDataLine;
      --lines;

      while ( dataLine + lineLength <= end ) {


        if ( dataLine[ lineLength - 1 ] == '\n' ) {
          char* timestamp = 0;
          dataLine[ lineLength - 1 ] = '\0';

          timestamp = skipWords( dataLine, 3 );

          if ( AND2( timestamp, timestamp + 13 < dataLine + lineLength ) ) {
            timestamp[ 4 ] = asDigit( timestamp[ 5 ] ); /* ddd */
            timestamp[ 5 ] = asDigit( timestamp[ 6 ] );
            timestamp[ 6 ] = asDigit( timestamp[ 7 ] );
            timestamp[ 7 ] = asDigit( timestamp[ 9 ] ); /* hh */
            timestamp[ 8 ] = asDigit( timestamp[ 10 ] );
            timestamp[ 9 ] = '0';
            timestamp[ 10 ] = '0';
            timestamp[ 11 ] = ' ';
            timestamp[ 12 ] = ' ';

            /* Count if line is within the given timestamps/domain: */

            subsetLineCount +=
              lineInSubset( dataLine,
                            data->firstTimestamp, data->lastTimestamp,
                           (const double (*)[2]) arguments->domain, 0, 0, 0 );
          }

          dataLine[ lineLength - 1 ] = '\n';
        } else {
          break;
        }

        dataLine += lineLength;
      }

      /* If there were lines in the subset, make a 2nd pass and copy them: */

      if ( subsetLineCount ) {
        data->lines = NEW( Line, subsetLineCount );

        if ( data->lines ) {
          const int is_wind = ! strcmp( data->variableName, "wind" );
          Line* destination = data->lines;
          char* dataLine = firstDataLine;

          while ( dataLine + lineLength <= end ) {
            Integer timestamp = 0;
            double longitude = 0.0;
            double latitude = 0.0;
            dataLine[ lineLength - 1 ] = '\0';

            if ( lineInSubset( dataLine,
                               data->firstTimestamp, data->lastTimestamp,
                               (const double (*)[2]) arguments->domain,
                               &timestamp, &longitude, &latitude ) ) {
              char* word = skipWords( dataLine, 4 );
              int ok = word != 0;

              if ( word ) {
                Station* const station = &destination->station;
                station->id = atoI( dataLine );
                station->longitude = longitude;
                station->latitude = latitude;
                destination->timestamp = timestamp;
                destination->value = atof( word );

                if ( is_wind ) {
                  word = skipWords( word, 1 );

                  if ( word ) {
                    destination->value2 = atof( word );
                  } else {
                    ok = 0;
                  }
                }

                if ( ok ) {
                  word = skipWords( word, 1 );
                  ok = word != 0;

                  if ( ok ) { /* Read note: */
                    strncpy( station->note, word,
                            sizeof station->note / sizeof station->note[0] - 1);
                    ++destination;
                    data->lineCount += 1;
                  }
                }
              }
            }

            dataLine[ lineLength - 1 ] = '\n';
            dataLine += lineLength;
          }

          data->ok = 1;
        }
      }
    }
  }

  /*
   * Free memory for fileData now since it is not needed by subsequent
   * routines and its (significant) memory can be reused by subsequent
   * allocations. (Avoid increasing process memory 'high-water mark'.)
   */

  FREE( data->fileData );
  data->fileDataLength = 0;

  data->ok = AND2( data->lineCount > 0, data->lines != 0 );

  if ( data->lineCount == 0 ) {
    fprintf( stderr,
            "\a\n\nThere are no data lines within the specified subset.\n\n");
  }

  POST4( data->fileData == 0,
         data->fileDataLength == 0,
         IMPLIES( data->ok,
           AND9( data->variableName[ 0 ] && data->units[ 0 ],
                 data->lines,
                 data->lineCount > 0,
                 data->lines[ 0 ].station.id > 0,
                 data->lines[ data->lineCount - 1 ].station.id > 0,
                 IN_RANGE( data->lines[ 0 ].station.longitude,
                           arguments->domain[ LONGITUDE ][ MINIMUM ],
                           arguments->domain[ LONGITUDE ][ MAXIMUM ] ),
                 IN_RANGE( data->lines[ 0 ].station.latitude,
                           arguments->domain[ LATITUDE ][ MINIMUM ],
                           arguments->domain[ LATITUDE ][ MAXIMUM ] ),
                 IN_RANGE( data->lines[ data->lineCount - 1].station.longitude,
                           arguments->domain[ LONGITUDE ][ MINIMUM ],
                           arguments->domain[ LONGITUDE ][ MAXIMUM ] ),
                 IN_RANGE( data->lines[ data->lineCount - 1].station.latitude,
                           arguments->domain[ LATITUDE ][ MINIMUM ],
                           arguments->domain[ LATITUDE ][ MAXIMUM ] ) ) ),
         IMPLIES( data->ok,
           AND6( data->lines[ 0 ].station.note,
                 data->lines[ 0 ].station.note[ 0 ],
                 data->lines[ 0 ].station.note[ NOTE_LENGTH ] == '\0',
                 data->lines[ data->lineCount - 1].station.note,
                 data->lines[ data->lineCount - 1].station.note[ 0 ],
                 data->lines[ data->lineCount - 1].station.note[ NOTE_LENGTH]
                   == '\0' ) ) );
}



/******************************************************************************
PURPOSE: lineInSubset - Is the line of data within the timetamp/domain?
INPUTS:  char* dataLine                Line of the data file to examine.
         Integer     firstTimestamp    First timestamp of subset.
         Integer     lastTimestamp     Last  timestamp of subset.
         const double domain[ 2 ][ 2 ] LonLat domain of subset.
OUTPUTS: Integer* timestamp            If not 0 and result is 1 then set it
                                       to the value parsed from the dataLine.
         double*   longitude           If not 0 and result is 1 then set it.
         double*   latitude            If not 0 and result is 1 then set it.
RETURNS: int 1 if dataLine is within the specified timestamps/domain, else 0.
******************************************************************************/

static int lineInSubset( char* dataLine,
                         Integer firstTimestamp, Integer lastTimestamp,
                         const double domain[ 2 ][ 2 ],
                         Integer* timestamp,
                         double* longitude, double* latitude ) {

  PRE8( dataLine, *dataLine,
        isValidTimestamp( firstTimestamp ), isValidTimestamp( lastTimestamp ),
        firstTimestamp <= lastTimestamp,
        domain, isValidDomain( domain ),
        IMPLIES_ELSE( timestamp,
                      AND2( longitude, latitude ),
                      AND2( ! longitude, ! latitude ) ) );

  int result = 0;
  char* word = skipWords( dataLine, 3 ); /* Skip STATION LONGITUDE LATITUDE */

  if ( word ) {
    const Integer dataTimestamp = atoI( word );

    if ( AND2( dataTimestamp > 0,
               IN_RANGE( dataTimestamp, firstTimestamp, lastTimestamp ) ) ) {
      double dataLongitude = 0.0;
      double dataLatitude = 0.0;
      const int ok =
        sscanf( dataLine, "%*s %lf %lf", &dataLatitude, &dataLongitude ) == 2;

      if ( AND3( ok,
                 IN_RANGE( dataLongitude,
                           domain[ LONGITUDE ][ MINIMUM ],
                           domain[ LONGITUDE ][ MAXIMUM ] ),
                 IN_RANGE( dataLatitude,
                           domain[ LATITUDE ][ MINIMUM ],
                           domain[ LATITUDE ][ MAXIMUM ] ) ) ) {

        if ( timestamp ) {
          *timestamp = dataTimestamp;
          *longitude = dataLongitude;
          *latitude  = dataLatitude;
        }

        result = 1;
      }
    }
  }

  POST2( IS_BOOL( result ),
         IMPLIES( AND2( result, timestamp ),
                  AND3( isValidTimestamp( *timestamp ),
                        IN_RANGE( *longitude,
                                  domain[ LONGITUDE ][ MINIMUM ],
                                  domain[ LONGITUDE ][ MAXIMUM ] ),
                        IN_RANGE( *latitude,
                                  domain[ LATITUDE ][ MINIMUM ],
                                  domain[ LATITUDE ][ MAXIMUM ] ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: uniqueStations - Create sorted array of unique stations in domain.
INPUTS:  Data*      data       lines, lineCount.
OUTPUTS: Data*      data       stations, stationCount.
                    ok         1 if successful, else 0.
******************************************************************************/

static void uniqueStations( Data* data ) {

  PRE6( data, data->ok, data->lines, data->lineCount > 0,
        data->stations == 0, data->stationCount == 0 );

  size_t line = 0;
  Integer previousStationId = 0;

  /*
   * The file data lines are grouped (but not sorted) by station
   * so they will be counted like uniq | sort | uniq:
   */

  for ( line = 0; line < data->lineCount; ++line ) {
    const Integer stationId = data->lines[ line ].station.id;

    if ( stationId != previousStationId ) {
      ++data->stationCount;
      previousStationId = stationId;
    }
  }

  CHECK ( data->stationCount );
  data->stations = NEW( Station, data->stationCount );

  if ( ! data->stations ) {
    data->stationCount = 0;
    data->ok = 0;
  } else { /* Now make a second pass and store the ids, lons, lats, notes: */
    size_t stationIndex = 0;
    size_t uniqueCount = 0;
    previousStationId = 0;

    for ( line = 0; line < data->lineCount; ++line ) {
      const Station* const station = &data->lines[ line ].station;
      const Integer stationId = station->id;
      CHECK( stationId > 0 );

      if ( stationId != previousStationId ) {
        CHECK( stationIndex < data->stationCount );
        data->stations[ stationIndex ].id        = stationId;
        data->stations[ stationIndex ].longitude = station->longitude;
        data->stations[ stationIndex ].latitude  = station->latitude;
        strncpy( data->stations[ stationIndex ].note, station->note,
                 NOTE_LENGTH + 1 );
        ++stationIndex;
        previousStationId = stationId;
      }
    }

    CHECK( stationIndex == data->stationCount );

    /* Sort the stations by id: */

    qsort( data->stations, data->stationCount, sizeof data->stations[ 0 ],
           stationComparer );

    /* Second uniq to avoid counting adjacent duplicates: */

    previousStationId = 0;

    for ( stationIndex = 0; stationIndex < data->stationCount;
          ++stationIndex ) {
      const Integer stationId = data->stations[ stationIndex ].id;

      if ( stationId != previousStationId ) {
        ++uniqueCount;
        previousStationId = stationId;
      }
    }

    /* Copy reduced array of unique stations: */

    if ( uniqueCount < data->stationCount ) {
      Station* uniqueStations = NEW( Station, uniqueCount );

      if ( uniqueStations ) {
        size_t uniqueStationIndex = 0;
        previousStationId = 0;

        for ( stationIndex = 0; stationIndex < data->stationCount;
              ++stationIndex ) {
          const Integer stationId = data->stations[ stationIndex ].id;

          if ( stationId != previousStationId ) {
            CHECK( uniqueStationIndex < uniqueCount );
            uniqueStations[ uniqueStationIndex ] =
              data->stations[ stationIndex ];
            ++uniqueStationIndex;
            previousStationId = stationId;
          }
        }

        data->stationCount = uniqueCount;
        FREE( data->stations );
        data->stations = uniqueStations;
        uniqueStations = 0;
      }
    }
  }

  POST( IMPLIES_ELSE( data->ok,
          AND8( data->stationCount > 0,
                data->stations,
                data->stations[ 0 ].id > 0,
                data->stations[ data->stationCount - 1 ].id > 0,
                IN_RANGE( data->stations[ 0 ].longitude, -180.0, 180.0 ),
                IN_RANGE( data->stations[ 0 ].latitude,   -90.0,  90.0 ),
             IN_RANGE( data->stations[ ( data->stationCount - 1 ) ].longitude,
                                        -180.0, 180.0 ),
             IN_RANGE( data->stations[ ( data->stationCount - 1 ) ].latitude,
                                        -90.0, 90.0 ) ),
          AND2( data->stationCount == 0, data->stations == 0 ) ) );
}



/******************************************************************************
PURPOSE: extractDataValues - Create array of data for subset.
INPUTS:  Data*      data       lines, lineCount, stations, stationCount
                               firstTimestamp, timesteps.
OUTPUTS: Data*      data       data[ timesteps ][ stations ].
                    ok         1 if successful, else 0.
                    timesteps  Divided by 24 if daily variable.
NOTES:   In cases where there is no data available for the given timestamp
         and/or station, missingValue will be stored.
******************************************************************************/

static void extractDataValues( Data* data ) {

  PRE9( data, data->ok, data->lines, data->lineCount > 0, data->stations,
        data->stationCount, data->firstTimestamp > 0, data->timesteps > 0,
        data->data == 0 );

  const int isWind = ! strcmp( data->variableName, "wind" );
  const int isDaily = strstr( data->variableName, "_daily" ) != 0;
  const size_t hoursPerTimestep = isDaily ? 24 : 1;
  const size_t timesteps0 = data->timesteps / hoursPerTimestep;
  const size_t timesteps  = timesteps0 > 0 ? timesteps0 : 1;
  data->timesteps = timesteps;

  {
    const size_t count = data->timesteps * data->stationCount;
    data->data = NEW( double, count + count * isWind );
    data->ok = data->data != 0;

    if ( data->data ) {
      double* value = data->data;
      double* value2 = isWind ? value + count : 0;
      size_t timestep = 0;
      Integer timestamp = data->firstTimestamp;

      do {
        size_t station = 0;

        do {
          const Integer stationId = data->stations[ station ].id;
          double v2 = -9999.0;
          *value++ = findValue( data, timestamp, stationId, &v2 );

          if ( isWind ) {
            *value2++ = v2;
          }

          ++station;
        } while ( station < data->stationCount );

        advanceTimestamp( hoursPerTimestep, &timestamp );
        ++timestep;
      } while ( timestep < data->timesteps );
    }
  }

  POST( data->ok == ( data->data != 0 ) );
}



/******************************************************************************
PURPOSE: writeHeader - Write ASCII header of subset to stdout.
INPUTS:  const Arguments*  arguments   description
         Data*             data        timesteps, stationCount, variableName.
******************************************************************************/

static void writeHeader( const Arguments* arguments, const Data* data ) {

  PRE3( isValidArguments( arguments ), isValidData( data ), data->ok );

  const int isWind = ! strcmp( data->variableName, "wind" );
  UTCTimestamp timestamp;
  toUTCTimestamp( data->firstTimestamp, timestamp );

  printf( "SITE 2.0\n" );
  printf( "%s\n", arguments->description );
  printf( "%s\n", timestamp );
  printf( "# data dimensions: timesteps stations\n" );
  printf( "%lu %lu\n", data->timesteps, data->stationCount );
  printf( "# Variable names:\n" );

  if ( isWind ) {
    printf( "wind_u wind_v\n" );
    printf( "# Variable units:\n" );
    printf( "m/s m/s\n" );
  } else {
    printf( "%s\n", data->variableName );
    printf( "# Variable units:\n" );
    printf( "%s\n", data->units );
  }

  printf( "# char notes[stations][80] and\n" );
  printf( "# MSB 64-bit integers ids[stations] and\n" );
  printf( "# IEEE-754 64-bit reals sites[stations][2=<longitude,latitude>]" );
  printf( " and\n");
  printf( "# IEEE-754 64-bit reals data[timesteps][stations]:\n");

  POST( data->ok );
}



/******************************************************************************
PURPOSE: writeXDR - Write XDR format output of subset to stdout.
INPUTS:  const Arguments*  arguments   description
         Data*             data        timesteps, stationCount, variableName,
                                       stations, data.
******************************************************************************/

static void writeXDR( const Arguments* arguments, Data* data ) {

  PRE3( isValidArguments( arguments ), isValidData( data ), data->ok );

  data->ok = 0;

  Integer* ids = NEW( Integer, data->stationCount );

  if ( ids ) {
    double* lonlats = NEW( double, data->stationCount * 2 );

    if ( lonlats ) {
      Integer* i = ids;
      double* ll = lonlats;
      const Station* stations = data->stations;
      size_t count = data->stationCount;
      data->ok = 1;
      writeHeader( arguments, data );

      do {
        *i++  = stations->id;
        *ll++ = stations->longitude;
        *ll++ = stations->latitude;
        printf( "%-79s\n", stations->note );
        ++stations;
      } while ( --count );

      count = data->stationCount;

#if IS_LITTLE_ENDIAN
      swapBytes8( ids, count );
#endif

      data->ok = fwrite( ids, sizeof *ids, count, stdout ) == count;

      if ( data->ok ) {
        count = data->stationCount * 2;

#if IS_LITTLE_ENDIAN
        swapBytes8( lonlats, count );
#endif

        data->ok = fwrite( lonlats, sizeof *lonlats, count, stdout ) == count;

        if ( data->ok ) {
          count = data->timesteps * data->stationCount;

          if ( ! strcmp( data->variableName, "wind" ) ) {
            count += count;
          }

#if IS_LITTLE_ENDIAN
          swapBytes8( data->data, count );
#endif

          data->ok =
            fwrite( data->data,  sizeof *data->data, count, stdout ) == count;

#if IS_LITTLE_ENDIAN
          swapBytes8( data->data, count );
#endif
        }
      }

      FREE( lonlats );
    }

    FREE( ids );
  }

  POST( IS_BOOL( data->ok ) );
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII spreadsheet (tab-separated text values)
         format output of subset to stdout.
INPUTS:  const Arguments*  arguments   description
         Data*             data        timesteps, stationCount, variableName,
                                       stations, data, firstTimestamp.
******************************************************************************/

static void writeASCII( const Arguments* arguments, const Data* data ) {

  PRE3( isValidArguments( arguments ), isValidData( data ), data->ok );

  const int isWind = ! strcmp( data->variableName, "wind" );
  const int isDaily = strstr( data->variableName, "_daily" ) != 0;
  const size_t hoursPerTimestep = isDaily ? 24 : 1;
  const double* values = data->data;
  const double* values2 =
    isWind ? values + data->stationCount * data->timesteps : 0;
  size_t timesteps    = data->timesteps;
  Integer yyyydddhhmm = data->firstTimestamp;
  const char* const headerStart =
    "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tSTATION(-)";
  const char* const headerFormat1 = "\t%s(%s)\tSITE_NAME\n";
  const char* const dataFormat1 =
    "%s\t%10.5f\t%10.5f\t%20lld\t%20.12e\t%44s\n";
  const char* const headerFormat2 = "\t%s(%s)\t%s(%s)\tSITE_NAME\n";
  const char* const dataFormat2 =
    "%s\t%10.5f\t%10.5f\t%20lld\t%20.12e\t%20.12e\t%44s\n";
  UTCTimestamp timestamp;
  toUTCTimestamp( yyyydddhhmm, timestamp );

  /* Write header row: */

  printf( "%s", headerStart );

  if ( isWind ) {
    printf( headerFormat2, "wind_u", "m/s", "wind_v", "m/s" );
  } else {
    printf( headerFormat1, data->variableName, data->units );
  }

  /* Write data rows: */

  do {
    size_t stationCount = data->stationCount;
    const Station* stations = data->stations;

    do {
      const Integer id       = stations->id;
      const double longitude = stations->longitude;
      const double latitude  = stations->latitude;
      const double value     = *values++;

      if ( isWind ) {
        const double value2 = *values2++;
        printf( dataFormat2,
                timestamp, longitude, latitude, id, value, value2,
                stations->note );

      } else {
        printf( dataFormat1,
                timestamp, longitude, latitude, id, value, stations->note );
      }
      ++stations;
    } while ( --stationCount );

    advanceTimestamp( hoursPerTimestep, &yyyydddhhmm );
    toUTCTimestamp( yyyydddhhmm, timestamp );
  } while ( --timesteps );

  POST( data->ok );
}



#ifndef NO_ASSERTIONS



/******************************************************************************
PURPOSE: isValidArguments - Is arguments in a valid state?
INPUTS:  const Arguments* arguments  Structure to examine.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int isValidArguments( const Arguments* arguments ) {
  int result =
    AND8( arguments,
          arguments->days > 0,
          arguments->fileNames,
          arguments->fileNames[ 0 ],
          arguments->fileNames[ arguments->days - 1 ],
          arguments->description,
          arguments->description[ 0 ],
          IS_VALID_OUTPUT_FORMAT( arguments->outputFormat ) );

  result = AND4( result,
                 isValidTimestamp( arguments->firstTimestamp ),
                 arguments->timesteps > 0,
                 isValidDomain( arguments->domain ) );

  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidData - Is data in a valid state?
INPUTS:  const Data* data  Structure to examine.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int isValidData( const Data* data ) {
  int result = data &&
    AND9( isValidTimestamp( data->firstTimestamp ),
          isValidTimestamp( data->lastTimestamp ),
          data->firstTimestamp <= data->lastTimestamp,
          data->timesteps > 0,
          data->lineCount > 0,
          data->stationCount > 0,
          IMPLIES( data->fileDataLength > 0, data->fileData ),
          data->variableName[ 0 ] && data->units[ 0 ],
          strlen( data->variableName ) <= NAME_LENGTH &&
          strlen( data->units )        <= NAME_LENGTH );

  result =
    AND9( result,
          IMPLIES( data->fileData, strlen( data->fileData ) < 100 ),
          data->lines,
          data->lines[ 0 ].station.id > 0,
          data->lines[ data->lineCount - 1 ].station.id > 0,
          IN_RANGE( data->lines[ 0 ].station.longitude, -180.0, 180.0 ),
          IN_RANGE( data->lines[ 0 ].station.latitude,   -90.0,  90.0 ),
          IN_RANGE( data->lines[ data->lineCount - 1 ].station.longitude,
                    -180.0, 180.0 ),
          IN_RANGE( data->lines[ data->lineCount - 1 ].station.latitude,
                     -90.0,  90.0 ) );

  result =
    AND9( result,
          data->stations,
          data->stations[ 0 ].id > 0,
          data->stations[ data->stationCount - 1 ].id > 0,
          IN_RANGE( data->stations[ 0 ].longitude, -180.0, 180.0 ),
          IN_RANGE( data->stations[ 0 ].latitude,   -90.0,  90.0 ),
          IN_RANGE( data->stations[ data->stationCount - 1 ].longitude,
                    -180.0, 180.0 ),
          IN_RANGE( data->stations[ data->stationCount - 1 ].latitude,
                     -90.0,  90.0 ),
          IS_BOOL( data->ok ) );

  result =
    AND7( result,
          data->stations[ 0 ].note,
          data->stations[ 0 ].note[ 0 ],
          data->stations[ 0 ].note[ NOTE_LENGTH ] == '\0',
          data->stations[ data->stationCount - 1 ].note,
          data->stations[ data->stationCount - 1 ].note[ 0 ],
          data->stations[ data->stationCount - 1 ].note[ NOTE_LENGTH] == '\0');

  POST( IS_BOOL( result ) );
  return result;

}



#endif /* ! defined( NO_ASSERTIONS ) */



/******************************************************************************
PURPOSE: isValidDomain - Is domain a valid longitude, latitude range?
INPUTS:  const double* domain[ 2 ][ 2 ]  Domain to examine.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int isValidDomain( const double domain[ 2 ][ 2 ] ) {
  PRE( domain );
  const int result =
    AND4( IN_RANGE( domain[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ),
          IN_RANGE( domain[ LONGITUDE ][ MAXIMUM ],
                    domain[ LONGITUDE ][ MINIMUM ], 180.0 ),
          IN_RANGE( domain[ LATITUDE  ][ MINIMUM ], -90.0, 90.0 ),
          IN_RANGE( domain[ LATITUDE  ][ MAXIMUM ],
                    domain[ LATITUDE  ][ MINIMUM ], 90.0 ) );
  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: stationComparer - Compare station ids and return -1, 0 or 1 when
         a < b or a == b or a > b.
INPUTS:  const void* a    First  station to compare.
         const void* a    Second station to compare.
RETURNS: int -1, 0 or 1 when a < b or a == b or a > b.
******************************************************************************/

static int stationComparer( const void* a, const void* b ) {
  PRE2( a, b );
  const Station* station1 = a;
  const Station* station2 = b;
  const Integer id1 = station1->id;
  const Integer id2 = station2->id;
  const int result = id1 < id2 ? -1 : id1 > id2 ? 1 : 0;
  CHECK2( id1 > 0, id2 > 0 );
  POST( IN4( result, -1, 0, 1 ) );
  return result;
}


/******************************************************************************
PURPOSE: lineComparer - Compare line timestamps/station ids and return -1 or 1
         when a < b or a > b.
INPUTS:  const void* a  First  line to compare.
         const void* b  Second line to compare.
RETURNS: int 1 if valid, else 0.
NOTES:   Duplicate lines are trapped as a defect.
         Lines look like:
         0123456789012345678901234567890123456789012345678901234567890
         000010301   48.9494  -57.9458 20041830330        3.0 ...
         350130021   31.7961 -106.5839 20042441730       -1.0 ...
         and comparisions are made on the timestamp (4th word) and then
         the stationId (1st word).
******************************************************************************/

static int lineComparer( const void* a, const void* b ) {
  PRE2( a, b );
  const Line* const line1 = a;
  const Line* const line2 = b;
  const Integer timestamp1 = line1->timestamp;
  const Integer timestamp2 = line2->timestamp;
  int result = 0;

  CHECK2( isValidTimestamp( timestamp1 ), isValidTimestamp( timestamp2 ) );

  if ( timestamp1 < timestamp2 ) {
    result = -1;
  } else if ( timestamp1 > timestamp2 ) {
    result = 1;
  } else { /* Equal timestamps so compare stationIds: */
    const Integer stationId1 = line1->station.id;
    const Integer stationId2 = line2->station.id;
    CHECK3( stationId1 > 0, stationId2 > 0, stationId1 != stationId2 );

    if ( stationId1 < stationId2 ) {
      result = -1;
    } else {
      result = 1;
    }
  }

  POST( OR2( result == -1, result == 1 ) ); /* No line/data duplicates! */

  return result;
}



/******************************************************************************
PURPOSE: findValue - Find the data value for the given timestamp and station.
INPUTS:  const Data* data     lines to search.
         Integer yyyydddhhmm  Timestamp to search for.
         int stationId        Station id to search for.
OUTPUTS: double* value2       2nd value if wind.
RETURNS: double found value or missingValue if not found.
******************************************************************************/

static double findValue( const Data* data, Integer yyyydddhhmm,
                         Integer stationId, double* value2 ) {

  PRE4( isValidData( data ), isValidTimestamp( yyyydddhhmm ), stationId > 0,
        value2 );

  const int isWind = ! strcmp( data->variableName, "wind" );
  const Line* const lines = data->lines;
  double result = missingValue;
  size_t index = data->lineCount;
  size_t lower = 0;
  size_t upper = data->lineCount;
  *value2 = missingValue;

  /* Binary search to find timestamp and stationId: */

  while ( lower < upper ) {
    const size_t middle = lower + ( upper - lower ) / 2;
    CHECK( IN_RANGE( middle, 0, data->lineCount - 1 ) );
    const Integer timestamp = lines[ middle ].timestamp;
    const Integer station   = lines[ middle ].station.id;

    if ( OR2( yyyydddhhmm < timestamp,
              AND2( yyyydddhhmm == timestamp, stationId < station ) ) ) {
      upper = middle;
    } else if ( OR2( yyyydddhhmm > timestamp,
                     AND2( yyyydddhhmm == timestamp, stationId > station ))) {
      lower = middle + 1;
    } else {
      CHECK2( yyyydddhhmm == timestamp, stationId == station );
      index = middle;
      lower = middle;
      upper = lower;
    }
  }

  if ( index != data->lineCount ) { /* Found value: */
    result = lines[ index ].value;

    if ( isWind ) {
      *value2 = lines[ index ].value2;
    }
  }

  if ( AND2( result < -98.0, strstr( data->variableName, "fire_" ) ) ) {
    result = missingValue;
  }

  return result;
}



/******************************************************************************
PURPOSE: allocate - Allocate and zero bytes of memory.
INPUTS:  size_t bytes Number of bytes to allocate.
RETURNS: void* address to sequence of zeroed bytes of memory, else 0 if failed.
******************************************************************************/

static void* allocate( size_t bytes ) {
  PRE( bytes > 0 );
  void* result = 0;
  int forceFailure = 0;
  CHECK( failureCountDown >= 0 );

  if ( failureCountDown > 0 ) { /* If we're counting down then... */
    --failureCountDown;
    forceFailure = failureCountDown == 0; /* Fail iff 0 is reached. */
  }

  if ( ! forceFailure ) {
    result = malloc( bytes );
  }

  if ( ! result ) {
    fprintf( stderr, "\a\n\nI'm sorry, can't allocate %lu bytes ", bytes );
    fprintf( stderr, "of memory to complete the requested action.\n" );
    perror( 0 );
  } else {
    memset( result, 0, bytes );
  }

  return result;
}



/******************************************************************************
PURPOSE: fileSize - Determine size of named file.
INPUTS:  const char* name  Name of file to examine.
RETURNS: size_t Size, in bytes, of named file, else 0 if failed.
******************************************************************************/

static size_t fileSize( const char* name ) {
  PRE2( name, *name );
  size_t result = 0;
  struct stat buf;

  if ( stat( name, &buf ) == -1 ) {
    fprintf( stderr, "\a\n\nFailed to determine size of file '%s'.\n", name );
    perror( 0 );
  } else {
    result = buf.st_size;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFile - Read named file into memory and return it as an allocated
         string.
INPUTS:  const char* name  Name of file to examine.
OUTPUTS: size_t* length    Length of string.
RETURNS: char* string contents of file, else 0 if failed.
******************************************************************************/

static char* readFile( const char* name, size_t* length ) {
  PRE3( name, *name, length );
  char* result = 0;
  *length = fileSize( name ) / sizeof (char);

  if ( *length > 0 ) {
    result = NEW( char, *length + 1 );

    if ( ! result ) {
      *length = 0;
    } else {
      FILE* file = fopen( name, "rb" );

      if ( file ) {
        const size_t itemsRead =
          fread( result, *length * sizeof (char), 1, file );

        if ( itemsRead != 1 ) {
          fprintf( stderr, "\a\n\nFailed to read entire file '%s'.\n", name );
          perror( 0 );
          FREE( result );
          *length = 0;
        } else {
          result[ *length ] = '\0'; /* Terminate string. */
        }

        fclose( file );
        file = 0;
      }
    }
  }

  POST( IMPLIES_ELSE( result, *length > 0, *length == 0 ) );
  return result;
}



#if IS_LITTLE_ENDIAN


#if 0

/******************************************************************************
PURPOSE: swapBytes4 - Swap byte order of each 4-byte word in the array.
INPUTS:  void*  array  Array of 4-byte words to swap.
         size_t count  Number of 4-byte words in array.
INPUTS:  void*  array  Array with bytes swapped.
******************************************************************************/

static void swapBytes4( void* array, size_t count ) {
  assert_static( sizeof (unsigned int) == 4 );
  CHECK( count > 0 );
  {
    unsigned int* word = array;

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

#endif


/******************************************************************************
PURPOSE: swapBytes8 - Swap byte order of each 8-byte word in the array.
INPUTS:  void*  array  Array of 8-byte words to swap.
         size_t count  Number of 8-byte words in array.
INPUTS:  void*  array  Array with bytes swapped.
******************************************************************************/

static void swapBytes8( void* array, size_t count ) {
  assert_static( sizeof (unsigned long long) == 8 );
  CHECK( count > 0 );
  {
    unsigned long long* word = array;

    do {
      const unsigned long long value = *word;
      const unsigned long long swapped =
        ( value & 0xff00000000000000LL ) >> 56 |
        ( value & 0x00ff000000000000LL ) >> 40 |
        ( value & 0x0000ff0000000000LL ) >> 24 |
        ( value & 0x000000ff00000000LL ) >>  8 |
        ( value & 0x00000000ff000000LL ) <<  8 |
        ( value & 0x0000000000ff0000LL ) << 24 |
        ( value & 0x000000000000ff00LL ) << 40 |
        ( value & 0x00000000000000ffLL ) << 56;
      *word = swapped;
      ++word;
      --count;
    } while ( count );
  }
}



#endif /* IS_LITTLE_ENDIAN. */



/******************************************************************************
PURPOSE: parseTimestamp - Parse string timestamp into its integer value.
INPUTS:  const char* string  String representation of timestamp.
OUTPUTS: Integer*    value   Value of string as yyyydddhh00.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static Integer parseTimestamp( const char* string, Integer* value ) {

  PRE2( string, value );

  const Integer yyyymmddhh = atoI( string );
  const int yyyymmdd       = yyyymmddhh / 100;
  const int hh             = yyyymmddhh % 100;
  const int result = AND2( IN_RANGE( hh, 0, 23 ), isValidDate( yyyymmdd ) );

  *value = 0;

  if ( ! result ) {
    fprintf( stderr, "\a\n\nInvalid timestamp specified '%s'.\n", string );
  } else {
    const Integer yyyyddd   = convertDate( yyyymmdd );
    const Integer yyyydddhh = yyyyddd * 100 + hh;
    *value = yyyydddhh * 100;
  }

  POST2( IS_BOOL( result ),
         IMPLIES_ELSE( result, isValidTimestamp( *value ), *value == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: toUTCTimestamp - Convert timestamp to ISO UTC string format.
INPUTS:  Integer yyyydddhhmm  Timestamp to convert.
OUTPUTS: UTCTimestamp string  ISO UTC string format for the timestamp.
******************************************************************************/

static void toUTCTimestamp( Integer yyyydddhhmm, UTCTimestamp string ) {
  PRE2( isValidTimestamp( yyyydddhhmm ), string );
  const int mm      = yyyydddhhmm % 100;
  const int hh      = yyyydddhhmm / 100 % 100;
  const int yyyyddd = yyyydddhhmm / 10000;
  const int yyyy    = yyyyddd / 1000;
  int mo = 0;
  int dd = 0;
  monthDay( yyyyddd, &mo, &dd );
  sprintf( string, "%04d-%02d-%02dT%02d:%02d:00-0000",
           yyyy, mo, dd, hh, mm );
  POST( strlen( string ) == UTC_TIMESTAMP_LENGTH );
}



#ifndef NO_ASSERTIONS



/******************************************************************************
PURPOSE: isValidTimestamp - Is the timestamp valid?
INPUTS:  Integer yyyydddhhmm The timestamp to examine.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int isValidTimestamp( Integer yyyydddhhmm ) {
  const int yyyy = yyyydddhhmm / 10000000;
  const int ddd  = yyyydddhhmm / 10000 % 1000;
  const int hh   = yyyydddhhmm / 100 % 100;
  const int mm   = yyyydddhhmm % 100;
  const int result =
    AND4( IN_RANGE( yyyy, 1900, 9999 ),
          IN_RANGE( ddd, 1, 365 + isLeapYear( yyyy ) ),
          IN_RANGE( hh, 0, 23 ),
          IN_RANGE( mm, 0, 59 ) );
  POST( IS_BOOL( result ) );
  return result;
}



#endif /* ! defined( NO_ASSERTIONS ) */



/******************************************************************************
PURPOSE: isValidDate - Is the date valid?
INPUTS:  Integer yyyymmdd The date to examine.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int isValidDate( Integer yyyymmdd ) {
  const int yyyy = yyyymmdd / 10000;
  const int mm   = yyyymmdd / 100 % 100;
  const int dd   = yyyymmdd % 100;
  const int result =
    AND3( IN_RANGE( yyyy, 1900, 9999 ),
          IN_RANGE( mm, 1, 12 ),
          IN_RANGE( dd, 1, daysPerMonth[ isLeapYear( yyyy ) ][ mm - 1 ] ) );
  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isLeapYear - Is yyyy a leap year (i.e., February has 29 days)?
INPUTS:  int yyyy  Year to examine.
RETURNS: int 1 if yyyy is a leap year, else 0.
******************************************************************************/

static int isLeapYear( int yyyy ) {
  PRE( yyyy > 1000 );
  const int result = yyyy % 4 == 0 && ( yyyy % 100 != 0 || yyyy % 400 == 0 );
  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: convertDate - Convert date from YYYYMMDD to YYYYDDD.
INPUTS:  int yyyymmdd  Year, month, day to convert.
RETURNS: int yyyyddd equivalent day.
******************************************************************************/

static int convertDate( int yyyymmdd ) {

  PRE3( yyyymmdd / 10000 > 1000,
        IN_RANGE(  yyyymmdd / 100 % 100, 1, 12 ),
        IN_RANGE(  yyyymmdd % 100, 1,
                   daysPerMonth[ isLeapYear( yyyymmdd / 10000 ) ]
                   [ yyyymmdd / 100 % 100 - 1] ) );

  const int yyyy = yyyymmdd / 10000;
  const int mm0  = yyyymmdd / 100 % 100 - 1;
  const int dd   = yyyymmdd % 100;
  const int leap = isLeapYear( yyyy );
  int result = 0;
  int month = 0;

  for ( month = 0; month < mm0; ++month ) {
    result += daysPerMonth[ leap ][ month ];
  }

  result += dd;
  result += yyyy * 1000;

  POST3( result / 1000 == yyyymmdd / 10000,
         result % 1000 > 0,
         result % 1000 <= 365 + isLeapYear( yyyymmdd / 10000 ) );
  return result;
}



/******************************************************************************
PURPOSE: monthDay - Extract month and day of month from YYYYDDD.
INPUTS:  int yyyyddd  Year and day of year to convert.
OUTPUTS: int* mm      Month [1, 12].
         int* dd      Day of month [1, daysPerMonth[ leap ][ mm - 1 ] ].
******************************************************************************/

static void monthDay( int yyyyddd, int* mm, int* dd ) {

  PRE4( yyyyddd / 1000 > 1000,
        IN_RANGE(  yyyyddd % 1000, 1, 365 + isLeapYear( yyyyddd / 1000 ) ),
        mm, dd );

  const int yyyy = yyyyddd / 1000;
  const int leap = isLeapYear( yyyy );
  int ddd  = yyyyddd % 1000;
  int month = 0;

  for ( month = 0;
        AND2( month < 12, ddd > daysPerMonth[ leap ][ month ] );
        ++month ) {
    ddd -= daysPerMonth[ leap ][ month ];
  }

  *mm = month + 1;
  *dd = ddd;
  POST2( IN_RANGE( *mm, 1, 12 ),
         IN_RANGE( *dd, 1,
                   daysPerMonth[ isLeapYear( yyyyddd / 1000 ) ][ *mm - 1 ] ));
}



/******************************************************************************
PURPOSE: advanceTimestamp - Advance timestamp by specified hours.
INPUTS:  const Integer  hours  Number of hours to advance timestamp by.
         Integer* yyyydddhhmm  Timestamp to increment.
OUTPUTS: Integer* yyyydddhhmm  Inremented timestamp.
******************************************************************************/

static void advanceTimestamp( const Integer hours, Integer* yyyydddhhmm ) {
  PRE3( hours > 0, yyyydddhhmm, isValidTimestamp( *yyyydddhhmm ) );
  Integer hour = 0;

  for ( hour = 0; hour < hours; ++hour ) {
    incrementTimestamp( yyyydddhhmm );
  }

  POST( isValidTimestamp( *yyyydddhhmm ) );
}



/******************************************************************************
PURPOSE: incrementTimestamp - Increment timestamp by one hour.
INPUTS:  Integer* yyyydddhhmm  Timestamp to increment.
OUTPUTS: Integer* yyyydddhhmm  Inremented timestamp.
******************************************************************************/

static void incrementTimestamp( Integer* yyyydddhhmm ) {
  PRE2( yyyydddhhmm, isValidTimestamp( *yyyydddhhmm ) );
  const Integer mm = *yyyydddhhmm % 100;
  Integer hh = *yyyydddhhmm / 100 % 100;
  ++hh;

  if ( hh < 24 ) { /* Just update the hour: */
    *yyyydddhhmm = *yyyydddhhmm / 10000 * 10000 + hh * 100 + mm;
  } else { /* Reset the hour and increment the day: */
    Integer yyyy = *yyyydddhhmm / 10000000;
    Integer ddd  = *yyyydddhhmm / 10000 % 1000;
    hh = 0;
    ++ddd;

    if ( AND2( ddd > 365, ddd > 365 + isLeapYear( yyyy ) ) ) {
      ddd = 1; /* First day,    */
      ++yyyy;  /* of next year. */
    }

    *yyyydddhhmm = yyyy * 10000000 + ddd * 10000 + hh * 00 + mm;
  }

  POST( isValidTimestamp( *yyyydddhhmm ) );
}



/******************************************************************************
PURPOSE: commasToSpaces - Change commas (and control-Ms) to spaces.
INPUTS:  char* string  String to modify.
OUTPUTS: char* string  String without commas.
RETURNS: size_t line count.
******************************************************************************/

static size_t commasToSpaces( char* string ) {
  PRE( string );
  size_t result = 0;
  char* s = string;

  while ( *s ) {
    const char c = *s;

    if ( c == ',' || c == '\r' ) {
      *s = ' ';
    } else if ( c == '\n' ) {
      ++result;
    }

    ++s;
  }

  return result;
}



/******************************************************************************
PURPOSE: indexOfString - Index of string in strings[] or -1 if not present.
INPUTS:  const char* string           String to search for.
         const char* const strings[]  Strings to search.
         int count                    Size of strings[].
RETURNS: int index of string in strings[], else -1 if not present.
******************************************************************************/

static int indexOfString( const char* string, const char* const strings[],
                          int count ) {

  PRE8( string, *string, strings, strings[ 0 ], *strings[ 0 ],
        count > 0, strings[ count - 1 ], *strings[ count - 1 ] );

  int result = -1;
  int index = 0;

  do {

    if ( strcmp( string, strings[ index ] ) == 0 ) {
      result = index;
      index = count;
    }

    ++index;
  } while ( index < count );

  POST( OR2( result == -1,
             AND2( result >= 0,
                   strcmp( string, strings[ result ] ) == 0 ) ) );

  return result;
}



/******************************************************************************
PURPOSE: linesMatch - Do strings match up to end of line/string?
INPUTS:  const char* line1  1st string to compare.
         const char* line2  2nd string to compare.
RETURNS: int 1 if matched, else 0.
******************************************************************************/

static int linesMatch( const char* line1, const char* line2 ) {
  PRE2( line1, line2 );
  int result = line1 == line2;

  if ( ! result ) {

    while ( AND4( *line1, *line2, *line1 != '\n', *line2 != '\n' ) ) {
      ++line1;
      ++line2;
    }

    result = OR2( AND2( *line1 == '\n', *line2 == '\n' ),
                  AND2( *line1 == '\0', *line2 == '\0' ) );
  }

  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: skipLine - Skip to end of line.
INPUTS:  const char* line     String to parse.
OUTPUTS: size_t*     skipped  Number of characters skipped, including '\n'.
RETURNS: const char* pointer into line after first '\n' or else '\0'.
******************************************************************************/

static const char* skipLine( const char* line, size_t* skipped ) {
  PRE2( line, skipped );
  const char* result = line;
  *skipped = 0;

  while ( AND2( *result, *result != '\n' ) ) {
    ++*skipped;
    ++result;
  }

  if ( *result == '\n' ) {
    ++*skipped;
    ++result;
  }

  POST( result );
  return result;
}



/******************************************************************************
PURPOSE: skipWords - Skip over a given number of words.
INPUTS:  char* string  String to scan.
         Integer count       Number of words to skip.
RETURNS: const char*  if skipped, else 0.
******************************************************************************/

static char* skipWords( char* string, int count ) {

  PRE2( string, count > 0 );

  char* result = string;
  int counter = 0;

  do {

    /* Skip any blanks before a possible word: */

    while ( isspace( *result ) ) {
      ++result;
    }

    if ( *result ) {
      ++counter; /* Count the word. */

      /* Skip the counted word: */

      do {
        ++result;
      } while ( AND2( *result, ! isspace( *result ) ) );
    }

  } while ( AND2( counter < count, *result ) );

  if ( *result ) {

    /* Skip any blanks after the last word: */

    while ( isspace( *result ) ) {
      ++result;
    }
  }

  if ( OR2( *result == '\0', counter != count ) ) {
    result = 0;
  }

  return result;
}



