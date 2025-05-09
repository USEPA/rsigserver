
/******************************************************************************
PURPOSE: http_connection.c - Routines for opening HTTP GET URL socket
         connections and reading ASCII and binary data from them.
NOTES:   Uses curl library to create sockets. Does not support HTTPS.
         Compile with -DNO_LIBCURL=1 to disable (replace with failing stubs)
         initialize/open/close routines but still allow
         read_http_connection_line() and read_http_connection_array()
         to be used with a given stream.
         http://curl.haxx.se/libcurl/c/curl_easy_setopt.html
HISTORY: 2012-06-06 plessel.todd@epa.gov
******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For FILE, stderr, fprintf(), snprintf(), fdopen().*/
#include <stdlib.h> /* For atexit(), malloc(), free(). */
#include <string.h> /* For memset(), strncat(), strchr(). */
#include <ctype.h>  /* For isprint(), isspace(). */

#ifndef _WIN32
#include <fcntl.h>      /* For O_NONBLOCK, F_GETFL, F_SETFL, fcntl(). */
#include <sys/select.h> /* For fd_set, struct timeval, FD_SET(), select(). */
#endif

#ifdef _WIN32
/* To try to use libcurl, comment-out the next line: */
#define NO_LIBCURL 1
#define CURL_STATICLIB 1
/* Disable for now due to undefined reference to _imp__curl_global_init, etc.*/
#endif

#ifndef NO_LIBCURL
#include <curl/curl.h>  /* For CURL*, curl_*(). */

/* Declare 'internal' CURL routine to use in lieu of select(): */

extern int Curl_socket_check( curl_socket_t readfd0,
                              curl_socket_t readfd1,
                              curl_socket_t writefd,
                              long timeout_ms );
#else

/* Define stubs. Client must use popen(curl.exe): */

typedef void CURL;
typedef int CURLcode;
typedef int curl_socket_t;
static void curl_global_cleanup( void ) {} /* Passed to atexit(). */
enum {
  CURLE_OK,
  CURL_GLOBAL_ALL, CURLINFO_LASTSOCKET, CURLOPT_URL, CURLOPT_CONNECT_ONLY,
  CURLOPT_TCP_KEEPALIVE, CURLOPT_VERBOSE,
  CURL_SOCKET_BAD, CURL_CSELECT_IN, CURL_CSELECT_OUT
};
#define curl_global_init( unused_ ) (-1)
#define curl_easy_init( unused_ ) 0
static void curl_easy_cleanup( CURL* unused_ ) {}
#define curl_easy_setopt( unused_1, unused_2, unused_3 ) (-1)
#define curl_easy_perform( unused_ ) (-1)
#define curl_easy_strerror( unused_ ) \
  "UNIMPLEMENTED: libcurl cannot be used with _WIN32. " \
  "Use CYGWIN (omit -mno-cygwin) or use popen( curl.exe ) instead."
#define curl_easy_send( unused_1, unused_2, unused_3, unused_4 ) (-1)
#define curl_easy_getinfo( unused_1, unused_2, unused_3 ) (-1)
#define Curl_socket_check( unused_1, unused_2, unused_3, unused_4 ) 0
#endif

#include <Assertions.h>      /* For PRE*(), POST*(), CHECK*(), DEBUG(). */
#include <Failure.h>         /* For failureMessage(). */
#include <http_connection.h> /* For public interface. */

/*================================== MACROS =================================*/

/*
 * Is the platform big-endian (MSB: most significant byte first) or
 * little-endian (LSB: least significant byte first)?
 */

#ifndef IS_LITTLE_ENDIAN
#if ( \
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
#define IS_LITTLE_ENDIAN 1
#else
#define IS_LITTLE_ENDIAN 0
#endif
#endif

/*================================== TYPES ==================================*/

enum { READ_MODE, WRITE_MODE };

/*========================== FORWARD DECLARATIONS ===========================*/

static int initialized = 0;
static int set_options( CURL* curl, const char* url );
static const char* skip_hostname( const char* url );
static int send_http_get_request( CURL* curl, const char* url );
static FILE* skip_http_header( CURL* curl );
static int set_socket_timeout( CURL* curl, int mode, int seconds );
static int set_socket_blocking( CURL* curl );
static int wait_on_socket(curl_socket_t socket, int mode, int timeout_seconds);
static void reverse_2byte_words_if_little_endian( size_t count, void* array );
static void reverse_4byte_words_if_little_endian( size_t count, void* array );
static void reverse_8byte_words_if_little_endian( size_t count, void* array );
static char* encode_spaces_and_percents( const char* string );
static int is_text( const char* string );

/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: initialize_http_connections - Initialize state before connecting.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES: Call once per program execution (i.e., this routine is not thread safe).
******************************************************************************/

int initialize_http_connections( void ) {
  PRE0( ! is_initialized_http_connections() );
  const CURLcode error = curl_global_init( CURL_GLOBAL_ALL );
  const int result = ! error;

  if ( error ) {
    failureMessage( "Failed because %s.", curl_easy_strerror( error ) );
  } else {
    initialized = 1;
    atexit( curl_global_cleanup );
  }

  POST0( result == is_initialized_http_connections() );
  return result;
}



/******************************************************************************
PURPOSE: is_initialized_http_connections - Called initialize_http_connections?
RETURNS: int 1 if initialized_http_connections() was called.
******************************************************************************/

int is_initialized_http_connections( void ) {
  const int result = initialized;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: open_http_connection - Open a HTTP/S GET URL connection to read from.
INPUTS:  const char* url Full HTTP/S URL.
         int timeout     Maximum number of seconds to wait to receive data
                         before giving up in failure or 0 to wait indefinitely.
OUTPUTS: FILE** stream   Stream to read contents from.
RETURNS: void* connection object to close when done reading from stream.
******************************************************************************/

void* open_http_connection( const char* url, int timeout, FILE** stream ) {
  PRE06( is_initialized_http_connections(),
         url, *url,
         OR2( strstr( url, "http://"  ) == url,
              strstr( url, "https://" ) == url ),
         timeout >= 0, stream );
  CURL* curl = curl_easy_init();
  *stream = 0;

  if ( curl ) {

    if ( set_options( curl, url ) ) {

      if ( set_socket_timeout( curl, WRITE_MODE, timeout ) ) {
        char* encoded_url = encode_spaces_and_percents( url );

        if ( encoded_url ) {

          if ( send_http_get_request( curl, encoded_url ) ) {

            if ( set_socket_timeout( curl, READ_MODE, timeout ) ) {

              if ( set_socket_blocking( curl ) ) {
                *stream = skip_http_header( curl );
              }
            }
          }

          free( encoded_url ), encoded_url = 0;
        }
      }
    }
  }

  if ( ! *stream ) {

    if ( curl ) {
      curl_easy_cleanup( curl ), curl = 0;
    }
  }

  POST0( is_initialized_http_connections() );
  return curl;
}



/******************************************************************************
PURPOSE: close_http_connection - Close connection.
INPUTS:  void* http_connection  Connection to close.
OUTPUTS: void* http_connection  Closed connection.
******************************************************************************/

void close_http_connection( void* http_connection ) {
  PRE0( is_initialized_http_connections() );

  if ( http_connection ) {
    curl_easy_cleanup( http_connection ), http_connection = 0;
  }

  POST02( is_initialized_http_connections(), http_connection == 0 );
}



/******************************************************************************
PURPOSE: read_http_connection_line - Read a line from a stream associated with
         an http connection.
INPUTS:  FILE* stream       Connection stream to read from.
         size_t size        Size of line.
OUTPUTS: char line[ size ]  Line read from stream.
RETURNS: int 1 if line read, else 0.
******************************************************************************/

int read_http_connection_line( FILE* stream, size_t size, char line[] ) {
  PRE03( stream, size, line );
  int result = 0;
  size_t bytes_read = 0;
  memset( line, 0, size * sizeof *line );

  do {
    fgets( line + bytes_read, size - bytes_read, stream );
    /* DEBUG( fprintf( stderr, "line = '%s'\n", line ); ) */
    bytes_read = strlen( line );
    CHECK( bytes_read < size );
  } while ( AND3( ! feof( stream ),
                  bytes_read < size - 1,
                  IMPLIES( bytes_read, line[ bytes_read - 1 ] != '\n' ) ) );

  if ( bytes_read > 0 ) {

    /* Replace '\r' (DOS control-M) character with ' ': */

    char* const c = strchr( line, '\r' );

    if ( c ) {
      *c = ' ';
    }

    result = is_text( line ); /* Check that only ASCII text characters. */
  }

  POST0( result == AND3( *line != '\0', !
                         strchr( line, '\r' ),
                         is_text( line ) ) );
  return result;

}



/******************************************************************************
PURPOSE: read_http_connection_array - Read an array of big-endian binary data
         from a stream associated with an http connection.
INPUTS:  FILE* stream      Connection stream to read from.
         size_t count      Number of words to read.
         size_t word_size  Bytes per word.
OUTPUTS: void* array       Data words read from stream.
RETURNS: int 1 if count words are read, else 0.
         if word_size is 2, 4 or 8 then word bytes will have been reversed on a
         little-endian platform.
******************************************************************************/

int read_http_connection_array( FILE* stream, size_t count, size_t word_size,
                                void* array ) {
  PRE04( stream, count, IN5( word_size, 1, 2, 4, 8 ), array );
  int result = 0;
  char* data = (char*) array;
  const size_t total_bytes_to_read = count * word_size;
  size_t total_bytes_read = 0;
  memset( data, 0, total_bytes_to_read );
  DEBUG( fprintf( stderr, "total_bytes_to_read = %lu\n", total_bytes_to_read);)

  do {
    const size_t bytes_remaining = total_bytes_to_read - total_bytes_read;
    const size_t bytes_read =
      fread( data + total_bytes_read, 1, bytes_remaining, stream );
    total_bytes_read += bytes_read;
    DEBUG( fprintf( stderr, "bytes_read = %lu\n", bytes_read ); )
    CHECK( total_bytes_read <= total_bytes_to_read );
  } while ( AND2( ! feof( stream ), total_bytes_read < total_bytes_to_read ) );

  result = total_bytes_read == total_bytes_to_read;

  if ( word_size == 2 ) {
    reverse_2byte_words_if_little_endian( count, array );
  } else if ( word_size == 4 ) {
    reverse_4byte_words_if_little_endian( count, array );
  } else if ( word_size == 8 ) {
    reverse_8byte_words_if_little_endian( count, array );
  }

  POST0( IS_BOOL( result ) );
  return result;

}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: set_options - Set configuration options before sending request.
INPUTS:  CURL* curl       Connection to use.
         const char* url  Full HTTP/S URL.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int set_options( CURL* curl, const char* url ) {

  PRE05( is_initialized_http_connections(),
         curl, url, *url,
         OR2( strstr( url, "http://"  ) == url,
              strstr( url, "https://" ) == url ) );

  int result = 0;
  CURLcode error = curl_easy_setopt( curl, CURLOPT_URL, url );

  if ( error == CURLE_OK ) {
    error = curl_easy_setopt( curl, CURLOPT_CONNECT_ONLY, 1 );

    if ( error == CURLE_OK ) {
      error = curl_easy_setopt( curl, CURLOPT_TCP_KEEPALIVE, 1 );
      DEBUG( fprintf( stderr, "KEEPALIVE = %d\n", error ); )
      error = CURLE_OK; /* Ignore error. */

      if ( error == CURLE_OK ) {
        error = curl_easy_setopt( curl, CURLOPT_TCP_KEEPINTVL, 60 );
        DEBUG( fprintf( stderr, "KEEPINTVL = %d\n", error ); )
        error = CURLE_OK; /* Ignore error. */

        if ( error == CURLE_OK ) {
          error = curl_easy_setopt( curl, CURLOPT_FOLLOWLOCATION, 1 );
          DEBUG( fprintf( stderr, "FOLLOWLOCATION = %d\n", error ); )
          error = CURLE_OK; /* Ignore error. */

          if ( error == CURLE_OK ) {
            DEBUG( error = curl_easy_setopt( curl, CURLOPT_VERBOSE, 1 ); )
            error = curl_easy_perform( curl );

            if ( error != CURLE_OK ) {

              /* Don't report failures to connect to internal host: */

              if ( ! AND2( strstr( url, "rtpmeta" ),
                           OR2( strstr( url, "REQUEST=GetVersion" ),
                                strstr( url, "TEST=1" ) ) ) ) {
                failureMessage( "Failed to connect to '%s'\n%s.\n",
                                url, curl_easy_strerror( error ) );
              }

            } else {
              result = 1;
            }
          }
        }
      }
    }
  }

  POST0( is_initialized_http_connections() );
  return result;
}



/******************************************************************************
PURPOSE: send_http_get_request - Send HTTP/S GET request.
INPUTS:  CURL* curl       Connection to use.
         const char* url  Full HTTP/S URL.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int send_http_get_request( CURL* curl, const char* url ) {

  PRE05( is_initialized_http_connections(),
         curl, url, *url,
         OR2( strstr( url, "http://"  ) == url,
              strstr( url, "https://" ) == url ) );

  int result = 0;
  const char* const query = skip_hostname( url );

  if ( query ) {
    enum { HTTP_GET_LENGTH = 4095 };
    char http_get_request[ HTTP_GET_LENGTH + 1 ] = "";
    const int isHttps = strstr( url, "https://" ) == url;
    const int start = 7 + isHttps;
    const size_t host_ip_name_length = query - ( url + start ); /* "http://" */
    const size_t query_length = strlen( query );
    memset( http_get_request, 0, sizeof http_get_request );
    CHECK( host_ip_name_length >= 5 );
    snprintf( http_get_request, HTTP_GET_LENGTH - query_length,
              "GET %s HTTP/1.0\nHost: ", query );
    CHECK( http_get_request[ HTTP_GET_LENGTH ] == '\0' );
    strncat( http_get_request, url + start, host_ip_name_length );
    CHECK( http_get_request[ HTTP_GET_LENGTH - 1 ] == '\0' );
    strcat( http_get_request, "\n\n" );
    CHECK( http_get_request[ HTTP_GET_LENGTH ] == '\0' );
    DEBUG( fprintf( stderr, "http_get_request = '%s'\n", http_get_request );)

    {
      const size_t http_get_request_length = strlen( http_get_request );
      size_t sent_length = 0;
      const CURLcode error =
        curl_easy_send( curl, http_get_request, http_get_request_length,
                        &sent_length );
      DEBUG( fprintf( stderr, "sent_length = %lu\n", sent_length ); )
      result = AND2( ! error, sent_length == http_get_request_length );

      if ( ! result ) {
        failureMessage( "Failed to send HTTP GET request because: %s.",
                        curl_easy_strerror( error ) );
      }
    }
  } else {
    failureMessage( "Missing query in URL '%s'.", url );
  }

  POST0( is_initialized_http_connections() );
  return result;
}



/******************************************************************************
PURPOSE: skip_hostname - Get pointer past host name of URL or 0 if none.
INPUTS:  const char* url   HTTP/S URL to examine.
RETURNS: const char* pointer into url after host name.
******************************************************************************/

static const char* skip_hostname( const char* url ) {
  PRE04( initialized, url, *url,
         OR2( strstr( url, "http://"  ) == url,
              strstr( url, "https://" ) == url ) );
  const char* result = 0;
  const char* slash = strchr( url, '/' );

  if ( slash ) {
    slash = strchr( slash + 1, '/' );

    if ( slash ) {
      slash = strchr( slash + 1, '/' );

      if ( slash ) {
        result = slash;
      }
    }
  }

  POST0( OR2( result == 0, *result == '/' ) );
  return result;
}



/******************************************************************************
PURPOSE: skip_http_header - Read and discard the HTTP/S header and return the
         subsequent message content size in bytes.
INPUTS:  CURL* curl  CURL object to read from.
RETURNS: FILE* stream to read contents from, else 0 and failure message.
******************************************************************************/

static FILE* skip_http_header( CURL* curl ) {
  PRE02( initialized, curl );
  FILE* result = 0;
  long temp = 0; /* CURL BUG: Compiling on Linux -O arg is long! */
  const CURLcode error = curl_easy_getinfo(curl, CURLINFO_LASTSOCKET, &temp );
  curl_socket_t socket = (curl_socket_t) temp;
  DEBUG( fprintf( stderr, "socket = %d\n", (int) socket ); )

  if ( error ) {
    failureMessage( "Failed to get socket because: %s\n",
                    curl_easy_strerror( error ) );
  } else {
    result = fdopen( socket, "rb" ); /* Don't close, else socket closed! */

    if ( ! result ) {
      failureMessage( "Failed to get stream for socket." );
    } else {
      enum { LINE_LENGTH = 255 };
      char line[ LINE_LENGTH + 1 ] = "";
      memset( line, 0, sizeof line );
      DEBUG( fprintf( stderr, "Received header:\n" ); )

      while ( AND2( fgets( line, LINE_LENGTH, result ),
                    ! IN3( *line, '\n', '\r' ) ) ) {
        DEBUG( fprintf( stderr, "%s", line ); )
        memset( line, 0, sizeof line );
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: set_socket_timeout - Set socket read/write timeout.
INPUTS:  CURL* curl            CURL object to set.
         int mode              READ_MODE or WRITE_MODE.
         int timeout_seconds   Seconds to timeout waiting.
RETURNS: int 1 if successful, else 0 and failure message is printed to stderr.
******************************************************************************/

static int set_socket_timeout( CURL* curl, int mode, int seconds ) {
  PRE04( initialized, curl, IN3( mode, READ_MODE, WRITE_MODE ), seconds >= 0 );
  int result = 0;
  long temp = 0; /* CURL BUG: Compiling on Linux -O arg is long! */
  const CURLcode error = curl_easy_getinfo(curl, CURLINFO_LASTSOCKET, &temp );
  curl_socket_t socket = (curl_socket_t) temp;
  DEBUG( fprintf( stderr, "socket = %d\n", (int) socket ); )

  if ( error ) {
    failureMessage( "Failed to get socket because: %s\n",
                    curl_easy_strerror( error ) );
  } else if ( ! wait_on_socket( socket, mode, seconds ) ) {
    failureMessage( "Failed waiting on socket." );
  } else {
    result = 1;
  }

  POST0( IS_BOOL( result ) );
  DEBUG( fprintf( stderr, "set_socket_timeout() returning %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: set_socket_blocking - Set socket to blocking (so freads wait for EOF).
INPUTS:  CURL* curl            CURL object to set.
         int mode              READ_MODE or WRITE_MODE.
         int timeout_seconds   Seconds to timeout waiting.
RETURNS: int 1 if successful, else 0 and failure message is printed to stderr.
******************************************************************************/

static int set_socket_blocking( CURL* curl ) {
  PRE02( initialized, curl );

#ifdef _WIN32
  const int result = curl != 0; /* No-op implementation of this routine. */
#else
  int result = 0;
  long temp = 0; /* CURL BUG: Compiling on Linux -O arg is long! */
  const CURLcode error = curl_easy_getinfo(curl, CURLINFO_LASTSOCKET, &temp );
  curl_socket_t socket = (curl_socket_t) temp;
  DEBUG( fprintf( stderr, "socket = %d\n", (int) socket ); )

  if ( error ) {
    failureMessage( "Failed to get socket because: %s\n",
                    curl_easy_strerror( error ) );
  } else {
    int flags = fcntl( socket, F_GETFL );

    if ( flags < 0 ) {
      failureMessage( "Failed to get socket flags." );
    } else {
      DEBUG( fprintf( stderr, "get socket %d flags = %x\n", socket, flags ); )
      flags ^= O_NONBLOCK;
      DEBUG( fprintf( stderr, "set socket %d flags = %x\n", socket, flags ); )

      if ( fcntl( socket, F_SETFL, flags ) < 0 ) {
        failureMessage( "Failed to set socket to non-blocking state." );
      } else {
        result = 1;
      }
    }
  }
#endif

  POST0( IS_BOOL( result ) );
  DEBUG( fprintf( stderr, "set_socket_blocking() returning %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: wait_on_socket - Check and wait to read from or write to a socket.
INPUTS:  curl_socket_t socket  Socket to wait on.
         int mode              READ_MODE or WRITE_MODE.
         int timeout_seconds   Seconds to timeout waiting or 0 if indefinitely.
RETURNS: int 1 if socket is ready, else 0 if timed-out or other error.
******************************************************************************/

static int wait_on_socket( curl_socket_t socket, int mode,
                           int timeout_seconds ) {
  PRE04( initialized, socket >= 0, IN3( mode, READ_MODE, WRITE_MODE ),
         timeout_seconds >= 0 );

  int result = 0;

#ifdef _WIN32

  /* Use 'internal' CURL routine: */

  const long timeout_ms = timeout_seconds ? timeout_seconds * 1000 : -1;

  if ( mode == READ_MODE ) {
    const int mask = Curl_socket_check(socket, 0, CURL_SOCKET_BAD, timeout_ms);
    result = AND2( mask > 0, mask & CURL_CSELECT_IN );
    DEBUG( fprintf( stderr, "socket = %d, timeout_ms = %ld, mask = %0x, "
                   "CURL_CSELECT_IN = %0x, result = %d\n", 
                   socket, timeout_ms, mask, CURL_CSELECT_IN, result ); )
  } else {
    const int mask = Curl_socket_check(socket, CURL_SOCKET_BAD, 0, timeout_ms);
    result = AND2( mask > 0, mask & CURL_CSELECT_OUT );
    DEBUG( fprintf( stderr, "socket = %d, timeout_ms = %ld, mask = %0x, "
                    "CURL_CSELECT_OUT = %0x, result = %d\n", 
                    socket, timeout_ms, mask, CURL_CSELECT_OUT, result ); )
  }

#else

  /* Use POSIX select() routine: */

  struct timeval timeout;
  fd_set read;
  fd_set write;
  fd_set error;
  FD_ZERO( &read );
  FD_ZERO( &write );
  FD_ZERO( &error );
  memset( &timeout, 0, sizeof timeout );
  timeout.tv_sec = timeout_seconds;

  if ( mode == READ_MODE ) {
    FD_SET( socket, &read );
  } else {
    FD_SET( socket, &write );
  }

  /*
   * select() blocks until either timeout expires or
   * a socket is ready to read from or write to
   * then returns the number of signalled/ready sockets or
   * 0 if timed-out or -1 if an error occurred.
   */

  DEBUG( fprintf( stderr, "calling select on socket %d mode %d timeout %d\n",
                  socket + 1, mode, timeout_seconds ); )

  result = select( socket + 1, &read, &write, &error,
                   timeout_seconds ? &timeout : 0 );

  DEBUG( fprintf( stderr, "select() result = %d\n", result ); )
  result = ( result == 1 );

#endif

  DEBUG( fprintf( stderr, "wait_on_socket() result = %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: reverse_2byte_words_if_little_endian - Reverse 2-byte words of array
         if little-endian platform.
INPUTS:  size_t count  Number of words in array.
         void* array   Array of words to reverse.
OUTPUTS: void* array   Array of reversed words.
******************************************************************************/

static void reverse_2byte_words_if_little_endian( size_t count, void* array ) {

#if IS_LITTLE_ENDIAN

  PRE02( count, array );
  assert_static( sizeof (unsigned short) == 2 );
  unsigned short* const a = (unsigned short*) array;
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    const unsigned short value = a[ index ];
    const unsigned short new_value =
      ( value & 0xff00 ) >> 8 |
      ( value & 0x00ff ) << 8;
    CHECK( IMPLIES( IN3( value, 0, 0xffff ), new_value == value ) );
    a[ index ] = new_value;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: reverse_4byte_words_if_little_endian - Reverse 4-byte words of array
         if little-endian platform.
INPUTS:  size_t count  Number of words in array.
         void* array   Array of words to reverse.
OUTPUTS: void* array   Array of reversed words.
******************************************************************************/

static void reverse_4byte_words_if_little_endian( size_t count, void* array ) {

#if IS_LITTLE_ENDIAN

  PRE02( count, array );
  assert_static( sizeof (unsigned int) == 4 );
  unsigned int* const a = (unsigned int*) array;
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    const unsigned int value = a[ index ];
    const unsigned int new_value =
      ( value & 0xff000000 ) >> 24 |
      ( value & 0x00ff0000 ) >>  8 |
      ( value & 0x0000ff00 ) <<  8 |
      ( value & 0x000000ff ) << 24;
    CHECK( IMPLIES( IN3( value, 0, 0xffffffff ), new_value == value ) );
    a[ index ] = new_value;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: reverse_8byte_words_if_little_endian - Reverse 8-byte words of array
         if little-endian platform.
INPUTS:  size_t count  Number of words in array.
         void* array   Array of words to reverse.
OUTPUTS: void* array   Array of reversed words.
******************************************************************************/

static void reverse_8byte_words_if_little_endian( size_t count, void* array ) {

#if IS_LITTLE_ENDIAN

  PRE02( count, array );
  assert_static( sizeof (unsigned long long) == 8 );
  unsigned long long* const a = (unsigned long long*) array;
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    const unsigned long long value = a[ index ];
    const unsigned long long new_value =
    ( value & 0xff00000000000000ULL ) >> 56 |
    ( value & 0x00ff000000000000ULL ) >> 40 |
    ( value & 0x0000ff0000000000ULL ) >> 24 |
    ( value & 0x000000ff00000000ULL ) >>  8 |
    ( value & 0x00000000ff000000ULL ) <<  8 |
    ( value & 0x0000000000ff0000ULL ) << 24 |
    ( value & 0x000000000000ff00ULL ) << 40 |
    ( value & 0x00000000000000ffULL ) << 56;
    CHECK( IMPLIES( IN3(value, 0, 0xffffffffffffffffULL), new_value == value));
    a[ index ] = new_value;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: encode_spaces_and_percents - Change spaces to %20 and %s to %25.
INPUTS:  const char* string  String to encode.
RETURNS: char* Copy of string with spaces changed to %20 and %s to %25.
NOTES:   Free when no longer needed.
******************************************************************************/

static char* encode_spaces_and_percents( const char* string ) {
  PRE0( string );
  char* result = 0;
  const int contains_key = strstr( string, "key=" ) != 0;

  if ( contains_key ) { /* Don't encode string since it will invalidate key. */
    result = strdup( string );
  } else {
    const size_t length = strlen( string );
    const char* s = string;
    size_t count = 0;

    do {
      const char* const space = strchr( s, ' ' );
      const char* const percent = strchr( s, '%' );

      if ( space ) {
        ++count;
        s = space + 1;
      } else if ( percent ) {
        ++count;
        s = percent + 1;
      } else {
        ++s;
      }

    } while ( *s );

    result = (char*) malloc( length + 1 + count + count );

    if ( result ) {
      char* r = result;

      for ( s = string; *s; ++s ) {

        if ( *s == ' ' ) {
          *r++ = '%';
          *r++ = '2';
          *r++ = '0';
        } else if ( *s == '%' ) {
          *r++ = '%';
          *r++ = '2';
          *r++ = '5';
        } else {
          *r++ = *s;
        }
      }

      *r = '\0';
    }
  }

  POST0( IMPLIES( result, strchr( result, ' ' ) == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: is_text - Is string non-empty and only contains isprint() characters?
INPUTS:  const char* string  String to check.
RETURNS: 1 if each character in string isprint(), else 0.
******************************************************************************/

static int is_text( const char* string ) {
  PRE0( string );
  int result = *string != '\0';

  for ( ; AND2( *string, result ); ++string ) {
    result = OR2( isprint( *string ), isspace( *string ) );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



#ifdef __cplusplus
}
#endif
