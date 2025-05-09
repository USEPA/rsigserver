
#ifndef HTTP_CONNECTION_H
#define HTTP_CONNECTION_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: http_connection.h - Routines for opening HTTP GET URL socket
         connections and reading ASCII and binary data from them.
NOTES:   Uses curl library to create sockets. Does not support HTTPS.
         Example usage:

#include <stdio.h>
#include <string.h>

#include <http_connection.h>

int main( void ) {
  const char* const url =
    "http://rtpmeta.epa.gov/cgi-bin/rsigserver?"
    "SERVICE=wcs&VERSION=1.0.0&REQUEST=GetCoverage&COVERAGE=uvnet.irradiance&"
    "TIME=1996-01-01T00:00:00Z/1996-01-02T23:59:59Z&BBOX=-118,33,-117,35&"
    "FORMAT=ascii";
  const int timeout = 300;
  int lines = 0;

  if ( initialize_http_connections() ) {
    FILE* stream = 0;
    void* connection = open_http_connection( url, timeout, &stream );

    if ( connection ) {
      char line[ 256 ] = "";
      memset( line, 0, sizeof line );

      fprintf( stderr, "Reading data:...\n" );

      while ( read_http_connection_line( stream, sizeof line / sizeof *line,
                                         line ) ) {
        char timestamp[ 24 + 1 ] = "";
        float longitude = 0.0;
        float latitude = 0.0;
        int station = 0;
        float data = 0.0;
        memset( timestamp, 0, sizeof timestamp );

        if ( lines == 0 ) {
          ++lines;
          printf( "%s\n", line );
        } else if ( sscanf( line, "%24s %f %f %d %f",
                            timestamp, &longitude, &latitude, &station,
                            &data ) == 5 ) {
          ++lines;
          printf( "%s\t%f\t%f\t%d\t%f\n",
                  timestamp, longitude, latitude, station, data );
        }
      }

      fprintf( stderr, "Read %d lines.\n", lines );
      close_http_connection( connection ), connection = 0;
    }
  }

  return lines < 2;
}

HISTORY: 2012-06-06 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h> /* For FILE. */

/*================================ FUNCTIONS ================================*/

extern int initialize_http_connections( void );
extern int is_initialized_http_connections( void );
extern void* open_http_connection(const char* url, int timeout, FILE** stream);
extern void close_http_connection( void* http_connection );
extern int read_http_connection_line( FILE* stream, size_t size, char line[] );
extern int read_http_connection_array( FILE* stream, size_t count,
                                       size_t word_size, void* array );

#ifdef __cplusplus
}
#endif

#endif /* HTTP_CONNECTION_H */


