
#ifndef STREAM_H
#define STREAM_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Stream.h - Declare streams for reading and writing to ASCII and
         portable binary (XDR/IEEE-754/MSB) files, pipes and sockets.
         Provides an efficient, portable (multi-platform, cross-language
         compatible), convenient and compatible alternative to fopen, popen
         and socket and associated fread/fwrite/xdr_vector calls.

NOTES:   Example usage:

         #include <Stream.h>

         int main( void ) {
           int ok = 0;
           Stream* stream = newPipeStream( "rsh -l cws sequoia ls", "r" );

           if ( stream ) {

             while ( stream->ok( stream ) && ! stream->isAtEnd( stream ) ) {
               char string[ 10 ];
               stream->readString( stream, string, 10 );

               if ( stream->ok( stream ) ) {
                 printf( "%s\n", string );
               }
             }

             ok = stream->ok( stream );
             FREE_OBJECT( stream );
           }

           return ! ok;
         }

         cc -n32 -mips3 -xansi -fullwarn -g \
            -I. -I/usr/include -I/usr/include/sys -I/usr/include/rpc \
            -I/usr/include/netinet -I../../../../include -I. \
            -o test_pipe test_pipe.c \
            -L/usr/lib32/mips3 -L/usr/lib32 -L../../../../lib/IRIX \
            -lStream.debug -lMemory.debug -lFailure.debug -lc

         test_pipe

         Uses a pipe to list (up to the first 9 letters of) names of files
         in the home directory on the remote host.

HISTORY: 1995/12, Todd Plessel, EPA/MMTSI, Created.
         1998/11, Todd Plessel, EPA/LM, Re-written using Integer,Real,CQSP.
         1999/11, Todd Plessel, EPA/LM, Re-written using embedded pointers.
STATUS:  reviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h> /* For FILE. */

#include <BasicNumerics.h> /* For Integer, Real.                   */
#include <Memory.h>        /* To export FREE_OBJECT() for clients. */

/*================================== TYPES ==================================*/

typedef struct StreamData StreamData;
typedef struct Stream Stream;

/* Constructors: */

extern Stream* newFileStream( const char* fileName, const char* mode );
extern Stream* newPipeStream( const char* command,  const char* mode );
extern Stream* newServerSocketStream( Integer port );
extern Stream* newClientSocketStream( Integer port, const char* host );

struct Stream {

  /* Commands: */

  void (*free)( Stream* self ); /* Destructor used via FREE_OBJECT() macro. */
  void (*flush)( Stream* self );

  void (*seekFromStart)(    Stream* self, Integer offset );
  void (*seekFromEnd)(      Stream* self, Integer offset );
  void (*seekFromCurrent)(  Stream* self, Integer offset );

  void (*readString)( Stream* self, char* s, Integer n );
  void (*readWord)(   Stream* self, char* s, Integer n );
  void (*readByte)(   Stream* self, void* x );

  void (*read8BitInteger)(    Stream* self, Integer* x );
  void (*read16BitInteger)(   Stream* self, Integer* x );
  void (*read32BitInteger)(   Stream* self, Integer* x );
  void (*read64BitInteger)(   Stream* self, Integer* x );

  void (*read32BitReal)(      Stream* self, Real*    x );
  void (*read64BitReal)(      Stream* self, Real*    x );

  void (*readBytes)(          Stream* self, void*    a, Integer n );
  void (*readUpToNBytes)(     Stream* self, void*    a, Integer n,
                              Integer* actualCount );

  void (*read8BitIntegers)(   Stream* self, Integer* a, Integer n );
  void (*read16BitIntegers)(  Stream* self, Integer* a, Integer n );
  void (*read32BitIntegers)(  Stream* self, Integer* a, Integer n );
  void (*read64BitIntegers)(  Stream* self, Integer* a, Integer n );

  void (*read32BitReals)(     Stream* self, Real*    a, Integer n );
  void (*read64BitReals)(     Stream* self, Real*    a, Integer n );


  void (*writeString)(        Stream* self, const char* s, ... );
  void (*writeByte)(          Stream* self, unsigned char x );

  void (*write8BitInteger)(   Stream* self, Integer x );
  void (*write16BitInteger)(  Stream* self, Integer x );
  void (*write32BitInteger)(  Stream* self, Integer x );
  void (*write64BitInteger)(  Stream* self, Integer x );

  void (*write32BitReal)(     Stream* self, Real    x );
  void (*write64BitReal)(     Stream* self, Real    x );

  void (*writeBytes)(         Stream* self, const void*    a, Integer n );

  void (*write8BitIntegers)(  Stream* self, const Integer* a, Integer n );
  void (*write16BitIntegers)( Stream* self, const Integer* a, Integer n );
  void (*write32BitIntegers)( Stream* self, const Integer* a, Integer n );
  void (*write64BitIntegers)( Stream* self, const Integer* a, Integer n );

  void (*write32BitReals)(    Stream* self, const Real*    a, Integer n );
  void (*write64BitReals)(    Stream* self, const Real*    a, Integer n );


  /* Queries: */

  Integer     (*invariant)(  const Stream* self ); /* Must always return 1. */
  Integer     (*ok)(         const Stream* self ); /* Did last command work? */
  Integer     (*isReadable)( const Stream* self );
  Integer     (*isWritable)( const Stream* self );
  Integer     (*isSeekable)( const Stream* self );
  Integer     (*isBlocking)( const Stream* self );
  Integer     (*isAtEnd)(    const Stream* self ); /* Not thread-safe. */

  Integer     (*offset)(     const Stream* self ); /* Not thread-safe? */
  Integer     (*size)(       const Stream* self );
  const char* (*name)(       const Stream* self );
  FILE*       (*file)(       const Stream* self ); /* HACK: Demeter! */
  Integer     (*descriptor)( const Stream* self ); /* HACK: Demeter! */


  /* Private: */

  StreamData* data;
};

#ifdef __cplusplus
}
#endif

#endif /* STREAM_H */

