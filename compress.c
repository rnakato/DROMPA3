/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "compress.h"
#include "filehandle.h"

#define BUFFSIZE 1024

void compress2gz(char *filename){
  FILE *IN=NULL, *OUT=NULL;
  char filename_gz[256];
  z_stream z;
  int status, flush;
  unsigned int count;
  unsigned char inbuf[BUFFSIZE], outbuf[BUFFSIZE];

  /* initinalization */
  z.zalloc = Z_NULL;
  z.zfree  = Z_NULL;
  z.opaque = Z_NULL;
  if(deflateInit2(&z, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 31, 8, Z_DEFAULT_STRATEGY) != Z_OK){   // gzip
    fprintf(stderr, "deflateInit error: %s\n", z.msg);
    exit(1);
  }
  z.avail_in = 0; 
  z.next_out = outbuf;
  z.avail_out = BUFFSIZE;
  flush = Z_NO_FLUSH;

  /* open stream */
  sprintf(filename_gz, "%s.gz", filename);
  IN = my_fopen(filename, FILE_MODE_READ);
  OUT = my_fopen(filename_gz, FILE_MODE_WRITE);

  /* compress */
  while(1){
    if(!z.avail_in){ 
      z.next_in = inbuf;  
      z.avail_in = fread(inbuf, 1, BUFFSIZE, IN); 
      if (z.avail_in < BUFFSIZE) flush = Z_FINISH;
    }
    status = deflate(&z, flush); 
    if(status == Z_STREAM_END) break;
    if(status != Z_OK){ 
      fprintf(stderr, "deflate error: %s\n", z.msg);
      exit(1);
    }
    if(!z.avail_out){ 
      if(fwrite(outbuf, 1, BUFFSIZE, OUT) != BUFFSIZE){
	fprintf(stderr, "Write error\n");
	exit(1);
      }
      z.next_out = outbuf; 
      z.avail_out = BUFFSIZE; 
    }
  }

  count = BUFFSIZE - z.avail_out;
  if(count){
    if(fwrite(outbuf, 1, count, OUT) != count){
      fprintf(stderr, "Write error\n");
      exit(1);
    }
  }

  /* free */
  if(deflateEnd(&z) != Z_OK){
    fprintf(stderr, "deflateEnd: %s\n", z.msg);
    exit(1);
  }
  fclose(IN);
  fclose(OUT);
  remove(filename);
  return;
}
