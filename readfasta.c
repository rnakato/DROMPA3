/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "readfasta.h"
#include "filehandle.h"
#include "alloc.h"
#include "macro.h"

/* Nが入っていた場合は-1を返す。 */
TYPE_FASTAGCARRAY *make_fastaGCarray(char *filename, long length, int flen4gc){
  int i,s,e, n=0;
  int state=0;
  char c;

  TYPE_FASTAGCARRAY *array = (TYPE_FASTAGCARRAY *)my_calloc(length, sizeof(TYPE_FASTAGCARRAY), "fastaGCarray");
  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  while((c = fgetc(IN)) != EOF){
    switch(state){
    case 0:
      if(c=='>') state=1;
      break;
    case 1:    /*header*/
      if(c=='\n') state=2;
      break;
    case 2:    /*body*/
      if(c=='>') goto final;
      else if(isalpha((int)c)){
	if(c=='G' || c=='C' || c=='g' || c=='c'){
	  s = max(0, n - flen4gc);
	  e = n;
	  for(i=s; i<e; i++){
	    if(array[i] != -1) array[i]++;
	  }
	}else if(c=='A' || c=='T' || c=='a' || c=='t'){
	  /* none */
	}else{  /* N and others */
	  s = max(0, n - flen4gc);
	  e = n;
	  for(i=s; i<e; i++) array[i] = -1;
	}

	n++;
	if(n > length){
	  fprintf(stderr, "ERROR:%s:length %ld < n %d\n",__FUNCTION__, length, n);
	  exit(0);
	}
      }
      break;
    }
  }

  final:
  fclose(IN);
  return array;
}
