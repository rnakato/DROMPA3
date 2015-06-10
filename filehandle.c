/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdlib.h>
#include <sys/stat.h>
#include "filehandle.h"

void isfile(char *file){
  struct stat st;
  if(stat(file, &st)){
    fprintf(stderr, "error: cannot open %s.\n", file);
    exit(0);
  }
  return;
}

void remove_file(char *filename){
  struct stat st;
  if(!stat(filename, &st)) remove(filename);
  return;
}

FILE *my_fopen(char *filename, File_Mode mode){
  FILE *IN=NULL;
  switch(mode){
  case FILE_MODE_READ:	IN = fopen(filename, "r"); break;
  case FILE_MODE_WRITE:	IN = fopen(filename, "w"); break;
  case FILE_MODE_WB:	IN = fopen(filename, "wb"); break;
  case FILE_MODE_A:	IN = fopen(filename, "a+"); break;
  default: fprintf(stderr,"[E] Invalid File_Mode: %d.\n", mode); exit(1);
  }
  if(!IN){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  return IN;
}

void my_system(char *command){
  int val = system(command);
  if(val) fprintf(stderr, "system: command error. %s\n", command);
}
