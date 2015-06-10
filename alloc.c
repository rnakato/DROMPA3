/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "alloc.h"

void *my_calloc(size_t n, size_t s, char *name){
  void *p;
  p = calloc(n,s);
  if(!p){
    fprintf(stderr,"[E]failed calloc: %s\n", name); 
    exit(1);
  }
  return p;
}

void *my_realloc(void *p, size_t s, char *name){
  p = realloc(p,s);
  if(!p){
    fprintf(stderr,"[E]failed realloc: %s\n", name); 
    exit(1);
  }
  return p;
}
