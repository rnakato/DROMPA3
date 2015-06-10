/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _ALLOC_H_
#define _ALLOC_H_

#include <stdio.h>
#include <stdlib.h>

#ifndef MEMTEST
#define MYFREE(p) {if(p){free(p); (p)=NULL;} }
#endif

void *my_calloc(size_t n, size_t s, char *name);
void *my_realloc(void *p, size_t s, char *name);

#endif /* _ALLOC_H_ */
