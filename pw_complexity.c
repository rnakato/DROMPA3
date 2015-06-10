/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "pw_complexity.h"

static void func(PwParam *p, Mapfile *mapfile, RefGenome *g, long nread, CompStats *cs){
  int i, chr, nread_chr;
  Strand strand;

  cs->nt_all = 0;
  cs->nt_nonred = 0;
  cs->nt_red = 0;

  double r = p->num4cmp/(double)nread;
  if(r>1){
    printf("Warning: number of reads is < %d million.\n", (int)(p->num4cmp/NUM_1M));
    cs->tv = 1;
  }
  double rnum = r*RAND_MAX;

  /* extract nread*p reads */
  for(chr=1; chr<g->chrnum; chr++){
    for(strand=0; strand<STRANDNUM; strand++){
      nread_chr = mapfile->chr[chr].seq[strand].n_read_infile;
      for(i=0; i<nread_chr; i++){
	if(rand() >= rnum) continue;
	cs->nt_all++;
	if(!mapfile->readarray[chr][strand].delete[i]) cs->nt_nonred++;
	else cs->nt_red++;
      }
    }
  }
  /* store the stats */
  cs->complexity = cs->nt_nonred/(double)cs->nt_all;

  return;
}

void calc_complexity(PwParam *p, Mapfile *mapfile, RefGenome *g){
  srand((unsigned)time(NULL));

  /* get nread */
  /* ignore peak regions */
  long nread = mapfile->genome->both.n_read_infile;
  //  long nread_nonred = mapfile->genome->both.n_read_nonred;

  func(p, mapfile, g, nread,        &(mapfile->cs_raw));
  //  func(p, mapfile, g, nread_nonred, &(mapfile->cs_nonred));
  
  return;
}

