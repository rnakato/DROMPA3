/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "seq.h"
#include "alloc.h"
#include "memtest.h"

const char str_strand[][STRANDNUM]={"+","-"};

RefGenome *refgenome_new(){
  RefGenome *p = (RefGenome *)my_calloc(1, sizeof(RefGenome), "RefGenome");
  return p;
}

void refgenome_delete(RefGenome *p){
  int i;
  for(i=1; i<p->chrnum; i++) MYFREE(p->chr[i].name);
  MYFREE(p->chr);
  MYFREE(p);
  return;
}

BedFile *bedfile_new(int chrnum){
  int chr;
  BedFile *p = (BedFile *)my_calloc(1, sizeof(BedFile), "BedFile");
  p->chr = (BedChr *)my_calloc(chrnum, sizeof(BedChr), "bedchr");
  for(chr=1; chr<chrnum; chr++){
    p->chr[chr].arraynum = BEDFILE_MAX;
    p->chr[chr].bed = (struct bed *)my_calloc(p->chr[chr].arraynum, sizeof(struct bed), "struct bed");
  }
  return p;
}

void bedfile_delete(BedFile *p, int chrnum){
  int chr;
  for(chr=1; chr<=chrnum; chr++) MYFREE(p->chr[chr].bed);
  MYFREE(p->chr);
  MYFREE(p);
  return;
}
