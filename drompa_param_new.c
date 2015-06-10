/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <math.h>
#include "drompa_param_new.h"
#include "alloc.h"
#include "memtest.h"

DrParam *drparam_new(){
  DrParam *p = (DrParam *)my_calloc(1, sizeof(DrParam), "DrParam");

  //  p->binsize = BINSIZE_DEFAULT;
  p->smoothing = SMOOTHING_DEFAULT;
  p->width4lambda = WIDTH4LAMBDA_DEFAULT;
  p->mpthre = THRE_LOW_MAPPABILITY;

  p->ntype = TYPE_RATIO_TOTALREAD;

  p->pthre_internal = -log10(PTHRE_INTERNAL_DEFAULT);
  p->pthre_enrich = -log10(PTHRE_ENRICH_DEFAULT);
  p->qthre = QTHRE_DEFAULT;

  return p;
}

void drparam_delete(DrParam *p){
  MYFREE(p);
  return;
}

static Samplefile *Samplefile_new(int chrnum){
  Samplefile *p = (Samplefile *)my_calloc(1, sizeof(Samplefile), "Samplefile");
  p->genome = (Library *)my_calloc(1,      sizeof(Library), "Library_genomoe");
  p->chr    = (Library *)my_calloc(chrnum, sizeof(Library), "Library_chr");
  return p;
}

static void Samplefile_delete(Samplefile *p){
  MYFREE(p->genome);
  MYFREE(p->chr);
  MYFREE(p);
  return;
}

static LibCompare *LibCompare_new(int chrnum){
  LibCompare *p = (LibCompare *)my_calloc(1, sizeof(LibCompare), "LibCompare");
  p->genome = (LCstats *)my_calloc(1,      sizeof(LCstats), "LCstats_genome");
  p->chr    = (LCstats *)my_calloc(chrnum, sizeof(LCstats), "LCstats_chr");
  return p;
}

static void LibCompare_delete(LibCompare *p){
  MYFREE(p->genome);
  MYFREE(p->chr);
  MYFREE(p);
  return;
}

SamplePair *SamplePair_new(int num, int chrnum){
  if(num<=0){
    fprintf(stderr, "ERROR:%s:num = %d.\n",__FUNCTION__, num); exit(1);
  }
  SamplePair *p = (SamplePair *)my_calloc(num, sizeof(SamplePair), "SamplePair");

  int i;
  for(i=0; i<num; i++){
    p[i].ChIP  = Samplefile_new(chrnum);
    p[i].Input = Samplefile_new(chrnum);
    p[i].comp  = LibCompare_new(chrnum);
    p[i].binnum = (int *)my_calloc(chrnum, sizeof(int), "sample->binnum");
    p[i].peak  = (Peak *)my_calloc(1, sizeof(Peak), "sample->peak");
    p[i].peak->bs = (struct bs *)my_calloc(PEAKNUM_DEFAULT, sizeof(struct bs), "peak");
    p[i].peak->arraynum = PEAKNUM_DEFAULT;
    p[i].copyC = p[i].copyI = p[i].copycomp = -1;
  }
  return p;
}

void SamplePair_delete(SamplePair *p, int num){
  int i;
  for(i=0; i<num; i++){
    if(p[i].copyC==-1) Samplefile_delete(p[i].ChIP);
    if(p[i].Input->argv && p[i].copyI==-1) Samplefile_delete(p[i].Input);
    if(p[i].copycomp==-1) LibCompare_delete(p[i].comp);
    MYFREE(p[i].peak->bs);
  }
  MYFREE(p);
  return;
}
