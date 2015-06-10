/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "pw_param_new.h"
#include "alloc.h"
#include "memtest.h"

PwParam *pwparam_new(){
  PwParam *p = (PwParam *)my_calloc(1, sizeof(PwParam), "PwParam");

  p->binsize = BINSIZE_DEFAULT;
  p->fraglen = FRAGMENT_LEN;
  p->max_fraglen = MAX_FRAGMENT_LEN;
  p->pcrfilter = true;

  p->num4rpm = NUM4RPM_DEFAULT;
  p->num4depth = NUM4DEPTH_DEFAULT;

  //  p->mpsize  = MPBLSIZE_DEFAULT;
  p->mpthre  = THRE_LOW_MAPPABILITY;
  p->num4cmp = NUM4CMP_DEFAULT*NUM_1M;

  p->flen4gc = FLEN4GC_DEFAULT;

  p->gcdepth = 1;
  return p;
}

Mapfile *mapfile_new(int chrnum, PwParam *pwparam){
  int chr;
  Strand strand;
  Mapfile *p = (Mapfile *)my_calloc(1, sizeof(Mapfile), "mapfile");
  p->genome     = (SeqStats *)my_calloc(1,        sizeof(SeqStats), "mapfile->genome");
  p->chr        = (SeqStats *)my_calloc(chrnum, sizeof(SeqStats), "mapfile->chr");
  p->wstats.genome = (WigStatsMember *)my_calloc(1,        sizeof(WigStatsMember), "mapfile->wstats.genome");
  p->wstats.chr    = (WigStatsMember *)my_calloc(chrnum, sizeof(WigStatsMember), "mapfile->wstats.chr");
  p->readarray = (Readarray **)my_calloc(chrnum, sizeof(Readarray *), "mapfile->readarray");

  /* Wigstats */
  for(chr=1; chr<chrnum; chr++){
    p->readarray[chr] = (Readarray *)my_calloc(2, sizeof(Readarray), "mapfile->readarray[chr]");
    for(strand=0; strand<STRANDNUM; strand++){
      p->readarray[chr][strand].F3 = (int *)my_calloc(READARRAY_NUM, sizeof(int), "mapfile->readarray[chr].F3");
      if(pwparam->rtype==READTYPE_PAIR) p->readarray[chr][strand].F5 = (int *)my_calloc(READARRAY_NUM, sizeof(int), "mapfile->readarray[chr].F5");
      p->readarray[chr][strand].weight = (int *)my_calloc(READARRAY_NUM, sizeof(int), "mapfile->readarray[chr].weight");
      p->readarray[chr][strand].delete = (bool *)my_calloc(READARRAY_NUM, sizeof(bool), "mapfile->readarray[chr].delete");
      if(pwparam->enrichfile) p->readarray[chr][strand].ignore = (bool *)my_calloc(READARRAY_NUM, sizeof(bool), "mapfile->readarray[chr].ignore");
      p->readarray[chr][strand].narray = READARRAY_NUM;
    }
  }
  return p;
}

void mapfile_delete(Mapfile *p, int chrnum){
  int chr;
  for(chr=1; chr<chrnum; chr++){
    MYFREE(p->wstats.chr[chr].darray_all);
    MYFREE(p->wstats.chr[chr].darray_bg);
  }
  MYFREE(p->readarray);
  MYFREE(p->wstats.genome->darray_all);
  MYFREE(p->wstats.genome->darray_bg);
  MYFREE(p->wstats.chr);
  MYFREE(p->wstats.genome);
  MYFREE(p->chr);
  MYFREE(p->genome);
  MYFREE(p);
  return;
}

void pwparam_delete(PwParam *p, RefGenome *g){
  if(p->bedfilename) bedfile_delete(p->enrichfile, g->chrnum);
  MYFREE(p);
  return;
}
