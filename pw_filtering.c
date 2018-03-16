/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "alloc.h"
#include "pw_filtering.h"
#include "macro.h"

#define PRIME 137

typedef struct{
  char *key;
  int val;
} Hash;

static void filtering_eachchr_single(PwParam *p, Mapfile *mapfile, int chr, Strand strand);
static void filtering_eachchr_pair(PwParam *p, Mapfile *mapfile, int chr);
static Hash *alloc_hash(int hashsize);
static void delete_hash(Hash *hashtable, int hashsize);
static int primes(int max);

void check_redundant_reads(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int chr, threshold;
  Strand strand;

  if(p->thre_filter) threshold = p->thre_filter;
  else{
    threshold = max(1, (int)(mapfile->genome->both.n_readname *10/(double)g->genome->len_mpbl));
  }
  printf("\nChecking redundant reads: redundancy threshold >%d\n", threshold);
  mapfile->threshold4filtering = threshold;

  mapfile->cs_raw.nt_all = 0;
  mapfile->cs_raw.nt_nonred = 0;
  mapfile->cs_raw.nt_red = 0;  // calculate library complexity
  double r = p->num4cmp/(double)mapfile->genome->both.n_read_infile;
  
  if(r>1){
    printf("Warning: number of reads is < %d million.\n", (int)(p->num4cmp/NUM_1M));
    mapfile->cs_raw.tv = 1;
  }
  p->r4cmp = r*RAND_MAX;

  for(chr=1; chr<g->chrnum; chr++){
    FLUSH("%s..", g->chr[chr].name);
    if(p->rtype==READTYPE_SINGLE){
      for(strand=0; strand<STRANDNUM; strand++) filtering_eachchr_single(p, mapfile, chr, strand);
    }else{
      filtering_eachchr_pair(p, mapfile, chr);
    }
  }
  mapfile->cs_raw.complexity = mapfile->cs_raw.nt_nonred/(double)mapfile->cs_raw.nt_all;

  printf("done.\n");
  return;
}

static unsigned int calchashvalue(char *key, int hashsize){
  unsigned int h=0;
  do{ h = ((h<<4) + *key - 0x20) % hashsize; }while(*++key);
  return h;
}

static bool hashfunc(Hash *hashtable, char *key, int thre, int hashsize){
  int h = calchashvalue(key, hashsize);
  int cnt = hashsize;
  bool delete = false;
  while(cnt--){
    if(!hashtable[h].val){
      hashtable[h].key = strdup(key);
      hashtable[h].val++;
      goto final;
    }else if(!strcmp(hashtable[h].key, key)){
      hashtable[h].val++;
      if(hashtable[h].val > thre) delete = true;
      goto final;
    }
    h += PRIME;
    if(h >= hashsize) h = h % hashsize;
  }
  printf("Error: Hash Table Full.\n");
  exit(0);
 final:
  return delete;
}

static void filtering_eachchr_single(PwParam *p, Mapfile *mapfile, int chr, Strand strand){
  int i;
  struct seq *seq = &(mapfile->chr[chr].seq[strand]);
  Readarray *read = &(mapfile->readarray[chr][strand]);
  char str[256];

  int ntemp = max(seq->n_read_infile*2, NUM_1M);
  int hashsize = primes(ntemp);
  Hash *hashtable = alloc_hash(hashsize);
  
  for(i=0; i<seq->n_read_infile; i++){
    sprintf(str, "%d", read->F3[i]);
    read->delete[i] = hashfunc(hashtable, str, mapfile->threshold4filtering, hashsize);
    if(read->delete[i]) seq->n_read_red++;
    else seq->n_read_nonred++;
  }
  delete_hash(hashtable, hashsize);

  LOG("%s, n_read_infile: %ld, n_read_nonred: %ld, n_read_red: %ld\n",__FUNCTION__, seq->n_read_infile, seq->n_read_nonred, seq->n_read_red);

  // calcuate library complexity
  hashtable = alloc_hash(hashsize);
  for(i=0; i<seq->n_read_infile; i++){
    if(rand() >= p->r4cmp) continue;
    sprintf(str, "%d", read->F3[i]);
    mapfile->cs_raw.nt_all++;
    if(hashfunc(hashtable, str, mapfile->threshold4filtering, hashsize)) mapfile->cs_raw.nt_red++;
    else mapfile->cs_raw.nt_nonred++;
  }
  delete_hash(hashtable, hashsize);

  return;
}

static void filtering_eachchr_pair(PwParam *p, Mapfile *mapfile, int chr){
  int i, F3=-1, F5=-1;
  long nread;
  Strand strand;
  struct seq *seq = NULL;
  Readarray *read = NULL;
  char str[256];

  int ntemp = max(mapfile->chr[chr].both.n_read_infile*2, NUM_1M);
  int hashsize = primes(ntemp);
  Hash *hashtable = alloc_hash(hashsize);

  for(strand=0; strand<STRANDNUM; strand++){
    seq = &(mapfile->chr[chr].seq[strand]);
    nread = seq->n_read_infile;
    read = &(mapfile->readarray[chr][strand]);
    for(i=0; i<nread; i++){
      F3 = read->F3[i];
      F5 = read->F5[i];
      sprintf(str, "%d-%d", min(F3, F5), max(F3, F5));
      read->delete[i] = hashfunc(hashtable, str, mapfile->threshold4filtering, hashsize);
      if(read->delete[i]) seq->n_read_red++;
      else seq->n_read_nonred++;
    }
  }
  delete_hash(hashtable, hashsize);
  
  // calcuate library complexity
  hashtable = alloc_hash(hashsize);
  for(strand=0; strand<STRANDNUM; strand++){
    seq = &(mapfile->chr[chr].seq[strand]);
    nread = seq->n_read_infile;
    read = &(mapfile->readarray[chr][strand]);
    for(i=0; i<nread; i++){
      if(rand() >= p->r4cmp) continue;
      mapfile->cs_raw.nt_all++;
      F3 = read->F3[i];
      F5 = read->F5[i];
      sprintf(str, "%d-%d", min(F3, F5), max(F3, F5));
      if(hashfunc(hashtable, str, mapfile->threshold4filtering, hashsize)) mapfile->cs_raw.nt_red++;
      else mapfile->cs_raw.nt_nonred++;
    }
  }
  delete_hash(hashtable, hashsize);
  
  return;
}

static Hash *alloc_hash(int hashsize){
  Hash *hashtable = (Hash *)my_calloc(hashsize, sizeof(Hash), "hashtable");
  return hashtable;
}

static void delete_hash(Hash *hashtable, int hashsize){
  int i;
  for(i=0; i<hashsize; i++) MYFREE(hashtable[i].key);
  MYFREE(hashtable);
  return;
}

static int primes(int max){ 
  int m, i;
  int mhalf = max/2;

  for(m=max; m>=2; m--){
    for(i=3; i<mhalf; i+=2){
      if(!(m%i)) break;
    }
    if(i>=mhalf) break;
  }
  return m;
}
