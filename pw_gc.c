/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <string.h>
#include "pw_gc.h"
#include "alloc.h"
#include "readfile.h"
#include "stringp.h"
#include "filehandle.h"
#include "readfasta.h"
#include "macro.h"

#define TYPE_READARRAY unsigned short
#define FRAG_IGNORE 5
#define GCDIST_THRE 1e-5
#define GCDEPTH_THRE 1e-3

static void make_GCdist(PwParam *p, Mapfile *mapfile, RefGenome *g);
static int *makeGCarray_genome(TYPE_FASTAGCARRAY *fastaGCarray, char *ar, long chrlen, int flen4gc);
static int *makeGCarray_read(Mapfile *mapfile, TYPE_FASTAGCARRAY *fastaGCarray, char *ar, int chr, long chrlen, int flen4gc);
static void weight_read(PwParam *p, Mapfile *mapfile, RefGenome *g);
static void output_GCdist(PwParam *p, Mapfile *mapfile, RefGenome *g);
static char *make_fragarray(PwParam *p, RefGenome *g, int chr);

void GCnorm(PwParam *p, Mapfile *mapfile, RefGenome *g){
  printf("\nNormalize with GC distribution...\n");

  /* estimate the GC distribution of the sample */
  make_GCdist(p, mapfile, g);
  /* weight each read */
  weight_read(p, mapfile, g);

  return;
}

static void make_GCdist(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int i;
  /* use longest chromosome as reference genome */
  int chrref = g->chrmax;
  printf("chromosome for GC distribution: %s\n",g->chr[chrref].name);

  /* make flagarray which contains the position to be considered */
  char *ar = make_fragarray(p, g, chrref);

  /* make fastaGCarray */
  char *fstfile = alloc_str_new(p->genomefile, 50);
  sprintf(fstfile, "%s/%s.fa", p->genomefile, g->chr[chrref].name);

  TYPE_FASTAGCARRAY *fastaGCarray = make_fastaGCarray(fstfile, g->chr[chrref].len, p->flen4gc);
  MYFREE(fstfile);

  /* make GCarray for genome */
  g->GCdist = makeGCarray_genome(fastaGCarray, ar, g->chr[chrref].len, p->flen4gc);
  /* make GCarray for reads */
  mapfile->GCdist = makeGCarray_read(mapfile, fastaGCarray, ar, chrref, g->chr[chrref].len, p->flen4gc);

  MYFREE(fastaGCarray);
  MYFREE(ar);

  int max=0;
  for(i=0; i<=p->flen4gc; i++){
    g->sum_GCdist       += g->GCdist[i];
    mapfile->sum_GCdist += mapfile->GCdist[i];
    if(max < mapfile->GCdist[i]){
      max = mapfile->GCdist[i];
      mapfile->maxGC = i;
    }
  }
  mapfile->GCweight = (double *)my_calloc(p->flen4gc +1, sizeof(double), "mapfile->GCweight");

  /* output GCdist.xls */
  output_GCdist(p, mapfile, g);

  return;
}

static int *makeGCarray_genome(TYPE_FASTAGCARRAY *fastaGCarray, char *ar, long chrlen, int flen4gc){
  int i, gc;
  int *array = (int *)my_calloc(flen4gc+1, sizeof(int), "fragarray");
  long max = chrlen - FRAG_IGNORE - flen4gc;

  for(i= FRAG_IGNORE + flen4gc; i<max; i++){
    if(ar[i]){
      gc = fastaGCarray[i];
      if(gc != -1) array[gc]++;
    }
  }
  return array;
}

static int define_posi4gc(Strand strand, Readarray *read, int i, int flen4gc, long chrlen){
  long posi;
  if(strand==STRAND_PLUS) posi = min(read->F3[i] + FRAG_IGNORE, chrlen -1);
  else                    posi = max(read->F3[i] - FRAG_IGNORE - flen4gc, 0);
  return posi;
}

static int *makeGCarray_read(Mapfile *mapfile, TYPE_FASTAGCARRAY *fastaGCarray, char *ar, int chr, long chrlen, int flen4gc){
  long i, nread;
  int gc, posi;
  Strand strand;
  struct seq *seq = NULL;
  Readarray *read = NULL;
  int *array = (int *)my_calloc(flen4gc+1, sizeof(int), "fragarray");

  for(strand=0; strand<STRANDNUM; strand++){
    seq   = &(mapfile->chr[chr].seq[strand]);
    read  = &(mapfile->readarray[chr][strand]);
    nread = seq->n_read_infile;
    for(i=0; i<nread; i++){
      if(read->delete[i]) continue;
      posi = define_posi4gc(strand, read, i, flen4gc, chrlen);
      if(!ar[posi] || !ar[posi + flen4gc]) continue;
      gc = fastaGCarray[posi];
      if(gc != -1) array[gc]++;
    }
  }
  return array;
}

static void weight_read(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int i, chr;
  int gc, posi;
  long nread;
  Strand strand;
  struct seq *seq = NULL;
  Readarray *read = NULL;
  TYPE_FASTAGCARRAY *fastaGCarray=NULL;

  FLUSH("add weight to reads...");
  char *fstfile = alloc_str_new(p->genomefile, 100);
  
  for(chr=1; chr<g->chrnum; chr++){
    FLUSH("%s...", g->chr[chr].name);
    sprintf(fstfile, "%s/%s.fa", p->genomefile, g->chr[chr].name);

    fastaGCarray = make_fastaGCarray(fstfile, g->chr[chr].len, p->flen4gc);

    for(strand=0; strand<STRANDNUM; strand++){
      seq   = &(mapfile->chr[chr].seq[strand]);
      read  = &(mapfile->readarray[chr][strand]);
      nread = seq->n_read_infile;
      for(i=0; i<nread; i++){
	posi = define_posi4gc(strand, read, i, p->flen4gc, g->chr[chr].len);
	gc = fastaGCarray[posi];
	if(gc != -1) read->weight[i] *= mapfile->GCweight[gc];
	seq->n_read_afterGC += read->weight[i];
      }
      mapfile->chr[chr].both.n_read_afterGC += mapfile->chr[chr].seq[strand].n_read_afterGC;
    }
    MYFREE(fastaGCarray);
  }
  MYFREE(fstfile);
  printf("done.\n");

  /* add_stats_to_genome_readafterGC */
  for(chr=1; chr<g->chrnum; chr++){
    for(strand=0; strand<STRANDNUM; strand++){
      mapfile->genome->seq[strand].n_read_afterGC += mapfile->chr[chr].seq[strand].n_read_afterGC;
    }
    mapfile->genome->both.n_read_afterGC += mapfile->chr[chr].both.n_read_afterGC;
  }
  return;
}

static void output_GCdist(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int i;
  double r_genome=0, r_read=0, r_depth=0;
  char *filename = alloc_str_new(p->output_dir, strlen(p->output_prefix) +100);
  sprintf(filename, "%s/%s.GCdist.xls", p->output_dir, p->output_prefix);

  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);
  fprintf(OUT, "GC\tgenome prop\treads prop\tdepth\tweight\n");
  for(i=0; i<=p->flen4gc; i++){
    r_genome = g->GCdist[i]/(double)g->sum_GCdist;
    r_read = mapfile->GCdist[i]/(double)mapfile->sum_GCdist;
    r_depth = g->GCdist[i] ? mapfile->GCdist[i]/(double)g->GCdist[i]: 0;

    //    fprintf(OUT, "%d\t%d\t%f\t%d\t%f\t%f\t", i, g->GCdist[i], r_genome, mapfile->GCdist[i], r_read, r_depth);
    fprintf(OUT, "%d\t%f\t%f\t%f\t", i, r_genome, r_read, r_depth);
    if(!r_read || r_genome < GCDIST_THRE) mapfile->GCweight[i] = 0;
    else{
      if(p->gcdepth && r_depth < GCDEPTH_THRE) mapfile->GCweight[i] = 1;
      else mapfile->GCweight[i] = r_genome/r_read;
    }
    fprintf(OUT, "%f\n", mapfile->GCweight[i]);
  }
  fclose(OUT);
  printf("fragment distribution is output to %s.\n", filename);

  MYFREE(filename);
  return;
}


static char *make_fragarray(PwParam *p, RefGenome *g, int chr){
  long i;
  int j;
  char *filename = NULL;
  char *ar = NULL;

  // mappability
  if(p->mpfile){
    filename = alloc_str_new(p->mpfile, 100);
    sprintf(filename, "%s_%s_binary.txt", p->mpbinaryfile, g->chr[chr].name);
    ar = read_mosaics_binary(filename, g->chr[chr].len);
  }else{
    ar = (char *)my_calloc(g->chr[chr].len, sizeof(char), "ar");
    for(i=0; i<g->chr[chr].len; i++) ar[i]=1;
  }
  
  // ignore bed regions
  int nbed, beds, bede;
  if(p->bedfilename){
    nbed = p->enrichfile->chr[chr].num;
    for(i=0; i<nbed; i++){
      beds = p->enrichfile->chr[chr].bed[i].s;
      bede = p->enrichfile->chr[chr].bed[i].e;
      for(j=beds; j<=bede; j++) ar[j] = 0;
    }
  }
  return ar;
}

