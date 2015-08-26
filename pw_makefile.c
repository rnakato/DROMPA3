/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "alloc.h"
#include "common.h"
#include "pw_makefile.h"
#include "readfile.h"
#include "filehandle.h"
#include "stringp.h"
#include "macro.h"
#include "pw_estimate.h"
#include "outputfile.h"

#define PRINTWARNING_W(w) fprintf(stderr, "Warning: Read scaling weight = %.2f. Too much scaling up will bring noisy results.\n", w)

static void makewig_chr(PwParam *p, Mapfile *mapfile, RefGenome *g, int chr);
static TYPE_WIGARRAY *make_wigarray(PwParam *p, Mapfile *mapfile, RefGenome *g, int chr, bool free);
static void normalize_wigarray_to_rpkm(PwParam *p, Mapfile *mapfile, TYPE_WIGARRAY **wigarrayp, RefGenome *g, int chr);
static void get_start_end_of_read(PwParam *p, Mapfile *mapfile, int *s, int *e, int chr, int i, Strand strand, const char *func);

static int comparray(const void *c1, const void *c2){
  TYPE_WIGARRAY n1 = *(TYPE_WIGARRAY *)c1;
  TYPE_WIGARRAY n2 = *(TYPE_WIGARRAY *)c2;
  if(n1<n2) return -1;
  else if(n1==n2) return 0;
  else return 1;
}

static TYPE_WIGARRAY calc_qnt(TYPE_WIGARRAY *array, int binnum, double per){
  int i, num=0;
  for(i=0; i<binnum; i++){
    if(array[i]) num++;
  }

  TYPE_WIGARRAY *sortarray = (TYPE_WIGARRAY *)my_calloc(num, sizeof(TYPE_WIGARRAY), "sortarray");
  num=0;
  for(i=0; i<binnum; i++){
    if(array[i]) sortarray[num++] = array[i];
  }

  qsort(sortarray, num, sizeof(TYPE_WIGARRAY), comparray);
  TYPE_WIGARRAY qnt = sortarray[(int)(num*per)];
  MYFREE(sortarray);

  return qnt;
}

static int define_wstats_thre(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int chr = g->chrmax;
  int binnum = p->binnum_chr[chr];
  TYPE_WIGARRAY *wigarray = make_wigarray(p, mapfile, g, chr, false);


  int i, thre=0, num;
  do{
    num=0;
    thre += N_ARRAY_NB;
    for(i=0; i<binnum; i++){
      if(WIGARRAY2VALUE(wigarray[i]) <= thre ) num++;
    }
    //    printf("thre=%d, num=%d, binnum=%d\n",thre, num, binnum);
  }while(num < binnum/5);

  mapfile->wstats.num95 = calc_qnt(wigarray, binnum, 0.95);

  MYFREE(wigarray);
  return thre;
}

void makewig(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int chr;
  printf("Convert read data to array: \n");

  mapfile->wstats.thre = define_wstats_thre(p, mapfile, g);
  mapfile->wstats.n_darray = max(NUM_DARRAY, mapfile->wstats.thre);
  mapfile->wstats.genome->darray_all = (int *)my_calloc(mapfile->wstats.n_darray +1, sizeof(int), "wstats.genome->darray_all");
  mapfile->wstats.genome->darray_bg  = (int *)my_calloc(mapfile->wstats.n_darray +1, sizeof(int), "wstats.genome->darray_bg");
  for(chr=1; chr<g->chrnum; chr++){
    mapfile->wstats.chr[chr].darray_all = (int *)my_calloc(mapfile->wstats.n_darray +1, sizeof(int), "wstats.chr[chr]->darray_all");
    mapfile->wstats.chr[chr].darray_bg  = (int *)my_calloc(mapfile->wstats.n_darray +1, sizeof(int), "wstats.chr[chr]->darray_bg");
  }

  for(chr=1; chr<g->chrnum; chr++) makewig_chr(p, mapfile, g, chr);
  printf("done.\n");
  return;
}

static void add_wigstats_bg(Mapfile *p, TYPE_WIGARRAY *array, int chr, int binnum){
  int i, v, num=0;
  TYPE_WIGARRAY val;
  double ave=0;

  for(i=0; i<binnum; i++){
    val = array[i];
    if(!val){
      p->wstats.chr[chr].darray_bg[0]++;
      continue;
    }
    if(val >= p->wstats.num95) continue;
    v = WIGARRAY2VALUE(val);

    if(v >= p->wstats.n_darray - 1) p->wstats.chr[chr].darray_bg[p->wstats.n_darray]++;
    else p->wstats.chr[chr].darray_bg[v]++;

    ave += v;
    num++;

    if(p->wstats.chr[chr].max < val) p->wstats.chr[chr].max = val;
  }

  p->wstats.genome->ave += ave;
  ave /= num;

  /* 分散の計算 */
  double var=0;
  for(i=0; i<binnum; i++){
    if(!array[i] || array[i] >= p->wstats.num95) continue;
    v = WIGARRAY2VALUE(array[i]);
    var += (v - ave)*(v - ave);
  }
  var /= num -1;

  double mu = ave/var;
  if(mu>=1) mu = 0.9;
  if(mu<=0) mu = 0.1; 
  double n = ave * mu /(1 - mu);
  LOG("ave=%f, var=%f, p=%f, n=%f, num=%d\n", ave, var, mu, n, num);

  p->wstats.chr[chr].ave = ave;
  p->wstats.chr[chr].var = var;
  p->wstats.chr[chr].nb_p = mu;
  p->wstats.chr[chr].nb_n = n;

  /* add to ws_genome */
  for(i=0; i<=p->wstats.n_darray; i++) p->wstats.genome->darray_bg[i] += p->wstats.chr[chr].darray_bg[i];

  return;
}

static void add_wigstats_all(Mapfile *p, TYPE_WIGARRAY *array, int chr, int binnum){
  int i, v;
  TYPE_WIGARRAY val;

  for(i=0; i<binnum; i++){
    val = array[i];
    if(!val){
      p->wstats.chr[chr].darray_all[0]++;
      continue;
    }
    v = WIGARRAY2VALUE(val);

    if(v >= p->wstats.n_darray - 1) p->wstats.chr[chr].darray_all[p->wstats.n_darray]++;
    else p->wstats.chr[chr].darray_all[v]++;

    if(p->wstats.chr[chr].max < val) p->wstats.chr[chr].max = val;
  }

  // add to ws_genome
  if(p->wstats.genome->max < p->wstats.chr[chr].max) p->wstats.genome->max = p->wstats.chr[chr].max;
  for(i=0; i<=p->wstats.n_darray; i++) p->wstats.genome->darray_all[i] += p->wstats.chr[chr].darray_all[i];

  return;
}

static void makewig_chr(PwParam *p, Mapfile *mapfile, RefGenome *g, int chr){
  int i, binnum = p->binnum_chr[chr];
  TYPE_WIGARRAY *wigarray = make_wigarray(p, mapfile, g, chr, true);
  FLUSH("%s..", g->chr[chr].name);

#ifdef SHOW_WIGFILE
  int i;
  for(i=0; i<binnum; i++) fprintf(stderr,"raw: %s\t%d\t%.3f\n",g->chr[chr].name, i*p->binsize, WIGARRAY2VALUE(wigarray[i]));
#endif

  /* read mpbl file */
  double r;
  TYPE_WIGARRAY *mparray=NULL;
  char *mpfilename=NULL;
  if(p->mpfile){
    mpfilename = alloc_str_new(p->mpfile, 100);
    sprintf(mpfilename, "%s_%s_bin%d.txt", p->mpfile, g->chr[chr].name, p->binsize);
    LOG("mpfile: %s\n", mpfilename);
    mparray = read_mosaics_binfile(mpfilename, binnum, p->binsize);
  }else{
    mparray = (TYPE_WIGARRAY *)my_calloc(binnum, sizeof(TYPE_WIGARRAY), "mparray");
    TYPE_WIGARRAY d = VALUE2WIGARRAY(1);
    for(i=0; i<binnum; i++) mparray[i] = d;
  }
  /* copy wigarray to bgarray */
  int num=0;
  TYPE_WIGARRAY mpthretemp = VALUE2WIGARRAY(p->mpthre);
  TYPE_WIGARRAY *bgarray = (TYPE_WIGARRAY *)my_calloc(binnum, sizeof(TYPE_WIGARRAY), "bgarray");
  
  for(i=0; i<binnum; i++){
    if(WIGARRAY2VALUE(wigarray[i]) <= mapfile->wstats.thre && mparray[i] >= mpthretemp) bgarray[num++] = wigarray[i];
  }
  mapfile->wstats.genome->num += mapfile->wstats.chr[chr].num = num;

  /* make wigarray stats */
  add_wigstats_bg(mapfile, wigarray, chr, num);
  MYFREE(bgarray);
  add_wigstats_all(mapfile, wigarray, chr, binnum);
  
  //  if(chr==1) pw_estimate_zinb_opt(mapfile, wigarray, mparray, chr, binnum);

  /*  mpbl normalization */
  for(i=0; i<binnum; i++){
    r = WIGARRAY2VALUE(mparray[i]);
    //      printf("%d\t%f\t%d\t", i*p->binsize, r, wigarray[i]);
    if(r < p->mpthre) continue;      // mpblが低すぎる場合は補正・カウントしない
    wigarray[i] /= r;
    //      printf("%d\n", wigarray[i]);
  }

  MYFREE(mparray);
  if(p->mpfile) MYFREE(mpfilename);
  
  /* total read normalization */
  normalize_wigarray_to_rpkm(p, mapfile, &wigarray, g, chr);

#ifdef SHOW_WIGFILE
  for(i=0; i<binnum; i++) fprintf(stderr, "normalized: %s\t%d\t%.3f\n", g->chr[chr].name, i*p->binsize, WIGARRAY2VALUE(wigarray[i]));
#endif

  /* output data */
  output_bindata(p->output_dir, p->output_prefix, g, wigarray, p->gtfile, p->binsize, p->binnum_chr[chr], chr, p->wtype);

  MYFREE(wigarray);
  return;
}

static void define_s_e(PwParam *p, Mapfile *mapfile, RefGenome *g, int *s, int *e, double *w, int chr, Strand strand, int i){
  get_start_end_of_read(p, mapfile, s, e, chr, i, strand, __FUNCTION__);

  *w = INT2WEIGHT(mapfile->readarray[chr][strand].weight[i]);
  if(*w <0) fprintf(stderr, "Warning: weight of read[%s][strand%s] is %f.\n",g->chr[chr].name, str_strand[strand], *w);
  
  if(p->usereadcenter){    /* consider only center region of fragments */
    *s = (*s + *e - p->usereadcenter)/2;
    *e = (*s + *e + p->usereadcenter)/2;
  }
  if(*e >= g->chr[chr].len){  // end posion > chrlen. (*s < 0) is omitted. 
  }
  *s = max(0, *s);
  *e = min(*e, g->chr[chr].len -1);

  return;
}

static void free_mapfile_chr(PwParam *pwparam, Mapfile *p, int chr){
  Strand strand;
  for(strand=0; strand<STRANDNUM; strand++){
    MYFREE(p->readarray[chr][strand].F3);
    if(pwparam->rtype==READTYPE_PAIR) MYFREE(p->readarray[chr][strand].F5);
    MYFREE(p->readarray[chr][strand].weight);
    MYFREE(p->readarray[chr][strand].delete);
    if(pwparam->bedfilename) MYFREE(p->readarray[chr][strand].ignore);
  }
  MYFREE(p->readarray[chr]);
  return;
}

static TYPE_WIGARRAY *make_wigarray(PwParam *p, Mapfile *mapfile, RefGenome *g, int chr, bool free){
  int i,j, s, e, sbin, ebin, readnum;
  double w;
  TYPE_WIGARRAY *wigarray = (TYPE_WIGARRAY *)my_calloc(p->binnum_chr[chr], sizeof(TYPE_WIGARRAY), "wigarray");
  Strand strand;

  for(strand=0; strand<STRANDNUM; strand++){
    readnum = mapfile->chr[chr].seq[strand].n_read_infile;
    for(i=0; i<readnum; i++){
      if(mapfile->readarray[chr][strand].delete[i]) continue;
      define_s_e(p, mapfile, g, &s, &e, &w, chr, strand, i);  /* define start and end position and weight of a fragment */
      sbin = s/p->binsize;
      ebin = e/p->binsize;
      for(j=sbin; j<=ebin; j++) wigarray[j] += VALUE2WIGARRAY(w);
    }
  }

  /* free mapfile*/
  if(free) free_mapfile_chr(p, mapfile, chr);

  return wigarray;
}

static void normalize_wigarray_to_rpkm(PwParam *p, Mapfile *mapfile, TYPE_WIGARRAY **wigarrayp, RefGenome *g, int chr){
  int i;
  double w=0, dn=0, nm=0;
  TYPE_WIGARRAY *wigarray = *wigarrayp;

  switch(p->ntype){
  case NORMTYPE_NONE: return; /* skip this normalization */
  case NORMTYPE_GENOME_READ:
    dn = mapfile->genome->both.n_read_nonred;
    w = dn ? p->num4rpm/dn: 0;
    if(chr==g->chrnum-1){
      printf("\ngenomic read number = %lu, after=%d, w=%.3f\n", mapfile->genome->both.n_read_nonred, p->num4rpm, w);
      if(w>2) PRINTWARNING_W(w);
    }
    break;
  case NORMTYPE_GENOME_DEPTH:
    dn = mapfile->genome->depth;
    w = dn ? p->num4depth/dn: 0;
    if(chr==g->chrnum-1){
      printf("\ngenomic depth = %.2f, after=%.2f, w=%.3f\n", mapfile->genome->depth, p->num4depth, w);
      if(w>2) PRINTWARNING_W(w);
    }
    break;
  case NORMTYPE_CHROM_READ:
    nm = p->num4rpm * (g->chr[chr].len_mpbl/(double)g->genome->len_mpbl);
    dn = mapfile->chr[chr].both.n_read_nonred;
    w = dn ? nm/dn: 0;
    printf("read number = %lu, after=%.1f, w=%.3f\n", mapfile->chr[chr].both.n_read_nonred, nm, w);
    if(w>2) PRINTWARNING_W(w);
    break;
  case NORMTYPE_CHROM_DEPTH:
    dn = mapfile->chr[chr].depth;
    w = dn ? p->num4depth/dn: 0;
    printf("depth = %.2f, after=%.2f, w=%.3f\n", mapfile->chr[chr].depth, p->num4depth, w);
    if(w>2) PRINTWARNING_W(w);
    break;
  }

  for(i=0; i<p->binnum_chr[chr]; i++){
    if(wigarray[i]) wigarray[i] *= w;
  }

  mapfile->chr[chr].w = w;
  mapfile->chr[chr].both.n_read_rpkm = mapfile->chr[chr].both.n_read_nonred * w;
  if(p->ntype==NORMTYPE_GENOME_READ || p->ntype==NORMTYPE_GENOME_DEPTH){
    mapfile->genome->w = w;
    mapfile->genome->both.n_read_rpkm = mapfile->genome->both.n_read_nonred * w;
  }else{
    mapfile->genome->both.n_read_rpkm += mapfile->chr[chr].both.n_read_rpkm;
  }
  
  return;
}


void calc_FRiP(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int i,j, chr;
  long nread;
  Strand strand;
  struct seq *seq;
  int s,e;
  char *array = NULL;

  /* count reads */
  for(chr=1; chr<g->chrnum; chr++){
    array = makearray_inbed(&(p->enrichfile->chr[chr]), (int)g->chr[chr].len);
    if(!array) continue;

    for(strand=0; strand<STRANDNUM; strand++){
      seq = &(mapfile->chr[chr].seq[strand]);
      nread = seq->n_read_infile;
      for(i=0; i<nread; i++){
	if(p->pcrfilter && mapfile->readarray[chr][strand].delete[i]) continue;

	/* readの両端をs,eに格納 */
	get_start_end_of_read(p, mapfile, &s, &e, chr, i, strand, __FUNCTION__);
	//	printf("%s, %d, %d-%d, %d, %d, %s\n", g->chr[chr].name, i, s, e, mapfile->readarray[chr][strand].F3[i], mapfile->readarray[chr][strand].F5[i], str_strand[strand]);
	for(j=s; j<=e; j++){
	  if(array[j]){
	    mapfile->readarray[chr][strand].ignore[i] = true;
	    mapfile->chr[chr].n_read_inbed++;
	    break;
	  }
	}
      }
    }
    /* add to genome */
    mapfile->genome->n_read_inbed += mapfile->chr[chr].n_read_inbed;
    MYFREE(array);
  }

  /* calc FRiP */
  for(chr=1; chr<g->chrnum; chr++){
    if(p->pcrfilter) nread = mapfile->chr[chr].both.n_read_nonred;
    else             nread = mapfile->chr[chr].both.n_read_infile;
    mapfile->chr[chr].FRiP = mapfile->chr[chr].n_read_inbed / (double)nread;
  }
  if(p->pcrfilter) nread = mapfile->genome->both.n_read_nonred;
  else             nread = mapfile->genome->both.n_read_infile;
  mapfile->genome->FRiP  = mapfile->genome->n_read_inbed / (double)nread;
  
  return;
}

static void get_start_end_of_read(PwParam *p, Mapfile *mapfile, int *s, int *e, int chr, int i, Strand strand, const char *func){
  if(p->rtype==READTYPE_SINGLE){
    if(strand==STRAND_PLUS){
      *s = mapfile->readarray[chr][strand].F3[i];
      *e = *s + p->fraglen;
    }else{
      *e = mapfile->readarray[chr][strand].F3[i];
      *s = *e - p->fraglen;
    }
  }else{   // paired_end
    if(strand==STRAND_PLUS){
      *s = mapfile->readarray[chr][strand].F3[i];
      *e = mapfile->readarray[chr][strand].F5[i];
    }else{
      *s = mapfile->readarray[chr][strand].F5[i];
      *e = mapfile->readarray[chr][strand].F3[i];
    }
  }
  /*if(*s > *e){
    fprintf(stderr, "%s error: invalud i%d strand%s s-e (%d)-(%d)\n", func, i, str_strand[strand],*s,*e);
    fprintf(stderr, "chr%d, i%d, F3 %d, F5 %d, strand %d\n", chr,i,mapfile->readarray[chr][strand].F3[i], mapfile->readarray[chr][strand].F5[i], strand);
    exit(0);
    }*/
  return;
}
