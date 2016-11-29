/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "drompa_param_new.h"
#include "drompa_readfile.h"
#include "dd_init.h"
#include "dd_readannotation.h"
#include "dp_call.h"
#include "dd_stroke.h"
#include "dd_profile.h"
#include "dd_heatmap.h"
#include "dd_otherfunc.h"
#include "stringp.h"
#include "alloc.h"
#include "macro.h"
#include "filehandle.h"
#include "readfile.h"

static void drompafunc(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g);
static void drompa4chr(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g, int chr);
static void mergepdf_all(DrParam *p, DDParam *d, RefGenome *g);
static DDParam *ddparam_new();
static void ddparam_delete(DDParam *p);

int main(int argc, char **argv){
#ifdef CLOCK
  clock_t t1 = clock();
#endif

  DrParam   *drparam = drparam_new();
  DDParam   *ddparam = ddparam_new();
  RefGenome *refgenome = refgenome_new();
  SamplePair *sample=NULL;
  dd_argv_init(argc, argv, drparam, ddparam, &sample, refgenome); /* sample構造体は中で確保される */

  /* calculate the global parameter */
  dr_calc_global_param(drparam, sample, refgenome);
  ddparam->pageheight = calc_pageheight(drparam, ddparam);

#ifdef CLOCK

  clock_t t2 = clock();
  print_time(t1,t2);
#endif

  drompafunc(drparam, ddparam, sample, refgenome);

#ifdef CLOCK
  clock_t t3 = clock();
  print_time(t2,t3);
#endif

  SamplePair_delete(sample, drparam->samplenum);
  ddparam_delete(ddparam);
  drparam_delete(drparam);

  /*--- dbg_malloc (for degug) ---*/
#ifdef MEMTEST
  dbg_print_alloc_count(stdout);
  //  dbg_print_all_alloc_block(stdout);
  dbg_print_memory_max(stdout);
#endif
  return 0;
}

static void drompafunc(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g){
  if(p->ftype == FTYPE_GOVERLOOK) genome_overlook(p, d, g);
  else if(p->ftype==FTYPE_PROFILE){
    d->filename_profile = alloc_str_new(p->headname, 10);
    sprintf(d->filename_profile, "%s.R", p->headname);
    remove_file(d->filename_profile);
  }else if(p->ftype==FTYPE_HEATMAP){
    draw_heatmap(p, d, sample, g);
    return;
  }else if(p->ftype==FTYPE_FRIP){
    dd_counttags(p, d, g, sample);
    return;
  }else if(p->ftype==FTYPE_COMPARE_INTENSITY){
    dd_compare_intensity(p, d, g, sample);
  }else if(p->ftype==FTYPE_COMPARE_GENEBODY){
    dd_compare_genebody(p, d, g, sample);
  }
  if(p->ftype==FTYPE_TR){
    dd_travelling_ratio(p, d, g, sample);
    return;
  }

  int chr;
  for(chr=1; chr<g->chrnum; chr++) drompa4chr(p, d, sample, g, chr);
  
  if(d->makefig && !d->png) mergepdf_all(p, d, g);
  if(p->ftype==FTYPE_PROFILE && !d->pdetail) show_profile(p, d, sample);
  return;
}


static void read_annotation(DrParam *p, DDParam *d, RefGenome *g, int chr){
  if(d->GC.argv) read_graph(&(d->GC), g, chr, "GC%",   GC_MMIN, GC_MMAX, GTYPE_SINGLE);
  if(d->GD.argv) read_graph(&(d->GD), g, chr, "genes", GD_MMIN, GD_MMAX, GTYPE_SINGLE);
  if(p->ftype==FTYPE_GV) return;

  if(d->gene.argv || d->arsfile || d->terfile) d->gene.gene = (Gene *)my_calloc(STRUCT_GENE_MAX, sizeof(Gene), "gene");
  if(d->gene.argv)   read_gene(&(d->gene), g, chr, d->gftype);
  if(d->arsfile)     read_ARS_OriDB(d->arsfile, &(d->gene), chr);
  if(d->terfile)     read_TER_scer(d->terfile, &(d->gene), chr);
  if(d->repeat.argv) read_repeat(d, g, chr);

  if(d->genefile_argv) read_genefile(d, d->drawregion, g, chr);
  return;
}

static char *make_peakarray(Peak *peak, int binsize, int binnum, int chr){
  int i,j,s,e;
  char *array = (char *)my_calloc(binnum, sizeof(char), "char");
  for(i=0; i<peak->num; i++){
    if(peak->bs[i].chr != chr) continue;
    s = peak->bs[i].start / binsize;
    e = peak->bs[i].end   / binsize;
    for(j=s; j<=e; j++) array[j] = 1;
  }
  return array;
}

static void drompa4chr(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g, int chr){
  int i;
  char *prefix=NULL;
  char *filename=NULL;

  if(!p->includeYM && (!strcmp(g->chr[chr].name, "chrY") || !strcmp(g->chr[chr].name, "chrM") || !strcmp(g->chr[chr].name, "chrMT"))) return;
  if(d->chronly && d->chronly != chr) return;

  if(d->drawregion_argv && !d->drawregion->chr[chr].num) return;
  /* read annotation files */
  read_annotation(p, d, g, chr);
  //  printf("num=%d\n", d->drawregion->chr[chr].num);
  if(d->genefile_argv && !d->drawregion->chr[chr].num) return;
  
  FLUSH("%s: ", g->chr[chr].name);  
  
  if(p->ftype==FTYPE_PD){
    for(i=0; i<d->pdnum; i++){
      if(d->pd_prop) read_graph(&(d->PD[i]), g, chr, d->PD[i].name, 0, 0, GTYPE_MULTI_PROP);
      else           read_graph(&(d->PD[i]), g, chr, d->PD[i].name, 0, 0, GTYPE_MULTI);
    }
  }else{
    /* read wigfile */
    for(i=0; i<p->samplenum; i++){
      //      printf("sample%d: reading ChIP file...", i+1);
      FLUSH("sample%d..", i+1);
      if(sample[i].copyC==-1) dr_read_wigdata(p, &(sample[i]), sample[i].ChIP, g, chr);
      if(sample[i].Input->argv && sample[i].copyI==-1){
	//	printf("reading Input file...\n");
	dr_read_wigdata(p, &(sample[i]), sample[i].Input, g, chr);
      }//else printf("\n");
      if(sample[i].peak_argv) sample[i].peakarray = make_peakarray(sample[i].peak, sample->binsize, sample->binnum[chr], chr);
    }
    printf("\n");
  }
  if(p->ftype==FTYPE_PROFILE) add_profile(p, d, g, sample, chr);
  if(d->makefig){
    /* read gapfile */
    if(p->gapfile){
      filename = alloc_str_new(p->gapfile, 100);
      sprintf(filename, "%s_%s_bin%d.txt", p->gapfile, g->chr[chr].name, sample->binsize);
      LOG("gapfile: %s\n", filename);
      d->gaparray = read_mosaics_binfile(filename, sample->binnum[chr], sample->binsize);
      MYFREE(filename);
    }
    if(p->mpfile){
      filename = alloc_str_new(p->mpfile, 100);
      sprintf(filename, "%s_%s_bin%d.txt", p->mpfile, g->chr[chr].name, sample->binsize);
      LOG("mpfile: %s\n", filename);
      d->mparray = read_mosaics_binfile(filename, sample->binnum[chr], sample->binsize);
      MYFREE(filename);
    }
    prefix = alloc_str_new(p->headname, 100);
    printf("drawing...\n");
    sprintf(prefix, "%s_%s", p->headname, g->chr[chr].name);
    draw(p, d, sample, g, chr, prefix);
    sprintf(d->command_mergepdf, "%s%s.pdf ", d->command_mergepdf, prefix);
    MYFREE(prefix);
  }
  
  /* delete wigfile */
  if(p->ftype!=FTYPE_PD){
    for(i=0; i<p->samplenum; i++){
      if(sample[i].copyC==-1) dr_delete_wigdata(sample[i].ChIP);
      if(sample[i].Input->argv && sample[i].copyI==-1) dr_delete_wigdata(sample[i].Input);
      if(sample[i].peak_argv) MYFREE(sample[i].peakarray);
    }
  }
  
  return;
}

static void mergepdf_all(DrParam *p, DDParam *d, RefGenome *g){
  int chr;
  char *filename = alloc_str_new(p->headname, 20);
  printf("\nmerging pdf files...\n");

  sprintf(filename, "%s.pdf", p->headname);
  remove_file(filename);
  /* cpdf */
  sprintf(d->command_mergepdf, "%s-o %s.pdf", d->command_mergepdf, p->headname);
  my_system(d->command_mergepdf);
  /* rm */
  if(d->rmchr){
    for(chr=1; chr<g->chrnum; chr++){
      sprintf(filename, "%s_%s.pdf", p->headname, g->chr[chr].name);
      remove(filename);
    }
  }
  MYFREE(filename);
  return;
}

static DDParam *ddparam_new(){
  DDParam *p = (DDParam *)my_calloc(1, sizeof(DDParam), "DDParam");
  p->GC.wsize = GCSIZE_DEFAULT;
  p->GD.wsize = GDSIZE_DEFAULT;
  p->pdsize   = PDSIZE_DEFAULT;
  p->linenum_per_page = 1;
  p->barnum = BARNUM_DEFAULT;
  p->ystep = YSTEP_DEFAULT;
  p->visualize_ctag = 1;
  p->backcolors  = 1;
  p->stroke_ymem = 1;
  p->stroke_ylab = 1;
  p->width_per_line = LS_DEFAULT;
  p->ntype = 1;

  /* readfile */
  p->genefile_len = 50000; // 50 kbp

  /* profile and heatmap */
  p->showse = true;
  p->compwidth = 2500;
  p->ptype = TSS;

  p->tssthre = 1;

  return p;
}

static void ddparam_delete(DDParam *p){
  MYFREE(p->command_mergepdf);
  MYFREE(p->PD);
  MYFREE(p);
  return;
}
