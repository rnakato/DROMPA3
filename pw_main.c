/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "readfile.h"
#include "filehandle.h"
#include "pw_filtering.h"
#include "pw_param_new.h"
#include "pw_init.h"
#include "pw_readmapfile.h"
#include "pw_makefile.h"
#include "pw_gc.h"
#include "stringp.h"
#include "macro.h"
#include "pw_estimate.h"
#include "pw_ccp.h"

static void add_read_red_to_genome(Mapfile *mapfile, RefGenome *g, bool);
static void calc_depth(PwParam *p, Mapfile *mapfile, RefGenome *g);
static void output_wigstats(PwParam *p, Mapfile *mapfile, RefGenome *g);
static void output_stats(PwParam *p, Mapfile *mapfile, RefGenome *g);

int main(int argc, char *argv[]){
#ifdef CLOCK
  clock_t t1 = clock();
#endif

  PwParam *p = pwparam_new();
  RefGenome *g = refgenome_new();
  pw_argv_init(argc, argv, p, g);
  
#ifdef DEBUG
  printf("\ngenome, len:%ld, len_mpbl:%ld, mpbl:%f\n", g->genome->len, g->genome->len_mpbl, g->genome->p_mpbl);
  int chr;
  for(chr=0; chr<g->chrnum; chr++){
    printf("i%d, %s\tlen:%ld, len_mpbl:%ld, mpbl:%f\n", chr, g->chr[chr].name, g->chr[chr].len, g->chr[chr].len_mpbl, g->chr[chr].p_mpbl);
  }
  print_sizeof(Mapfile, "Mapfile");
  print_sizeof(SeqStats, "SeqStats");
#endif

#ifdef CLOCK
  clock_t t2 = clock();
  print_time(t1,t2);
#endif
  /* read mapfile */
  Mapfile *mapfile = read_mapfile(p, g);
  
#ifdef CLOCK
  clock_t t3 = clock();
  print_time(t2,t3);
#endif

  if(!mapfile->genome->both.n_read_infile) {
    fprintf(stderr, "ERROR: there is no read in input file.\n");
    exit(0);
  }
  
  /* PCR bias filtering and ignore enrichregions */
  if(p->pcrfilter) check_redundant_reads(p, mapfile, g);
  add_read_red_to_genome(mapfile, g, p->pcrfilter);

  if(p->ccp) pw_ccp(p, mapfile, g->chrmax, (int)g->chr[g->chrmax].len);

#ifdef CLOCK
  clock_t t4 = clock();
  print_time(t3,t4);
#endif

  if(p->bedfilename) calc_FRiP(p, mapfile, g);

  /* calculate depth */
  calc_depth(p, mapfile, g);

#ifdef CLOCK
  clock_t t5 = clock();
  print_time(t4,t5);
#endif

  /* GC normalization */
  if(p->genomefile) GCnorm(p, mapfile, g);

  /* make and output wigdata */
  makewig(p, mapfile, g);

  /* output wigarray_stats, 
     calculate genome coverage */
  output_wigstats(p, mapfile, g);

#ifdef CLOCK
  clock_t t6 = clock();
  print_time(t5,t6);
#endif

  /* output stats */
  output_stats(p, mapfile, g);

  /* End processing */
  mapfile_delete(mapfile, g->chrnum);
  pwparam_delete(p, g);
  refgenome_delete(g);
  
  /*--- dbg_malloc (for degug) ---*/
#ifdef MEMTEST
  dbg_print_alloc_count(stdout);
  //  dbg_print_all_alloc_block(stdout);
  dbg_print_memory_max(stdout);
#endif
  return 0;
}

static void output_wigstats(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int i,chr;
  char *outputfile = alloc_str_new(p->output_dir, strlen(p->output_prefix) +100);
  sprintf(outputfile, "%s/%s.%d.binarray_dist.xls", p->output_dir, p->output_prefix, p->binsize);
  FILE *OUT = my_fopen(outputfile, FILE_MODE_WRITE);

  pw_estimate_nb_param(mapfile, g);
  /* for Poisson */
  mapfile->wstats.genome->ave /= mapfile->wstats.genome->num;

  fprintf(OUT, "Poisson: lambda = %f\n", mapfile->wstats.genome->ave);
  fprintf(OUT, "Negative binomial: p=%f, n=%f, p0=%f\n", mapfile->wstats.genome->nb_p, mapfile->wstats.genome->nb_n, mapfile->wstats.genome->nb_p0);
  fprintf(OUT, "<genome>\n");
  fprintf(OUT, "read number\tAll regions\tprop all\tBG regions\tprop bg\tNB simulated (%s)\tNB simulated (genome)\tPoisson simulated\n", g->chr[g->chrmax].name);
  for(i=0; i<mapfile->wstats.n_darray; i++){
    fprintf(OUT, "%d\t%d\t%f\t", i, mapfile->wstats.genome->darray_all[i], mapfile->wstats.genome->darray_all[i]/(double)p->binnum_genome);
    fprintf(OUT, "%d\t%f\t",        mapfile->wstats.genome->darray_bg[i],  mapfile->wstats.genome->darray_bg[i] /(double)mapfile->wstats.genome->num);
    fprintf(OUT, "%f\t", pw_get_negative_binomial(i, mapfile->wstats.chr[g->chrmax].nb_p, mapfile->wstats.chr[g->chrmax].nb_n));
    fprintf(OUT, "%f\t", pw_get_zeroinflated_negative_binomial(i, mapfile->wstats.genome->nb_p, mapfile->wstats.genome->nb_n, mapfile->wstats.genome->nb_p0));
    fprintf(OUT, "%f\n", pw_get_poisson(i, mapfile->wstats.genome->ave));
  }

  /* calculate genome coverage */
  for(chr=1; chr<g->chrnum; chr++){
    int s=0;
    for(i=0;i<=mapfile->wstats.n_darray;i++){
      //      printf("%d, %d\n",i,mapfile->wstats.chr[chr].darray_bg[i]);
      s += mapfile->wstats.chr[chr].darray_bg[i];
    }
    g->chr[chr].gcov = (mapfile->wstats.chr[chr].num - mapfile->wstats.chr[chr].darray_bg[0])/(double)mapfile->wstats.chr[chr].num;
    //    printf("chr=%d, s=%d, mapfile->wstats.chr[chr].num=%d, 0=%d, gcov=%f (%d/%d)\n",chr,s,mapfile->wstats.chr[chr].num, mapfile->wstats.chr[chr].darray_all[0], g->chr[chr].gcov, (mapfile->wstats.chr[chr].num - mapfile->wstats.chr[chr].darray_bg[0]), mapfile->wstats.chr[chr].num);
  }
  g->genome->gcov = (mapfile->wstats.genome->num - mapfile->wstats.genome->darray_bg[0])/(double)mapfile->wstats.genome->num;

  fclose(OUT);
  return;
}

static void add_read_red_to_genome(Mapfile *mapfile, RefGenome *refgenome, bool pcrfilter){
  int chr;
  Strand strand;
  for(chr=1; chr<refgenome->chrnum; chr++){
    for(strand=0; strand<STRANDNUM; strand++){

      if(!pcrfilter) mapfile->chr[chr].seq[strand].n_read_nonred = mapfile->chr[chr].seq[strand].n_read_infile;
      
      mapfile->chr[chr].both.n_read_nonred += mapfile->chr[chr].seq[strand].n_read_nonred;
      mapfile->chr[chr].both.n_read_red    += mapfile->chr[chr].seq[strand].n_read_red;
      mapfile->genome->seq[strand].n_read_nonred += mapfile->chr[chr].seq[strand].n_read_nonred;
      mapfile->genome->seq[strand].n_read_red    += mapfile->chr[chr].seq[strand].n_read_red;
      LOG("%s: chr%d%s: n_read_infile: %ld\tn_readname: %.1Lf\tn_nonred: %ld\tn_red: %ld\n",__FUNCTION__,chr, str_strand[strand], mapfile->chr[chr].seq[strand].n_read_infile, mapfile->chr[chr].seq[strand].n_readname, mapfile->chr[chr].seq[strand].n_read_nonred, mapfile->chr[chr].seq[strand].n_read_red);
    }
    mapfile->genome->both.n_read_nonred += mapfile->chr[chr].both.n_read_nonred;
    mapfile->genome->both.n_read_red += mapfile->chr[chr].both.n_read_red;
      LOG("%s: chr%dboth: n_read_infile: %ld\tn_readname: %.1Lf\tn_nonred: %ld\tn_red: %ld\n",__FUNCTION__,chr, mapfile->chr[chr].both.n_read_infile, mapfile->chr[chr].both.n_readname, mapfile->chr[chr].both.n_read_nonred, mapfile->chr[chr].both.n_read_red);
  }
  LOG("%s: genomeboth: n_read_infile: %ld\tn_readname: %.1Lf\tn_nonred: %ld\tn_red: %ld\n",__FUNCTION__,mapfile->genome->both.n_read_infile, mapfile->genome->both.n_readname, mapfile->genome->both.n_read_nonred, mapfile->genome->both.n_read_red);
}

static void calc_depth(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int chr;
  for(chr=1; chr<g->chrnum; chr++){
    mapfile->chr[chr].depth = g->chr[chr].len_mpbl ? mapfile->chr[chr].both.n_read_nonred*p->fraglen/(double)g->chr[chr].len_mpbl: 0;
    LOG("%s: chr%dboth: depth: %f\n",__FUNCTION__,chr, mapfile->chr[chr].depth);
  }
  mapfile->genome->depth = g->genome->len_mpbl ? mapfile->genome->both.n_read_nonred*p->fraglen/(double)g->genome->len_mpbl: 0;
  LOG("%s: genomeboth: depth: %f\n",__FUNCTION__,mapfile->genome->depth);
  
  return;
}

static void print_SeqStats(FILE *OUT, PwParam *pwparam, SeqStats *p, long nread, Fastadata *fd){
  /* genome data */
  fprintf(OUT, "%s\t", fd->name);
  fprintf(OUT, "%ld\t%ld\t%.3f\t", fd->len, fd->len_mpbl, fd->p_mpbl);
  /* total reads*/  fprintf(OUT, "%ld\t%ld\t%ld\t%.1f%%\t", p->both.n_read_infile, p->seq[STRAND_PLUS].n_read_infile, p->seq[STRAND_MINUS].n_read_infile, p->both.n_read_infile*100/(double)nread);
  /* nonredundant reads */
  fprintf(OUT, "%ld (%.1f%%)\t", p->both.n_read_nonred, p->both.n_read_nonred*100/(double)p->both.n_read_infile);
  fprintf(OUT, "%ld (%.1f%%)\t", p->seq[STRAND_PLUS].n_read_nonred,  p->seq[STRAND_PLUS].n_read_nonred*100/(double)p->seq[STRAND_PLUS].n_read_infile);
  fprintf(OUT, "%ld (%.1f%%)\t", p->seq[STRAND_MINUS].n_read_nonred, p->seq[STRAND_MINUS].n_read_nonred*100/(double)p->seq[STRAND_MINUS].n_read_infile);
  /* redundant reads */
  fprintf(OUT, "%ld (%.1f%%)\t", p->both.n_read_red, p->both.n_read_red*100/(double)p->both.n_read_infile);
  fprintf(OUT, "%ld (%.1f%%)\t", p->seq[STRAND_PLUS].n_read_red, p->seq[STRAND_PLUS].n_read_red*100/(double)p->seq[STRAND_PLUS].n_read_infile);
  fprintf(OUT, "%ld (%.1f%%)\t", p->seq[STRAND_MINUS].n_read_red, p->seq[STRAND_MINUS].n_read_red*100/(double)p->seq[STRAND_MINUS].n_read_infile);
  /* reads after GCnorm */
  if(pwparam->genomefile){
    fprintf(OUT, "%.1f (%.2f)\t", p->both.n_read_afterGC, p->both.n_read_afterGC/(double)p->both.n_read_infile);
    fprintf(OUT, "%.1f (%.2f)\t", p->seq[STRAND_PLUS].n_read_afterGC, p->seq[STRAND_PLUS].n_read_afterGC/(double)p->seq[STRAND_PLUS].n_read_infile);
    fprintf(OUT, "%.1f (%.2f)\t", p->seq[STRAND_MINUS].n_read_afterGC, p->seq[STRAND_MINUS].n_read_afterGC/(double)p->seq[STRAND_MINUS].n_read_infile);
  }
  /* read depth */
  fprintf(OUT, "%.2f\t", p->depth);
  /* weight */
  if(p->w) fprintf(OUT, "%.2f\t", p->w); else fprintf(OUT, " - \t");
  /* read number after rpkm normalization */
  if(pwparam->ntype==NORMTYPE_NONE) fprintf(OUT, "%ld\t", p->both.n_read_nonred); else fprintf(OUT, "%ld\t", p->both.n_read_rpkm);
  /* genome coverage */
  fprintf(OUT, "%.2f\t", fd->gcov);
  /* FRiP */
  if(pwparam->bedfilename) fprintf(OUT, "%.4f\t", p->FRiP);

  return;
}

static void output_stats(PwParam *p, Mapfile *mapfile, RefGenome *g){
  int chr;
  char strtemp1[MAXLEN_INSCOMMA], strtemp2[MAXLEN_INSCOMMA];

  char *filename = alloc_str_new(p->output_dir, strlen(p->output_prefix) +100);
  sprintf(filename, "%s/%s.%d.xls", p->output_dir, p->output_prefix, p->binsize);
  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);

  fprintf(OUT, "parse2wig version %s\n", VERSION);
  fprintf(OUT, "Input file:\t\"%s\"\n", p->inputfile);
  if(p->pcrfilter) fprintf(OUT, "Redundancy threshold: >%d\n", mapfile->threshold4filtering);
  /* complexity */
  insComma(mapfile->cs_raw.nt_nonred, strtemp1);
  insComma(mapfile->cs_raw.nt_all, strtemp2);

  if(p->pcrfilter) {
    if(mapfile->cs_raw.tv) fprintf(OUT, "Library complexity: (%.3f) (%s / %s)\n", mapfile->cs_raw.complexity, strtemp1, strtemp2);
    else fprintf(OUT, "Library complexity: %.3f (%s / %s)\n", mapfile->cs_raw.complexity, strtemp1, strtemp2);
  }
  if(p->genomefile) fprintf(OUT, "GC summit: %d\n", mapfile->maxGC);
  fprintf(OUT, "Poisson: lambda = %f\n", mapfile->wstats.genome->ave);
  fprintf(OUT, "Negative binomial: p=%f, n=%f, p0=%f\n", mapfile->wstats.genome->nb_p, mapfile->wstats.genome->nb_n, mapfile->wstats.genome->nb_p0);

  /* SeqStats */
  fprintf(OUT, "\tlength\tmappable base\tmappability\t");
  fprintf(OUT, "total reads\t\t\t\t");
  fprintf(OUT, "nonredundant reads\t\t\t");
  fprintf(OUT, "redundant reads\t\t\t");
  if(p->genomefile) fprintf(OUT, "reads (GCnormed)\t\t\t");
  fprintf(OUT, "read depth\t");
  fprintf(OUT, "scaling weight\t");
  fprintf(OUT, "normalized read number\t");
  fprintf(OUT, "genome coverage\t");
  if(p->bedfilename) fprintf(OUT, "FRiP\t");
  fprintf(OUT, "\n");
  fprintf(OUT, "\t\t\t\t");
  fprintf(OUT, "both\tforward\treverse\t %% genome\t");
  fprintf(OUT, "both\tforward\treverse\t");
  fprintf(OUT, "both\tforward\treverse\t");
  if(p->genomefile) fprintf(OUT, "both\tforward\treverse\t");
  fprintf(OUT, "\n");

  /*** genome ***/
  print_SeqStats(OUT, p, mapfile->genome, mapfile->genome->both.n_read_infile, g->genome);
  fprintf(OUT, "\n");
  /*** chromosome ***/
  for(chr=1; chr<g->chrnum; chr++){
    print_SeqStats(OUT, p, &(mapfile->chr[chr]), mapfile->genome->both.n_read_infile, &(g->chr[chr]));
    fprintf(OUT, "\n");
  }

  fclose(OUT);
  printf("stats is output in %s.\n", filename);
  MYFREE(filename);

  return;
}

