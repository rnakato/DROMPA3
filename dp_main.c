/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "drompa_gv.h"
#include "drompa_param_new.h"
#include "drompa_readfile.h"
#include "dp_init.h"
#include "dp_call.h"
#include "stringp.h"
#include "filehandle.h"
#include "alloc.h"
#include "macro.h"

static void drompa4chr(DrParam *p, SamplePair *sample, RefGenome *g, int chr);
static void OutputPeaks(DrParam *p, SamplePair *sample, RefGenome *g);

/* enrichment, pvalueのwigファイルを出力できるように拡張*/

int main(int argc, char **argv){

#ifdef CLOCK
  clock_t t1 = clock();
#endif

  DrParam *p = drparam_new();
  RefGenome *g = refgenome_new();
  SamplePair *sample = NULL;
  dp_argv_init(argc, argv, p, &sample, g); /* sample構造体は中で確保される */

  /* calculate the global parameter */
  dr_calc_global_param(p, sample, g);

#ifdef CLOCK
  clock_t t2 = clock();
  print_time(t1,t2);
#endif

  /* for chromosome */
  FLUSH("\nPeak-calling: ");
  int chr;
  for(chr=1; chr<g->chrnum; chr++) drompa4chr(p, sample, g, chr);
  printf("done.\n");

#ifdef CLOCK
  clock_t t3 = clock();
  print_time(t2,t3);
#endif

  printf("output peaklist...");
  OutputPeaks(p, sample, g);
  printf("done.\n");

  SamplePair_delete(sample, 1);
  drparam_delete(p);

  /*--- dbg_malloc (for degug) ---*/
#ifdef MEMTEST
  dbg_print_alloc_count(stdout);
  //  dbg_print_all_alloc_block(stdout);
  dbg_print_memory_max(stdout);
#endif
  return 0;
}

static void drompa4chr(DrParam *p, SamplePair *sample, RefGenome *g, int chr){
  if(!p->includeYM && (!strcmp(g->chr[chr].name, "chrY") || !strcmp(g->chr[chr].name, "chrM") || !strcmp(g->chr[chr].name, "chrMT"))) return;

  FLUSH("%s..", g->chr[chr].name);  
  dr_read_wigdata(p, sample, sample->ChIP, g, chr);
  if(sample->Input->argv){
    dr_read_wigdata(p, sample, sample->Input, g, chr);
  }
  dp_call(p, sample, g, chr);

  dr_delete_wigdata(sample->ChIP);
  if(sample->Input->argv) dr_delete_wigdata(sample->ChIP);
  return;
}

static int comp(const void *c1, const void *c2){
  struct bs n1 = *(struct bs *)c1;
  struct bs n2 = *(struct bs *)c2;
  if(n1.p_enr > n2.p_enr) return -1;
  else if(n1.p_enr == n2.p_enr) return 0;
  else return 1;
}

static void calc_qvalue_BH(Peak **peak_ref, char *Input_argv){
  int i;
  Peak *peak = *peak_ref;
  qsort(peak->bs, peak->num, sizeof(struct bs), comp);
  for(i=0; i<peak->num; i++){
    if(Input_argv){
      peak->bs[i].qvalue = pow(10, -peak->bs[i].p_enr)   * peak->num /(double)(i+1);
      //      printf("%f %f %f %f\n",peak->bs[i].qvalue, peak->bs[i].p_enr, pow(10, -peak->bs[i].p_enr) , peak->num /(double)(i+1));
    }else           peak->bs[i].qvalue = pow(10, -peak->bs[i].p_inter) * peak->num /(double)(i+1);
  }

  double min = peak->bs[peak->num-1].qvalue;
  for(i=peak->num-2; i>=0; i--){
    if(peak->bs[i].qvalue > min) peak->bs[i].qvalue = min;
    if(peak->bs[i].qvalue < min) min = peak->bs[i].qvalue;
  }
  *peak_ref = peak;
  return;
}

/* 正しいbed形式で出力できるオプションを作る */
static void OutputPeaks(DrParam *p, SamplePair *sample, RefGenome *g){
  int i, s, e, maxposi;
  struct bs *bs = NULL;
  calc_qvalue_BH(&(sample->peak), sample->Input->argv);

  char *filename = alloc_str_new(p->headname, 10);
  char *filename_bed = alloc_str_new(p->headname, 10);
  sprintf(filename, "%s.xls", p->headname);
  sprintf(filename_bed, "%s.bed", p->headname);
  LOG("outputfilename:%s, %s\n", filename, filename_bed);
  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);
  FILE *OUT_BED = my_fopen(filename_bed, FILE_MODE_WRITE);

  fprintf(OUT, "chromosome\tstart\tend\tpeak summit\tpeak width\tmax IP\t-log10(pvalue (IP internal))\t");
  if(sample->Input->argv) fprintf(OUT, "-log10(pvalue IP/Input)\t");
  fprintf(OUT, "FDR\n");

  int id4bed=1;
  for(i=0, bs=sample->peak->bs; i<sample->peak->num; i++, bs++){
    if(bs->qvalue < p->qthre){
      s = bs->start    * sample->binsize;
      e = (bs->end +1) * sample->binsize -1;
      maxposi = bs->maxposi * sample->binsize + sample->binsize /2;
      fprintf(OUT, "%s\t%d\t%d\t%d\t%d\t%.1f\t%.2f\t", g->chr[bs->chr].name, s, e, maxposi, e-s, bs->maxIP, bs->p_inter);
      if(sample->Input->argv) fprintf(OUT, "%.2f\t%.5f\n", bs->p_enr, bs->qvalue);
      else fprintf(OUT, "%.5f\n", bs->qvalue);
      fprintf(OUT_BED, "%s\t%d\t%d\t%s_%d\t", g->chr[bs->chr].name, s, e, p->headname, id4bed++);
      if(sample->Input->argv) fprintf(OUT_BED, "%.2f\n", bs->p_enr);
      else fprintf(OUT_BED, "%.5f\n", bs->p_inter);
    }
  }
  fclose(OUT);
  fclose(OUT_BED);
  MYFREE(filename);
  MYFREE(filename_bed);
}
