/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "dp_init.h"
#include "drompa_readfile.h"
#include "readfile.h"
#include "alloc.h"
#include "argument.h"
#include "common.h"
#include "macro.h"
#include "drompa_usage.h"
#include "stringp.h"

#define PRINT_ERROR(args...) do{ printf(args); goto err; }while(0)

static void init_dump(DrParam *p, SamplePair *sample);

static void print_version(){
  printf("drompa_peakcall version %s\n", VERSION);
  exit(0);
}
static void print_usage_base(){
  fprintf(stderr, "Usage: drompa_peakcall [--version] <command>\n");
  fprintf(stderr, "Command: PC_SHARP    peak-calling (for sharp mode)\n");
  fprintf(stderr, "         PC_BROAD    peak-calling (for broad mode)\n");
  fprintf(stderr, "         PC_E        peak-calling (enrichment ratio)\n");
  exit(0);
}

static void check_argv1(DrParam *p, char *argv, int *binsize){
  if(!strcmp(argv, "-h") || !strcmp(argv, "--help")) print_usage_base();
  if(!strcmp(argv, "--version")) print_version();
  else if(!strcmp(argv, "PC_SHARP")){
    p->ftype = FTYPE_PEAKCALL_SHARP;
  }else if(!strcmp(argv, "PC_BROAD")){
    p->ftype = FTYPE_PEAKCALL_BROAD;
    *binsize = BINSIZE_BROAD;
    p->smoothing = SMOOTHING_BROAD;
  }else if(!strcmp(argv, "PC_ENRICH")){
    p->ftype = FTYPE_PEAKCALL_E;
    p->enrichthre = ENRICHTHRE_DEFAULT;
  }else{
    fprintf(stderr, "Please specify command.\n\n");
    print_usage_base();
  }
  return;
}

void dp_argv_init(int argc, char **argv, DrParam *p, SamplePair **sample, RefGenome *g){
  int i, chr, nline;
  Elem clm[ELEM_NUM];
  if(argc<=1) print_usage_base();

  int binsize = BINSIZE_DEFAULT;
  check_argv1(p, argv[1], &binsize);
  argv++; argc--;

  char **str_argv = (char **)my_calloc(1024, sizeof(char *), "str_argv");
  char **str_bed = (char **)my_calloc(1024, sizeof(char *), "str_bed");
  const Argument args[]={
    {"--version", ARGUMENT_TYPE_FUNCTION, print_version,        NULL},
    {"-h"    ,    ARGUMENT_TYPE_FUNCTION, print_usage_dp,       NULL},
    {"--help",    ARGUMENT_TYPE_FUNCTION, print_usage_dp,       NULL},
    {"-gt",       ARGUMENT_TYPE_STRING,   &p->gtfile,           NULL},
    {"-p",        ARGUMENT_TYPE_STRING,   &p->headname,         NULL},
    {"-i",        ARGUMENT_TYPE_STRING_MULTI, &(str_argv), &(p->samplenum)},
    {"-binsize",  ARGUMENT_TYPE_INTEGAR,  &binsize,             NULL},
    {"-mp",       ARGUMENT_TYPE_STRING,   &p->mpfile,           NULL},
    {"-mpthre",   ARGUMENT_TYPE_FLOAT,    &p->mpthre,           NULL},
    {"-includeYM",ARGUMENT_TYPE_FLAG_ON,  &p->includeYM,        NULL},
    {"-if",       ARGUMENT_TYPE_INTEGAR,  &p->itype,            NULL},
    {"-sm",       ARGUMENT_TYPE_INTEGAR,  &p->smoothing,        NULL},
    {"-width4lmd",ARGUMENT_TYPE_INTEGAR,  &p->width4lambda,     NULL},
    {"-ipm",      ARGUMENT_TYPE_FLOAT,    &p->IPmaxthre,        NULL},
    {"-ethre",    ARGUMENT_TYPE_FLOAT,    &p->enrichthre,       NULL},
    {"-pthre_internal",ARGUMENT_TYPE_LOG10, &p->pthre_internal, NULL},
    {"-pthre_enrich",  ARGUMENT_TYPE_LOG10, &p->pthre_enrich,   NULL},
    {"-qthre",    ARGUMENT_TYPE_FLOAT,    &p->qthre,            NULL},
    {"-norm",     ARGUMENT_TYPE_INTEGAR,  &p->ntype,            NULL},
    {"-outputwig",ARGUMENT_TYPE_INTEGAR,  &p->outputwig,        NULL},
    {"-owtype",   ARGUMENT_TYPE_INTEGAR,  &p->wtype,            NULL},
    {"-odir",     ARGUMENT_TYPE_STRING,   &p->output_dir,       NULL},
    {"-ignore",   ARGUMENT_TYPE_STRING_MULTI, &(str_bed), &(p->n_igregion)},
    {NULL,        ARGUMENT_TYPE_NONE,     NULL,                 NULL},
  };
  argument_read(&argc, argv, args);

  /* checkparam */
  if(print_error_peakcall(p,g)) goto err;
  if(p->outputwig){
    if(!p->output_dir) p->output_dir = strdup("drompadir");
    mkdir(p->output_dir, 0755);
  }
  if(!p->samplenum) PRINT_ERROR("please specify sample file.\n");
  if(p->samplenum>1) PRINT_ERROR("drompa_peakcall accepts only one ChIP-input pair.\n");

  (*sample) = scan_samplestr(p, str_argv, NULL, g->chrnum);
  (*sample)->binsize = binsize;
  for(chr=1; chr<g->chrnum; chr++){
    (*sample)->binnum[chr] = g->chr[chr].len/(*sample)->binsize +1;
    if(g->chr[chr].len < (*sample)->binsize) fprintf(stderr, "Warning: length of %s (%ld) is shorter than binsize (%d).\n", g->chr[chr].name, g->chr[chr].len, (*sample)->binsize);
  }
  if(p->ftype == FTYPE_PEAKCALL_E && !(*sample)->Input->argv) PRINT_ERROR("PC_ENRICH requires Input.\n");

  /* bedfile */
  if(p->n_igregion){
    p->igregion = (BedFile **)my_calloc(p->n_igregion, sizeof(BedFile *), "p->igregion");
    for(i=0; i<p->n_igregion; i++){
      nline = ParseLine_arbit(str_bed[i], clm, ',');
      if(nline >2) PRINT_ERROR("error: -ignore %d has ',' more than 2: %s\n", i+1, str_bed[i]);
      p->igregion[i] = read_bedfile(clm[0].str, g);
      p->igregion[i]->argv = strdup(clm[0].str);
    }
  }

  init_dump(p, *sample);
  scan_mappability_chrtotal(p->mpfile, g);

  MYFREE(str_argv);
  MYFREE(str_bed);
  return;
 err:
  print_usage_dp();
}

static void init_dump(DrParam *p, SamplePair *sample){
  int i;
  char str_wigfiletype[][20]={"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH"};
  char *str_norm[]={ "OFF", "TOTALREAD", "NCIS"};

  printf("\n======== %s ========\n", str_ftype[p->ftype]);
  printf("drompa_peakcall version %s\n", VERSION);
  printf("   output prefix: %s\n", p->headname);
  printf("   genometable: %s\n", p->gtfile);
  printf("   ChIP: %s\n", sample->ChIP->argv);
  if(sample->Input->argv) printf("   Input: %s\n", sample->Input->argv);
  else printf("   Input: NONE\n");
  printf("   Input format: %s\n", str_wigfiletype[p->itype]);
  printf("   binsize: %d\n", sample->binsize);
  if(p->smoothing) printf("   smoothing width: %d bp\n", p->smoothing);
  if(p->mpfile) printf("Mappability file: %s\n", p->mpfile);
  printf("   Peak intensity threshold: %.2f\n", p->IPmaxthre);
  printf("   ChIP/Input normalization: %s\n", str_norm[p->ntype]);
  printf("   Enrichment threshold: %.2f\n", p->enrichthre);
  printf("   p-value threshold (internal, -log10): %.2e\n", p->pthre_internal);
  printf("   p-value threshold (enrichment, -log10): %.2e\n", p->pthre_enrich);
  printf("   q-value threshold: %.2e\n", p->qthre);
  if(p->n_igregion){
    printf("   ignore list:");
    for(i=0; i<p->n_igregion; i++) printf(" %s", p->igregion[i]->argv);
    printf("\n");
  }
  printf("======================================\n");
  return;
}
