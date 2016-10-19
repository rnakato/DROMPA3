/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dd_init.h"
#include "dd_gv.h"
#include "drompa_readfile.h"
#include "readfile.h"
#include "argument.h"
#include "common.h"
#include "macro.h"
#include "alloc.h"
#include "stringp.h"
#include "filehandle.h"
#include "drompa_usage.h"

typedef struct{
  char **str_argv, **str_argv_overlay;
  char **str_pd;
  char **str_bed;
  char **str_inter;
} StructInit;

typedef struct{
  int binsize, binsize_overlay;
  double scale_tag, scale_tag_overlay;
  double scale_ratio, scale_ratio_overlay;
  double scale_pvalue, scale_pvalue_overlay;
} SampleParam;

#define PRINT_ERROR(args...) do{ fprintf(stderr, args); goto err; }while(0)

static void check_argv1(DrParam *p, DDParam *d, char *argv, SampleParam *sp);
static void check_ddparam(DrParam *p, DDParam *d, SamplePair **sample, StructInit *st, RefGenome *g, SampleParam *sp);
static void init_dump(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample);
static void check_sample_copy(DrParam *p, SamplePair **sample);
static void print_version();
static void print_usage_base();

static StructInit *StructInit_new(int num){
  if(num<=0){ printf("error:%s: num = %d.\n",__FUNCTION__, num); exit(0);}
  StructInit *p = (StructInit *)my_calloc(num, sizeof(StructInit), "StructInit");
  p->str_argv   = (char **)my_calloc(num, sizeof(char *), "str_argv");
  p->str_argv_overlay = (char **)my_calloc(num, sizeof(char *), "str_argv_overlay");
  p->str_pd     = (char **)my_calloc(num, sizeof(char *), "str_pd");
  p->str_bed    = (char **)my_calloc(num, sizeof(char *), "str_bed");
  p->str_inter  = (char **)my_calloc(num, sizeof(char *), "str_inter");
  return p;
}

static void StructInit_delete(StructInit *p){
  MYFREE(p->str_argv);
  MYFREE(p->str_argv_overlay);
  MYFREE(p->str_pd);
  MYFREE(p->str_bed);
  MYFREE(p->str_inter);
  MYFREE(p);
  return;
}

void dd_argv_init(int argc, char **argv, DrParam *p, DDParam *d, SamplePair **sample, RefGenome *g){
  if(argc<=1) print_usage_base();
  
  SampleParam sp;
  sp.binsize      = sp.binsize_overlay      = BINSIZE_DEFAULT;
  sp.scale_tag    = sp.scale_tag_overlay    = SCALE_TAG_DEFAULT;
  sp.scale_ratio  = sp.scale_ratio_overlay  = SCALE_RATIO_DEFAULT;
  sp.scale_pvalue = sp.scale_pvalue_overlay = SCALE_PVALUE_DEFAULT;

  check_argv1(p, d, argv[1], &sp);
  argv++; argc--;
  if(argc<=1) print_error_dd(p->ftype);

  StructInit *st = StructInit_new(1024);
  const Argument args[]={
    {"--version", ARGUMENT_TYPE_FUNCTION, print_version,        NULL},
    {"-h"    ,    ARGUMENT_TYPE_FUNCTION, print_usage_base,     NULL},
    {"--help",    ARGUMENT_TYPE_FUNCTION, print_usage_base,     NULL},
    {"-gt",       ARGUMENT_TYPE_STRING,   &p->gtfile,           NULL},
    {"-p",        ARGUMENT_TYPE_STRING,   &p->headname,         NULL},
    {"-i",        ARGUMENT_TYPE_STRING_MULTI, &(st->str_argv), &(p->samplenum_1st)},
    {"-ioverlay", ARGUMENT_TYPE_STRING_MULTI, &(st->str_argv_overlay), &(p->samplenum_overlay)},
    {"-pd",       ARGUMENT_TYPE_STRING_MULTI, &(st->str_pd),   &(d->pdnum)},
    {"-binsize",  ARGUMENT_TYPE_INTEGAR,  &sp.binsize,          NULL},
    {"-binsize2", ARGUMENT_TYPE_INTEGAR,  &sp.binsize_overlay,  NULL},
    {"-mp",       ARGUMENT_TYPE_STRING,   &p->mpfile,           NULL},
    {"-mpthre",   ARGUMENT_TYPE_FLOAT,    &p->mpthre,           NULL},
    {"-gap",      ARGUMENT_TYPE_STRING,   &p->gapfile,          NULL},
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
    /* DDParam */
    {"-nosig",    ARGUMENT_TYPE_FLAG_OFF, &d->do_peakcall,      NULL},
    {"-chr",      ARGUMENT_TYPE_INTEGAR,  &d->chronly,          NULL},
    {"-bed",      ARGUMENT_TYPE_STRING_MULTI, &(st->str_bed), &(d->bednum)},
    {"-inter",    ARGUMENT_TYPE_STRING_MULTI, &(st->str_inter), &(d->internum)},
    {"-GC",       ARGUMENT_TYPE_STRING,   &d->GC.argv,          NULL},
    {"-gcsize",   ARGUMENT_TYPE_INTEGAR,  &d->GC.wsize,         NULL},
    {"-GD",       ARGUMENT_TYPE_STRING,   &d->GD.argv,          NULL},
    {"-gdsize",   ARGUMENT_TYPE_INTEGAR,  &d->GD.wsize,         NULL},
    {"-gene",     ARGUMENT_TYPE_STRING,   &d->gene.argv,        NULL},
    {"-gftype",   ARGUMENT_TYPE_INTEGAR,  &d->gftype,           NULL},
    {"-ars",      ARGUMENT_TYPE_STRING,   &d->arsfile,          NULL},
    {"-ter",      ARGUMENT_TYPE_STRING,   &d->terfile,          NULL},
    {"-showars",  ARGUMENT_TYPE_FLAG_ON,  &d->show_ars,         NULL},
    {"-repeat",   ARGUMENT_TYPE_STRING,   &d->repeat.argv,      NULL},
    {"-pdsize",   ARGUMENT_TYPE_INTEGAR,  &d->pdsize,           NULL},
    {"-png",      ARGUMENT_TYPE_FLAG_ON,  &d->png,              NULL},
    {"-rmchr",    ARGUMENT_TYPE_FLAG_ON,  &d->rmchr,            NULL},
    {"-ls",       ARGUMENT_TYPE_INTEGAR,  &d->width_per_line,   NULL},
    {"-lpp",      ARGUMENT_TYPE_INTEGAR,  &d->linenum_per_page, NULL},
    {"-r",        ARGUMENT_TYPE_STRING,   &d->drawregion_argv,  NULL},
    {"-genefile", ARGUMENT_TYPE_STRING,   &d->genefile_argv,  NULL},
    {"-len_genefile", ARGUMENT_TYPE_INTEGAR, &d->genefile_len,NULL},
    {"-scale_tag",    ARGUMENT_TYPE_FLOAT, &sp.scale_tag,     NULL},
    {"-scale_ratio",  ARGUMENT_TYPE_FLOAT, &sp.scale_ratio,   NULL},
    {"-scale_pvalue", ARGUMENT_TYPE_FLOAT, &sp.scale_pvalue,  NULL},
    {"-scale_tag2",   ARGUMENT_TYPE_FLOAT, &sp.scale_tag_overlay,   NULL},
    {"-scale_ratio2", ARGUMENT_TYPE_FLOAT, &sp.scale_ratio_overlay, NULL},
    {"-scale_pvalue2",ARGUMENT_TYPE_FLOAT, &sp.scale_pvalue_overlay,NULL},
    {"-show_ctag", ARGUMENT_TYPE_INTEGAR, &d->visualize_ctag,   NULL},
    {"-show_itag", ARGUMENT_TYPE_INTEGAR, &d->visualize_itag,   NULL},
    {"-showratio", ARGUMENT_TYPE_INTEGAR, &d->visualize_ratio,  NULL},
    {"-showpinter",  ARGUMENT_TYPE_INTEGAR, &d->visualize_p_inter,  NULL},
    {"-showpenrich", ARGUMENT_TYPE_INTEGAR, &d->visualize_p_enrich, NULL},
    {"-bn",     ARGUMENT_TYPE_INTEGAR, &d->barnum,      NULL},
    {"-ystep",  ARGUMENT_TYPE_FLOAT,   &d->ystep,       NULL},
    {"-offbg",  ARGUMENT_TYPE_FLAG_OFF, &d->backcolors, NULL},
    {"-offymem",ARGUMENT_TYPE_FLAG_OFF, &d->stroke_ymem,NULL},
    {"-offylab",ARGUMENT_TYPE_FLAG_OFF, &d->stroke_ylab,NULL},
    {"-viz",    ARGUMENT_TYPE_INTEGAR, &d->viz,        NULL},
    /* CG */
    {"-cgthre", ARGUMENT_TYPE_FLOAT, &d->cgthre,        NULL},
    /* PD */
    {"-prop",   ARGUMENT_TYPE_FLAG_ON,  &d->pd_prop,    NULL},
    /* TR */
    {"-tssthre", ARGUMENT_TYPE_FLOAT, &d->tssthre,      NULL},
    /* profile and heatmap */
    {"-ptype", ARGUMENT_TYPE_INTEGAR, &d->ptype,     NULL},
    {"-stype", ARGUMENT_TYPE_INTEGAR, &d->stype,     NULL},
    {"-ntype", ARGUMENT_TYPE_INTEGAR, &d->ntype,     NULL},
    {"-cw",    ARGUMENT_TYPE_INTEGAR, &d->compwidth, NULL},
    {"-offse", ARGUMENT_TYPE_FLAG_OFF, &d->showse,   NULL},
    {"-hmsort", ARGUMENT_TYPE_INTEGAR, &d->hmsort,   NULL},
    {"-sortgbody", ARGUMENT_TYPE_FLAG_ON, &d->sortgbody, NULL},
    {"-pdetail",   ARGUMENT_TYPE_FLAG_ON, &d->pdetail,   NULL},
    {NULL, ARGUMENT_TYPE_NONE, NULL, NULL},
  };
  argument_read(&argc, argv, args);
  check_ddparam(p, d, sample, st, g, &sp);

  init_dump(p, d, g, *sample);

  StructInit_delete(st);
  return;
}

static void check_argv1(DrParam *p, DDParam *d, char *argv, SampleParam *sp){
  if(!strcmp(argv, "-h") || !strcmp(argv, "--help")) print_usage_base();
  if(!strcmp(argv, "--version")) print_version();
  else if(!strcmp(argv, "PC_SHARP")){
    p->ftype = FTYPE_PEAKCALL_SHARP;
    d->do_peakcall = 1;
    d->makefig = 1;
  }else if(!strcmp(argv, "PC_BROAD")){
    p->ftype = FTYPE_PEAKCALL_BROAD;
    d->do_peakcall = 1;
    d->makefig = 1;
    p->smoothing = SMOOTHING_BROAD;
    d->width_per_line = LS_DEFAULT*2;
    sp->binsize     = sp->binsize_overlay     = BINSIZE_BROAD;
    sp->scale_tag   = sp->scale_tag_overlay   = SCALE_TAG_BROAD;
    sp->scale_ratio = sp->scale_ratio_overlay = SCALE_RATIO_BROAD;
  }else if(!strcmp(argv, "PC_ENRICH")){
    p->ftype = FTYPE_PEAKCALL_E;
    d->do_peakcall = 1;
    p->enrichthre = ENRICHTHRE_DEFAULT;
    d->visualize_ctag = 0;
    d->visualize_ratio = 1;
    d->makefig = 1;
  }else if(!strcmp(argv, "GV")){
    p->ftype = FTYPE_GV;
    d->gw = true;
    d->makefig = 1;
    d->visualize_ctag = 0;
    d->visualize_ratio = 1;
    d->rmchr = 1;
    sp->binsize = sp->binsize_overlay = BINSIZE_DEFAULT_GV;
    sp->scale_ratio = sp->scale_ratio_overlay = 1;
  }else if(!strcmp(argv, "PD")){
    p->ftype = FTYPE_PD;
    d->gw = true;
    d->makefig = 1;
    d->rmchr = 1;
  }else if(!strcmp(argv, "FRIP"))     p->ftype = FTYPE_FRIP;
  else if(!strcmp(argv, "CI"))        p->ftype = FTYPE_COMPARE_INTENSITY;
  else if(!strcmp(argv, "CG"))        p->ftype = FTYPE_COMPARE_GENEBODY;
  else if(!strcmp(argv, "GOVERLOOK")) p->ftype = FTYPE_GOVERLOOK;
  else if(!strcmp(argv, "PROFILE"))   p->ftype = FTYPE_PROFILE;
  else if(!strcmp(argv, "HEATMAP"))   p->ftype = FTYPE_HEATMAP;
  else if(!strcmp(argv, "TR"))        p->ftype = FTYPE_TR;
  else{
    fprintf(stderr, "Please specify command.\n\n");
    print_usage_base();
  }
  return;
}

static Graph *scan_pdstr(DrParam *p, char **str, int pdnum, int pdsize){
  int i, nline;
  Elem clm[ELEM_NUM];
  Graph *graph = (Graph *)my_calloc(pdnum, sizeof(Graph), "d->PD");
  for(i=0; i<pdnum; i++){
    nline = ParseLine_arbit(str[i], clm, ',');
    if(nline >2){
      fprintf(stderr, "error: sample %d has ',' more than 2: %s\n", i+1, str[i]);
      print_error_dd(p->ftype);
    }
    LOG("%d PD:%s name:%s\n", i+1, clm[0].str, clm[1].str);
    if(!strcmp(clm[0].str, "")){
      fprintf(stderr, "please specify ChIP sample %d: %s.\n", i+1, str[i]);
      print_error_dd(p->ftype);
    }else                      graph[i].argv = strdup(clm[0].str);
    if(strcmp(clm[1].str, "")) graph[i].name = strdup(clm[1].str);
    else                       graph[i].name = graph[i].argv;
    graph[i].wsize = pdsize;
  }
  return graph;
}

static void check_ddparam(DrParam *p, DDParam *d, SamplePair **sample, StructInit *st, RefGenome *g, SampleParam *sp){
  int i,j,chr, nline, lenmax=0;
  Elem clm[ELEM_NUM];

  if(print_error_peakcall(p,g)) goto err;
  if(g->genome->len > THRE_GENOMELEN4WG) d->large_genome = true;

  if((p->ftype==FTYPE_COMPARE_GENEBODY || p->ftype == FTYPE_TR) && !d->gene.argv) PRINT_ERROR("error: %s requires -gene option.\n", str_ftype[p->ftype]);
  if((p->ftype==FTYPE_FRIP || p->ftype==FTYPE_COMPARE_INTENSITY) && !d->bednum)   PRINT_ERROR("error: %s requires -bed option.\n", str_ftype[p->ftype]);
  if(p->ftype==FTYPE_GOVERLOOK && !range(d->bednum, 1, 3)) PRINT_ERROR("error: please specify -bed option (up to 3.)\n");
  if(p->ftype==FTYPE_PROFILE || p->ftype==FTYPE_HEATMAP){
    if(!range(d->ptype, 1, PTYPENUM-1))        PRINT_ERROR("error: Invalid input ptype: %d.\n",       d->ptype);
    if(d->ptype==BEDSITES && !d->bednum)       PRINT_ERROR("error: ptype %d requires -bed option.\n", d->ptype);
    if(d->ptype!=BEDSITES && !d->gene.argv)    PRINT_ERROR("error: ptype %d requires -gene option.\n",d->ptype);
    if(!range(d->stype, READDIST, PVALUEDIST)) PRINT_ERROR("error: Invalid input stype: %d.\n",       d->stype);
    if(!range(d->ntype, 0, 1))                 PRINT_ERROR("error: Invalid input ntype: %d.\n",       d->ntype);
  }
  if(!range(d->visualize_ratio, 0, 2))         PRINT_ERROR("error: Invalid input showratio: %d.\n",   d->visualize_ratio);


  if(p->ftype == FTYPE_GV || p->ftype == FTYPE_PD){
    d->drawregion_argv=NULL;
    d->genefile_argv=NULL;
  }
  if(d->makefig){
    if(d->genefile_argv){
      d->drawregion_argv=NULL;
      isfile(d->genefile_argv);
      d->drawregion = bedfile_new(g->chrnum);
      //      show_bedfile(d->drawregion, g->chrnum);
    }else if(d->drawregion_argv){
      isfile(d->drawregion_argv);
      d->drawregion = read_bedfile(d->drawregion_argv, g);
      //      show_bedfile(d->drawregion, g->chrnum);
    }
    if(d->width_per_line<=0) PRINT_ERROR("error: -ls= %d.\n", d->width_per_line);
    else d->width_per_line *= NUM_1K;
    if(d->gw==true){
      for(i=1; i<g->chrnum; i++){
	if(lenmax < g->chr[i].len) lenmax = g->chr[i].len;
      }
      d->width_per_line = lenmax +1;
      d->gene.argv = NULL;
      d->repeat.argv = NULL;
      d->arsfile = NULL;
    }
    if(d->internum){ // interaction files
      d->inter = (InteractionSet *)my_calloc(d->internum, sizeof(InteractionSet), "d->inter");
      for(i=0; i<d->internum; i++){
	nline = ParseLine_arbit(st->str_inter[i], clm, ',');
	if(nline >2) PRINT_ERROR("error: -inter %d has ',' more than 2: %s\n", i+1, st->str_inter[i]);
	isfile(clm[0].str);
	read_interactionfile(&(d->inter[i]), clm[0].str, g);
	if(nline==2) d->inter[i].name = strdup(clm[1].str);
	else d->inter[i].name = "Interaction";
      }
    }
    d->command_mergepdf = (char *)my_calloc((strlen(p->headname)+100)*g->chrnum, sizeof(char), "command_margepdf");
    sprintf(d->command_mergepdf, "cpdf ");
  }

  if(d->bednum){ // bedfiles
    if(p->ftype==FTYPE_PROFILE || p->ftype==FTYPE_HEATMAP) d->peak = (Peak **)my_calloc(d->bednum, sizeof(Peak *), "d->peak");
    else d->bed = (BedFile **)my_calloc(d->bednum, sizeof(BedFile *), "d->bed");
    for(i=0; i<d->bednum; i++){
      nline = ParseLine_arbit(st->str_bed[i], clm, ',');
      if(nline >2) PRINT_ERROR("error: bed %d has ',' more than 2: %s\n", i+1, st->str_bed[i]);
      isfile(clm[0].str);
      if(p->ftype==FTYPE_PROFILE || p->ftype==FTYPE_HEATMAP){
	d->peak[i] = read_peakfile(clm[0].str, g);
	d->peak[i]->argv = strdup(clm[0].str);
	if(nline==2) d->peak[i]->name = strdup(clm[1].str);
	else d->peak[i]->name = d->peak[i]->argv;
      }else{
	d->bed[i] = read_bedfile(clm[0].str, g);
	d->bed[i]->argv = strdup(clm[0].str);
	if(nline==2) d->bed[i]->name = strdup(clm[1].str);
	else d->bed[i]->name = d->bed[i]->argv;
	//      show_bedfile(d->bed[i], g->chrnum);
      }
    }
  }
  if(p->ftype==FTYPE_PD){
    if(!d->pdnum) PRINT_ERROR("error: please specify sample file.\n");
    else d->PD = scan_pdstr(p, st->str_pd, d->pdnum, d->pdsize);
  }else if(p->ftype != FTYPE_GOVERLOOK){
    if(!p->samplenum_1st) PRINT_ERROR("error: please specify sample file.\n");
    else{
      p->samplenum = p->samplenum_1st + p->samplenum_overlay;
      (*sample) = scan_samplestr(p, st->str_argv, st->str_argv_overlay, g->chrnum);
    }
    check_sample_copy(p, sample);
    for(i=0; i<p->samplenum; i++){
      if((*sample)[i].peak_argv){
	isfile((*sample)[i].peak_argv);
	(*sample)[i].peak = read_peakfile((*sample)[i].peak_argv, g);
      }
      if((*sample)[i].binsize <= 0){
	if(i<p->samplenum_1st) (*sample)[i].binsize = sp->binsize;
	else (*sample)[i].binsize = sp->binsize_overlay;
      }
      for(chr=1; chr<g->chrnum; chr++){
	(*sample)[i].binnum[chr] = g->chr[chr].len/(*sample)[i].binsize +1;
	if(g->chr[chr].len < (*sample)[i].binsize) fprintf(stderr, "Warning: length of %s (%ld) is shorter than binsize (%d).\n", g->chr[chr].name, g->chr[chr].len, (*sample)[i].binsize);
      }
      if((*sample)[i].scale_tag    <= 0){
	if(i<p->samplenum_1st) (*sample)[i].scale_tag = sp->scale_tag;
	else (*sample)[i].scale_tag = sp->scale_tag_overlay;
      }
      if((*sample)[i].scale_ratio  <= 0){
	if(i<p->samplenum_1st) (*sample)[i].scale_ratio = sp->scale_ratio;
	else (*sample)[i].scale_ratio = sp->scale_ratio_overlay;
      }
      if((*sample)[i].scale_pvalue <= 0){
	if(i<p->samplenum_1st) (*sample)[i].scale_pvalue = sp->scale_pvalue;
	else (*sample)[i].scale_pvalue = sp->scale_pvalue_overlay;
      }
    }
  }
  if(p->ftype==FTYPE_PROFILE || p->ftype==FTYPE_HEATMAP){
    d->cwbin = d->compwidth/(*sample)[0].binsize;
    int num;
    if(d->ptype==GENE100) num = 300; //GENEBLOCKNUM*3;
    else num = d->cwbin*2+1;
    for(i=0; i<p->samplenum; i++){
      (*sample)[i].profile.IP = (double *)my_calloc(num, sizeof(double), "profile.IP");
      if(d->stype == ENRICHDIST) (*sample)[i].profile.Input = (double *)my_calloc(num, sizeof(double), "profile.Input");
      if(d->showse){
	(*sample)[i].profile.SEarray = (double *)my_calloc(num, sizeof(double), "profile.SE");
	(*sample)[i].profile.SE      = (double **)my_calloc(num, sizeof(double *), "profile.SE");
	d->SEnum = PEAKNUM_DEFAULT;
	for(j=0; j<num; j++) (*sample)[i].profile.SE[j] = (double *)my_calloc(d->SEnum, sizeof(double), "profile.SE[j]");
      }
    }
  }
  
  return;
 err:
  print_error_dd(p->ftype);
}

static void init_dump(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample){
  int i;
  char str_bool[][20] = {"OFF", "ON"};
  char str_wigfiletype[][20] = {"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH"};
  char str_gftype[][20] = {"refFlat", "Ensembl", "gtf", "SGD"};
  char str_Inputread[][20] = {"OFF", "ALL", "1st"};
  char str_ratio[][20] = {"OFF", "liner", "log"};
  char str_norm[][20] = {"OFF", "TOTALREAD", "NCIS"};
  char str_stype[][60]={ "ChIP read", "Enrichment ratio", "Enrichment P-value"};
  char str_ptype[][20]={ "NONE", "TSS", "TTS", "GENE100", "SPECIFIEDSITES"};
  char str_ntype[][20]={ "WHOLE GENOME", "TARGET REGIONS ONLY"};

  printf("\n======== %s ========\n", str_ftype[p->ftype]);
  printf("drompa_draw version %s\n", VERSION);
  printf("   output prefix: %s\n", p->headname);
  printf("   genometable: %s\n", p->gtfile);
  printf("   Input format: %s\n", str_wigfiletype[p->itype]);
  if(p->smoothing) printf("   smoothing width: %d bp\n", p->smoothing);
  if(p->mpfile) printf("Mappability file: %s\n", p->mpfile);
  printf("   ChIP/Input normalization: %s\n", str_norm[p->ntype]);
  if(p->ftype==FTYPE_PEAKCALL_SHARP || p->ftype==FTYPE_PEAKCALL_BROAD || p->ftype==FTYPE_PEAKCALL_E){
    printf("   Peak intensity threshold: %.2f\n", p->IPmaxthre);
    printf("   Enrichment threshold: %.2f\n", p->enrichthre);
    printf("   p-value threshold (internal, -log10): %.2e\n", p->pthre_internal);
    printf("   p-value threshold (enrichment, -log10): %.2e\n", p->pthre_enrich);
    printf("   q-value threshold: %.2e\n", p->qthre);
  }
  if(p->ftype==FTYPE_PROFILE || p->ftype==FTYPE_HEATMAP){
    printf("   show type: %s\n", str_stype[d->stype]);
    printf("   profile type: %s\n", str_ptype[d->ptype]);
    printf("   profile normalization: %s\n", str_ntype[d->ntype]);
  }
  if(p->ftype != FTYPE_GOVERLOOK){
    printf("\nSamples\n");
    for(i=0; i<p->samplenum_1st; i++){
      printf("   ChIP%d: %s\tname: %s\tbinsize:%d\n", i+1, sample[i].ChIP->argv, sample[i].linename, sample[i].binsize);
      if(sample[i].Input->argv) printf("   Input%d: %s\n", i+1, sample[i].Input->argv);
      if(sample[i].peak_argv) printf("      peak list: %s\n", sample[i].peak_argv);
    }
    for(; i<p->samplenum; i++){
      printf("   Overlayed ChIP%d: %s\tname: %s\tbinsize:%d\n", i+1, sample[i].ChIP->argv, sample[i].linename, sample[i].binsize);
      if(sample[i].Input->argv) printf("   Input%d: %s\n", i-p->samplenum_1st+1, sample[i].Input->argv);
      if(sample[i].peak_argv) printf("      peak list: %s\n", sample[i].peak_argv);
    }
  }
  if(p->ftype==FTYPE_PD){
    for(i=0; i<d->pdnum; i++) printf("   IP%d: %s\tname: %s\n", i+1, d->PD[i].argv, d->PD[i].name);
    printf("   pdsize: %d\n", d->pdsize);
  }
  printf("\nAnnotations:\n");
  if(d->gene.argv)   printf("   gene file: %s, format: %s\n", d->gene.argv, str_gftype[d->gftype]);
  if(d->arsfile)     printf("   ARS file: %s\n",  d->arsfile);
  if(d->terfile)     printf("   TER file: %s\n",  d->terfile);
  if(d->repeat.argv) printf("   repeat file: %s\n", d->repeat.argv);
  if(d->GC.argv)     printf("   GCcontents file: %s\n", d->GC.argv);
  if(d->GD.argv)     printf("   gene density file: %s\n", d->GD.argv);
  if(d->internum){
    for(i=0; i<d->internum; i++) printf("   interaction file %d: %s\n",i+1, d->inter[i].argv);
  }
  if(d->drawregion_argv) printf("   region file: %s\n", d->drawregion_argv);
  for(i=0; i<d->bednum; i++){
    if(p->ftype==FTYPE_PROFILE || p->ftype==FTYPE_HEATMAP) printf("   bedfile%d: %s\n", i+1, d->peak[i]->argv);
    else printf("   bedfile%d: %s, name: %s\n", i+1, d->bed[i]->argv, d->bed[i]->name);
  }
  if(d->chronly) printf("output %s only.\n", g->chr[d->chronly].name);

  if(d->makefig){
    printf("\nFigure parameter:\n");
    printf("   display read: ChIP %s, Input %s\n", str_bool[d->visualize_ctag], str_Inputread[d->visualize_itag]);
    printf("   display enrichment: %s\n", str_ratio[d->visualize_ratio]);
    printf("   display pvalue (internal): %s\n", str_bool[d->visualize_p_inter]);
    printf("   display pvalue (ChIP/Input): %s\n", str_bool[d->visualize_p_enrich]);
    printf("   background color: %s\n", str_bool[d->backcolors]);
    printf("   Y label: %s\n", str_bool[d->stroke_ylab]);
    printf("   Y memory: %s\n", str_bool[d->stroke_ymem]);
  }

  printf("======================================\n");
  return;
}

static void check_sample_copy(DrParam *p, SamplePair **sample){
  int i,j;
  char *argC, *argI;
  for(i=0; i<p->samplenum; i++){
    argC = (*sample)[i].ChIP->argv;
    argI = (*sample)[i].Input->argv;
    for(j=0; j<i; j++){
      if(!strcmp(argC, (*sample)[j].ChIP->argv)){
	(*sample)[i].ChIP = (*sample)[j].ChIP;
	(*sample)[i].copyC = j;
	break;
      }
    }
    if(argI){
      for(j=0; j<i; j++){
	if((*sample)[j].Input->argv && !strcmp(argI, (*sample)[j].Input->argv)){
	  (*sample)[i].Input = (*sample)[j].Input;
	  (*sample)[i].copyI = j;
	  break;
	}
      }
      for(j=0; j<i; j++){
	if(!strcmp(argC, (*sample)[j].ChIP->argv) && !strcmp(argI, (*sample)[j].Input->argv)){
	  (*sample)[i].comp = (*sample)[j].comp;
	  (*sample)[i].copycomp = j;
	  break;
	}
      }
    }
  }
#ifdef DEBUG
  for(i=0; i<p->samplenum; i++){
    LOG("sample %d ChIP:%08x Input:%08x comp:%08x\n", i, (int)(*sample)[i].ChIP, (int)(*sample)[i].Input, (int)(*sample)[i].comp);
  }
#endif
  return;
}

static void print_usage_base(){
  fprintf(stderr, "Usage: drompa_draw [--version] <command>\n");
  fprintf(stderr, "Command: PC_SHARP    peak-calling (for sharp mode)\n");
  fprintf(stderr, "         PC_BROAD    peak-calling (for broad mode)\n");
  fprintf(stderr, "         PC_ENRICH   peak-calling (enrichment ratio)\n");
  fprintf(stderr, "         GV          global-view visualization\n");
  fprintf(stderr, "         PD          peak density\n");
  fprintf(stderr, "         FRIP        accumulate read counts in bed regions specified\n");
  fprintf(stderr, "         CI          compare peak-intensity between two samples\n");
  fprintf(stderr, "         CG          output ChIP-reads in each gene body\n");
  fprintf(stderr, "         GOVERLOOK   genome-wide overlook of peak positions\n");
  fprintf(stderr, "         PROFILE     make R script of averaged read density\n");
  fprintf(stderr, "         HEATMAP     make heatmap of multiple samples\n");
  fprintf(stderr, "         TR          calculate the travelling ratio (pausing index) for each gene\n\n");
  exit(0);
}

static void print_version(){
  printf("drompa_draw version %s\n", VERSION);
  exit(0);
}
