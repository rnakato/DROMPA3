/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "pw_init.h"
#include "argument.h"
#include "readfile.h"
#include "filehandle.h"
#include "macro.h"
#include "alloc.h"

static void check_pwparam(PwParam *p, RefGenome *g, char *, char *);
static void init_dump(PwParam *p);

static void print_usage(){
  fprintf(stderr, "parse2wig version %s\n", VERSION);
  fprintf(stderr, "Usage: parse2wig [option] {-pair} -i <inputfile> -o <output> -gt <genome_table>\n\n");
  fprintf(stderr, "       <inputfile>\tMapping file. Multiple files are allowed (separated by ',')\n");
  fprintf(stderr, "       <output>\tPrefix of output files\n");
  fprintf(stderr, "       <genome_table>\tTab-delimited file describing the name and length of each chromosome\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "       -f {SAM|BAM|BOWTIE}: format of input file (default:SAM) \n");
  fprintf(stderr, "       -odir: output directory name (default: 'parse2wigdir') \n");
  fprintf(stderr, "       -binsize: bin size (default: %d bp) \n", BINSIZE_DEFAULT);
  fprintf(stderr, "       -rcenter <int>: consider <int> bp around the center of fragment \n");
  fprintf(stderr, "       -of <int>: output format (default: 0)\n");
  fprintf(stderr, "           0: binary (.bin) \n");
  fprintf(stderr, "           1: compressed wig (.wig.gz) \n");
  fprintf(stderr, "           2: uncompressed wig (.wig) \n");
  fprintf(stderr, "           3: bedGraph (.bedGraph) \n");
  fprintf(stderr, "           4: bigWig (.bw) \n");
  fprintf(stderr, "       -num4cmp <int>: read number for calculating complexity (default: %d M) \n", NUM4CMP_DEFAULT);
  fprintf(stderr, "       -ccp: make cross-correlation profile (default: off)\n");
  fprintf(stderr, "       For single-end file:\n");
  fprintf(stderr, "           -flen: expected fragment length (default: %d bp) \n", FRAGMENT_LEN);
  fprintf(stderr, "            (Automatically calculated for paired-end mode)\n");
  fprintf(stderr, "       For paired-end file:\n");
  fprintf(stderr, "           -pair: add when the input file is paired-end \n");
  fprintf(stderr, "           -maxins: maximum fragment length (default: %d bp) \n", MAX_FRAGMENT_LEN);
  fprintf(stderr, "       To take peaks into consideration:\n");
  fprintf(stderr, "           -bed <bedfile>: specify the BED file of enriched regions (e.g., peak regions)\n\n");
  fprintf(stderr, "       For normalization:\n");
  fprintf(stderr, "           Total read normalization: \n");
  fprintf(stderr, "              -n {NONE|GR|GD|CR|CD} (default:NONE)\n");
  fprintf(stderr, "                 NONE;\tnot normalize\n");
  fprintf(stderr, "                 GR;\tfor whole genome, read number\n");
  fprintf(stderr, "                 GD;\tfor whole genome, read depth\n");
  fprintf(stderr, "                 CR;\tfor each chromosome, read number\n");
  fprintf(stderr, "                 CD;\tfor each chromosome, read depth\n");
  fprintf(stderr, "              -np <int>:\tread number after normalization (default:%d (%d million))\n", NUM4RPM_DEFAULT, NUM4RPM_DEFAULT/NUM_1M);
  fprintf(stderr, "              -nd <double>:\tread depth after normalization (default:%.1f)\n", NUM4DEPTH_DEFAULT);
  fprintf(stderr, "           PCR bias: \n");
  fprintf(stderr, "              -nofilter: do not filter PCR bias \n");
  fprintf(stderr, "              -thre_pb: PCRbias threshold (default: more than max(1 read, 10 times greater than genome average)) \n");
  fprintf(stderr, "           Mappability: \n");
  fprintf(stderr, "              -mp <mappability file>: specify when normalizing reads with mappability \n");
  fprintf(stderr, "              -mpthre <double>: threshold of low mappability (default: < %.1f)\n", THRE_LOW_MAPPABILITY);
  fprintf(stderr, "           GC content: (require genome sequence, take large time and memory) \n");
  fprintf(stderr, "              -GC <genome>: specify the reference genome sequence used for read mapping\n");
  fprintf(stderr, "              -flen4gc <int>: fragment length for calculation of GC distribution\n");
  fprintf(stderr, "              -gcdepthoff: do not consider depth of GC contents\n");
  fprintf(stderr, "\n");
  exit(0);
}

static void print_version(){
  printf("parse2wig version %s\n", VERSION);
  exit(0);
}

void pw_argv_init(int argc, char *argv[], PwParam *p, RefGenome *g){
  if(argc<=1) print_usage();
  
  char *ftype=NULL, *ntype=NULL;
  const Argument args[]={
    {"--version", ARGUMENT_TYPE_FUNCTION, print_version,     NULL},
    {"-h"    ,    ARGUMENT_TYPE_FUNCTION, print_usage,       NULL},
    {"--help",    ARGUMENT_TYPE_FUNCTION, print_usage,       NULL},
    {"-i",        ARGUMENT_TYPE_STRING,   &p->inputfile,     NULL},
    {"-o",        ARGUMENT_TYPE_STRING,   &p->output_prefix, NULL},
    {"-gt",       ARGUMENT_TYPE_STRING,   &p->gtfile,        NULL},
    {"-odir",     ARGUMENT_TYPE_STRING,   &p->output_dir,    NULL},
    {"-pair",     ARGUMENT_TYPE_FLAG_ON,  &p->rtype,         NULL},
    {"-ccp",      ARGUMENT_TYPE_FLAG_ON,  &p->ccp,           NULL},
    {"-opair",    ARGUMENT_TYPE_FLAG_ON,  &p->out_readpair,  NULL},
    {"-binsize",  ARGUMENT_TYPE_INTEGAR,  &p->binsize,       NULL},
    {"-flen",     ARGUMENT_TYPE_INTEGAR,  &p->fraglen,       NULL},
    {"-maxins",   ARGUMENT_TYPE_INTEGAR,  &p->max_fraglen,   NULL},
    {"-nofilter", ARGUMENT_TYPE_FLAG_OFF, &p->pcrfilter,     NULL},
    {"-thre_pb",  ARGUMENT_TYPE_INTEGAR,  &p->thre_filter,   NULL},
    {"-n",        ARGUMENT_TYPE_STRING,   &ntype,            NULL},
    {"-np",       ARGUMENT_TYPE_INTEGAR,  &p->num4rpm,       NULL},
    {"-nd",       ARGUMENT_TYPE_FLOAT,    &p->num4depth,     NULL},
    {"-mp",       ARGUMENT_TYPE_STRING,   &p->mpfile,        NULL},
    {"-mpthre",   ARGUMENT_TYPE_FLOAT,    &p->mpthre,        NULL},
    {"-mpbin",    ARGUMENT_TYPE_STRING,   &p->mpbinaryfile,  NULL},
    {"-f",        ARGUMENT_TYPE_STRING,   &ftype,            NULL},
    {"-of",       ARGUMENT_TYPE_INTEGAR,  &p->wtype,         NULL},
    {"-rcenter",  ARGUMENT_TYPE_INTEGAR,  &p->usereadcenter, NULL},
    {"-bed",      ARGUMENT_TYPE_STRING,   &p->bedfilename,   NULL},
    {"-num4cmp",  ARGUMENT_TYPE_INTEGAR,  &p->num4cmp,       NULL},
    {"-GC",       ARGUMENT_TYPE_STRING,   &p->genomefile,    NULL},
    {"-flen4gc",  ARGUMENT_TYPE_INTEGAR,  &p->flen4gc,       NULL},
    {"-gcdepthoff", ARGUMENT_TYPE_FLAG_OFF,  &p->gcdepth,    NULL},
    {NULL,        ARGUMENT_TYPE_NONE,     NULL,              NULL},
  };
  argument_read(&argc, argv, args);

  check_pwparam(p, g, ftype, ntype);
  init_dump(p);
  scan_mappability_chrtotal(p->mpfile, g);

  return;
}

static void check_pwparam(PwParam *p, RefGenome *g, char *ftype, char *ntype){
  int chr;
  if(!p->inputfile){ fprintf(stderr, "please specify inputfile.\n"); goto err;}
  if(!p->output_dir) p->output_dir = strdup("parse2wigdir");
  if(!p->output_prefix){ fprintf(stderr, "please specify output prefix.\n"); goto err;}
  if(p->binsize <=0){ fprintf(stderr, "error: binsize = %d.\n", p->binsize); goto err;}
  if(!p->gtfile){ fprintf(stderr, "please specify genome_table.\n");         goto err;}
  else parse_genometable(p->gtfile, g);
  p->binnum_chr = (int *)my_calloc(g->chrnum, sizeof(int), "p->binnum_chr");
  for(chr=1; chr<g->chrnum; chr++){
    p->binnum_genome += p->binnum_chr[chr] = g->chr[chr].len/p->binsize +1;
  }
  if(p->rtype==READTYPE_SINGLE && p->fraglen <=0){   fprintf(stderr, "error: fragment length = %d.\n", p->fraglen);             goto err;}
  if(p->rtype==READTYPE_PAIR && p->max_fraglen <=0){ fprintf(stderr, "error: maximum fragment length = %d.\n", p->max_fraglen); goto err;}
  if(!range(p->wtype, 0, PWFILETYPENUM-1)){          fprintf(stderr, "error: Invalid wigfile type = %d.\n", p->wtype);          goto err;}
  if(p->pcrfilter && p->thre_filter <0) {            fprintf(stderr, "error: -thre_pb = %d.\n", p->thre_filter); goto err;}
  if(p->bedfilename){
    isfile(p->bedfilename);
    p->enrichfile = read_bedfile(p->bedfilename, g);
    //    show_bedfile(p->enrichfile, g->chrnum);
  }
  if(p->genomefile){
    if(!p->mpbinaryfile){ fprintf(stderr, "error: -GC option requires -mpbin option.\n"); goto err;}
    isfile(p->genomefile);
  }

  if(ftype){
    if(!strcmp(ftype, "SAM"))         p->ftype = FILETYPE_SAM;
    else if(!strcmp(ftype, "BAM"))    p->ftype = FILETYPE_BAM;
    else if(!strcmp(ftype, "BOWTIE")) p->ftype = FILETYPE_BOWTIE;
    else{ fprintf(stderr, "Invalid input filetype: %s\n", ftype); exit(0);}
  }

  if(ntype){
    if(!strcmp(ntype, "NONE"))    p->ntype = NORMTYPE_NONE;
    else if(!strcmp(ntype, "GR")) p->ntype = NORMTYPE_GENOME_READ;
    else if(!strcmp(ntype, "GD")) p->ntype = NORMTYPE_GENOME_DEPTH;
    else if(!strcmp(ntype, "CR")) p->ntype = NORMTYPE_CHROM_READ;
    else if(!strcmp(ntype, "CD")) p->ntype = NORMTYPE_CHROM_DEPTH;
    else{ fprintf(stderr, "Invalid normalize type: %s\n", ntype); exit(0);}
  }
  if(p->ntype==NORMTYPE_GENOME_READ || p->ntype==NORMTYPE_CHROM_READ){
    if(p->num4rpm <=0){ printf("error: -np = %d.\n", p->num4rpm); goto err;}
  }
  if(p->ntype==NORMTYPE_GENOME_DEPTH || p->ntype==NORMTYPE_CHROM_DEPTH){
    if(p->num4depth <=0){ printf("error: -nd = %.2f.\n", p->num4depth); goto err;}
  }

  if(p->mpfile && p->mpthre <=0){
    printf("error: threshold for mappability = %.2f.\n", p->mpthre);
    goto err;
  }

  mkdir(p->output_dir, 0755);

  return;
 err:
  print_usage();
}

static void init_dump(PwParam *p){
  char str_Inputfiletype[][10]={"SAM", "BAM", "BOWTIE"};
  char str_end[][16]={"SINGLE_END", "PAIRED_END"};
  char str_bool[][10]={"OFF", "ON"};
  char str_wigfiletype[][20]={"BINARY", "COMPRESSED WIG", "WIG", "BEDGRAPH", "BIGWIG"};
  char str_normtype[][20]={"NONE", "GENOME_READ", "GENOME_DEPTH", "CHROM_READ", "CHROM_DEPTH"};
  
  printf("\n======================================\n");
  printf("parse2wig version %s\n\n", VERSION);
  printf("Input file: %s\n", p->inputfile);
  printf("\tFormat: %s, %s\n", str_Inputfiletype[p->ftype], str_end[p->rtype]);
  printf("Output file: %s/%s\n", p->output_dir, p->output_prefix);
  printf("\tFormat: %s\n", str_wigfiletype[p->wtype]);
  printf("Genome-table file: %s\n", p->gtfile);
  printf("Binsize: %d bp\n", p->binsize);
  if(p->rtype==READTYPE_SINGLE) printf("Fragment length: %d\n", p->fraglen);
  else printf("Maximum fragment length: %d\n", p->max_fraglen);
  printf("PCR bias filtering: %s\n", str_bool[p->pcrfilter]);
  printf("Read number for library complexity: %.1f M\n", p->num4cmp/(double)NUM_1M);
  if(p->pcrfilter && p->thre_filter) printf("PCR bias threshold: >%d\n", p->thre_filter);
  if(p->bedfilename) printf("Enriched regions file: %s (%d sites)\n", p->bedfilename, p->enrichfile->num);
  printf("Cross-correlation profile: %s\n", str_bool[p->ccp]);
  
  printf("\nTotal read normalization: %s", str_normtype[p->ntype]);
  if(p->ntype==NORMTYPE_GENOME_READ || p->ntype==NORMTYPE_CHROM_READ){
    printf("\tnormed read: %.1f M for genome\n", p->num4rpm/(double)NUM_1M);
  }
  if(p->ntype==NORMTYPE_GENOME_DEPTH || p->ntype==NORMTYPE_CHROM_DEPTH){
    printf("\tnormed depth: %.2f.\n", p->num4depth);
  }
  else printf("\n");
  if(p->mpfile){
    printf("Mappability normalization:\n");
    printf("\tfile prefix: %s\n", p->mpfile);
    printf("\tLow mappablitiy threshold: %.2f\n", p->mpthre);
  }
  if(p->genomefile){
    printf("Correcting GC bias:\n");
    printf("\tChromosome directory: %s\n", p->genomefile);
    printf("\tmappability binary prefix: %s\n", p->mpbinaryfile);
    printf("\tLength for GC distribution: %d\n", p->flen4gc);
    printf("\tConsidering gc depth: %s\n", str_bool[p->gcdepth]);
  }
  printf("======================================\n");
  return;
}
