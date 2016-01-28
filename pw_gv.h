/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_GV_H_
#define _PW_GV_H_

#include <stdbool.h>
#include "common.h"
#include "seq.h"

/* default parameter */
#define FRAGMENT_LEN 150
#define MAX_FRAGMENT_LEN 500
#define NUM4RPM_DEFAULT 20000000
#define NUM4DEPTH_DEFAULT 1.0
#define NUM4CMP_DEFAULT 10
#define FLEN4GC_DEFAULT 120

#define DIST_READLEN_MAX 200
#define DIST_FRAGLEN_MAX 1000

#define NUM_DARRAY 100
#define READARRAY_NUM 50000 

/* weight */
#define WEIGHT2INT(v) ((v) * 10000.0)
#define INT2WEIGHT(v) ((v) * (1.0/10000.0))

typedef enum {FILETYPE_SAM, FILETYPE_BAM, FILETYPE_BOWTIE, FILETYPE_TAGALIGN} Inputfiletype;
typedef enum {READTYPE_SINGLE, READTYPE_PAIR} Readtype;
typedef enum{
  NORMTYPE_NONE,
  NORMTYPE_GENOME_READ,
  NORMTYPE_GENOME_DEPTH,
  NORMTYPE_CHROM_READ,
  NORMTYPE_CHROM_DEPTH,
} Normtype;

typedef struct{
  char *inputfile;
  char *output_prefix;
  char *output_dir;
  char *gtfile;
  char *mpfile;

  Inputfiletype ftype;
  Readtype rtype;
  PWfile_Type wtype; 
  Normtype ntype;

  int binsize;
  int *binnum_chr, binnum_genome;
  int fraglen;
  int max_fraglen;
  
  int num4rpm;
  double num4depth;

  int num4cmp;
  double r4cmp;

  bool pcrfilter;
  int thre_filter;
  double mpthre;
  bool usereadcenter;

  /* bed file */
  char *bedfilename;
  BedFile *enrichfile;

  /* for GC normalization */
  char *genomefile;
  char *mpbinaryfile;
  int flen4gc;
  int gcdepth;

  /* cross-correlation*/
  bool ccp;
  bool out_readpair;

} PwParam;

typedef struct{
  int *F3;
  int *F5;
  int *weight;
  bool *delete; // for filtering redundant reads
  bool *ignore; // for ignoring peak regions
  //  char *lcmp;   // for calculating complexity

  int narray;
} Readarray;

typedef struct{
  int dist_readlen_F3[DIST_READLEN_MAX];
  int dist_readlen_F5[DIST_READLEN_MAX];
  int dist_fraglen[DIST_FRAGLEN_MAX +1]; /* +1: over DIST_FRAGLEN_MAX */
} FragStats;

typedef struct{
  TYPE_WIGARRAY max;
  int *darray_all, *darray_bg;
  int num;
  double ave, var;
  double nb_p, nb_n, nb_p0;  /* for negative binomial model */
} WigStatsMember;

typedef struct{
  int n_darray;
  int thre;
  TYPE_WIGARRAY num95;
  WigStatsMember *genome, *chr;
} WigStats;

struct seq{
  long n_read_infile;     /* allow multiread (-kn) */
  long double n_readname; /* number of unique readname */
  long n_read_nonred;     /* number of nonredundant reads (use as "reads number") */
  long n_read_rpkm;       /* number of reads after normalization (n_read_nonred * rpkm weight) */
  long n_read_red;        /* number of redundant reads (filtered as PCR bias) */
  double n_read_afterGC;  /* number of reads after GC normalization */
};

typedef struct{
  struct seq seq[STRANDNUM];
  struct seq both;
  double depth;
  double w;
  
  /* FRiP */
  long n_read_inbed;
  double FRiP;
} SeqStats;

typedef struct{
  int nt_all, nt_nonred, nt_red; 
  double complexity;
  int tv;
} CompStats;

typedef struct{
  SeqStats *genome, *chr;
  Readarray **readarray;
  FragStats fstats;
  WigStats wstats;
  CompStats cs_raw, cs_nonred;
  int threshold4filtering;

  /* for GC*/
  int *GCdist;
  int maxGC;
  int sum_GCdist;
  double *GCweight;
} Mapfile;

#endif /* _PW_GV_H_ */
