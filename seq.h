/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
/* 一般的なデータ形式を格納する構造体の宣言 */
#ifndef _SEQ_H_
#define _SEQ_H_

#define BEDFILE_MAX 1000

typedef enum{
  STRAND_PLUS,
  STRAND_MINUS,
  STRANDNUM
} Strand;

extern const char str_strand[][STRANDNUM];

typedef struct{
  int start;
  int end;
} SE;

typedef struct{
  char *name;
  long len, len_mpbl;
  double p_mpbl;  /* mappability */
  double gcov;    /* genome coverage for bin */
} Fastadata;

typedef struct{
  Fastadata *genome;
  Fastadata *chr;
  int chrnum;
  int chrmax;
  int *GCdist;
  long sum_GCdist;
} RefGenome;

typedef struct{
  int r;
  int g;
  int b;
} Rgbval;

/* bed format */
struct bed{
  /* required fields */
  int s,e;
  /* optional fields */
  char *name;
  double score;
  Strand strand;
  int ts, te;     /* thickStart, thickEnd (used by UCSC drawing code) */
  Rgbval itemRgb; /* RGB colour value */
};

typedef struct{
  struct bed *bed;
  int num;
  int arraynum;
} BedChr;

typedef struct{
  BedChr *chr;
  int nline;
  int num;
  long len_total;
  char *name;
  char *argv;
} BedFile;


RefGenome *refgenome_new();
void refgenome_delete(RefGenome *);
BedFile *bedfile_new(int);
void bedfile_delete(BedFile *p, int);

#endif /* _SEQ_H_ */
