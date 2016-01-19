/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
/* drompaとparse2wigで共有するパラメータを格納 */
#ifndef _COMMON_H_
#define _COMMON_H_

#ifdef CLOCK
#include <time.h>
#endif

#define VERSION "3.2.0"
#define BINSIZE_DEFAULT 100
#define THRE_LOW_MAPPABILITY 0.3

#define NUM_1K 1000
#define NUM_1M 1000000
#define NUM_10M 10000000
#define NUM_100M 100000000

/* wigfile */
#define TYPE_WIGARRAY int
#define VALUE2WIGARRAY(v) ((v) * 1000.0)
#define WIGARRAY2VALUE(v) ((v) * (1.0/1000.0))

typedef enum{
  TYPE_BINARY,
  TYPE_COMPRESSWIG,
  TYPE_UNCOMPRESSWIG,
  TYPE_BEDGRAPH,
  TYPE_BIGWIG,
  PWFILETYPENUM
} PWfile_Type;

struct bs{
  int chr;
  int start;
  int end;
  int maxposi;
  double maxIP;
  double enrich;
  double p_inter, p_enr;
  double qvalue;
};

typedef struct{
  char *argv;
  char *name;
  struct bs *bs;
  int num;
  int arraynum;
  int num_nonzero;
} Peak;

typedef struct{
  struct bs head, tail;
  double pval;
} Interaction;

typedef struct{
  char *argv;
  char *name;
  Interaction *site;
  int num;
  int arraynum;
} InteractionSet;


#endif /* _COMMON_H_ */
