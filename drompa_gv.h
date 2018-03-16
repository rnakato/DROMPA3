/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DROMPA_GV_H_
#define _DROMPA_GV_H_
#include <stdbool.h>
#include "common.h"
#include "seq.h"

#define BINSIZE_DEFAULT 100
#define BINSIZE_BROAD 1000
#define BINSIZE_DEFAULT_GV 100000
#define SMOOTHING_DEFAULT 500
#define SMOOTHING_BROAD 10000
#define WIDTH4LAMBDA_DEFAULT 100000
#define PEAKNUM_DEFAULT 10000

#define GCSIZE_DEFAULT 500000
#define GDSIZE_DEFAULT 500000

#define ENRICHTHRE_DEFAULT 2
#define PTHRE_INTERNAL_DEFAULT 1e-4
#define PTHRE_ENRICH_DEFAULT 1e-4
#define QTHRE_DEFAULT 1e-3

typedef enum{
  FTYPE_PEAKCALL_SHARP,
  FTYPE_PEAKCALL_BROAD,
  FTYPE_PEAKCALL_E,
  FTYPE_GV,
  FTYPE_PD,
  FTYPE_FRIP,
  FTYPE_COMPARE_INTENSITY,
  FTYPE_MULTICI,
  FTYPE_COMPARE_GENEBODY,
  FTYPE_GOVERLOOK,
  FTYPE_PROFILE,
  FTYPE_HEATMAP,
  FTYPE_TR
} Function_Type;

typedef enum{
  OWTYPE_NONE,
  OWTYPE_RATIO,
  OWTYPE_P_INTER,
  OWTYPE_P_ENRICH,
} OUTPUTWIG_Type;

typedef enum{
  TYPE_RATIO_NONE,
  TYPE_RATIO_TOTALREAD,
  TYPE_RATIO_NCIS,
  CHIPNORMALIZETYPENUM
} ChIPNormalize_Type;

typedef struct{
  int samplenum, samplenum_1st, samplenum_overlay;
  char *headname;
  int smoothing;
  int width4lambda;
  int flen;

  Function_Type ftype;
  PWfile_Type itype;
  ChIPNormalize_Type ntype;

  char *gtfile;
  char *mpfile;
  double mpthre;
  char *gapfile;

  bool includeYM;
  int outputwig;
  PWfile_Type wtype;
  int showzero;

  char *output_dir;

  double IPmaxthre;
  double enrichthre;
  double pthre_internal, pthre_enrich;
  double qthre;

  /*ignore region*/
  BedFile **igregion;
  int n_igregion;

} DrParam;

typedef struct{
  double ratio;
} LCstats;

typedef struct _LibCompare{
  LCstats *genome, *chr;
} LibCompare;

typedef struct{
  long nread;
} Library;

typedef struct _Samplefile{
  char *argv;
  Library *genome, *chr;
  TYPE_WIGARRAY *data;

  double lambda;
  double nb_p, nb_n, nb_p0;

} Samplefile;

typedef struct{
  double *IP;
  double *Input;
  double **SE, *SEarray;
  double IPsum, Inputsum;
} Profile;

typedef struct{
  Samplefile *ChIP, *Input;
  LibCompare *comp;
  int copyC, copyI, copycomp;
  int overlay;

  Peak *peak;
  char *peak_argv;
  char *peakarray;
  char *linename;
  int binsize;
  int *binnum;

  double scale_tag;
  double scale_ratio;
  double scale_pvalue;

  Profile profile;

} SamplePair;

extern char str_ftype[][17];

#endif /* _DROMPA_GV_H_ */
