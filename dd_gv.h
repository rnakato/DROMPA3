/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DD_GV_H_
#define _DD_GV_H_
#include <stdbool.h>
#include "seq.h"

#define STRUCT_GENE_MAX 10000
#define PDSIZE_DEFAULT 100000
#define LS_DEFAULT 1000
#define SCALE_TAG_DEFAULT 30
#define SCALE_RATIO_DEFAULT 5
#define SCALE_PVALUE_DEFAULT 5
#define SCALE_TAG_BROAD 50
#define SCALE_RATIO_BROAD 3

// for Graph
#define GC_MMIN 20
#define GC_MMAX 70
#define GD_MMIN 0
#define GD_MMAX 40

// page config
#define OFFSET_X 190
#define OFFSET_Y 50
#define WIDTH_DRAW 820
#define PAGEWIDTH  1088

#define FONTS_STYLE "Arial" //"sans-serif"

// for profile and heatmap

// for GV and PD
#define THRE_GENOMELEN4WG 20000000 // 20Mbp

// for page
#define BARNUM_DEFAULT 2
#define YSTEP_DEFAULT 20

typedef enum{
  ZERO,
  TSS,
  TTS,
  GENE100,
  BEDSITES,
  PTYPENUM
} Profile_Type;

typedef enum{
  READDIST,
  ENRICHDIST,
  PVALUEDIST
} ShowDist_Type;

typedef enum{
  GTYPE_SINGLE,
  GTYPE_MULTI,
  GTYPE_MULTI_PROP
} Graph_Type;

typedef enum{
  LTYPE_PVALUE_INTER,
  LTYPE_PVALUE_ENRICH,
  LTYPE_RATIO,
  LTYPE_RATIO_GV,
  LTYPE_CHIP,
  LTYPE_INPUT
} LINE_Type;

typedef struct{
  char *argv;
  char *name;
  int wsize;
  double *array;
  int arraynum;

  double mmin, mmax;
} Graph;

typedef enum{
  CODING,
  NONCODING,
  MIRNA,
  PSEUDO,
  PROCESS,
  OTHERS,
  ARS,
  TER,
  CENTROMERE,
  rRNA,
  LTR
} GeneAnnoType;

typedef struct{
  char name[128];
  char ID[128];
  int chr;
  char dir;
  int start;
  int end;
  int exonnum;
  SE *exon;
  GeneAnnoType genetype;
  bool delete;
} Gene;

typedef enum{
  GFTYPE_REFFLAT,
  GFTYPE_ENSEMBL,
  GFTYPE_GTF,
  GFTYPE_SGD
} GeneFile_Type;

typedef struct{
  char *argv;
  Gene *gene;
  int num, numall;
  GeneFile_Type type;
} GeneSet;

typedef enum{
  RepeatMasker,
  SimpleRepeats,
  Microsatellite,
  RM_SINE,
  RM_LINE,
  RM_LTR,
  RM_DNA,
  RM_Simple,
  RM_Low_Com,
  RM_Satellite,
  RM_RNA,
  RM_Other,
  REPEATTYPENUM
} RepeatType;

typedef struct{
  int start;
  int end;
  char dir;
  RepeatType rtype;
  char name[32];
  char class[32];
} Repeat;

typedef struct{
  char *argv;
  Repeat *repeat;
  int num;
} RepeatSet;

typedef struct{
  /* figure */
  int chronly;
  int makefig;
  char *command_mergepdf;
  int width_per_line;
  int linenum_per_page;
  int visualize_p_inter;
  int visualize_p_enrich;
  int visualize_ratio;
  int visualize_ctag;
  int visualize_itag;
  int barnum;
  int do_peakcall;
  int viz;
  int backcolors;
  int stroke_ymem;
  int stroke_ylab;
  int cl_orange;

  /* annotation*/
  Graph GC;
  Graph GD;
  Graph *PD;
  int pdnum;
  GeneSet gene;
  GeneFile_Type gftype;
  char *arsfile;
  char *terfile;
  bool *show_ars;
  RepeatSet repeat;
  int pdsize;
  int png;
  int rmchr;
  char *drawregion_argv;
  BedFile *drawregion;
  BedFile **bed;
  int bednum;
  InteractionSet *inter;
  int internum;

  /* page_info */
  int pageheight;
  double ystep;
  bool large_genome;
  bool gw;

  /* CG */
  double cgthre;
  /* PD */
  bool pd_prop;
  /* TR */
  double tssthre;

  /* profile and heatmap */
  char *filename_profile;
  Peak **peak;
  int ptype;
  int stype;
  int ntype;
  int compwidth;
  int cwbin;
  int ntotal_profile;
  int ntotal_skip;
  bool showse;
  int SEnum;
  int hmsort;
  bool sortgbody;
  bool pdetail;

  TYPE_WIGARRAY *mparray, *gaparray;
} DDParam;

#endif /* _DD_GV_H_ */
