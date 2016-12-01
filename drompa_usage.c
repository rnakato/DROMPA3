/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include "drompa_usage.h"
#include "dd_gv.h"
#include "readfile.h"
#include "macro.h"

char str_ftype[][16]={"PC_SHARP", "PC_BROAD", "PC_ENRICH", "GV", "PD", "FRIP", "CI", "CG", "GOVERLOOK", "PROFILE", "HEATMAP", "TR"};

static void print_gt_and_output(){
  fprintf(stderr, "       <output>: Name of output files\n");
  fprintf(stderr, "       <genometable>: Tab-delimited file describing the name and length of each chromosome\n");
}

static void print_usage_if_binsize(int binsize){
  fprintf(stderr, "       -if <int>: input file format (default: 0)\n");
  fprintf(stderr, "            0: binary file (.bin)\n");
  fprintf(stderr, "            1: compressed wig file (.wig.gz)\n");
  fprintf(stderr, "            2: uncompressed wig file (.wig)\n");
  fprintf(stderr, "            3: bedGraph (.bedGraph)\n");
  fprintf(stderr, "       -binsize <int>  : bin size (default: %d bp)\n", binsize);
}

static void print_usage_io(){
  fprintf(stderr, "       -norm <int>     : normalization between ChIP and Input\n");
  fprintf(stderr, "            0; not normalize\n");
  fprintf(stderr, "            1; with total read number (default)\n");
  fprintf(stderr, "            2; with NCIS method\n");
  fprintf(stderr, "       -includeYM      : output peaks of chromosome Y and M (default: off)\n");
  fprintf(stderr, "       -sm <int> : smooth bin with <int> bp window (default: %d bp)\n", SMOOTHING_DEFAULT);
}

static void print_mpbl(){
  fprintf(stderr, "       -mp <mappability file>: for normalization of mappability\n");
  fprintf(stderr, "       -mpthre <double>: threshold of low mappability (default: %.1f)\n", THRE_LOW_MAPPABILITY);
  fprintf(stderr, "       -gap <gap file>: shade highly gapped regions\n");
}

static void print_usage_threshold(){
  fprintf(stderr, "\n   threshold:\n");
  fprintf(stderr, "       -pthre_enrich   <double>: p-value threshold for ChIP internal (Poisson, default < %.1e)\n", PTHRE_INTERNAL_DEFAULT);
  fprintf(stderr, "       -pthre_internal <double>: p-value threshold for ChIP/Input enrichment (binomial, default < %.1e)\n", PTHRE_ENRICH_DEFAULT);
  fprintf(stderr, "       -ethre <double>: IP/Input enrichment threshold (default: %d)\n", ENRICHTHRE_DEFAULT);
  fprintf(stderr, "       -ipm <double>: read intensity threshold of peak summit (default: 0)\n");
}


void print_usage_dp(){
  //  fprintf(stderr, "DROMPA version %s\n", VERSION);
  fprintf(stderr, "Usage: drompa_peakcall COMMAND [options] -p <output> -gt <genometable> -i <ChIP>,<Input>\n");
  fprintf(stderr, "       <ChIP>\tChIP sample\n");
  fprintf(stderr, "       <Input>\tInput (control) sample\n");
  print_gt_and_output();
  fprintf(stderr, "\nOptions:\n");
  print_usage_if_binsize(BINSIZE_DEFAULT);
  fprintf(stderr, "       -width4lmd <int> : width for calculating local lambda (default: %d bp)\n", WIDTH4LAMBDA_DEFAULT);
  print_usage_io();
  fprintf(stderr, "       -outputwig <int> : output bin data\n");
  fprintf(stderr, "            0; ChIP/Input ratio\n");
  fprintf(stderr, "            1; output ChIP-internal p-value\n");
  fprintf(stderr, "            2; ChIP/Input enrichment p-value\n");
  fprintf(stderr, "       -owtype <int>: output format (default: 0)\n");
  fprintf(stderr, "           0: binary (.bin) \n");
  fprintf(stderr, "           1: compressed wig (.wig.gz) \n");
  fprintf(stderr, "           2: uncompressed wig (.wig) \n");
  fprintf(stderr, "           3: bedGraph (.bedGraph) \n");
  fprintf(stderr, "           4: bigWig (.bw) \n");
  fprintf(stderr, "       -odir: output directory name (default: 'drompadir') \n");
  print_usage_threshold();
  fprintf(stderr, "       -qthre <double>: q-value threshold (Bonferroni-Hochberg, default < %.3f)\n", QTHRE_DEFAULT);
  fprintf(stderr, "\n   annotations:\n");
  print_mpbl();
  fprintf(stderr, "\n");
  exit(0);
}

static void print_graph(){
  fprintf(stderr, "       -GC <filename>: draw GCcontents graph\n");
  fprintf(stderr, "       -gcsize <int> : window size for -GC (default: %d kb)\n", GCSIZE_DEFAULT/NUM_1K);
  fprintf(stderr, "       -GD <filename>: draw gene density (number of genes per each window)\n");
  fprintf(stderr, "       -gdsize <int> : window size for -GD (default: %d kb)\n", GDSIZE_DEFAULT/NUM_1K);
  return;
}
static void print_gene_anno(){
  fprintf(stderr, "       -gene <gene.txt>: gene annotation file\n");
  fprintf(stderr, "       -gftype <int>   : format of <gene.txt>\n");
  fprintf(stderr, "            0: RefFlat (default)\n");
  fprintf(stderr, "            1: Ensembl\n");
  fprintf(stderr, "            2: GTF (for S. pombe)\n");
  fprintf(stderr, "            3: SGD (for S. cerevisiae)\n");
  return;
}

static void print_usage_annotations(){
  print_gene_anno();
  fprintf(stderr, "       -ars <ARS.txt>: ARS list (for yeast)\n");
  fprintf(stderr, "       -showars: display ARS only (do not display genes)\n");
  fprintf(stderr, "       -ter <TER.txt>: TER list (for S.cerevisiae)\n");
  fprintf(stderr, "       -bed <bedfile>,<name>: specify the bedfile and name (<name> can be omited)\n");
  fprintf(stderr, "       -repeat <repeat.txt>: display repeat annotation (RepeatMasker format)\n");
}

static void print_usage_drawing1(){
  fprintf(stderr, "       -nosig: do not call peaks\n");
  fprintf(stderr, "       -png: output with png format (Note: output each page separately)\n");
  fprintf(stderr, "       -ls <int>: length per one line (kb) (default: %d kb)\n", LS_DEFAULT);
  fprintf(stderr, "       -lpp <int>: number of lines per one page (default: 1)\n");
  fprintf(stderr, "       -chr <int>: output only a chromosome specified\n");
  fprintf(stderr, "       -rmchr: remove chromosome-separated pdf files\n");
  fprintf(stderr, "       -r <regionfile>: specify the regions to visualize\n");
  fprintf(stderr, "       -genefile <genelist>: specify the gene loci to visualize\n");
  fprintf(stderr, "       -len_genefile <int>: extended length for each gene locus (default: 50000 bp)\n\n");
}

static void print_usage_drawing2(){
  fprintf(stderr, "       -bn <int>: number of separations for y-axis (default: %d)\n", BARNUM_DEFAULT);
  fprintf(stderr, "       -ystep <double>: height of read line (default: %d)\n", YSTEP_DEFAULT);
  fprintf(stderr, "       -show_ctag <int>: display ChIP-read lines   (0:off 1:on)\n");
  fprintf(stderr, "       -show_itag <int>: display Input-read lines  (0:off 1:all 2:first one)\n");
  fprintf(stderr, "       -showratio <int>: display ChIP/Input ratio (0:off 1:liner scale 2:logscale)\n");
  fprintf(stderr, "       -showpinter <int>:  display logp for ChIP internal (negative binomial test)\n");
  fprintf(stderr, "       -showpenrich <int>: display logp for ChIP/Input enrichment (binomial test)\n");
  fprintf(stderr, "       -scale_tag <double>:    scale of read line\n");
  fprintf(stderr, "       -sct<int> <double>:     scale of <int>st read line\n");
  fprintf(stderr, "       -scale_ratio <double>:  scale of ratio line\n");
  fprintf(stderr, "       -scale_pvalue <double>: scale of pvalue line\n");
  fprintf(stderr, "       -scr<int> <double>:     scale of <int>st ratio line\n");
  fprintf(stderr, "       -viz <int>: color of read profile\n");
  fprintf(stderr, "            0(default): normal color\n");
  fprintf(stderr, "            1: semitransparent color\n");
  fprintf(stderr, "       -offylab:   delete Y label\n");
  fprintf(stderr, "       -offymem:   delete Y memory\n");
  fprintf(stderr, "       -offbg:     delete background color of read lines\n\n");
}

static void print_usage_CG_POL2(Function_Type ftype){
  fprintf(stderr, "Usage: drompa_draw %s -p <output> -gt <genometable> -i <ChIP>,,<name> [-i <ChIP>,,<name> ...]\n", str_ftype[ftype]);
  print_gt_and_output();
  fprintf(stderr, "       -i: specify ChIP data and name (separated by ',,', <name> can be omitted)\n\n");
  fprintf(stderr, "Options:\n");
  print_usage_if_binsize(BINSIZE_DEFAULT);
  print_gene_anno();
  if(ftype == FTYPE_COMPARE_GENEBODY) fprintf(stderr, "       -cgthre <double>: minimum threshold per kbp\n");
  fprintf(stderr, "       -outputYM       : output peaks of chromosome Y and M (default: ignore)\n\n");
}

static void print_usage_CI_FRIP(Function_Type ftype){
  fprintf(stderr, "Usage: drompa_draw %s -p <output> -gt <genometable> -i <ChIP>,,<name> -i <ChIP>,,<name> -bed <bedfile>\n", str_ftype[ftype]);
  print_gt_and_output();
  fprintf(stderr, "       -i: specify ChIP data and name (separated by ',,', <name> can be omitted)\n");
  fprintf(stderr, "       -bed: region file (e.g., peak list)\n\n");
  fprintf(stderr, "Options:\n");
  print_usage_if_binsize(BINSIZE_DEFAULT);
  fprintf(stderr, "\n");
}

static void print_usage_PD(){
  fprintf(stderr, "Usage: drompa_draw PD [options] -p <output> -gt <genometable> -pd <pdfile>,<name> [-pd <pdfile>,<name> ...]\n");
  print_gt_and_output();
  fprintf(stderr, "       -pd: specify peak density file and name (separated by ',' <name> can be omited)\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "       -pdsize: windowsize of peak density (default: %d kbp)\n", PDSIZE_DEFAULT/NUM_1K);
  print_graph();
  fprintf(stderr, "\n");
}

static void print_usage_GOVERLOOK(){
  fprintf(stderr, "Usage: drompa_draw GOVERLOOK -p <output> -gt <genometable> -bed <bedfile>,<name> [-bed <bedfile>,<name> ...]\n");
  print_gt_and_output();
  fprintf(stderr, "       -bed: specify the peak file and name (separated by ',' <name> can be omited)\n\n");
}


void print_usage_dd(Function_Type ftype){
  fprintf(stderr, "Usage: drompa_draw %s [options] -p <output> -gt <genometable> -i <ChIP>,<Input>,<name> [-i <ChIP>,<Input>,<name> ...]\n", str_ftype[ftype]);
  print_gt_and_output();
  fprintf(stderr, "       -i: specify ChIP data, Input data and name of ChIP sample (separated by ',', <Input> and <name> can be omitted)\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\n   Input and output:\n");
  int binsize;
  if(ftype==FTYPE_GV) binsize = BINSIZE_DEFAULT_GV; else binsize = BINSIZE_DEFAULT;
  print_usage_if_binsize(binsize);
  if(ftype==FTYPE_PEAKCALL_SHARP || ftype==FTYPE_PEAKCALL_BROAD || ftype==FTYPE_PEAKCALL_E){
    print_usage_io();
    //    print_mpbl();
    print_usage_threshold();
    fprintf(stderr, "\n   annotations:\n");
    print_usage_annotations();
    print_graph();
    fprintf(stderr, "\n   drawing parameter:\n");
    print_usage_drawing1();
    print_usage_drawing2();
  }else if(ftype==FTYPE_GV){
    fprintf(stderr, "       -bs<int1> <int2>: change bin size of <int1>st sample\n");
    print_usage_io();
    fprintf(stderr, "\n   annotations:\n");
    print_graph();
    fprintf(stderr, "\n   drawing parameter:\n");
    print_usage_drawing2();
  }else if(ftype==FTYPE_PROFILE || ftype==FTYPE_HEATMAP){
    print_gene_anno();
    fprintf(stderr, "       -bed <peakfile>:  peaklist (multiple files allowed)\n\n");
    fprintf(stderr, "       -cw <int>: width from the center (default: 2.5kb)\n");
    fprintf(stderr, "       -norm <int>     : normalization between ChIP and Input\n");
    fprintf(stderr, "            0; not normalize\n");
    fprintf(stderr, "            1; with total read number (default)\n");
    fprintf(stderr, "            2; with NCIS method\n");
    fprintf(stderr, "       -ptype <int>: 1; around TSS, 2; around TES, 3; divide gene into 100 subregions 4; around peak sites\n");
    fprintf(stderr, "       -stype <int>: show type (0; ChIP read (default) 1; ChIP/Input enrichment)\n");
    if(ftype==FTYPE_PROFILE){
      fprintf(stderr, "       -ntype <int>: normalization type (0; not normalized (default) 1; normalized by number of read mapped to target regions)\n");
      fprintf(stderr, "       -offse: omit the standard error in profile.\n\n");
    }
    if(ftype==FTYPE_HEATMAP){
      fprintf(stderr, "       -scale_tag    <double>: maxvalue for -stype0\n");
      fprintf(stderr, "       -scale_ratio  <double>: maxvalue for -stype1\n");
      fprintf(stderr, "       -scale_pvalue <double>: maxvalue for -stype2\n\n");
      fprintf(stderr, "       -hmsort <int>: sort lists for <int>th sample \n\n");
    }
  }
}

void print_error_dd(Function_Type ftype){
  if(ftype==FTYPE_PD) print_usage_PD();
  else if(ftype==FTYPE_FRIP || ftype==FTYPE_COMPARE_INTENSITY) print_usage_CI_FRIP(ftype);
  else if(ftype==FTYPE_COMPARE_GENEBODY || ftype==FTYPE_TR)  print_usage_CG_POL2(ftype);
  else if(ftype==FTYPE_GOVERLOOK)         print_usage_GOVERLOOK();
  else print_usage_dd(ftype);
  exit(0);
}

#define PRINT_ERROR(args...) do{ printf(args); return 1; }while(0)

int print_error_peakcall(DrParam *p, RefGenome *g){
  if(!p->headname)    PRINT_ERROR("please specify output prefix.\n");
  if(p->mpfile && p->mpthre <=0) PRINT_ERROR("error: threshold for mappability = %.2f.\n", p->mpthre);
  if(!p->gtfile)      PRINT_ERROR("please specify genome_table.\n");
  else parse_genometable(p->gtfile, g);
  if(!range(p->itype, 0, PWFILETYPENUM-2))        PRINT_ERROR("error: Invalid input filetype: %d.\n", p->itype);
  if(!range(p->ntype, 0, CHIPNORMALIZETYPENUM-1)) PRINT_ERROR("error: Invalid input filetype: %d.\n", p->ntype);
  if(p->smoothing <0) PRINT_ERROR("error: smoothing size = %d.\n", p->smoothing);
  if(p->IPmaxthre <0) PRINT_ERROR("error: -ipm = %.2f < 0.\n", p->IPmaxthre);

  return 0;
}
