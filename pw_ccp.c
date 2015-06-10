/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <string.h>
#include <math.h>
#include "pw_ccp.h"
#include "alloc.h"
#include "filehandle.h"
#include "stringp.h"

static TYPE_WIGARRAY calc_qnt(TYPE_WIGARRAY *array, int chrlen, double per);
static double calc_corr(TYPE_WIGARRAY* x, TYPE_WIGARRAY* y, int);
static int comp_wig(const void *c1, const void *c2);

static TYPE_WIGARRAY *chrarray_new(Mapfile *mapfile, int chr, int chrlen, Strand strand){
  int i;
  Readarray *read = &(mapfile->readarray[chr][strand]);
  int readnum = mapfile->chr[chr].seq[strand].n_read_infile;
  TYPE_WIGARRAY *array = (TYPE_WIGARRAY *)my_calloc(chrlen, sizeof(TYPE_WIGARRAY), "chrarray");
  for(i=0; i<readnum; i++){
    if(read->delete[i]) continue;
    array[read->F3[i]]++;
  }
  return array;
}

void pw_ccp(PwParam *p, Mapfile *mapfile, int chr, int chrlen){
  int i;

  printf("Making cross-correlation profile...\n");
  TYPE_WIGARRAY *plus  = chrarray_new(mapfile, chr, chrlen, STRAND_PLUS);
  TYPE_WIGARRAY *minus = chrarray_new(mapfile, chr, chrlen, STRAND_MINUS);

  TYPE_WIGARRAY qnt99 = calc_qnt(plus, chrlen, 0.99);
  for(i=0; i<chrlen; i++){
    if(plus[i]  > qnt99) plus[i]  = qnt99;
    if(minus[i] > qnt99) minus[i] = qnt99;
  }

  char *outputfile = alloc_str_new(p->output_dir, strlen(p->output_prefix) +100);
  sprintf(outputfile, "%s/%s.ccp.xls", p->output_dir, p->output_prefix);
  FILE *OUT = my_fopen(outputfile, FILE_MODE_WRITE);
  fprintf(OUT, "strand-shift\tcross-correlation\n");
  double cc=0;
  int start=1000;
  int num=chrlen-2500;
  for(i=-500; i<1500; i+=5){
    cc = calc_corr(plus + start, minus + start + i, num);
    fprintf(OUT, "%d\t%f\n", i, cc);
  }
  fclose(OUT);
  printf("Output to %s.\n",outputfile);

  MYFREE(outputfile);
  MYFREE(plus);
  MYFREE(minus);
  return;
}

static TYPE_WIGARRAY calc_qnt(TYPE_WIGARRAY *array, int chrlen, double per){
  int i, num=0;
  TYPE_WIGARRAY qnt=0;
  TYPE_WIGARRAY *sortarray = (TYPE_WIGARRAY *)my_calloc(chrlen, sizeof(TYPE_WIGARRAY), "sortarray");
  for(i=0; i<chrlen; i++){
    if(array[i]) sortarray[num++] = array[i];
  }
  qsort(sortarray, num, sizeof(TYPE_WIGARRAY), comp_wig);
  qnt = sortarray[(int)(num*per)];
  MYFREE(sortarray);
  return qnt;
}

static double calc_mean(int n, TYPE_WIGARRAY *a){
  double mean = 0;
  int i;
  for(i=0; i<n; i++) mean += a[i];
  mean /= (double)n;
  return mean;
}

static double calc_corr(TYPE_WIGARRAY *x, TYPE_WIGARRAY *y, int num){
  int i;
  double dx, dy;

  static double mx=-1, xx=0;
  if(mx==-1){
    mx = calc_mean(num, x);
    for(i=0; i<num; i++){
      dx = x[i] - mx;
      xx += dx * dx;
    }
    xx = sqrt(xx);
  }

  double yy=0, xy=0;
  static int sumy=-1;
  if(sumy==-1){
    sumy=0;
    for(i=0; i<num; i++) sumy += y[i];
  }else{
    for(i=0; i<5; i++){
      sumy -= *(y-i);
      sumy += y[num-i-1];
    }
  }
  double my = sumy/(double)num;
  for(i=0; i<num; i++){
    dx = x[i] - mx;
    dy = y[i] - my;
    yy += dy * dy;
    xy += dx * dy;
  }
  yy = sqrt(yy);

  return(xy / (xx*yy));
}

static int comp_wig(const void *c1, const void *c2){
  TYPE_WIGARRAY n1 = *(TYPE_WIGARRAY *)c1;
  TYPE_WIGARRAY n2 = *(TYPE_WIGARRAY *)c2;
  if(n1<n2) return -1;
  else if(n1==n2) return 0;
  else return 1;
}
