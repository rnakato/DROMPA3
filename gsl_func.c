/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <math.h>
#include <gsl/gsl_cdf.h>
#include "gsl_func.h"

/*static double calc_locallambda(DrParam *p, Samplefile *a, int i, int binnum){
  int num = (p->width4lambda/p->binsize)>>1;
  int s = max(i - num, 0);
  int e = min(i + num, binnum);
  int j, d=0;
  for(j=s; j<=e; j++) d += a->data[j];
  return WIGARRAY2VALUE(d/(double)(e-s+1));
  }*/

double binomial_test(double n1_ref, double n2_ref, double ratio){
  int n1, n2;
  double p=0.5;  // null model
  double pvalue;
  if(ratio > 1){  /* 大きい方を小さい方に合わせる */
    n1 = (int)ceil(n1_ref/ratio); // rounded up
    n2 = (int)ceil(n2_ref);
  }else{
    n1 = (int)ceil(n1_ref);
    n2 = (int)ceil(n2_ref*ratio);
  }
  if((n1 < n2) || (!n1 && !n2)) return 0;
  
  pvalue = gsl_cdf_binomial_Q(n1, p, n1+n2);
  if(!pvalue) pvalue = 1e-314;
  pvalue = -log10(pvalue);
  return pvalue;
}

double zero_inflated_binomial_test(double k, double p, double n){
  double r, pval;
  // static double pmin=1;
  if(!k) pval =0;
  else {
    pval = gsl_cdf_negative_binomial_Q(k, p, n);
    if(!pval) pval = 1e-314;
    //    else if(pmin > pval) pmin = pval;
    //printf("k=%f, p=%f, n=%f, pval=%.3e, pmin=%.3e\n",k, p, n, pval, pmin);
  }
  if(!pval) r = 0; else r = -log10(pval);
  return r;
}
