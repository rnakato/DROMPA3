/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DP_CALL_H_
#define _DP_CALL_H_

#include "drompa_gv.h"
#include "seq.h"

void dp_call(DrParam *p, SamplePair *sample, RefGenome *g, int chr);
void calc_ratio_and_pvalues(double *ratio, double *pval_inter, double *pval_enr, SamplePair *sample, int i);
int judge_significant(DrParam *p, double p_enr, double p_inter, double ratio_i, double iptag, char *Input_argv);

#endif /* _DP_CALL_H_ */
