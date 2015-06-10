/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _GSL_FUNC_H_
#define _GSL_FUNC_H_

#include "drompa_gv.h"

double calc_ratio(Samplefile *ChIP, Samplefile *Input, int i, double ratio);
double binomial_test(double n1_ref, double n2_ref, double ratio);
double zero_inflated_binomial_test(double k, double p, double n);

#endif /* _GSL_FUNC_H_ */
