/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_ESTIMATE_H_
#define _PW_ESTIMATE_H_

#include "pw_gv.h"
#include "seq.h"

#define N_ARRAY_NB 10

void pw_estimate_nb_param(Mapfile *mapfile, RefGenome *g);

double pw_get_negative_binomial(int k, double p, double n);
double pw_get_zeroinflated_negative_binomial(int, double, double, double);
double pw_get_poisson(int k, int lambda);
void pw_estimate_zinb_opt(Mapfile *mapfile, TYPE_WIGARRAY *wigarray, TYPE_WIGARRAY *mparray, int chr, int binnum);

#endif /* _PW_ESTIMATE_H_ */
