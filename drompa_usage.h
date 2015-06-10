/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
/* drompaとparse2wigで共有するパラメータを格納 */
#ifndef _DROMPA_USAGE_H_
#define _DROMPA_USAGE_H_

#include "drompa_gv.h"
#include "seq.h"

void print_usage_dp();
void print_usage_dd(Function_Type ftype);
void print_error_dd(Function_Type ftype);
int print_error_peakcall(DrParam *p, RefGenome *g);

#endif /* _DROMPA_USAGE_H_ */
