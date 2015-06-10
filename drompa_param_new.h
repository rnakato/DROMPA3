/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _DR_PARAM_NEW_H_
#define _DR_PARAM_NEW_H_

#include "drompa_gv.h"

DrParam *drparam_new();
void drparam_delete(DrParam *p);
SamplePair *SamplePair_new(int num, int chrnum);
void SamplePair_delete(SamplePair *p, int num);

#endif /* _DR_PARAM_NEW_H_ */
