/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_PARAM_NEW_H_
#define _PW_PARAM_NEW_H_

#include "pw_gv.h"

PwParam *pwparam_new();
Mapfile *mapfile_new(int chrnum, PwParam *p);
void mapfile_delete(Mapfile *mapfile, int);
void pwparam_delete(PwParam *p, RefGenome *);
#endif /* _PW_PARAM_NEW_H_ */
