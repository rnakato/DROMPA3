/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _PW_MAKEFILE_H_
#define _PW_MAKEFILE_H_

#include "seq.h"
#include "pw_gv.h"

void makewig(PwParam *p, Mapfile *mapfile, RefGenome *g);
void calc_FRiP(PwParam *p, Mapfile *mapfile, RefGenome *g);
#endif /* _PW_MAKEFILE_H_ */
