/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _READFILE_H_
#define _READFILE_H_

#include "seq.h"
#include "common.h"

char *read_mosaics_binary(char *filename, long chrlen);
TYPE_WIGARRAY *read_mosaics_binfile(char *filename, int, int);
void parse_genometable(char *, RefGenome *);
void scan_mappability_chrtotal(char *, RefGenome *);
int changechr_str2int(char *, RefGenome *);
BedFile *read_bedfile(char *, RefGenome *);
void show_bedfile(BedFile *p, int);
char *makearray_inbed(BedChr *p, int arraysize);
Peak *read_peakfile(char *filename, RefGenome *g);
void read_interactionfile(InteractionSet *inter, char *argv, RefGenome *g);

#endif /* _READFILE_H_ */
