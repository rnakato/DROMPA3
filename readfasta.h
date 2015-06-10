/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _READFASTA_H_
#define _READFASTA_H_

#define TYPE_FASTAGCARRAY short

TYPE_FASTAGCARRAY *make_fastaGCarray(char *filename, long length, int flen4gc);

#endif /* _READFASTA_H_ */
