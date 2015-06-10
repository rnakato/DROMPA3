/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _FILEHANDLE_H_
#define _FILEHANDLE_H_

#include <stdio.h>

typedef enum{
  FILE_MODE_READ,
  FILE_MODE_WRITE,
  FILE_MODE_WB,
  FILE_MODE_A
} File_Mode;

void isfile(char *file);
FILE *my_fopen(char *filename, File_Mode mode);
void my_system(char *command);
void remove_file(char *filename);

#endif /* _FILEHANDLE_H_ */
