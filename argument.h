/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _ARGUMENT_H_
#define _ARGUMENT_H_

#include <stdbool.h>

typedef enum{
  ARGUMENT_TYPE_NONE,
  ARGUMENT_TYPE_FUNCTION,
  ARGUMENT_TYPE_FLAG_ON,
  ARGUMENT_TYPE_FLAG_OFF,
  ARGUMENT_TYPE_INTEGAR,
  ARGUMENT_TYPE_INTEGAR_MULTI,
  ARGUMENT_TYPE_FLOAT,
  ARGUMENT_TYPE_STRING,
  ARGUMENT_TYPE_STRING_MULTI,
  ARGUMENT_TYPE_LOG10,
} Argument_Type;

#define ARGUMENT_FLAG_ON 1
#define ARGUMENT_FLAG_OFF 0

typedef struct{
  const char *name;
  Argument_Type type;
  void *value;
  int *num;
} Argument;

void argument_read(int *argcp, char *argv[], const Argument args[]);

#endif /* _ARGUMENT_H_ */
