/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#ifndef _STRINGP_H_
#define _STRINGP_H_

#define STR_LEN 100000
#define ELEM_NUM 256
#define MAXLEN_INSCOMMA 20

typedef struct{
  char str[10240];
} Elem;

void chomp(char *);
int ParseLine(char *, Elem clm[]);
int ParseLine_arbit(char *, Elem clm[], char);
char *delimit_str(char *str, char token);
char *alloc_str_new(char *str, int add);
void insComma(int num, char str[MAXLEN_INSCOMMA]);

#endif /* _STRINGP_H_ */
