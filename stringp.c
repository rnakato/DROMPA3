/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */

#include <string.h>
#include "alloc.h"
#include "stringp.h"
#include "memtest.h"

void chomp(char *str){
  char *p = strchr(str, '\n');
  if(p) p[0]='\0';
  return;
}

int ParseLine(char *str, Elem clm[]){
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len, sizeof(char), "ParseLine");
  for(i=0; i<=len; i++){
    if(str[i]=='\0' || str[i]=='\n'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      MYFREE(strtemp);
      return ++num;
    }
    if(str[i]=='\t'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      num++; 
      if(num >= ELEM_NUM){
	fprintf(stderr, "error: too many columns: %s", str);
	exit(0);
      }
      j=0;
    }else{
      strtemp[j]=str[i];
      j++;
    }
  }
  MYFREE(strtemp);
  return num;
}

int ParseLine_arbit(char *str, Elem clm[], char token){
  int i, j=0, num=0, len=strlen(str);
  char *strtemp = (char *)my_calloc(len+1, sizeof(char), "strtemp");
  for(i=0; i<=len; i++){
    if(str[i]=='\0'){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      MYFREE(strtemp);
      return ++num;
    }
    if(str[i]==token){
      strtemp[j]='\0';
      strcpy(clm[num].str, strtemp);
      num++; 
      j=0;
    }else{
      strtemp[j]=str[i];
      j++;
    }
  }
  MYFREE(strtemp);
  return num;
}

char *delimit_str(char *str, char token){
  int i, len=strlen(str);
  for(i=0; i<=len; i++){
    if(str[i]==token){
      str[i]='\0';
      break;
    }
  }
  return str;
}

char *alloc_str_new(char *str, int add){
  char *p=NULL;
  int length = strlen(str);
  //  printf("str=%s, len=%d, len+add=%d\n", str, length, length + add);
  p = (char *)my_calloc(length + add, sizeof(char), "alloc_str_new");
  return p;
}

void insComma(int num, char str[MAXLEN_INSCOMMA]){
  int i, j, keta, temp;
  int minus = 0;
  
  /* 下位から順に数字を取り出し3桁の区切りに,を入れる */
  i = keta = 0;
  if(num <0){
    num = -num;
    minus = 1;
  }
  do{
    str[i++] = num % 10 + '0';  // <-- 今回のクイズ
    keta++;
    num /= 10;
    if(keta %3 == 0 && num != 0) str[i++] = ',';
  }while(num != 0 && i < MAXLEN_INSCOMMA);
  if(minus) str[i++] = '-';
  str[i] = '\0';

  /* 逆順にコピーする */
  j = i-1;
  for(i=0; j>i; i++,j--) {
    temp = str[i];
    str[i] = str[j];
    str[j] = temp;
  }
  return;
}
