/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include "argument.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "macro.h"

static const char *delete_arg(int *argcp, char *argv[], int i){
  const char *ret;
  if(i >= *argcp) return NULL;
  ret = argv[i];
  LOG("argv=%s\n", argv[i]);
  (*argcp)--;
  for(; i< *argcp; i++) argv[i] = argv[i+1];
  argv[i] = NULL;
  return (ret);
}

static int check_argv_dup(char *argv, char *argvarray[], int nargv, const Argument *args){
  int i;
  for(i=0; i<nargv; i++){
    if(!strcmp(argv, argvarray[i]) && !args->num){
      fprintf(stderr, "ERROR: option %s specified multiple times.\n", argv);
      exit(1);
    }
  }
  return 0;
}

static const Argument *check_arg(int argc, char *argv[], int i, const Argument *args){
  static char *argvarray[256];
  static int nargv=0;
  if(i >= argc) return NULL;
  for(; args->name; args++){
    if(!strcmp(argv[i], args->name)){
      check_argv_dup(argv[i], argvarray, nargv, args);
      argvarray[nargv++] = argv[i];
      return args;
    }
  }
  return NULL;
}

void argument_read(int *argcp, char *argv[], const Argument argument[]){
  int i,j;
  const char *p=NULL;
  const Argument *arg;

  for(i=1; i <*argcp; i++){
    while((arg = check_arg(*argcp, argv, i, argument))!=NULL){
      delete_arg(argcp, argv, i);
      switch(arg->type){
      case ARGUMENT_TYPE_NONE:
	/* NONE */
	break;
      case ARGUMENT_TYPE_FUNCTION:
	((void (*)())(arg->value))(); 
	break;
      case ARGUMENT_TYPE_FLAG_ON:
	*((int *)(arg->value)) = ARGUMENT_FLAG_ON;
	break;
      case ARGUMENT_TYPE_FLAG_OFF: 
	*((int *)(arg->value)) = ARGUMENT_FLAG_OFF;
	break;
      case ARGUMENT_TYPE_INTEGAR:
	p = delete_arg(argcp, argv, i);
	if(p){
	  if(p[0]=='-') goto err_novalue;
	  for(j=0; p[j] && isdigit((int)p[j]); j++);
	  if(p[j]!='\0') goto err_nodigit;
	  *((int *)(arg->value)) = atoi(p);
	}
	break;
      case ARGUMENT_TYPE_INTEGAR_MULTI:
	p = delete_arg(argcp, argv, i);
	if(p){
	  if(p[0]=='-') goto err_novalue;
	  (*((int **)(arg->value)))[(*arg->num)++] = atoi(p);
	}
	break;
      case ARGUMENT_TYPE_FLOAT:
	p = delete_arg(argcp, argv, i);
	if(p){
	  if(p[0]=='-') goto err_novalue;
	  for(j=0; p[j] && (isdigit((int)p[j]) || p[j]=='.'); j++);
	  if(p[j]!='\0') goto err_nodigit;
	  *((double *)(arg->value)) = atof(p);
	}
	break;
      case ARGUMENT_TYPE_STRING:
	p = delete_arg(argcp, argv, i);
	if(p){
	  if(p[0]=='-') goto err_novalue;
	  *((const char **)(arg->value)) = p;
	}
	break;
      case ARGUMENT_TYPE_STRING_MULTI:
	p = delete_arg(argcp, argv, i);
	if(p){
	  if(p[0]=='-') goto err_novalue;
	  (*((const char ***)(arg->value)))[(*arg->num)++] = p;
	}
	break;
      case ARGUMENT_TYPE_LOG10:
	p = delete_arg(argcp, argv, i);
	if(p){
	  if(p[0]=='-') goto err_novalue;
	  for(j=0; p[j] && (isdigit((int)p[j]) || p[j]=='.'); j++);
	  if(p[j]!='\0') goto err_nodigit;
	  *((double *)(arg->value)) = -log10(atof(p));
	}
	break;
      }
    }
  }

  int on=0;
  while(*argcp > 1){
    fprintf(stderr, "ERROR: Invalid Argument: %s\n", argv[1]);
    delete_arg(argcp, argv, 1);
    on=1;
  }
  if(on) exit(1);

  return;
 err_novalue:
    fprintf(stderr, "ERROR: no value specified for %s.\n", arg->name);
    exit(1);
 err_nodigit:
    fprintf(stderr, "ERROR: Invalid value for %s: \"%s\".\n", arg->name, p);
    exit(1);
}
