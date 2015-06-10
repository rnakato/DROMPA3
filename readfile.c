/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "readfile.h"
#include "filehandle.h"
#include "stringp.h"
#include "macro.h"
#include "memtest.h"

#define READFILE_PEAKNUM 10000

char *read_mosaics_binary(char *filename, long chrlen){
  isfile(filename);
  //  printf("filename: %s\n",filename);
  int num=0;
  char c;
  char *array = (char *)my_calloc(chrlen, sizeof(char), "mosaicsbinaryarray");

  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  while((c = fgetc(IN)) != EOF){
    if(c==' ') continue;
    //    printf("c = %c\n",c);
    if(c=='1') array[num]++;
    num++;
    if(num >= chrlen-1) break;
  }
  fclose(IN);

  return array;
}

TYPE_WIGARRAY *read_mosaics_binfile(char *filename, int arraynum, int binsize){
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  isfile(filename);

  int num;
  TYPE_WIGARRAY *array = (TYPE_WIGARRAY *)my_calloc(arraynum, sizeof(TYPE_WIGARRAY), "mosaicsbinarray");

  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);

    num = atoi(clm[0].str)/binsize;
    if(num < 0 || num > arraynum){
      fprintf(stderr, "\nERROR:%s: invalid arraynum: %d, arraynum:%d\n",__FUNCTION__, num, arraynum);
      exit(0);
    }
    array[num] = VALUE2WIGARRAY(atof(clm[1].str));
  }
  fclose(IN);
  MYFREE(str);
  
  return array;
}

void parse_genometable(char *argv, RefGenome *g){
  int num=1, nummax=100;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];

  isfile(argv);

  g->genome = (Fastadata *)my_calloc(1, sizeof(Fastadata), "g->genome");
  g->genome->name = strdup("Whole genome");
  g->chr = (Fastadata *)my_calloc(nummax, sizeof(Fastadata), "g->chr");

  FILE *IN = my_fopen(argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    g->chr[num].name = strdup(clm[0].str);
    g->chr[num].len = atoi(clm[1].str);
    g->genome->len += g->chr[num].len;
    num++;
    if(num>=nummax){
      nummax += 100;
      g->chr = (Fastadata *)my_realloc(g->chr, nummax*sizeof(Fastadata), "g->chr");
    }
  }
  fclose(IN);
  MYFREE(str);
  
  g->chrnum = num;
  return;
}

void scan_mappability_chrtotal(char *mpfile, RefGenome *g){
  int chr;
  char *filename=NULL;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  FILE *IN;

  if(!mpfile){   // mappabilityÇégÇÌÇ»Ç¢èÍçáÇÕÉQÉmÉÄëSëÃÇ™mappableÇ…
    for(chr=1; chr<g->chrnum; chr++){
      g->genome->len_mpbl += g->chr[chr].len_mpbl = g->chr[chr].len;
      g->chr[chr].p_mpbl = g->chr[chr].len ? g->chr[chr].len_mpbl/(double)g->chr[chr].len: 0;
    }
  }else{
    filename = alloc_str_new(mpfile, 50);
    sprintf(filename, "%s_genome.txt", mpfile);
    IN = my_fopen(filename, FILE_MODE_READ);
    while((fgets(str, STR_LEN, IN))!=NULL){
      if(str[0]=='\n') continue;
      chomp(str);
      ParseLine(str, clm);
      chr = changechr_str2int(clm[0].str, g);
      g->genome->len_mpbl += g->chr[chr].len_mpbl = atoi(clm[1].str);
      g->chr[chr].p_mpbl = g->chr[chr].len ? g->chr[chr].len_mpbl/(double)g->chr[chr].len: 0;
    }
    fclose(IN);
    MYFREE(filename);
  }
  
  g->genome->p_mpbl = g->genome->len ? g->genome->len_mpbl/(double)g->genome->len: 0;

  /* ç≈Ç‡mpblîzóÒí∑ÇÃí∑Ç¢êıêFëÃÇäiî[ */
  int lenmax=0;
  for(chr=1; chr<g->chrnum; chr++){
    if(lenmax < g->chr[chr].len_mpbl){
      lenmax = g->chr[chr].len_mpbl;
      g->chrmax = chr;
    }
  }

  MYFREE(str);
  return;
}

static char *checkchrname(char *name){
  char *p = strstr(name, "chr");
  if(p) return p+3;
  else return name;
}

int changechr_str2int(char *str, RefGenome *g){
  int i, chr=0;
  for(i=1; i<g->chrnum; i++){
    if(!strcmp(checkchrname(str), checkchrname(g->chr[i].name))) chr = i;
  }
  return chr;
}

BedFile *read_bedfile(char *filename, RefGenome *g){
  BedFile *p = bedfile_new(g->chrnum);
  int chr, num, arraynum;
  int nline, nline_rgb;

  Elem clm[ELEM_NUM], clmrgb[3];
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n' || str[0]=='#') continue;
    chomp(str);
    /* browser line */
    if(strstr(str, "browser")) continue;
    /* track line */
    if(strstr(str, "track")) continue;
    /* data line */

    nline = ParseLine(str, clm);
    if(!strcmp(clm[0].str, "chromosome")) continue;
    if(!p->nline) p->nline = nline;
    if(nline<3){
      fprintf(stderr, "error: column number is: %d < 3\n", nline);
      fprintf(stderr, "line is: %s\n", str);
      exit(1);
    }

    /* required line */
    chr = changechr_str2int(clm[0].str, g);
    if(!chr) continue;
    num = p->chr[chr].num;
    arraynum = p->chr[chr].arraynum;
    p->chr[chr].bed[num].s = max(0, atoi(clm[1].str));
    p->chr[chr].bed[num].e = min(atoi(clm[2].str), g->chr[chr].len);
    p->len_total += p->chr[chr].bed[num].e - p->chr[chr].bed[num].s;

    if(p->chr[chr].bed[num].s > p->chr[chr].bed[num].e){
      fprintf(stderr, "Warning: start position > end position in bedfile: %s\n", filename);
    }

    /* optional line */
    p->chr[chr].bed[num].itemRgb.r = -1;
    if(nline >3) p->chr[chr].bed[num].name  = strdup(clm[3].str);
    if(nline >4) p->chr[chr].bed[num].score = atof(clm[4].str);
    if(nline >5){
      if(!strcmp(clm[5].str, "+")) p->chr[chr].bed[num].strand = STRAND_PLUS;
      else p->chr[chr].bed[num].strand = STRAND_MINUS;
    }
    if(nline >6) p->chr[chr].bed[num].ts = atoi(clm[6].str);
    if(nline >7) p->chr[chr].bed[num].te = atoi(clm[7].str);
    if(nline >8){
      nline_rgb = ParseLine_arbit(clm[8].str, clmrgb, ',');
      if(nline_rgb >3){
	fprintf(stderr, "error: itemRgb column is: %s", clm[8].str);
      }
      p->chr[chr].bed[num].itemRgb.r = atoi(clmrgb[0].str);
      p->chr[chr].bed[num].itemRgb.g = atoi(clmrgb[1].str);
      p->chr[chr].bed[num].itemRgb.b = atoi(clmrgb[2].str);
    }
    //    printf("%s %d %d %d\n", p->chr[chr].bed[num].name, p->chr[chr].bed[num].itemRgb.r, p->chr[chr].bed[num].itemRgb.g,p->chr[chr].bed[num].itemRgb.b);

    num++;
    if(num >= arraynum){
      arraynum += BEDFILE_MAX;
      p->chr[chr].bed = (struct bed *)my_realloc(p->chr[chr].bed, sizeof(struct bed)*arraynum, "struct bed");
    }
    p->chr[chr].arraynum = arraynum;
    p->chr[chr].num = num;
  }
  fclose(IN);
  MYFREE(str);

  for(chr=1; chr<g->chrnum; chr++) p->num += p->chr[chr].num;

  return p;
}

void show_bedfile(BedFile *p, int chrnum){
  int i, chr, nline;
  for(chr=1; chr<chrnum; chr++){
    for(i=0; i<p->chr[chr].num; i++){
      printf("chr%d\t%d\t%d\t", chr, p->chr[chr].bed[i].s, p->chr[chr].bed[i].e);
      nline = p->nline - 3;
      if(nline--) printf("%s\t", p->chr[chr].bed[i].name);
      if(nline--) printf("%f\t", p->chr[chr].bed[i].score);
      if(nline--) printf("%s\t", str_strand[p->chr[chr].bed[i].strand]);
      if(nline--) printf("%d\t", p->chr[chr].bed[i].ts);
      if(nline--) printf("%d\t", p->chr[chr].bed[i].te);
      if(nline--) printf("%d,%d,%d\t", p->chr[chr].bed[i].itemRgb.r, p->chr[chr].bed[i].itemRgb.g, p->chr[chr].bed[i].itemRgb.b);
      printf("\n");
    }
  }
  return;
}

char *makearray_inbed(BedChr *p, int arraysize){
  int i,j, beds, bede;
  if(!p->num) return NULL;
  int nbed = p->num;
  char *array = (char *)my_calloc(arraysize, sizeof(char), "array_inbed");
  for(i=0; i<nbed; i++){
    beds = p->bed[i].s;
    bede = p->bed[i].e;
    for(j=beds; j<=bede; j++) array[j] = 1;
  }
  return array;
}

Peak *read_peakfile(char *filename, RefGenome *g){
  int chr;
  int nline;
  Elem clm[ELEM_NUM];
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");

  Peak *p = (Peak *)my_calloc(1, sizeof(Peak), "Peak");
  p->arraynum = READFILE_PEAKNUM;
  p->bs = (struct bs *)my_calloc(p->arraynum, sizeof(struct bs), "peak->bs");

  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n' || str[0]=='#') continue;
    chomp(str);

    nline = ParseLine(str, clm);
    if(!strcmp(clm[0].str, "chromosome")) continue;
    if(nline<3){
      fprintf(stderr, "error: column number is: %d < 3\n", nline);
      fprintf(stderr, "line is: %s\n", str);
      exit(1);
    }

    /* required line */
    chr = changechr_str2int(clm[0].str, g);
    if(!chr) continue;
    p->bs[p->num].chr = chr;
    p->bs[p->num].start = max(0, atoi(clm[1].str));
    p->bs[p->num].end = min(atoi(clm[2].str), g->chr[chr].len);
    if(p->bs[p->num].start > p->bs[p->num].end){
      fprintf(stderr, "Warning: start position > end position in bedfile: %s\n", filename);
    }
    if(clm[3].str){
      if(atoi(clm[3].str)< p->bs[p->num].start || atoi(clm[3].str) > p->bs[p->num].end) p->bs[p->num].maxposi = (p->bs[p->num].start + p->bs[p->num].end)/2;
      else p->bs[p->num].maxposi = atoi(clm[3].str);
    }
    /* ingore rest columns */
    p->num++;
    if(p->num >= p->arraynum){
      p->arraynum += READFILE_PEAKNUM;
      p->bs = (struct bs *)my_realloc(p->bs, sizeof(struct bs)*p->arraynum, "peak->bs");
    }
  }
  fclose(IN);
  MYFREE(str);

  return p;
}

void read_interactionfile(InteractionSet *inter, char *argv, RefGenome *g){
  int num=0, nline;
  Elem clm[ELEM_NUM];
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  inter->arraynum = READFILE_PEAKNUM;
  inter->site = (Interaction *)my_calloc(inter->arraynum, sizeof(Interaction), "inter->site");
  inter->argv = strdup(argv);

  FILE *IN = my_fopen(argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);

    nline = ParseLine(str, clm);
    if(nline<6){
      //      fprintf(stderr, "Warning: column number is: %d < 6\n", nline);
      continue;
    }

    inter->site[num].head.chr   = changechr_str2int(clm[0].str, g);
    inter->site[num].head.start = atoi(clm[1].str);
    inter->site[num].head.end   = atoi(clm[2].str);
    inter->site[num].tail.chr   = changechr_str2int(clm[3].str, g);
    inter->site[num].tail.start = atoi(clm[4].str);
    inter->site[num].tail.end   = atoi(clm[5].str);
    if(!inter->site[num].head.chr || !inter->site[num].tail.chr) continue;
    if(inter->site[num].head.start > inter->site[num].head.end || inter->site[num].tail.start > inter->site[num].tail.end) continue;
    if(nline>=7) inter->site[num].pval = atof(clm[6].str);

    num++;
    if(num >= inter->arraynum){
      inter->arraynum += READFILE_PEAKNUM;
      inter->site = (Interaction *)my_realloc(inter->site, inter->arraynum*sizeof(Interaction), "inter->site");
    }
  }
  fclose(IN);
  MYFREE(str);
  inter->num = num;

  /*  int i;
  for(i=0;i<num;i++) printf("%d %d %d %d %d %d %.2e\n", inter->site[i].head.chr, inter->site[i].head.start, inter->site[i].head.end,inter->site[i].tail.chr, inter->site[i].tail.start, inter->site[i].tail.end, inter->site[i].pval);
  exit(0);*/

  return;
}
