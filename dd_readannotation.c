/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <string.h>
#include "dd_readannotation.h"
#include "alloc.h"
#include "filehandle.h"
#include "readfile.h"
#include "stringp.h"
#include "macro.h"

#define EXON_MAX 2560
#define STRUCT_REPEAT_MAX 1000
#define MEM_MAX_DEFAULT 20

static void read_refFlat(GeneSet *p, RefGenome *g, int chr);
static void read_ENSEMBL(GeneSet *p, RefGenome *g, int chr);
static void read_gtf(GeneSet *p, RefGenome *g, int chr);
static void read_gene_SGD(GeneSet *p, int chr);

void read_graph(Graph *graph, RefGenome *g, int chr, char *name, double min, double max, Graph_Type type){
  int i,start;
  char filename[256];
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  double maxtemp=0;

  graph->name = strdup(name);
  graph->arraynum = g->chr[chr].len/graph->wsize +1;
  graph->array = (double *)my_calloc(graph->arraynum, sizeof(double), "graph_array");

  if(type==GTYPE_SINGLE) sprintf(filename, "%s/%s-bs%d", graph->argv, g->chr[chr].name, graph->wsize);
  else                   sprintf(filename, "%s.%d.%s.xls", graph->argv, graph->wsize, g->chr[chr].name);

  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    start = atoi(clm[0].str);
    if(start % graph->wsize){
      printf("%d %d\n",start, graph->wsize);
      fprintf(stderr,"[E]graph: invalid start position or binsize: <%s>\n", graph->argv); exit(1);
    }
    graph->array[start/graph->wsize] = atof(clm[1].str);
    if(maxtemp < atof(clm[1].str)) maxtemp = atof(clm[1].str);
  }
  fclose(IN);

  if(type==GTYPE_MULTI_PROP){
    double sum=0;
    graph->mmax = 0;
    for(i=0; i<graph->arraynum; i++) sum += graph->array[i];
    for(i=0; i<graph->arraynum; i++){
      graph->array[i] *= 100/sum;
      if(graph->mmax < graph->array[i]) graph->mmax = graph->array[i];
    }
  }else{
    if(max){
      graph->mmin = min;
      graph->mmax = max;
    }else{
      graph->mmin = 0;
      graph->mmax = max(MEM_MAX_DEFAULT, maxtemp);
    }
  }
  MYFREE(str);
  return;
}

void read_gene(GeneSet *p, RefGenome *g, int chr, GeneFile_Type gftype){
  if(!chr){
    fprintf(stderr, "error %s: chr=%d\n", __FUNCTION__, chr);
    return;
  }

  switch(gftype){
  case GFTYPE_REFFLAT: read_refFlat(p, g, chr); break;
  case GFTYPE_ENSEMBL: read_ENSEMBL(p, g, chr); break;
  case GFTYPE_GTF:     read_gtf(p, g, chr);     break;
  case GFTYPE_SGD:     read_gene_SGD(p, chr);   break;
  default: break;
  }
  p->numall += p->num;
  return;
}

static void read_refFlat(GeneSet *p, RefGenome *g, int chr){
  int i,j, num=0;
  int nummax = STRUCT_GENE_MAX;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  Elem *exonclm = (Elem *)my_calloc(EXON_MAX, sizeof(Elem), "exonclm");

  FILE *IN = my_fopen(p->argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);

    if(changechr_str2int(clm[2].str, g) != chr) continue;
    p->gene[num].chr = chr;
    strcpy(p->gene[num].ID, clm[1].str);
    if(strcmp(clm[0].str, "")) strcpy(p->gene[num].name, clm[0].str);
    else strcpy(p->gene[num].name, clm[1].str);
    if(strstr(clm[1].str, "NM")) p->gene[num].genetype = CODING;
    else p->gene[num].genetype = NONCODING;
    if(!strcmp(clm[3].str,"+")) p->gene[num].dir = 1; else p->gene[num].dir = -1;
    p->gene[num].start   = atoi(clm[4].str);
    p->gene[num].end     = atoi(clm[5].str);
    p->gene[num].exonnum = atoi(clm[8].str);
    p->gene[num].exon    = (SE *)my_calloc(p->gene[num].exonnum, sizeof(SE), "exon");
    ParseLine_arbit(clm[9].str, exonclm, ',');
    for(i=0; i<p->gene[num].exonnum; i++) p->gene[num].exon[i].start = atoi(exonclm[i].str);
    ParseLine_arbit(clm[10].str, exonclm, ',');
    for(i=0; i<p->gene[num].exonnum; i++) p->gene[num].exon[i].end = atoi(exonclm[i].str);
    num++;
    if(num>=nummax){
      nummax += STRUCT_GENE_MAX;
      p->gene = (Gene *)my_realloc(p->gene, nummax*sizeof(Gene), "gene");
    }
  }

  char dir;
  for(i=0; i<num; i++){
    dir = p->gene[i].dir;
    for(j=0; j<i; j++){
      if((dir == 1  && p->gene[j].dir == 1  && p->gene[j].start == p->gene[i].start) ||
	 (dir == -1 && p->gene[j].dir == -1 && p->gene[j].end   == p->gene[i].end)){
	p->gene[i].delete = 1;
	continue;
      }
    }
  }
  fclose(IN);

  p->num = num;

  MYFREE(str);
  MYFREE(exonclm);
  return;
}

static void read_ENSEMBL(GeneSet *p, RefGenome *g, int chr){
  int i, num=0;
  int nummax = STRUCT_GENE_MAX;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM], temp[2];
  Elem *exonclm = (Elem *)my_calloc(EXON_MAX, sizeof(Elem), "exonclm");

  FILE *IN = my_fopen(p->argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    if(changechr_str2int(clm[5].str, g) != chr) continue;
    p->gene[num].chr = chr;
    strcpy(p->gene[num].name, clm[0].str);

    if(!strcmp(clm[1].str, "protein_coding"))            p->gene[num].genetype = CODING;
    else if(!strcmp(clm[1].str, "processed_transcript")) p->gene[num].genetype = PROCESS;
    else if(!strcmp(clm[1].str, "pseudogene"))           p->gene[num].genetype = PSEUDO;
    else if(!strcmp(clm[1].str, "miRNA"))                p->gene[num].genetype = MIRNA;
    else if(strstr(clm[1].str, "RNA"))                   p->gene[num].genetype = NONCODING;
    else                                                 p->gene[num].genetype = OTHERS;
    p->gene[num].dir = atoi(clm[2].str);
    p->gene[num].start  = atoi(clm[3].str);
    p->gene[num].end    = atoi(clm[4].str);
    p->gene[num].exonnum = atoi(clm[7].str);
    p->gene[num].exon = (SE *)my_calloc(p->gene[num].exonnum, sizeof(SE), "exon");
    ParseLine_arbit(clm[8].str, exonclm, ',');
    for(i=0; i<p->gene[num].exonnum; i++){
      ParseLine_arbit(exonclm[i].str, temp, '-');
      p->gene[num].exon[i].start = atoi(temp[0].str);
      p->gene[num].exon[i].end   = atoi(temp[1].str);
    }
    strcpy(p->gene[num].ID, clm[9].str);
    num++;
    if(num>=nummax){
      nummax += STRUCT_GENE_MAX;
      p->gene = (Gene *)my_realloc(p->gene, nummax*sizeof(Gene), "gene");
    }
  }
  fclose(IN);

  p->num = num;

  MYFREE(str);
  MYFREE(exonclm);
  return;
}

static void read_gtf(GeneSet *p, RefGenome *g, int chr){
  int num=0;
  int nummax=STRUCT_GENE_MAX;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  char *p1=NULL, *p2=NULL;

  FILE *IN = my_fopen(p->argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    if(changechr_str2int(clm[0].str, g) != chr) continue;
    p->gene[num].chr = chr;
    if(strcmp(clm[2].str, "exon")) continue;
    if(!strcmp(clm[1].str, "protein_coding"))  p->gene[num].genetype = CODING;
    else if(!strcmp(clm[1].str, "pseudogene")) p->gene[num].genetype = PSEUDO;
    else if(!strcmp(clm[1].str, "rRNA"))       p->gene[num].genetype = rRNA;
    else                                       p->gene[num].genetype = NONCODING;
    p->gene[num].start = atoi(clm[3].str);
    p->gene[num].end   = atoi(clm[4].str);
    if(!strcmp(clm[6].str, "+")) p->gene[num].dir = 1; else p->gene[num].dir = -1;
    p1 = strstr(clm[8].str, "gene_name");
    p2 = strstr(p1+11, "\"");
    strncpy(p->gene[num].name, p1+11, p2-p1-11);
    num++;
    if(p->num>=nummax){
      nummax += STRUCT_GENE_MAX;
      p->gene = (Gene *)my_realloc(p->gene, nummax*sizeof(Gene), "gene");
    }
  }
  fclose(IN);

  p->num = num;

  MYFREE(str);
  return;
}

static void read_gene_SGD(GeneSet *p, int chr){
  int num=0, chrtemp=0, nummax=STRUCT_GENE_MAX;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];

  FILE *IN = my_fopen(p->argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    if(!strcmp(clm[8].str, "2-micron")) chrtemp=18;
    else if(atoi(clm[8].str)) chrtemp = atoi(clm[8].str);
    if(chrtemp != chr) continue;
    p->gene[num].chr = chr;

    if(!strcmp(clm[1].str, "ARS")) continue;
    else if(!strcmp(clm[1].str, "centromere")){
      sprintf(p->gene[num].name, "CEN_chr%d", chr);
      p->gene[num].genetype = CENTROMERE;
    }else if(!strcmp(clm[1].str, "ORF")){ 
      if(clm[4].str[0]!='\0') strcpy(p->gene[num].name, clm[4].str);
      else strcpy(p->gene[num].name, clm[3].str);
      p->gene[num].genetype = CODING;
    }else if(!strcmp(clm[1].str, "tRNA")){
      if(clm[4].str[0]!='\0') strcpy(p->gene[num].name, clm[4].str);
      else strcpy(p->gene[num].name, clm[3].str);
      p->gene[num].genetype = NONCODING;
    }else if(!strcmp(clm[1].str, "rRNA")){
      strcpy(p->gene[num].name, "rRNA");
      p->gene[num].genetype = rRNA;
    }else continue;

    if(!strcmp(clm[11].str, "C")){
      p->gene[num].start = atoi(clm[10].str);
      p->gene[num].end   = atoi(clm[9].str);
      p->gene[num].dir = -1;
    }else{
      p->gene[num].start = atoi(clm[9].str);
      p->gene[num].end   = atoi(clm[10].str);
      p->gene[num].dir = 1;
    }
    if(!strcmp(clm[1].str,"ARS") || !strcmp(clm[1].str,"centromere")) p->gene[num].dir = 0;
    num++;
    if(num>=nummax){
      nummax += STRUCT_GENE_MAX;
      p->gene = (Gene *)my_realloc(p->gene, nummax*sizeof(Gene), "gene");
    }
  }
  fclose(IN);

  p->num = num;

  MYFREE(str);
  return;
}

void read_ARS_OriDB(char *argv, GeneSet *p, int chr){
  int nummax = STRUCT_GENE_MAX;
  printf("reading ARS file...\n");
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  
  FILE *IN = my_fopen(argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    if(atoi(clm[4].str) != chr) continue;
    p->gene[p->num].chr = chr;
    if(strstr(clm[2].str,"ARS")) strcpy(p->gene[p->num].name, clm[2].str);
    else sprintf(p->gene[p->num].name, "ARS_%s", clm[2].str);
    p->gene[p->num].start = atoi(clm[5].str);
    p->gene[p->num].end   = atoi(clm[6].str);
    p->gene[p->num].dir   = 0;
    p->gene[p->num].genetype = ARS;
    p->num++;
    if(p->num>=nummax){
      nummax += STRUCT_GENE_MAX;
      p->gene = (Gene *)my_realloc(p->gene, nummax*sizeof(Gene), "gene");
    }
  }
  fclose(IN);

  MYFREE(str);
  return;
}

void read_TER_scer(char *argv, GeneSet *p, int chr){
  int nummax=STRUCT_GENE_MAX;
  printf("reading TER file...\n");
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  
  FILE *IN = my_fopen(argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    if(atoi(clm[1].str) != chr) continue;
    p->gene[p->num].chr = chr;
    strcpy(p->gene[p->num].name, clm[0].str);
    p->gene[p->num].start = atoi(clm[2].str);
    p->gene[p->num].end   = atoi(clm[3].str);
    p->gene[p->num].dir   = 0;
    p->gene[p->num].genetype = TER;
    p->num++;
    if(p->num>=nummax){
      nummax += STRUCT_GENE_MAX;
      p->gene = (Gene *)my_realloc(p->gene, nummax*sizeof(Gene), "gene");
    }
  }
  fclose(IN);
  MYFREE(str);
  return;
}

void read_repeat(DDParam *d, RefGenome *g, int chr){
  int on=0, num=0, nummax = STRUCT_REPEAT_MAX;
  RepeatType type;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  d->repeat.repeat = (Repeat *)my_calloc(nummax, sizeof(Repeat), "repeat");

  printf("reading repeat file...\n");
  FILE *IN = my_fopen(d->repeat.argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    ParseLine(str, clm);
    if(changechr_str2int(clm[0].str, g)!=chr){
      if(!on) continue; 
      else break;
    }
    on=1;
    d->repeat.repeat[num].start = atoi(clm[1].str);
    d->repeat.repeat[num].end   = atoi(clm[2].str);
    if(clm[3].str[0]=='+') d->repeat.repeat[num].dir = 1;
    else d->repeat.repeat[num].dir = -1;
    if(!strcmp(clm[5].str, "SINE"))               type = RM_SINE;
    else if(strstr(clm[5].str, "LINE"))           type = RM_LINE;
    else if(strstr(clm[5].str, "LTR"))            type = RM_LTR;
    else if(strstr(clm[5].str, "DNA"))            type = RM_DNA;
    else if(strstr(clm[5].str, "Simple_repeat"))  type = RM_Simple;
    else if(strstr(clm[5].str, "Low_complexity")) type = RM_Low_Com;
    else if(strstr(clm[5].str, "Satellite"))      type = RM_Satellite;
    else if(strstr(clm[5].str, "RNA"))            type = RM_RNA;
    else type = RM_Other;
    d->repeat.repeat[num].rtype = type;
    strcpy(d->repeat.repeat[num].name, clm[4].str);
    strcpy(d->repeat.repeat[num].class, clm[5].str);
    num++;
    if(num >= nummax){
      nummax += STRUCT_REPEAT_MAX;
      d->repeat.repeat = (Repeat *)my_realloc(d->repeat.repeat, nummax*sizeof(Repeat), "repeat");
    }
  }
  fclose(IN);

  d->repeat.num = num;

  MYFREE(str);
  return;
}

void read_genefile(DDParam *d, BedFile *p, RefGenome *g, int chr){
  int i, s,e,num=0, on=0;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  GeneSet *gene=&(d->gene);
  
  FILE *IN = my_fopen(d->genefile_argv, FILE_MODE_READ);
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n' || str[0]=='#') continue;
    chomp(str);
    
    on=0;
    s=0; e=0;
    for(i=0; i<gene->num; i++){
      if(!strcmp(gene->gene[i].name, str) && !gene->gene[i].delete){
	printf("str %s, name=%s, %d, %d, %d\n",str, gene->gene[i].name, gene->gene[i].dir, gene->gene[i].start, gene->gene[i].end);
	if(!s && !e){
	  s = gene->gene[i].start;
	  e = gene->gene[i].end;
	  on++;
	}else{
	  s = min(s, gene->gene[i].start);
	  e = max(e, gene->gene[i].end);
	}
      }
      LOG("\n\n\nnum=%d, %d, %d\n",num, s, e);
    }
    if(on){
      p->chr[chr].bed[num].s = max(0, s - d->genefile_len);
      p->chr[chr].bed[num].e = min(e + d->genefile_len, g->chr[chr].len -1);
      num++;
    }
  }
  p->chr[chr].num = num;
  //  printf("a%d, b%d\n", p->chr[chr].num, num);
  fclose(IN);
  MYFREE(str);

  return;
}
