/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <string.h>
#include <zlib.h>
#include "drompa_param_new.h"
#include "drompa_readfile.h"
#include "alloc.h"
#include "stringp.h"
#include "filehandle.h"
#include "readfile.h"
#include "macro.h"

#define ARRAYNUM 50

static void read_parse2wig_stats(Samplefile *p, RefGenome *g, int binsize);

/* 1:ChIP   2:Input   3:name   4:peaklist   5:binsize
   6:scale_tag   7:scale_ratio   8:scale_pvalue */
SamplePair *scan_samplestr(DrParam *p, char **str, int chrnum){
  int i, nline;
  SamplePair *sample = SamplePair_new(p->samplenum, chrnum);
  for(i=0; i<p->samplenum; i++){
    Elem *clm = (Elem *)my_calloc(ELEM_NUM, sizeof(Elem), "clm");
    nline = ParseLine_arbit(str[i], clm, ',');
    if(nline >8){
      fprintf(stderr, "error: sample %d has ',' more than 8: %s\n", i+1, str[i]);
      exit(1);
    }
    LOG("%d ChIP:%s Input:%s name:%s peak:%s\n", i+1, clm[0].str, clm[1].str, clm[2].str, clm[3].str);
    if(!strcmp(clm[0].str, "")){
      fprintf(stderr, "please specify ChIP sample %d: %s.\n", i+1, str[i]);
      exit(1);
    }else                      sample[i].ChIP->argv  = strdup(clm[0].str);
    if(strcmp(clm[1].str, "")) sample[i].Input->argv = strdup(clm[1].str);
    if(strcmp(clm[2].str, "")) sample[i].linename    = strdup(clm[2].str);
    else                       sample[i].linename    = sample[i].ChIP->argv;
    if(strcmp(clm[3].str, "")) sample[i].peak_argv   = strdup(clm[3].str);
    if(strcmp(clm[4].str, "")) sample[i].binsize     = atoi(clm[4].str);
    if(strcmp(clm[5].str, "")) sample[i].scale_tag   = atof(clm[5].str);
    if(strcmp(clm[6].str, "")) sample[i].scale_ratio = atof(clm[6].str);
    if(strcmp(clm[7].str, "")) sample[i].scale_pvalue= atof(clm[7].str);
    LOG("sample %d ChIP:%s Input:%s name:%s peak:%s binsize:%d scale_tag %f scale_ratio %f scale_pvalue %f\n", i, sample[i].ChIP->argv, sample[i].Input->argv, sample[i].linename, sample[i].peak_argv, sample[i].binsize, sample[i].scale_tag, sample[i].scale_ratio, sample[i].scale_pvalue);
    MYFREE(clm);
  }
  return sample;
}

static double definestep(DrParam *p, SamplePair *sample, RefGenome *g){
  double ave=0;
  int i, num=0;
  int chr = g->chrmax;
  int binnum = sample->binnum[chr];
  printf("chr=%d\n",chr);

  dr_read_wigdata(p, sample, sample->ChIP, g, chr);
  dr_read_wigdata(p, sample, sample->Input, g, chr);

  for(i=0; i<binnum; i++){
    if(!sample->ChIP->data[i] || !sample->Input->data[i]) continue;
    ave += sample->ChIP->data[i] + sample->Input->data[i];
    num++;
  }
  ave = ave/(double)num;

  int s, cnt=1;
  double step=0;
  while(1){
    s = ave*2*cnt/ARRAYNUM;
    if(!s) cnt *= 10;
    else{ step = s/(double)cnt; break;}
  }
  //  printf("%f, %f\n", ave, step);

  dr_delete_wigdata(sample->ChIP);
  dr_delete_wigdata(sample->Input);

  return step;
}

static void calc_NCIS(DrParam *p, SamplePair *sample, RefGenome *g){
  int i, chr, binnum;
  double step;
  double ratioarray[ARRAYNUM], IParray[ARRAYNUM], Inputarray[ARRAYNUM];
  int num, n[ARRAYNUM];
  for(i=0;i<ARRAYNUM;i++){
    ratioarray[i]=0;
    IParray[i]=0;
    Inputarray[i]=0;
    n[i]=0;
  }

  step = definestep(p, sample, g);

  for(chr=1; chr<g->chrnum; chr++){
    binnum = sample->binnum[chr];

    dr_read_wigdata(p, sample, sample->ChIP, g, chr);
    dr_read_wigdata(p, sample, sample->Input, g, chr);

    for(i=0; i<binnum; i++){
      if(!sample->ChIP->data[i] || !sample->Input->data[i]) continue;
      num = min((sample->ChIP->data[i] + sample->Input->data[i])/step, ARRAYNUM-1);
      ratioarray[num] += sample->ChIP->data[i]/(double)sample->Input->data[i];
      IParray[num]    += sample->ChIP->data[i];
      Inputarray[num] += sample->Input->data[i];
      n[num]++;
      //      printf("plot %s\t%f\t%f\n", g->chr[chr].name, WIGARRAY2VALUE(sample->ChIP->data[i]), WIGARRAY2VALUE(sample->Input->data[i]));
    }
    dr_delete_wigdata(sample->ChIP);
    dr_delete_wigdata(sample->Input);
  }

  int binthre=0;
  for(i=0;i<ARRAYNUM;i++) binthre += n[i];
  binthre *= 0.75; // NCIS threshold
  int bintemp=0;

  double r=0, r_pre=0;
  int tw=0, bwt=0;
  for(i=0; i<ARRAYNUM; i++){
    bintemp += n[i];
    if(!bwt && bintemp>binthre) bwt = i;
    r = n[i] ? ratioarray[i]/n[i]: 0;
    if(bintemp > binthre && r_pre && !tw && r > r_pre) tw=i;
    LOG("%d\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%d\n",i, i*step, r, tw, IParray[i]/n[i], Inputarray[i]/n[i], n[i]);
    r_pre = r;
  }
  if(!tw) tw = bwt;
  double IPsum=0, Inputsum=0;
  for(i=0; i<tw; i++){
    IPsum += IParray[i];
    Inputsum += Inputarray[i];
  }

  sample->comp->genome->ratio = Inputsum ? IPsum/Inputsum: 1;
  LOG("NCIS: tw %d, bwt %d, IPsum %f, Inputsum %f, ratio %f\n", tw, bwt, IPsum, Inputsum, sample->comp->genome->ratio);
  return;
}

void dr_calc_global_param(DrParam *p, SamplePair *sample, RefGenome *g){
  int i, chr;
  /* read_parse2wig_stats */
  for(i=0; i<p->samplenum; i++){
    if(sample[i].copyC==-1) read_parse2wig_stats(sample[i].ChIP, g, sample[i].binsize);
    if(!sample[i].Input->argv) continue;

    if(sample[i].copyI==-1) read_parse2wig_stats(sample[i].Input, g, sample[i].binsize);   
    if(sample[i].copycomp==-1){
      if(p->ntype==TYPE_RATIO_NONE){
	sample[i].comp->genome->ratio = 1;
	for(chr=1; chr<g->chrnum; chr++) sample[i].comp->chr[chr].ratio = 1;
      }else if(p->ntype==TYPE_RATIO_TOTALREAD){
	sample[i].comp->genome->ratio = sample[i].Input->genome->nread ? sample[i].ChIP->genome->nread/(double)sample[i].Input->genome->nread: 0;
	for(chr=1; chr<g->chrnum; chr++){
	  sample[i].comp->chr[chr].ratio = sample[i].Input->chr[chr].nread ? sample[i].ChIP->chr[chr].nread/(double)sample[i].Input->chr[chr].nread: 0;
	  LOG("sample%d: %s: ChIP read=%ld, Input read=%ld, ChIP/Input = %.2f\n", i+1, g->chr[chr].name, sample[i].ChIP->chr[chr].nread, sample[i].Input->chr[chr].nread, sample[i].comp->chr[chr].ratio);
	}
      }else if(p->ntype==TYPE_RATIO_NCIS){
	calc_NCIS(p, &(sample[i]), g);
      }
      printf("sample%d: genome: ChIP read=%ld, Input read=%ld, ChIP/Input = %.2f\n", i+1, sample[i].ChIP->genome->nread, sample[i].Input->genome->nread, sample[i].comp->genome->ratio);
    }
  }
  return;
}

static void read_parse2wig_stats(Samplefile *p, RefGenome *g, int binsize){

  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  char *filename = alloc_str_new(p->argv, 50);
  sprintf(filename, "%s.%d.xls", p->argv, binsize);
  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  char *tp;

  int chr;
  int on=0;
  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    if(!on){
      if(strstr(str, "Poisson: λ=")){
	//	printf("%s\n", str+12);
	p->lambda = atof(str+12);
      }else if(strstr(str, "Negative binomial:")){
	//	printf("%f\n", atof(str+21));
	p->nb_p = atof(str+21);
	//printf("%f\n", atof(tp+2));
	tp = strstr(str, "n=");
	p->nb_n = atof(tp+2);
	tp = strstr(str, "p0=");
	p->nb_p0 = atof(tp+3);
      }else if(strstr(str, "Whole genome")){
	chomp(str);
	ParseLine(str, clm);
	//	p->genome->nread = atoi(delimit_str(clm[8].str, ' '));
	p->genome->nread = atoi(clm[16].str);
	on=1;
      }
    }else{
      chomp(str);
      ParseLine(str, clm);
      chr = changechr_str2int(clm[0].str, g);
      if(!chr){
	fprintf(stderr, "Warning:%s: %s does not appear in genome_table.\n", __FUNCTION__, clm[0].str);
	continue;
      }
      p->chr[chr].nread = atoi(delimit_str(clm[8].str, ' '));
    }
  }
  fclose(IN);

  //  printf("%f %f %f %f\n",p->lambda, p->nb_p, p->nb_n, p->nb_p0);

  MYFREE(filename);
  MYFREE(str);
  return;
}

static void read_wig(TYPE_WIGARRAY **data, void *IN, int binsize, PWfile_Type itype){
  char *status=NULL;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  int start, on=0;
  char *tp=NULL;
  
  while(1){
    if(itype==TYPE_COMPRESSWIG) status = gzgets(IN, str, STR_LEN);
    else status = fgets(str, STR_LEN, IN);
    if(!status) break;
    if(str[0]=='\n') continue;
    if(!on && !strstr(str, "variableStep")) continue;
    else if(strstr(str, "variableStep")){
      tp = strtok(str, "\t");   /* スペースとタブを両方認識できるようにstrtokを使う */
      tp = strtok(NULL, "\t");
      tp = strtok(NULL, "\t");
      if(atoi(tp+5) != binsize){    // span of wigfile != binsize
	fprintf(stderr,"ERROR:%s: invalid binsize : %d.\n",__FUNCTION__, atoi(tp+5)); exit(1);
      }
      on=1;
      continue;
    }
    chomp(str);
    ParseLine(str, clm);
    start = atoi(clm[0].str) -1;   // wigファイルは1から始まる
    if(start%binsize){
      fprintf(stderr,"ERROR:%s: invalid start position: %d / %d\n",__FUNCTION__, start, binsize); exit(1);
    }
    (*data)[start/binsize] = VALUE2WIGARRAY(atof(clm[1].str));
  }
  MYFREE(str);
  return;
}

static void read_wiggz(TYPE_WIGARRAY **data, char *prefix, int binsize, PWfile_Type itype){
  char *filename = alloc_str_new(prefix, 50);
  sprintf(filename, "%s.%d.wig.gz", prefix, binsize);

  struct gzFile_s *gzIN=NULL;
  if((gzIN = gzopen(filename, "r"))==NULL){
    fprintf(stderr,"[E] Cannot open <%s>.\n", filename); 
    exit(1);
  }
  read_wig(data, gzIN, binsize, itype);
  gzclose(gzIN);
  MYFREE(filename);
  return;
}

static void read_rawwig(TYPE_WIGARRAY **data, char *prefix, int binsize, PWfile_Type itype){
  char *filename = alloc_str_new(prefix, 50);
  sprintf(filename, "%s.%d.wig", prefix, binsize);
  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  read_wig(data, IN, binsize, itype);
  fclose(IN);
  MYFREE(filename);
  return;
}

static void read_binary(TYPE_WIGARRAY **data, char *prefix, int binsize, int datasize){
  char *filename = alloc_str_new(prefix, 50);
  sprintf(filename, "%s.%d.bin", prefix, binsize);
  
  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  if(fread((*data), sizeof(TYPE_WIGARRAY)*datasize, 1, IN)!=1){
    if(ferror(IN)!=0) printf("ferror\n");
    if(feof(IN)!=0) printf("feof\n");
    fprintf(stderr,"[E]fread error:%s\n", filename);
    exit(1);
  }
  fclose(IN);

  MYFREE(filename);
  return;
}

static void read_bedGraph(TYPE_WIGARRAY **data, char *prefix, int chrref, int binsize, RefGenome *g){
  char *filename = alloc_str_new(prefix, 50);
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];

  sprintf(filename, "%s.%d.bedGraph", prefix, binsize);
  int chr, on=0, start;
  FILE *IN = my_fopen(filename, FILE_MODE_READ);
  while(fgets(str, STR_LEN, IN)){
    if(str[0]=='\n') continue;
    if(!on && !strstr(str, "track")) continue;
    else if(strstr(str, "track")){
      on=1;
      continue;
    }
    chomp(str);
    ParseLine_arbit(str, clm, ' ');
    chr = changechr_str2int(clm[0].str, g);
    if(!chr){
      fprintf(stderr,"ERROR:%s: invalid chr: %s\n",__FUNCTION__, clm[0].str); exit(1);
    }
    if(chr != chrref) continue;
    start = atoi(clm[1].str);
    if(start%binsize){
      fprintf(stderr,"ERROR:%s: invalid start position: %d / %d\n",__FUNCTION__, start, binsize); exit(1);
    }
    (*data)[start/binsize] = VALUE2WIGARRAY(atof(clm[3].str));
  }
  
  fclose(IN);
  MYFREE(str);
  MYFREE(filename);
  return;
}

static void smooth_tags(TYPE_WIGARRAY **data, int smoothing, int binsize, int datasize){
  int i,j;
  int smooth_numsum = smoothing/binsize;
  if(smooth_numsum <=1) return;
  if(!(smooth_numsum%2)) smooth_numsum++; // odd num
  int smooth_num = smooth_numsum/2;
  int ave, temp[smooth_num];

  for(i=0; i<smooth_num; i++) temp[i] = (*data)[i];
  for(i=smooth_num; i<datasize-smooth_num; i++){
    ave=0;
    for(j=0; j<smooth_num; j++) ave += temp[j];
    for(j=i; j<=i+smooth_num; j++) ave += (*data)[j];
    for(j=1; j<smooth_num; j++) temp[j-1] = temp[j];
    temp[smooth_num-1] = (*data)[i];
    (*data)[i] = ave/(double)smooth_numsum;
  }
  return;
}

void dr_read_wigdata(DrParam *p, SamplePair *sample, Samplefile *s, RefGenome *g, int chr){
  int binsize = sample->binsize;
  int binnum = sample->binnum[chr];
  s->data = (TYPE_WIGARRAY *)my_calloc(binnum, sizeof(TYPE_WIGARRAY), "sample->data");

  char *prefix = alloc_str_new(s->argv, strlen(g->chr[chr].name) +10);
  sprintf(prefix, "%s_%s", s->argv, g->chr[chr].name);
  if(p->itype==TYPE_BINARY)             read_binary(&(s->data), prefix, binsize, binnum);
  else if(p->itype==TYPE_COMPRESSWIG)   read_wiggz(&(s->data), prefix, binsize, p->itype);
  else if(p->itype==TYPE_UNCOMPRESSWIG) read_rawwig(&(s->data), prefix, binsize, p->itype);
  else if(p->itype==TYPE_BEDGRAPH)      read_bedGraph(&(s->data), s->argv, chr, binsize, g);

  /* smoothing*/
  if(p->smoothing) smooth_tags(&(s->data), p->smoothing, binsize, binnum);

  MYFREE(prefix);
  return;
}

void dr_delete_wigdata(Samplefile *s){
  MYFREE(s->data);
  return;
}
