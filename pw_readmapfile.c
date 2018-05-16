/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <string.h>
#include <zlib.h>
#include "pw_readmapfile.h"
#include "pw_param_new.h"
#include "alloc.h"
#include "readfile.h"
#include "stringp.h"
#include "filehandle.h"
#include "macro.h"

typedef struct{
  char *name;
  int chr_F3, chr_F5;
  int F3, F5;
  Strand strand;
  int fraglen;
  int readlen_F3, readlen_F5;
  short num_multimapped;
} FragmentData;

static void do_parse(PwParam *p, char *inputfile, Mapfile *mapfile, RefGenome *);
static void parse_sam(PwParam *p, Mapfile *mapfile, RefGenome *g, char *inputfile);
static void parse_bowtie(PwParam *p, Mapfile *mapfile, RefGenome *g, char *inputfile);
static void parse_tagAlign(PwParam *p, Mapfile *mapfile, RefGenome *g, char *inputfile);

static void add_SeqStats_to_genome(Mapfile *mapfile, RefGenome *g){
  int chr;
  Strand strand;
  for(chr=1; chr<g->chrnum; chr++){
    for(strand=0;strand<STRANDNUM;strand++){
      mapfile->genome->seq[strand].n_read_infile += mapfile->chr[chr].seq[strand].n_read_infile;
      mapfile->genome->seq[strand].n_readname    += mapfile->chr[chr].seq[strand].n_readname;
      mapfile->chr[chr].both.n_read_infile       += mapfile->chr[chr].seq[strand].n_read_infile;
      mapfile->chr[chr].both.n_readname          += mapfile->chr[chr].seq[strand].n_readname;
      LOG("%s: chr%d%s: n_read_infile: %ld\tn_readname: %.1Lf\n",__FUNCTION__,chr, str_strand[strand], mapfile->chr[chr].seq[strand].n_read_infile, mapfile->chr[chr].seq[strand].n_readname);
    }
    mapfile->genome->both.n_read_infile += mapfile->chr[chr].both.n_read_infile;
    mapfile->genome->both.n_readname    += mapfile->chr[chr].both.n_readname;
    LOG("%s: chr%dboth: n_read_infile: %ld\tn_readname: %.1Lf\n",__FUNCTION__,chr, mapfile->chr[chr].both.n_read_infile, mapfile->chr[chr].both.n_readname);
  }
  for(strand=0; strand<STRANDNUM; strand++){
    LOG("%s: genome%s: n_read_infile: %ld\tn_readname: %.1Lf\n",__FUNCTION__, str_strand[strand], mapfile->genome->seq[strand].n_read_infile, mapfile->genome->seq[strand].n_readname);
  }
  LOG("%s: genomeboth: n_read_infile: %ld\tn_readname: %.1Lf\n",__FUNCTION__, mapfile->genome->both.n_read_infile, mapfile->genome->both.n_readname);
  return;
}

static void output_read_fragment_distribution(PwParam *p, Mapfile *mapfile){
  FILE *OUT = NULL;
  int i;
  char *outputfile = alloc_str_new(p->output_dir, strlen(p->output_prefix) +100);
  
  /* read length distribution */
  sprintf(outputfile, "%s/%s.readlength_dist.xls", p->output_dir, p->output_prefix);
  OUT = my_fopen(outputfile, FILE_MODE_WRITE);
  fprintf(OUT, "F3 distritbution\n");
  fprintf(OUT, "length\tread number\tproportion\n");
  for(i=0; i<DIST_READLEN_MAX; i++){
    if(mapfile->fstats.dist_readlen_F3[i]) fprintf(OUT, "%d\t%d\t%.3f\n", i, mapfile->fstats.dist_readlen_F3[i], mapfile->fstats.dist_readlen_F3[i]/(double)mapfile->genome->both.n_read_infile);
  }
  if(p->rtype==READTYPE_PAIR){
    fprintf(OUT, "\n\nF5 distritbution\n");
    fprintf(OUT, "length\tread number\tproportion\n");
    for(i=0; i<DIST_READLEN_MAX; i++){
      if(mapfile->fstats.dist_readlen_F5[i]) fprintf(OUT, "%d\t%d\t%.3f\n", i, mapfile->fstats.dist_readlen_F5[i], mapfile->fstats.dist_readlen_F5[i]/(double)mapfile->genome->both.n_read_infile);
    }
  }
  fclose(OUT);

  /* fragment length distribution */
  int flen_max=0;
  if(p->rtype==READTYPE_PAIR){
    sprintf(outputfile, "%s/%s.fragmentlength_dist.xls", p->output_dir, p->output_prefix);
    OUT = my_fopen(outputfile, FILE_MODE_WRITE);
    fprintf(OUT, "length\tread number\tproportion\n");
    for(i=0; i<DIST_FRAGLEN_MAX; i++){
      if(mapfile->fstats.dist_fraglen[i]) fprintf(OUT, "%d\t%d\t%.3f\n", i, mapfile->fstats.dist_fraglen[i], mapfile->fstats.dist_fraglen[i]/(double)mapfile->genome->both.n_read_infile);
      if(flen_max < mapfile->fstats.dist_fraglen[i]){
	flen_max = mapfile->fstats.dist_fraglen[i];
	p->fraglen = i;
      }
    }
    if(mapfile->fstats.dist_fraglen[i]) fprintf(OUT, ">%d\t%d\t%.3f\n", DIST_FRAGLEN_MAX, mapfile->fstats.dist_fraglen[i], mapfile->fstats.dist_fraglen[i]/(double)mapfile->genome->both.n_read_infile);
    fprintf(OUT, "estimated fragment length: %d\n", p->fraglen);
    fclose(OUT);
  }

  MYFREE(outputfile);
  return;
}

Mapfile *read_mapfile(PwParam *p, RefGenome *g){
  Mapfile *mapfile = mapfile_new(g->chrnum, p);
  char *inputfile = strdup(p->inputfile);
  char *tp = strtok(inputfile, "," );
  do_parse(p, tp, mapfile, g);
  while(tp){
    tp = strtok(NULL, ",");
    if(tp) do_parse(p, tp, mapfile, g);
  }
  MYFREE(inputfile);

  add_SeqStats_to_genome(mapfile, g);
  output_read_fragment_distribution(p, mapfile);  /* output distributions of read length and fragment length */

  return mapfile;
}

static void do_parse(PwParam *p, char *inputfile, Mapfile *mapfile, RefGenome *g){
  isfile(inputfile);
  printf("Parsing %s...\n",inputfile);
  switch(p->ftype){
  case FILETYPE_BAM:    /* same as SAM */
  case FILETYPE_SAM:    parse_sam(p, mapfile, g, inputfile);    break;
  case FILETYPE_BOWTIE: parse_bowtie(p, mapfile, g, inputfile); break;
  case FILETYPE_TAGALIGN: parse_tagAlign(p, mapfile, g, inputfile); break;
  }
  printf("done.\n");
}

static void init_frag(FragmentData *p){
  MYFREE(p->name);
  p->chr_F3 = p->chr_F5 = 0;
  p->F3 = p->F5 = -1;
  p->strand  = -1;
  p->fraglen = -1;
  p->readlen_F3 = p->readlen_F5 = -1;
  p->num_multimapped = 1;
}

static int define_readstartposition(int posi, int readlen, Strand strand, Inputfiletype ftype){
  int s=0;
  if(strand==STRAND_PLUS) s = posi; else s = posi + readlen;
  if(ftype==FILETYPE_SAM || ftype==FILETYPE_BAM) s--; /* SAM format uses 1-based position */
  return s;
}

/* frag構造体からMapfileにデータを変換 
 * readlenを考慮したstart, endをMapfileに格納 
 * SAMfileのposition-1もここで行う */
static void add_fragment_to_readarray(PwParam *pwparam, Mapfile *p, FragmentData *p_frag){
  int chr = p_frag->chr_F3;
  Strand strand = p_frag->strand, strandF5=0;
  int num = p->chr[chr].seq[strand].n_read_infile;
  int narray = p->readarray[chr][strand].narray;
  
  //  printf("%s\tchr%d, strand:%s, F3len:%d, F5len:%d, fraglen:%d, n_read_infile: %ld\tn_readname: %.1Lf\n",p_frag->name, chr, str_strand[strand], p_frag->readlen_F3, p_frag->readlen_F5, p_frag->fraglen, p->chr[chr].seq[strand].n_read_infile, p->chr[chr].seq[strand].n_readname);

  p->readarray[chr][strand].F3[num] = define_readstartposition(p_frag->F3, p_frag->readlen_F3, strand, pwparam->ftype);
  if(pwparam->rtype==READTYPE_PAIR){
    if(strand==STRAND_PLUS) strandF5 = STRAND_MINUS; else strandF5 = STRAND_PLUS;
    p->readarray[chr][strand].F5[num] = define_readstartposition(p_frag->F5, p_frag->readlen_F5, strandF5, pwparam->ftype);
  }
  p->readarray[chr][strand].weight[num] = WEIGHT2INT(1/(double)p_frag->num_multimapped);
  //  printf("F3:%d, num_multimapped:%d\n", p->readarray[chr][strand].F3[num], p->readarray[chr][strand].weight[num]);
  //  printf("F3:%d, F5:%d, num_multimapped:%d, strand%s\n", p->readarray[chr][strand].F3[num], p->readarray[chr][strand].F5[num], p->readarray[chr][strand].weight[num]);
  num++;
  p->chr[chr].seq[strand].n_read_infile = num;
  p->chr[chr].seq[strand].n_readname += 1/(double)p_frag->num_multimapped;
  
  if(num >= narray){
    narray += READARRAY_NUM;
    p->readarray[chr][strand].F3 = (int *)my_realloc(p->readarray[chr][strand].F3, sizeof(int)*(narray), "mapfile->readarray[chr].F3");
    if(pwparam->rtype==READTYPE_PAIR) p->readarray[chr][strand].F5 = (int *)my_realloc(p->readarray[chr][strand].F5, sizeof(int)*(narray), "mapfile->readarray[chr].F5");
    p->readarray[chr][strand].weight = (int *)my_realloc(p->readarray[chr][strand].weight, sizeof(int)*(narray), "mapfile->readarray[chr].weight");
    p->readarray[chr][strand].delete = (bool *)my_realloc(p->readarray[chr][strand].delete, sizeof(bool)*(narray), "mapfile->readarray[chr].delete");
    if(pwparam->enrichfile) p->readarray[chr][strand].ignore = (bool *)my_realloc(p->readarray[chr][strand].ignore, sizeof(bool)*(narray), "mapfile->readarray[chr].ignore");
  }

  /* fragment length distribution */
  if(p_frag->fraglen <0){
    fprintf(stderr, "ERROR: fragment length = %d < 0.\n", p_frag->fraglen);
    exit(1);
  }
  if(p_frag->fraglen < DIST_FRAGLEN_MAX) p->fstats.dist_fraglen[p_frag->fraglen]++;
  else p->fstats.dist_fraglen[DIST_FRAGLEN_MAX]++;

  p->readarray[chr][strand].narray = narray;
  return;
}

#ifdef READSV
static int check_sv(int sv, PwParam *p, char *name){
#else
static int check_sv(int sv, PwParam *p){
#endif
#ifdef READSV
  LOG("name:%s\n",name);
  // for paired-end
  LOG("   the read is paired in sequencing: %d\n",sv&1);
  LOG("   the read is mapped in a proper pair: %d\n",sv&2);
  LOG("   the query sequence itself is unmapped: %d\n",sv&4);
  LOG("   the mate is unmapped: %d\n",sv&8);
  LOG("   strand of the query (1 for reverse): %d\n",sv&16);
  LOG("   strand of the mate: %d\n",sv&32);
  LOG("   the read is the first read(F3) in a pair: %d\n",sv&64);
  LOG("   the read is the second read(F5) in a pair: %d\n",sv&128);
  LOG("   the alignment is not primary: %d\n",sv&256);
  LOG("   the read fails platform/vendor quality checks: %d\n",sv&512);
  LOG("   the read is either a PCR or an optical duplicate: %d\n",sv&1024);
  /*  LOG("   template having multiple segments in sequencing: %d\n",sv&1);
  LOG("   each segment properly aligned according to the aligner: %d\n",sv&2);
  LOG("   segment unmapped: %d\n",sv&4);
  LOG("   next segment in the template unmapped: %d\n",sv&8);
  LOG("   SEQ being reverse complemented: %d\n",sv&16);
  LOG("   SEQ of the next segment in the template being reversed: %d\n",sv&32);
  LOG("   the first segment in the template: %d\n",sv&64);
  LOG("   the last segment in the template: %d\n",sv&128);
  LOG("   secondary alignment: %d\n",sv&256);
  LOG("   not passing quality controls: %d\n",sv&512);
  LOG("   PCR or optical duplicate: %d\n",sv&1024);
  LOG("   supplementary alignment: %d\n",sv&2048);
  */
#endif

  // unmapped reads
  if(sv&4) goto err;
  // low quality reads
  if(sv&512 || sv&1024) goto err;
  if(p->rtype==READTYPE_PAIR){
    // unproper pair
    if(!(sv&2)) goto err;
    // unmatched pairs and interchromosomal pairs
    if(sv&8) goto err;
    // read pair mapped in same strand (for paired-end)
    if((sv&16 && sv&32) || (!(sv&16) && !(sv&32))) goto err;
  }
  
  return 0;
 err:
  return 1;
}

static void parse_sam(PwParam *p, Mapfile *mapfile, RefGenome *g, char *inputfile){
  FILE *IN=NULL;
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  int i, sv=0, nclm=0;
  char *command = alloc_str_new(inputfile, 100);
  FragmentData *p_frag = (FragmentData *)my_calloc(1, sizeof(FragmentData), "p_frag");
  
  /** SAM or BAM **/
  if(p->ftype==FILETYPE_SAM) IN = my_fopen(inputfile, FILE_MODE_READ);
  else if(p->ftype==FILETYPE_BAM){
    sprintf(command, "samtools view -h %s", inputfile);
    if(!(IN = popen(command, "r"))){
      fprintf(stderr,"error: cannot read %s.\n", inputfile);
      exit(0);
    }
  }
  
  /** Read each line **/
  while((fgets(str, STR_LEN, IN))!=NULL){ 
    if(str[0]=='\n' || str[0]=='@') continue;
    chomp(str);
    //    printf("%s\n",str);
    nclm = ParseLine(str, clm);
    sv = atoi(clm[1].str); // bitwise FLAG
#ifdef READSV
    if(check_sv(sv, p, clm[0].str)) continue;
#else
    if(check_sv(sv, p)) continue;
#endif

    init_frag(p_frag);
    p_frag->chr_F3 = changechr_str2int(clm[2].str, g);
    p_frag->name = strdup(clm[0].str);
    if(!p_frag->chr_F3){ printf("invalid chr name: %s\n", clm[2].str); exit(0);}
    for(i=0; i<nclm; i++){ if(!strncmp(clm[i].str, "NH:i:", 5)) p_frag->num_multimapped = atoi(clm[i].str+5);}
    
    if(p->rtype==READTYPE_SINGLE){  // single_end
      if(sv&16) p_frag->strand = STRAND_MINUS; else p_frag->strand = STRAND_PLUS;
      p_frag->readlen_F3 = strlen(clm[9].str);
      if(atoi(clm[8].str)) fprintf(stderr, "Warning: parsing paired-end file as single-end.\n");
      p_frag->F3      = atoi(clm[3].str);
      p_frag->fraglen = p->fraglen;
      add_fragment_to_readarray(p, mapfile, p_frag);
      mapfile->fstats.dist_readlen_F3[p_frag->readlen_F3]++;
    }else{                          // paired_end
      if(strcmp(clm[6].str, "="))              continue; // skip unmatched pairs and interchromosomal pairs
      if(atoi(clm[3].str) == atoi(clm[7].str)) continue; // skip single read
      if(sv&64){      // F3 read
	//if(!p_frag->chr_F5){ printf("invalid chr name: %s\n", clm[6].str); exit(0);}
	if(sv&16) p_frag->strand = STRAND_MINUS; else p_frag->strand = STRAND_PLUS;
	p_frag->readlen_F3 = strlen(clm[9].str);
	p_frag->F3      = atoi(clm[3].str);
	p_frag->F5      = atoi(clm[7].str);
	p_frag->fraglen = abs(atoi(clm[8].str));
	if(p_frag->fraglen <= p->max_fraglen && p_frag->fraglen > 0) add_fragment_to_readarray(p, mapfile, p_frag);
	mapfile->fstats.dist_readlen_F3[p_frag->readlen_F3]++;
      }else{          // F5 read
	p_frag->readlen_F5 = strlen(clm[9].str);
	mapfile->fstats.dist_readlen_F5[p_frag->readlen_F5]++;
      }
    }
  }

  fclose(IN);
  MYFREE(command);
  MYFREE(str);
  MYFREE(p_frag);
  return;
}

static void parse_bowtie(PwParam *p, Mapfile *mapfile, RefGenome *g, char *inputfile){
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  Elem clm[ELEM_NUM];
  int nclm=0;
  FragmentData *p_frag = (FragmentData *)my_calloc(1, sizeof(FragmentData), "p_frag");
  FILE *IN = my_fopen(inputfile, FILE_MODE_READ);
  char *tp=NULL, *tptemp=NULL;

  while((fgets(str, STR_LEN, IN))!=NULL){
    if(str[0]=='\n') continue;
    chomp(str);
    //    printf("%s\n",str);
    nclm = ParseLine(str, clm);
    if(nclm < 8){
      fprintf(stderr, "please use a bowtie file in which lines are not suppressed.\n");
      exit(0);
    }

    if(!tp) init_frag(p_frag); else MYFREE(p_frag->name);
    p_frag->name = strdup(clm[0].str);
    p_frag->num_multimapped = atoi(clm[6].str) +1;

    if(p->rtype==READTYPE_SINGLE){  // single_end
      p_frag->chr_F3 = changechr_str2int(clm[2].str, g);
      if(!p_frag->chr_F3){ printf("invalid chr name: %s\n", clm[2].str); exit(0);}
      p_frag->readlen_F3 = strlen(clm[4].str);
      if(strstr(clm[0].str, "/2")) fprintf(stderr, "Warning: parsing paired-end file as single-end.\n");
      p_frag->F3      = atoi(clm[3].str);
      p_frag->fraglen = p->fraglen;
      if(clm[1].str[0] == '+') p_frag->strand = STRAND_PLUS; else p_frag->strand = STRAND_MINUS;
      add_fragment_to_readarray(p, mapfile, p_frag);
      mapfile->fstats.dist_readlen_F3[p_frag->readlen_F3]++;
    }else{                           // paired_end
      if(!tp) tp = delimit_str(p_frag->name, ' ');
      else{
	tptemp = delimit_str(p_frag->name, ' ');
	if(strcmp(tp, tptemp)){
	  fprintf(stderr, "ERROR: Invalid read pair. %s - %s\n",tp, tptemp);
	  exit(1);
	}
	tp=NULL;
      }
      if(strstr(clm[0].str, "/1")){  // F3 read
	p_frag->chr_F3     = changechr_str2int(clm[2].str, g);
	if(!p_frag->chr_F3){ printf("invalid chr name: %s\n", clm[2].str); exit(0);}
	p_frag->readlen_F3 = strlen(clm[4].str);
	p_frag->F3         = atoi(clm[3].str);
	if(clm[1].str[0] == '+') p_frag->strand = STRAND_PLUS; else p_frag->strand = STRAND_MINUS;
	mapfile->fstats.dist_readlen_F3[p_frag->readlen_F3]++;
      }else{                         // F5 read
	p_frag->chr_F5     = changechr_str2int(clm[2].str, g);
	if(!p_frag->chr_F5){ printf("invalid chr name: %s\n", clm[2].str); exit(0);}
	p_frag->readlen_F5 = strlen(clm[4].str);
	p_frag->F5         = atoi(clm[3].str);
	mapfile->fstats.dist_readlen_F5[p_frag->readlen_F5]++;
      }
      if(tptemp){
	if(p_frag->strand == STRAND_PLUS) p_frag->fraglen = p_frag->F5 + p_frag->readlen_F5 - p_frag->F3;
	else p_frag->fraglen = p_frag->F3 + p_frag->readlen_F3 - p_frag->F5;
	if(p_frag->fraglen <= p->max_fraglen && p_frag->fraglen > 0) add_fragment_to_readarray(p, mapfile, p_frag);
	tptemp=NULL;
      }
    }
  }
  fclose(IN);
  MYFREE(str);
  MYFREE(p_frag);

  return;
}

static void parse_tagAlign(PwParam *p, Mapfile *mapfile, RefGenome *g, char *inputfile){
  char *str = (char *)my_calloc(STR_LEN, sizeof(char), "str");
  char *c=NULL;
  Elem clm[ELEM_NUM];
  int nclm=0;
  FragmentData *p_frag = (FragmentData *)my_calloc(1, sizeof(FragmentData), "p_frag");
  int zipped=0;
  if(strstr(inputfile, ".gz")) zipped=1;
  FILE *IN=NULL;
  struct gzFile_s *gzIN=NULL;
  if(zipped){
    if((gzIN = gzopen(inputfile, "r"))==NULL){
      fprintf(stderr,"[E] Cannot open <%s>.\n", inputfile); 
      exit(1);
    }
  }else{
    IN = my_fopen(inputfile, FILE_MODE_READ);
  }

  while(1) {
    if (zipped) c = gzgets(gzIN, str, STR_LEN);
    else        c = fgets(str, STR_LEN, IN);
    if(!c) break;

    if(str[0]=='\n') continue;
    chomp(str);
    nclm = ParseLine(str, clm);
    //    printf("%s %d\n",str, nclm);
    if(nclm < 6){
      fprintf(stderr, "please use tagAlign (BED3+3) file.\n");
      exit(0);
    }
    p_frag->num_multimapped = 1;

    if(p->rtype==READTYPE_SINGLE){  // single_end
      p_frag->chr_F3 = changechr_str2int(clm[0].str, g);
      if(!p_frag->chr_F3){ printf("invalid chr name: %s\n", clm[2].str); exit(0);}
      p_frag->F3      = atoi(clm[1].str);
      p_frag->readlen_F3 = abs(atoi(clm[2].str) - atoi(clm[1].str));
      p_frag->fraglen = p->fraglen;
      if(clm[5].str[0] == '+') p_frag->strand = STRAND_PLUS; else p_frag->strand = STRAND_MINUS;
      //      printf("%d %d %d %d\n", p_frag->F3, p_frag->readlen_F3, p_frag->fraglen, p_frag->strand);
      add_fragment_to_readarray(p, mapfile, p_frag);
      mapfile->fstats.dist_readlen_F3[p_frag->readlen_F3]++;
    }else{   // paired_end
      fprintf(stderr, "paired-end tagAlign format is not supported.\n");
      exit(0);
    }
  }
  if(zipped){
    gzclose(gzIN);
  }else{
    fclose(IN);
  }
  MYFREE(str);
  MYFREE(p_frag);

  return;
}

