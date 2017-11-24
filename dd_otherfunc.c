/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <string.h>
#include <math.h>
#include "dd_otherfunc.h"
#include "drompa_readfile.h"
#include "dd_readannotation.h"
#include "gsl_func.h"
#include "alloc.h"
#include "filehandle.h"
#include "stringp.h"
#include "macro.h"

#define LENMINUS4TR 30
#define LENPLUS4TR 300
#define NUM4HIST 5000
#define DEV4HIST 100

typedef struct{
  double tag_TSS, ave_TSS;
  double tag_body, ave_body;
  double TR;
} StPIndex;

void dd_counttags(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample){
  int i,j,k,s,e,chr;
  double sum, counttag_sum[p->samplenum];
  for(k=0; k<p->samplenum; k++) counttag_sum[k]=0;
  char *filename = alloc_str_new(p->headname, 20);
  sprintf(filename, "%s.xls", p->headname);
  remove_file(filename);

  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);
  fprintf(OUT, "chr\tstart\tend\t");
  for(i=0; i<p->samplenum; i++){
    fprintf(OUT, "%s", sample[i].linename);
    if(i==p->samplenum-1) fprintf(OUT, "\n"); else fprintf(OUT, "\t");
  }

  BedChr *bed=NULL;
  for(chr=1; chr<g->chrnum; chr++){
    if(!p->includeYM && (!strcmp(g->chr[chr].name, "chrY") || !strcmp(g->chr[chr].name, "chrM") || !strcmp(g->chr[chr].name, "chrMT"))) continue;
    FLUSH("%s..", g->chr[chr].name);
    for(i=0; i<p->samplenum; i++){
      if(sample[i].copyC==-1) dr_read_wigdata(p, &(sample[i]), sample[i].ChIP, g, chr);
    }
    bed = &(d->bed[0]->chr[chr]);
    for(i=0; i<bed->num; i++){
      //      fprintf(OUT, "%s\t%d\t%d\t", g->chr[chr].name, bed->bed[i].s, bed->bed[i].e);
      for(k=0; k<p->samplenum; k++){
	s = bed->bed[i].s / sample[k].binsize;
	e = bed->bed[i].e / sample[k].binsize;
	sum=0;
	for(j=s; j<=e; j++) sum += WIGARRAY2VALUE(sample[k].ChIP->data[j]);
	counttag_sum[k] += sum;
	//	fprintf(OUT, "%.1f", sum);
	//if(k==p->samplenum-1) fprintf(OUT, "\n"); else fprintf(OUT, "\t");
      }
    }
    for(i=0; i<p->samplenum; i++){
      if(sample[i].copyC==-1) dr_delete_wigdata(sample[i].ChIP);
    }
  }
  fprintf(OUT, "total tags\t\t\t");
  for(i=0; i<p->samplenum; i++){
    fprintf(OUT, "%.1f (%.1f%%)", counttag_sum[i], counttag_sum[i]*100/(double)sample[i].ChIP->genome->nread);
    if(i==p->samplenum-1) fprintf(OUT, "\n"); else fprintf(OUT, "\t");
  }
  fprintf(OUT, "bed regions total (bp)\t%ld\n", d->bed[0]->len_total);

  printf("done.\n");
  fclose(OUT);
  MYFREE(filename);
  return;
}

void dd_compare_genebody(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample)
{
  int i,j,k,l,s,e, chr, len;
  double IPtag[p->samplenum], exontag[p->samplenum], introntag[p->samplenum];
  char *filename = alloc_str_new(p->headname, 20);
  sprintf(filename, "%s.xls", p->headname);
  remove_file(filename);
  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);
  Gene *gene;

  fprintf(OUT, "name\tchromosome\tstart\tend\tlength\t");
  for(i=0; i<p->samplenum; i++){
    fprintf(OUT, "%s\tper kbp\texon\tintron\t intron/exon", sample[i].linename);
    if(i==p->samplenum-1) fprintf(OUT, "\n"); else fprintf(OUT, "\t");
  }
  for(chr=1; chr<g->chrnum; chr++){
    if(!p->includeYM && (!strcmp(g->chr[chr].name, "chrY") || !strcmp(g->chr[chr].name, "chrM") || !strcmp(g->chr[chr].name, "chrMT"))) continue;
    FLUSH("%s..", g->chr[chr].name);
    d->gene.num=0;
    d->gene.gene = (Gene *)my_calloc(STRUCT_GENE_MAX, sizeof(Gene), "gene");
    read_gene(&(d->gene), g, chr, d->gftype);
    for(i=0; i<p->samplenum; i++){
      if(sample[i].copyC==-1) dr_read_wigdata(p, &(sample[i]), sample[i].ChIP, g, chr);
    }
    gene = d->gene.gene;

    for(i=0; i<d->gene.num; i++){
      for(j=0; j<p->samplenum; j++){
	IPtag[j]=0;
	exontag[j]=0;
	introntag[j]=0;
      }
      /* gene body */
      s = gene[i].start / sample[0].binsize;
      e = gene[i].end   / sample[0].binsize;
      len = gene[i].end - gene[i].start;
      for(k=0; k<p->samplenum; k++){
	for(j=s; j<=e; j++) IPtag[k] += WIGARRAY2VALUE(sample[k].ChIP->data[j]);
      }
      /* exon, intron */
      for(l=0; l< gene[i].exonnum; l++){
	for(k=0; k<p->samplenum; k++){
	  s = gene[i].exon[l].start / sample[k].binsize;  /*exon*/
	  e = gene[i].exon[l].end   / sample[k].binsize;
	  for(j=s; j<=e; j++) exontag[k] += WIGARRAY2VALUE(sample[k].ChIP->data[j]);
	}
	if(l < gene[i].exonnum -1){   /*intron*/
	  for(k=0; k<p->samplenum; k++){
	  s = gene[i].exon[l].end     / sample[k].binsize;
	  e = gene[i].exon[l+1].start / sample[k].binsize;
	    for(j=s; j<=e; j++) introntag[k] += WIGARRAY2VALUE(sample[k].ChIP->data[j]);
	  }
	}
      }
      int on=0;
      for(k=0; k<p->samplenum; k++){
	if(IPtag[k]/(double)len*NUM_1K < d->cgthre) on=1;
      }
      if(on) continue;
      fprintf(OUT, "%s\t%s\t%d\t%d\t%d\t",gene[i].name, g->chr[chr].name, gene[i].start, gene[i].end, len);
      for(j=0; j<p->samplenum; j++){
	fprintf(OUT, "%.2f\t%.2f", IPtag[j], IPtag[j]/len*NUM_1K);
	if(d->gene.gene[i].exonnum){
	  fprintf(OUT, "\t%.2f\t%.2f\t%.2f", exontag[j], introntag[j], exontag[j] ? introntag[j]/exontag[j]: 0);
	}else fprintf(OUT, "\t\t\t");
	if(j==p->samplenum-1) fprintf(OUT, "\n"); else fprintf(OUT, "\t");
      }
    }
    for(i=0; i<p->samplenum; i++){
      if(sample[i].copyC==-1) dr_delete_wigdata(sample[i].ChIP);
    }

    MYFREE(gene);
  }
  printf("done.\n");
  fclose(OUT);
  MYFREE(filename);
  return;
}

void dd_compare_intensity(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample)
{
  int i,j,s,e,n,chr;
  double IPtag1=0, IPtag2=0,A=0,M=0, max, min, pval_enr;
  char *filename = alloc_str_new(p->headname, 20);
  sprintf(filename, "%s.xls", p->headname);
  remove_file(filename);
  BedChr *bed=NULL;

  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);
  fprintf(OUT, "chromosome\tstart\tend\tpeak width\t%s\t%s\tlog average(A)\tlog enrichment (M)\t-log10(p)\n", sample[0].linename, sample[1].linename);

  for(chr=1; chr<g->chrnum; chr++){
    if(!p->includeYM && (!strcmp(g->chr[chr].name, "chrY") || !strcmp(g->chr[chr].name, "chrM") || !strcmp(g->chr[chr].name, "chrMT"))) continue;

    printf("%s..", g->chr[chr].name);
    int nbedsites=0;
    for (j=0; j<d->bednum; j++) nbedsites += d->bed[j]->chr[chr].num;
    if (!nbedsites) {
      printf("(no site)\n");
      continue;
    }
    for (i=0; i<2; i++) {
      printf("load sample%d..", i+1);
      if(sample[i].copyC==-1) dr_read_wigdata(p, &(sample[i]), sample[i].ChIP, g, chr);
    }
    for(j=0; j<d->bednum; j++){
      bed = &(d->bed[j]->chr[chr]);
      for(i=0; i<bed->num; i++){
	s = bed->bed[i].s / sample[0].binsize;
	e = bed->bed[i].e / sample[0].binsize;
	n = e - s + 1;
	IPtag1=0; IPtag2=0;
	for(; s<=e; s++){
	  IPtag1 += WIGARRAY2VALUE(sample[0].ChIP->data[s]);
	  IPtag2 += WIGARRAY2VALUE(sample[1].ChIP->data[s]);
	}
	IPtag1 /= n;
	IPtag2 /= n;
	A = log2((IPtag1 +1)*(IPtag2+1))/2;
	M = log2(IPtag1 +1) - log2(IPtag2+1);
	max = max(IPtag1, IPtag2);
	min = min(IPtag1, IPtag2);
	pval_enr = binomial_test(max, min, 1);
	fprintf(OUT, "%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", g->chr[chr].name, bed->bed[i].s, bed->bed[i].e, bed->bed[i].e-bed->bed[i].s, IPtag1, IPtag2, A, M, pval_enr);
      }
    }
    for(i=0; i<2; i++){
      if(sample[i].copyC==-1) dr_delete_wigdata(sample[i].ChIP);
    }
    printf("\n");
  }
  printf("done.\n");
  fclose(OUT);
  MYFREE(filename);
  return;
}

void dd_multici(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample)
{
  int i,j,k,l,s,e,n,chr;
  char *filename = alloc_str_new(p->headname, 20);
  sprintf(filename, "%s.xls", p->headname);
  remove_file(filename);
  BedChr *bed=NULL;

  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);

  fprintf(OUT, "chromosome\tstart\tend\tpeak width");
  for (i=0; i<p->samplenum; i++) fprintf(OUT, "\t%s", sample[i].linename);
  fprintf(OUT, "\n");

  for (chr=1; chr<g->chrnum; chr++) {
    if (!p->includeYM && (!strcmp(g->chr[chr].name, "chrY") || !strcmp(g->chr[chr].name, "chrM") || !strcmp(g->chr[chr].name, "chrMT"))) continue;
    printf("%s..", g->chr[chr].name);
    int nbedsites=0;
    for (j=0; j<d->bednum; j++) nbedsites += d->bed[j]->chr[chr].num;
    if (!nbedsites) {
      printf("(no site)\n");
      continue;
    }
    for (i=0; i<p->samplenum; i++) {
      printf("load sample%d..", i+1);
      if(sample[i].copyC==-1) dr_read_wigdata(p, &(sample[i]), sample[i].ChIP, g, chr);
    }
    for (j=0; j<d->bednum; j++) {
      bed = &(d->bed[j]->chr[chr]);
      for (i=0; i<bed->num; i++) {
	//	printf("%d %d\n", i, bed->num);
	s = bed->bed[i].s / sample[0].binsize;
	e = bed->bed[i].e / sample[0].binsize;
	n = e - s + 1;
	//	printf("s%d e%d n%d\n", s, e, n);
	fprintf(OUT, "%s\t%d\t%d\t%d", g->chr[chr].name, bed->bed[i].s, bed->bed[i].e, bed->bed[i].e-bed->bed[i].s);
	double tag=0;
	for (k=0; k<p->samplenum; k++) {
	  for (l=s; l<=e; l++) tag += WIGARRAY2VALUE(sample[k].ChIP->data[l]);
	  fprintf(OUT, "\t%.2f", tag/n);
	}
	fprintf(OUT, "\n");
      }
    }
    for (i=0; i<p->samplenum; i++) {
      if (sample[i].copyC==-1) dr_delete_wigdata(sample[i].ChIP);
    }
    printf("\n");
  }
  
  printf("done.\n");
  fclose(OUT);
  MYFREE(filename);
  return;
}

static void draw_TRgraph(DrParam *p, SamplePair *sample, int hist[][NUM4HIST], int ngene){
  int i,j;
  int value[p->samplenum];
  for(j=0; j<p->samplenum; j++) value[j]=0;
  char *filename = alloc_str_new(p->headname, 20);
  sprintf(filename, "%s.fig.xls", p->headname);
  remove_file(filename);

  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);
  for(i=0; i<p->samplenum; i++) fprintf(OUT, "\t%s", sample[i].linename);
  fprintf(OUT, "\n");
  for(i=0; i<NUM4HIST; i++){
    fprintf(OUT, "%.2f", i/(double)DEV4HIST);
    for(j=0; j<p->samplenum; j++){
      value[j] += hist[j][i];
      fprintf(OUT, "\t%.2f", value[j]/(double)ngene);
    }
    fprintf(OUT, "\n"); 
  }
  fclose(OUT);
  MYFREE(filename);
  return;
}

static void init_pi(StPIndex *pi){
  pi->tag_TSS = pi->ave_TSS = 0;
  pi->tag_body = pi->ave_body = 0;
  pi->TR = 0;
  return;
}

void dd_travelling_ratio(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample){
  int i,j,k, chr, num, ngene=0, s_tss, e_tss, s_body, e_body;
  int on;
  int binsize = sample[0].binsize;
  char *filename = alloc_str_new(p->headname, 20);
  sprintf(filename, "%s.xls", p->headname);
  remove_file(filename);
  StPIndex pi[p->samplenum];
  
  Gene *gene=NULL;
  int hist[p->samplenum][NUM4HIST];
  for(i=0; i<p->samplenum; i++){
    for(j=0; j<NUM4HIST; j++) hist[i][j]=0;
  }
  
  FILE *OUT = my_fopen(filename, FILE_MODE_WRITE);
  fprintf(OUT, "\t\t\t");
  for(i=0; i<p->samplenum; i++) fprintf(OUT, "%s\t\t\t\t\t", sample[i].linename);
  fprintf(OUT, "\n");

  fprintf(OUT, "ID\tname\tchromosome");
  for(i=0; i<p->samplenum; i++) fprintf(OUT, "\tTSS\tave\tgenebody\tave\tTSS/body (TR)");
  fprintf(OUT, "\n");
  
  for(chr=1; chr<g->chrnum; chr++){
    if(!p->includeYM && (!strcmp(g->chr[chr].name, "chrY") || !strcmp(g->chr[chr].name, "chrM") || !strcmp(g->chr[chr].name, "chrMT"))) continue;
    FLUSH("%s..", g->chr[chr].name);
    if(d->gene.argv){
      d->gene.gene = (Gene *)my_calloc(STRUCT_GENE_MAX, sizeof(Gene), "gene");
      read_gene(&(d->gene), g, chr, d->gftype);
    }
    for(i=0; i<p->samplenum; i++){
      if(sample[i].copyC==-1) dr_read_wigdata(p, &(sample[i]), sample[i].ChIP, g, chr);
    }
    
    gene = d->gene.gene;
    for(i=0; i<d->gene.num; i++){
      if(gene[i].dir==1){
	s_tss = (gene[i].start - LENMINUS4TR)/ binsize;
	e_tss = (gene[i].start + LENPLUS4TR) / binsize;
	s_body = e_tss;
	e_body = gene[i].end / binsize;
      }else{
	s_tss = (gene[i].end - LENPLUS4TR) / binsize;
	e_tss = (gene[i].end + LENMINUS4TR)/ binsize;
	s_body = gene[i].start / binsize;
	e_body = s_tss;
      }
      if(s_body >= e_body || s_tss >= e_tss) continue;
      on=0;
      for(k=0; k<p->samplenum; k++){
	init_pi(&(pi[k]));
	for(j=s_tss; j<=e_tss; j++) pi[k].tag_TSS += WIGARRAY2VALUE(sample[k].ChIP->data[j]);
	pi[k].ave_TSS = pi[k].tag_TSS/(e_tss - s_tss +1);
	for(j=s_body; j<=e_body; j++) pi[k].tag_body += WIGARRAY2VALUE(sample[k].ChIP->data[j]);
	pi[k].ave_body = pi[k].tag_body/(e_body - s_body +1);
	pi[k].TR = pi[k].ave_body ? pi[k].ave_TSS/pi[k].ave_body: 0;
	//pi[k].TR = (pi[k].tag_TSS+1) / (pi[k].tag_body+1);   /* stalling score (ZINBA) */
	if(pi[k].ave_TSS >d->tssthre && pi[k].ave_body) on=1;
      }
      if(!on) continue;
      fprintf(OUT, "%s\t%s\t%s", gene[i].ID, gene[i].name, g->chr[chr].name);
      for(k=0; k<p->samplenum; k++){
	fprintf(OUT, "\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", pi[k].tag_TSS, pi[k].ave_TSS, pi[k].tag_body, pi[k].ave_body, pi[k].TR);
	num = min(pi[k].TR*DEV4HIST, NUM4HIST-1);
	hist[k][num]++;
      }
      fprintf(OUT, "\n");
      ngene++;
    }
    for(i=0; i<p->samplenum; i++){
      if(sample[i].copyC==-1) dr_delete_wigdata(sample[i].ChIP);
    }
    MYFREE(gene);
    d->gene.num=0;
  }
  fclose(OUT);
  MYFREE(filename);

  draw_TRgraph(p, sample, hist, ngene);
  printf("done.\n");
  return;
}
