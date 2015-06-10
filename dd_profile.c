/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdio.h>
#include <math.h>
#include "dd_profile.h"
#include "alloc.h"
#include "filehandle.h"
#include "stringp.h"
#include "macro.h"
#include "color.h"

#define GENEBLOCKNUM 100
static double **pdarray_new(int snum, int arraynum);
static void delete_pdarray(double **p, int snum);

static void func(DDParam *d, SamplePair *sample, int samplenum, int i, int posi){
  int k;
  for(k=0; k<samplenum; k++){
    sample[k].profile.IP[i] += sample[k].ChIP->data[posi];
    sample[k].profile.IPsum += sample[k].ChIP->data[posi];
    if(d->showse) sample[k].profile.SE[i][d->ntotal_profile] = sample[k].ChIP->data[posi];
    if(d->stype == ENRICHDIST){
      sample[k].profile.Input[i] += sample[k].Input->data[posi];
      sample[k].profile.Inputsum += sample[k].Input->data[posi];
    }
  }
  return;
}

static void output_detail(DDParam *d, SamplePair *sample, char *sitename, int samplenum, int arraynum, char *chrname){
  int i,k;
  FILE *OUT = my_fopen(d->filename_profile, FILE_MODE_A);
  for(k=0; k<samplenum; k++){
    fprintf(OUT, "%s\t%s\t%s", chrname, sitename, sample[k].linename);
    for(i=0; i<arraynum; i++){
      if(d->stype == ENRICHDIST) fprintf(OUT, "\t%.2f", sample[k].profile.Input[i] ? sample[k].profile.IP[i]/(double)sample[k].profile.Input[i]: 0);
      else fprintf(OUT, "\t%.2f", WIGARRAY2VALUE(sample[k].profile.IP[i]));
    }
    fprintf(OUT, "\n");
  }
  fclose(OUT);
  return;
}

static void add_around_sites(DDParam *d, RefGenome *g, SamplePair *sample, int samplenum, int wigposi, char *sitename, int binsize, int binnum, int cwbin, int chr, int dir){
  int i,j,k;
  wigposi /= binsize;
  if(wigposi - cwbin < 0 || wigposi + cwbin >= binnum){
    d->ntotal_skip++;
    return;
  }
  for(i=0; i<cwbin*2+1; i++) func(d, sample, samplenum, i, wigposi - (cwbin-i)*dir);
  if(d->pdetail) output_detail(d, sample, sitename, samplenum, cwbin*2+1, g->chr[chr].name);
  d->ntotal_profile++;

  if(d->showse && d->ntotal_profile>= d->SEnum){
    d->SEnum += PEAKNUM_DEFAULT;
    for(k=0; k<samplenum; k++){
      for(j=0; j<cwbin*2+1; j++){
	sample[k].profile.SE[j] = (double *)my_realloc(sample[k].profile.SE[j], d->SEnum*sizeof(double), "sample[k].profile.SE[j]");
      }
    }
  }
  return;
}

void add_profile(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample, int chr){
  int i,j,k, wigposi;
  int samplenum = p->samplenum;
  int binsize = sample[0].binsize;
  int cwbin = d->cwbin;

  if(d->ptype==GENE100){  // gene100
    int len, s=0, e=0, arraynum=GENEBLOCKNUM*3;
    double len100=0;
    double **pdarray = pdarray_new(samplenum, arraynum);
    Gene *gene = d->gene.gene;

    for(i=0; i<d->gene.num; i++){
      len = gene[i].end - gene[i].start;
      if(gene[i].end + len > g->chr[chr].len || gene[i].start - len < 0){
	d->ntotal_skip++;
	continue;
      }
      len100 = len/(double)GENEBLOCKNUM;
      for(j=0; j<arraynum; j++){
	if(gene[i].dir==1){
	  s = (gene[i].start - len + len100 *j)    / binsize;
	  e = (gene[i].start - len + len100 *(j+1))/ binsize;
	}else{
	  s = (gene[i].end + len - len100 * (j+1))/ binsize;
	  e = (gene[i].end + len - len100 * j)    / binsize;
	}
	for(k=s; k<=e; k++) func(d, sample, samplenum, j, k);
      }
      if(d->pdetail) output_detail(d, sample, gene[i].name, samplenum, arraynum, g->chr[chr].name);
      d->ntotal_profile++;

      if(d->showse && d->ntotal_profile>= d->SEnum){
	d->SEnum += PEAKNUM_DEFAULT;
	for(k=0; k<samplenum; k++){
	  for(j=0; j<arraynum; j++){
	    sample[k].profile.SE[j] = (double *)my_realloc(sample[k].profile.SE[j], d->SEnum*sizeof(double), "sample[k].profile.SE[j]");
	  }
	}
      }

    }
    delete_pdarray(pdarray, samplenum);
  }
  else if(d->ptype==BEDSITES){   // bedfile
    Peak *peak=NULL;
    char sitename[300];
    for(i=0; i<d->bednum; i++){
      peak = d->peak[i];
      for(j=0; j<peak->num; j++){
	if(peak->bs[j].chr != chr) continue;
	wigposi = peak->bs[j].maxposi;
	sprintf(sitename, "%d", wigposi);
	add_around_sites(d, g, sample, samplenum, wigposi, sitename, binsize, sample->binnum[chr], cwbin, chr, 1);
      }
    }
  }
  else{    // TSS,TTS
    Gene *gene = d->gene.gene;
    for(i=0; i<d->gene.num; i++){
      if(gene[i].dir==1){
	if(d->ptype==TSS) wigposi = gene[i].start;
	else              wigposi = gene[i].end;
      }else{
	if(d->ptype==TSS) wigposi = gene[i].end;
	else              wigposi = gene[i].start;
      }
      add_around_sites(d, g, sample, samplenum, wigposi, gene[i].name, binsize, sample->binnum[chr], cwbin, chr, gene[i].dir);
    }
  }
  return;
}

static void output_Rfile(DrParam *p, DDParam *d, SamplePair *sample, int num, int samplenum){
  int i,j;
  double color[]={CLR_RED, CLR_BLUE, CLR_GREEN, CLR_BLACK, CLR_PINK, CLR_DARKORANGE, CLR_PURPLE, CLR_GRAY3, CLR_OLIVE, CLR_YELLOW3, CLR_LIGHTCYAN, CLR_LIGHTCORAL, CLR_SALMON, CLR_GREEN2, CLR_BLUE3, CLR_PURPLE2};
  char *posistr[]={"",
		   "Distance from TSS (bp)",
		   "Distance from TTS (bp)",
		   "%% gene length from TSS",
		   "Distance from the peak summit (bp)"};

  FILE *OUT = my_fopen(d->filename_profile, FILE_MODE_WRITE);
  fprintf(OUT, "# bsnum_allowed=%d, bsnum_skipped=%d\n", d->ntotal_profile, d->ntotal_skip);
  for(i=0; i<samplenum; i++){
    fprintf(OUT, "p%d <- c(", i+1);
    for(j=0; j<num; j++){
      fprintf(OUT, "%.3f", sample[i].profile.IP[j]);
      if(j != num-1) fprintf(OUT, ","); else fprintf(OUT, ")\n");
    }

    if(d->showse){
      fprintf(OUT, "p%d_SE <- c(", i+1);
      for(j=0; j<num; j++){
	fprintf(OUT, "%.3f", 1.96 * sample[i].profile.SEarray[j]);  // 95% CI
	if(j != num-1) fprintf(OUT, ","); else fprintf(OUT, ")\n");
      }
      fprintf(OUT, "p%d_upper <- p%d + p%d_SE\n", i+1, i+1, i+1);
      fprintf(OUT, "p%d_lower <- p%d - p%d_SE\n", i+1, i+1, i+1);
    }
  }
  if(d->ptype != GENE100) fprintf(OUT, "x <- seq.int(from=-%d,to=%d,by=%d)\n", d->compwidth, d->compwidth, sample[0].binsize);
  else                    fprintf(OUT, "x <- seq.int(from=-%d,to=%d,by=1)\n", GENEBLOCKNUM, GENEBLOCKNUM*2-1);
  fprintf(OUT, "ymax <- ceiling(max(c(");
  for(i=0; i<samplenum; i++){
    fprintf(OUT, "max(p%d)", i+1);
    if(i<samplenum-1) fprintf(OUT, ","); else fprintf(OUT, ")))\n");
  }
  fprintf(OUT, "pdf(\"%s.pdf\",6,6)\n", p->headname);

  double min=1;
  for(i=0; i<samplenum; i++){
    for(j=0; j<num; j++){
      if(min > sample[i].profile.IP[j]) min = sample[i].profile.IP[j];
    }
  }
  min = (int)(min*10)/10.0;
  if(!min) min = 0.1;

  /* datalines */
  fprintf(OUT, "plot(x,p%d,type=\"l\",col=rgb(%.3f,%.3f,%.3f),xlab=\"%s\",ylab=\"read density\", ",1, color[0], color[1], color[2], posistr[d->ptype]);
  if(d->ptype == GENE100) fprintf(OUT, "log=\"y\", ylim=c(%.1f, ymax))\n", min);
  else                    fprintf(OUT, "ylim=c(0, ymax))\n");  

  for(i=1; i<samplenum; i++) fprintf(OUT, "lines(x,p%d,col=rgb(%.3f,%.3f,%.3f))\n", i+1, color[i*3], color[i*3+1], color[i*3+2]);
  for(i=0; i<samplenum; i++){
    fprintf(OUT, "polygon(c(x, rev(x)), c(p%d_lower, rev(p%d_upper)), col=rgb(%.3f,%.3f,%.3f,0.3), border=NA)\n", i+1, i+1, color[i*3], color[i*3+1], color[i*3+2]);
  }

  /* legend */
  fprintf(OUT, "legend(\"bottomleft\",c(");
  for(i=0; i<samplenum; i++){
    fprintf(OUT, "\"%s\"", sample[i].linename);
    if(i<samplenum-1) fprintf(OUT, ","); else fprintf(OUT, "),lty=c(");
  }
  for(i=0; i<samplenum; i++){
    fprintf(OUT, "1");
    if(i<samplenum-1) fprintf(OUT, ","); else fprintf(OUT, "),col=c(");
  }
  for(i=0; i<samplenum; i++){
    fprintf(OUT, "rgb(%.3f,%.3f,%.3f)", color[i*3], color[i*3+1], color[i*3+2]);
    if(i<samplenum-1) fprintf(OUT, ","); else fprintf(OUT, "), lwd=1.5)\n");
  }
  fclose(OUT);
  return;
}

static double calc_SE(double *array, int num){
  int i;
  if(!num) return 0;
  double ave=0;
  for(i=0; i<num; i++){
    ave += array[i];
  }
  ave = WIGARRAY2VALUE(ave)/num;
  double var=0, diff;
  for(i=0; i<num; i++){
    diff = WIGARRAY2VALUE(array[i]) - ave;
    var += diff * diff;
  }
  var /= num-1; 
  double se = sqrt(var/num);
  return se;
}

/* profiles and standard errors (condidence interval) */
void show_profile(DrParam *p, DDParam *d, SamplePair *sample){
  int i,j, num;
  double r=0;
  int cwbin = d->cwbin;
  int samplenum = p->samplenum;

  double s=0;
  for(i=0; i<samplenum; i++){
    if(!i) s = sample[i].profile.IPsum;
    else if(s > sample[i].profile.IPsum) s = sample[i].profile.IPsum;
  }
  if(d->ptype == GENE100) num = GENEBLOCKNUM*3; else num = cwbin*2+1;

  for(i=0; i<samplenum; i++){
    if(d->stype == READDIST){
      if(!d->ntype) r = 1 /(double)(d->ntotal_profile);
      else          r = s / sample[i].profile.IPsum /(double)(d->ntotal_profile);
      for(j=0; j<num; j++){
	sample[i].profile.IP[j] = WIGARRAY2VALUE(sample[i].profile.IP[j]) * r;
      }
    }
    else if(d->stype == ENRICHDIST){
      for(j=0; j<num; j++){ 
	sample[i].profile.IP[j] = CALCRATIO(sample[i].profile.IP[j], sample[i].profile.Input[j], sample->comp->genome->ratio);
      }
    }
  }

  /* SE */
  if(d->showse){
    for(i=0; i<samplenum; i++){
      for(j=0; j<num; j++){ 
	sample[i].profile.SEarray[j] = calc_SE(sample[i].profile.SE[j], d->ntotal_profile);
      }
    }
  }

  output_Rfile(p, d, sample, num, samplenum);

  /* R command */
  char *command = alloc_str_new(d->filename_profile, 50);
  sprintf(command, "R --vanilla < %s", d->filename_profile);
  printf("\nExecute command: %s\n", command);
  my_system(command);
  MYFREE(command);
  return;
}


static double **pdarray_new(int snum, int arraynum){
  int i;
  double **p = (double **)my_calloc(snum, sizeof(double *), "pdarray *");
  for(i=0; i<snum; i++){
    p[i] = (double *)my_calloc(arraynum, sizeof(double), "pdarray");
  }
  return p;
}

static void delete_pdarray(double **p, int snum){
  int i;
  for(i=0; i<snum; i++){ MYFREE(p[i]);}
  MYFREE(p);
  return;
}
