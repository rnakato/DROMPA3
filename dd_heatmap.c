/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of DROMPA sources.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cairo.h>
#include <cairo-pdf.h>
#include <gtk/gtk.h>
#include "dd_heatmap.h"
#include "drompa_readfile.h"
#include "dd_readannotation.h"
#include "gsl_func.h"
#include "alloc.h"
#include "filehandle.h"
#include "macro.h"
#include "color.h"

#define GENEBLOCKNUM 10
#define FONTSIZE 20
#define MERGIN_BETWEEN_SAMPLE 30

typedef struct{
  gdouble *array;
  gint id;
} HMline;

typedef struct{
  HMline *line;
} HMarray;

gint binwidth;
gint sort_posi;
gdouble lat;            // height of each line
gdouble dot;            // width of each bin
gdouble wd4sample;      // width of each sample

#define rel_xline(cr, x1, y1, xlen) do{		\
    cairo_move_to(cr, x1, y1);			\
    cairo_rel_line_to(cr, xlen, 0); }while(0)
#define rel_yline(cr, x1, y1, ylen) do{		\
    cairo_move_to(cr, x1, y1);			\
    cairo_rel_line_to(cr, 0, ylen); }while(0)

static void draw_page(DrParam *p, DDParam *d, cairo_t *cr, SamplePair *sample, gint num, HMarray *hmarray, gdouble height_draw);
static void getvalue_eachsample(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample, HMarray *hmarray, gint sampleid, Gene *gene, gint num);
static void draw_eachsample(DDParam *d, cairo_t *cr, SamplePair *sample, HMarray *hmarray, gint sampleid, gint num, char *name, gdouble x, gint on);
static void show_indicatior(DDParam *d, cairo_t *cr, SamplePair *sample, gdouble height_draw);
static void stroke_ylabel(DDParam *d, cairo_t *cr, gint num);
static void define_color(cairo_t *cr, gdouble r);
static void showtext_cr(cairo_t *cr, gdouble x, gdouble y, gchar *str, gint fontsize);
static HMarray *HMarray_new(gint samplenum, gint sitenum, gint width);
static void HMarray_delete(HMarray *p, gint samplenum, gint sitenum);

static int HMlinecomp(const void *c1, const void *c2){
  HMline n1 = *(HMline *)c1;
  HMline n2 = *(HMline *)c2;
  if(n1.array[sort_posi] > n2.array[sort_posi]) return -1;
  else if(n1.array[sort_posi] == n2.array[sort_posi]) return 0;
  else return 1;
}

static int HMlinetotalcomp(const void *c1, const void *c2){
  HMline n1 = *(HMline *)c1;
  HMline n2 = *(HMline *)c2;
  int i;
  double sum1=0, sum2=0;
  for(i=GENEBLOCKNUM; i<GENEBLOCKNUM*2; i++){  // from TSS to TTS
    sum1 += n1.array[i];
    sum2 += n2.array[i];
  }
  if(sum1 > sum2) return -1;
  else if(sum1 == sum2) return 0;
  else return 1;
}

static gint get_siteset(DDParam *d, RefGenome *g, Gene **gene){
  gint i, chr, num=0, j=0;
  GeneSet *geneset=NULL;
  if(d->ptype == BEDSITES){
    for(i=0; i<d->bednum; i++) num += d->peak[i]->num;
  }else{
    geneset = (GeneSet *)my_calloc(g->chrnum, sizeof(GeneSet), "GeneSet");
    for(chr=1; chr<g->chrnum; chr++){
      geneset[chr].argv = d->gene.argv;
      geneset[chr].gene = (Gene *)my_calloc(STRUCT_GENE_MAX, sizeof(Gene), "gene");
      read_gene(&(geneset[chr]), g, chr, d->gftype);
      num += geneset[chr].num;
    }
    if(num){
      *gene = (Gene *)my_calloc(num, sizeof(Gene), "gene");
      for(chr=1; chr<g->chrnum; chr++){
	for(i=0; i< geneset[chr].num; i++) (*gene)[j++] = geneset[chr].gene[i];
	MYFREE(geneset[chr].gene);
      }
    }
  }
  if(!num){
    fprintf(stderr, "error: no target site.\n");
    exit(0);
  }

  return num;
}

static void getlat(DDParam *d, gint num){
  if(num < 100) lat = 20;
  else if(num < 1000) lat = 4;
  else if(num < 5000) lat = 0.8;
  else if(num < 20000) lat = 0.2;
  else if(num < 100000) lat = 0.04;
  else lat = 0.004;
  if(d->png) lat = max(0.1, lat);
  return;
}

static void copy_line(HMline *to, HMline *from){
  (*to).array = (*from).array;
  (*to).id = (*from).id;
  return;
}

static void sort_eachbedfile(DDParam *d, HMarray *p){
  gint i,j, diff=0, diff_nonzero=0;
  HMline *line;
  Peak *peak;
  for(i=0; i<d->bednum; i++){
    peak = d->peak[i];
    for(j=0; j<peak->num; j++){
      if(p->line[j+diff].array[sort_posi]) peak->num_nonzero++;
    }

    line = (HMline *)my_calloc(peak->num, sizeof(HMline), "hmline");
    // for(j=0; j<peak->num; j++) line[j] = p->line[j+diff];
        for(j=0; j<peak->num; j++) copy_line(&(line[j]), &(p->line[j+diff]));
    qsort(line, peak->num, sizeof(HMline), HMlinecomp);
    //for(j=0; j<peak->num_nonzero; j++) p->line[j+diff_nonzero] = line[j];
    for(j=0; j<peak->num_nonzero; j++) copy_line(&(p->line[j+diff_nonzero]), &(line[j])); 
    MYFREE(line);

    diff += peak->num;
    diff_nonzero += peak->num_nonzero;
  }
  return;
}

void draw_heatmap(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g){
  gint i, num;
  gchar filename[128];
  char *format;
  if(!d->png) format="pdf"; else format="png";
  sprintf(filename, "%s.heatmap.%s", p->headname, format);
  remove_file(filename);

  // global variant
  if(d->ptype != GENE100){
    dot = 4;
    binwidth  = d->cwbin*2;
    sort_posi = d->cwbin;
  }else{
    dot = 80/(gdouble)GENEBLOCKNUM;
    binwidth  = GENEBLOCKNUM*3;
    sort_posi = GENEBLOCKNUM;
  }
  wd4sample = binwidth * dot;

  /* hmarray */
  Gene *gene=NULL;
  gint sitenum = get_siteset(d, g, &gene);
  LOG("%s sitenum=%d\n",__FUNCTION__,sitenum);
  HMarray *hmarray = HMarray_new(p->samplenum, sitenum, binwidth);
  for(i=0; i<p->samplenum; i++) getvalue_eachsample(p, d, g, &(sample[i]), &(hmarray[i]), i, gene, sitenum);

  /* sort hmarray */
  gint sitenum_nonzero = 0;
  if(d->hmsort){
    if(d->ptype == BEDSITES){
      sort_eachbedfile(d, &(hmarray[d->hmsort-1]));
      for(i=0; i<d->bednum; i++) sitenum_nonzero += d->peak[i]->num_nonzero;
    }else{
      if(d->ptype == GENE100 && d->sortgbody) qsort(hmarray[d->hmsort-1].line, sitenum, sizeof(HMline), HMlinetotalcomp);
      else qsort(hmarray[d->hmsort-1].line, sitenum, sizeof(HMline), HMlinecomp);
      for(i=0; i<sitenum; i++){
	if(hmarray[d->hmsort-1].line[i].array[sort_posi]) sitenum_nonzero++;
      }
    }
  }
  if(!d->hmsort) num = sitenum;
  else           num = sitenum_nonzero;
  //  printf("num=%d, %d %d\n",num,sitenum, sitenum_nonzero);
  getlat(d, num);
  //  for(i=0; i<num; i++) printf("%d, %f %d\n", i, hmarray[d->hmsort-1].line[i].array[sort_posi], hmarray[d->hmsort].line[i].id);
  
  gdouble width_draw = wd4sample * p->samplenum + MERGIN_BETWEEN_SAMPLE * (p->samplenum-1);
  gdouble width_sheet = width_draw + OFFSET_X *2;
  gdouble height_draw = num * lat;
  gdouble height_sheet = max(100, height_draw + OFFSET_Y *2);
  
  // drawing
  cairo_surface_t *surface;
  if(!d->png) surface = cairo_pdf_surface_create(filename, width_sheet, height_sheet);
  else        surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width_sheet, height_sheet);

  FLUSH("Drawing...");
  cairo_t *cr = cairo_create(surface);
#ifdef DEBUG
  rel_xline(cr, 0, OFFSET_Y, width_sheet);
  rel_yline(cr, OFFSET_X, 0, height_sheet);
#endif

  draw_page(p, d, cr, sample, num, hmarray, height_draw);

  FLUSH("finishing...");
  cairo_destroy(cr);
  if(d->png) cairo_surface_write_to_png(surface, filename);
  cairo_surface_destroy(surface);

  printf("done.\n");
  HMarray_delete(hmarray, p->samplenum, num);
  return;
}

static void draw_page(DrParam *p, DDParam *d, cairo_t *cr, SamplePair *sample, gint num, HMarray *hmarray, gdouble height_draw){
  gint i;
  // background
  cairo_set_source_rgb(cr, CLR_WHITE);
  cairo_paint(cr);

  // ylabel
  stroke_ylabel(d, cr, num);

  gdouble x = OFFSET_X;
  gint on=0;
  for(i=0; i<p->samplenum; i++){
    draw_eachsample(d, cr, sample, hmarray, i, num, sample[i].linename, x, on);

    if(!on && strlen(sample[i].linename) > 30) on=1;
    else on=0;
    x += wd4sample + MERGIN_BETWEEN_SAMPLE;
  }

  // indicator
  show_indicatior(d, cr, sample, height_draw);

  return;
}

static gdouble getvalue(DDParam *d, TYPE_WIGARRAY *ChIP, TYPE_WIGARRAY *Input, gdouble ratio, gint posi){
  gdouble val=0;
  switch(d->stype){
  case READDIST:   val = WIGARRAY2VALUE(ChIP[posi]); break;
  case ENRICHDIST: val = CALCRATIO(ChIP[posi], Input[posi], ratio); break;
  case PVALUEDIST: val = binomial_test(WIGARRAY2VALUE(ChIP[posi]), WIGARRAY2VALUE(Input[posi]), ratio); break;
  }
  return val;
}

static void get_gene100(DDParam *d, SamplePair *sample, HMarray *hmarray, Gene *gene, gint chr, gint binsize, gint binnum, gint num){
  gint i,j,l, s=0, e=0;
  gdouble len, len100;
  TYPE_WIGARRAY arrayIP[binwidth], arrayInput[binwidth];
  for(i=0; i<num; i++){
    if(gene[i].chr != chr) continue;
    for(j=0; j<binwidth; j++){
      arrayIP[j]=0; arrayInput[j]=0;
    }
    len = gene[i].end - gene[i].start;
    len100 = len/GENEBLOCKNUM;
    for(j=0; j<binwidth; j++){
      if(gene[i].dir==1){
	s = (gene[i].start - len + len100 * j)    / binsize;
	e = (gene[i].start - len + len100 * (j+1))/ binsize;
      }else{
	s = (gene[i].end + len - len100 * (j+1))/ binsize;
	e = (gene[i].end + len - len100 * j)    / binsize;
      }
      if(s<0 || e>=binnum) continue;
      for(l=s; l<=e; l++){
	arrayIP[j] += sample->ChIP->data[l];
	if(d->stype != READDIST) arrayInput[j] += sample->Input->data[l];
      }
      hmarray->line[i].array[j] = getvalue(d, arrayIP, arrayInput, sample->comp->genome->ratio, j);
    }
  }
  return;
}

static void func(DDParam *d, SamplePair *sample, gdouble **array, gint wigposi, gint cwbin, gint binnum, gint dir){
  gint i;
  if(wigposi - cwbin < 0 || wigposi + cwbin >= binnum) return;
  for(i=-cwbin; i<cwbin; i++){
    (*array)[i+cwbin] = getvalue(d, sample->ChIP->data, sample->Input->data, sample->comp->genome->ratio, wigposi + i*dir);
  }
  return;
}

static void get_aroundpeak(DDParam *d, SamplePair *sample, HMarray *hmarray, gint chr, gint binsize, gint binnum){
  gint bid;
  gint i,wigposi,num=0;
  Peak *peak=NULL;
  for(bid=0; bid<d->bednum; bid++){
    peak = d->peak[bid];
    LOG("%s bid%d num%d\n",__FUNCTION__,bid,peak->num);
    for(i=0; i<peak->num; i++, num++){
      if(peak->bs[i].chr != chr) continue;
      wigposi = peak->bs[i].maxposi/binsize;
      func(d, sample, &(hmarray->line[num].array), wigposi, d->cwbin, binnum, 1);
    }
  }
  return;
}

static void get_aroundgene(DDParam *d, SamplePair *sample, HMarray *hmarray, Gene *gene, gint chr, gint binsize, gint binnum, gint num){
  gint i, wigposi;
  for(i=0; i<num; i++){
    if(gene[i].chr != chr) continue;
    if(gene[i].dir==1){
      if(d->ptype==TSS) wigposi = gene[i].start;
      else              wigposi = gene[i].end;
    }else{
      if(d->ptype==TSS) wigposi = gene[i].end;
      else              wigposi = gene[i].start;
    }
    wigposi /= binsize;
    func(d, sample, &(hmarray->line[i].array), wigposi, d->cwbin, binnum, gene->dir);
  }
  return;
}

static void getvalue_eachsample(DrParam *p, DDParam *d, RefGenome *g, SamplePair *sample, HMarray *hmarray, gint sampleid, Gene *gene, gint num){
  gint chr;

  FLUSH("sample %d: ", sampleid+1);
  for(chr=1; chr<g->chrnum; chr++){
    FLUSH("%s..", g->chr[chr].name);
    dr_read_wigdata(p, sample, sample->ChIP, g, chr);
    if(d->stype != READDIST){
      if(!sample->Input->argv){
	fprintf(stderr, "error: no Input specified for sample %d.\n", sampleid+1);
	exit(0);
      }
      dr_read_wigdata(p, sample, sample->Input, g, chr);
    }

    if(d->ptype == BEDSITES)                    get_aroundpeak(d, sample, hmarray, chr, sample->binsize, sample->binnum[chr]);
    else if(d->ptype == TSS || d->ptype == TTS) get_aroundgene(d, sample, hmarray, gene, chr, sample->binsize, sample->binnum[chr], num);
    else                                        get_gene100(d, sample, hmarray, gene, chr, sample->binsize, sample->binnum[chr], num);
   
    dr_delete_wigdata(sample->ChIP);
    if(d->stype != READDIST) dr_delete_wigdata(sample->Input);
  }

  printf("done.\n");
  return;
}

static void stroke_dot(DDParam *d, cairo_t *cr, SamplePair *sample, gdouble val, gdouble x, gdouble y){
  gdouble r=0;
  switch(d->stype){
  case READDIST:
    r = min(1, val/sample->scale_tag);
    break;
  case ENRICHDIST:
    val = val ? log10(val) : 0;
    r = max(-1, min(1, val/sample->scale_ratio));
    break;
  case PVALUEDIST:
    r = min(1, val/sample->scale_pvalue);
    break;
  }
  if(!r) return;

  define_color(cr, r);
  rel_xline(cr, x, y, dot);
  cairo_stroke(cr);
  return;
}

static void stroke_hmarray(DDParam *d, cairo_t *cr, SamplePair *sample, gdouble *array, gdouble x_ref, gdouble *y){
  gint i;
  gdouble x;
  for(i=0, x=x_ref; i<binwidth; i++, x += dot) stroke_dot(d, cr, sample, array[i], x, *y);
  *y += lat;
  return;
}

static void draw_eachsample(DDParam *d, cairo_t *cr, SamplePair *sample, HMarray *hmarray, gint sampleid, gint num, char *name, gdouble x, gint on){
  gint i,id;
  gdouble y=OFFSET_Y;
  // data
  cairo_set_line_width(cr, lat);

  if(d->ptype == BEDSITES){
    if(!d->hmsort){
      for(i=0; i<num; i++){
	stroke_hmarray(d, cr, sample, hmarray[sampleid].line[i].array, x, &y);
	//	if(i==1800)printf("\naaa %d %f\n",i,hmarray[sampleid].line[i].array[10]);
      }
    }else{
      for(i=0; i<num; i++){
	if(sampleid==d->hmsort-1) id = i;
	else        	          id = hmarray[d->hmsort-1].line[i].id;
	stroke_hmarray(d, cr, sample, hmarray[sampleid].line[id].array, x, &y);
      }
    }
  }else{
    for(i=0; i<num; i++){
      if(!d->hmsort) stroke_hmarray(d, cr, sample, hmarray[sampleid].line[i].array, x, &y);
      else{
	if(sampleid==d->hmsort-1) id = i;
	else  	                  id = hmarray[d->hmsort-1].line[i].id;
	stroke_hmarray(d, cr, sample, hmarray[sampleid].line[id].array, x, &y);
      }
    }
  }
  // frame and xlabel
  cairo_set_line_width(cr, 1.5);
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_rectangle(cr, x, OFFSET_Y, wd4sample, lat*num);
  cairo_stroke(cr);
  gdouble xname = max(x,x + wd4sample/2 - FONTSIZE * strlen(name)/2);
  showtext_cr(cr, xname, OFFSET_Y -5 - on * FONTSIZE, name, FONTSIZE);
  // stroke line between each bedfile
  if(d->ptype == BEDSITES){
    gint n=0;
    y = OFFSET_Y;
    for(i=0; i<d->bednum-1; i++){
      if(!d->hmsort) n = d->peak[i]->num;
      else           n = d->peak[i]->num_nonzero;
      y += lat * n;
      rel_xline(cr, x, y, wd4sample);
      cairo_stroke(cr);
    }
  }
  // xaxis scale
  char xscale[256];
  y = OFFSET_Y +lat*num +FONTSIZE +2;
  if(d->ptype != GENE100){
    sprintf(xscale, "-%.1fkb", d->compwidth/(double)NUM_1K);
    showtext_cr(cr, x, y, xscale, FONTSIZE);
    sprintf(xscale, "%.1fkb", d->compwidth/(double)NUM_1K);
    showtext_cr(cr, x +wd4sample - FONTSIZE * strlen(xscale)/2, y, xscale, FONTSIZE);
  }else{
    sprintf(xscale, "-100");
    showtext_cr(cr, x - 12, y, xscale, 16);
    sprintf(xscale, "0");
    showtext_cr(cr, x +GENEBLOCKNUM*dot - 2, y, xscale, 16);
    sprintf(xscale, "100");
    showtext_cr(cr, x +GENEBLOCKNUM*2*dot - 12, y, xscale, 16);
    sprintf(xscale, "200");
    showtext_cr(cr, x +wd4sample - 19, y, xscale, 16);
    sprintf(xscale, "%% length from TSS");
    showtext_cr(cr, x +7*dot, y+20, xscale, 16);
  }

  return;
}

static void show_indicatior(DDParam *d, cairo_t *cr, SamplePair *sample, gdouble height_draw){
  gdouble lw = 10, scale=0;
  char str[32];
  gdouble x = OFFSET_X - 60;
  gdouble y = min(OFFSET_Y + 20, OFFSET_Y + height_draw - 60);
  cairo_set_line_width(cr, lw);
  cairo_set_source_rgba(cr, CLR_BLACK, 1);

  if(d->stype==READDIST){
    sprintf(str, "ChIP reads");
    scale = sample->scale_tag;
  }else if(d->stype==ENRICHDIST){
    sprintf(str, "LogRatio");
    scale = sample->scale_ratio;
  }else if(d->stype==PVALUEDIST){
    sprintf(str, "P-value");
    scale = sample->scale_pvalue;
  }
  showtext_cr(cr, x-30, y-10, str, 14);

  gdouble r = 1;
  gint len=12;
  gdouble height=0;
  while(1){
    define_color(cr, r);
    rel_xline(cr, x, y, len);
    cairo_stroke(cr);
    sprintf(str, "%.1f", r*scale);
    cairo_set_source_rgba(cr, CLR_BLACK, 1);
    if(r==-1 || r==0 || r==1) showtext_cr(cr, x-35, y+lw/2, str, 12);
    r = (int)(r*10-1)/(double)10;  // add 0.1
    y += lw;
    height += lw;
    if((d->stype==READDIST || d->stype==PVALUEDIST) && r<0) break;
    if((d->stype==ENRICHDIST) && r<-1) break;
  }

  cairo_set_line_width(cr, 1);
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_rectangle(cr, x, y-height-lw/2, len, height);
  cairo_stroke(cr);

  return;
}

static void stroke_ylabel(DDParam *d, cairo_t *cr, gint num){
  gint i;
  gdouble height; 
  char str[128];
  gint x = OFFSET_X - 10;
  gint y = -OFFSET_Y - FONTSIZE * strlen(str)/2;
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_select_font_face(cr, FONTS_STYLE, CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, FONTSIZE);
  cairo_rotate(cr, -M_PI / 2);

  if(d->ptype == BEDSITES){
    for(i=0; i<d->bednum; i++){
      if(!d->hmsort) height = lat * d->peak[i]->num;
      else           height = lat * d->peak[i]->num_nonzero;
      sprintf(str, "%s", d->peak[i]->name);
      cairo_move_to(cr, y - height/2, x);  // rotateしているのでxとyが逆になる。yはマイナス
      cairo_show_text(cr, str);
      cairo_stroke(cr);
      y -= height;
    }
  }else{
    sprintf(str, "Genes");
    cairo_move_to(cr, y - lat*num/2, x);
    cairo_show_text(cr, str);
    cairo_stroke(cr);
  }
  cairo_rotate(cr, M_PI / 2);
  return;
}

static void define_color(cairo_t *cr, gdouble r){
  if(r<0) cairo_set_source_rgba(cr, 1+r*0.8, 1+r*0.8, 1,       1);
  else    cairo_set_source_rgba(cr, 1,       1-r*0.8, 1-r*0.8, 1);
  return;
}

static void showtext_cr(cairo_t *cr, gdouble x, gdouble y, gchar *str, gint fontsize){
  cairo_move_to(cr, x, y);
  cairo_select_font_face(cr, FONTS_STYLE, CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, fontsize);
  cairo_show_text(cr, str);
  cairo_stroke(cr);
  return;
}

static HMarray *HMarray_new(gint samplenum, gint sitenum, gint width){
  gint i,j;
  HMarray *p = (HMarray *)my_calloc(samplenum, sizeof(HMarray), "hmarray");
  for(i=0; i<samplenum; i++){
    p[i].line = (HMline *)my_calloc(sitenum, sizeof(HMline), "hmarray array*");
    for(j=0; j<sitenum; j++){
      p[i].line[j].array = (gdouble *)my_calloc(width, sizeof(gdouble), "hmarray array");
      p[i].line[j].id = j;
    }
  }
  return p;
}

static void HMarray_delete(HMarray *p, gint samplenum, gint sitenum){
  gint i,j;
  for(i=0; i<samplenum; i++){
    for(j=0; j<sitenum; j++) MYFREE(p[i].line[j].array);
    MYFREE(p[i].line);
  }
  MYFREE(p);
  return;
}
