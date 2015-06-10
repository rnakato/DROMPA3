/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp> This file
 * is a part of DROMPA sources.
 */
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <cairo.h>
#include <cairo-pdf.h>
#include <cairo-ps.h>
#include <cairo-svg.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include "dd_stroke.h"
#include "common.h"
#include "alloc.h"
#include "stringp.h"
#include "dp_call.h"
#include "filehandle.h"
#include "macro.h"
#include "color.h"
#include "gsl_func.h"

#define BOXHEIGHT_GRAPH 80
#define BOXHEIGHT_INTERACTION 50
#define BOXHEIGHT_PD 40
#define BOXHEIGHT_GENEBOX_EXON 130
#define BOXHEIGHT_GENEBOX_NOEXON 80
#define LINEHEIGHT_BED 10
#define HEIGHT_REPEAT 10
#define MERGIN_BETWEEN_GRAPH_DATA 15
#define MERGIN_BETWEEN_DATA 6
#define MERGIN_BETWEEN_LINE 30
#define MERGIN_BETWEEN_GOVERLOOK 60
#define MERGIN_FOR_REPEAT 15
#define MERGIN_FOR_BED 16
#define WIDTH_REPEAT 5.5
#define WIDTH_BED 6

#define MEMNUM_GC 10
#define MEMNUM_PD 4

#define GAP_THRE 0.9

enum{
  LINE_CHIP,
  LINE_INPUT,
  LINE_RATIO,
  LINE_PVALUE_INTER,
  LINE_PVALUE_ENRICH
};

typedef struct{
  gdouble x1, x2, xcen, xwid;
  gdouble x_name;
  gint y_name, ylen;
  gint ybar;
  gint on_plus, on_minus;
} Dposi;

/* global variable */
gint yaxis_now;
gdouble dot_per_bp;

#define BP2XAXIS(v) ((v) * dot_per_bp + OFFSET_X)
#define rel_xline(cr, x1, y1, xlen) do{		\
    cairo_move_to(cr, x1, (int)y1);		\
    cairo_rel_line_to(cr, xlen, 0); }while(0)
#define rel_yline(cr, x1, y1, ylen) do{		\
    cairo_move_to(cr, x1, (int)y1);		\
    cairo_rel_line_to(cr, 0, (int)ylen); }while(0)

static void draw_page(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g, cairo_t *cr, int page_curr, int num_line, int num_page, int pstart, int pend, int region_no, int chr);
static void draw_dataline(DrParam *p, DDParam *d, SamplePair *sample, cairo_t *cr, int xstart, int xend);
static void stroke_xaxis(DDParam *d, cairo_t *cr, gint lstart, gint lend);
static void stroke_xaxis_num(DDParam *d, cairo_t *cr, gint lstart, gint lend, gint yaxis, gint fontsize);
static void stroke_ARS(DDParam *d, cairo_t *cr, gint xstart, gint xend);
static void define_posi(Dposi *p, gint start, gint end, gint dir, gchar *name, gint xstart, gint yaxis, gint dif, gint gene_posi_cnt);
static void draw_annotation(DDParam *d, cairo_t *cr, gint xstart, gint xend);
static void draw_graph(DDParam *d, cairo_t *cr, Graph *graph, gint xstart, gint xend, gint memnum, gint boxheight, bool color, bool xaxis);
static void draw_bedfile(DDParam *d, cairo_t *cr, gint xstart, gint xend, gint chr);
static void draw_repeat(DDParam *d, cairo_t *cr, gint xstart, gint xend);
static void print_pagetitle(cairo_t *cr, char *chrname, gint page_curr, gint num_page, gint region_no);
static void show_colors(cairo_t *cr, gint x, gint *ycen, gchar *label, gdouble r, gdouble g, gdouble b);
static void stroke_rectangle(cairo_t *cr, gint x1, gint x2, gint y1, gint height);
static void fill_rectangle(cairo_t *cr, gint x1, gint x2, gint y1, gint height, gdouble r, gdouble g, gdouble b, double alpha);
static void showtext_cr(cairo_t *cr, gdouble x, gdouble y, gchar *str, gint fontsize);

void merge_pdf(char *prefix){
  char *command = alloc_str_new(prefix, strlen(prefix) +50);
  /* pdftk */
  sprintf(command, "pdftk %s_*.pdf cat output %s.pdf", prefix, prefix);
  my_system(command);
  /* rm */
  sprintf(command, "rm %s_*.pdf", prefix);
  my_system(command);
  MYFREE(command);
  return;
}

int calc_pageheight(DrParam *p, DDParam *d){
  gint height=OFFSET_Y*2, height_sample=0, height_lpp=0, height_dataline=0;
  gint lineheight = d->ystep * d->barnum + MERGIN_BETWEEN_DATA;

  /* GC and GD */
  if(d->GC.argv) height += BOXHEIGHT_GRAPH + MERGIN_BETWEEN_GRAPH_DATA;
  if(d->GD.argv) height += BOXHEIGHT_GRAPH + MERGIN_BETWEEN_GRAPH_DATA;
  if(p->ftype==FTYPE_PD){  // PD
    height += BOXHEIGHT_PD *d->pdnum + MERGIN_BETWEEN_GRAPH_DATA *(d->pdnum-1);
  }else{
    /* annotation */
    if(d->gftype == GFTYPE_REFFLAT || d->gftype == GFTYPE_ENSEMBL) height_lpp += BOXHEIGHT_GENEBOX_EXON;
    else height_lpp += BOXHEIGHT_GENEBOX_NOEXON;
    height_lpp += MERGIN_BETWEEN_DATA;
    height_lpp += BOXHEIGHT_INTERACTION * d->internum;
    /* dataline */
    if(d->visualize_p_inter)  height_sample += lineheight;
    if(d->visualize_p_enrich) height_sample += lineheight;
    if(d->visualize_ratio)    height_sample += lineheight;
    if(d->visualize_ctag)     height_sample += lineheight;
    if(d->visualize_itag==1)  height_sample += lineheight;

    height_dataline += height_sample * p->samplenum;
    if(d->visualize_itag==2) height_dataline += lineheight;

    height_lpp += height_dataline;
    /* bedline */
    if(d->bednum){
      height_lpp += MERGIN_FOR_BED;
      height_lpp += LINEHEIGHT_BED * d->bednum;
    }
    if(d->repeat.argv) height_lpp += MERGIN_FOR_REPEAT + HEIGHT_REPEAT * (REPEATTYPENUM-3);
    height += height_lpp * d->linenum_per_page + MERGIN_BETWEEN_LINE * (d->linenum_per_page-1);
  }

  //  printf("%d %d %d %d %d\n", height, height_lpp, height_dataline, height_sample, lineheight);

  return height;
}

void draw_region(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g, int s, int e, char *prefix, int region_no, int chr){
  int i;
  int num_line = (e-s-1) / d->width_per_line +1;
  int num_page = (num_line -1) / d->linenum_per_page +1;
  char *format;
  if(!d->png) format="pdf"; else format="png";
  char *filename = alloc_str_new(prefix, 128);
  sprintf(filename, "%s.%s", prefix, format);
  remove_file(filename);

  /* initialize global variables */
  dot_per_bp = WIDTH_DRAW/(double)d->width_per_line;

  cairo_surface_t *surface;
  cairo_t *cr;
  for(i=0; i<num_page; i++){
    if(num_page > 1) sprintf(filename, "%s_%.4d_%.4d.%s", prefix, region_no +1, i+1, format);
    else             sprintf(filename, "%s_%.4d.%s",      prefix, region_no +1,      format);

    FLUSH("   page %5d/%5d/%5d\r", i+1, num_page, region_no+1);

    if(!d->png) surface = cairo_pdf_surface_create(filename, PAGEWIDTH, d->pageheight);
    else        surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, PAGEWIDTH, d->pageheight);
    cr = cairo_create(surface);
    draw_page(p, d, sample, g, cr, i, num_line, num_page, s, e, region_no, chr);
    cairo_destroy(cr);
    if(d->png) cairo_surface_write_to_png(surface, filename);
    cairo_surface_destroy(surface);
  }
  printf("\n");
  MYFREE(filename);

  char *tempstr=NULL;
  if(num_page > 1 && !d->png){
    tempstr = alloc_str_new(prefix, 20);
    sprintf(tempstr, "%s_%.4d", prefix, region_no +1);
    merge_pdf(tempstr);
    MYFREE(tempstr);
  }
  return;
}

static void draw_interaction(cairo_t *cr, InteractionSet *inter, gint xstart, gint xend, gint chr){
  gint i;
  gint boxheight = BOXHEIGHT_INTERACTION - 5;
  gint xcen_head, xcen_tail;
  gdouble xc, radius;
  gdouble r;
  Interaction *site;

  /* keys */
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  showtext_cr(cr, 40, yaxis_now + boxheight/2, inter->name, 12);

  yaxis_now += boxheight;
  cairo_set_line_width(cr, 2);
  for(i=0; i<inter->num; i++){
    site = &(inter->site[i]);
    if(site->head.chr != chr && site->tail.chr != chr) continue;
    else if(site->head.chr != chr || site->tail.chr != chr) cairo_set_source_rgba(cr, CLR_GREEN4, 0.8);  // inter-chromosomal
    else cairo_set_source_rgba(cr, CLR_RED, 0.8);  // intra-chromosomal
    xcen_head = xcen_tail = -1;
    if(overlap(site->head.start, site->head.end, xstart, xend)) xcen_head = (site->head.start + site->head.end)/2 - xstart;
    if(overlap(site->tail.start, site->tail.end, xstart, xend)) xcen_tail = (site->tail.start + site->tail.end)/2 - xstart;
    if(site->head.chr == chr && site->tail.chr == chr && xcen_head >= 0 && xcen_tail >= 0){
      xc = BP2XAXIS((xcen_head+xcen_tail)/2);
      radius = (xcen_tail - xcen_head)/2 * dot_per_bp;
      r = max(radius/boxheight, 1);
      cairo_scale(cr,1,1/r);
      cairo_arc(cr, xc, yaxis_now*r, radius, M_PI, 2*M_PI);
      cairo_stroke(cr);
      cairo_scale(cr,1,r);
    }
    if((site->head.chr == chr && xcen_head > 0) && (site->tail.chr != chr || xcen_tail < 0)){
      xc = BP2XAXIS(xend - xstart);
      radius = (xend - xcen_head - xstart) * dot_per_bp;
      r = max(radius/boxheight, 1);
      cairo_scale(cr,1,1/r);
      cairo_arc(cr, xc, yaxis_now*r, radius, M_PI, 1.5*M_PI);
      cairo_stroke(cr);
      cairo_scale(cr,1,r);
    }
    if((site->head.chr != chr || xcen_head < 0) && (site->tail.chr == chr && xcen_tail > 0)){
      xc = OFFSET_X;
      radius = xcen_tail * dot_per_bp;
      r = max(radius/boxheight, 1);
      cairo_scale(cr,1,1/r);
      cairo_arc(cr, xc, yaxis_now*r, radius, 1.5*M_PI, 2*M_PI);
      cairo_stroke(cr);
      cairo_scale(cr,1,r);
    }
  }
  cairo_stroke(cr);
  yaxis_now += 5;
  return;
}

static void draw_page(DrParam *p, DDParam *d, SamplePair *sample, RefGenome *g, cairo_t *cr, int page_curr, int num_line, int num_page, int pstart, int pend, int region_no, int chr){
  int i,j;
  int xstart, xend;
  int lstart, lend;   /* no of start and end lines for this page */
  lstart = page_curr * d->linenum_per_page;
  if(page_curr == num_page-1) lend = num_line;
  else lend = (page_curr+1) * d->linenum_per_page;

  yaxis_now = OFFSET_Y;  // initialize yaxis_now
  if(d->backcolors){
    cairo_set_source_rgb(cr, CLR_WHITE);
    cairo_paint(cr);
  }

  for(i=lstart; i<lend; i++){
    xstart = pstart + i * d->width_per_line;
    if(i==num_line-1) xend = pend;
    else xend = pstart + (i+1) * d->width_per_line -1;
    if(xstart >= xend) continue;

    if(d->GC.argv){
      draw_graph(d, cr, &(d->GC), xstart, xend, MEMNUM_GC, BOXHEIGHT_GRAPH, false, true);
      yaxis_now += MERGIN_BETWEEN_GRAPH_DATA;
    }
    if(d->GD.argv){
      draw_graph(d, cr, &(d->GD), xstart, xend, MEMNUM_GC, BOXHEIGHT_GRAPH, false, true);
      yaxis_now += MERGIN_BETWEEN_GRAPH_DATA;
    }

    if(p->ftype==FTYPE_PD){
      for(j=0; j<d->pdnum; j++){
	draw_graph(d, cr, &(d->PD[j]), xstart, xend, MEMNUM_PD, BOXHEIGHT_PD, true, false);
	if(j!=d->pdnum-1) yaxis_now += MERGIN_BETWEEN_GRAPH_DATA;
      }
      stroke_xaxis_num(d, cr, xstart, xend, yaxis_now, 9);
    }else{
      if(d->gene.argv || d->arsfile || d->terfile) draw_annotation(d, cr, xstart, xend);
      if(d->internum){
	for(j=0; j<d->internum; j++) draw_interaction(cr, &(d->inter[j]), xstart, xend, chr);
      }
      draw_dataline(p, d, sample, cr, xstart, xend);
      if(d->bednum) draw_bedfile(d, cr, xstart, xend, chr);
      if(d->repeat.argv) draw_repeat(d, cr, xstart, xend);
    }
    yaxis_now += MERGIN_BETWEEN_LINE;
  }
  print_pagetitle(cr, g->chr[chr].name, page_curr, num_page, region_no);
  return;
}

static void peakcall(DrParam *p, cairo_t *cr, SamplePair *sample, gint i){
  double ratio=0, pval_inter=0, pval_enr=0;
  calc_ratio_and_pvalues(&ratio, &pval_inter, &pval_enr, sample, i);
  if(judge_significant(p, pval_enr, pval_inter, ratio, WIGARRAY2VALUE(sample->ChIP->data[i]), sample->Input->argv)) cairo_set_source_rgba(cr, CLR_PINK, 1);
  return;
}

static gdouble define_value_and_color(DrParam *p, DDParam *d, cairo_t *cr, SamplePair *sample, gint i, LINE_Type type){
  gdouble value=0, ratio, data;
  switch(type){
  case LTYPE_CHIP:
    value = WIGARRAY2VALUE(sample->ChIP->data[i]) /sample->scale_tag;

    cairo_set_source_rgba(cr, CLR_GREEN3, 1);

    if(sample->peakarray && sample->peakarray[i]) cairo_set_source_rgba(cr, CLR_PINK, 1);
    if(!sample->peakarray) peakcall(p, cr, sample, i);
    
    break;
  case LTYPE_INPUT:
    value = WIGARRAY2VALUE(sample->Input->data[i]) /sample->scale_tag;
    cairo_set_source_rgba(cr, CLR_BLUE, 1); 
    break;
  case LTYPE_RATIO_GV: /* same as LTYPE_RATIO */
  case LTYPE_RATIO:
    ratio = CALCRATIO(sample->ChIP->data[i], sample->Input->data[i], sample->comp->genome->ratio);
    if(!ratio) value = 0;
    else{
      if(d->visualize_ratio==2) value = log2(ratio) / sample->scale_ratio;
      else value = ratio / sample->scale_ratio; 
    }
    if(type==LTYPE_RATIO_GV){
      if(d->visualize_ratio==1 && value > 1) cairo_set_source_rgba(cr, CLR_PINK, 1);
      else if(d->visualize_ratio==2 && value > 0) cairo_set_source_rgba(cr, CLR_PINK, 1);
      else cairo_set_source_rgba(cr, CLR_GRAY3, 1);
      break;
    }
    if(p->ftype!=FTYPE_PEAKCALL_E){
      if(ratio > p->enrichthre) cairo_set_source_rgba(cr, CLR_PINK, 1);
      else cairo_set_source_rgba(cr, CLR_GRAY, 1);
    }else cairo_set_source_rgba(cr, CLR_ORANGE, 1);
    if(d->visualize_ratio==2 && value < 0) cairo_set_source_rgba(cr, CLR_GRAY3, 1);
    break;
  case LTYPE_PVALUE_INTER:
    data = zero_inflated_binomial_test(WIGARRAY2VALUE(sample->ChIP->data[i]), sample->ChIP->nb_p, sample->ChIP->nb_n);
    value = data /sample->scale_pvalue;
    if(data > p->pthre_internal) cairo_set_source_rgba(cr, CLR_PINK, 1);
    else cairo_set_source_rgba(cr, CLR_GRAY, 1);
    break;
  case LTYPE_PVALUE_ENRICH:
    data = binomial_test(WIGARRAY2VALUE(sample->ChIP->data[i]), WIGARRAY2VALUE(sample->Input->data[i]), sample->comp->genome->ratio);
    value = data /sample->scale_pvalue;
    if(data > p->pthre_enrich) cairo_set_source_rgba(cr, CLR_PINK, 1);
    else cairo_set_source_rgba(cr, CLR_GRAY, 1);
    break;
  default:
    fprintf(stderr, "error: invalid stroke value.\n");
    exit(1);
  }
  return value;
}

static void stroke_bindata(DrParam *p, DDParam *d, cairo_t *cr, SamplePair *sample, gint xstart, gint xend, LINE_Type type){
  gint i,binsize = sample->binsize;
  gint sbin = xstart/binsize;
  gint ebin = xend/binsize;
  gdouble xcen, len, value=0;
  gint thin = min(d->width_per_line/(1000*binsize), 20);
  gint barnum_minus = d->barnum/2;
  gint barnum_plus = d->barnum - barnum_minus;
  gdouble diff = binsize * dot_per_bp;
 
  xcen = BP2XAXIS(binsize/2);  // initial position
  if(thin > 1) cairo_set_line_width(cr, (binsize*dot_per_bp*thin));
  else cairo_set_line_width(cr, (binsize*dot_per_bp));

  for(i=sbin; i<ebin; i++, xcen += diff){
    if(thin > 1 && i%thin) continue;

    if(p->gapfile && d->gaparray[i] >= GAP_THRE){
      cairo_set_source_rgba(cr, CLR_BLACK, 0.3);
      rel_yline(cr, xcen, (gint)yaxis_now - d->ystep * d->barnum, (gint)(d->ystep * d->barnum));
      cairo_stroke(cr);
    }else if(p->mpfile && d->mparray[i] < p->mpthre){
      cairo_set_source_rgba(cr, CLR_PURPLE, 0.3);
      rel_yline(cr, xcen, (gint)yaxis_now - d->ystep * d->barnum, (gint)(d->ystep * d->barnum));
      cairo_stroke(cr);
    }
    value = define_value_and_color(p, d, cr, sample, i, type);
    if(!value) continue;

    if(d->visualize_ratio==2){
      if(value>0) len = -d->ystep * (min(value, barnum_plus));
      else        len = -d->ystep * (max(value, -barnum_minus));
      rel_yline(cr, xcen, (gint)yaxis_now - d->ystep * barnum_minus, (gint)len);
    }else{
      len = -d->ystep * (min(value, d->barnum));
      rel_yline(cr, xcen, (gint)yaxis_now, (gint)len);
    }
    cairo_stroke(cr);
  }

  /*  if(p_opt->localmax){ // localmax
    cairo_set_line_width(cr, 1.0);
    cairo_set_source_rgba(cr, CLR_BLUE, 1);
    for(j=0; j<p->ChIP.peak.num; j++){
      if(p->ChIP.peak.bs[j].chr != chr) continue;
      if(p->ChIP.peak.bs[j].maxposi*binsize>= xstart && p->ChIP.peak.bs[j].maxposi*binsize < xend){
	xcen = BP2XAXIS(p->ChIP.peak.bs[j].maxposi*binsize - xstart +1);
	rel_yline(cr, xcen, d->yaxis - d->ystep * d->barnum, d->ystep * d->barnum);
	cairo_stroke(cr);
      }
    }
    }*/

  return;
}

static void stroke_dataframe(DDParam *d, cairo_t *cr, SamplePair *sample, gint xstart, gint xend, LINE_Type type){
  gint i;
  gdouble width = (xend-xstart+1) * dot_per_bp;
  gdouble height = d->ystep * d->barnum;
  gdouble y;
  gchar str[128];

  /* frame */
  cairo_set_line_width(cr, 0.4);
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  rel_xline(cr, OFFSET_X, yaxis_now, width);
  rel_yline(cr, OFFSET_X, yaxis_now-height, height);
  cairo_stroke(cr);
  if(d->backcolors) fill_rectangle(cr, xstart, xend, yaxis_now-height, height, CLR_BLUE3, 0.1);
  
  /* y memory */
  if(d->backcolors){
    cairo_set_source_rgba(cr, CLR_BLACK, 0.5);
    for(i=0, y=yaxis_now; i<d->barnum; i++, y-=d->ystep) rel_xline(cr, OFFSET_X, y, width);
    cairo_stroke(cr);
  }

  /* y key */
  gint barnum_minus = d->barnum/2;
  gdouble x, scale=0;
  if(type==LTYPE_CHIP || type==LTYPE_INPUT) scale = sample->scale_tag;
  else if(type==LTYPE_RATIO || type==LTYPE_RATIO_GV) scale = sample->scale_ratio;
  else if(type==LTYPE_PVALUE_ENRICH || type==LTYPE_PVALUE_INTER) scale = sample->scale_pvalue;
  else{ fprintf(stderr, "error:%s: invalid type.\n", __FUNCTION__); exit(1);}

  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  x = OFFSET_X + width + 7;
  for(i=0; i<=d->barnum; i++){
    if(d->visualize_ratio==2){
      if(i < barnum_minus) sprintf(str, "1/%d", (int)pow(2, (barnum_minus-i)*scale));
      else sprintf(str, "%d", (int)pow(2, (i - barnum_minus)*scale));
    }else  sprintf(str, "%.1f", (gdouble)(i*scale));
    y = yaxis_now - i*(d->ystep - 1.5);
    showtext_cr(cr, x, y, str, 9);
  }

  /* keys */
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  switch(type){
  case LTYPE_CHIP: sprintf(str, "%s", sample->linename); break;
  case LTYPE_INPUT: sprintf(str, "Input"); break;
  case LTYPE_RATIO_GV: /* same as LTYPE_RATIO */ 
  case LTYPE_RATIO: 
    if(d->visualize_ctag) sprintf(str, "IP/Input");
    else sprintf(str, "%s", sample->linename);
    break;
  case LTYPE_PVALUE_INTER:
    if(d->visualize_ctag || d->visualize_ratio) sprintf(str, "pval (ChIP internal)");
    else sprintf(str, "%s", sample->linename);
    break;
  case LTYPE_PVALUE_ENRICH: 
    if(d->visualize_ctag || d->visualize_ratio) sprintf(str, "pval (IP/Input)");
    else sprintf(str, "%s", sample->linename);
    break;
  }
  x = 40;
  y = yaxis_now - d->ystep * d->barnum/2; 
  showtext_cr(cr, x, y, str, 12);

  return;
}

static void stroke_readdist(DrParam *p, DDParam *d, cairo_t *cr, SamplePair *sample, gint xstart, gint xend, LINE_Type type){
  //  printf("type=%d\n", type);
  yaxis_now += d->ystep * d->barnum + MERGIN_BETWEEN_DATA;
  stroke_bindata(p, d, cr, sample, xstart, xend, type);
  stroke_dataframe(d, cr, sample, xstart, xend, type);
  stroke_xaxis(d, cr, xstart, xend);
  return;
}

static void draw_dataline(DrParam *p, DDParam *d, SamplePair *sample, cairo_t *cr, int xstart, int xend){
  int i;
  for(i=0; i<p->samplenum; i++){
    if(d->visualize_p_inter)  stroke_readdist(p, d, cr, &(sample[i]), xstart, xend, LTYPE_PVALUE_INTER);
    if(d->visualize_p_enrich) stroke_readdist(p, d, cr, &(sample[i]), xstart, xend, LTYPE_PVALUE_ENRICH);
    if(d->visualize_ratio){
      if(p->ftype==FTYPE_GV)  stroke_readdist(p, d, cr, &(sample[i]), xstart, xend, LTYPE_RATIO_GV);
      else                    stroke_readdist(p, d, cr, &(sample[i]), xstart, xend, LTYPE_RATIO);
    }if(d->visualize_ctag)    stroke_readdist(p, d, cr, &(sample[i]), xstart, xend, LTYPE_CHIP);
    if(d->visualize_itag==1)  stroke_readdist(p, d, cr, &(sample[i]), xstart, xend, LTYPE_INPUT);
  }
  if(d->visualize_itag==2) stroke_readdist(p, d, cr, &(sample[0]), xstart, xend, LTYPE_INPUT);

  stroke_xaxis_num(d, cr, xstart, xend, yaxis_now, 9);
  return;
}

static gint setline(gint lstart, gint interval){
  gint posi = lstart-1;
  if(!posi%interval) return posi;
  posi = posi/interval +1;
  posi *= interval;
  return posi;
}

static gint set_interval_large(DDParam *d){
  gint interval;
  if(d->gw==true){
    if(d->large_genome==false) interval = 100*NUM_1K;  // 100kbp
    else interval = 10*NUM_1M;  // 10Mbp
  }else{
    interval = d->width_per_line/10;
  }
  return interval;
}

static void stroke_xaxis(DDParam *d, cairo_t *cr, gint lstart, gint lend){
  gint i;
  gdouble x;
  gint interval_large = set_interval_large(d);
  gint interval = interval_large/10;

  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  for(i=setline(lstart, interval); i<=lend; i+=interval){
    x = BP2XAXIS(i-lstart);
    if(!(i%interval_large)){
      cairo_set_line_width(cr, 1);
      rel_yline(cr, x, yaxis_now-4, 8);
    }else{
      cairo_set_line_width(cr, 0.5);
      rel_yline(cr, x, yaxis_now-1.5, 3);
    }
    cairo_stroke(cr);
  }
  return;
}

static void stroke_xaxis_num(DDParam *d, cairo_t *cr, gint lstart, gint lend, gint yaxis, gint fontsize){
  gint i, mega, kilo;
  gdouble x;
  gchar str[128];
  gint interval_large = set_interval_large(d);

  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  for(i=setline(lstart, interval_large); i<=lend; i+=interval_large){
    x = BP2XAXIS(i-lstart);
    if(d->gw==true){
      if(d->large_genome==false) sprintf(str, "%dk", i/NUM_1K);
      else sprintf(str, "%dM", i/NUM_1M);
      showtext_cr(cr, x - 3*strlen(str), yaxis+10, str, fontsize);
    }else{
      mega = i/NUM_1M;
      kilo = (i%NUM_1M)/NUM_1K;
      if(d->width_per_line > 10*NUM_1K) sprintf(str, "%.3fM", i/(double)NUM_1M);
      else if(d->width_per_line > 10){
	if(mega) sprintf(str, "%d,%.3fK", mega, (i%NUM_1M)/(double)NUM_1K);
	else sprintf(str, "%.3fK", i/(double)NUM_1K);
      }else{
	if(mega) sprintf(str, "%d,%d,%d", mega, kilo, i%NUM_1K);
	else if(kilo) sprintf(str, "%d,%d", kilo, i%NUM_1K);
	else sprintf(str, "%d", i);
      }
      showtext_cr(cr, x - 3*strlen(str), yaxis+10, str, fontsize);
    }
  }
  return;
}

static void stroke_ARS(DDParam *d, cairo_t *cr, gint xstart, gint xend){
  gint i, ars_on=0;
  Gene *gene = d->gene.gene;
  Dposi *posi = (Dposi *)my_calloc(1, sizeof(Dposi), "define_posi");

  cairo_set_line_width(cr, 0.3);
  for(i=0; i<d->gene.num; i++){
    if(!overlap(gene[i].start, gene[i].end, xstart, xend)) continue;

    define_posi(posi, gene[i].start, gene[i].end, gene[i].dir, gene[i].name, xstart, yaxis_now, 6, 1);
    if(strstr(gene[i].name, "ARS")){
      cairo_set_source_rgba(cr, CLR_RED, 1);
      rel_yline(cr, posi->xcen, yaxis_now -2, posi->ylen +14 - 8 * ars_on);
      showtext_cr(cr, posi->x_name, posi->y_name +8 - 8*ars_on, gene[i].name, 8);
      if(ars_on==2) ars_on=0; else ars_on++;
    }else if(strstr(gene[i].name, "CEN")){
      cairo_set_source_rgba(cr, CLR_GREEN, 1);
      rel_yline(cr, posi->xcen, yaxis_now -2, posi->ylen -2);
      showtext_cr(cr, posi->x_name, posi->y_name -8, gene[i].name, 8);
      //      rel_yline(cr, posi->xcen, yaxis_now -2, posi->ylen -14);
      // showtext_cr(cr, posi->x_name, posi->y_name-20, gene[i].name, 8);
    }else if(strstr(gene[i].name, "TER")){
      cairo_set_source_rgba(cr, CLR_OLIVE, 1);
      rel_yline(cr, posi->xcen, yaxis_now -2, posi->ylen -5);
      showtext_cr(cr, posi->x_name, posi->y_name-11, gene[i].name, 7);
    }else continue;
    gene[i].delete=1;
  }
  
  MYFREE(posi);
  return;
}

static void def_xposi(Dposi *p, gint start, gint end, gint xstart, gchar *name){
  p->x1 = BP2XAXIS(start - xstart +1);
  p->x2 = BP2XAXIS(end   - xstart +1);
  p->xcen = (p->x1 + p->x2) /2;
  p->xwid = p->x2 - p->x1;

  p->x_name = p->xcen - 3.25*strlen(name) + 6;
  if(p->x_name < 0) p->x_name = 10;
  else if(p->x_name > PAGEWIDTH) p->x_name = PAGEWIDTH - 50;

  return;
}

static void define_posi(Dposi *p, gint start, gint end, gint dir, gchar *name, gint xstart, gint yaxis, gint dif, gint gene_posi_cnt){
  def_xposi(p, start, end, xstart, name);

  p->ybar = yaxis - dif*dir;
  if(dir==1){
    p->y_name = p->ybar -5 - p->on_minus*6;
    if(p->on_minus==gene_posi_cnt){ p->on_minus =0;}
    else{ (p->on_minus)++;}
  }else if(dir == -1){
    p->y_name = p->ybar +9 + p->on_plus*6;
    if(p->on_plus==gene_posi_cnt){ p->on_plus=0;}
    else{ (p->on_plus)++;}
  }else{
    p->y_name = yaxis -22;
  }
  p->ylen = p->y_name - yaxis_now;
  return;
}

static void show_geneanno(DDParam *d, cairo_t *cr){
  gint x = 50;
  gint y = yaxis_now-20;

  cairo_set_line_width(cr, 2.5);
  show_colors(cr, x, &y, "coding", CLR_BLUE);
  show_colors(cr, x, &y, "noncoding", CLR_GREEN);
  if(d->gftype == GFTYPE_ENSEMBL){
    show_colors(cr, x, &y, "processed transcript", CLR_ORANGE);
    show_colors(cr, x, &y, "microRNA", CLR_PINK);
    show_colors(cr, x, &y, "pseudo", CLR_GRAY2);
    show_colors(cr, x, &y, "others", CLR_BLACK);
  }
  else if(d->gftype == GFTYPE_SGD){
    show_colors(cr, x, &y, "rRNA", CLR_BLACK);
    show_colors(cr, x, &y, "LTR", CLR_PURPLE);
  }
  return;
}

static void stroke_gene_yeast(DDParam *d, cairo_t *cr, gint xstart, gint xend){
  gint i;
  Gene  *gene = d->gene.gene;
  Dposi *posi = (Dposi *)my_calloc(1, sizeof(Dposi), "posi");

  show_geneanno(d, cr);

  cairo_set_line_width(cr, 1.5);
  for(i=0; i<d->gene.num; i++){
    if(gene[i].delete) continue;
    if(!overlap(gene[i].start, gene[i].end, xstart, xend)) continue;

    define_posi(posi, gene[i].start, gene[i].end, gene[i].dir, gene[i].name, xstart, yaxis_now, 6, 1);
    rel_xline(cr, posi->x1, posi->ybar, posi->xwid);
    cairo_stroke(cr);

    switch(gene[i].genetype){
    case CENTROMERE:
      cairo_set_source_rgba(cr, CLR_GREEN, 1);
      rel_yline(cr, posi->xcen, yaxis_now -2, posi->ylen);
      showtext_cr(cr, posi->x_name, posi->y_name-6, gene[i].name, 8);
      break;
    case ARS:
      cairo_set_source_rgba(cr, CLR_RED, 1);
      rel_yline(cr, posi->xcen, yaxis_now -2, posi->ylen +2);
      showtext_cr(cr, posi->x_name, posi->y_name-4, gene[i].name, 7);
      break;
    case TER:
      cairo_set_source_rgba(cr, CLR_OLIVE, 1);
      rel_yline(cr, posi->xcen, yaxis_now -2, posi->ylen -5);
      showtext_cr(cr, posi->x_name, posi->y_name-11, gene[i].name, 7);
      break;
    case rRNA:
      cairo_set_source_rgba(cr, CLR_BLACK, 1);
      showtext_cr(cr, posi->x_name, posi->y_name, gene[i].name, 6);
      break;
    case LTR:
      cairo_set_source_rgba(cr, CLR_PURPLE, 1);
      showtext_cr(cr, posi->x_name, posi->y_name, gene[i].name, 6);
      break;
    case NONCODING:
      cairo_set_source_rgba(cr, CLR_GREEN, 1);
      showtext_cr(cr, posi->x_name, posi->y_name, gene[i].name, 6);
      break;
    default:
      cairo_set_source_rgba(cr, CLR_BLUE, 1);
      showtext_cr(cr, posi->x_name, posi->y_name, gene[i].name, 6);
      break;
    }
  }
  MYFREE(posi);
  return;
}

static void stroke_gene_vertebrate(DDParam *d, cairo_t *cr, gint xstart, gint xend){
  gint i, j, on_plus=0, on_minus=0;
  gdouble x, xlen;
  gint ybar=0, y_name=0;
  Gene *gene = d->gene.gene;
  Dposi *posi = (Dposi *)my_calloc(1, sizeof(Dposi), "posi");

  show_geneanno(d, cr);

  for(i=0; i<d->gene.num; i++){
    if(gene[i].delete) continue;
    if(!overlap(gene[i].start, gene[i].end, xstart, xend)) continue;

    if(gene[i].dir == 1){  // + strand
      ybar = yaxis_now -6 - on_minus * 13;
      y_name = ybar -6;
      if(on_minus==3) on_minus=0; else on_minus++;
    }else{                 // - strand
      ybar = yaxis_now +12 + on_plus * 13;
      y_name = ybar +11;
      if(on_plus==3) on_plus=0; else on_plus++;
    }

    if(gene[i].genetype==CODING)         cairo_set_source_rgba(cr, CLR_BLUE, 1);
    else if(gene[i].genetype==NONCODING) cairo_set_source_rgba(cr, CLR_GREEN, 1);
    else if(gene[i].genetype==MIRNA)     cairo_set_source_rgba(cr, CLR_PINK, 1);
    else if(gene[i].genetype==PROCESS)   cairo_set_source_rgba(cr, CLR_ORANGE, 1);
    else if(gene[i].genetype==PSEUDO)    cairo_set_source_rgba(cr, CLR_GRAY2, 1);
    else                                 cairo_set_source_rgba(cr, CLR_BLACK, 1);

    def_xposi(posi, gene[i].start, gene[i].end, xstart, gene[i].name);
    /* gene body */
    cairo_set_line_width(cr, 1.5);
    rel_yline(cr, posi->x1, ybar-4, 8);
    rel_yline(cr, posi->x2, ybar-4, 8);
    cairo_stroke(cr);
    cairo_set_line_width(cr, 3);
    rel_xline(cr, posi->x1, ybar, posi->xwid);
    cairo_stroke(cr);
    /* exon */
    cairo_set_line_width(cr, 6);
    for(j=0; j<gene[i].exonnum; j++){
      x = BP2XAXIS(gene[i].exon[j].start - xstart +1);
      xlen = (gene[i].exon[j].end - gene[i].exon[j].start) * dot_per_bp;
      rel_xline(cr, x, ybar, xlen);
      cairo_stroke(cr);
    }
    /* name */
    cairo_set_source_rgba(cr, CLR_BLACK, 1);
    showtext_cr(cr, posi->x_name, y_name, gene[i].name, 10);
  }
  MYFREE(posi);
  return;
}

static void draw_annotation(DDParam *d, cairo_t *cr, gint xstart, gint xend){
  int boxheight;
  if(d->gftype == GFTYPE_REFFLAT || d->gftype == GFTYPE_ENSEMBL) boxheight = BOXHEIGHT_GENEBOX_EXON;
  else boxheight = BOXHEIGHT_GENEBOX_NOEXON;

  gdouble y1 = yaxis_now;
  yaxis_now += boxheight/2;

  if(d->show_ars){
    stroke_ARS(d, cr, xstart, xend);
    showtext_cr(cr, 70, yaxis_now, "ARS", 12);
  }else{
    if(d->arsfile) stroke_ARS(d, cr, xstart, xend);
    if(d->gftype == GFTYPE_SGD || d->gftype == GFTYPE_GTF) stroke_gene_yeast(d, cr, xstart, xend);
    else stroke_gene_vertebrate(d, cr, xstart, xend);
  }
  /* frame */
  stroke_rectangle(cr, xstart, xend, y1, boxheight);
  if(d->backcolors) fill_rectangle(cr, xstart, xend, y1, boxheight, CLR_YELLOW2, 0.1);

  /* genome line */
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_set_line_width(cr, 1.5);
  rel_xline(cr, OFFSET_X, yaxis_now, (xend - xstart +1) * dot_per_bp);
  cairo_stroke(cr);

  /* memory */
  stroke_xaxis(d, cr, xstart, xend);
  stroke_xaxis_num(d, cr, xstart, xend, yaxis_now, 7);

  yaxis_now += boxheight/2 + MERGIN_BETWEEN_DATA;
  return;
}

static gdouble defy_graph(gdouble value, gint boxheight, gdouble p_low, gdouble p_high){
  gdouble y;
  y = yaxis_now + boxheight/2 - boxheight*(value - p_low)/(p_high - p_low);
  return y;
}

static void check_bottom(cairo_t *cr, gdouble x_pre, gdouble y_pre, gdouble x_cen, gdouble y_cen, gint bottom){
  if(y_pre > bottom && y_cen > bottom) return;

  gdouble x1=x_pre, x2=x_cen, y1=y_pre, y2=y_cen;
  gdouble xlen, ylen, ydiff;
  if(y_pre > bottom || y_cen > bottom){
    xlen = abs(x2-x1);
    ylen = abs(y2-y1);
    if(y_pre > bottom){
      ydiff = abs(y_cen-bottom);
      x1 = x2 - xlen*(ydiff/ylen);
      y1 = bottom;
    }else{
      ydiff = abs(y_pre-bottom);
      x2 = x1 - xlen*(ydiff/ylen);
      y2 = bottom;
    }
  }
  cairo_move_to(cr, x1, y1);
  cairo_line_to(cr, x2, y2);
  cairo_stroke(cr);
  return;
}

static void draw_graph(DDParam *d, cairo_t *cr, Graph *graph, gint xstart, gint xend, gint memnum, gint boxheight, bool color, bool xaxis){
  gint i;
  gchar str[16];
  gint s = xstart/graph->wsize, e = xend/graph->wsize +1;
  gdouble xcen, xpre, ycen, ypre;
  gdouble width = (xend-xstart+1) * dot_per_bp;
  gdouble diff = graph->wsize * dot_per_bp;

  yaxis_now += boxheight/2;

  /* graph line */
  if(color==false) cairo_set_source_rgba(cr, CLR_GREEN, 1); else cairo_set_source_rgba(cr, CLR_BLUE, 1);
  cairo_set_line_width(cr, 0.6);
  xpre = BP2XAXIS(0);
  xcen = BP2XAXIS(0.5*graph->wsize);
  ypre = defy_graph(graph->array[s], boxheight, graph->mmin, graph->mmax);
  for(i=s; i<e; i++, xcen += diff){
    ycen = defy_graph(graph->array[i], boxheight, graph->mmin, graph->mmax);
    check_bottom(cr, xpre, ypre, xcen, ycen, yaxis_now + boxheight/2 + 10);
    xpre = xcen;
    ypre = ycen;
  }
  /* keys */
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  showtext_cr(cr, OFFSET_X - 5*strlen(graph->name)-55, yaxis_now, graph->name, 13);

  yaxis_now += boxheight/2;

  /* axis */
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_set_line_width(cr, 0.4);
  rel_yline(cr, OFFSET_X, yaxis_now - boxheight, boxheight);
  cairo_stroke(cr);
  cairo_set_line_width(cr, 1.5);
  rel_xline(cr, OFFSET_X, yaxis_now, width);
  cairo_stroke(cr);
  stroke_xaxis(d, cr, xstart, xend);

  /* memory */
  gdouble x, y;
  gdouble mem = (graph->mmax - graph->mmin)/memnum;
  //  gdouble memnum = (graph->mmax - graph->mmin)/graph->mem;
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_set_line_width(cr, 0.5);
  for(i=0; i<=memnum; i++){
    if(mem <1) sprintf(str, "%.2f", graph->mmin + i*mem);
    else       sprintf(str, "%d", (int)(graph->mmin + i*mem));
    x = OFFSET_X - 5*strlen(str) - 7;
    y = yaxis_now - i*boxheight/memnum;
    showtext_cr(cr, x, y+2, str, 9);
    rel_xline(cr, OFFSET_X-2, y, 2);
    cairo_stroke(cr);
  }
  if(xaxis==true) stroke_xaxis_num(d, cr, xstart, xend, yaxis_now, 9);

  return;
}

static void draw_bedfile(DDParam *d, cairo_t *cr, gint xstart, gint xend, gint chr){
  gint i,j;
  gdouble x1, xlen;
  gint lineheight_half = LINEHEIGHT_BED/2;
  BedChr *bed=NULL;
  
  yaxis_now += MERGIN_FOR_BED;

  cairo_set_line_width(cr, WIDTH_BED);
  for(j=0; j<d->bednum; j++){
    bed = &(d->bed[j]->chr[chr]);
    yaxis_now += lineheight_half;

    /* color */
    if(j%2) cairo_set_source_rgba(cr, CLR_GRAY4, 1);
    else cairo_set_source_rgba(cr, CLR_GREEN, 1);
    for(i=0; i<bed->num; i++){
      if(overlap(bed->bed[i].s, bed->bed[i].e, xstart, xend)){
	if(bed->bed[i].itemRgb.r != -1){
	  // color
	  cairo_set_source_rgba(cr, bed->bed[i].itemRgb.r/(double)255, bed->bed[i].itemRgb.g/(double)255, bed->bed[i].itemRgb.b/(double)255, 0.6);
	}
	x1 = BP2XAXIS(bed->bed[i].s - xstart +1);
	xlen = (bed->bed[i].e - bed->bed[i].s) * dot_per_bp;
	rel_xline(cr, x1, yaxis_now, xlen);
	cairo_stroke(cr);
	if(bed->bed[i].itemRgb.r != -1){
	  // name
	  cairo_set_source_rgba(cr, CLR_BLACK, 1);
	  showtext_cr(cr, x1 + xlen/2, yaxis_now +3, bed->bed[i].name, 7);
	}
      }
    }
    cairo_set_source_rgba(cr, CLR_BLACK, 1);
    showtext_cr(cr, 80, yaxis_now + lineheight_half, d->bed[j]->name, 11);
    
    yaxis_now += lineheight_half;
  }

  /* frame */
  gint boxheight = LINEHEIGHT_BED *d->bednum;
  stroke_rectangle(cr, xstart, xend, yaxis_now - boxheight, boxheight);
  stroke_xaxis(d, cr, xstart, xend);
  stroke_xaxis_num(d, cr, xstart, xend, yaxis_now, 9);

  return;
}

static void stroke_repeat(cairo_t *cr, RepeatSet *rset, gint xstart, gint xend, gint height, RepeatType type, gchar *label){
  gint i;
  gint y1 = yaxis_now, ycen = y1 + height/2;
  Dposi *posi = (Dposi *)my_calloc(1, sizeof(Dposi), "define_posi");
  Repeat *repeat = rset->repeat;
  gint num = rset->num;

  for(i=0; i<num; i++){
    if(repeat[i].rtype != type) continue;
    if(!overlap(repeat[i].start, repeat[i].end, xstart, xend)) continue;
    /* color */
    if(type%2) cairo_set_source_rgba(cr, CLR_GRAY4, 1);
    else cairo_set_source_rgba(cr, CLR_PURPLE2, 1);

    if(!strcmp(repeat[i].class, "centr")) cairo_set_source_rgba(cr, CLR_GREEN, 1);
    else if(!strcmp(repeat[i].class, "telo")) cairo_set_source_rgba(cr, CLR_BLUE, 1);
    else if(!strcmp(repeat[i].class, "acro")) cairo_set_source_rgba(cr, CLR_RED, 1);
    else if(!strcmp(repeat[i].name, "HSATII")) cairo_set_source_rgba(cr, CLR_YELLOW, 1);
    else if(!strcmp(repeat[i].name, "(CATTC)n")) cairo_set_source_rgba(cr, CLR_YELLOW2, 1);
    else if(!strcmp(repeat[i].name, "(GAATG)n")) cairo_set_source_rgba(cr, CLR_OLIVE, 1);
    if(!strcmp(repeat[i].name, "ALR/Alpha")) cairo_set_source_rgba(cr, CLR_ORANGE, 1);
    
    define_posi(posi, repeat[i].start, repeat[i].end, repeat[i].dir, repeat[i].name, xstart, yaxis_now, 6, 1);
    cairo_set_line_width(cr, WIDTH_REPEAT);
    rel_xline(cr, posi->x1, ycen, posi->xwid);
    cairo_stroke(cr);
  }
  /* keys */
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  showtext_cr(cr, OFFSET_X - 80, ycen+4, label, 10);
  MYFREE(posi);

  yaxis_now += height;
  return;
}

static void draw_repeat(DDParam *d, cairo_t *cr, gint xstart, gint xend){
  char *repeat_type_str[]={"RepeatMasker", "SimpleRepeats", "Microsatellite", "SINE", "LINE", "LTR", "DNA", "Simple", "Low_Complexity", "Satellite", "RNA", "Other"};

  yaxis_now += MERGIN_FOR_REPEAT;
  gint y1 = yaxis_now;

  gint lineheight = HEIGHT_REPEAT;
  RepeatType type;
  gint boxheight = 0;
  for(type=RM_SINE; type<=RM_Other; type++){
    stroke_repeat(cr, &(d->repeat), xstart, xend, lineheight, type, repeat_type_str[type]);
    boxheight += lineheight;
  }

  stroke_rectangle(cr, xstart, xend, y1, boxheight);

  /* show_colors */
  gint x=20, ycen = y1 + lineheight/2;
  cairo_set_line_width(cr, 5.5);
  show_colors(cr, x, &ycen, "ALR/Alpha", CLR_ORANGE);
  show_colors(cr, x, &ycen, "centr-other", CLR_GREEN);
  show_colors(cr, x, &ycen, "telo", CLR_BLUE);
  show_colors(cr, x, &ycen, "acro", CLR_RED);
  show_colors(cr, x, &ycen, "HSATII", CLR_YELLOW);
  show_colors(cr, x, &ycen, "(CATTC)n", CLR_YELLOW2);
  show_colors(cr, x, &ycen, "(GAATG)n", CLR_OLIVE);
  return;
}

static void print_pagetitle(cairo_t *cr, char *chrname, gint page_curr, gint num_page, gint region_no){
  gchar str[1280];
  if(num_page>1) sprintf(str, "%s_%d_%d", chrname, region_no+1, (page_curr+1));
  else sprintf(str, "%s_%d", chrname, region_no+1);

  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  showtext_cr(cr, 50, 30, str, 16);
  return;
}

static void show_colors(cairo_t *cr, gint x, gint *ycen, gchar *label, gdouble r, gdouble g, gdouble b){
  cairo_set_source_rgba(cr, r,g,b, 1);
  rel_xline(cr, x, *ycen, 20);
  cairo_stroke(cr);
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  showtext_cr(cr, x+24, *ycen+4, label, 10);
  *ycen += 15;
  return;
}

static void showtext_cr(cairo_t *cr, gdouble x, gdouble y, gchar *str, gint fontsize){
  cairo_move_to(cr, x, y);
  cairo_select_font_face(cr, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, fontsize);
  cairo_show_text(cr, str);
  cairo_stroke(cr);
  return;
}

static void stroke_rectangle(cairo_t *cr, gint x1, gint x2, gint y1, gint height){
  cairo_set_line_width(cr, 0.4);
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  cairo_rectangle(cr, BP2XAXIS(0), y1, (x2 - x1) * dot_per_bp, height);
  cairo_stroke(cr);
  return;
}

static void fill_rectangle(cairo_t *cr, gint x1, gint x2, gint y1, gint height, gdouble r, gdouble g, gdouble b, double alpha){
  cairo_set_source_rgba(cr, r,g,b, alpha);
  cairo_rectangle(cr, BP2XAXIS(0), y1, (x2 - x1) * dot_per_bp, height);
  cairo_fill(cr);
  return;
}


static void go_show_colors(cairo_t *cr, gint x, gint *y, gchar *label){
  rel_xline(cr, x, *y, 20);
  cairo_stroke(cr);
  cairo_set_source_rgba(cr, CLR_BLACK, 1);
  showtext_cr(cr, x+30, *y+4, label, 14);
  *y += 15;
}

void genome_overlook(DrParam *p, DDParam *d, RefGenome *g){
  cairo_surface_t *surface;
  cairo_t *cr;
  gint i,k, chr;
  gdouble x, len;
  gint lat = 15;
  gchar *filename = alloc_str_new(p->headname, 20);
  sprintf(filename, "%s.pdf", p->headname);
  remove_file(filename);

  dot_per_bp = WIDTH_DRAW/(gdouble)(g->chr[g->chrmax].len +1);
  gint height = OFFSET_Y*2 + 15*d->bednum + (lat + MERGIN_BETWEEN_GOVERLOOK) * g->chrnum;
  surface = cairo_pdf_surface_create(filename, PAGEWIDTH, height);
  cr = cairo_create(surface);

  gint yaxis = OFFSET_Y;
  cairo_set_line_width(cr, 3.5);
  for(k=0; k<d->bednum; k++){
    switch(k){
    case 0: cairo_set_source_rgba(cr, CLR_RED,   1); break;
    case 1: cairo_set_source_rgba(cr, CLR_GREEN, 1); break;
    case 2: cairo_set_source_rgba(cr, CLR_BLUE,  2); break;
    }
    go_show_colors(cr, OFFSET_X, &yaxis, d->bed[k]->name);
  }
  yaxis += lat;

  for(chr=1; chr<g->chrnum; chr++){
    stroke_rectangle(cr, 0, (int)g->chr[chr].len, yaxis-8, 16);
    showtext_cr(cr, OFFSET_X-50, yaxis+5, g->chr[chr].name, 12);

    cairo_set_line_width(cr, lat);
    for(k=0; k<d->bednum; k++){
      switch(k){
      case 0: cairo_set_source_rgba(cr, CLR_RED,   1); break;
      case 1: cairo_set_source_rgba(cr, CLR_GREEN, 1); break;
      case 2: cairo_set_source_rgba(cr, CLR_BLUE,  1); break;
      }
      for(i=0; i<d->bed[k]->chr[chr].num; i++){
	x = BP2XAXIS(d->bed[k]->chr[chr].bed[i].s);
	len = (d->bed[k]->chr[chr].bed[i].e - d->bed[k]->chr[chr].bed[i].s) * dot_per_bp;
	rel_xline(cr, x, yaxis-lat + lat*k, len);
	cairo_stroke(cr);
      }
    }
    yaxis += MERGIN_BETWEEN_GOVERLOOK;
  }

  printf("%s is created.\n", filename);
  cairo_destroy(cr);
  cairo_surface_destroy(surface);
  exit(0);
}
