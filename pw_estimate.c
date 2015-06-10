#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>
#include "pw_estimate.h"
#include "macro.h"

typedef struct{
  TYPE_WIGARRAY *wig;
  TYPE_WIGARRAY *mp;
  int binnum;
} ParStr;

int thre;

static void func_iteration(gsl_multimin_fminimizer *s);
//static double f_nb(const gsl_vector *v, void *params);
static double f_zinb_reg(const gsl_vector *v, void *params);
static double f_zinb_const(const gsl_vector *v, void *params);

static gsl_multimin_fminimizer *gsl_multimin_fminimizer_new(size_t ndim){
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;  /* ネルダーとミードのシンプレックス法 */
  gsl_multimin_fminimizer            *s = gsl_multimin_fminimizer_alloc(T, ndim); 
  return s;
}

static gsl_vector *gsl_vector_new(int ndim, double init){
  gsl_vector *p = gsl_vector_alloc(ndim);
  gsl_vector_set_all(p, init);
  return p;
}

double pw_get_negative_binomial(int k, double p, double n){
  return gsl_ran_negative_binomial_pdf(k, p, n);
}

double pw_get_zeroinflated_negative_binomial(int k, double p, double n, double p0){
  double r;
  if(!k){
    r = p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n);
  }else{
    r = (1 - p0) * gsl_ran_negative_binomial_pdf(k, p, n);
  }
  return r;
}

double pw_get_poisson(int k, int lambda){
  return gsl_ran_poisson_pdf(k, lambda);
}

void pw_estimate_nb_param(Mapfile *mapfile, RefGenome *g){
  int i;
  size_t ndim=3;

  /* initialization */
  thre = mapfile->wstats.thre;
  double par[thre];
  for(i=0; i<thre; i++){
    if(!mapfile->wstats.genome->darray_bg[i]) par[i] = 0;
    //    else par[i] = -log(mapfile->wstats.genome->darray_bg[i] /(double)mapfile->wstats.genome->num);
    else par[i] = mapfile->wstats.genome->darray_bg[i] /(double)mapfile->wstats.genome->num;
    LOG("par[%d]=%f, darray_bg=%d, num=%d\n", i, par[i], mapfile->wstats.genome->darray_bg[i], mapfile->wstats.genome->num);
  }
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_new(ndim);
  gsl_vector *x = gsl_vector_alloc(ndim);
  gsl_vector_set(x, 0, mapfile->wstats.chr[g->chrmax].nb_p);
  gsl_vector_set(x, 1, mapfile->wstats.chr[g->chrmax].nb_n);
  gsl_vector_set(x, 2, 0.02); /* p0 */
  gsl_vector *ss = gsl_vector_new(ndim, 0.1); /* step size */

  gsl_multimin_function minex_func;
  minex_func.n = ndim;
  minex_func.f = &f_zinb_const;
  minex_func.params = (void *)&par;
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  func_iteration(s);

  mapfile->wstats.genome->nb_p  = gsl_vector_get(s->x, 0);
  if(mapfile->wstats.genome->nb_p<=0) mapfile->wstats.genome->nb_p=0.01;
  if(mapfile->wstats.genome->nb_p>=1) mapfile->wstats.genome->nb_p=0.99;
  mapfile->wstats.genome->nb_n  = gsl_vector_get(s->x, 1);
  mapfile->wstats.genome->nb_p0 = gsl_vector_get(s->x, 2);
  if(mapfile->wstats.genome->nb_p0<0) mapfile->wstats.genome->nb_p0=0.0;
  if(mapfile->wstats.genome->nb_p0>1) mapfile->wstats.genome->nb_p0=1.0;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  return;
}

void pw_estimate_zinb_opt(Mapfile *mapfile, TYPE_WIGARRAY *wigarray, TYPE_WIGARRAY *mparray, int chr, int binnum){
  size_t ndim=3;

  /* initialization */
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_new(ndim);

  gsl_vector *x = gsl_vector_alloc(ndim);
  gsl_vector_set(x, 0, mapfile->wstats.chr[chr].nb_p);
  gsl_vector_set(x, 1, mapfile->wstats.chr[chr].nb_n);
  gsl_vector_set(x, 2, 0.2); /* for Beta function */
  gsl_vector *ss = gsl_vector_new(ndim, 0.1); /* step size */

  ParStr par;
  par.wig = wigarray;
  par.mp  = mparray;
  par.binnum = binnum;
  gsl_multimin_function minex_func;
  minex_func.n = ndim;
  minex_func.f = &f_zinb_reg;
  minex_func.params = (void *)&par;

  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  func_iteration(s);

  mapfile->wstats.genome->nb_p = gsl_vector_get(s->x, 0);
  if(mapfile->wstats.genome->nb_p<0) mapfile->wstats.genome->nb_p=0;
  if(mapfile->wstats.genome->nb_p>1) mapfile->wstats.genome->nb_p=1.0;
  mapfile->wstats.genome->nb_n = gsl_vector_get(s->x, 1);
  mapfile->wstats.genome->nb_p0 = gsl_vector_get(s->x, 2);
  if(mapfile->wstats.genome->nb_p0<0) mapfile->wstats.genome->nb_p0=0;
  if(mapfile->wstats.genome->nb_p0>1) mapfile->wstats.genome->nb_p0=1.0;

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  
  exit(0);

  return;
}

static void func_iteration(gsl_multimin_fminimizer *s){
  size_t iter = 0;
  int status;
  double size;
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);  /* 繰り返し計算を1回行う。予期しない問題が発生した場合はエラーコードを返す。 */
    if(status) break;
    size = gsl_multimin_fminimizer_size(s);       /* sのその時点での最小点の最良推定値を返す */
    status = gsl_multimin_test_size(size, 1e-3);  /* sizeが閾値(1e-3)より小さければGSL_SUCCESS を、そうでなければGSL_CONTINUEを返す。 */
    if(status == GSL_SUCCESS){ LOG("converged to minimum at\n");}
    LOG("%5zd p=%.5f n=%.5f p0=%.5f f() = %10.5f size = %f\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), gsl_vector_get(s->x, 2), gsl_multimin_fminimizer_minimum(s), size);
  }while(status == GSL_CONTINUE && iter < 1000);

  return;
}

static double f_zinb_reg(const gsl_vector *v, void *params){
  int i, binnum;
  double p, n, m, p0, r, fxy=0;
  ParStr *par = (ParStr *)params;
  p = gsl_vector_get(v, 0);
  n = gsl_vector_get(v, 1);
  m = gsl_vector_get(v, 2);
  printf("zinb_reg p=%f, n=%f, m=%f\n",p, n, m);

  binnum = par->binnum;
  TYPE_WIGARRAY *wig = par->wig;
  TYPE_WIGARRAY *mp = par->mp;
  for(i=0; i<binnum; i++){  /* wig[i]が得られる確率を最大化 */
    if(!mp[i]) p0=1; else p0 = gsl_sf_beta(WIGARRAY2VALUE(mp[i]), m);
    printf("%d p0=%f\n", i, p0);
    if(!wig[i]){
      r = p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n);
    }else{
      r = (1 - p0) * gsl_ran_negative_binomial_pdf(WIGARRAY2VALUE(wig[i]), p, n);
    }
    fxy += log(r);
  }
  printf("fxy=%f\n",fxy);
  return fxy;
}

static double f_zinb_const(const gsl_vector *v, void *params){
  int i;
  double p, n, r, fxy=0, p0;
  double *par = (double *)params;
  p = gsl_vector_get(v, 0);
  if(p<=0) p=0.01;
  if(p>=1) p=0.99;
  n = gsl_vector_get(v, 1);
  p0 = gsl_vector_get(v, 2);
  if(p0<0) p0=0;
  if(p0>1) p0=1.0;
  LOG("zinb_const p=%f, n=%f, p0=%f\n", p, n, p0);
  
  for(i=0; i<thre; i++){
    /*    if(!i) r = -log(p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n));
	  else   r = -log(     (1 - p0) * gsl_ran_negative_binomial_pdf(i, p, n));*/
    if(!i) r = p0 + (1 - p0) * gsl_ran_negative_binomial_pdf(0, p, n);
    else   r =      (1 - p0) * gsl_ran_negative_binomial_pdf(i, p, n);
    //    printf("%d: %f - %f\n", i, par[i] , r);
    fxy += (par[i] - r)*(par[i] - r);
  }
  return fxy;
}

/*static double f_nb(const gsl_vector *v, void *params){
  int i;
  double p, n, r, fxy=0;
  double *par = (double *)params;
  p = gsl_vector_get(v, 0);
  n = gsl_vector_get(v, 1);
  for(i=0; i<N_ARRAY; i++){
    r = -log(gsl_ran_negative_binomial_pdf(i, p, n));
    fxy += (par[i] - r)*(par[i] - r);
  }
  return fxy;
  }*/
