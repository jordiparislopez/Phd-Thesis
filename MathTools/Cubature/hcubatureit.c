typedef int (*integrandit) (unsigned ndim, const double *x, void *,
	unsigned fdim, double *fval, double *values);

	/* a vector integrand of a vector of npt points: x[i*ndim + j] is the
	j-th coordinate of the i-th point, and the k-th function evaluation
	for the i-th point is returned in fval[i*fdim + k].  Return 0 on success
	or nonzero to terminate the integration. */
	typedef int (*integrandit_v) (unsigned ndim, size_t npt,
		const double *x, void *,
		unsigned fdim, double *fval, double* values);

		/* Different ways of measuring the absolute and relative error when
		we have multiple integrands, given a vector e of error estimates
		in the individual components of a vector v of integrands.  These
		are all equivalent when there is only a single integrand. */



/* vectorized wrapper around non-vectorized integrands */
typedef struct fv_datait_s { integrandit f; void *fdata; } fv_datait;
static int fvit(unsigned ndim, size_t npt,
	      const double *x, void *d_,
	      unsigned fdim, double *fval, double* values)
{
     fv_datait *d = (fv_datait *) d_;
     integrandit f = d->f;
     void *fdata = d->fdata;
     unsigned i;
     /* printf("npt = %u\n", npt); */
     for (i = 0; i < npt; ++i)
	  if (f(ndim, x + i*ndim, fdata, fdim, fval + i*fdim, values))
	       return FAILURE;
     return SUCCESS;
}

struct ruleit_s; /* forward declaration */

typedef int (*evalErrorit_func)(struct ruleit_s *r,
			      unsigned fdim, integrandit_v f, void *fdata,
			      unsigned nR, region *R);
typedef void (*destroyit_func)(struct ruleit_s *r);




typedef struct ruleit_s {
     unsigned dim, fdim;         /* the dimensionality & number of functions */
     unsigned num_points;       /* number of evaluation points */
     unsigned num_regions; /* max number of regions evaluated at once */
     double *pts; /* points to eval: num_regions * num_points * dim */
     double *vals; /* num_regions * num_points * fdim */
     evalErrorit_func evalErrorit;
     destroyit_func destroyit;
} ruleit;

static void destroyit_region(region *R)
{
     destroy_hypercube(&R->h);
     free(R->ee);
     R->ee = 0;
}

static void destroyit_ruleit(ruleit *r)
{
     if (r) {
		if (r->destroyit) r->destroyit(r);
	  free(r->pts);
	  free(r);
     }
}


static int allocit_rule_pts(ruleit *r, unsigned num_regions)
{
     if (num_regions > r->num_regions) {
	  free(r->pts);
	  r->pts = r->vals = NULL;
	  r->num_regions = 0;
	  num_regions *= 2; /* allocate extra so that
			       repeatedly calling alloc_rule_pts with
			       growing num_regions only needs
			       a logarithmic number of allocations */
	  r->pts = (double *) malloc(sizeof(double) *
				     (num_regions
				      * r->num_points * (r->dim + r->fdim)));
	  if (r->fdim + r->dim > 0 && !r->pts) return FAILURE;
	  r->vals = r->pts + num_regions * r->num_points * r->dim;
	  r->num_regions = num_regions;
     }
     return SUCCESS;
}

static ruleit *make_ruleit(size_t sz, /* >= sizeof(rule) */
		       unsigned dim, unsigned fdim, unsigned num_points,
		       evalErrorit_func evalErrorit, destroyit_func destroyit)
{
     ruleit *r;

     if (sz < sizeof(ruleit)) return NULL;
     r = (ruleit *) malloc(sz);
     if (!r) return NULL;
     r->pts = r->vals = NULL;
     r->num_regions = 0;
     r->dim = dim; r->fdim = fdim; r->num_points = num_points;
     r->evalErrorit = evalErrorit;
     r->destroyit = destroyit;
     return r;
}


static ruleit *make_rule75genzmalikit(unsigned dim, unsigned fdim)
{
     rule75genzmalik *r;

     if (dim < 2) return NULL; /* this rule does not support 1d integrals */

     /* Because of the use of a bit-field in evalR_Rfs, we are limited
	to be < 32 dimensions (or however many bits are in unsigned).
	This is not a practical limitation...long before you reach
	32 dimensions, the Genz-Malik cubature becomes excruciatingly
	slow and is superseded by other methods (e.g. Monte-Carlo). */
     if (dim >= sizeof(unsigned) * 8) return NULL;

     r = (rule75genzmalik *) make_rule(sizeof(rule75genzmalik),
				       dim, fdim,
				       num0_0(dim) + 2 * numR0_0fs(dim)
				       + numRR0_0fs(dim) + numR_Rfs(dim),
				       rule75genzmalik_evalError,
				       destroy_rule75genzmalik);
     if (!r) return NULL;

     r->weight1 = (real(12824 - 9120 * to_int(dim) + 400 * isqr(to_int(dim)))
		   / real(19683));
     r->weight3 = real(1820 - 400 * to_int(dim)) / real(19683);
     r->weight5 = real(6859) / real(19683) / real(1U << dim);
     r->weightE1 = (real(729 - 950 * to_int(dim) + 50 * isqr(to_int(dim)))
		    / real(729));
     r->weightE3 = real(265 - 100 * to_int(dim)) / real(1458);

     r->p = (double *) malloc(sizeof(double) * dim * 3);
     if (!r->p) { destroyit_ruleit((ruleit *) r); return NULL; }
     r->widthLambda = r->p + dim;
     r->widthLambda2 = r->p + 2 * dim;

     return (ruleit *) r;
}

/***************************************************************************/
/* 1d 15-point Gaussian quadrature rule, based on qk15.c and qk.c in
   GNU GSL (which in turn is based on QUADPACK). */
   static int rule15gaussit_evalErrorit(ruleit *r,
	   unsigned fdim, integrandit_v f, void *fdata,
	   unsigned nR, region *R)
	   {
		   /* Gauss quadrature weights and kronrod quadrature abscissae and
		   weights as evaluated with 80 decimal digit arithmetic by
		   L. W. Fullerton, Bell Labs, Nov. 1981. */
		   const unsigned n = 8;
		   const double xgk[8] = {  /* abscissae of the 15-point kronrod rule */
			   0.991455371120812639206854697526329,
			   0.949107912342758524526189684047851,
			   0.864864423359769072789712788640926,
			   0.741531185599394439863864773280788,
			   0.586087235467691130294144838258730,
			   0.405845151377397166906606412076961,
			   0.207784955007898467600689403773245,
			   0.000000000000000000000000000000000
			   /* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule.
			   xgk[0], xgk[2], ... to optimally extend the 7-point gauss rule */
		   };
		   static const double wg[4] = {  /* weights of the 7-point gauss rule */
			   0.129484966168869693270611432679082,
			   0.279705391489276667901467771423780,
			   0.381830050505118944950369775488975,
			   0.417959183673469387755102040816327
		   };
		   static const double wgk[8] = { /* weights of the 15-point kronrod rule */
			   0.022935322010529224963732008058970,
			   0.063092092629978553290700663189204,
			   0.104790010322250183839876322541518,
			   0.140653259715525918745189590510238,
			   0.169004726639267902826583426598550,
			   0.190350578064785409913256402421014,
			   0.204432940075298892414161999234649,
			   0.209482141084727828012999174891714
		   };
		   unsigned j, k, iR;
		   size_t npts = 0;
		   double *pts, *vals;

		   if (allocit_rule_pts(r, nR)) return FAILURE;
		   pts = r->pts; vals = r->vals;

		   for (iR = 0; iR < nR; ++iR) {
			   const double center = R[iR].h.data[0];
			   const double halfwidth = R[iR].h.data[1];

			   pts[npts++] = center;

			   for (j = 0; j < (n - 1) / 2; ++j) {
				   int j2 = 2*j + 1;
				   double w = halfwidth * xgk[j2];
				   pts[npts++] = center - w;
				   pts[npts++] = center + w;
			   }
			   for (j = 0; j < n/2; ++j) {
				   int j2 = 2*j;
				   double w = halfwidth * xgk[j2];
				   pts[npts++] = center - w;
				   pts[npts++] = center + w;
			   }

			   R[iR].splitDim = 0; /* no choice but to divide 0th dimension */
		   }

		   if (f(1, npts, pts, fdata, fdim, vals,vals))
		   return FAILURE;

		   for (k = 0; k < fdim; ++k) {
			   const double *vk = vals + k;
			   for (iR = 0; iR < nR; ++iR) {
				   const double halfwidth = R[iR].h.data[1];
				   double result_gauss = vk[0] * wg[n/2 - 1];
				   double result_kronrod = vk[0] * wgk[n - 1];
				   double result_abs = fabs(result_kronrod);
				   double result_asc, mean, err;

				   /* accumulate integrals */
				   npts = 1;
				   for (j = 0; j < (n - 1) / 2; ++j) {
					   int j2 = 2*j + 1;
					   double v = vk[fdim*npts] + vk[fdim*npts+fdim];
					   result_gauss += wg[j] * v;
					   result_kronrod += wgk[j2] * v;
					   result_abs += wgk[j2] * (fabs(vk[fdim*npts])
					   + fabs(vk[fdim*npts+fdim]));
					   npts += 2;
				   }
				   for (j = 0; j < n/2; ++j) {
					   int j2 = 2*j;
					   result_kronrod += wgk[j2] * (vk[fdim*npts]
						   + vk[fdim*npts+fdim]);
						   result_abs += wgk[j2] * (fabs(vk[fdim*npts])
						   + fabs(vk[fdim*npts+fdim]));
						   npts += 2;
					   }

					   /* integration result */
					   R[iR].ee[k].val = result_kronrod * halfwidth;

					   /* error estimate
					   (from GSL, probably dates back to QUADPACK
					   ... not completely clear to me why we don't just use
					   fabs(result_kronrod - result_gauss) * halfwidth */
					   mean = result_kronrod * 0.5;
					   result_asc = wgk[n - 1] * fabs(vk[0] - mean);
					   npts = 1;
					   for (j = 0; j < (n - 1) / 2; ++j) {
						   int j2 = 2*j + 1;
						   result_asc += wgk[j2] * (fabs(vk[fdim*npts]-mean)
						   + fabs(vk[fdim*npts+fdim]-mean));
						   npts += 2;
					   }
					   for (j = 0; j < n/2; ++j) {
						   int j2 = 2*j;
						   result_asc += wgk[j2] * (fabs(vk[fdim*npts]-mean)
						   + fabs(vk[fdim*npts+fdim]-mean));
						   npts += 2;
					   }
					   err = fabs(result_kronrod - result_gauss) * halfwidth;
					   result_abs *= halfwidth;
					   result_asc *= halfwidth;
					   if (result_asc != 0 && err != 0) {
						   double scale = pow((200 * err / result_asc), 1.5);
						   err = (scale < 1) ? result_asc * scale : result_asc;
					   }
					   if (result_abs > DBL_MIN / (50 * DBL_EPSILON)) {
						   double min_err = 50 * DBL_EPSILON * result_abs;
						   if (min_err > err) err = min_err;
					   }
					   R[iR].ee[k].err = err;

					   /* increment vk to point to next batch of results */
					   vk += 15*fdim;
				   }
			   }
			   return SUCCESS;
		   }

static ruleit *make_rule15gaussit(unsigned dim, unsigned fdim)
{
    if (dim != 1) return NULL; /* this rule is only for 1d integrals */

    return make_ruleit(sizeof(ruleit), dim, fdim, 15,
		      rule15gaussit_evalErrorit, 0);
}






static int eval_regionsit(unsigned nR, region *R,
			integrandit_v f, void *fdata, ruleit *r)
{
     unsigned iR;
     if (nR == 0) return SUCCESS; /* nothing to evaluate */
     if (r->evalErrorit(r, R->fdim, f, fdata, nR, R)) return FAILURE;
     for (iR = 0; iR < nR; ++iR)
	  R[iR].errmax = errMax(R->fdim, R[iR].ee);
     return SUCCESS;
}






static int rulecubatureiterative(ruleit *r, unsigned fdim,
			integrandit_v f, void *fdata,
			const hypercube *h,
			size_t maxEval,
			double reqAbsError, double reqRelError,
			error_norm norm,
			double *val, double *err, int parallel)
			{
     size_t numEval = 0;
     heap regions;
     unsigned i, j;
     region *R = NULL; /* array of regions to evaluate */
     size_t nR_alloc = 0;
     esterr *ee = NULL;

     if (fdim <= 1) norm = ERROR_INDIVIDUAL; /* norm is irrelevant */
     if (norm < 0 || norm > ERROR_LINF) return FAILURE; /* invalid norm */

     regions = heap_alloc(1, fdim);
     if (!regions.ee || !regions.items) goto bad;

     ee = (esterr *) malloc(sizeof(esterr) * fdim);
     if (!ee) goto bad;

     nR_alloc = 2;
     R = (region *) malloc(sizeof(region) * nR_alloc);
     if (!R) goto bad;
     R[0] = make_region(h, fdim);
     if (!R[0].ee
	 || eval_regionsit(1, R, f, fdata, r)
	 || heap_push(&regions, R[0]))
	       goto bad;
     numEval += r->num_points;

     while (numEval < maxEval || !maxEval) {
	  if (converged(fdim, regions.ee, reqAbsError, reqRelError, norm))
	       break;

	  if (parallel) {

	       size_t nR = 0;
	       for (j = 0; j < fdim; ++j) ee[j] = regions.ee[j];
	       do {
		    if (nR + 2 > nR_alloc) {
			 nR_alloc = (nR + 2) * 2;
			 R = (region *) realloc(R, nR_alloc * sizeof(region));
			 if (!R) goto bad;
		    }
		    R[nR] = heap_pop(&regions);
		    for (j = 0; j < fdim; ++j) ee[j].err -= R[nR].ee[j].err;
		    if (cut_region(R+nR, R+nR+1)) goto bad;
		    numEval += r->num_points * 2;
		    nR += 2;
		    if (converged(fdim, ee, reqAbsError, reqRelError, norm))
			 break; /* other regions have small errs */
	       } while (regions.n > 0 && (numEval < maxEval || !maxEval));
	       if (eval_regionsit(nR, R, f, fdata, r)
		   || heap_push_many(&regions, nR, R))
		    goto bad;
	  }
	  else { /* minimize number of function evaluations */
	       R[0] = heap_pop(&regions); /* get worst region */
	       if (cut_region(R, R+1)
		   || eval_regionsit(2, R, f, fdata, r)
		   || heap_push_many(&regions, 2, R))
		    goto bad;
	       numEval += r->num_points * 2;
	  }
     }

     /* re-sum integral and errors */

	 for (j = 0; j < fdim; ++j) val[j] = err[j] = 0;
     for (i = 0; i < regions.n; ++i) {
	  for (j = 0; j < fdim; ++j) {
	       val[j] += regions.items[i].ee[j].val;
	       err[j] += regions.items[i].ee[j].err;
	  }
	  destroy_region(&regions.items[i]);
     }
     /* printf("regions.nalloc = %d\n", regions.nalloc); */
     free(ee);
     heap_free(&regions);
     free(R);

     return SUCCESS;

	 bad:
     free(ee);
     heap_free(&regions);
     free(R);
     return FAILURE;
}

static int cubatureiterative(unsigned fdim, integrandit_v f, void *fdata,
		    unsigned dim, const double *xmin, const double *xmax,
		    size_t maxEval, double reqAbsError, double reqRelError,
		    error_norm norm,
		    double *val, double *err, int parallel)
{
     ruleit *r;
     hypercube h;
     int status;
     unsigned i;

     if (fdim == 0) /* nothing to do */ return SUCCESS;
     if (dim == 0) { /* trivial integration */
	  if (f(0, 1, xmin, fdata, fdim, val,val)) return FAILURE;
	  for (i = 0; i < fdim; ++i) err[i] = 0;
	  return SUCCESS;
     }
     r = dim == 1 ? make_rule15gaussit(dim, fdim)
 	          : make_rule75genzmalikit(dim, fdim);
     if (!r) {
	  for (i = 0; i < fdim; ++i) {
	       val[i] = 0;
	       err[i] = HUGE_VAL;
	  }
	  return FAILURE;
     }

	 h = make_hypercube_range(dim, xmin, xmax);
     status = !h.data ? FAILURE
	  : rulecubatureiterative(r, fdim, f, fdata, &h,
				maxEval, reqAbsError, reqRelError, norm,
				val, err, parallel);
				     destroy_hypercube(&h);

     destroyit_ruleit(r);

     return status;
}



int hcubatureiterative(unsigned fdim, integrandit f, void *fdata,
	      unsigned dim, const double *xmin, const double *xmax,
	      size_t maxEval, double reqAbsError, double reqRelError,
	      error_norm norm,
	      double *val, double *err)
{
     int i,ret;
     fv_datait d;

     if (fdim == 0) return SUCCESS; /* nothing to do */

     d.f = f; d.fdata = fdata;

	 for(i=0; i<fdim; i++){val[i]=43;}
     ret = cubatureiterative(fdim, fvit, &d, dim, xmin, xmax,
		    maxEval, reqAbsError, reqRelError, norm, val, err, 0);

			cout << val[0] << endl;


		//  ret = cubatureiterative(fdim, fvit, &d, dim, xmin, xmax,
		// 	   maxEval, reqAbsError, reqRelError, norm, val, err, 0);

     return ret;
}
