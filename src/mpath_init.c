#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
 *   The following name(s) appear with different usages
 *     e.g., with different numbers of arguments:
 *
 *         penGLM
 *
 *           This needs to be resolved in the tables and any declarations.
 *           */

/* FIXME: 
 *    Check these declarations against the C/Fortran source code.
 *    */

/* .Fortran calls */
extern void F77_NAME(checkconvergence)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(compute_h)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(copymatrix)(void *, void *, void *, void *);
extern void F77_NAME(deveval)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(find_activeset)(void *, void *, void *, void *, void *);
extern void F77_NAME(gfunc)(void *, void *, void *, void *, void *);
extern void F77_NAME(glmlink)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(glmreg_fit_fortran)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(glmregnb_fortran)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gradient)(void *, void *, void *, void *, void *);
extern void F77_NAME(linkinv)(void *, void *, void *, void *);
extern void F77_NAME(lmax_zipath)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void*);
extern void F77_NAME(loglikfor)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(loss)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nclreg_ad)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nclreg_fortran)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nonconvexloss)(void *, void *, void *, void *);
extern void F77_NAME(outloop)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(outprod)(void *, void *, void *, void *, void *);
extern void F77_NAME(penglm)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pnorm_fortran)(void *);
extern void F77_NAME(pred)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(theta_ml)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zeval)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ziloss)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zi_onelambda)(void *, void *, void *, void *, void *, void*, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zipath_active)(void *, void *, void*, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zipath_nonactive)(void *, void *, void*, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
	    {"checkconvergence",           (DL_FUNC) &F77_NAME(checkconvergence),    8},
	    {"compute_h",                  (DL_FUNC) &F77_NAME(compute_h),           7},
	    {"copymatrix",                 (DL_FUNC) &F77_NAME(copymatrix),          4},
	    {"deveval",                    (DL_FUNC) &F77_NAME(deveval),             7},
	    {"find_activeset",             (DL_FUNC) &F77_NAME(find_activeset),      5},
	    {"gfunc",                      (DL_FUNC) &F77_NAME(gfunc),               5},
	    {"glmlink",                    (DL_FUNC) &F77_NAME(glmlink),             6},
	    {"glmreg_fit_fortran",         (DL_FUNC) &F77_NAME(glmreg_fit_fortran), 29},
	    {"glmregnb_fortran",           (DL_FUNC) &F77_NAME(glmregnb_fortran),   31},
	    {"gradient",                   (DL_FUNC) &F77_NAME(gradient),            5},
	    {"linkinv",                    (DL_FUNC) &F77_NAME(linkinv),             4},
	    {"lmax_zipath",                (DL_FUNC) &F77_NAME(lmax_zipath),        18},
	    {"loglikfor",                  (DL_FUNC) &F77_NAME(loglikfor),           7},
	    {"loss",                       (DL_FUNC) &F77_NAME(loss),                7},
	    {"nclreg_ad",                  (DL_FUNC) &F77_NAME(nclreg_ad),          35},
	    {"nclreg_fortran",             (DL_FUNC) &F77_NAME(nclreg_fortran),     36},
	    {"nonconvesloss",              (DL_FUNC) &F77_NAME(nonconvexloss),       4},
	    {"outloop",                    (DL_FUNC) &F77_NAME(outloop),            36},
	    {"outprod",                    (DL_FUNC) &F77_NAME(outprod),             5},
	    {"penglm",                     (DL_FUNC) &F77_NAME(penglm),              7},
	    {"pnorm_fortran",              (DL_FUNC) &F77_NAME(pnorm_fortran),       1},
	    {"pred",                       (DL_FUNC) &F77_NAME(pred),                9},
	    {"theta_ml",                   (DL_FUNC) &F77_NAME(theta_ml),            8},
	    {"zeval",                      (DL_FUNC) &F77_NAME(zeval),               7},
	    {"ziloss",                     (DL_FUNC) &F77_NAME(ziloss),             10},
	    {"zi_onelambda",               (DL_FUNC) &F77_NAME(zi_onelambda),       43},
	    {"zipath_active",              (DL_FUNC) &F77_NAME(zipath_active),      41},
	    {"zipath_nonactive",           (DL_FUNC) &F77_NAME(zipath_nonactive),   41},
					        {NULL, NULL, 0}
};

void R_init_mpath(DllInfo *dll)
{
	    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
	        R_useDynamicSymbols(dll, FALSE);
}

