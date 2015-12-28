#include "utils.h"
//#define DEBUG 1

/* Pointers to not-yet-instantiated structures. These will be referred to globally,
 * so I'll make a cache of them here so that any conflicts can be resolved.
 */

pval_finder * pcomp=NULL;
double const * theta_ptr=NULL;
double (*fp)(int, double*) = NULL;
double* theta_old=NULL, * temp_theta=NULL;
double* old_sums=NULL, * old_total=NULL;
double * constraints=NULL;

/* nmmin-compatible functions for optimisation. They need to have the same 
 * face in order to allow the use of function pointers to switch the function
 * in the barrier computation. */

double get_pval(int n, double* par) {
	pcomp->compute_pval(par);
    return -std::log(pcomp->pval);	
}

double get_sumsq(int n, double* par) {
	double stuff=0, diff;
	for (int i=0; i<n; ++i) { 
		diff = theta_ptr[i] - par[i];
		stuff += diff * diff;
	}
	pcomp->summate(par); // Need to do this so barrier function has correct totals to work on.
	return stuff;		
}


/* Computes the function + barrier value, given the old and new theta values. 
 * This barrier function constrains the proportions to obvious values (between
 * 0 and a sum of unity) and the less obvious (constraints on the partial sums
 * to enforce a positive beta).
 */

double barrier_fn(int n, double* par, void* ex) { 
	const double mu=0.0001;
    double bval=0, fnval=(*fp)(n, par);
    int i;
    for (i=0; i<n; ++i) { bval+=std::log(par[i])*theta_old[i] - par[i]; }
    bval += (1 - *old_total)*std::log(1 - pcomp->total) + pcomp->total;

    if (!ISNA(constraints[0])) { // Protects against NA constraints when current p-value is 1.
        for (i=0; i<(pcomp->nexp); ++i) { 
            bval += (constraints[i] - old_sums[i])*std::log(constraints[i] - pcomp->total_per_exp[i]) + pcomp->total_per_exp[i]; 
        }
    }
    if (!R_FINITE(bval)) { return R_PosInf; }
    return fnval - bval*mu; 
}


/* Constrained optimisation, with multiple outer iterations to refine the log
 * barrier.  We use Nelder-Mead as BFGS is a bit more fragile regarding how it
 * deals with NA values.  Even a pure R implementation of BFGS with constrOptim
 * falls apart, so it's nothing to do with the code.
 */

const int outer_it=100, maxit=500; 
const double abstol=R_NegInf, reltol=0.00000001, outer_tol = 0.00001, nm_alpha=1, nm_beta=0.5, nm_gamma=2;
#ifdef DEBUG
int tracer=0;
#else
const int tracer=0;
#endif

double constr_optim(double* current_theta) { 
	int fail=0, fncount=0;
	
	// Setting up inputs. 
	const int& nprops=pcomp->nprops;
	int ne=0, prop;
	for (prop=0; prop<nprops; ++prop) { temp_theta[prop]=theta_old[prop]=current_theta[prop]; }
	double oldtrueval=(*fp)(nprops, theta_old), trueval=oldtrueval;
	for (ne=0; ne<(pcomp->nexp); ++ne) { old_sums[ne]=pcomp->total_per_exp[ne]; }
	*old_total=pcomp->total;
	double oldbarval=barrier_fn(nprops, theta_old, NULL), barval=oldbarval; 

	for (int oit=0; oit<outer_it; ++oit){
		// Nelder-Mead; Need to use temp_theta, as nmmin changes both temp_theta and current_theta (surprise!).
		nmmin(nprops, temp_theta, current_theta, &barval, barrier_fn, &fail, abstol, reltol, 
				NULL, nm_alpha, nm_beta, nm_gamma, tracer, &fncount, maxit);
		for (prop=0; prop<nprops; ++prop) { temp_theta[prop]=theta_old[prop]=current_theta[prop]; }
		trueval=(*fp)(nprops, current_theta);
		for (ne=0; ne<(pcomp->nexp); ++ne) { old_sums[ne]=pcomp->total_per_exp[ne]; }
		*old_total=pcomp->total;

		if (! (R_FINITE(barval) && R_FINITE(oldbarval)) ) {
			throw std::runtime_error("Nelder-Mead optimisation returned a NA value from non-NA start"); 
		} else if (std::abs(barval-oldbarval) < (std::abs(barval)+0.001)*outer_tol) {
			break;					
		}
		if (trueval > oldtrueval) { break; }
		oldbarval = barval;
		oldtrueval = trueval;
	}
	return trueval;
}

/**************
 * Main loop. *
 **************/

extern "C" {

SEXP R_load_stuff(SEXP alpha, SEXP mprop, SEXP theta, SEXP combos, SEXP order) try {

	// Checking the p-values, and the ordering thereof.
	if (!isNumeric(alpha)) { throw std::runtime_error("alpha must be a double-precision scalar"); }
	const int ngenes=LENGTH(alpha);
	const double* aptr=REAL(alpha);
    if (!isInteger(order)) { throw std::runtime_error("ordering must be an integer vector"); }
	if (LENGTH(order)!=ngenes) { throw std::runtime_error("length of ordering vector should be equal to the number of tests"); }
	const int* optr=INTEGER(order);

	// Checking the M values.
	if (!isNumeric(mprop)) { throw std::runtime_error("m-proportions must be a double-precision vector"); }
	const int nexp=LENGTH(mprop)/ngenes;
	double const** mptrs=(double const**)R_alloc(nexp, sizeof(double*));
	mptrs[0]=REAL(mprop);
	for (int i=1; i<nexp; ++i) { mptrs[i]=mptrs[i-1]+ngenes; }

	// Checking our starting proportions.
	if (!isNumeric(theta)) { throw std::runtime_error("theta proportions must be a double-precision vector"); }
	const int nprops=LENGTH(theta);
	theta_ptr=REAL(theta);

	// Checking the combinations.
	if (!isInteger(combos)) { throw std::runtime_error("combination matrix must be a double-precision matrix"); }
	if (LENGTH(combos)!=nprops*nexp) { throw std::runtime_error("combination matrix has invalid dimensions"); }
	pval_finder px(nexp, nprops, INTEGER(combos));
	pcomp=&px;
	
	/********************************************************
 	 * checking whether our starting proportions are okay. 
 	 ********************************************************/
	
	bool initok=true;
	double* ref_sum_per_exp=(double*)R_alloc(nexp, sizeof(double));
	{ 
		double total=0;
		for (int i=0; i<nprops; ++i) {
			total += theta_ptr[i];
			if (theta_ptr[i] <= 0 || theta_ptr[i]>=1) { 
				initok=false;
				break;
			}
		}
		if (total >= 1) { initok=false; }
		px.summate(theta_ptr);
		for (int j=0; j<nexp; ++j) { ref_sum_per_exp[j]=px.total_per_exp[j]; }
	}
	double* good_enough=(double*)R_alloc(nprops, sizeof(double));
	for (int prop=0; prop<nprops; ++prop) { good_enough[prop]=theta_ptr[prop]; }
	
	// Initializing other odds and ends.
	theta_old=(double*)R_alloc(nprops, sizeof(double));
	temp_theta=(double*)R_alloc(nprops, sizeof(double));
	old_sums=(double*)R_alloc(nexp+1, sizeof(double));
	old_total=old_sums+nexp;
	constraints=(double*)R_alloc(nexp, sizeof(double));

	SEXP output=PROTECT(allocVector(REALSXP, ngenes));
try {	
	double* opptr=REAL(output);
	double* current_theta=(double*)R_alloc(nprops, sizeof(double));

	/********************************************************
 	 * Searching through all genes.
 	 ********************************************************/

#ifdef DEBUG
	int tester=62; // Debugging for this gene.
	--tester;
#endif

	for (int odex=0; odex<ngenes; ++odex) {
		const int& gene=optr[odex];
		const double& cur_alpha=aptr[gene];
		pcomp->set_alpha(cur_alpha);
	
#ifdef DEBUG
		if (gene==tester) { 
			Rprintf("Gene: %i (%i) %.10f\n", gene, odex, cur_alpha);
		}
#endif		

		/* Checking if the theta estimates are already within the feasible
 		 * region. We skip the checks for 0 and 1'ness, because they're
 		 * obvious. We also compute the constraints with some protection 
 		 * against undefined values when cur_alpha is 1. Also, mptrs must 
 		 * be positive which avoids impossible constraints.
 		 */
		bool notokay=false;
		for (int exp=0; exp<nexp; ++exp) { 
			pcomp->mvals[exp]=mptrs[exp][gene];
			if (1 - cur_alpha > 0.00000001) {
				constraints[exp]=mptrs[exp][gene]/(1-cur_alpha); 
				if (!notokay && initok && ref_sum_per_exp[exp] >= constraints[exp]) { notokay=true; }
			} else {  constraints[exp]=R_NaReal; }
		}

		// Dummy variables, so barrier_fn doesn't get uninitialized values for R_FINITE check.
		(*old_total)=0.5;
		for (int exp=0; exp<nexp; ++exp) { old_sums[exp]=0.5; }

		/* Finding the closest point in the feasible region to the estimated theta.
		 * We pick a point in the feasible region and we minimize the sum of squared
		 * differences to the estimated theta. This should work as the space is
		 * convex i.e. there should only be one  maximum to a point outside the space.
		 * We try to hot-start it from the previous gene, if possible. The ordering
		 * is such that the constraints should expand, so hot-starting should cut time.
		 */
		if (notokay || !initok) {
			fp=get_sumsq;
			if (good_enough[0] > 0 && R_FINITE(barrier_fn(nprops, good_enough, NULL))) { 
				for (int prop=0; prop < nprops; ++prop) { current_theta[prop]=good_enough[prop]; }
			} else {
				double minval=1;
				for (int exp=0; exp<nexp; ++exp) { 
					if (constraints[exp] < minval) { minval=constraints[exp]; }
				}
				minval /= nprops + 1;
				for (int prop=0; prop<nprops; ++prop) { current_theta[prop]=minval; }
			}
#ifdef DEBUG			
			if (tester==gene) { 
				Rprintf("Old proportions are:");
				for (int prop=0; prop<nprops; ++prop) { Rprintf("\t%.10f", current_theta[prop]); }
				Rprintf("\n");
				tracer=1;
			}
#endif
			constr_optim(current_theta);			
			for (int prop=0; prop<nprops; ++prop) { good_enough[prop]=current_theta[prop]; }
		} else {
			for (int prop=0; prop<nprops; ++prop) { current_theta[prop]=theta_ptr[prop]; }
			good_enough[0]=-1;
		}

#ifdef DEBUG
		if (tester==gene) { 
			Rprintf("New proportions are:");
			for (int prop=0; prop<nprops; ++prop) { Rprintf("\t%.10f", current_theta[prop]); }
			Rprintf("\n");
			tracer=0;
		}
#endif

		// Using the discovered point for constrained optimization to max out the p-value.
	  	fp=get_pval;
	 	opptr[gene]=std::exp(-constr_optim(current_theta));

#ifdef DEBUG
		if (tester==gene) {
			Rprintf("P-value is %.10f\n", opptr[gene]);
			Rprintf("Final proportions are:");
			for (int prop=0; prop<nprops; ++prop) { Rprintf("\t%.10f", current_theta[prop]); }
			Rprintf("\n");
			tracer=0;
		}
#endif		
	}
} catch (std::exception& e){
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

}

/**********************
 * Old BFGS stuff. ****
 **********************/

//void deriv_sumsq(int n, double* par, double* gr) {
//	for (int i=0; i<n; ++i) { gr[i]=2*(par[i]-theta_ptr[i]); }	
//	return;
//}

//double *beta_deriv=NULL, *pval_coeff=NULL;
//void deriv_pval(int n, double* par, double* gr) {
//	int i, j, k;
//	double temp;
//
//	// Computing the 1-beta values for each experiment, and their derivatives.
//	for (i=0; i<(*nexperiments); ++i) {
//		temp=0;
//		for (j=0; j < n; ++j) {
//			if (combo_ptrs[j][i]) { temp+=par[j]; }
//		}
//		beta_ptr[i]=(1-mval_ptr[i]-*curalpha*temp)/(1-temp);
//		beta_deriv[i]=beta_ptr[i]/(1-temp); 
//	}
//
//	// Running through and computing the p-value coefficients.
//	double pval=0, total=0;
//	for (j=0; j<n; ++j) {
//		temp=par[j];
//		for (i=0; i<(*nexperiments); ++i) { 
//			if (combo_ptrs[j][i]) { temp*=(*curalpha); }
//			else { temp*=beta_ptr[i]; }
//		}
//		pval += temp;
//		pval_coeff[j] = temp;
//		total += par[j];
//	}
//
//	// Running through and using everything for differentiation.
//	for (j=0; j<n; ++j) { gr[j] = pval_coeff[j]/par[j]; }
//	for (i=0; i<(*nexperiments); ++i) {
//		for (j=0; j < n; ++j) {
//			if (combo_ptrs[j][i]) { 
//				for (k=0; k < n; ++k) { 
//					if (!combo_ptrs[k][i]) { gr[j]+=pval_coeff[k]/beta_ptr[i]*beta_deriv[i]; }
//				}
//			}
//		}
//	}
//
//	// Realising that we're dealing with a negative logP-value.
//	for (j=0; j<n; ++j) { 
//		gr[j]/=-pval;
//		gr[j]+=1/total;
//	}
//	return;
//}

//void (*gfp)(int, double*, double*)=NULL;
//void deriv_barrier(int n, double* par, double* gr, void* ex) {
//	(*gfp)(n, par, gr);
//    double bval=0, temp=0, temp2=0;
//    int i, j; 
//    for (i=0; i<n; ++i) { 
//		gr[i] -= mu*(theta_old[i] / par[i] - 1);
//        temp+=par[i];
//        temp2+=theta_old[i];
//    }
//	bval=mu*(1-(1-temp2)/(1-temp));
//	for (i=0; i<n; ++i) { gr[i] -= bval; }
//
//    if (!ISNA(constraints[0])) { 
//        for (i=0; i<(*nexperiments); ++i) {
//            temp=temp2=0;
//            for (j=0; j < n; ++j) {
//                if (combo_ptrs[j][i]) { 
//                    temp+=par[j]; 
//                    temp2+=theta_old[j];
//                }
//            }
//			bval = mu*(1-(constraints[i]-temp2)/(constraints[i]-temp));
//			for (j=0; j<n; ++j) {
//				if (combo_ptrs[j][i]) { gr[j] -= bval; }
//			}
//        }
//    }
//    return;
//}

