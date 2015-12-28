#include "utils.h"

truth_keeper::truth_keeper(const int& nr, const int& nc, const int* ptr) : by_row(nr), by_column(nc), inv_by_column(nc), inv_by_row(nr) {
	int base=0, c, r;
	for (c=0; c<nc; ++c) {
		for (r=0; r<nr; ++r, ++base) {
			if (ptr[base]) { by_column[c].push_back(r); }
			else { inv_by_column[c].push_back(r); }
		}
	}		
	for (r=0; r<nr; ++r) {
		base=r;
		for (c=0; c<nc; ++c, base+=nr) {
			if (ptr[base]) { by_row[r].push_back(c); }
			else { inv_by_row[r].push_back(c); }
		}
	}
	return;
}

pval_finder::pval_finder(const int& ne, const int& np, const int* xptr) : total_per_exp(ne), nexp(ne), nprops(np), mvals(ne), beta(ne), alphas(ne+1,1),
	design(ne, np, xptr) {};

void pval_finder::set_alpha(const double& a) 
/* Sets up the type I error rate threshold, and computes
 * the powers thereof to save time later.
 */
{
	alphas[1]=a;
	for (i=2; i<=nexp; ++i) { alphas[i]=alphas[i-1]*a; }
	return;
}

void pval_finder::summate(const double* props) 
/* Computing the total proportions for each experiment (i.e., total truth nulls,
 * for use in computing the type II error rate) as well as the total sum of 
 * proportions for all proportions.
 */
{
	for (i=0; i<nexp; ++i) {
		double& cursum=(total_per_exp[i]=0);
		const std::deque<int>& current=design.by_row[i];
		for (j=0; j < current.size(); ++j) { cursum+=props[current[j]]; }
	}
	total=0;
	for (j=0; j<nprops; ++j) { total += props[j]; }
	return;
}

void pval_finder::compute_pval(const double* props) 
/* Computes the IUT p-value. Summation is necessary to obtain updated totals,
 * depending on what is entering as 'props'. Note that the beta deque actually 
 * holds 1-beta.
 */
{
	summate(props);
	for (i=0; i<nexp; ++i) { beta[i]=(1-mvals[i]-alphas[1]*total_per_exp[i])/(1-total_per_exp[i]); }

	pval=0;
	for (j=0; j<nprops; ++j) {
		temp=props[j]*alphas[design.by_column[j].size()];
		const std::deque<int>& current=design.inv_by_column[j];
		for (i=0; i<current.size(); ++i) { temp*=beta[current[i]]; }
		pval += temp;
	}
	pval /= total;
	return;
}
