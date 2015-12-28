#ifndef GIUT_H
#define GIUT_H

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Applic.h"
#include <stdexcept>
#include <cmath>
#include <deque>

/* Assorted variables which mean we don't have to check for 1'ness
 * in the combination matrix every time. Rather, the indices relevant
 * to each row or column are returned. This should shave off some 
 * time during the iteration.
 */ 

struct truth_keeper {
	truth_keeper(const int&, const int &, const int*);
	std::deque<std::deque<int> > by_row, by_column;
	std::deque<std::deque<int> > inv_by_row, inv_by_column;
};

/* A structure that computes the p-value for each set of incoming proportions.
 * It records the alpha and the m-values and computes the constraints. It also
 * records the beta values, the sum of proportions for each experiment, the
 * total sum of proportions.
 */

struct pval_finder {
	pval_finder(const int&, const int&, const int*);
	void set_alpha(const double&); 
	void summate(const double*);
	void compute_pval (const double*);

	const int nprops, nexp;
	const truth_keeper design;
	
	double total;
	std::deque<double> total_per_exp;
	
	std::deque<double> mvals, beta;
	std::deque<double> alphas;

	int i, j;
	double pval, temp;
};

#endif
