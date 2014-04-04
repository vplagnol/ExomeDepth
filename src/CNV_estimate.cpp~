#include<iostream>
#include<vector>
#include<string>
#include <cmath>
#include "gsl_sf_gamma.h"

#include <Rinternals.h>
#include <cstdlib>
#include <sstream>
#include <map>

using namespace std;
 
extern "C" {
  SEXP CNV_estimate (const SEXP expected_a, const SEXP total_a, const SEXP observed_a);
  SEXP get_loglike_matrix (const SEXP phi_a, const SEXP expected_a, const SEXP total_a, const SEXP observed_a, const SEXP mixture_a);
  SEXP dbetabin_ab_vp (const SEXP x_a, const SEXP size_a, const SEXP shape1_a, const SEXP shape2_a);
}


SEXP dbetabin_ab_vp (const SEXP x_a, const SEXP size_a, const SEXP shape1_a, const SEXP shape2_a) {
  const double * observed = REAL( x_a );
  const double total = REAL( size_a )[0];
  const double a1 = REAL( shape1_a )[0];
  const double a2 = REAL( shape2_a )[0];

  const unsigned int nx = length( x_a );


  SEXP ab;
  PROTECT(ab = allocVector(REALSXP, nx));
  for (unsigned int i = 0; i != nx; i++) {
    double mans = gsl_sf_lnbeta (a1 + observed[i], a2 + total - observed[i])- gsl_sf_lnbeta (a1, a2);
    REAL(ab)[i] = mans;
  };

  UNPROTECT(1);
  return (ab);

}



double myprob ( const double expected, const double sd_lambda, const int total, const int observed) {
  double a1 = expected*expected*(1-expected)/(sd_lambda*sd_lambda) - expected;
  double a2 = (1-expected)/expected*a1;


  return (gsl_sf_lnbeta (a1 + observed, a2 + total - observed)- gsl_sf_lnbeta (a1, a2));
}

SEXP get_loglike_matrix (const SEXP phi_a, const SEXP expected_a, const SEXP total_a, const SEXP observed_a, const SEXP mixture_a) {

  const double * expected = REAL( expected_a );
  const int * total = INTEGER( total_a );
  const int * observed = INTEGER( observed_a );
  const unsigned int nCNVs = length( total_a );
  const double * phi = REAL( phi_a );
  const double mixture = REAL( mixture_a )[ 0 ];

  if (mixture != 1) {Rprintf("As a warning (this could be normal), the mixture coefficient is %f\n", mixture);}
  

  //----------------------------------------
  double odds_del = 1 - 0.5*mixture;
  double odds_dup = 1 + 0.5*mixture;

  SEXP ans;
  PROTECT(ans = allocMatrix(REALSXP, nCNVs, 3));  //matrix of probabilities for 1,2,3 cn
  double * rans = REAL(ans);
  for (unsigned int cnv = 0; cnv != nCNVs; cnv++) {

    double bestSd = sqrt(phi[ cnv ]*expected[ cnv ]*(1.-expected[ cnv ]));
    
    rans[cnv + nCNVs*0] =   myprob ( expected[ cnv ]*odds_del/ ( expected[ cnv ]*odds_del + 1 - expected[ cnv ] ), bestSd, total[ cnv ], observed[ cnv]);
    rans[cnv + nCNVs*1] =   myprob ( expected[ cnv ], bestSd, total[ cnv ], observed[ cnv]);
    rans[cnv + nCNVs*2] =   myprob ( expected[ cnv ]*odds_dup/ ( expected[ cnv ]*odds_dup + 1 - expected[ cnv ] )  , bestSd, total[ cnv ], observed[ cnv]);
  }
  UNPROTECT(1);
  
  return (ans);
}


SEXP CNV_estimate (const SEXP expected_a, const SEXP total_a, const SEXP observed_a) {


  
  const double * expected = REAL( expected_a );
  const int * total = INTEGER( total_a );
  const int * observed = INTEGER( observed_a );

  //------------ steps for the lambda parameters search
  double my_expected = expected[0];
  double max_sdE = my_expected/5.;
  double min_sdE = 0.001;
  double step = (max_sdE - min_sdE)/20.;


  //-------------
  unsigned int nCNVs = length( total_a );
  Rprintf("Number of CNVs: %d\n", nCNVs);

  double maxLL = - HUGE_VAL, bestSd = 0.;
  
  for (double sd_lambda = min_sdE; sd_lambda < max_sdE; sd_lambda += step) {

    double logLike = 0;
    for (unsigned int i = 0; i != nCNVs; i++) logLike += myprob ( expected[ i ], sd_lambda, total[ i ], observed[ i ]);
    
    if (logLike > maxLL) {
      bestSd = sd_lambda;
      maxLL = logLike;
    }
    
  }
  
  Rprintf("Best fitting sd_lambda: %f\n", bestSd);

  //-------------------------
  SEXP ans, myList, sdLambda;
  PROTECT( myList = allocVector(VECSXP, 2));

  //-------- first I return the estimated sd for lambda
  PROTECT(sdLambda = allocVector(REALSXP,1));           
  double * csdLambda = REAL(sdLambda);
  csdLambda[0] = bestSd;
  SET_VECTOR_ELT(myList, 0, sdLambda);

  
  //---------- Now the matric of likelihood
  PROTECT(ans = allocMatrix(REALSXP, nCNVs, 3));  //matrix of probabilities for 1,2,3 cn
  double * rans = REAL(ans);
  for (unsigned int cnv = 0; cnv != nCNVs; cnv++) {

    rans[cnv + nCNVs*0] =   myprob ( expected[ cnv ]*0.5/ ( expected[ cnv ]*0.5 +  1- expected[ cnv ] ), bestSd, total[ cnv ], observed[ cnv]);
    rans[cnv + nCNVs*1] =   myprob ( expected[ cnv ], bestSd, total[ cnv ], observed[ cnv]);
    rans[cnv + nCNVs*2] =   myprob ( expected[ cnv ]*1.5/ ( expected[ cnv ]*1.5 +  1- expected[ cnv ] )  , bestSd, total[ cnv ], observed[ cnv]);
  }
  SET_VECTOR_ELT(myList, 1, ans);

  UNPROTECT(3);
  
  return (myList);

}





//                        y | lambda ~ Binomial(n, lambda)                
//  P(lambda) = lambda^{a1 - 1} * (1 - lambda)^{a2 - 1} / B(a1, a2)   
//                         E[lambda] = a1 / (a1 + a2)                          
//              Var[lambda] = a1 * a2 / [(a1 + a2 + 1) * (a1 + a2)^2]          
//     The marginal beta-binomial distribution is:
//              P(y) = C(n, y) * B(a1 + y, a2 + n - y) / B(a1, a2)      
//  The function uses the parameterization p = a1 / (a1 + a2) = h(X b)
//     = h(eta) and phi = 1 / (a1 + a2 + 1), where h is the inverse of
//     the link function (logit or complementary log-log),
//                   E[y] = n * p                               
//               Var[y] = n * p * (1 - p) * [1 + (n - 1) * phi]     


//beta distribution: var(lambda) = p(1-p)/(1+a1 + a2)
//alpha = 


//Function: int gsl_sf_lnchoose_e (unsigned int n, unsigned int m, gsl_sf_result * result)
//These routines compute the logarithm of n choose m. This is equivalent to the sum \log(n!) - \log(m!) - \log((n-m)!).

