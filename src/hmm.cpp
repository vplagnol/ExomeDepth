#include<iostream>
#include<vector>
#include<string>
#include <cmath>

#include <Rinternals.h>
#include <cstdlib>
#include <sstream>

using namespace std;

extern "C" {
  SEXP C_hmm (const SEXP nstates, const SEXP nobs, const SEXP transitions, const SEXP probabilities, const SEXP positions, const SEXP expectedLength);
}



SEXP C_hmm (const SEXP nstates, const SEXP nobs, const SEXP transitions, const SEXP probabilities, const SEXP positions, const SEXP expectedLength) {


  const int nstates_c = *INTEGER (nstates);
  const int nobs_c = *INTEGER (nobs);
  const double * trans_c = REAL(transitions);  //reads first column first
  const double * proba_c = REAL(probabilities);

  
  const int *locations = INTEGER (positions);
  const double Expected = *REAL(expectedLength);

  vector<double> trans (3, 0.);
  double dist, dist_effect;

  SEXP myList = NULL;  //output variables
  SEXP final, calls_R;  //output variables

  // check that we have 3 states
  if (nstates_c != 3) {
    Rprintf("ERROR: The code must assume 3 states\n");
    return(myList);
  }
  
  vector<vector<double> > viterbi_prob;
  vector<vector< int> > from_where;
  
  //---------------- starting point in 0 state
  vector<int> temp_i(nstates_c,-1);
  from_where.push_back( temp_i );

  vector<double> temp(nstates_c,0.);
  temp[0] = 0.;
  for (int i = 1; i < nstates_c; i++) temp[i] = -HUGE_VAL;
  viterbi_prob.push_back ( temp );

  //-----------------------------------------------------------------------------
  //now run the HMM
  
  //suppose 10 observations, then nobs_c = 10, i goes to 9, we should be fine
  for (int i = 1; i != nobs_c; i++) {     
    viterbi_prob.push_back(temp);
    from_where.push_back( temp_i);

    dist = double(locations[ i ]) - double(locations[ i - 1 ]);
    //cout<<locations[i]<<"\t"<<locations[ i -1 ]<<"   ----   "<<dist<<"   ------------ "<<Expected<<endl;
    dist_effect = exp(-dist/Expected);


    for (int j = 0; j != nstates_c; j++) {         //for each state
      viterbi_prob[i][j] = - HUGE_VAL;

      // transition probabilities, note that we assume 3 states
      // note that the code is set up such that if dist is very large
      // dist_effect is very small and we are more likely to trqnsition
      // into a normal copy number state
      trans[0] = trans_c[j*nstates_c];
      trans[1] = dist_effect*trans_c[j*nstates_c + 1 ] + (1.0 - dist_effect)*trans_c[j*nstates_c];
      trans[2] = dist_effect*trans_c[j*nstates_c + 2 ] + (1.0 - dist_effect)*trans_c[j*nstates_c];
      
      for (int k = 0; k != nstates_c; k++) {   //state one is coming from: k means we come from k
	double newp = proba_c[j*nobs_c + i] + viterbi_prob[i - 1][k] + log(trans[ k ]);

	if (newp >  viterbi_prob[i][j] ) {
	  viterbi_prob[i][j] = newp;
	  from_where[i][j] = k;
	}
      }
      
      if (proba_c[j*nobs_c + i] == -HUGE_VAL) {from_where[i][j] = 0;}  //weak stuff, but it deals with me setting some likelihoods to -Inf
    
    }
  }

  //Rprintf("Done with the first step of the HMM, now running the trace back\n");

  //-------------------------------------------- trace back
  vector<int> trace_back ( nobs_c, 0);
  trace_back[ nobs_c - 1] = 0;  //assumes we finish at 0- odd but OK

  for (int i = 1; i != nobs_c; i++) {
    trace_back[ nobs_c - i - 1] = from_where[ nobs_c - i ] [trace_back[ nobs_c - i ] ];
  }


  //-----------------------------------------------------------------------------
  vector<vector<double> > summary;
  vector<double> temp_s(4, 0.);
  double start = -1., end = -1., nexons = 0;

  int current = 0;
  for (int i = 1; i != nobs_c; i++) {

    if (trace_back[i - 1 ] != trace_back[i]) {
      if (current == 0) start = i; 
      if (current != 0) {
	end = i-1;
	temp_s[0] = start + 1;  //because R starts at 1 not 0
	temp_s[1] = end + 1;  //because R starts at 1 not 0
	temp_s[2] = current;
	temp_s[3] = nexons;
	summary.push_back(temp_s);
	nexons = 0;
      }
    }

    if ( trace_back[i] != 0) nexons++;
    current =  trace_back[i];
  }

  int ncalls = summary.size();
  //Rprintf("Total number of calls: %d\n", ncalls);


  //-----------------------------------------------------------------------------
  PROTECT( myList = allocVector(VECSXP, 2));
  PROTECT(final = allocVector(REALSXP,nobs_c)); 
  PROTECT(calls_R = allocMatrix(REALSXP, summary.size(), 4));



  double * cfinal = REAL(final);
  for (int i = 0; i != nobs_c; i++) cfinal[ i ] = trace_back[ i ];
  SET_VECTOR_ELT(myList, 0, final);
  

  double * Rcalls = REAL(calls_R);
  for (int i = 0; i != ncalls; i++) {
    for (int j = 0; j != 4; j++) {
      Rcalls[ ncalls*j + i ] = summary[i][j];
    }
  }

  //SEXP dimnames;
  //PROTECT(dimnames = allocVector(VECSXP, 4));
  //SET_VECTOR_ELT(dimnames, 0, "t1");
  //SET_VECTOR_ELT(dimnames, 1, "t2");
  //SET_VECTOR_ELT(dimnames, 2, "t3");
  //SET_VECTOR_ELT(dimnames, 3, "t4");
  //setAttrib(calls_R, R_DimNamesSymbol, dimnames);
     

  SET_VECTOR_ELT(myList, 1, calls_R);


  UNPROTECT(3);
  return myList;


}
  
