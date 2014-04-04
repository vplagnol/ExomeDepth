#include<iostream>
#include<vector>
#include<string>
#include <cmath>

#include <Rinternals.h>
#include <cstdlib>
#include <sstream>

using namespace std;

extern "C" {
  SEXP C_hmm (const SEXP nstates, const SEXP nobs, const SEXP transitions, const SEXP probabilities, const SEXP positions);
}



SEXP C_hmm (const SEXP nstates, const SEXP nobs, const SEXP transitions, const SEXP probabilities, const SEXP positions) {


  const int nstates_c = *INTEGER (nstates);
  const int nobs_c = *INTEGER (nobs);
  const double * trans_c = REAL(transitions);  //reads first column first
  const double * proba_c = REAL(probabilities);
  
  Rprintf("Number of hidden states: %d\n",nstates_c);
  Rprintf("Number of data points: %d\n", nobs_c);

  SEXP myList, final, calls_R;  //output variables
  
  vector<vector<double> > viterbi_prob;
  vector<vector< int> > from_where;
  
  //---------------- starting point in 0 state
  Rprintf("Initializing the HMM\n");
  vector<int> temp_i(nstates_c,-1);
  from_where.push_back( temp_i );

  vector<double> temp(nstates_c,0.);
  temp[0] = 0.;
  for (int i = 1; i < nstates_c; i++) temp[i] = -HUGE_VAL;
  viterbi_prob.push_back ( temp );

  //-----------------------------------------------------------------------------
  //now run the HMM

  for (int i = 1; i != nobs_c; i++) {     
    viterbi_prob.push_back(temp);
    from_where.push_back( temp_i);

    for (int j = 0; j != nstates_c; j++) {         //for each state
      viterbi_prob[i][j] = - HUGE_VAL;

      for (int k = 0; k != nstates_c; k++) {   //state one is coming from
	//double newp = loglikelihood[i][j] + viterbi_prob[i - 1][k] + log(trans_matrix[k][j]);
	double newp = proba_c[j*nobs_c + i] + viterbi_prob[i - 1][k] + log(trans_c[j*nstates_c + k]);
	
	if (newp >  viterbi_prob[i][j] ) {
	  viterbi_prob[i][j] = newp;
	  from_where[i][j] = k;
	}
      }
      
      if (proba_c[j*nobs_c + i] == -HUGE_VAL) {from_where[i][j] = 0;}  //weak stuff, but it deals with me setting some likelihoods to -Inf
    
      //-----------
      if (from_where[i][j] == -1) {  //catching bugs
	for (int k = 0; k != nstates_c; k++) { 
	  Rprintf("bug"); return (myList);
	}
	Rprintf("Error: Value equal to -1\n"); 
	return myList;
      }
      //-------------


    }
  }

  Rprintf("Done with the first step of the HMM, now running the trace back\n");

  //-------------------------------------------- trace back
  vector<int> trace_back ( nobs_c, 0);
  trace_back[ nobs_c - 1] = 0;

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
  Rprintf("Total number of calls: %d\n", ncalls);


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
  
