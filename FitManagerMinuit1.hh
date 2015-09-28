#ifndef FITMANAGER_MINUIT1_HH
#define FITMANAGER_MINUIT1_HH

#include "TMinuit.hh" 
extern PdfBase* pdfPointer; 
extern int numPars; 

#include "TRandom3.hh"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

#ifdef OMP_ON
#pragma omp threadprivate (numPars)
#pragma omp threadprivate (pdfPointer)
#endif

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag); 

class FitManager { 
public:
  FitManager (PdfBase* dat);
  ~FitManager ();
  void setMaxCalls (double mxc) {overrideCallLimit = mxc;}
  void setupMinuit ();
  void runMigrad (); 
  void fit (); 
  TMinuit* getMinuitObject () {return minuit;} 
  void getMinuitValues () const;
  void setRandMinuitValues ( const int nSamples) ;
  void loadSample(const int iSample);
  TMinuit* minuit; 
private:
  double overrideCallLimit; 
  TRandom3 rnd;
  MatrixXd* sqrtCov;
  vector<VectorXd> samples;
};

#endif 
