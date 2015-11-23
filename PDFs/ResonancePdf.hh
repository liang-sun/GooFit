#ifndef RESONANCE_PDF_HH
#define RESONANCE_PDF_HH

#include "GooPdf.hh" 
#include "devcomplex.hh" 
#include <climits>

#define MAXNKNOBS 1000

enum ResPdfType {RES_RBW = 0, RES_LASS, RES_GS, RES_FLATTE, RES_GAUSS, RES_SPLINE, NONRES};

typedef devcomplex<fptype> (*resonance_function_ptr) (fptype, fptype, fptype, unsigned int*); 

class ResonancePdf : public GooPdf {
  // Service class intended to hold parametrisations of
  // resonances on Dalitz plots. Don't try to use this
  // as a standalone PDF! It should only be used as a
  // component in one of the friend classes. It extends
  // GooPdf so as to take advantage of the 
  // infrastructure, but will crash if used on its own. 

  friend class TddpPdf;
  friend class DalitzPlotPdf; 
  friend class IncoherentSumPdf; 
public:
  // Constructor for regular BW 
  ResonancePdf (string name, 
			  Variable* ar, 
			  Variable* ai, 
			  Variable* mass, 
			  Variable* width, 
			  unsigned int sp, 
			  unsigned int cyc,
              ResPdfType rpt = RES_RBW,
			  const bool symmDP = false
              ); 

  // Gounaris-Sakurai
  ResonancePdf (string name, 
			  Variable* ar, 
			  Variable* ai, 
			  unsigned int sp, 
			  Variable* mass, 
			  Variable* width, 
			  unsigned int cyc, 
              const bool symmDP = false); 
 
  // LASS constructor
  ResonancePdf (string name,
              Variable* ar,
              Variable* ai,
			  Variable* mass,
			  unsigned int sp,
              Variable* width,
              unsigned int cyc,const bool symmDP = false);
  

  // Nonresonant constructor
  ResonancePdf (string name, 
			  Variable* ar, 
			  Variable* ai);  

  // Gaussian constructor
  ResonancePdf (string name,
			  Variable* ar, 
			  Variable* ai,
			  Variable* mean, 
			  Variable* sigma,
			  unsigned int cyc,
              ResPdfType rpt = RES_GAUSS,
              const bool symmDP = false);

  // Flatte constructor (arXiv:1505.01710)
  ResonancePdf (string name,
			  Variable* ar, 
			  Variable* ai,
			  Variable* mean, 
			  Variable* g1,
			  Variable* rg2og1,
			  unsigned int cyc,
              ResPdfType rpt = RES_FLATTE,
              const bool symmDP = false );
  
  // Cubic spline method, a list of bin limits and a list of knobs
  ResonancePdf (string name,
          Variable* ar,
          Variable* ai,
          vector<fptype>& HH_bin_limits,
          vector<Variable*>& pwa_coefs_reals,
          vector<Variable*>& pwa_coefs_imags,
          unsigned int cyc, 
          ResPdfType rpt = RES_SPLINE,
          const bool symmDP = false);

  __host__ void storeParameters () const;// Do something special for PWA

private:
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}

  Variable* amp_real;
  Variable* amp_imag;
  const ResPdfType RPT;
  fptype* host_constants;
  /*
  Variable* mass;
  Variable* width;
  unsigned int spin;
  unsigned int cyclic_index;
  unsigned int eval_type;
  unsigned int resonance_type; 
  */ 
};

#endif
