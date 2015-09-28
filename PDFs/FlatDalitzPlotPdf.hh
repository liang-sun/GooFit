#ifndef FLATDALITZPLOT_PDF_HH
#define FLATDALITZPLOT_PDF_HH

#include "GooPdf.hh" 
#include "DalitzPlotHelpers.hh" 
#include "devcomplex.hh"

  
class FlatDalitzPlotPdf : public GooPdf {
public:
  FlatDalitzPlotPdf (std::string n, Variable* m12, Variable* m13, DecayInfo* decay);
  // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the 
  // coherent sum. The caching method requires that it be done this way or the ProdPdf
  // normalisation will get *really* confused and give wrong answers. 

  __host__ virtual fptype normalise () const;
protected:

private:
  DecayInfo* decayInfo; 
  Variable* _m12;
  Variable* _m13; 

};

#endif

