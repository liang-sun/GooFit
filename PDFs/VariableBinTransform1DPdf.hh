#ifndef VARBINTRANSFORM_PDF_HH
#define VARBINTRANSFORM_PDF_HH

#include "GooPdf.hh" 

// Transforms ND coordinates into a single bin number. 
class VariableBinTransform1DPdf : public GooPdf {
public:
  VariableBinTransform1DPdf (std::string n, Variable* _x, vector<fptype> binlimits/*, vector<fptype> binSizes, vector<int> numBins*/); 

private:

};

#endif
