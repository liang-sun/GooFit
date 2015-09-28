#include "FlatDalitzPlotPdf.hh"
#include <complex>
using std::complex; 


EXEC_TARGET fptype device_FlatDalitzPlot (fptype* evt, fptype* p, unsigned int* indices) {
  fptype motherMass = functorConstants[indices[1] + 0];
  fptype daug1Mass  = functorConstants[indices[1] + 1];
  fptype daug2Mass  = functorConstants[indices[1] + 2];
  fptype daug3Mass  = functorConstants[indices[1] + 3];

  fptype m12 = evt[indices[2 + indices[0]]];
  fptype m13 = evt[indices[3 + indices[0]]];

  fptype ret = 1.;
  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) ret = 0;
  //printf("%f\t%f\t%f\t%f\t%f\n", m12, m13, motherMass, daug1Mass,ret);
  return ret;

}

MEM_DEVICE device_function_ptr ptr_to_FlatDalitzPlot = device_FlatDalitzPlot; 

__host__ FlatDalitzPlotPdf::FlatDalitzPlotPdf (std::string n, 
							   Variable* m12, 
							   Variable* m13, 
							   DecayInfo* decay 
							   )
  : GooPdf(0, n) 
  , decayInfo(decay)
  , _m12(m12)
  , _m13(m13)
{
  registerObservable(_m12);
  registerObservable(_m13);

  fptype decayConstants[5];
  
  std::vector<unsigned int> pindices;
  pindices.push_back(registerConstants(5)); 
  decayConstants[0] = decayInfo->motherMass;
  decayConstants[1] = decayInfo->daug1Mass;
  decayConstants[2] = decayInfo->daug2Mass;
  decayConstants[3] = decayInfo->daug3Mass;
  decayConstants[4] = decayInfo->meson_radius;
  MEMCPY_TO_SYMBOL(functorConstants, decayConstants, 5*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);  

  GET_FUNCTION_ADDR(ptr_to_FlatDalitzPlot);
  initialise(pindices);

}

__host__ fptype FlatDalitzPlotPdf::normalise () const {
  recursiveSetNormalisation(1); // Not going to normalise efficiency, 
  // so set normalisation factor to 1 so it doesn't get multiplied by zero. 
  // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency, 
  // don't get zeroes through multiplying by the normFactor. 
  MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice);

  fptype ret = (_m12->upperlimit - _m12->lowerlimit)*(_m13->upperlimit - _m13->lowerlimit); // That complex number is a square, so it's fully real
  ret *= 5.58340972222222232e-01;

  host_normalisation[parameters] = 1.0/ret;
  return (fptype) ret; 
}
