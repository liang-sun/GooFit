#include "MixingTimeResolution_Aux.hh"
#include "GooPdf.hh" 

MixingTimeResolution::MixingTimeResolution () { dtimeMin = 1; dtimeMax = -1; }
MixingTimeResolution::~MixingTimeResolution () {}

void MixingTimeResolution::initIndex (void* dev_fcn_ptr) {
  resFunctionIdx = GooPdf::findFunctionIdx(dev_fcn_ptr); 
}

void MixingTimeResolution::setRange(Variable* dtvar){
    if (dtvar == NULL) return;
    dtimeMin = dtvar->lowerlimit; dtimeMax = dtvar->upperlimit; 
}
