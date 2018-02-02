#include "ResonancePdf.hh" 

MEM_DEVICE fptype cDeriatives[2*MAXNKNOBS];
//MEM_DEVICE fptype sharedCache[1024*5*MAXNKNOBS];

EXEC_TARGET fptype twoBodyCMmom (fptype rMassSq, fptype d1m, fptype d2m) {
  // For A -> B + C, calculate momentum of B and C in rest frame of A. 
  // PDG 38.16.

  fptype kin1 = 1 - POW(d1m+d2m, 2) / rMassSq;
  if (kin1 >= 0) kin1 = SQRT(kin1);
  else kin1 = 1;
  fptype kin2 = 1 - POW(d1m-d2m, 2) / rMassSq;
  if (kin2 >= 0) kin2 = SQRT(kin2);
  else kin2 = 1; 

  return 0.5*SQRT(rMassSq)*kin1*kin2; 
}

EXEC_TARGET fptype twoBodyCMmom (fptype rMassSq, fptype d1m, fptype d2m, fptype mR) {
  // For A -> B + C, calculate momentum of B and C in rest frame of A. 
  // PDG 38.16.

/*  fptype kin1 = rMassSq - POW(d1m+d2m, 2) ;
  if (kin1 >= 0) kin1 = SQRT(kin1);
  else kin1 = 1;
  fptype kin2 = rMassSq - POW(d1m-d2m, 2) ;
  if (kin2 >= 0) kin2 = SQRT(kin2);
  else kin2 = 1; 

  return 0.5*mR*kin1*kin2/rMassSq; //Equiv. to the original form if mR*mR == rMassSq*/
    fptype x = rMassSq;
    fptype y = d1m*d1m;
    fptype z = d2m*d2m;
	double l = (x - y - z)*(x - y - z) - 4*y*z;
/*	if (l<0){ 
		cout << "lambda - Error - lambda < 0" <<" mr = "<<SQRT(x)<< endl;
                exit(-1);      
	}*/
	return SQRT(l)/(2*mR);    
}

EXEC_TARGET fptype dampingFactorSquare (fptype cmmom, int spin, fptype mRadius) {
  fptype square = mRadius*mRadius*cmmom*cmmom;
  fptype dfsq = 1 + square; // This accounts for spin 1
  if (2 == spin) dfsq += 8 + 2*square + square*square; // Coefficients are 9, 3, 1.   

  // Spin 3 and up not accounted for. 
  return dfsq; 
}

EXEC_TARGET fptype spinFactor (unsigned int spin, fptype motherMass, fptype daug1Mass, fptype daug2Mass, fptype daug3Mass, fptype m12, fptype m13, fptype m23, unsigned int cyclic_index, fptype m02 = 0) {
  if (0 == spin) return 1; // Should not cause branching since every thread evaluates the same resonance at the same time. 
  /*
  // Copied from BdkDMixDalitzAmp
   
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass)); 
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug3Mass)); 
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass)); 
    
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m12 : m12)); 
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m23 : m13)); 
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 

  // The above, collapsed into single tests where possible. 
  fptype _mA = (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass);
  fptype _mB = (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass); 
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass)); 

  fptype _mAC = (PAIR_23 == cyclic_index ? m13 : m23);
  fptype _mBC = (PAIR_12 == cyclic_index ? m13 : m12);
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 
  */

  // Copied from EvtDalitzReso, with assumption that pairAng convention matches pipipi0 from EvtD0mixDalitz.
  // Again, all threads should get the same branch. 
  fptype _mA = (PAIR_12 == cyclic_index ? daug1Mass : (PAIR_13 == cyclic_index ? daug3Mass : daug2Mass));
  fptype _mB = (PAIR_12 == cyclic_index ? daug2Mass : (PAIR_13 == cyclic_index ? daug1Mass : daug3Mass));
  fptype _mC = (PAIR_12 == cyclic_index ? daug3Mass : (PAIR_13 == cyclic_index ? daug2Mass : daug1Mass));
  fptype _mAC = (PAIR_12 == cyclic_index ? m13 : (PAIR_13 == cyclic_index ? m23 : m12)); 
  fptype _mBC = (PAIR_12 == cyclic_index ? m23 : (PAIR_13 == cyclic_index ? m12 : m13)); 
  fptype _mAB = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23)); 

  if (m02 <=0) m02 = _mAB;
  fptype massFactor = 1.0/m02;
  fptype sFactor = -2; 
  sFactor *= ((_mBC - _mAC) + (massFactor*(motherMass*motherMass - _mC*_mC)*(_mA*_mA-_mB*_mB)));
  if (2 == spin) {
    sFactor *= sFactor; 
    fptype extraterm = ((_mAB-(2*motherMass*motherMass)-(2*_mC*_mC))+massFactor*pow((motherMass*motherMass-_mC*_mC),2));
    extraterm *= ((_mAB-(2*_mA*_mA)-(2*_mB*_mB))+massFactor*pow((_mA*_mA-_mB*_mB),2));
    extraterm /= 3;
    sFactor -= -4*extraterm;
  }
  return sFactor; 
}

EXEC_TARGET devcomplex<fptype> plainBW (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  fptype meson_radius           = functorConstants[indices[1]+4];

  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int spin             = indices[4];
  unsigned int cyclic_index     = indices[5];
  unsigned int doSwap           = indices[6]; 

  devcomplex<fptype> ret(0., 0.);
  fptype resmassSq  = resmass*resmass;
  for (int i=0;i<1+doSwap;i++){

  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype rMass = SQRT(rMassSq);  
  fptype frFactor = 1;
  fptype fdFactor = 1;

  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, 
					    (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), 
					    (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass),
                        rMass);
  fptype nominalDaughterMoms = twoBodyCMmom(resmassSq, 
					    (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), 
					    (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius); 
    // Form_Factor_Mother_Decay
    fptype measureDaughterMoms2 = twoBodyCMmom(motherMass*motherMass, 
                        rMass, 
//					    (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), 
					    (PAIR_12 == cyclic_index ? daug3Mass : daug2Mass),
                        rMass);
//                        motherMass);
    fptype nominalDaughterMoms2 = twoBodyCMmom(motherMass*motherMass, 
                        resmass, 
//					    (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), 
					    (PAIR_12 == cyclic_index ? daug3Mass : daug2Mass),
                        resmass);
//                        motherMass);

    fdFactor =  dampingFactorSquare(nominalDaughterMoms2, spin, 5.);
    fdFactor /= dampingFactorSquare(measureDaughterMoms2, spin, 5.); 
  }  
 
  // RBW evaluation
  fptype A = (resmassSq - rMassSq); 
  fptype B = resmassSq*reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor / SQRT(rMassSq);
  fptype C = 1.0 / (A*A + B*B); 
  devcomplex<fptype> retur(A*C, B*C); // Dropping F_D=1

  retur *= SQRT(frFactor*fdFactor); 
  fptype spinF = spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index); 
//  fptype spinF = spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index, resmass); 
  retur *= spinF; 
  ret += retur;
  if (doSwap) {fptype swpmass = m12; m12 = m13; m13 = swpmass;}
  }

  return ret; 
}

EXEC_TARGET devcomplex<fptype> gaussian (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  // indices[1] is unused constant index, for consistency with other function types. 
  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int cyclic_index     = indices[4]; 
  unsigned int doSwap           = indices[5]; 

  fptype ret(0.);
  for (int i=0;i<1+doSwap;i++){

  // Notice SQRT - this function uses mass, not mass-squared like the other resonance types. 
  fptype massToUse = SQRT(PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  massToUse -= resmass;
  massToUse /= reswidth;
  massToUse *= massToUse;
  fptype retur = EXP(-0.5*massToUse); 

  // Ignore factor 1/SQRT(2pi). 
  retur /= reswidth;
  ret += retur;
  if (doSwap) {fptype swpmass = m12; m12 = m13; m13 = swpmass;}
  }

  return devcomplex<fptype>(ret, 0); 
}

EXEC_TARGET fptype hFun (double s, double daug2Mass, double daug3Mass) {
  // Last helper function
  const fptype _pi = 3.14159265359;
  double sm   = daug2Mass + daug3Mass;
  double SQRTs = SQRT(s);
  double k_s = twoBodyCMmom(s, daug2Mass, daug3Mass);

  double val = ((2/_pi) * (k_s/SQRTs) * log( (SQRTs + 2*k_s)/(sm)));

  return val;
}

EXEC_TARGET fptype dh_dsFun (double s, double daug2Mass, double daug3Mass) {
  // Yet another helper function
  const fptype _pi = 3.14159265359;
  double k_s = twoBodyCMmom(s, daug2Mass, daug3Mass);
  
  double val = (hFun(s, daug2Mass, daug3Mass) * (1.0/(8.0*pow(k_s, 2)) - 1.0/(2.0 * s)) + 1.0/(2.0* _pi*s));
  return val;
}


EXEC_TARGET fptype dFun (double s, double daug2Mass, double daug3Mass) {
  // Helper function used in Gronau-Sakurai
  const fptype _pi = 3.14159265359;
  double sm   = daug2Mass + daug3Mass;
  double sm24 = sm*sm/4.0;
  double m    = SQRT(s);
  double k_m2 = twoBodyCMmom(s, daug2Mass, daug3Mass);
 
  double val = 3.0/_pi * sm24/pow(k_m2, 2) * log((m + 2*k_m2)/sm) + m/(2*_pi*k_m2) - sm24*m/(_pi * pow(k_m2, 3));
  return val;
}

EXEC_TARGET fptype fsFun (double s, double m2, double gam, double daug2Mass, double daug3Mass) {
  // Another G-S helper function
   
  double k_s   = twoBodyCMmom(s,  daug2Mass, daug3Mass);
  double k_Am2 = twoBodyCMmom(m2, daug2Mass, daug3Mass);
   
  double f     = gam * m2 / POW(k_Am2, 3);
  f           *= (POW(k_s, 2) * (hFun(s, daug2Mass, daug3Mass) - hFun(m2, daug2Mass, daug3Mass)) + (m2 - s) * pow(k_Am2, 2) * dh_dsFun(m2, daug2Mass, daug3Mass));
 
  return f;
}

EXEC_TARGET devcomplex<fptype> gouSak (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  fptype meson_radius           = functorConstants[indices[1]+4];

  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int spin             = indices[4];
  unsigned int cyclic_index     = indices[5];
  unsigned int doSwap           = indices[6]; 

  resmass *= resmass; 
  devcomplex<fptype> ret(0., 0.);
  for (int i=0;i<1+doSwap;i++){

  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype frFactor = 1;

  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2). 
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_12 == cyclic_index ? daug2Mass : daug3Mass));

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius); 
  }
  
  // Implement Gou-Sak:

  fptype D = (1.0 + dFun(resmass, daug2Mass, daug3Mass) * reswidth/SQRT(resmass));
  fptype E = resmass - rMassSq + fsFun(rMassSq, resmass, reswidth, daug2Mass, daug3Mass);
  fptype F = SQRT(resmass) * reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor;

  D       /= (E*E + F*F);
  devcomplex<fptype> retur(D*E, D*F); // Dropping F_D=1
  retur *= SQRT(frFactor);
  retur *= spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);
  ret += retur;
  if (doSwap) {fptype swpmass = m12; m12 = m13; m13 = swpmass;}
  }
  return ret; 
}


EXEC_TARGET devcomplex<fptype> lass (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  fptype motherMass             = functorConstants[indices[1]+0];
  fptype daug1Mass              = functorConstants[indices[1]+1];
  fptype daug2Mass              = functorConstants[indices[1]+2];
  fptype daug3Mass              = functorConstants[indices[1]+3];
  fptype meson_radius           = functorConstants[indices[1]+4];

  fptype resmass                = cudaArray[indices[2]];
  fptype reswidth               = cudaArray[indices[3]];
  unsigned int spin             = indices[4];
  unsigned int cyclic_index     = indices[5];
  unsigned int doSwap           = indices[6]; 

  fptype _a    = 0.22357;
  fptype _r    = -15.042;
  fptype _R    = 1; // ?
  fptype _phiR = 1.10644;
  fptype _B    = 0.614463;
  fptype _phiB = -0.0981907;
  resmass *= resmass;

  devcomplex<fptype> ret(0., 0.);
  for (int i=0;i<1+doSwap;i++){

  fptype rMassSq = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
  fptype frFactor = 1;

  // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <-> dm2).
  
  fptype measureDaughterMoms = twoBodyCMmom(rMassSq, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_23 == cyclic_index ? daug3Mass : daug2Mass));
  fptype nominalDaughterMoms = twoBodyCMmom(resmass, (PAIR_23 == cyclic_index ? daug2Mass : daug1Mass), (PAIR_23 == cyclic_index ? daug3Mass : daug2Mass));

  if (0 != spin) {
    frFactor =  dampingFactorSquare(nominalDaughterMoms, spin, meson_radius);
    frFactor /= dampingFactorSquare(measureDaughterMoms, spin, meson_radius);
  }

  //Implement LASS:
  /*
  fptype s = kinematics(m12, m13, _trackinfo[i]);
  fptype q = twoBodyCMmom(s, _trackinfo[i]);
  fptype m0  = _massRes[i]->getValFast();
  fptype _g0 = _gammaRes[i]->getValFast();
  int spin   = _spinRes[i];
  fptype g = runningWidthFast(s, m0, _g0, spin, _trackinfo[i], FrEval(s, m0, _trackinfo[i], spin));
  */

  fptype q = measureDaughterMoms;
  fptype g = reswidth * POW(measureDaughterMoms / nominalDaughterMoms, 2.0*spin + 1) * frFactor / SQRT(rMassSq);

  // background phase motion
  fptype cot_deltaB = (1.0 / (_a*q)) + 0.5*_r*q;
  fptype qcot_deltaB = (1.0 / _a) + 0.5*_r*q*q;

  // calculate resonant part
  devcomplex<fptype> expi2deltaB = devcomplex<fptype>(qcot_deltaB,q)/devcomplex<fptype>(qcot_deltaB,-q);
  devcomplex<fptype>  resT = devcomplex<fptype>(cos(_phiR+2*_phiB),sin(_phiR+2*_phiB))*_R;

  devcomplex<fptype> prop = devcomplex<fptype>(1, 0)/devcomplex<fptype>(resmass-rMassSq, SQRT(resmass)*g);
  // resT *= prop*m0*_g0*m0/twoBodyCMmom(m0*m0, _trackinfo[i])*expi2deltaB;
  resT *= prop*(resmass*reswidth/nominalDaughterMoms)*expi2deltaB;

  // calculate bkg part
  resT += devcomplex<fptype>(cos(_phiB),sin(_phiB))*_B*(cos(_phiB)+cot_deltaB*sin(_phiB))*SQRT(rMassSq)/devcomplex<fptype>(qcot_deltaB,-q);

  resT *= SQRT(frFactor);
  resT *= spinFactor(spin, motherMass, daug1Mass, daug2Mass, daug3Mass, m12, m13, m23, cyclic_index);
  ret += resT;
  if (doSwap) {fptype swpmass = m12; m12 = m13; m13 = swpmass;}
  }
  
  return ret;
}

EXEC_TARGET devcomplex<fptype> flatte(fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  // indices[1] is unused constant index, for consistency with other function types. 
  fptype resmass                = cudaArray[indices[2]];
  fptype g1                     = cudaArray[indices[3]];
  fptype g2                     =  (cudaArray[indices[4]])*g1;
  unsigned int cyclic_index     = indices[5]; 
  unsigned int doSwap           = indices[6]; 

  fptype pipmass = 0.13957018;
  fptype pi0mass = 0.1349766;
  fptype kpmass = 0.493677;
  fptype k0mass = 0.497614;

  fptype twopimasssq = 4*pipmass*pipmass;
  fptype twopi0masssq = 4*pi0mass*pi0mass;
  fptype twokmasssq = 4*kpmass*kpmass;
  fptype twok0masssq = 4*k0mass*k0mass;

  devcomplex<fptype> ret(0., 0.);
  for (int i=0;i<1+doSwap;i++){

  fptype rhopipi_real = 0, rhopipi_imag = 0;
  fptype rhokk_real = 0, rhokk_imag = 0;

  fptype s = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));

  if (s>=twopimasssq) rhopipi_real += (2./3)*SQRT(1-twopimasssq/s);//Above pi+pi- threshold
  else rhopipi_imag += (2./3)*SQRT(-1+twopimasssq/s);
  if (s>=twopi0masssq) rhopipi_real += (1./3)*SQRT(1-twopi0masssq/s);//Above pi0pi0 threshold
  else rhopipi_imag += (1./3)*SQRT(-1+twopi0masssq/s);
  if (s>=twokmasssq) rhokk_real += 0.5*SQRT(1-twokmasssq/s);//Above K+K- threshold
  else rhokk_imag += 0.5*SQRT(-1+twokmasssq/s);
  if (s>=twok0masssq) rhokk_real += 0.5*SQRT(1-twok0masssq/s);//Above K0K0 threshold
  else rhokk_imag += 0.5*SQRT(-1+twok0masssq/s);
  fptype A = (resmass*resmass - s) + resmass*(rhopipi_imag*g1+rhokk_imag*g2); 
  fptype B = resmass*(rhopipi_real*g1+rhokk_real*g2);
  fptype C = 1.0 / (A*A + B*B); 
  devcomplex<fptype> retur(A*C, B*C);
  ret += retur;
  if (doSwap) {fptype swpmass = m12; m12 = m13; m13 = swpmass;}
  }

  return ret; 
}

/*EXEC_TARGET*/
__host__ bool Complex_Derivative (unsigned int nKnobs, const fptype* x, const devcomplex<fptype>* y, devcomplex<fptype>* y2) {
    int i,k;
    unsigned int n = nKnobs;
    devcomplex<fptype> *u = new devcomplex<fptype>[n];
    fptype sig,p,qn,un;
    devcomplex<fptype> yp1 = 2.*(y[1] - y[0]) / (x[1] - x[0]);
    devcomplex<fptype> ypn = 2.*(y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);
    assert (y2!=0) ; // Must be a valid pointer

    /* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/
    if(yp1.real > 0.99e30) {
        y2[0].real = 0.;
        u[0].real = 0.;
	}
	else{
		y2[0].real=-0.5;
		u[0].real=(3.0/(x[1]-x[0]))*((y[1].real-y[0].real)/(x[1]-x[0])-yp1.real);
	}
    if(yp1.imag > 0.99e30) {
        y2[0].imag = 0.;
        u[0].imag = 0.;
	}
	else{
		y2[0].imag=-0.5;
		u[0].imag=(3.0/(x[1]-x[0]))*((y[1].imag-y[0].imag)/(x[1]-x[0])-yp1.imag);
	}
	
	
/* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors*/
	
	for(i=1;i<n-1;i++){
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1].real+2.0;
		y2[i].real=(sig-1.0)/p;
		u[i].real=(y[i+1].real-y[i].real)/(x[i+1]-x[i]) - (y[i].real-y[i-1].real)/(x[i]-x[i-1]);
		u[i].real=(6.0*u[i].real/(x[i+1]-x[i-1])-sig*u[i-1].real)/p;
		p=sig*y2[i-1].imag+2.0;
		y2[i].imag=(sig-1.0)/p;
		u[i].imag=(y[i+1].imag-y[i].imag)/(x[i+1]-x[i]) - (y[i].imag-y[i-1].imag)/(x[i]-x[i-1]);
		u[i].imag=(6.0*u[i].imag/(x[i+1]-x[i-1])-sig*u[i-1].imag)/p;
	}
	
	/* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/
	
	if(ypn.real > 0.99e30) {
        qn = 0.;
        un = 0.;
	}
	else{
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn.real-(y[n-1].real-y[n-2].real)/(x[n-1]-x[n-2]));
	}
	y2[n-1].real=(un-qn*u[n-2].real)/(qn*y2[n-2].real+1.0);
	if(ypn.imag > 0.99e30) {
        qn = 0.;
        un = 0.;
	}
	else{
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn.imag-(y[n-1].imag-y[n-2].imag)/(x[n-1]-x[n-2]));
	}
	y2[n-1].imag=(un-qn*u[n-2].imag)/(qn*y2[n-2].imag+1.0);
	
	/* This is the backsubstitution loop of the tridiagonal algorithm */
	
	for(k=n-2;k>=0;k--) {
        y2[k].real=y2[k].real*y2[k+1].real+u[k].real;
        y2[k].imag=y2[k].imag*y2[k+1].imag+u[k].imag;
    }
    delete [] u;
    return true;
}

EXEC_TARGET devcomplex<fptype> cubicspline (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  devcomplex<fptype> ret(0,0);
  unsigned int cyclic_index     = indices[2]; 
  unsigned int doSwap           = indices[3]; 
  const unsigned int nKnobs = indices[4];
  unsigned int idx = 5; // Next index
  unsigned int i = 0;
  const unsigned int pwa_coefs_idx  = idx;
  idx += 2*nKnobs;
  const fptype* mKKlimits = &(functorConstants[indices[idx]]);
  fptype mAB = m12, mAC = m13, mBC = m23;
  switch (cyclic_index) {
    case PAIR_13:  mAB = m13; mAC = m12; break;
    case PAIR_23:  mAB = m23; mAC = m12; mBC = m13; break;
  }
  int khiAB = 0, khiAC = 0;
  fptype dmKK, aa, bb, aa3, bb3;
  unsigned int timestorun = 1+doSwap;
  for (;khiAB<nKnobs;){
      if (mAB<mKKlimits[khiAB]) break;
      khiAB ++;
  }
//  while (khiAB<nKnobs&&mAB>mKKlimits[++khiAB]);
  if (khiAB <=0 || khiAB == nKnobs)  timestorun = 0; 
  for (;khiAC<nKnobs;){
      if (mAC<mKKlimits[khiAC]) break;
      khiAC ++;
  }
//  while (khiAC<nKnobs&&mAC>mKKlimits[++khiAC]);
  if ((khiAC <=0 || khiAC == nKnobs)&&timestorun)  timestorun = 1; //Only run once even if set to be symmetric

  for (i=0;i<timestorun;i++){
  unsigned int kloAB = khiAB -1;//, kloAC = khiAC -1;
  unsigned int twokloAB = kloAB + kloAB;
  unsigned int twokhiAB = khiAB + khiAB;
  fptype pwa_coefs_real_kloAB = cudaArray[indices[pwa_coefs_idx+twokloAB]];
  fptype pwa_coefs_real_khiAB = cudaArray[indices[pwa_coefs_idx+twokhiAB]];
  fptype pwa_coefs_imag_kloAB = cudaArray[indices[pwa_coefs_idx+twokloAB+1]];
  fptype pwa_coefs_imag_khiAB = cudaArray[indices[pwa_coefs_idx+twokhiAB+1]];
  fptype pwa_coefs_prime_real_kloAB = cDeriatives[twokloAB];
  fptype pwa_coefs_prime_real_khiAB = cDeriatives[twokhiAB];
  fptype pwa_coefs_prime_imag_kloAB = cDeriatives[twokloAB+1];
  fptype pwa_coefs_prime_imag_khiAB = cDeriatives[twokhiAB+1];
//  printf("m12: %f: %f %f %f %f %f %f %d %d %d\n", mAB, mKKlimits[0], mKKlimits[nKnobs-1], pwa_coefs_real_khiAB, pwa_coefs_imag_khiAB, pwa_coefs_prime_real_khiAB, pwa_coefs_prime_imag_khiAB, khiAB, khiAC, timestorun );

  dmKK = mKKlimits[khiAB] - mKKlimits[kloAB];
  aa = ( mKKlimits[khiAB] - mAB )/dmKK;
  bb = 1 - aa;
  aa3 = aa * aa * aa; bb3 = bb * bb * bb;
//  ret += aa * pwa_coefs[kloAB] + bb * pwa_coefs[khiAB] + ((aa3 - aa)*pwa_coefs_prime[kloAB] + (bb3 - bb) * pwa_coefs_prime[khiAB]) * (dmKK*dmKK)/6.0;
  ret.real += aa * pwa_coefs_real_kloAB + bb * pwa_coefs_real_khiAB + ((aa3 - aa)*pwa_coefs_prime_real_kloAB + (bb3 - bb) * pwa_coefs_prime_real_khiAB) * (dmKK*dmKK)/6.0;
  ret.imag += aa * pwa_coefs_imag_kloAB + bb * pwa_coefs_imag_khiAB + ((aa3 - aa)*pwa_coefs_prime_imag_kloAB + (bb3 - bb) * pwa_coefs_prime_imag_khiAB) * (dmKK*dmKK)/6.0;
  khiAB = khiAC;  mAB = mAC;
  }
  return ret;
}

EXEC_TARGET devcomplex<fptype> nonres (fptype m12, fptype m13, fptype m23, unsigned int* indices) {
  return devcomplex<fptype>(1, 0); 
}


EXEC_TARGET void getAmplitudeCoefficients (devcomplex<fptype> a1, devcomplex<fptype> a2, fptype& a1sq, fptype& a2sq, fptype& a1a2real, fptype& a1a2imag) {
  // Returns A_1^2, A_2^2, real and imaginary parts of A_1A_2^*
  a1sq = a1.abs2();
  a2sq = a2.abs2();
  a1 *= conj(a2);
  a1a2real = a1.real;
  a1a2imag = a1.imag; 
}

MEM_DEVICE resonance_function_ptr ptr_to_RBW = plainBW;
MEM_DEVICE resonance_function_ptr ptr_to_GOUSAK = gouSak; 
MEM_DEVICE resonance_function_ptr ptr_to_GAUSSIAN = gaussian;
MEM_DEVICE resonance_function_ptr ptr_to_NONRES = nonres;
MEM_DEVICE resonance_function_ptr ptr_to_LASS = lass;
MEM_DEVICE resonance_function_ptr ptr_to_FLATTE = flatte;
MEM_DEVICE resonance_function_ptr ptr_to_SPLINE = cubicspline;

ResonancePdf::ResonancePdf (string name, 
						Variable* ar, 
						Variable* ai, 
						Variable* mass, 
						Variable* width, 
						unsigned int sp, 
						unsigned int cyc, ResPdfType rpt, const bool symmDP) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
  , RPT(rpt)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Making room for index of decay-related constants. Assumption:
  // These are mother mass and three daughter masses in that order.
  // They will be registered by the object that uses this resonance,
  // which will tell this object where to find them by calling setConstantIndex. 

  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width)); 
  pindices.push_back(sp);
  pindices.push_back(cyc); 
  pindices.push_back((unsigned int)symmDP); 

  switch( rpt ){
      case RES_GS: 
          GET_FUNCTION_ADDR(ptr_to_GOUSAK);
          break;
      case RES_LASS: 
          GET_FUNCTION_ADDR(ptr_to_LASS);
          break;
      default:
          GET_FUNCTION_ADDR(ptr_to_RBW);
  }
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name, 
						Variable* ar, 
						Variable* ai, 
						unsigned int sp, 
						Variable* mass, 
						Variable* width, 
						unsigned int cyc, const bool symmDP) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
  , RPT(RES_GS)
{
  // Same as BW except for function pointed to. 
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width)); 
  pindices.push_back(sp);
  pindices.push_back(cyc); 
  pindices.push_back((unsigned int)symmDP); 

  GET_FUNCTION_ADDR(ptr_to_GOUSAK);
  initialise(pindices); 
} 
 
   
ResonancePdf::ResonancePdf (string name,
                                                Variable* ar,
                                                Variable* ai,
						Variable* mass,
                                                unsigned int sp,
                                                Variable* width,
                                                unsigned int cyc,
                                                const bool symmDP)
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
  , RPT(RES_LASS)
{
  // Same as BW except for function pointed to.
  vector<unsigned int> pindices;
  pindices.push_back(0);
  pindices.push_back(registerParameter(mass));
  pindices.push_back(registerParameter(width));
  pindices.push_back(sp);
  pindices.push_back(cyc);
  pindices.push_back((unsigned int)symmDP);

  GET_FUNCTION_ADDR(ptr_to_LASS);
  initialise(pindices);
}


ResonancePdf::ResonancePdf (string name, 
						Variable* ar, 
						Variable* ai ) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
  , RPT(NONRES)
{
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Dummy index for constants - won't use it, but calling 
  // functions can't know that and will call setConstantIndex anyway. 
  GET_FUNCTION_ADDR(ptr_to_NONRES);
  initialise(pindices); 
}

ResonancePdf::ResonancePdf (string name,
						Variable* ar, 
						Variable* ai,
						Variable* mean, 
						Variable* sigma,
						unsigned int cyc, ResPdfType rpt, const bool symmDP ) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
  , RPT(rpt)
{
  assert(rpt == RES_GAUSS);
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  // Dummy index for constants - won't use it, but calling 
  // functions can't know that and will call setConstantIndex anyway. 
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma)); 
  pindices.push_back(cyc); 
  pindices.push_back((unsigned int)symmDP); 

  GET_FUNCTION_ADDR(ptr_to_GAUSSIAN);
  initialise(pindices); 

}

ResonancePdf::ResonancePdf (string name,
						Variable* ar, 
						Variable* ai,
						Variable* mean, 
						Variable* g1,
						Variable* rg2og1,
						unsigned int cyc, ResPdfType rpt, const bool symmDP) 
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
  , RPT(rpt)
{
  assert(rpt == RES_FLATTE);
  vector<unsigned int> pindices; 
  pindices.push_back(0); 
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(g1)); 
  pindices.push_back(registerParameter(rg2og1)); 
  pindices.push_back(cyc); 

  GET_FUNCTION_ADDR(ptr_to_FLATTE);
  initialise(pindices); 

}
   
ResonancePdf::ResonancePdf (string name,
                                                Variable* ar,
                                                Variable* ai,
                                                vector<fptype>& HH_bin_limits,
                                                vector<Variable*>& pwa_coefs_reals,
                                                vector<Variable*>& pwa_coefs_imags,
                                                unsigned int cyc, ResPdfType rpt, const bool symmDP)
  : GooPdf(0, name)
  , amp_real(ar)
  , amp_imag(ai)
  , RPT(rpt)
{
  assert(rpt == RES_SPLINE);
  vector<unsigned int> pindices;
  const unsigned int nKnobs = HH_bin_limits.size();
  host_constants = new fptype[nKnobs];
  pindices.push_back(0);
  pindices.push_back(cyc);
  pindices.push_back((unsigned int)symmDP);
  pindices.push_back(nKnobs);
  for (int i=0;i<pwa_coefs_reals.size();i++){
      host_constants[i] = HH_bin_limits[i];
      pindices.push_back(registerParameter(pwa_coefs_reals[i]));
      pindices.push_back(registerParameter(pwa_coefs_imags[i]));
  }
  pindices.push_back(registerConstants(nKnobs));
  MEMCPY_TO_SYMBOL(functorConstants, host_constants, nKnobs*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice); 

  GET_FUNCTION_ADDR(ptr_to_SPLINE);

  initialise(pindices);

/*  cachedDerivatives = new DEVICE_VECTOR<devcomplex<fptype> >(nKnobs);
  void* dummy = thrust::raw_pointer_cast(cachedDerivatives->data());
  MEMCPY_TO_SYMBOL(PWA_COEF_PRIME, &dummy, sizeof(devcomplex<fptype>*), 0, cudaMemcpyHostToDevice);*/
//  delete [] host_constants;
}

__host__ void ResonancePdf::storeParameters () const {
    PdfBase::storeParameters ();
    if (RPT == RES_SPLINE){
        parCont params;
        getParameters(params); 
        const unsigned nKnobs = params.size()/2;
        devcomplex<fptype>* y = new devcomplex<fptype>[nKnobs];
        devcomplex<fptype>* y2 = new devcomplex<fptype>[nKnobs];
        unsigned int i = 0; 
        for (parIter v = params.begin(); v != params.end(); ++v,++i) {
             unsigned int idx = i/2;
             fptype value = host_params[(*v)->index];
             if (i%2 == 0) y[idx].real = value;
             else y[idx].imag = value;
        }
        Complex_Derivative (nKnobs, host_constants, y, y2);
        fptype* y2_flat = new fptype[2*nKnobs];
        for (i=0;i<nKnobs;i++) { y2_flat[i+i] = y2[i].real; y2_flat[i+i+1] = y2[i].imag; }
        MEMCPY_TO_SYMBOL(cDeriatives, y2_flat, 2*nKnobs*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
        delete [] y;
        delete [] y2;
        delete [] y2_flat;
    }
}
