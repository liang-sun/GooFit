// ROOT stuff
#include "TRandom.hh"
#include "TCanvas.h" 
#include "TFile.h" 
#include "TH1F.h" 
#include "TH2F.h" 
#include "TStyle.h" 
#include "TRandom3.hh" 
#include "TLegend.h" 
#include "TText.h" 
#include "TLine.h" 

// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include "Variable.hh" 
#include "PolynomialPdf.hh" 
#include "DalitzPlotPdf.hh" 
#include "DalitzVetoPdf.hh" 
#include "ResonancePdf.hh" 
#include "AddPdf.hh"
#include "ProdPdf.hh"
#include "GooPdf.hh" 
#include "FitManager.hh" 
#include "UnbinnedDataSet.hh"
#include "EventWeightedAddPdf.hh"
#include "SmoothHistogramPdf.hh"

using namespace std;

TCanvas* foo; 
TCanvas* foodal; 
timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 
UnbinnedDataSet* data = 0; 
TH2F* weightHistogram = 0; 
TH2F* underlyingBins = 0; 

Variable* m12 = 0;
Variable* m13 = 0;
Variable* eventNumber = 0; 
Variable* wSig0 = 0; 
bool fitMasses = false; 
Variable* fixedRhoMass  = new Variable("rho_mass", 0.7758, 0.01, 0.7, 0.8);
Variable* fixedRhoWidth = new Variable("rho_width", 0.1503, 0.01, 0.1, 0.2); 

const fptype _mD0 = 1.86484; 
const fptype _mD02 = _mD0 *_mD0;
const fptype _mD02inv = 1./_mD02; 
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;

// Constants used in more than one PDF component. 
Variable* motherM = new Variable("motherM", _mD0);
Variable* chargeM = new Variable("chargeM", piPlusMass);
Variable* neutrlM = new Variable("neutrlM", piZeroMass);
Variable* massSum = new Variable("massSum", _mD0*_mD0 + 2*piPlusMass*piPlusMass + piZeroMass*piZeroMass); // = 3.53481 
Variable* veto_min = new Variable("veto_min", 0.475*0.475);
Variable* veto_max = new Variable("veto_min", 0.505*0.505);
Variable* constantOne = new Variable("constantOne", 1); 
Variable* constantZero = new Variable("constantZero", 0); 
  
std::vector<PdfBase*> comps;

GooPdf* kzero_veto = 0; 
char strbuffer[1000]; 
double mesonRad  = 1.5;
DalitzPlotPdf* signalDalitz; 
bool doEffSwap = false;
bool saveEffPlot = true;

fptype cpuGetM23 (fptype massPZ, fptype massPM) {
  return (_mD02 + piZeroMass*piZeroMass + piPlusMass*piPlusMass + piPlusMass*piPlusMass - massPZ - massPM); 
}

bool cpuDalitz (fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3) {
  if (m12 < POW(dm1 + dm2, 2)) return false; // This m12 cannot exist, it's less than the square of the (1,2) particle mass.
  if (m12 > POW(bigM - dm3, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter. 
  
  // Calculate energies of 1 and 3 particles in m12 rest frame. 
  fptype e1star = 0.5 * (m12 - dm2*dm2 + dm1*dm1) / SQRT(m12); 
  fptype e3star = 0.5 * (bigM*bigM - m12 - dm3*dm3) / SQRT(m12); 

  // Bounds for m13 at this value of m12.
  fptype minimum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) + SQRT(e3star*e3star - dm3*dm3), 2);
  if (m13 < minimum) return false;
  fptype maximum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) - SQRT(e3star*e3star - dm3*dm3), 2);
  if (m13 > maximum) return false;

  return true; 
}

void getToyData (std::string toyFileName, double sigweight = 0.9) {
  TH2F dalitzplot("dalitzplot", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit); 
  std::vector<Variable*> vars;
  vars.push_back(m12);
  vars.push_back(m13);
  vars.push_back(eventNumber); 
  vars.push_back(wSig0); 
  data = new UnbinnedDataSet(vars); 

  std::ifstream reader;
  reader.open(toyFileName.c_str()); 
  std::string buffer;
  while (!reader.eof()) {
    reader >> buffer;
    if (buffer == "====") break; 
    std::cout << buffer; 
  }

  double dummy = 0; 
  wSig0->value = sigweight;
  while (!reader.eof()) {
    reader >> dummy;
    reader >> dummy;      // m23, m(pi+ pi-), called m12 in processToyRoot convention. 
    reader >> m12->value; // Already swapped according to D* charge. m12 = m(pi+pi0)
    reader >> m13->value;

    // Errors on Dalitz variables
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 

    reader >> dummy; // Decay time
    reader >> dummy; // sigma_t

    reader >> dummy; // Md0
    reader >> dummy; // deltaM
    reader >> dummy; // ProbSig
    reader >> dummy; // Dst charge
    reader >> dummy; // Run
    reader >> dummy; // Event
    reader >> dummy; // Signal and four bkg fractions. 
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 
    reader >> dummy; 

    // EXERCISE 1 (preliminary): Impose an artificial reconstruction efficiency
    // by throwing out events with a probability linear in m12. 
    // NB! This exercise continues below. 

    // EXERCISE 2: Instead of the above efficiency, impose a 
    // K0 veto, by throwing out events with 0.475 < m23 < 0.505. 
    double M_23 = cpuGetM23(m12->value, m13->value);
    if (M_23>veto_min->value && M_23<veto_max->value) continue;

    // EXERCISE 3: Use both the above. 

    eventNumber->value = data->getNumEvents(); 
    data->addEvent(); 

    dalitzplot.Fill(m12->value, m13->value); 
  }
  reader.close();
  
  TRandom3 donram(0); 
  int nsig = data->getNumEvents();
  // Generate background events based on a flat distribution across DP plane
  for (int ib = 0; ib < nsig*(1-sigweight)/sigweight; ib++){
    do{
    m12->value = donram.Uniform(m12->lowerlimit, m12->upperlimit);
    m13->value = donram.Uniform(m13->lowerlimit, m13->upperlimit);
    }while(!cpuDalitz(m12->value, m13->value, _mD0, piZeroMass, piPlusMass, piPlusMass));
    double M_23 = cpuGetM23(m12->value, m13->value);
    if (M_23>veto_min->value && M_23<veto_max->value) {ib--; continue;}
    eventNumber->value = data->getNumEvents(); 
    data->addEvent(); 
    dalitzplot.Fill(m12->value, m13->value); 
  }

  dalitzplot.SetStats(false); 
  dalitzplot.Draw("colz");
  foodal->SaveAs("dalitzplot.png"); 
}

//Taking care of bin blocks around the DP boundary
void createWeightHistogram () {
  weightHistogram = new TH2F("weightHistogram", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit);
  weightHistogram->SetStats(false); 
  double step12 = (m12->upperlimit - m12->lowerlimit) / m12->numbins;
  double step13 = (m13->upperlimit - m13->lowerlimit) / m13->numbins;

  for (int i = 1; i <= m12->numbins; ++i) {
    for (int j = 1; j < m13->numbins; ++j) {
      double maxCount = 0;
      double count = 0; 
      for (double currM12 = m12->lowerlimit + step12*(i-1) + 0.05*step12; currM12 < m12->lowerlimit + step12*i; currM12 += 0.1*step12) {
	for (double currM13 = m13->lowerlimit + step13*(j-1) + 0.05*step13; currM13 < m13->lowerlimit + step13*j; currM13 += 0.1*step13) {
	  maxCount++;
	  if (!cpuDalitz(currM12, currM13, _mD0, piZeroMass, piPlusMass, piPlusMass)) continue;
	  count++; 
	}
      }
      if (0.1 > maxCount) continue;
      count /= maxCount;
      weightHistogram->SetBinContent(i, j, count); 
    }
  }
  //Histogram used to draw eff. plot
  underlyingBins = new TH2F("underlyingBins", "", 
			    m12->numbins, m12->lowerlimit, m12->upperlimit, 
			    m13->numbins, m13->lowerlimit, m13->upperlimit);
  underlyingBins->SetStats(false); 
}

GooPdf* makeEfficiencyPdf () {
  vector<Variable*> lvars;
  lvars.push_back(m12); 
  lvars.push_back(m13);  
  BinnedDataSet* binEffData = new BinnedDataSet(lvars); 
  createWeightHistogram();
  // Now testing your efficiency data by uniformly generating m12,m13 values 
  TRandom3 donram(0); 
  for (int i = 0; i < 1e5; i++){
    do{
    m12->value = donram.Uniform(m12->lowerlimit, m12->upperlimit);
    m13->value = donram.Uniform(m13->lowerlimit, m13->upperlimit);
    }while(!cpuDalitz(m12->value, m13->value, _mD0, piZeroMass, piPlusMass, piPlusMass));
    //Weight will not be one if the physics boundary crosses the bin square.
    double weight = weightHistogram->GetBinContent(weightHistogram->FindBin(m12->value, m13->value));
    binEffData->addWeightedEvent(weight);
    if (underlyingBins) underlyingBins->Fill(m12->value, m13->value, weight);
    // Imposing the requirement on efficiency symmetry along m12=m13 line
      if (doEffSwap){
      double swapmass = m12->value; m12->value = m13->value; m13->value = swapmass;
      weight = weightHistogram->GetBinContent(weightHistogram->FindBin(m12->value, m13->value));
      binEffData->addWeightedEvent(weight);
      if (underlyingBins) underlyingBins->Fill(m12->value, m13->value, weight);
      swapmass = m12->value; m12->value = m13->value; m13->value = swapmass;   
      }
  }
  if (saveEffPlot) {
    foodal->cd();
    underlyingBins->Draw("colz"); 
    foodal->SaveAs("plots/efficiency_bins.png");
    foodal->SetLogz(true);
    foodal->SaveAs("plots/efficiency_bins_log.png");
    foo->cd(); 
  }

  //Variable* effSmoothing = new Variable("effSmoothing", 1.0, 0.1, 0, 1.25); 
  Variable* effSmoothing = new Variable("effSmoothing", 0);   
  SmoothHistogramPdf* ret = new SmoothHistogramPdf("efficiency", binEffData, effSmoothing); 

  return ret; 
}

GooPdf* makeKzeroVeto () {
  if (kzero_veto) return kzero_veto; 

  VetoInfo* kVetoInfo = new VetoInfo();
  kVetoInfo->cyclic_index = PAIR_23; 
  kVetoInfo->minimum = veto_min;
  kVetoInfo->maximum = veto_max;
  vector<VetoInfo*> vetos; vetos.push_back(kVetoInfo); 
  kzero_veto = new DalitzVetoPdf("kzero_veto", m12, m13, motherM, neutrlM, chargeM, chargeM, vetos); 
  return kzero_veto;
}

GooPdf* makeFlatBkgDalitzPdf() {
  VetoInfo* kVetoInfo = new VetoInfo();
  kVetoInfo->cyclic_index = PAIR_23; 
  kVetoInfo->minimum = veto_min;
  kVetoInfo->maximum = veto_max;
  vector<VetoInfo*> vetos; vetos.push_back(kVetoInfo); 
  GooPdf *ret = new DalitzVetoPdf("flatbkgPdf", m12, m13, motherM, neutrlM, chargeM, chargeM, vetos); 
  return ret;
}

DalitzPlotPdf* makeSignalPdf (GooPdf* eff = 0) {
  DecayInfo* dtop0pp = new DecayInfo();
  dtop0pp->motherMass  = _mD0; 
  dtop0pp->daug1Mass  = piZeroMass;
  dtop0pp->daug2Mass  = piPlusMass;
  dtop0pp->daug3Mass  = piPlusMass;
  dtop0pp->meson_radius  = 1.5; 
 
  ResonancePdf* rhop  = new ResonancePdf("rhop",
							     new Variable("rhop_amp_real", 1),
							     new Variable("rhop_amp_imag", 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_12);


  bool fixAmps = false;

  ResonancePdf* rhom  = new ResonancePdf("rhom", 
							     fixAmps ? new Variable("rhom_amp_real", 0.714) : 
							     new Variable("rhom_amp_real",  0.714, 0.001, 0, 0),
							     fixAmps ? new Variable("rhom_amp_imag", -0.025) :
							     new Variable("rhom_amp_imag", -0.025, 0.1, 0, 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_13);

  ResonancePdf* rho0  = new ResonancePdf("rho0", 
							     fixAmps ? new Variable("rho0_amp_real", 0.565) : 
							     new Variable("rho0_amp_real", 0.565, 0.001, 0, 0),
							     fixAmps ? new Variable("rho0_amp_imag", 0.164) :
							     new Variable("rho0_amp_imag", 0.164, 0.1, 0, 0),
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_23);

  Variable* sharedMass = new Variable("rhop_1450_mass", 1.465, 0.01, 1.0, 2.0);
  Variable* shareWidth = new Variable("rhop_1450_width", 0.400, 0.01, 0.01, 5.0); 

  ResonancePdf* rhop_1450  = new ResonancePdf("rhop_1450", 
								  fixAmps ? new Variable("rhop_1450_amp_real", -0.174) : 
								  new Variable("rhop_1450_amp_real", -0.174, 0.001, 0, 0),
								  fixAmps ? new Variable("rhop_1450_amp_imag", -0.117) :
								  new Variable("rhop_1450_amp_imag", -0.117, 0.1, 0, 0),
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_12);

  ResonancePdf* rho0_1450  = new ResonancePdf("rho0_1450", 
								  fixAmps ? new Variable("rho0_1450_amp_real", 0.325) : 
								  new Variable("rho0_1450_amp_real", 0.325, 0.001, 0, 0),
								  fixAmps ? new Variable("rho0_1450_amp_imag", 0.057) : 
								  new Variable("rho0_1450_amp_imag", 0.057, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_23);

  ResonancePdf* rhom_1450  = new ResonancePdf("rhom_1450", 
								  fixAmps ? new Variable("rhom_1450_amp_real", 0.788) : 
								  new Variable("rhom_1450_amp_real", 0.788, 0.001, 0, 0),
								  fixAmps ? new Variable("rhom_1450_amp_imag", 0.226) : 
								  new Variable("rhom_1450_amp_imag", 0.226, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_13);

  sharedMass = new Variable("rhop_1700_mass",  1.720, 0.01, 1.6, 1.9);
  shareWidth = new Variable("rhop_1700_width", 0.250, 0.01, 0.1, 1.0); 

  
  ResonancePdf* rhop_1700  = new ResonancePdf("rhop_1700", 
								  fixAmps ? new Variable("rhop_1700_amp_real", 2.151) : 
								  new Variable("rhop_1700_amp_real",  2.151, 0.001, 0, 0),
								  fixAmps ? new Variable("rhop_1700_amp_imag", -0.658) : 
								  new Variable("rhop_1700_amp_imag", -0.658, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_12);
  
  ResonancePdf* rho0_1700  = new ResonancePdf("rho0_1700", 
								  fixAmps ? new Variable("rho0_1700_amp_real",  2.400) : 
								  new Variable("rho0_1700_amp_real",  2.400, 0.001, 0, 0),
								  fixAmps ? new Variable("rho0_1700_amp_imag", -0.734) : 
								  new Variable("rho0_1700_amp_imag", -0.734, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_23);
  
  ResonancePdf* rhom_1700  = new ResonancePdf("rhom_1700", 
								  fixAmps ? new Variable("rhom_1700_amp_real",  1.286) : 
								  new Variable("rhom_1700_amp_real",  1.286, 0.001, 0, 0),
								  fixAmps ? new Variable("rhom_1700_amp_imag", -1.532) : 
								  new Variable("rhom_1700_amp_imag", -1.532, 0.1, 0, 0),  
								  sharedMass,
								  shareWidth,
								  1,
								  PAIR_13);
  
  ResonancePdf* f0_980  = new ResonancePdf("f0_980", 
							       fixAmps ? new Variable("f0_980_amp_real",  0.008 * (-_mD02)) : 
							       new Variable("f0_980_amp_real",  0.008 * (-_mD02), 0.001, 0, 0),
							       fixAmps ? new Variable("f0_980_amp_imag", -0.013 * (-_mD02)) : 
							       new Variable("f0_980_amp_imag", -0.013 * (-_mD02), 0.1, 0, 0),  
							       new Variable("f0_980_mass",     0.980, 0.01, 0.8, 1.2),
							       new Variable("f0_980_width",    0.044, 0.001, 0.001, 0.08),
							       (unsigned int)0,
							       PAIR_23);
  
  ResonancePdf* f0_1370  = new ResonancePdf("f0_1370", 
								fixAmps ? new Variable("f0_1370_amp_real", -0.058 * (-_mD02)) : 
								new Variable("f0_1370_amp_real", -0.058 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1370_amp_imag",  0.026 * (-_mD02)) : 
								new Variable("f0_1370_amp_imag",  0.026 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1370_mass",     1.434, 0.01, 1.2, 1.6),
								new Variable("f0_1370_width",    0.173, 0.01, 0.01, 0.4),
							    (unsigned int)0,
								PAIR_23);
  
  ResonancePdf* f0_1500  = new ResonancePdf("f0_1500", 
								fixAmps ? new Variable("f0_1500_amp_real", 0.057 * (-_mD02)) : 
								new Variable("f0_1500_amp_real", 0.057 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1500_amp_imag", 0.012 * (-_mD02)) : 
								new Variable("f0_1500_amp_imag", 0.012 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1500_mass",     1.507, 0.01, 1.3, 1.7),
								new Variable("f0_1500_width",    0.109, 0.01, 0.01, 0.3),
							    (unsigned int)0,
								PAIR_23);
  
  ResonancePdf* f0_1710  = new ResonancePdf("f0_1710", 
								fixAmps ? new Variable("f0_1710_amp_real", 0.070 * (-_mD02)) : 
								new Variable("f0_1710_amp_real", 0.070 * (-_mD02), 0.001, 0, 0),
								fixAmps ? new Variable("f0_1710_amp_imag", 0.087 * (-_mD02)) : 
								new Variable("f0_1710_amp_imag", 0.087 * (-_mD02), 0.1, 0, 0),  
								new Variable("f0_1710_mass",     1.714, 0.01, 1.5, 2.9), 
								new Variable("f0_1710_width",    0.140, 0.01, 0.01, 0.5),
								(unsigned int)0,
								PAIR_23);
  
  ResonancePdf* f2_1270  = new ResonancePdf("f2_1270", 
								fixAmps ? new Variable("f2_1270_amp_real", -1.027 * (-_mD02inv)) : 
								new Variable("f2_1270_amp_real", -1.027 * (-_mD02inv), 0.001, 0, 0),
								fixAmps ? new Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv)) : 
								new Variable("f2_1270_amp_imag", -0.162 * (-_mD02inv), 0.1, 0, 0),  
								new Variable("f2_1270_mass",     1.2754, 0.01, 1.0, 1.5),
								new Variable("f2_1270_width",    0.1851, 0.01, 0.01, 0.4),
								2,
								PAIR_23);
  
  ResonancePdf* f0_600  = new ResonancePdf("f0_600", 
							       fixAmps ? new Variable("f0_600_amp_real", 0.068 * (-_mD02)) : 
							       new Variable("f0_600_amp_real", 0.068 * (-_mD02), 0.001, 0, 0),
							       fixAmps ? new Variable("f0_600_amp_imag", 0.010 * (-_mD02)) : 
							       new Variable("f0_600_amp_imag", 0.010 * (-_mD02), 0.1, 0, 0),  
							       new Variable("f0_600_mass",     0.500, 0.01, 0.3, 0.7),
							       new Variable("f0_600_width",    0.400, 0.01, 0.2, 0.6), 
							       (unsigned int)0,
							       PAIR_23);
  
  ResonancePdf* nonr  = new ResonancePdf("nonr",
							     fixAmps ? new Variable("nonr_amp_real", 0.5595 * (-1)) : 
							     new Variable("nonr_amp_real", 0.5595 * (-1),   0.001, 0, 0),
							     fixAmps ? new Variable("nonr_amp_imag", -0.108761 * (-1)) : 
							     new Variable("nonr_amp_imag", -0.108761* (-1), 0.1, 0, 0)); 

  dtop0pp->resonances.push_back(nonr); 
  dtop0pp->resonances.push_back(rhop);
  dtop0pp->resonances.push_back(rho0); 
  dtop0pp->resonances.push_back(rhom); 
  dtop0pp->resonances.push_back(rhop_1450); 
  dtop0pp->resonances.push_back(rho0_1450); 
  dtop0pp->resonances.push_back(rhom_1450); 
  dtop0pp->resonances.push_back(rhop_1700); 
  dtop0pp->resonances.push_back(rho0_1700); 
  dtop0pp->resonances.push_back(rhom_1700); 
  dtop0pp->resonances.push_back(f0_980); 
  dtop0pp->resonances.push_back(f0_1370); 
  dtop0pp->resonances.push_back(f0_1500); 
  dtop0pp->resonances.push_back(f0_1710); 
  dtop0pp->resonances.push_back(f2_1270); 
  dtop0pp->resonances.push_back(f0_600); 

  if (!fitMasses) {
    for (vector<ResonancePdf*>::iterator res = dtop0pp->resonances.begin(); res != dtop0pp->resonances.end(); ++res) {
      (*res)->setParameterConstantness(true); 
    }
  }

  if (!eff) {
    // By default create a constant efficiency. 
    vector<Variable*> offsets;
    vector<Variable*> observables;
    vector<Variable*> coefficients; 

    observables.push_back(m12);
    observables.push_back(m13);
    offsets.push_back(constantZero);
    offsets.push_back(constantZero);
    coefficients.push_back(constantOne); 
    eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
  }
  comps.clear();
  comps.push_back(eff);
  if (!kzero_veto) makeKzeroVeto();
  comps.push_back(kzero_veto);
  ProdPdf* effWithVeto = new ProdPdf("effWithVeto", comps);

  return new DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, effWithVeto);
}

void drawFitPlotsWithPulls(TH1* hd, TH1* ht, TH1* hb, string plotdir){
    const char* hname = hd->GetName();
    char obsname[10];
    for (int i=0;;i++) {
        if (hname[i]=='_') obsname[i] = '\0';
        else obsname[i] = hname[i];
        if (obsname[i] == '\0') break;
    }
    ht->Scale(hd->Integral()/ht->Integral());
    hb->Scale(hd->Integral()/ht->Integral());
    foo->cd(); 
    foo->Clear();
    hd->Draw("ep");
    ht->Draw("lsame");
    hb->SetLineStyle(kDashed);
    hb->Draw("lsame");
    sprintf(strbuffer, "%s/%s_fit.C", plotdir.c_str(), obsname);
    foo->SaveAs(strbuffer);
    sprintf(strbuffer, "%s/%s_fit.pdf", plotdir.c_str(), obsname);
    foo->SaveAs(strbuffer);
/*    sprintf(strbuffer, "%s/%s_fit_log.pdf", plotdir.c_str(), obsname);
    foo->SaveAs(strbuffer);*/
}

void makeToyDalitzPdfPlots (GooPdf* overallSignal, string plotdir = "plots") {
  TH1F m12_dat_hist("m12_dat_hist", "", m12->numbins, m12->lowerlimit, m12->upperlimit);
  m12_dat_hist.SetStats(false); 
  m12_dat_hist.SetMarkerStyle(8); 
  m12_dat_hist.SetMarkerSize(1.2);
  m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
  sprintf(strbuffer, "Events / %.1f MeV", 1e3*m12_dat_hist.GetBinWidth(1));
  m12_dat_hist.GetYaxis()->SetTitle(strbuffer); 
  TH1F m12_pdf_hist("m12_pdf_hist", "", m12->numbins, m12->lowerlimit, m12->upperlimit);
  m12_pdf_hist.SetStats(false); 
  m12_pdf_hist.SetLineColor(kBlue); 
  m12_pdf_hist.SetLineWidth(3); 
  TH1* m12_pdf_hist_bkg = (TH1*)m12_pdf_hist.Clone("m12_pdf_hist_bkg");
  TH1F m13_dat_hist("m13_dat_hist", "", m13->numbins, m13->lowerlimit, m13->upperlimit);
  m13_dat_hist.SetStats(false); 
  m13_dat_hist.SetMarkerStyle(8); 
  m13_dat_hist.SetMarkerSize(1.2);
  m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
  sprintf(strbuffer, "Events / %.1f MeV", 1e3*m13_dat_hist.GetBinWidth(1));
  m13_dat_hist.GetYaxis()->SetTitle(strbuffer); 
  TH1F m13_pdf_hist("m13_pdf_hist", "", m13->numbins, m13->lowerlimit, m13->upperlimit);
  m13_pdf_hist.SetStats(false); 
  m13_pdf_hist.SetLineColor(kBlue); 
  m13_pdf_hist.SetLineWidth(3); 
  TH1* m13_pdf_hist_bkg = (TH1*)m13_pdf_hist.Clone("m13_pdf_hist_bkg");
  TH1F m23_dat_hist("m23_dat_hist", "", m13->numbins, m13->lowerlimit, m13->upperlimit);
  m23_dat_hist.SetStats(false); 
  m23_dat_hist.SetMarkerStyle(8); 
  m23_dat_hist.SetMarkerSize(1.2);
  m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{-}) [GeV]");
  sprintf(strbuffer, "Events / %.1f MeV", 1e3*m13_dat_hist.GetBinWidth(1));
  m23_dat_hist.GetYaxis()->SetTitle(strbuffer); 
  TH1F m23_pdf_hist("m23_pdf_hist", "", m13->numbins, m13->lowerlimit, m13->upperlimit);
  m23_pdf_hist.SetStats(false); 
  m23_pdf_hist.SetLineColor(kBlue); 
  m23_pdf_hist.SetLineWidth(3); 
  TH1* m23_pdf_hist_bkg = (TH1*)m23_pdf_hist.Clone("m23_pdf_hist_bkg");
  TH2F dalitzpp0_dat_hist("dalitzpp0_dat_hist", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit);
  dalitzpp0_dat_hist.SetStats(false); 
  dalitzpp0_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV]");
  dalitzpp0_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV]");
  TH2F dalitzpp0_pdf_hist("dalitzpp0_pdf_hist", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit);
/*  dalitzpp0_pdf_hist.GetXaxis()->SetTitle("m^{2}(K^{-} #pi^{0}) [GeV^{2}]");
  dalitzpp0_pdf_hist.GetYaxis()->SetTitle("m^{2}(K^{-} #pi^{+}) [GeV^{2}]");*/
  dalitzpp0_pdf_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{0}) [GeV^{2}]");
  dalitzpp0_pdf_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{0}) [GeV^{2}]");
  dalitzpp0_pdf_hist.SetStats(false); 
  double totalPdf = 0; 
  double totalDat = 0; 
  double totalSigProb = 0;
  double totalBGProb = 0;
  int evtCounter = 0; 
  double totalPdf_bg = 0;
  for (unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
    double data_m12 = data->getValue(m12, evt);
    m12_dat_hist.Fill(data_m12); 
    double data_m13 = data->getValue(m13, evt);
    m13_dat_hist.Fill(data_m13); 
    dalitzpp0_dat_hist.Fill(data_m12, data_m13);
    m23_dat_hist.Fill(cpuGetM23(data_m12, data_m13)); 
    totalSigProb += data->getValue(wSig0, evt);    
    totalBGProb += 1 - data->getValue(wSig0, evt);    
    totalDat++; 
  }
  wSig0->value = totalSigProb / totalDat;
    std::vector<Variable*> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber); 
    vars.push_back(wSig0);
    UnbinnedDataSet currData(vars); 

  for (int i = 0; i < m12->numbins; ++i) {
      m12->value = m12->lowerlimit + (m12->upperlimit - m12->lowerlimit)*(i + 0.5) / m12->numbins; 
      for (int j = 0; j < m13->numbins; ++j) {
          m13->value = m13->lowerlimit + (m13->upperlimit - m13->lowerlimit)*(j + 0.5) / m13->numbins; 
          if (!cpuDalitz(m12->value, m13->value, _mD0, piPlusMass, piPlusMass, piZeroMass)) continue;
          eventNumber->value = evtCounter; 
          evtCounter++;
          currData.addEvent(); 
      }
  }
  overallSignal->setData(&currData);
  signalDalitz->setDataSize(currData.getNumEvents(),4); 
  std::vector<std::vector<double> > pdfValues;
  overallSignal->getCompProbsAtDataPoints(pdfValues);
  for (unsigned int j = 0; j < pdfValues[0].size(); ++j) {
	double currm12 = currData.getValue(m12, j);
	double currm13 = currData.getValue(m13, j);

      dalitzpp0_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);
      m12_pdf_hist.Fill(currm12, pdfValues[0][j]);
      m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
      m23_pdf_hist.Fill(cpuGetM23(currm12, currm13), pdfValues[0][j]); 
      // NB: 1st for signal, and 2nd for bkg.
      m12_pdf_hist_bkg->Fill(currm12, pdfValues[2][j]);
      m13_pdf_hist_bkg->Fill(currm13, pdfValues[2][j]);
      m23_pdf_hist_bkg->Fill(cpuGetM23(currm12, currm13), pdfValues[2][j]); 
      totalPdf     += pdfValues[0][j]; 
      totalPdf_bg     += pdfValues[2][j]; 
  }
/*  m12_pdf_hist_bkg->Scale(totalBGProb/totalPdf_bg);
  m13_pdf_hist_bkg->Scale(totalBGProb/totalPdf_bg);
  m23_pdf_hist_bkg->Scale(totalBGProb/totalPdf_bg);*/
  foodal->cd(); 
  foodal->SetLogz(false);
  dalitzpp0_pdf_hist.Draw("colz");
  foodal->SaveAs((plotdir + "/dalitzpp0_pdf.png").c_str());
/*  m12_pdf_hist.Draw("");
  foodal->SaveAs((plotdir + "/m12_pdf_hist.png").c_str());
  m13_pdf_hist.Draw("");
  foodal->SaveAs((plotdir + "/m13_pdf_hist.png").c_str());
  if (!data) return;*/
  dalitzpp0_dat_hist.Draw("colz");
  foodal->SaveAs((plotdir + "/dalitzpp0_dat.png").c_str());

  drawFitPlotsWithPulls(&m12_dat_hist, &m12_pdf_hist, m12_pdf_hist_bkg, plotdir);
  drawFitPlotsWithPulls(&m13_dat_hist, &m13_pdf_hist, m13_pdf_hist_bkg, plotdir);
  drawFitPlotsWithPulls(&m23_dat_hist, &m23_pdf_hist, m23_pdf_hist_bkg, plotdir);
  delete m12_pdf_hist_bkg;
  delete m13_pdf_hist_bkg;
  delete m23_pdf_hist_bkg;
}

void runToyFit (std::string toyFileName) {
  m12 = new Variable("m12", 0, 3);
  m13 = new Variable("m13", 0, 3); 
  m12->numbins = 240;
  m13->numbins = 240;
  eventNumber = new Variable("eventNumber", 0, INT_MAX);
  wSig0 = new Variable("wSig0", 0, 1);
  getToyData(toyFileName);

  // EXERCISE 1 (real part): Create a PolynomialPdf which models
  // the efficiency you imposed in the preliminary, and use it in constructing
  // the signal PDF. 

  // EXERCISE 2: Create a K0 veto function and use it as the efficiency. 

  // EXERCISE 3: Make the efficiency a product of the two functions
  // from the previous exercises.

  int oldBins1 = m12->numbins;
  int oldBins2 = m13->numbins;
  m12->numbins = 60;
  m13->numbins = 60;  //Use different choice of binning for efficiency
  GooPdf* eff = makeEfficiencyPdf();
  m12->numbins = oldBins1;
  m13->numbins = oldBins2;   
  signalDalitz = makeSignalPdf(eff); 
  signalDalitz->setDataSize(data->getNumEvents(),4); 
  GooPdf* bkgFlatPdf = makeFlatBkgDalitzPdf();
  std::vector<Variable*> evtWeights;
  evtWeights.push_back(wSig0);
  comps.clear();
  comps.push_back(signalDalitz);
  comps.push_back(bkgFlatPdf);
  EventWeightedAddPdf* totPdf = new EventWeightedAddPdf("total", evtWeights, comps);
  totPdf->setData(data); 
  FitManager datapdf(totPdf); 
  
  gettimeofday(&startTime, NULL);
  startCPU = times(&startProc);
  datapdf.fit(); 
  stopCPU = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  makeToyDalitzPdfPlots(totPdf);   
}

int main (int argc, char** argv) {
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetFuncWidth(1);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(1);
  gStyle->SetPalette(1, 0);
  foo = new TCanvas(); 
  foodal = new TCanvas(); 
  foodal->Size(10, 10);

  cudaSetDevice(0);
  runToyFit(argv[1]);

  // Print total minimization time
  double myCPU = stopCPU - startCPU;
  double totalCPU = myCPU; 

  timersub(&stopTime, &startTime, &totalTime);
  std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
  std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl; 
  std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl; 
  myCPU = stopProc.tms_utime - startProc.tms_utime;
  std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;

  return 0; 
}
