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
#include "TTree.h" 

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

using namespace std;

TCanvas* foo; 
TCanvas* foodal; 
timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 
UnbinnedDataSet* data = 0; 

Variable* m12 = 0;
Variable* m13 = 0;
Variable* eventNumber = 0; 
bool fitMasses = false; 
Variable* fixedRhoMass  = new Variable("rho_mass", 0.77549, 0.01, 0.7, 0.8);
Variable* fixedRhoWidth = new Variable("rho_width", 0.1491, 0.01, 0.1, 0.2); 

Variable* veto_min = new Variable("veto_min", 0.475*0.475);
Variable* veto_max = new Variable("veto_min", 0.505*0.505);

const fptype _mDp = 1.86962;
const fptype piPlusMass = 0.13957018;
const fptype piZeroMass = 0.1349766;

// D+ -> pi- (1) pi+ (2) pi+ (3)

const fptype D1Mass = piPlusMass;//pi+
const fptype D2Mass = piPlusMass;//pi
const fptype D3Mass = piPlusMass;
const fptype D1Mass2 = D1Mass*D1Mass;
const fptype D2Mass2 = D2Mass*D2Mass;
const fptype D3Mass2 = D3Mass*D3Mass;
const fptype MMass = _mDp;
const fptype MMass2 = MMass*MMass;
const fptype MMass2inv = 1./MMass2; 

// Constants used in more than one PDF component. 
Variable* motherM = new Variable("motherM", MMass);
Variable* dau1M = new Variable("dau1M", D1Mass);
Variable* dau2M = new Variable("dau2M", D2Mass);
Variable* dau3M = new Variable("dau3M", D3Mass);
Variable* massSum = new Variable("massSum", MMass2 + D1Mass2+D2Mass2+D3Mass2); // = 3.53481 
Variable* constantOne = new Variable("constantOne", 1); 
Variable* constantZero = new Variable("constantZero", 0); 
  
std::vector<PdfBase*> comps;

GooPdf* kzero_veto = 0; 
char strbuffer[1000]; 
double mesonRad  = 1.5;
DalitzPlotPdf* signalDalitz; 

void makeToyDalitzData (GooPdf* overallSignal, const int iSeed = 0, string datadir = ".", const int nTotal = 1e5 ) ;
DalitzPlotPdf* makeSignalPdf (GooPdf* eff = 0) ;

fptype cpuGetM23 (fptype massPZ, fptype massPM) {
  return (massSum->value - massPZ - massPM); 
}

bool cpuDalitz (fptype m_12, fptype m_13, fptype bigM = MMass, fptype dm1 = D1Mass, fptype dm2 = D2Mass, fptype dm3 = D3Mass) {
//  if (m_12 > m_13) return false; // Only the upper corner is considered
  if (m_12 < POW(dm1 + dm2, 2)) return false; // This m_12 cannot exist, it's less than the square of the (1,2) particle mass.
  if (m_12 > POW(bigM - dm3, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter. 
  
  // Calculate energies of 1 and 3 particles in m_12 rest frame. 
  fptype e1star = 0.5 * (m_12 - dm2*dm2 + dm1*dm1) / SQRT(m_12); 
  fptype e3star = 0.5 * (bigM*bigM - m_12 - dm3*dm3) / SQRT(m_12); 

  // Bounds for m_13 at this value of m_12.
  fptype minimum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) + SQRT(e3star*e3star - dm3*dm3), 2);
  if (m_13 < minimum) return false;
  fptype maximum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) - SQRT(e3star*e3star - dm3*dm3), 2);
  if (m_13 > maximum) return false;

  return true; 
}

void makeToyDalitzData (GooPdf* overallSignal, const int iSeed, string datadir, const int nTotal ) {
    std::vector<Variable*> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber); 
    data = new UnbinnedDataSet(vars);
    UnbinnedDataSet currData(vars); 
  std::vector<std::vector<double> > pdfValues;
    int ncount = 0;
    TRandom3 donram(iSeed); 
  for (int i = 0; i < m12->numbins; ++i) {
      m12->value = m12->lowerlimit + (m12->upperlimit - m12->lowerlimit)*(i + 0.5) / m12->numbins; 
      for (int j = 0; j < m13->numbins; ++j) {
          m13->value = m13->lowerlimit + (m13->upperlimit - m13->lowerlimit)*(j + 0.5) / m13->numbins; 
          if (!cpuDalitz(m12->value, m13->value)) continue;
          eventNumber->value = ncount;
          ncount++;
          currData.addEvent(); 
      }
    }
    signalDalitz->setDataSize(currData.getNumEvents());
    overallSignal->setData(&currData);
//    std::vector<std::vector<double> > pdfValues;
    overallSignal->getCompProbsAtDataPoints(pdfValues);
    TH2F dalitzpp_dat_hist("dalitzpp_dat_hist", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit);
    dalitzpp_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-}#pi^{+}) [GeV^{2}]");
    dalitzpp_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV^{2}]");
    ncount = 0;
    ofstream writer;
    sprintf(strbuffer, "%s/dalitz_D3pi_toyMC_%03d.txt", datadir.c_str(), iSeed);
    writer.open(strbuffer);
    vector<double> fIntegral;
    fIntegral.push_back(pdfValues[0][0]);
    Int_t ncells = pdfValues[0].size();
    for (unsigned int j = 1; j < ncells; ++j) {
        fIntegral.push_back(pdfValues[0][j]+fIntegral[j-1]);
    }
    for (unsigned int j = 0; j < ncells; ++j)  fIntegral[j] /= fIntegral[ncells-1];
    ncount = 0;
    int nEvents = donram.Poisson(nTotal);
    for (int iEvt = 0;iEvt<nEvents;iEvt++){
        double r = donram.Rndm();
        //Binary search for fIntegral[cell-1] < r < fIntegral[cell]
        int lo = 0, hi = ncells-1, mid;
        while(lo <= hi){
            mid = lo + (hi-lo)/2;
            if( r<=fIntegral[mid]&&(mid==0||r>fIntegral[mid-1])) break;
            else if (r > fIntegral[mid] ) lo = mid+1;
            else hi = mid-1;
        }
        int j = mid;
        double currm12 = currData.getValue(m12, j);
        currm12 += (m12->upperlimit - m12->lowerlimit)*(donram.Rndm() - 0.5) / m12->numbins;
        double currm13 = currData.getValue(m13, j);
        currm13 += (m13->upperlimit - m13->lowerlimit)*(donram.Rndm() - 0.5) / m13->numbins;
        eventNumber->value = ncount++;
        dalitzpp_dat_hist.Fill(currm12, currm13);
        data->addEvent();
        writer << ncount-1 << '\t'<<currm12 << '\t'<<currm13<<std::endl;
    }
    writer.close(); 
    std::cout<<"Entries generated: "<<data->numEvents()<<std::endl;
    foodal->cd(); 
    foodal->SetLogz(false);
    dalitzpp_dat_hist.Draw("colz");
    foodal->SaveAs("dalitzpp_dat_temp.png");
}

void runToyGeneration(int numFile = 0){
  m12   = new Variable("m12",   0.0, 3.0);
  m12->numbins = 1500;
//  m12   = new Variable("m12",   0.4, 3.0);
  m13   = new Variable("m13",   0.0, 3.0); 
  m13->numbins = 1500;
  eventNumber = new Variable("eventNumber", 0, INT_MAX);
  signalDalitz = makeSignalPdf(); 
  vector<PdfBase*> comps;
  comps.clear(); 
  comps.push_back(signalDalitz);
//  comps.push_back(sig0_jsugg); 
  std::cout << "Creating overall PDF\n"; 
  ProdPdf* overallSignal = new ProdPdf("overallSignal", comps);
  gettimeofday(&startTime, NULL);
  startCPU = times(&startProc);
//  makeToyDalitzData (signalDalitz);
  makeToyDalitzData (overallSignal, numFile);
  stopCPU = times(&stopProc);
  gettimeofday(&stopTime, NULL);
}

void getToyData (std::string toyFileName) {
  TH2F dalitzplot("dalitzplot", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit); 
  std::vector<Variable*> vars;
  vars.push_back(m12);
  vars.push_back(m13);
  vars.push_back(eventNumber); 
  data = new UnbinnedDataSet(vars); 
//  const int MAXEVT = 1e4;

  const string suffix = ".root";
  if (toyFileName.rfind(suffix)+suffix.length() == toyFileName.length()){
      std::cout<<"Reading file "<<toyFileName<<std::endl;
      TFile*f = TFile::Open(toyFileName.c_str());
      TTree*t = (TTree*)f->Get("DecayTree");
      std::cout<<"Entries: "<<t->GetEntries()<<std::endl;
      assert(t);
      double m2_12, m2_13;
      t->SetBranchAddress("s12", &m2_12);
      t->SetBranchAddress("s13", &m2_13);
      for (int i=0;i<t->GetEntries()/*&&i<MAXEVT*/;i++){
          t->GetEntry(i);
          if (i%10) continue;// Only accept 10% of data
          m12->value = m2_12;
          m13->value = m2_13;
          eventNumber->value = data->getNumEvents(); 
          data->addEvent(); 
          dalitzplot.Fill(m12->value, m13->value); 
      }
      f->Close();
  }
  else{
  std::ifstream reader;
  reader.open(toyFileName.c_str()); 
  std::string buffer;
  while (!reader.eof()) {
    reader >> buffer;
    if (buffer == "====") break; 
    std::cout << buffer; 
  }

  double dummy = 0; 
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

    eventNumber->value = data->getNumEvents(); 
    data->addEvent(); 

    dalitzplot.Fill(m12->value, m13->value); 
  }}

  dalitzplot.SetStats(false); 
  dalitzplot.Draw("colz");
  foodal->SaveAs("dalitzplot.png"); 
}

GooPdf* makeKzeroVeto () {
  if (kzero_veto) return kzero_veto; 

  VetoInfo* kVetoInfo = new VetoInfo();
  kVetoInfo->cyclic_index = PAIR_12; 
  kVetoInfo->minimum = veto_min;
  kVetoInfo->maximum = veto_max;
  VetoInfo* kVetoInfo2 = new VetoInfo();
  kVetoInfo2->cyclic_index = PAIR_13; 
  kVetoInfo2->minimum = veto_min;
  kVetoInfo2->maximum = veto_max;
  vector<VetoInfo*> vetos; 
  vetos.push_back(kVetoInfo); 
  vetos.push_back(kVetoInfo2); 
  kzero_veto = new DalitzVetoPdf("kzero_veto", m12, m13, motherM, dau1M, dau2M, dau3M, vetos); 
  return kzero_veto;
}

DalitzPlotPdf* makeSignalPdf (GooPdf* eff ) {
  DecayInfo* dtop0pp = new DecayInfo();
  dtop0pp->motherMass  = MMass; 
  dtop0pp->daug1Mass  = D1Mass;
  dtop0pp->daug2Mass  = D2Mass;
  dtop0pp->daug3Mass  = D3Mass;
  dtop0pp->meson_radius  = 1.5; 
 
  Variable *rho0_amp_real = new Variable("rho0_amp_real", 1);
  Variable *rho0_amp_imag = new Variable("rho0_amp_imag", 0);  
  ResonancePdf* rho0  = new ResonancePdf("rho0",
							     rho0_amp_real,
							     rho0_amp_imag,          
							     fixedRhoMass,
							     fixedRhoWidth,
							     1,
							     PAIR_12, RES_RBW, true);


  bool fixAmps = false;

  Variable *f0_980_amp_real = new Variable("f0_980_amp_real", -.29999,  0.001, 0, 0);
  Variable *f0_980_amp_imag = new Variable("f0_980_amp_imag", -0.2, 0.001, 0, 0);  
  Variable* f0980Mass  = new Variable("f0_980_mass", 0.990, 0.01, 0.8, 1.2);
  f0980Mass->fixed = true;
  Variable* f0980Width = new Variable("f0_980_width", 0.070, 0.001, 0.001, 0.28); 
  f0980Width->value = 0.165; f0980Width->fixed = true;
  Variable* f0_980_rg2og1 = new Variable("f0_980_rg2og1",    4.21, 0.1, 1e-3, 10);
  f0_980_rg2og1->fixed = true;
  f0_980_amp_real->fixed = f0_980_amp_imag->fixed = fixAmps;
  ResonancePdf* f0_980  = new ResonancePdf("f0_980", 
                                   f0_980_amp_real,
                                   f0_980_amp_imag,
                                   f0980Mass, 
                                   f0980Width,                                  
							       //(unsigned int)0,
                                   f0_980_rg2og1,
							       PAIR_12, 
                                   //RES_RBW, 
                                   RES_FLATTE,
                                   true);

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
								  PAIR_12, RES_RBW, true);


  ResonancePdf* f2_1270  = new ResonancePdf("f2_1270", 
								fixAmps ? new Variable("f2_1270_amp_real", -1.027 * (-MMass2inv)) : 
								new Variable("f2_1270_amp_real", -1.027 * (-MMass2inv), 0.001, 0, 0),
								fixAmps ? new Variable("f2_1270_amp_imag", -0.162 * (-MMass2inv)) : 
								new Variable("f2_1270_amp_imag", -0.162 * (-MMass2inv), 0.1, 0, 0),  
								new Variable("f2_1270_mass",     1.2751, 0.01, 1.0, 1.5),
								new Variable("f2_1270_width",    0.1851, 0.01, 0.01, 0.4),
								2,
								PAIR_12, RES_RBW, true);
  
  ResonancePdf* f0_500  = new ResonancePdf("f0_500", 
							       fixAmps ? new Variable("f0_500_amp_real", 0.068 * (-MMass2)) : 
							       new Variable("f0_500_amp_real", 0.068 * (-MMass2), 0.001, 0, 0),
							       fixAmps ? new Variable("f0_500_amp_imag", 0.010 * (-MMass2)) : 
							       new Variable("f0_500_amp_imag", 0.010 * (-MMass2), 0.1, 0, 0),  
							       new Variable("f0_500_mass",     0.500, 0.01, 0.3, 0.7),
							       new Variable("f0_500_width",    0.400, 0.01, 0.2, 0.6), 
							       (unsigned int)0,
							       PAIR_12, RES_RBW, true);
  
  ResonancePdf* nonr  = new ResonancePdf("nonr",
							     fixAmps ? new Variable("nonr_amp_real", 0.5595 * (-1)) : 
							     new Variable("nonr_amp_real", 0.5595 * (-1),   0.001, 0, 0),
							     fixAmps ? new Variable("nonr_amp_imag", -0.108761 * (-1)) : 
							     new Variable("nonr_amp_imag", -0.108761* (-1), 0.1, 0, 0)); 


  dtop0pp->resonances.push_back(rho0);
  dtop0pp->resonances.push_back(f0_980);
  dtop0pp->resonances.push_back(f2_1270);
  dtop0pp->resonances.push_back(rhop_1450);
  dtop0pp->resonances.push_back(f0_500);
  dtop0pp->resonances.push_back(nonr); 

  if (!fitMasses) {
    for (vector<ResonancePdf*>::iterator res = dtop0pp->resonances.begin(); res != dtop0pp->resonances.end(); ++res) {
      const std::string resname = (*res)->getName();
//      if (resname != "f0_500")
      if (resname != "f0_980")
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
//  comps.push_back(kzero_veto);
  ProdPdf* effWithVeto = new ProdPdf("effWithVeto", comps);

  return new DalitzPlotPdf("signalPDF", m12, m13, eventNumber, dtop0pp, effWithVeto);
}

void drawFitPlotsWithPulls(TH1* hd, TH1* ht, string plotdir){
    const char* hname = hd->GetName();
    char obsname[10];
    for (int i=0;;i++) {
        if (hname[i]=='_') obsname[i] = '\0';
        else obsname[i] = hname[i];
        if (obsname[i] == '\0') break;
    }
    ht->Scale(hd->Integral()/ht->Integral());
    foo->cd(); 
    foo->Clear();
    ht->Draw("c");
    hd->Draw("epsame");
    sprintf(strbuffer, "%s/%s_fit.png", plotdir.c_str(), obsname);
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
  m12_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
  sprintf(strbuffer, "Events / %.1f MeV", 1e3*m12_dat_hist.GetBinWidth(1));
  m12_dat_hist.GetYaxis()->SetTitle(strbuffer); 
  TH1F m12_pdf_hist("m12_pdf_hist", "", m12->numbins, m12->lowerlimit, m12->upperlimit);
  m12_pdf_hist.SetStats(false); 
  m12_pdf_hist.SetLineColor(kBlue); 
  m12_pdf_hist.SetLineWidth(3); 
  TH1F m13_dat_hist("m13_dat_hist", "", m13->numbins, m13->lowerlimit, m13->upperlimit);
  m13_dat_hist.SetStats(false); 
  m13_dat_hist.SetMarkerStyle(8); 
  m13_dat_hist.SetMarkerSize(1.2);
  m13_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
  sprintf(strbuffer, "Events / %.1f MeV", 1e3*m13_dat_hist.GetBinWidth(1));
  m13_dat_hist.GetYaxis()->SetTitle(strbuffer); 
  TH1F m13_pdf_hist("m13_pdf_hist", "", m13->numbins, m13->lowerlimit, m13->upperlimit);
  m13_pdf_hist.SetStats(false); 
  m13_pdf_hist.SetLineColor(kBlue); 
  m13_pdf_hist.SetLineWidth(3); 
  TH1F m23_dat_hist("m23_dat_hist", "", m13->numbins, m13->lowerlimit, m13->upperlimit);
  m23_dat_hist.SetStats(false); 
  m23_dat_hist.SetMarkerStyle(8); 
  m23_dat_hist.SetMarkerSize(1.2);
  m23_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{+} #pi^{+}) [GeV]");
  sprintf(strbuffer, "Events / %.1f MeV", 1e3*m13_dat_hist.GetBinWidth(1));
  m23_dat_hist.GetYaxis()->SetTitle(strbuffer); 
  TH1F m23_pdf_hist("m23_pdf_hist", "", m13->numbins, m13->lowerlimit, m13->upperlimit);
  m23_pdf_hist.SetStats(false); 
  m23_pdf_hist.SetLineColor(kBlue); 
  m23_pdf_hist.SetLineWidth(3); 
  double totalPdf = 0; 
  double totalDat = 0; 
  TH2F dalitzpp_dat_hist("dalitzpp_dat_hist", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit);
  dalitzpp_dat_hist.SetStats(false); 
  dalitzpp_dat_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
  dalitzpp_dat_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV]");
  TH2F dalitzpp_pdf_hist("dalitzpp_pdf_hist", "", m12->numbins, m12->lowerlimit, m12->upperlimit, m13->numbins, m13->lowerlimit, m13->upperlimit);
/*  dalitzpp_pdf_hist.GetXaxis()->SetTitle("m^{2}(K^{-} #pi^{0}) [GeV^{2}]");
  dalitzpp_pdf_hist.GetYaxis()->SetTitle("m^{2}(K^{-} #pi^{+}) [GeV^{2}]");*/
  dalitzpp_pdf_hist.GetXaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV^{2}]");
  dalitzpp_pdf_hist.GetYaxis()->SetTitle("m^{2}(#pi^{-} #pi^{+}) [GeV^{2}]");
  dalitzpp_pdf_hist.SetStats(false); 
    std::vector<Variable*> vars;
    vars.push_back(m12);
    vars.push_back(m13);
    vars.push_back(eventNumber); 
    UnbinnedDataSet currData(vars); 
  int evtCounter = 0; 

  for (int i = 0; i < m12->numbins; ++i) {
      m12->value = m12->lowerlimit + (m12->upperlimit - m12->lowerlimit)*(i + 0.5) / m12->numbins; 
      for (int j = 0; j < m13->numbins; ++j) {
          m13->value = m13->lowerlimit + (m13->upperlimit - m13->lowerlimit)*(j + 0.5) / m13->numbins; 
          if (!cpuDalitz(m12->value, m13->value)) continue;
          eventNumber->value = evtCounter; 
          evtCounter++;
          currData.addEvent(); 
      }
  }
  overallSignal->setData(&currData);
  signalDalitz->setDataSize(currData.getNumEvents()); 
  std::vector<std::vector<double> > pdfValues;
  overallSignal->getCompProbsAtDataPoints(pdfValues);
  for (unsigned int j = 0; j < pdfValues[0].size(); ++j) {
	double currm12 = currData.getValue(m12, j);
	double currm13 = currData.getValue(m13, j);

      dalitzpp_pdf_hist.Fill(currm12, currm13, pdfValues[0][j]);
      m12_pdf_hist.Fill(currm12, pdfValues[0][j]);
      m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
      m23_pdf_hist.Fill(cpuGetM23(currm12, currm13), pdfValues[0][j]); 
      totalPdf     += pdfValues[0][j]; 
  }
  foodal->cd(); 
  foodal->SetLogz(false);
  dalitzpp_pdf_hist.Draw("colz");
  foodal->SaveAs((plotdir + "/dalitzpp_pdf.png").c_str());
/*  m12_pdf_hist.Draw("");
  foodal->SaveAs((plotdir + "/m12_pdf_hist.png").c_str());
  m13_pdf_hist.Draw("");
  foodal->SaveAs((plotdir + "/m13_pdf_hist.png").c_str());
  if (!data) return;*/
  for (unsigned int evt = 0; evt < data->getNumEvents(); ++evt) {
    double data_m12 = data->getValue(m12, evt);
    m12_dat_hist.Fill(data_m12); 
    double data_m13 = data->getValue(m13, evt);
    m13_dat_hist.Fill(data_m13); 
    dalitzpp_dat_hist.Fill(data_m12, data_m13);
    m23_dat_hist.Fill(cpuGetM23(data_m12, data_m13)); 
    totalDat++; 
  }
  dalitzpp_dat_hist.Draw("colz");
  foodal->SaveAs((plotdir + "/dalitzpp_dat.png").c_str());

  drawFitPlotsWithPulls(&m12_dat_hist, &m12_pdf_hist, plotdir);
  drawFitPlotsWithPulls(&m13_dat_hist, &m13_pdf_hist, plotdir);
  drawFitPlotsWithPulls(&m23_dat_hist, &m23_pdf_hist, plotdir);
}

void runToyFit (std::string toyFileName) {
  m12 = new Variable("m12", 0, 3);
  m13 = new Variable("m13", 0, 3); 
  m12->numbins = 300;
  m13->numbins = 300;
  eventNumber = new Variable("eventNumber", 0, INT_MAX);
  getToyData(toyFileName);

  // EXERCISE 1 (real part): Create a PolynomialPdf which models
  // the efficiency you imposed in the preliminary, and use it in constructing
  // the signal PDF. 

  // EXERCISE 2: Create a K0 veto function and use it as the efficiency. 

  // EXERCISE 3: Make the efficiency a product of the two functions
  // from the previous exercises.

  signalDalitz = makeSignalPdf(); 
  comps.clear();
  comps.push_back(signalDalitz);
  ProdPdf* overallSignal = new ProdPdf("overallSignal", comps);
  overallSignal->setData(data); 
  signalDalitz->setDataSize(data->getNumEvents()); 
  FitManager datapdf(overallSignal); 
  
  gettimeofday(&startTime, NULL);
  startCPU = times(&startProc);
  datapdf.fit(); 
  stopCPU = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  datapdf.getMinuitValues();

  //Get the fractions w/ uncertainties
  vector<double> fracList;
  signalDalitz->getFractions(fracList);
/*  const int nRes = fracList.size();
  vector <float> fractions[nRes];
  float mean[nRes];
  float rms[nRes];
  for (int ii=0;ii<nRes;ii++) mean[ii] = rms[ii] = 0;
  for (int ii=0;ii<nSamples;ii++){
      datapdf.loadSample(ii);
      signalDalitz->getFractions(fracList);
      for (int jj=0;jj<nRes; jj++) {
          fractions[jj].push_back(fracList[jj]);
          mean[jj] += fracList[jj];
          rms[jj] += fracList[jj]*fracList[jj];
      }
  }
  TH1F* hFracs[nRes];
  TFile * froot = new TFile("fractionHists.root", "recreate");
  for (int ii=0;ii<nRes;ii++) {
      mean[ii] /= nSamples;
      rms[ii] = sqrt(rms[ii]/nSamples-mean[ii]*mean[ii]);
      sprintf(strbuffer, "hfrac_res%d", ii);
      hFracs[ii] = new TH1F(strbuffer, "", 100, mean[ii]-4*rms[ii], mean[ii]+4*rms[ii]);
      for (int jj=0;jj<nSamples;jj++)
          hFracs[ii]->Fill(fractions[ii][jj]);
      hFracs[ii]->Write();
  }
  froot->Close();*/

  makeToyDalitzPdfPlots(overallSignal);   
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
  if (argc ==2) runToyFit(argv[1]);
  else{
  int fitToRun = atoi(argv[1]);
  int ifile = -1;
  switch (fitToRun) {
  case 0:
      ifile = atoi(argv[2]);
      sprintf(strbuffer, "dalitz_D3pi_toyMC_%03d.txt", ifile);
      runToyFit(strbuffer);  break;
  case 1: runToyGeneration(atoi(argv[2]));  break; 
  default: break; 
  }
  }

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
