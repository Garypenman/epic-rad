//need to use reconstructed variables for all particles
//rather than set index by number, need to use lambda method like is done for positron.

#include "ePICReaction.h"
#include "ParticleCreator.h"
#include "ePICParticleCreator.h"
#include "Indicing.h"
#include "Histogrammer.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TMath.h>

void ProcessMCMatchedY(){
  using namespace rad::names::data_type; //for Rec(), Truth()
  
  gBenchmark->Start("df total");
  
  rad::config::ePICReaction epic{"events","/w/work5/home/garyp/eic/Farm/Y4260/recon/*.root"};
  epic.SetBeamsFromMC(); //for this file 0=ebeam 1=pbeam
  
  //epic.AliasColumns();
  //epic.AliasMC();
  //epic.AliasColumnsAndMC();
  epic.AliasColumnsAndMatchWithMC();
  
  //these work but only if MCScatteredElectrons and MCScatteredProtons only have 1 particle in them??
  //apparently this isnt always the case. Benching them for now.
  epic.setScatElectron(rad::indice::UseAsID(0), {"MCScatteredElectrons_objIdx.index"});
  epic.setParticleIndex("pprime",rad::indice::UseAsID(0),{"MCScatteredProtons_objIdx.index"},2212);

  
  //particle creator
  rad::epic::ePICParticleCreator epic_particles{epic};
  
  //SUBTRACT 2 FROM WHATEVER THE INDEXING IS IN THE 
  //MCPARTICLES/HEPMC3 LIST, SINCE BEAM PARTICLES ARE
  //REMOVED FROM THE LIST
  //scattered electron
  //epic.setScatElectronIndex(0); 
  //epic_particles.LowQ2Electron();
  epic_particles.MCMatchedLowQ2Electron();
  
  //recoil proton (the baryon)
  //epic.setParticleIndex("pprime",1);
  //epic_particles.RomanPotProton();
  epic_particles.MCMatchedRomanPotProton();
  epic_particles.MCMatchedB0Proton();
  
  epic.setParticleIndex("pim",2,-211);
  epic.setParticleIndex("pip",3,211);
  
  epic.setParticleIndex("ele",4,11);
  epic.setParticleIndex("pos",5,-11);
  
  epic.Particles().Sum("Jpsi",{"ele","pos"});
  epic.Particles().Sum("Y",{"pim","pip","Jpsi"});
  
  epic.setMesonParticles({"Jpsi","pip","pim"});
  
  
  epic.Particles().Miss("calc_pprime",{rad::names::ScatEle().data(),"Y"});
  
  //set this if not detecting proton
  // i.e. semi inclusive measurement
  //epic.setBaryonParticles({"calc_pprime"});
  
  //set this if we are goign to detect the proton directly
  //for now RP and B0 dont seem to have the acceptance
  //at 18x275 atleast!
  epic.setBaryonParticles({"pprime"});
  
  
  //must call this after all particles are configured
  epic.makeParticleMap();
  
  
  //option filtering of reconstructed tracks
  //epic.Filter("el_OK==1&&po_OK==1","partFilter");
  epic.Filter("rec_pmag[scat_ele]>0.1","pmag_scat_ele_filt");
  //epic.Filter("rec_pmag[pprime]>0.1","pmag_pprime_filt");
  
  //////////////////////////////////////////////////////////
  // Now define calculated variables
  // Note reconstructed variables will have rec_ prepended
  // truth variables will have tru_ prepended
  //////////////////////////////////////////////////////////

  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"Welec","{scat_ele}");
  rad::rdf::Mass(epic,"Whad","{Y,pprime}");
  rad::rdf::Mass(epic,"JMass","{Jpsi}");
  rad::rdf::Mass(epic,"YMass","{Y}");
  rad::rdf::Mass(epic,"MissNMass","{calc_pprime}");
  rad::rdf::Q2(epic,"Q2");

  //t distribution, column name
  rad::rdf::TTop(epic,"t_top");
  rad::rdf::TBot(epic,"t_bot");
  rad::rdf::TPrimeTop(epic,"tp_top");
  rad::rdf::TPrimeBot(epic,"tp_bot");
  
  //CM production angles
  rad::rdf::CMAngles(epic,"CM");

  //exlusivity
  rad::rdf::MissMass(epic,"MissMass","{scat_ele,pprime,Y}");
  rad::rdf::MissP(epic,"MissP_Meson","{scat_ele,Y}");
  rad::rdf::MissPt(epic,"MissPt_Meson","{scat_ele,Y}");
  rad::rdf::MissPz(epic,"MissPz_Meson","{scat_ele,Y}");
  rad::rdf::MissTheta(epic,"MissTheta_Meson","{scat_ele,Y}");

  ///////////////////////////////////////////////////////////
  //Define histograms
  ///////////////////////////////////////////////////////////
  rad::histo::Histogrammer histo{"set1",epic};
  // we can create many histograms by splitting events into
  // bins, where the same histogram is produced for the given bin
  // e.g. create 10 bins in tru_W between 4 and 54 GeV 
  //histo.Splitter().AddRegularDimension(Truth()+"W", rad::histo::RegularSplits(10,4,14) );
  //can add as many split dimensions as we like
  //histo.Splitter().AddRegularDimension("xxx", rad::histo::RegularSplits(nbins,low,high) );
  histo.Init({Rec(),Truth()});//will create histograms for rec and truth

  histo.Create<TH1D,double>({"hQ2",";Q^{2} [GeV^{2}]",100,0,2.},{"Q2"});
  histo.Create<TH1D,double>({"hWelec",";W (electro miss mass) [GeV/c^{2}]",100,0,200.},{"Welec"});
  histo.Create<TH1D,double>({"hWhad",";W (hadro final state) [GeV/c^{2}]",100,0,200.},{"Whad"});
  histo.Create<TH1D,double>({"hJMass",";M(e-,e+) [GeV/c^{2}]",100,2.,5.},{"JMass"});
  histo.Create<TH1D,double>({"hYMass",";M(e-,e+, #pi^{+},#pi^{-}) [GeV/c^{2}]",100,2.,5.},{"YMass"});
  histo.Create<TH1D,double>({"hMissMass",";M_{miss} [GeV/c^{2}]",1000,-10,10},{"MissMass"});
  
  histo.Create<TH1D,double>({"httop",";t(top vertex) [GeV^{2}]",100,-1,5},{"t_top"});
  histo.Create<TH1D,double>({"htbot",";t(bottom vertex) [GeV^{2}]",100,-1,5},{"t_bot"});
  histo.Create<TH1D,double>({"htptop",";t'(top vertex) [GeV^{2}]",100,-1,5},{"tp_top"});
  histo.Create<TH1D,double>({"htpbot",";t'(bottom vertex) [GeV^{2}]",100,-1,5},{"tp_bot"});
  
  histo.Create<TH1D,double>({"hcthCM",";cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"hphCM",";#phi_{CM}",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  
  histo.Create<TH1D,double>({"hmissP",";p_{miss}(e',Y)",1000,-100,100},{"MissP_Meson"});
  histo.Create<TH1D,double>({"hmissPt",";p_{t,miss}(e',Y)",100,0,10},{"MissPt_Meson"});
  histo.Create<TH1D,double>({"hmissPz",";p_{z,miss}(e',Y)",1000,-100,100},{"MissPz_Meson"});
  histo.Create<TH1D,double>({"hmissTheta",";#theta_{miss}(e',Y)",100,-TMath::Pi(),TMath::Pi()},{"MissTheta_Meson"});
 
  //particle momenta
  histo.Create<TH1D,double>({"hpmag_elec",";p_{e'} [GeV/c]",100,-1,20},{"pmag[scat_ele]"});
  histo.Create<TH1D,double>({"heta_elec",";#eta_{e'} ",100,-10,10},{"eta[scat_ele]"});
  histo.Create<TH1D,double>({"htheta_elec",";#theta_{e'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"theta[scat_ele]"});
  histo.Create<TH1D,double>({"hphi_elec",";#phi_{e'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"phi[scat_ele]"});
  
  histo.Create<TH1D,double>({"hpmag_pprime",";p_{p'} [GeV/c]",100,-1,300},{"pmag[pprime]"});
  histo.Create<TH1D,double>({"heta_pprime",";#eta_{p'} ",100,-10,10},{"eta[pprime]"});
  histo.Create<TH1D,double>({"htheta_pprime",";#theta_{p'} [rad]",100,0,0.1},{"theta[pprime]"});
  histo.Create<TH1D,double>({"hphi_pprime",";#phi_{p'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"phi[pprime]"});
  
  //finally, book lazy snapshot before processing
  //benchmark here will be zero if lazy snapshot
  //which now works and will be booked till trigger
  gBenchmark->Start("snapshot");
  epic.BookLazySnapshot("MCMatchedY.root");
  //epic.ImmediateSnapshot("MCMatchedY.root");
  gBenchmark->Stop("snapshot");
  gBenchmark->Print("snapshot");

  gBenchmark->Start("processing");
  ///////////////////////////////////////////////////////////
  //Draw histograms
  ///////////////////////////////////////////////////////////
  
  
  TCanvas *c00 = new TCanvas();
  c00->Divide(3,2);
  c00->cd(1)->SetLogy();
  histo.DrawSame("hQ2",gPad);
  c00->cd(2)->SetLogy();
  histo.DrawSame("hWelec",gPad);
  c00->cd(3)->SetLogy();
  histo.DrawSame("hWhad",gPad);
  c00->cd(4)->SetLogy();
  histo.DrawSame("hJMass",gPad);
  c00->cd(5)->SetLogy();
  histo.DrawSame("hYMass",gPad);
  c00->cd(6)->SetLogy();
  histo.DrawSame("hMissMass",gPad);

  TCanvas *c01 = new TCanvas();
  c01->Divide(2,2);
  c01->cd(1)->SetLogy();
  histo.DrawSame("httop",gPad);
  c01->cd(2)->SetLogy();
  histo.DrawSame("htbot",gPad);
  c01->cd(3)->SetLogy();
  histo.DrawSame("htptop",gPad);
  c01->cd(4)->SetLogy();
  histo.DrawSame("htpbot",gPad);
  
  TCanvas *c02 = new TCanvas();
  c02->Divide(2,2);
  c02->cd(1)->SetLogy();
  histo.DrawSame("hmissP",gPad);
  c02->cd(2)->SetLogy();
  histo.DrawSame("hmissPt",gPad);
  c02->cd(3)->SetLogy();
  histo.DrawSame("hmissPz",gPad);
  c02->cd(4)->SetLogy();
  histo.DrawSame("hmissTheta",gPad);
  
  //reco elec kin
  TCanvas *c03 = new TCanvas();
  c03->Divide(2,2);
  c03->cd(1)->SetLogy();
  histo.DrawSame("hpmag_elec",gPad);
  c03->cd(2)->SetLogy();
  histo.DrawSame("heta_elec",gPad);
  c03->cd(3)->SetLogy();
  histo.DrawSame("htheta_elec",gPad);
  c03->cd(4)->SetLogy();
  histo.DrawSame("hphi_elec",gPad);
  
  //reco proton kin
  TCanvas *c04 = new TCanvas();
  c04->Divide(2,2);
  c04->cd(1)->SetLogy();
  histo.DrawSame("hpmag_pprime",gPad);
  c04->cd(2)->SetLogy();
  histo.DrawSame("heta_pprime",gPad);
  c04->cd(3)->SetLogy();
  histo.DrawSame("htheta_pprime",gPad);
  c04->cd(4)->SetLogy();
  histo.DrawSame("hphi_pprime",gPad);
  
  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");

  //save all histograms to file
  histo.File("MCMatchedY_hists.root");

  gBenchmark->Stop("df total");
  gBenchmark->Print("df total");
  
  
}
