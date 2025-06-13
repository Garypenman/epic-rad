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
#include "gammaN_2_Spin0Spin0SpinHalfRDF.h"
#include <TBenchmark.h>
#include <TCanvas.h>
#include <TMath.h>


void ProcessMCMatched_DVCS(){
  using namespace rad::names::data_type; //for Rec(), Truth()

  gBenchmark->Start("df");

   rad::config::ePICReaction epic{"events","~/EIC/data/DVCS/DVCS_10x100.root"};
   epic.SetBeamsFromMC(0,3); //for this file 0=ebeam 1=pbeam
   epic.AliasColumnsAndMatchWithMC();
  
  //Assign particles names and indices
  //indicing comes from ordering in hepmc file
   //epic.setBeamIonIndex(3);
   // epic.setBeamElectronIndex(0);
   epic.setScatElectronIndex(0);
  //give final state hadrons names,
  //if we give a PDG code it will generate el_OK branches etc
  //el_OK = 1 if electron reconstructed with right PDG
  epic.setParticleIndex("photon",1);
  epic.setParticleIndex("proton",2);

  //particle creator
  rad::epic::ePICParticleCreator epic_particles{epic};
  //hack in proton from far forward particles
  epic_particles.MCMatchedFarForwardProton("proton");

  //Group particles into top and bottom vertices
  //aka Meson and Baryon components
  //this is required for calculating reaction kinematics
  //e.g. t distributions
  epic.setBaryonParticles({"proton"});
  epic.setMesonParticles({"photon"});

  //must call this after all particles are configured
  epic.makeParticleMap();
  
  //masses column name, {+ve particles}, {-ve particles}
  rad::rdf::MissMass(epic,"W","{scat_ele}");
  rad::rdf::Q2(epic,"Q2");
  rad::rdf::Mass(epic,"Whad","{photon,proton}");

  //t distribution, column name
  rad::rdf::TBot(epic,"tb_pn");
  rad::rdf::TPrimeBot(epic,"tbp_pn");
  rad::rdf::TTop(epic,"tt_gZ");
  rad::rdf::TPrimeTop(epic,"ttp_gZ");

  //CM production angles
  rad::rdf::CMAngles(epic,"CM");
  rad::rdf::CosThetaProtonRest(epic,"CosTheta_PRest");
  rad::rdf::PhiProtonRest(epic,"Phi_PRest");

  //exlusivity
  rad::rdf::MissMass(epic,"MissMass","{scat_ele,proton,photon}");
  rad::rdf::MissP(epic,"MissP","{scat_ele,proton,photon}");
  rad::rdf::MissE(epic,"MissE","{scat_ele,proton,photon}");
  rad::rdf::MissPt(epic,"MissPt","{scat_ele,proton,photon}");
  rad::rdf::MissPz(epic,"MissPz","{scat_ele,proton,photon}");
  rad::rdf::MissTheta(epic,"MissTheta","{scat_ele,proton,photon}");

  rad::rdf::MissMass(epic,"MissMassPhoto","{scat_ele,photon}");
  rad::rdf::MissP(epic,"MissPPhoto","{scat_ele,photon}");
  rad::rdf::MissE(epic,"MissEPhoto","{scat_ele,photon}");
  rad::rdf::MissPt(epic,"MissPtPhoto","{scat_ele,photon}");
  rad::rdf::MissPz(epic,"MissPzPhoto","{scat_ele,photon}");
  rad::rdf::MissTheta(epic,"MissThetaPhoto","{scat_ele,photon}");
  

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
  histo.Create<TH1D,double>({"hWelec",";W (electro miss mass) [GeV/c^{2}]",100,0,200.},{"W"});
  histo.Create<TH1D,double>({"hWhad",";W (hadro final state) [GeV/c^{2}]",100,0,200.},{"Whad"});

  histo.Create<TH1D,double>({"httop",";t(top vertex) [GeV^{2}]",100,-1,5},{"t_top"});
  histo.Create<TH1D,double>({"htbot",";t(bottom vertex) [GeV^{2}]",100,-1,5},{"t_bot"});
  histo.Create<TH1D,double>({"htptop",";t'(top vertex) [GeV^{2}]",100,-1,5},{"tp_top"});
  histo.Create<TH1D,double>({"htpbot",";t'(bottom vertex) [GeV^{2}]",100,-1,5},{"tp_bot"});
  
  histo.Create<TH1D,double>({"hcthCM",";cos(#theta_{CM})",100,-1,1},{"CM_CosTheta"});
  histo.Create<TH1D,double>({"hphCM",";#phi_{CM}",100,-TMath::Pi(),TMath::Pi()},{"CM_Phi"});
  
  histo.Create<TH1D,double>({"hMissMassPhoto",";M_{miss} [GeV/c^{2}]",1000,-10,10},{"MissMassPhoto"});
  histo.Create<TH1D,double>({"hmissPPhoto",";p_{miss}(e',Y)",1000,-200,200},{"MissPPhoto"});
  histo.Create<TH1D,double>({"hmissEPhoto",";E_{miss}(e',Y)",1000,-200,200},{"MissEPhoto"});
  histo.Create<TH1D,double>({"hmissPtPhoto",";p_{t,miss}(e',Y)",100,0,10},{"MissPtPhoto"});
  histo.Create<TH1D,double>({"hmissPzPhoto",";p_{z,miss}(e',Y)",1000,-100,100},{"MissPzPhoto"});
  histo.Create<TH1D,double>({"hmissThetaPhoto",";#theta_{miss}(e',Y)",100,-TMath::Pi(),TMath::Pi()},{"MissThetaPhoto"});
 
  histo.Create<TH1D,double>({"hMissMass",";M_{miss} [GeV/c^{2}]",1000,-10,10},{"MissMass"});
  histo.Create<TH1D,double>({"hmissP",";p_{miss}(e',Y)",1000,-100,100},{"MissP"});
  histo.Create<TH1D,double>({"hmissE",";E_{miss}(e',Y)",1000,-100,100},{"MissE"});
  histo.Create<TH1D,double>({"hmissPt",";p_{t,miss}(e',Y)",100,0,10},{"MissPt"});
  histo.Create<TH1D,double>({"hmissPz",";p_{z,miss}(e',Y)",1000,-100,100},{"MissPz"});
  histo.Create<TH1D,double>({"hmissTheta",";#theta_{miss}(e',Y)",100,-TMath::Pi(),TMath::Pi()},{"MissTheta"});
 
  //particle momenta
  histo.Create<TH1D,double>({"hpmag_elec",";p_{e'} [GeV/c]",100,-1,20},{"pmag[scat_ele]"});
  histo.Create<TH1D,double>({"heta_elec",";#eta_{e'} ",100,-10,10},{"eta[scat_ele]"});
  histo.Create<TH1D,double>({"htheta_elec",";#theta_{e'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"theta[scat_ele]"});
  histo.Create<TH1D,double>({"hphi_elec",";#phi_{e'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"phi[scat_ele]"});
  
  histo.Create<TH1D,double>({"hpmag_proton",";p_{p'} [GeV/c]",100,-1,300},{"pmag[proton]"});
  histo.Create<TH1D,double>({"heta_proton",";#eta_{p'} ",100,-10,10},{"eta[proton]"});
  histo.Create<TH1D,double>({"htheta_proton",";#theta_{p'} [rad]",1000,0,1},{"theta[proton]"});
  histo.Create<TH1D,double>({"hphi_proton",";#phi_{p'} [rad]",100,-TMath::Pi(),TMath::Pi()},{"phi[proton]"});
  histo.Create<TH2D,double,double>({"htheta_protonVmissp",";#theta_{p'} versus missing proton[rad]",100,0,0.01,100,0,0.01},{"theta[proton]","MissThetaPhoto"});
 
  //if I want a ROOT ttree
  epic.BookLazySnapshot("DVCS.root");
  gBenchmark->Start("processing");
			
  //save all histograms to file
  histo.File("DVCS_hists.root");

  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");
  gBenchmark->Stop("df");
  gBenchmark->Print("df");
  
}
