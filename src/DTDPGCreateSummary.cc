/*
 *  See header file for a description of this class.
 *
 *  $Date: 2013/07/04 17:16:16 $
 * revised U.G. 2007/07/24
 * Modified M.C Fouz. 2008/04/22 
 *          Updated for version 200 
 *          include the run & sector number on gif files
 *          2013: Local Trigger efficiency
 * Modified M.C Fouz. 2016/10/11 
 *          Updated for new DQM TM. Theta plots only 2 bins (Low & High) 
 *             (Before there were 8 for phi and theta)
 * Modified M.C Fouz. 2018/01/24
 *          Disabling the call to the RO utility due to the
 *          new uROS code, and forcing to read all wheels and sectors
 *          New program to be implemented in future to skeep Wh & S not
 *          present in data if needed
 *
 */


#include <UserCode/DTDPGAnalysis/src/DTDPGCreateSummary.h>
#include "UserCode/DTDPGAnalysis/src/DTROUtility.h"

// Framework
#include <FWCore/Utilities/interface/Exception.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EventSetup.h>

// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TPostScript.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TFolder.h"
#include "TText.h"
#include "boost/filesystem.hpp"

using namespace edm;
using namespace std;
using namespace boost::filesystem;

DTDPGCreateSummary::DTDPGCreateSummary(const edm::ParameterSet& ps) : myFile(0){

  LogVerbatim("DTDPGSummary") << "[DTDPGCreateSummary]: Constructor";

  myParameters = ps;
  myRunNumber = ps.getUntrackedParameter<int>("runNumber",0);
  ProcessDDUTrigger = ps.getUntrackedParameter<bool>("IncludeDDUTrigger","False");
  
  //  TCanvas c1("c1","",600,780);
  myCanvas = new TCanvas("c1", "",201,81,999,699);
  myCanvas->SetFillColor(0);

  gStyle->SetPalette(1);
  gStyle->SetStatColor(0);
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.29);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.30);
  gStyle->SetTitleH(0.1);
  gStyle->SetTitleW(0.97);


  gStyle->SetTitleFillColor(10);
}

DTDPGCreateSummary::~DTDPGCreateSummary(){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateSummary]: Destructor ";

}


void DTDPGCreateSummary::beginJob(){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateSummary]: BeginJob";
}


void DTDPGCreateSummary::beginRun(const edm::Run& run,const edm::EventSetup& context){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateSummary]: BeginRun";
  context.get<MuonGeometryRecord>().get(myMuonGeom);

}

void DTDPGCreateSummary::analyze(const edm::Event& e, const edm::EventSetup& context){
  
}


void DTDPGCreateSummary::endJob(){

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateSummary]: endJob";

  myFile = new TFile(myParameters.getUntrackedParameter<string>("rootFileName", "DTEfficiencyTest.root").c_str(), "READ");

  char cMainFolder[30];sprintf(cMainFolder,"DQMData/Run %d", myRunNumber);
  TFolder * mfolder=(TFolder*)myFile->Get(cMainFolder);
  if(!mfolder) {
    throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateSummary]: Folder = " << cMainFolder << " does not exist!!! Check the run number" << endl;     
  }
  
  myMainFolder.append(cMainFolder);
  myMainFolder.append("/DT/Run summary/");
  
  // FOR ROS DATA....
  //DTROUtility roUtil(myFile,mainFold);
  //map<int,std::vector<int> > mySectors = roUtil.getSectorMap();
  // FOR uROS DATA.... Forcing to look to all the sectors
     map<int,std::vector<int> > mySectors;
     vector<int> rosInTheFED;
     for(int isec=1;isec<15;isec++)
       rosInTheFED.push_back(isec);

     for (int iw=-2;iw<3;iw++)
       mySectors[iw] = rosInTheFED;
  // END FOR uROS DATA  Forcing to look to all the sectors

  edm::LogVerbatim ("DTDPGSummary") << "[DTDPGCreateSummary]: Creating summaries!";
 
  std::map<int,std::vector<int> >::const_iterator whIt  = mySectors.begin();
  std::map<int,std::vector<int> >::const_iterator whEnd = mySectors.end();
  for(;whIt!=whEnd;++whIt) {
    int wh = (*whIt).first;
    std::vector<int>::const_iterator secIt  = (*whIt).second.begin();
    std::vector<int>::const_iterator secEnd = (*whIt).second.end();
    for(;secIt!=secEnd;++secIt) {
      int sec = (*secIt);
      if (myParameters.getUntrackedParameter<bool>("DigiHistos", false))               { createDigiPlots(wh,sec); }
      if (myParameters.getUntrackedParameter<bool>("RecoHistos", false))               { createRecoPlots(wh,sec); }
      if (myParameters.getUntrackedParameter<bool>("ResoHistos", false))               { createResoPlots(wh,sec); }
//      if (myParameters.getUntrackedParameter<bool>("TimeSegHistos", false))            { createTimeSegmentPlots(wh,sec); }
      if (myParameters.getUntrackedParameter<bool>("EfficiencyHistos", false))         { createEfficiencyPlots(wh,sec); }
      if (myParameters.getUntrackedParameter<bool>("TestPulsesHistos", false))         { createTestPulsePlots(wh,sec); }
      if (myParameters.getUntrackedParameter<bool>("TriggerHistos", false) && sec<=12) { createTriggerPlots(wh,sec); }
    }
  }

  delete myCanvas;

  if(myFile) {
    myFile->Close();
    delete myFile;
  }
}

void DTDPGCreateSummary::createGifFile(string fileName, int wh, int sec) {

  stringstream gifBase, gifTag;
  gifBase << "Run" << myRunNumber << "/Wheel" << showpos << wh << "/Sect" << noshowpos << setfill('0')<< setw(2) << sec;
  try {
    create_directories(gifBase.str());
  } catch(const std::exception & ex) {
    throw cms::Exception("DTDPGCreateSummaryError")<< "[DTDPGCreateSummary]: Exception " << ex.what() << " thrown creating " << gifBase.str() << " directory" << endl;
  }
  
  gifTag << "_r" << myRunNumber << "_W" << wh << "_S" << sec;
  string gifFile = gifBase.str() + "/" + fileName + gifTag.str() + ".gif";
  myCanvas->Update();
  myCanvas->Print(gifFile.c_str());

}
  
  
void DTDPGCreateSummary::createDigiPlots(int wh, int sec) {
  bool  checkDraw=false;
  bool  checkDraw2=false;
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.29);
  //gStyle->SetStatY(0.98);
  gStyle->SetStatY(0.89);
  gStyle->SetStatH(0.30);
  gStyle->SetTitleH(0.1);
  gStyle->SetTitleW(0.97);
  gStyle->SetTitleBorderSize(2);



  gStyle->SetOptTitle();
  gStyle->SetOptStat();
  myCanvas->Clear();
  myCanvas->Divide(3,4);
  stringstream wheel;   wheel   << wh;
  stringstream sector;  sector  << sec;
  string digiFolder = myMainFolder + "01-Digi/Wheel" + wheel.str();
    
  // TimeBox Histos
  //vector<DTChamber*>::const_iterator chDigiIt = myMuonGeom->chambers().begin();
  //vector<DTChamber*>::const_iterator chDigiEnd = myMuonGeom->chambers().end();
  vector<const DTChamber*>::const_iterator chDigiIt = myMuonGeom->chambers().begin();
  vector<const DTChamber*>::const_iterator chDigiEnd = myMuonGeom->chambers().end();

  for (; chDigiIt != chDigiEnd; ++chDigiIt) {
    DTChamberId ch = (*chDigiIt)->id();
    if (ch.sector() == sec && ch.wheel()==wh) {
      stringstream station; station << ch.station();

      vector<const DTSuperLayer*>::const_iterator slIt = (*chDigiIt)->superLayers().begin(); 
      vector<const DTSuperLayer*>::const_iterator slEnd = (*chDigiIt)->superLayers().end();

      for(; slIt != slEnd; ++slIt) {
	DTSuperLayerId sl = (*slIt)->id();

	stringstream superLayer; superLayer << sl.superlayer();
	
	string histoName = digiFolder +  "/Sector" + sector.str() + "/Station" + station.str() +"/TimeBox_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str()+ "_SL" + superLayer.str();
	TH1F *histoTimeBox = (TH1F*) myFile->Get(histoName.c_str());

	if(histoTimeBox){
        checkDraw=true;
	  int pad = (ch.station() - 1)*3 + sl.superlayer();
	  myCanvas->cd(pad);
	  gStyle->SetStatW(0.19);
	  gStyle->SetOptStat(10);
	  histoTimeBox->GetXaxis()->SetLabelSize(0.08);
	  histoTimeBox->Draw();
	}
      }
    }
  }

  if(checkDraw) createGifFile("Timeboxes_filtered",wh,sec);    

  // Occupancy in Time Histos
  gStyle->SetOptStat();
  gStyle->SetStatW(0.29);
  myCanvas->Clear();
  myCanvas->Divide(3,4);
  TH1F *layerOccupancy[3][4];
  memset(layerOccupancy,0,12*sizeof(TH1F*));
    
  chDigiIt = myMuonGeom->chambers().begin();
  chDigiEnd = myMuonGeom->chambers().end();
  for (; chDigiIt != chDigiEnd; ++chDigiIt) {      
    DTChamberId ch = (*chDigiIt)->id();
    if (ch.sector() == sec && ch.wheel() == wh ) {
      stringstream station; station << ch.station();
      //bool found = false;  // Not being used and in 62X gives erors for being set but not used


      string histoNameOccupancy = digiFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/OccupancyInTimeHits_perCh_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str();     
      TH2F *histoOccupancy = (TH2F*) myFile -> Get(histoNameOccupancy.c_str());
      if(histoOccupancy){
      checkDraw2=true;
      //found = true;  // Not being used and in 62X gives erors for being set but not used
	  
	vector<const DTSuperLayer*>::const_iterator slIt = (*chDigiIt)->superLayers().begin(); 
	vector<const DTSuperLayer*>::const_iterator slEnd = (*chDigiIt)->superLayers().end();
      
	for(; slIt != slEnd; ++slIt) {
	  DTSuperLayerId sl = (*slIt)->id();
	  vector<const DTLayer*>::const_iterator lIt = (*slIt)->layers().begin(); 
	  vector<const DTLayer*>::const_iterator lEnd = (*slIt)->layers().end();
	  stringstream superLayer; superLayer << sl.superlayer();

	  for(; lIt != lEnd; ++lIt) {
	    DTLayerId layerId = (*lIt)->id();
           
	    // Create the new 1D histos
	    stringstream layer; layer << layerId.layer();
	    string histoTitleLayer= "In time occupancy: Station " + station.str() + "/SL" + superLayer.str();
	    string histoNameLayer = "Occupancy_MB" + station.str() + "SL" + superLayer.str() + "L" + layer.str();
	    layerOccupancy[sl.superlayer()-1][layerId.layer()-1]= new TH1F(histoNameLayer.c_str(),histoTitleLayer.c_str(),100, 0.,100.);

	    // Fill the 1D Histos
	    int iby=(sl.superlayer()-1)*4+layerId.layer();
	    int nbx= histoOccupancy->GetNbinsX() ;
	    float totentr=0;
	    for(int ibx=1;ibx<nbx;ibx++) {
	      float yval=histoOccupancy->GetBinContent(ibx,iby);
	      layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->SetBinContent(ibx,yval);
	      totentr+=yval;
	    }
	    layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->SetEntries(totentr);
	      
	    //Draw the histos
	    int Nstation = ch.station();
	    if (Nstation > 10) Nstation = 4;
	    int pad = (Nstation-1)*3 + sl.superlayer() ;
	    myCanvas->cd(pad);
	    layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->SetMinimum(0.);
	    layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->SetStats(0);
	    layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->GetXaxis()->SetLabelSize(0.08);
	    layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->GetYaxis()->SetLabelSize(0.06);
	    layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->SetLineColor( layerId.layer() );
	    if( layerId.layer() == 1 ) layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->Draw();
	    else layerOccupancy[sl.superlayer()-1][layerId.layer()-1]->Draw("same");
	  }  
	}  
      }  
    }
  }
    
  if(checkDraw2) createGifFile("Occupancies",wh,sec);

  for (int i=0;i<3;++i) {
    for (int j=0;j<4;++j) {
      if (layerOccupancy[i][j]) {
	delete layerOccupancy[i][j];
      }
    }
  }
    
}
  

void DTDPGCreateSummary::createTriggerPlots(int wh, int sec) {
  bool checkDraw=false;
  bool checkDraw2=false;
  bool checkDraw3=false;
  bool checkDraw4=false;
  bool checkDraw5=false;
  bool checkDraw6=false;
  bool checkDraw7=false;
  
  stringstream wheel; wheel << wh;
  stringstream sector; sector << sec;

  char Whname[5][20]={"Wm2","Wm1","W0","W1","W2"};
  stringstream wheelb;  
  wheelb << Whname[wh+2];

 
  std::vector<std::string> trigHwSources;
  //trigHwSources.push_back("DCC");
  trigHwSources.push_back("TM");
  if(ProcessDDUTrigger )trigHwSources.push_back("DDU");

  vector<string>::const_iterator hwSrcIt  = trigHwSources.begin();
  vector<string>::const_iterator hwSrcEnd = trigHwSources.end();

  for (; hwSrcIt != hwSrcEnd; ++hwSrcIt) {  
    // Trigger quality (Phi & Theta)
    myCanvas->Clear();
    myCanvas->Divide(2,4);

    //vector<DTChamber*>::const_iterator chTriggerIt = myMuonGeom->chambers().begin();
    //vector<DTChamber*>::const_iterator chTriggerEnd = myMuonGeom->chambers().end();
    vector<const DTChamber*>::const_iterator chTriggerIt = myMuonGeom->chambers().begin();
    vector<const DTChamber*>::const_iterator chTriggerEnd = myMuonGeom->chambers().end();
    string triggerFolder;
      //if ((*hwSrcIt)=="DCC") {
      if ((*hwSrcIt)=="TM") {
	//triggerFolder = myMainFolder + "03-LocalTrigger-DCC/Wheel" + wheel.str();
	triggerFolder = myMainFolder + "03-LocalTrigger-TM/Wheel" + wheel.str();
      }
      else {
	triggerFolder = myMainFolder + "04-LocalTrigger-DDU/Wheel" + wheel.str();
      }
    
    for (; chTriggerIt != chTriggerEnd; ++chTriggerIt) {
      DTChamberId ch = (*chTriggerIt)->id();
      if (ch.sector() == sec  && ch.wheel() == wh ) {      
	stringstream station; station << ch.station();
    //Changed to move to TM IN (06/10/2016)
    //string histoNameBXQual = triggerFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/LocalTriggerPhi/" + (*hwSrcIt) + "_BXvsQual_W" +  wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
	string histoNameBXQual = triggerFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/LocalTriggerPhiIn/" + (*hwSrcIt) + "_BXvsQual_In_W" +  wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
	TH1F *histoBXQual = (TH1F*) myFile -> Get(histoNameBXQual.c_str());
	if(histoBXQual){
        checkDraw=true;
	  myCanvas->cd( 2*ch.station() - 1 );
	  histoBXQual->GetXaxis()->SetTitle("Quality Flag, Phi trigger");
	  histoBXQual->GetYaxis()->SetTitle("BX id .");
	  histoBXQual->GetYaxis()->SetTitleSize(0.1);
	  histoBXQual->GetYaxis()->SetTitleOffset(0.4);
	  histoBXQual->GetYaxis()->SetRangeUser(-9.0,50.1);
	  histoBXQual->GetXaxis()->SetTitleSize(0.08);
	  histoBXQual->SetStats( 0 );
        histoBXQual->SetFillColor( 0 );
	  histoBXQual->GetXaxis()->SetLabelSize(0.12);
	  histoBXQual->GetYaxis()->SetLabelSize(0.08);
	  histoBXQual->Draw("colz");
	}
	histoNameBXQual = triggerFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/LocalTriggerTheta/" + (*hwSrcIt)+ "_ThetaBXvsQual_W" +  wheel.str() + "_Sec" + sector.str() + "_St" + station.str();    
	TH1F *histoBXQualTh = (TH1F*) myFile -> Get(histoNameBXQual.c_str());
	if(histoBXQualTh){
        checkDraw=true;
	  myCanvas->cd( 2*(ch.station()-1) + 2 );
	  histoBXQualTh->GetXaxis()->SetTitle("Quality Flag, Theta trigger");
	  histoBXQualTh->GetYaxis()->SetTitle("BX id .");
	  histoBXQualTh->GetYaxis()->SetTitleSize(0.1);
	  histoBXQualTh->GetYaxis()->SetTitleOffset(0.4);
	  histoBXQualTh->GetXaxis()->SetTitleSize(0.08);
       // After changes on DQM for the TM, only 2 bins in the histogram
       // for(int ib=1;ib<8;ib++) histoBXQualTh->GetXaxis()->SetBinLabel(ib," ");
       // histoBXQualTh->GetXaxis()->SetBinLabel(2,"L");
	 // histoBXQualTh->GetXaxis()->SetBinLabel(4,"H");
        histoBXQualTh->GetXaxis()->SetBinLabel(1,"L");
	  histoBXQualTh->GetXaxis()->SetBinLabel(2,"H");

	  histoBXQualTh->SetStats( 0 );
        histoBXQualTh->SetFillColor( 0 );
	  histoBXQualTh->GetXaxis()->SetLabelSize(0.12);
	  histoBXQualTh->GetYaxis()->SetLabelSize(0.08);
	  histoBXQualTh->Draw("colz");
	}
      } 
    }

    if(checkDraw) createGifFile("TriggQuality"+(*hwSrcIt),wh,sec);


    // Trig efficiency 2d plots 
    myCanvas->Clear();
    myCanvas->Divide(3,4);
    string histo2DTag[3] = {"TrigEffPosvsAnglePhi","TrigEffPosvsAngleHHHLPhi","TrackPosvsAngle"};

    chTriggerIt = myMuonGeom->chambers().begin();
    chTriggerEnd = myMuonGeom->chambers().end();
    for (; chTriggerIt != chTriggerEnd; ++chTriggerIt) {
      DTChamberId ch = (*chTriggerIt)->id();
    
      if (ch.sector() == sec  && ch.wheel() == wh ) {      
	stringstream station; station << ch.station();
	
	for (int i=0;i<3;++i){
	  string histoName2DEff = triggerFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/Segment/"+ (*hwSrcIt) + "_" + histo2DTag[i] + "_W" +  wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
	  TH2F *histo2DEff = (TH2F*) myFile -> Get(histoName2DEff.c_str());
	  if(histo2DEff){
          checkDraw2=true;
	    myCanvas->cd((ch.station()*3)-2+i);
          histo2DEff->SetFillColor( 0 );
	    //histo2DEff->SetStats( 0 );
          bool printStats = i==2;
          histo2DEff->SetStats( printStats );
	    histo2DEff->GetXaxis()->SetTitle("angle [deg]");
	    histo2DEff->GetYaxis()->SetTitle("position [cm] .");
	    histo2DEff->GetYaxis()->SetTitleSize(0.08);
	    histo2DEff->GetYaxis()->SetTitleOffset(0.6);
	    histo2DEff->GetYaxis()->SetLabelSize(0.08);
	    histo2DEff->GetXaxis()->SetTitleSize(0.075);
	    histo2DEff->GetXaxis()->SetLabelSize(0.075);
	    histo2DEff->GetXaxis()->SetTitleOffset(0.75);
	    histo2DEff->SetMinimum(0.);
	    histo2DEff->SetMaximum(1.);
	    //const char* option = i ? "box" : "colz";
          const char* option = i==2 ? "box" : "colz";
	    histo2DEff->Draw(option);
	  }
	}
      }
    }      

    if(checkDraw2) createGifFile("TriggEfficiencyPosvsAngle"+(*hwSrcIt),wh,sec);
  }



  if(ProcessDDUTrigger )
  {
  // Local PhiDDU triggerPlots
    string  LocalEfftriggerFolder = myMainFolder + "DTOfflineAnalyzer/";
  // Phi=======================
    myCanvas->Clear();
    myCanvas->Divide(4,4);
    string histoLTrigTagPhi[4] = {"hDDU_LocalTrackPosandTrig","hDDU_LocalTrackPosHHHLandTrig","hDDU_LocalTrackPos","hDDU_LocalTrackPosHHHL"};
    string histoLTrigTitPhi[4] = {"Local Phi Trig Eff (DDU)","Local Phi Trig Eff-HHHL (DDU)","Local Phi Track Pos >3hits","Local Phi Track Pos >6hits"};

    //vector<DTChamber*>::const_iterator chTrigIt = myMuonGeom->chambers().begin();
    //vector<DTChamber*>::const_iterator chTrigEnd = myMuonGeom->chambers().end();
    vector<const DTChamber*>::const_iterator chTrigIt = myMuonGeom->chambers().begin();
    vector<const DTChamber*>::const_iterator chTrigEnd = myMuonGeom->chambers().end();

    chTrigIt = myMuonGeom->chambers().begin();
    chTrigEnd = myMuonGeom->chambers().end();
    for (; chTrigIt != chTrigEnd; ++chTrigIt) {
      DTChamberId ch = (*chTrigIt)->id();
    
      if (ch.sector() == sec  && ch.wheel() == wh ) {      
	stringstream station; station << ch.station();
	
	for (int i=0;i<4;++i){
	  string histoNameLTrigEff = LocalEfftriggerFolder +  histoLTrigTagPhi[i] + "_"  +  wheelb.str() +"_S" + sector.str() + "_MB" + station.str();
	  TH1F *histoLTDEff = (TH1F*) myFile -> Get(histoNameLTrigEff.c_str());


	  if(histoLTDEff){
          checkDraw3=true;
          if(i<2)
          {// Compute the Efficiency
	       string histoName2LT = LocalEfftriggerFolder +  histoLTrigTagPhi[i+2] + "_"  +  wheelb.str() +"_S" + sector.str() + "_MB" + station.str();
	       TH1F *histo2LT = (TH1F*) myFile -> Get(histoName2LT.c_str());

             if(histo2LT)
               if(histo2LT->GetEntries()>0) histoLTDEff->Divide(histo2LT);
          }
          string histoNewTit=histoLTrigTitPhi[i] + " W"  +  wheel.str() +"_S" + sector.str() + "_MB" + station.str();
          histoLTDEff->SetTitle(histoNewTit.c_str());

	    myCanvas->cd((ch.station()*4)-3+i);
          histoLTDEff->SetFillColor( 0 );
	    //histoLTDEff->SetStats( 0 );
          //bool printStats = i==2;
          bool printStats = i>1;
          histoLTDEff->SetStats( printStats );
	    histoLTDEff->GetXaxis()->SetTitle("position [cm] .");
	    if(i<2)histoLTDEff->GetYaxis()->SetTitle("Efficiency");
	     else histoLTDEff->GetYaxis()->SetTitle("Entries");
	    histoLTDEff->GetYaxis()->SetTitleSize(0.08);
	    //histoLTDEff->GetYaxis()->SetTitleOffset(0.6);
	    histoLTDEff->GetYaxis()->SetTitleOffset(1.2);
	    histoLTDEff->GetYaxis()->SetLabelSize(0.08);
	    histoLTDEff->GetXaxis()->SetTitleSize(0.075);
	    histoLTDEff->GetXaxis()->SetLabelSize(0.075);
	    histoLTDEff->GetXaxis()->SetTitleOffset(0.75);
	    histoLTDEff->SetMinimum(0.);
	    if (i<2) histoLTDEff->SetMaximum(1.1);
	    histoLTDEff->Draw();
	  }
	}
      }
    }      

    if(checkDraw3) createGifFile("TriggLocalEfficiency",wh,sec);

  // Local PhiDDU triggerPlots
  // Theta  ===================
    myCanvas->Clear();
    myCanvas->Divide(4,4);
    string histoLTrigTagTheta[4] = {"hDDU_LocalTrackThetaPosandTrig","hDDU_LocalTrackThetaPosHandTrig","hDDU_LocalTrackThetaPos","hDDU_LocalTrackThetaPosH"};
    string histoLTrigTitTheta[4] = {"Local Theta Trig Eff (DDU)","Local Theta Trig Eff - H (DDU)","Local Theta Track Pos","Local Theta Track Pos 4hits"};

    int check=0;
    chTrigIt = myMuonGeom->chambers().begin();
    chTrigEnd = myMuonGeom->chambers().end();
    for (; chTrigIt != chTrigEnd; ++chTrigIt) {
      DTChamberId ch = (*chTrigIt)->id();
    
      if (ch.sector() == sec  && ch.wheel() == wh ) {      
	stringstream station; station << ch.station();
	
	for (int i=0;i<4;++i){
	  string histoNameLTrigEff = LocalEfftriggerFolder +  histoLTrigTagTheta[i] + "_"  +  wheelb.str() +"_S" + sector.str() + "_MB" + station.str();
	  TH1F *histoLTDEff = (TH1F*) myFile -> Get(histoNameLTrigEff.c_str());

	  if(histoLTDEff){
        checkDraw4=true;
        check++;
          if(i<2)
          {// Compute the Efficiency
	       string histoName2LT = LocalEfftriggerFolder +  histoLTrigTagTheta[i+2] + "_"  +  wheelb.str() +"_S" + sector.str() + "_MB" + station.str();
	       TH1F *histo2LT = (TH1F*) myFile -> Get(histoName2LT.c_str());

             if(histo2LT) 
               if(histo2LT->GetEntries()>0) histoLTDEff->Divide(histo2LT);

          }

          string histoNewTit=histoLTrigTitTheta[i] + " W"  +  wheel.str() +"_S" + sector.str() + "_MB" + station.str();
          histoLTDEff->SetTitle(histoNewTit.c_str());

	    myCanvas->cd((ch.station()*4)-3+i);
          histoLTDEff->SetFillColor( 0 );
	    //histoLTDEff->SetStats( 0 );
          //bool printStats = i==2;
          bool printStats = i>2;
          histoLTDEff->SetStats( printStats );
	    histoLTDEff->GetXaxis()->SetTitle("position [cm] .");
	    if(i<2)histoLTDEff->GetYaxis()->SetTitle("Efficiency");
	     else histoLTDEff->GetYaxis()->SetTitle("Entries");
	    histoLTDEff->GetYaxis()->SetTitleSize(0.08);
	    //histoLTDEff->GetYaxis()->SetTitleOffset(0.6);
	    histoLTDEff->GetYaxis()->SetTitleOffset(1.2);
	    histoLTDEff->GetYaxis()->SetLabelSize(0.08);
	    histoLTDEff->GetXaxis()->SetTitleSize(0.075);
	    histoLTDEff->GetXaxis()->SetLabelSize(0.075);
	    histoLTDEff->GetXaxis()->SetTitleOffset(0.75);
	    histoLTDEff->SetMinimum(0.);
	    if(i<2)histoLTDEff->SetMaximum(1.1);
	    histoLTDEff->Draw();
	  }
	}
      }
    }      

    if(checkDraw4) createGifFile("TriggLocalThetaEfficiency",wh,sec);
  }// END if(ProcessDDUTrigger )


  // DCC Only Plots
  //string triggerFolder = myMainFolder + "03-LocalTrigger-DCC/Wheel" + wheel.str();
  
  // TM Only Plots
  string triggerFolder = myMainFolder + "03-LocalTrigger-TM/Wheel" + wheel.str();

  // Trigger LUTs   
  myCanvas->Clear();
  if (sec==4 || sec==10) {
    myCanvas->Divide(2,5);
  } else {
    myCanvas->Divide(2,4);
  }
  //vector<DTChamber*>::const_iterator chTriggerIt = myMuonGeom->chambers().begin();
  //vector<DTChamber*>::const_iterator chTriggerEnd = myMuonGeom->chambers().end();
  vector<const DTChamber*>::const_iterator chTriggerIt = myMuonGeom->chambers().begin();
  vector<const DTChamber*>::const_iterator chTriggerEnd = myMuonGeom->chambers().end();
  for (; chTriggerIt != chTriggerEnd; ++chTriggerIt) {
    DTChamberId ch = (*chTriggerIt)->id();
    int scsec = ch.sector()>12 ? (ch.sector()==13 ? 4 : 10) : ch.sector();
    if ( (scsec == sec  && ch.wheel() == wh)) {
      int istat = ch.sector()>12 ? 5 : ch.station();
      stringstream station; station << ch.station();
      stringstream chsector; chsector << ch.sector();
      //string histoNameLut = triggerFolder + "/Sector" + chsector.str() + "/Station" + station.str() + "/Segment/DCC_PhitkvsPhitrig_W" +  wheel.str() + "_Sec" + chsector.str() + "_St" + station.str();
      string histoNameLut = triggerFolder + "/Sector" + chsector.str() + "/Station" + station.str() + "/Segment/TM_PhitkvsPhitrig_W" +  wheel.str() + "_Sec" + chsector.str() + "_St" + station.str();
      TH2F *histoLut = (TH2F*) myFile -> Get(histoNameLut.c_str());
      if(histoLut){
      checkDraw5=true;
	myCanvas->cd( 2*istat - 1 );
      histoLut->SetFillColor(0);
	histoLut->GetXaxis()->SetTitle("Trigger local position [cm]");
	histoLut->GetYaxis()->SetTitle("Segment local position [cm]");
	histoLut->SetStats( 0 );
	histoLut->GetXaxis()->SetLabelSize(0.08);
	histoLut->GetYaxis()->SetLabelSize(0.08);
	histoLut->Draw("box");
      }

      //histoNameLut = triggerFolder + "/Sector" + chsector.str() + "/Station" + station.str() + "/Segment/DCC_PhibtkvsPhibtrig_W" +  wheel.str() + "_Sec" + chsector.str() + "_St" + station.str();
      histoNameLut = triggerFolder + "/Sector" + chsector.str() + "/Station" + station.str() + "/Segment/TM_PhibtkvsPhibtrig_W" +  wheel.str() + "_Sec" + chsector.str() + "_St" + station.str();
      histoLut = (TH2F*) myFile -> Get(histoNameLut.c_str());
      if(histoLut && istat!=3){
      checkDraw5=true;
	myCanvas->cd( 2*istat );
	histoLut->SetFillColor(0);
	histoLut->GetXaxis()->SetTitle("Trigger local direction [deg]");
	histoLut->GetYaxis()->SetTitle("Segment local direction [deg]");
	histoLut->SetStats( 0 );
	histoLut->GetXaxis()->SetLabelSize(0.08);
	histoLut->GetYaxis()->SetLabelSize(0.08);
	histoLut->Draw("box");  
      }
    }
  }

  //if(checkDraw5) createGifFile("TriggLUTsDCC",wh,sec);
  if(checkDraw5) createGifFile("TriggLUTsTM",wh,sec);

  // Trigger LUTs Residuals   
  myCanvas->Clear();
  if (sec==4 || sec==10) {
    myCanvas->Divide(2,5);
  } else {
    myCanvas->Divide(2,4);
  }

  gStyle->SetOptFit(1);

  chTriggerIt = myMuonGeom->chambers().begin();
  gStyle->SetStatW(0.15);
  chTriggerEnd = myMuonGeom->chambers().end();
  for (; chTriggerIt != chTriggerEnd; ++chTriggerIt) {
    DTChamberId ch = (*chTriggerIt)->id();
    int scsec = ch.sector()>12 ? (ch.sector()==13 ? 4 : 10) : ch.sector();
    if (scsec == sec  && ch.wheel() == wh ) {
      int istat = ch.sector()>12 ? 5 : ch.station();
      stringstream station; station << ch.station();
      stringstream chsector; chsector << ch.sector();
      //string histoNameLutRes = triggerFolder + "/Sector" + chsector.str() + "/Station" + station.str() + "/Segment/DCC_PhiResidual_W" +  wheel.str() + "_Sec" + chsector.str() + "_St" + station.str();
      string histoNameLutRes = triggerFolder + "/Sector" + chsector.str() + "/Station" + station.str() + "/Segment/TM_PhiResidual_W" +  wheel.str() + "_Sec" + chsector.str() + "_St" + station.str();
      TH1F *histoLutRes = (TH1F*) myFile -> Get(histoNameLutRes.c_str());
      if(histoLutRes){
      checkDraw6=true;
	myCanvas->cd( 2*istat - 1 );
	gStyle->SetOptStat(1111);
	histoLutRes->SetFillColor(0);
	histoLutRes->GetXaxis()->SetTitle("Trigger - Segment local position [cm]");
	histoLutRes->GetXaxis()->SetLabelSize(0.08);
	histoLutRes->GetYaxis()->SetLabelSize(0.08);
	histoLutRes->Draw();
      }

      //histoNameLutRes = triggerFolder + "/Sector" + chsector.str() + "/Station" + station.str() + "/Segment/DCC_PhibResidual_W" +  wheel.str() + "_Sec" + chsector.str() + "_St" + station.str();
      histoNameLutRes = triggerFolder + "/Sector" + chsector.str() + "/Station" + station.str() + "/Segment/TM_PhibResidual_W" +  wheel.str() + "_Sec" + chsector.str() + "_St" + station.str();
      histoLutRes = (TH1F*) myFile->Get(histoNameLutRes.c_str());
      if(histoLutRes && istat!=3){
      checkDraw6=true;
	myCanvas->cd( 2*istat );
	histoLutRes->SetFillColor(0);
	histoLutRes->GetXaxis()->SetTitle("Trigger - Segment local direction [deg]");
	histoLutRes->GetXaxis()->SetLabelSize(0.08);
	histoLutRes->GetYaxis()->SetLabelSize(0.08);
	histoLutRes->Draw();
      }
    }
  }

  //if(checkDraw6) createGifFile("TriggLUTsResidualDCC",wh,sec);
  if(checkDraw6) createGifFile("TriggLUTsResidualTM",wh,sec);
  gStyle->SetOptFit(0);


  // Quality vs Hits
  myCanvas->Clear();
  myCanvas->Divide(2,2);
  
  chTriggerIt = myMuonGeom->chambers().begin();
  chTriggerEnd = myMuonGeom->chambers().end();
  for (; chTriggerIt != chTriggerEnd; ++chTriggerIt) {
    DTChamberId ch = (*chTriggerIt)->id();
    if (ch.sector() == sec  && ch.wheel() == wh ) {      
      stringstream wheel; wheel << ch.wheel();
      stringstream station; station << ch.station();
      stringstream sector; sector << ch.sector();
    
      //string histoNameHits = triggerFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/Segment/DCC_HitstkvsQualtrig_W" +  wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
      string histoNameHits = triggerFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/Segment/TM_HitstkvsQualtrig_W" +  wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
      TH2F *histoHits = (TH2F*) myFile -> Get(histoNameHits.c_str());
      if(histoHits){
      checkDraw7=true;
	myCanvas->cd( ch.station() );
	histoHits->GetXaxis()->SetTitle("Phi trigger quality");
	histoHits->GetYaxis()->SetTitle("Hits associated to phi segment");
	histoHits->SetFillColor( 0 );
	histoHits->SetStats(0);
	histoHits->Draw("box");
      }
    } 
  }
    
  //if(checkDraw7) createGifFile("TriggQualvsHitsDCC",wh,sec);
  if(checkDraw7) createGifFile("TriggQualvsHitsTM",wh,sec);

}

  
void DTDPGCreateSummary::createRecoPlots(int wh, int sec) {
  bool checkDraw=false;

  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.29);
  //gStyle->SetStatY(0.98);
  gStyle->SetStatY(0.89);
  gStyle->SetStatH(0.30);
  gStyle->SetTitleH(0.1);
  gStyle->SetTitleW(0.97);
  gStyle->SetTitleBorderSize(2);

  myCanvas->Clear(); 
  myCanvas->Divide(2,2);
  stringstream wheel; wheel << wh;
  stringstream sector; sector << sec;
  string recoFolder = myMainFolder + "02-Segments/Wheel" + wheel.str();

  //vector<DTChamber*>::const_iterator chRecoIt = myMuonGeom->chambers().begin();
  //vector<DTChamber*>::const_iterator chRecoEnd = myMuonGeom->chambers().end();
  vector<const DTChamber*>::const_iterator chRecoIt = myMuonGeom->chambers().begin();
  vector<const DTChamber*>::const_iterator chRecoEnd = myMuonGeom->chambers().end();

  for (; chRecoIt != chRecoEnd; ++chRecoIt) {
    DTChamberId ch = (*chRecoIt)->id();
    if (ch.sector() == sec  && ch.wheel() == wh ) {      
      stringstream station; station << ch.station();
      string histoName4DSeg;
      if(ch.station() < 4)
	histoName4DSeg = recoFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/h4DSegmXvsYInCham_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str(); 
      if(ch.station() == 4)
	histoName4DSeg = recoFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/h4DSegmXInCham_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str(); 
      //TProfile *histo4DSeg = (TProfile*) myFile -> Get(histoName4DSeg.c_str());  // CB fix here it is not a profile
      histoName4DSeg = recoFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/h4DSegmNHits_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str(); 
      TH1F *histo4DSegNHits = (TH1F*) myFile -> Get(histoName4DSeg.c_str());
      if(histo4DSegNHits){
      checkDraw=true;
	int pad = ch.station();
	myCanvas->cd(pad);
	histo4DSegNHits->GetXaxis()->SetLabelSize(0.08);
	histo4DSegNHits->GetYaxis()->SetLabelSize(0.07);
	histo4DSegNHits->Draw();

      TText L1;
      L1.SetTextSize(0.1);
      L1.SetTextColor(4);
      char MBstation[10];
      sprintf(MBstation,"MB%d", ch.station());
      L1.DrawText(14.0,0.02,MBstation);


      }
    }
  }

  if(checkDraw) createGifFile("RecoSegments",wh,sec);
 
}


    
  
void DTDPGCreateSummary::createResoPlots(int wh, int sec) {

  bool checkDraw=false;
 
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.29);
  //gStyle->SetStatY(0.98);
  gStyle->SetStatY(0.89);
  gStyle->SetStatH(0.30);
  gStyle->SetTitleH(0.1);
  gStyle->SetTitleW(0.97);
  gStyle->SetTitleBorderSize(2);

  myCanvas->Clear();
  myCanvas->Divide(3,4);
  stringstream wheel; wheel << wh;
  stringstream sector; sector << sec;
  string resoFolder = myMainFolder + "02-Segments/Wheel" + wheel.str();
  //vector<DTChamber*>::const_iterator chResoIt = myMuonGeom->chambers().begin();
  //vector<DTChamber*>::const_iterator chResoEnd = myMuonGeom->chambers().end();
  vector<const DTChamber*>::const_iterator chResoIt = myMuonGeom->chambers().begin();
  vector<const DTChamber*>::const_iterator chResoEnd = myMuonGeom->chambers().end();

  for (; chResoIt != chResoEnd; ++chResoIt) {
    DTChamberId ch = (*chResoIt)->id();
    if (ch.sector() == sec  && ch.wheel() == wh ) {
      stringstream station; station << ch.station();
      vector<const DTSuperLayer*>::const_iterator slIt = (*chResoIt)->superLayers().begin(); 
      vector<const DTSuperLayer*>::const_iterator slEnd = (*chResoIt)->superLayers().end();
      
      for(; slIt != slEnd; ++slIt) {
	DTSuperLayerId sl = (*slIt)->id();
	stringstream superLayer; superLayer << sl.superlayer();
	
	string histoNameResDist = resoFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/hResDist_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str() + "_SL" + superLayer.str(); 
	TH1F *histoResDist = (TH1F*) myFile -> Get(histoNameResDist.c_str());
	if(histoResDist){
        checkDraw=true;
	  int pad = (ch.station() - 1)*3 + sl.superlayer();
	  myCanvas->cd(pad);
	  histoResDist->GetXaxis()->SetTitle(" Hit residuals (cm) ");
	  // histo->FitPanel(1);
	  // histo->Fit("f1","R");
	  // g1 = new TF1("g1","gaus",-0.15,0.15);
	  // histo->Fit("gaus");	  
	  histoResDist->GetXaxis()->SetLabelSize(0.08);
	  histoResDist->GetYaxis()->SetLabelSize(0.07);
	  histoResDist->Draw();
	}
      }
    }
  }

  if(checkDraw) createGifFile("HitResiduals",wh,sec);
    
}

/*
void DTDPGCreateSummary::createTimeSegmentPlots(int wh, int sec) {

  bool checkDraw=false; 
  gStyle->SetStatX(0.99);
  gStyle->SetStatW(0.29);
  //gStyle->SetStatY(0.98);
  gStyle->SetStatY(0.89);
  gStyle->SetStatH(0.30);
  gStyle->SetTitleH(0.1);
  gStyle->SetTitleW(0.97);
  gStyle->SetTitleBorderSize(2);

  myCanvas->Clear();
  myCanvas->Divide(2,2);

  stringstream wheel; wheel << wh;
  stringstream sector; sector << sec;
  string  TimeSegmentFolder = myMainFolder + "DTOfflineAnalyzer/";
  char Whname[5][20]={"Wm2","Wm1","W0","W1","W2"};

  //vector<DTChamber*>::const_iterator chTimeSegIt = myMuonGeom->chambers().begin();
  //vector<DTChamber*>::const_iterator chTimeSegEnd = myMuonGeom->chambers().end();
  vector<const DTChamber*>::const_iterator chTimeSegIt = myMuonGeom->chambers().begin();
  vector<const DTChamber*>::const_iterator chTimeSegEnd = myMuonGeom->chambers().end();
  for (; chTimeSegIt != chTimeSegEnd; ++chTimeSegIt) {
    DTChamberId ch = (*chTimeSegIt)->id();
    if (ch.sector() == sec  && ch.wheel() == wh ) {
      stringstream station; station << ch.station();

      stringstream wheelb;  
      wheelb << Whname[wh+2];
      string histoNameTimeSeg = TimeSegmentFolder + "hTimeSeg_W" + wheelb.str() +  "_S" + sector.str() + "_MB" + station.str();
      TH1F *histoTimeSeg = (TH1F*) myFile -> Get(histoNameTimeSeg.c_str());

      if(histoTimeSeg)
      {
         checkDraw=true; 
	   int pad = ch.station();
	   myCanvas->cd(pad);

         double peak = histoTimeSeg->GetBinCenter(histoTimeSeg->GetMaximumBin());
         double TimSegMean = 0;

         try {
          TF1 fTimSeg("fTimeSeg","gaus");
          histoTimeSeg->Fit(&fTimSeg,"CQO","",peak-30,peak+30);
          TimSegMean = fTimSeg.GetParameter(1);
          histoTimeSeg->Fit(&fTimSeg,"CQO","",TimSegMean-30,TimSegMean+30);
          TimSegMean = fTimSeg.GetParameter(1);
         } catch(const cms::Exception&) {
            edm::LogError("DTDPGSummary") << "[DTDPGCreateSummary]: Error fitting TimeSegment for Wheel " << wh
                                       <<" Sector " << sec << " Station " << station;
         }

	   gStyle->SetOptStat(1);
	   gStyle->SetOptFit(1);
         histoTimeSeg->GetXaxis()->SetLabelSize(0.08);
         histoTimeSeg->GetYaxis()->SetLabelSize(0.07);

         histoTimeSeg->Draw();
 
         TText L1;
         L1.SetTextSize(0.1);
         L1.SetTextColor(4);
         char MBstation[10];
         sprintf(MBstation,"MB%d", ch.station());
         L1.DrawText(-200.0,0.02,MBstation);
      }
    }
  }

  if(checkDraw) createGifFile("TimeSegment",wh,sec);
}
*/


  
void DTDPGCreateSummary::createEfficiencyPlots(int wh, int sec) {

  stringstream wheel; wheel << wh;
  stringstream sector; sector << sec;
  
  string efficiencyFolder = myMainFolder + "DTEfficiencyTask/Wheel" + wheel.str() ; 
  //vector<DTChamber*>::const_iterator chEffIt = myMuonGeom->chambers().begin();
  //vector<DTChamber*>::const_iterator chEffEnd = myMuonGeom->chambers().end();
  vector<const DTChamber*>::const_iterator chEffIt = myMuonGeom->chambers().begin();
  vector<const DTChamber*>::const_iterator chEffEnd = myMuonGeom->chambers().end();

  gStyle->SetTitleBorderSize(2);

  for (; chEffIt != chEffEnd; ++chEffIt) {
    myCanvas->Clear();
    myCanvas->Divide(4,3);
    DTChamberId ch = (*chEffIt)->id();
    if (ch.sector() == sec  && ch.wheel() == wh ) {
      stringstream station; station << ch.station();
      bool checkDraw=false;

      vector<const DTSuperLayer*>::const_iterator slIt = (*chEffIt)->superLayers().begin(); 
      vector<const DTSuperLayer*>::const_iterator slEnd = (*chEffIt)->superLayers().end();

      for(; slIt != slEnd; ++slIt) {
	DTSuperLayerId sl = (*slIt)->id();
	if( !(ch.station() == 4 &&  sl.superlayer() ==2)) {
	  vector<const DTLayer*>::const_iterator lIt = (*slIt)->layers().begin(); 
	  vector<const DTLayer*>::const_iterator lEnd = (*slIt)->layers().end();
	  stringstream superLayer; superLayer << sl.superlayer();
	  for(; lIt != lEnd; ++lIt) {
	    DTLayerId layerId = (*lIt)->id();
	    stringstream layer; layer << layerId.layer();
	    
          string histoName1 = efficiencyFolder + "/Station" + station.str() + "/Sector" + sector.str() 
               + "/SuperLayer" + superLayer.str() +  "/hRecSegmOccupancy_W" + wheel.str() + "_St" + station.str() 
               + "_Sec" + sector.str() + "_SL" + superLayer.str()+ "_L" + layer.str(); 
          string histoName2 = efficiencyFolder + "/Station" + station.str() + "/Sector" + sector.str() 
               + "/SuperLayer" + superLayer.str() +  "/hEffUnassOccupancy_W" + wheel.str() + "_St" + station.str() 
               + "_Sec" + sector.str() + "_SL" + superLayer.str()+ "_L" + layer.str(); 
          TH1F *histoAll = (TH1F*) myFile -> Get(histoName1.c_str());
          TH1F *histoEff = NULL ;
          if(myFile -> Get(histoName2.c_str()))
             histoEff = (TH1F*) myFile -> Get(histoName2.c_str())->Clone("histoEff");


	    if(histoEff){
            checkDraw=true;

            int nb1=histoAll->GetNbinsX();
            int nb2=histoEff->GetNbinsX();
            if(nb1 != nb2){ throw cms::Exception("DTDPGCreateSummaryError") << "[DTDPGCreateSummary]: ERROR : mistmach on number of bins" << endl;}
            for(int ibin=1;ibin<nb1;ibin++)
            {
              float xpos1=histoEff->GetBinContent(ibin);
              float xpos2=histoAll->GetBinContent(ibin);
              if(xpos2>0)histoEff->SetBinContent(ibin,xpos1/xpos2);
               else histoEff->SetBinContent(ibin,0.);
            }

	      int pad = (sl.superlayer() - 1)*4 + layerId.layer();
	      myCanvas->cd(pad);
	      histoEff->SetMinimum(0.);
	      histoEff->SetMaximum(1.1);
	      histoEff->GetXaxis()->SetLabelSize(0.07);
	      histoEff->GetYaxis()->SetLabelSize(0.07);
	      histoEff->SetStats(0);
	      histoEff->Draw();
	      // Look for DeadCells > Eff<2% to avoid some noise entries
	      TH1F *histoDead=(TH1F*)histoEff->Clone();
                                                                                                                     
	      string histoName = myMainFolder + "DTEfficiencyTask/Wheel" + wheel.str()  + "/Station" + station.str()
		+ "/Sector" + sector.str() + "/SuperLayer" + superLayer.str()
		+ "/hRecSegmOccupancy_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str() + "_SL" + superLayer.str() + "_L"
		+ layer.str();
	      // To check that there are tracks passing for cells and avoid not existing cells
	      TH1F *histoCheck = (TH1F*) myFile -> Get(histoName.c_str());
	      int nb=histoEff->GetNbinsX();
	      for(int ib=1;ib<nb+1;ib++)
		{
		  if(histoEff->GetBinContent(ib)<0.02 && histoCheck->GetBinContent(ib)>3)
		    histoDead->SetBinContent(ib,1.3);
		  else
		    histoDead->SetBinContent(ib,0.0);
                                                                                                                     
		}
	      histoDead->SetFillColor(5);
	      histoDead->Draw("histsame");
	      histoEff->Draw("same");
	    }
	  }
	}
      }
      if(checkDraw){
	stringstream name; name << "Cell_efficiency_MB" << ch.station();
	createGifFile(name.str(),wh,sec);
      }
    }
  }
    
  // Segment reconstruction efficiency vs X
  efficiencyFolder = myMainFolder + "01-DTChamberEfficiency/Wheel" + wheel.str();
  myCanvas->Clear();
  myCanvas->Divide(2,4);

  bool checkDraw=false;
  chEffIt = myMuonGeom->chambers().begin();
  chEffEnd = myMuonGeom->chambers().end();
  for (; chEffIt != chEffEnd; ++chEffIt) {
    DTChamberId ch = (*chEffIt)->id();
    if (ch.sector() == sec  && ch.wheel() == wh ) {      
      stringstream station; station << ch.station();
      string histoName = efficiencyFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/xEfficiency_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str();
      TH1F *histoEffX = (TH1F*) myFile -> Get(histoName.c_str());
      if(histoEffX){
      checkDraw=true;
	myCanvas->cd( (ch.station()-1)*2 +1 );
        string histoTitle=  "Segment reco eff.vs X: Station " + station.str();
	histoEffX->SetTitle(histoTitle.c_str());
	histoEffX->GetXaxis()->SetTitle(" X position in chamber (cm) ");
	histoEffX->GetXaxis()->SetLabelSize(0.08);
	histoEffX->GetYaxis()->SetLabelSize(0.08);
	histoEffX->Draw();
      }
    }
  }
    
  // Segment reco efficiency vs Y
  chEffIt = myMuonGeom->chambers().begin();
  chEffEnd = myMuonGeom->chambers().end();
  for (; chEffIt != chEffEnd; ++chEffIt) {
    DTChamberId ch = (*chEffIt)->id();
    if (ch.sector() == sec  && ch.wheel() == wh ) {      
      stringstream station; station << ch.station();
      string histoNameEffY = efficiencyFolder + "/Sector" + sector.str() + "/Station" + station.str() + "/yEfficiency_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str();
      TH1F *histoEffY = (TH1F*) myFile -> Get(histoNameEffY.c_str());
      if(histoEffY){
      checkDraw=true;
	myCanvas->cd( (ch.station()-1)*2 +2 );
        string histoTitle=  "Segment reco eff.vs Y: Station " + station.str();
	histoEffY->SetTitle(histoTitle.c_str());
	histoEffY->GetXaxis()->SetTitle(" Y position in chamber (cm) ");
	gStyle->SetOptStat(10);
	histoEffY->GetXaxis()->SetLabelSize(0.08);
	histoEffY->GetYaxis()->SetLabelSize(0.08);
	histoEffY->Draw();
      }
    }
  }

  if(checkDraw) createGifFile("SegmentEfficiency",wh,sec);

}

void DTDPGCreateSummary::createTestPulsePlots(int wh, int sec) {

  gStyle->SetOptStat();
  myCanvas->Clear();
  myCanvas->Divide(4,3);
  stringstream wheel; wheel << wh;
  stringstream sector; sector << sec;

  //vector<DTChamber*>::const_iterator chTPIt = myMuonGeom->chambers().begin();
  //vector<DTChamber*>::const_iterator chTPEnd = myMuonGeom->chambers().end();
  vector<const DTChamber*>::const_iterator chTPIt = myMuonGeom->chambers().begin();
  vector<const DTChamber*>::const_iterator chTPEnd = myMuonGeom->chambers().end();
  for (; chTPIt != chTPEnd; ++chTPIt) {
    DTChamberId ch = (*chTPIt)->id();
    if (ch.sector() == sec  && ch.wheel() == wh ) {
      stringstream station; station << ch.station();
      bool found = false;
      vector<const DTSuperLayer*>::const_iterator slIt = (*chTPIt)->superLayers().begin(); 
      vector<const DTSuperLayer*>::const_iterator slEnd = (*chTPIt)->superLayers().end();
      // Loop over the SLs
      for(; slIt != slEnd; ++slIt) {
	DTSuperLayerId sl = (*slIt)->id();
	vector<const DTLayer*>::const_iterator lIt = (*slIt)->layers().begin(); 
	vector<const DTLayer*>::const_iterator lEnd = (*slIt)->layers().end();
	stringstream superLayer; superLayer << sl.superlayer();
	// Loop over the Ls
	for(; lIt != lEnd; ++lIt) {
	  DTLayerId layerId = (*lIt)->id();
	  stringstream layer; layer << layerId.layer();
	  string testPulsesFolder = myMainFolder + "DTTestPulsesTask/Wheel" + wheel.str();
	  string histoName = testPulsesFolder + "/Station" + station.str() + "/Sector" + sector.str() + "/SuperLayer" + superLayer.str() +  "/TPProfile/TestPulses2D_W" + wheel.str() + "_St" + station.str() + "_Sec" + sector.str() + "_SL" + superLayer.str() + "_L" + layer.str();
	  TProfile *histoTP = (TProfile*) myFile->Get(histoName.c_str()); // CB fix here
	  if(histoTP){
	    found = true;
	    int pad = (sl.superlayer() - 1)*4 + layerId.layer();
	    myCanvas->cd(pad);
	    histoTP->Draw();
	  }
	}
      }
      if(found){
	myCanvas->Update();
      }
    }
  }
}

