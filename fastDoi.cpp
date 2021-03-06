// compile with
// g++ -o ../build/fastDoi fastDoi.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas


#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TFormula.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphDelaunay.h"
#include "TVector.h"
#include "TNamed.h"
#include "TPaveLabel.h"
#include "THStack.h"
#include "TFitResult.h"
#include "TMatrixD.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

// #include "./libraries/CrystalStructs.h"     // Crystal_t , detector_t
// #include "./libraries/Calibration.h"        // readTaggingData , readCalibration , setWandZcuts
// #include "./libraries/Utilities.h"          // read_directory , invert_a_matrix
// #include "./libraries/Extract.h"            // extractWithEMG , extractCTR , FindSmallestInterval

// forward declaration of usage info output
void usage();

// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}


//----------------//
//  MAIN PROGRAM  //
//----------------//
int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout	<< "Usage: " << std::endl << argv[0] ;
    usage();
    return 1;
  }

  std::stringstream streamCommand;
  for(int i=0 ; i < argc; i++)
  {
    streamCommand << argv[i] << " ";
  }

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);

  // default args
  std::string calibrationFileNames = "";
  std::string inputFolderName = "./";
  std::string inputFilePrefix = "";
  std::string outputFileName = "outputBareboneFile.root";
  float enMin = 4500.0;
  float enMax = 9000.0;
  float fxMin = 1.0;
  float fxMax = 2.4;
  float fyMin = 1.0;
  float fyMax = 1.7;
  float point = 0.5;
  float wmin = 0.45;
  float wmax = 0.85;

  // parse command line arguments
  static struct option longOptions[] =
  {
      { "folder", required_argument, 0, 0 },
      { "prefix", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "enMin", required_argument, 0, 0 },
      { "enMax", required_argument, 0, 0 },
      { "fxMin", required_argument, 0, 0 },
      { "fxMax", required_argument, 0, 0 },
      { "fyMin", required_argument, 0, 0 },
      { "fyMax", required_argument, 0, 0 },
      { "point", required_argument, 0, 0 },
      { "wmin", required_argument, 0, 0 },
      { "wmax", required_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "f:p:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		else if (c == 'f'){
      inputFolderName = (char *)optarg;
    }
    else if (c == 'p'){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 0){
      inputFolderName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      outputFileName = (char *)optarg;
    }

    else if (c == 0 && optionIndex == 3){
      enMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      enMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 5){
      fxMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      fxMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      fyMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      fyMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      point = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      wmin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      wmax = atof((char *)optarg);
    }


		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}



  // check if required are given and files actually exists
  // first, input given and not empty
  if(inputFilePrefix == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide the prefix of input files!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }





  //----------------------------------//
  // GET INPUT FILES(S)               //
  //----------------------------------//
  // read file in dir
  std::cout << std::endl;
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << "|         ANALYSIS FILES                 |" << std::endl;
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << std::endl;
  // get input files list
  std::vector<std::string> v;
  read_directory(inputFolderName, v);
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;
  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFilePrefix.size(),inputFilePrefix))
    {
      listInputFiles.push_back(inputFolderName + "/" + v[i]);
    }
  }
  // check if it's empty
  if(listInputFiles.size() == 0)
  {
    std::cout << std::endl;
    std::cout << "ERROR! Some input files do not exists! Aborting." << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }


  //----------------------------------------------------------//
  //  Get TChain from input TTree files                       //
  //----------------------------------------------------------//
  TChain* tree = new TChain("adc");  // create the input tchain and the analysis ttree
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    std::cout << "Adding file " << listInputFiles[i] << std::endl;
    tree->Add(listInputFiles[i].c_str());
  }
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << std::endl;
  std::vector<int> detector_channels;
  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  // fill a vector with the leaves names
  for(int i = 0 ; i < nLeaves ; i++)
  {
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  // count the entries that start with "ch"
  int numOfCh = 0;
  std::string ch_prefix("ch");
  std::string t_prefix("t");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
    {
      numOfCh++;
      detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
    }
  }
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;

  //set variables and branches
  ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous
  UShort_t      *charge;
  Float_t      *timeStamp;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      **bCharge;
  TBranch      **btimeStamp;
  charge = new UShort_t[numOfCh];
  timeStamp = new Float_t[numOfCh];
  bCharge = new TBranch*[numOfCh];
  btimeStamp = new TBranch*[numOfCh];
  // set branches for reading the input files
  tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    //empty the stringstreams
    std::stringstream sname;
    sname << "ch" << detector_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
    sname.str("");
    sname << "t" << detector_channels[i];
    tree->SetBranchAddress(sname.str().c_str(),&timeStamp[detector_channels[i]],&btimeStamp[detector_channels[i]]);
    sname.str("");
  }

  // sim data
  Float_t RealX;
  Float_t RealY;
  Float_t RealZ;
  TBranch *bRealX;
  TBranch *bRealY;
  TBranch *bRealZ;
  tree->SetBranchAddress("RealX", &RealX, &bRealX);
  tree->SetBranchAddress("RealY", &RealY, &bRealY);
  tree->SetBranchAddress("RealZ", &RealZ, &bRealZ);

  // std::cout << "Entries = " << tree->GetEntries() << std::endl;

  float xmppc[16] = {-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};
  float ymppc[16] = {-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8};

  TH2F *scatter = new TH2F("Flood Map","Flood Map",1000,-7,7,1000,-7,7);
  TH2F *scatterReal = new TH2F("Flood Map Real","Flood Map Real",1000,-7,7,1000,-7,7);
  TH2F *Rz_w = new TH2F("Rz_w","Rz_w",100,0,1,25,0,16);
  TH1F *spectrum = new TH1F("spectrum","spectrum",500,0,10000);
  TH1F *H_realZ = new TH1F("H_realZ","H_realZ",100,-8,8);
  TH1F *H_w = new TH1F("H_w","H_w",500,0,1);

  TH3I *scatter3D = new TH3I("Flood Map 3D","Flood Map 3D",100,-7,7,100,-7,7,100,0,1);

  TH1F *spectrum1Cry = new TH1F("spectrum1Cry","spectrum 1Cry",500,0,10000);
  TH2F *scatter1Cry = new TH2F("Flood Map1Cry","Flood Map 1Cry",1000,-7,7,1000,-7,7);
  TH2F *scatterReal1Cry = new TH2F("Flood Map Real 1Cry","Flood Map Real 1Cry",1000,-7,7,1000,-7,7);
  TH3I *scatter3D1Cry = new TH3I("Flood Map 3D 1Cry","Flood Map 3D 1Cry",100,-7,7,100,-7,7,100,0,1);

  TH1F *ctr = new TH1F("ctr","ctr",200,-2e-9,2e-9);
  TH1F *ctr_avg= new TH1F("ctr_avg","ctr_avg",200,-2e-9,2e-9);
  TH1F *ctr_doi= new TH1F("ctr_doi","ctr_doi",200,-2e-9,2e-9);

  TH1F *ctr_n[8];
  ctr_n[0] = new TH1F("ctr_n0","ctr_n0",200,0,15e-9);
  ctr_n[1] = new TH1F("ctr_n1","ctr_n1",200,0,15e-9);
  ctr_n[2] = new TH1F("ctr_n2","ctr_n2",200,0,15e-9);
  ctr_n[3] = new TH1F("ctr_n3","ctr_n3",200,0,15e-9);
  ctr_n[4] = new TH1F("ctr_n4","ctr_n4",200,0,15e-9);
  ctr_n[5] = new TH1F("ctr_n5","ctr_n5",200,0,15e-9);
  ctr_n[6] = new TH1F("ctr_n6","ctr_n6",200,0,15e-9);
  ctr_n[7] = new TH1F("ctr_n7","ctr_n7",200,0,15e-9);


  TH2F *ctrVSw = new TH2F("ctrVSw","ctrVSw",100,0,1,200,-2e-9,2e-9);
  TH2F *ctrVSw_n[8] ;
  ctrVSw_n[0] = new TH2F("ctrVSw_n0","ctrVSw_n0",100,0,1,200,0e-9,10e-9);
  ctrVSw_n[1] = new TH2F("ctrVSw_n1","ctrVSw_n1",100,0,1,200,0e-9,10e-9);
  ctrVSw_n[2] = new TH2F("ctrVSw_n2","ctrVSw_n2",100,0,1,200,0e-9,10e-9);
  ctrVSw_n[3] = new TH2F("ctrVSw_n3","ctrVSw_n3",100,0,1,200,0e-9,10e-9);
  ctrVSw_n[4] = new TH2F("ctrVSw_n4","ctrVSw_n4",100,0,1,200,0e-9,10e-9);
  ctrVSw_n[5] = new TH2F("ctrVSw_n5","ctrVSw_n5",100,0,1,200,0e-9,10e-9);
  ctrVSw_n[6] = new TH2F("ctrVSw_n6","ctrVSw_n6",100,0,1,200,0e-9,10e-9);
  ctrVSw_n[7] = new TH2F("ctrVSw_n7","ctrVSw_n7",100,0,1,200,0e-9,10e-9);
  int n_ids[8] = {5,6,7,9,11,13,14,15};

  TProfile *profile_n[8] ;
  // profile__n[0] = new TH2F("profile__n_n0","ctrVSw_n0",100,0,1,200,0e-9,10e-9);
  // profile__n[1] = new TH2F("profile__n_n1","ctrVSw_n1",100,0,1,200,0e-9,10e-9);
  // profile__n[2] = new TH2F("profile__n_n2","ctrVSw_n2",100,0,1,200,0e-9,10e-9);
  // profile__n[3] = new TH2F("profile__n_n3","ctrVSw_n3",100,0,1,200,0e-9,10e-9);
  // profile__n[4] = new TH2F("profile__n_n4","ctrVSw_n4",100,0,1,200,0e-9,10e-9);
  // profile__n[5] = new TH2F("profile__n_n5","ctrVSw_n5",100,0,1,200,0e-9,10e-9);
  // profile__n[6] = new TH2F("profile__n_n6","ctrVSw_n6",100,0,1,200,0e-9,10e-9);
  // profile__n[7] = new TH2F("profile__n_n7","ctrVSw_n07",100,0,1,200,0e-9,10e-9);

  TF1 *line_n[8] ;
  line_n[0] = new TF1("line_n","pol1",0,1);
  line_n[1] = new TF1("line_n","pol1",0,1);
  line_n[2] = new TF1("line_n","pol1",0,1);
  line_n[3] = new TF1("line_n","pol1",0,1);
  line_n[4] = new TF1("line_n","pol1",0,1);
  line_n[5] = new TF1("line_n","pol1",0,1);
  line_n[6] = new TF1("line_n","pol1",0,1);
  line_n[7] = new TF1("line_n","pol1",0,1);



  // std::stringstream var;
  //
  // for(int i = 0 ; i < 16 ; i++)
  // {
  //   var << "ch"
  // }


  // TList *formulasAnalysis = new TList();
  // std::vector<Crystal_t> crystal;
  //
  // for(unsigned int i = 0 ; i < calibrationFile.size() ; i++)
  // {
  //   readCalibration(calibrationFile[i],       // this calib file
  //                   tree,                     // input TChain (same for everyone)
  //                   formulasAnalysis,         // TList of all TTreeFormula
  //                   crystal);                 // structure of all crystals found in all calib lifes
  //
  //
  // }
  // // optionally set w and z limits, and write values into crystal struct
  // // setWandZcuts(crystal);
  //
  //
  // // list the crystals with calibration data found
  // std::cout << "Calibration data found for crystals: " << std::endl;
  // for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  // {
  //   if(crystal[i].accepted)
  //   {
  //     std::cout << crystal[i].number << std::endl;
  //   }
  // }
  //
  //
  //
  //
  // MAIN LOOP
  long long int counter = 0;
  // tree->SetNotify(formulasAnalysis);
  long long int neventAnalysis = tree->GetEntries();
  std::cout << "Total number of events = " << neventAnalysis << std::endl;
  long int goodEventsAnalysis = 0;
  long int counterAnalysis = 0;
  for (long long int i=0;i<neventAnalysis;i++)
  {

    tree->GetEvent(i);              //read complete accepted event in memory

    float totalCharge = 0;
    //find max charge
    // and calc sum
    float maxCharge = 0;
    int chID = -1;
    for(int iCh = 0 ; iCh < 16 ; iCh++)
    {
      if(charge[iCh] > maxCharge)
      {
        maxCharge = charge[iCh];
        chID = iCh;
      }
      totalCharge += charge[iCh];
    }

    // for(int iCh = 0 ; iCh < 16 ; iCh++)
    // {
    //   totalCharge += charge[iCh];
    // }

    float FloodX = 0;
    float FloodY = 0;
    float FloodZ = 0;

    for(int iCh = 0 ; iCh < 16 ; iCh++)
    {
      FloodX += charge[iCh]*xmppc[iCh];
      FloodY += charge[iCh]*ymppc[iCh];
    }

    FloodX = FloodX/totalCharge;
    FloodY = FloodY/totalCharge;
    FloodZ = maxCharge/totalCharge;

    spectrum->Fill(totalCharge);
    scatter->Fill(FloodX,FloodY);
    scatter3D->Fill(FloodX,FloodY,FloodZ);
    scatterReal->Fill(RealX,RealY);

    // select one "peak in 2d"
    if( (FloodX > fxMin ) && (FloodX < fxMax ) && (FloodY > fyMin ) && (FloodY < fyMax )  )
    {
      spectrum1Cry->Fill(totalCharge);
      scatter1Cry->Fill(FloodX,FloodY);
      scatter3D1Cry->Fill(FloodX,FloodY,FloodZ);
      scatterReal1Cry->Fill(RealX,RealY);
      // select photopeak
      if( (totalCharge > enMin ) && (totalCharge < enMax ) )
      {
        goodEventsAnalysis++;
        H_realZ->Fill(RealZ);
        H_w->Fill(FloodZ);
        Rz_w->Fill(FloodZ,(7.5 - RealZ) );
        ctr->Fill(timeStamp[chID]);
        ctrVSw->Fill(FloodZ,timeStamp[chID]);
        for(int iN = 0; iN < 8 ; iN++)
        {
          ctrVSw_n[iN]->Fill(FloodZ,timeStamp[n_ids[iN]]);
          ctr_n[iN]->Fill(timeStamp[n_ids[iN]]);
        }
      }
    }







    // LOOP ON ALL ACCEPTED CRYSTALS
    // I.E. ALL CRYSTALS IN THE CALIBRATION FILES
    // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    // {
    //   // do whatever you want
    //   // in particular, some tips:
    //   // 1) if(crystal[iCry].FormulaTagAnalysis->EvalInstance()){} --> this will select events in the photopeak of the external reference crystal, if there is one. Notice that each crystal can have its own FormulaTagAnalysis, since they can come from different calibration runs
    //   // 2) if(crystal[iCry].FormulaAnalysis->EvalInstance()) --> this will select events in a crystal (geometrical position in u-v-w and photopeak condition)
    //   // 3) to calculate w, if needed
    //   // float FloodZ = calculateFloodZ(charge,crystal[iCry]);
    //   // notice that the detectorSaturation values are written directly into the crystal struct, because they can be acquired in different conditions (but it has to be the same photodetector!)
    //
    //   //example here, selecting photopeak events in ext ref, photopeak events the crystal(s) found in calibration file(s), and calculating FloodZ for each
    //   if(crystal[iCry].accepted)
    //   {
    //     if(crystal[iCry].FormulaTagAnalysis->EvalInstance())
    //     {
    //       if(crystal[iCry].FormulaAnalysis->EvalInstance())
    //       {
    //         goodEventsAnalysis++;
    //         float FloodZ = calculateFloodZ(charge,crystal[iCry]);
    //         // std::cout << FloodZ << std::endl;
    //       }
    //     }
    //   }
    // }

    //LOOP COUNTER
    counterAnalysis++;
    int perc = ((100*counterAnalysis)/neventAnalysis);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }

  std::cout << "Good events = " << goodEventsAnalysis << std::endl;

  gStyle->SetOptStat(0);
  TCanvas *c_gr = new TCanvas("c_pfx","",200,10,700,500);
  c_gr->SetGrid();

  TProfile *pfx =  Rz_w->ProfileX("pfx");
  pfx->SetTitle("");
  pfx->GetXaxis()->SetTitle("w");
  pfx->GetYaxis()->SetTitle("z [mm]");

  std::cout << "value = "<<pfx->GetBinError(pfx->FindBin(point))*2.355 << std::endl;






  c_gr->cd();
  pfx->Draw();

  TProfile *pfx_central =  ctrVSw->ProfileX("pfx_central");
  pfx_central->SetTitle("");
  pfx_central->GetXaxis()->SetTitle("w");
  pfx_central->GetYaxis()->SetTitle("t [s]");

  TF1 * line_central = new TF1("line_n","pol1",0,1);
  pfx_central->Fit(line_central,"","",wmin,wmax);




  for(int iN = 0; iN < 8 ; iN++)
  {
    std::stringstream pname;
    pname << "pfx_n" << iN;
    profile_n[iN] = ctrVSw_n[iN]->ProfileX(pname.str().c_str());
    profile_n[iN]->Fit(line_n[iN],"","",wmin,wmax);

  }


  counter = 0;
  // tree->SetNotify(formulasAnalysis);
  // long long int neventAnalysis = tree->GetEntries();
  // std::cout << "Total number of events = " << neventAnalysis << std::endl;
  goodEventsAnalysis = 0;
  counterAnalysis = 0;
  for (long long int i=0;i<neventAnalysis;i++)
  {

    tree->GetEvent(i);              //read complete accepted event in memory

    float totalCharge = 0;
    //find max charge
    // and calc sum
    float maxCharge = 0;
    int chID = -1;
    for(int iCh = 0 ; iCh < 16 ; iCh++)
    {
      if(charge[iCh] > maxCharge)
      {
        maxCharge = charge[iCh];
        chID = iCh;
      }
      totalCharge += charge[iCh];
    }

    // for(int iCh = 0 ; iCh < 16 ; iCh++)
    // {
    //   totalCharge += charge[iCh];
    // }

    float FloodX = 0;
    float FloodY = 0;
    float FloodZ = 0;

    for(int iCh = 0 ; iCh < 16 ; iCh++)
    {
      FloodX += charge[iCh]*xmppc[iCh];
      FloodY += charge[iCh]*ymppc[iCh];
    }

    FloodX = FloodX/totalCharge;
    FloodY = FloodY/totalCharge;
    FloodZ = maxCharge/totalCharge;

    // spectrum->Fill(totalCharge);
    // scatter->Fill(FloodX,FloodY);
    // scatter3D->Fill(FloodX,FloodY,FloodZ);
    // scatterReal->Fill(RealX,RealY);

    // select one "peak in 2d"
    if( (FloodX > fxMin ) && (FloodX < fxMax ) && (FloodY > fyMin ) && (FloodY < fyMax )  )
    {
      // spectrum1Cry->Fill(totalCharge);
      // scatter1Cry->Fill(FloodX,FloodY);
      // scatter3D1Cry->Fill(FloodX,FloodY,FloodZ);
      // scatterReal1Cry->Fill(RealX,RealY);
      // select photopeak
      if(FloodZ > wmin && FloodZ < wmax)
      {
        float m_w = (wmin+wmax)/2.0;
        if( (totalCharge > enMin ) && (totalCharge < enMax ) )
        {
          // goodEventsAnalysis++;
          // H_realZ->Fill(RealZ);
          // H_w->Fill(FloodZ);
          // Rz_w->Fill(FloodZ,(7.5 - RealZ) );
          // ctr->Fill(timeStamp[chID]);
          // ctrVSw->Fill(FloodZ,timeStamp[chID]);

          float weight = 0;
          float corrected_t = 0;


          float ref_t = line_central->Eval(m_w);
          float t_shift = ref_t - line_central->Eval(FloodZ);


          weight += 1.0/(pow(ctr->GetRMS(),2));

          corrected_t += timeStamp[chID] * 1.0/(pow(ctr->GetRMS(),2));

          for(int iN = 0; iN < 8 ; iN++)
          {
            float corr_time = line_n[iN]->Eval(FloodZ);
            weight += 1.0/(pow(ctr_n[iN]->GetRMS(),2));
            corrected_t += (timeStamp[n_ids[iN]] - corr_time) * 1.0/(pow(ctr_n[iN]->GetRMS(),2));
          }
          corrected_t = corrected_t / weight;
          ctr_avg->Fill(corrected_t);
          ctr_doi->Fill(corrected_t + t_shift);
        }
      }

    }







    // LOOP ON ALL ACCEPTED CRYSTALS
    // I.E. ALL CRYSTALS IN THE CALIBRATION FILES
    // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    // {
    //   // do whatever you want
    //   // in particular, some tips:
    //   // 1) if(crystal[iCry].FormulaTagAnalysis->EvalInstance()){} --> this will select events in the photopeak of the external reference crystal, if there is one. Notice that each crystal can have its own FormulaTagAnalysis, since they can come from different calibration runs
    //   // 2) if(crystal[iCry].FormulaAnalysis->EvalInstance()) --> this will select events in a crystal (geometrical position in u-v-w and photopeak condition)
    //   // 3) to calculate w, if needed
    //   // float FloodZ = calculateFloodZ(charge,crystal[iCry]);
    //   // notice that the detectorSaturation values are written directly into the crystal struct, because they can be acquired in different conditions (but it has to be the same photodetector!)
    //
    //   //example here, selecting photopeak events in ext ref, photopeak events the crystal(s) found in calibration file(s), and calculating FloodZ for each
    //   if(crystal[iCry].accepted)
    //   {
    //     if(crystal[iCry].FormulaTagAnalysis->EvalInstance())
    //     {
    //       if(crystal[iCry].FormulaAnalysis->EvalInstance())
    //       {
    //         goodEventsAnalysis++;
    //         float FloodZ = calculateFloodZ(charge,crystal[iCry]);
    //         // std::cout << FloodZ << std::endl;
    //       }
    //     }
    //   }
    // }

    //LOOP COUNTER
    counterAnalysis++;
    int perc = ((100*counterAnalysis)/neventAnalysis);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }

  TF1 *g1 = new TF1("g1","gaus",-2e-9,2e-9);
  TF1 *g2 = new TF1("g2","gaus",-2e-9,2e-9);
  TF1 *g3 = new TF1("g3","gaus",-2e-9,2e-9);

  ctr->Fit(g1);
  ctr_avg->Fit(g2);
  ctr_doi->Fit(g3);

  std::cout << g1->GetParameter(2)*2.355*sqrt(2) << " "
            << g2->GetParameter(2)*2.355*sqrt(2) << " "
            << g3->GetParameter(2)*2.355*sqrt(2) << " "
            << std::endl;

  //
  // // sort crystals struct (can be useful)
  // std::sort(crystal.begin(), crystal.end(), compareByNumber);

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  // write whatever you want to save
  spectrum->Write();
  scatter->Write();
  scatter3D->Write();
  scatterReal->Write();


  spectrum1Cry   ->Write();
  scatter1Cry    ->Write();
  scatter3D1Cry  ->Write();
  scatterReal1Cry->Write();

  H_realZ->Write();
  H_w->Write();
  pfx->Write();
  c_gr->Write();
  Rz_w->Write();
  ctr->Write();
  ctr_avg->Write();
  ctr_doi->Write();
  ctrVSw->Write();
  pfx_central->Write();
  for(int iN = 0; iN < 8 ; iN++)
  {
    ctr_n[iN]->Write();
    ctrVSw_n[iN]->Write();
    profile_n[iN]->Write();
  }

  //save the full command line
  TNamed CommandNameD("Command",streamCommand.str().c_str());
  CommandNameD.Write();
  outputFile->Close();
  std::cout << "Results saved in file " << outputFileName << std::endl;

  return 0;
}
// end of main program




// feedback to user
void usage()
{
  std::cout << "\t"
            << "[-f|--folder] <folder> "
            << "\t"
            << "[-p|--prefix] <prefix> "
            << "\t"
            << "[-o|--output] <output> "
            << std::endl
            << "\t\t"
            << "<folder>         - path to folder were input files are located - default = ./ "
            << std::endl
            << "\t\t"
            << "<prefix>         - prefix of input TTree files "
            << std::endl
            << "\t\t"
            << "<output>         - output file name - default = outputBareboneFile.root "
            << std::endl
            << std::endl;
}
