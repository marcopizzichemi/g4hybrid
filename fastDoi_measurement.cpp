// compile with
// g++ -o ../build/fastDoi_measurement fastDoi_measurement.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas


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


#include "/home/marco/cernbox/Universita/NewClearPEM/Programs/ModuleCalibration/code/libraries/Extract.h"

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


// struct neighbour
// {
//   int channel;
//
// }


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
  int jumper = 2;
  float hardCut = 1; // cut of weird w part. = 1 mean no cut
  int useLine = 0;

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
    { "jumper", required_argument, 0, 0 },
    { "hardCut", required_argument, 0, 0 },
    { "useLine", no_argument, 0, 0 },
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
      jumper = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      hardCut = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      useLine = 1;
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
  // tree->SetBranchAddress("RealX", &RealX, &bRealX);
  // tree->SetBranchAddress("RealY", &RealY, &bRealY);
  // tree->SetBranchAddress("RealZ", &RealZ, &bRealZ);

  // std::cout << "Entries = " << tree->GetEntries() << std::endl;

  float xmppc[16] = {-1.6,-1.6,-4.8,-4.8,4.8 ,4.8 ,1.6 ,1.6 ,-4.8,-4.8,-1.6,-1.6,1.6 ,1.6 ,4.8 ,4.8};
  float ymppc[16] = {1.6, 4.8, 4.8, 1.6,1.6 ,4.8 ,4.8 ,1.6 ,-1.6,-4.8,-4.8,-1.6,-1.6,-4.8,-4.8,-1.6};

  // int channels[16] = {16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}; //j1
  int channels[16];
  int neighboursID[8] ;
  int hybridID;

  int channels_j1[16] = {16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}; //j1
  int channels_j2[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; //j2

  int neighboursID_j1[8] = {16,23,20,27,31,26,20,30}; //j1
  int neighboursID_j2[8] = {0,7,4,11,15,10,13,14}; //j2

  if(jumper == 1)
  {
    hybridID = 28;
    for(int i = 0 ; i < 16 ; i++)
    {
      channels[i] = channels_j1[i];
    }
    for(int i = 0 ; i < 8 ; i++)
    {
      neighboursID[i] = neighboursID_j1[i];
    }
  }
  else
  {
    if(jumper == 2)
    {
      hybridID = 12;
      for(int i = 0 ; i < 16 ; i++)
      {
        channels[i] = channels_j2[i];
      }
      for(int i = 0 ; i < 8 ; i++)
      {
        neighboursID[i] = neighboursID_j2[i];
      }
    }
    else
    {
      std::cout << "ERROR, invalid jumper!!! " << std::endl;
      return 1;
    }
  }

  float enLow = 14000;
  float enUp = 28000;



  // int hybridID = 28; //j1
  // int hybridID = 12; //j2


  // int neighboursID[8] = {16,23,20,27,31,26,20,30}; //j1
  // int neighboursID[8] = {0,7,4,11,15,10,13,14}; //j1

  TH1F *delay_n[8];
  delay_n[0] = new TH1F("delay_n0","delay_n0",25,-2e-9,10e-9);
  delay_n[1] = new TH1F("delay_n1","delay_n1",25,-2e-9,10e-9);
  delay_n[2] = new TH1F("delay_n2","delay_n2",25,-2e-9,10e-9);
  delay_n[3] = new TH1F("delay_n3","delay_n3",25,-2e-9,10e-9);
  delay_n[4] = new TH1F("delay_n4","delay_n4",25,-2e-9,10e-9);
  delay_n[5] = new TH1F("delay_n5","delay_n5",25,-2e-9,10e-9);
  delay_n[6] = new TH1F("delay_n6","delay_n6",25,-2e-9,10e-9);
  delay_n[7] = new TH1F("delay_n7","delay_n7",25,-2e-9,10e-9);

  TH2F *delayVSw_n[8];
  delayVSw_n[0] = new TH2F("delayVSw_n0","delayVSw_n0",100,0,1,25,-2e-9,10e-9);
  delayVSw_n[1] = new TH2F("delayVSw_n1","delayVSw_n1",100,0,1,25,-2e-9,10e-9);
  delayVSw_n[2] = new TH2F("delayVSw_n2","delayVSw_n2",100,0,1,25,-2e-9,10e-9);
  delayVSw_n[3] = new TH2F("delayVSw_n3","delayVSw_n3",100,0,1,25,-2e-9,10e-9);
  delayVSw_n[4] = new TH2F("delayVSw_n4","delayVSw_n4",100,0,1,25,-2e-9,10e-9);
  delayVSw_n[5] = new TH2F("delayVSw_n5","delayVSw_n5",100,0,1,25,-2e-9,10e-9);
  delayVSw_n[6] = new TH2F("delayVSw_n6","delayVSw_n6",100,0,1,25,-2e-9,10e-9);
  delayVSw_n[7] = new TH2F("delayVSw_n7","delayVSw_n7",100,0,1,25,-2e-9,10e-9);

  TH2F *delayVSw_corr_n[8];
  delayVSw_corr_n[0] = new TH2F("delayVSw_corr_n0","delayVSw_corr_n0",100,0,1,25,-2e-9,10e-9);
  delayVSw_corr_n[1] = new TH2F("delayVSw_corr_n1","delayVSw_corr_n1",100,0,1,25,-2e-9,10e-9);
  delayVSw_corr_n[2] = new TH2F("delayVSw_corr_n2","delayVSw_corr_n2",100,0,1,25,-2e-9,10e-9);
  delayVSw_corr_n[3] = new TH2F("delayVSw_corr_n3","delayVSw_corr_n3",100,0,1,25,-2e-9,10e-9);
  delayVSw_corr_n[4] = new TH2F("delayVSw_corr_n4","delayVSw_corr_n4",100,0,1,25,-2e-9,10e-9);
  delayVSw_corr_n[5] = new TH2F("delayVSw_corr_n5","delayVSw_corr_n5",100,0,1,25,-2e-9,10e-9);
  delayVSw_corr_n[6] = new TH2F("delayVSw_corr_n6","delayVSw_corr_n6",100,0,1,25,-2e-9,10e-9);
  delayVSw_corr_n[7] = new TH2F("delayVSw_corr_n7","delayVSw_corr_n7",100,0,1,25,-2e-9,10e-9);


  TProfile *pfx_n[8];
  TGraph *graph_n[8];
  TGraph *g_sigma_n[8];
  TF1 *lineReg_n[8];
  lineReg_n[0] = new TF1("lineReg_n0","pol1",0,1);
  lineReg_n[1] = new TF1("lineReg_n1","pol1",0,1);
  lineReg_n[2] = new TF1("lineReg_n2","pol1",0,1);
  lineReg_n[3] = new TF1("lineReg_n3","pol1",0,1);
  lineReg_n[4] = new TF1("lineReg_n4","pol1",0,1);
  lineReg_n[5] = new TF1("lineReg_n5","pol1",0,1);
  lineReg_n[6] = new TF1("lineReg_n6","pol1",0,1);
  lineReg_n[7] = new TF1("lineReg_n7","pol1",0,1);

  TH2F *ctrVSw = new TH2F("ctrVSw","ctrVSw",100,0,1,25,-2e-9,2e-9);
  TProfile *pfx;
  TF1 *lineReg =  new TF1("lineReg","pol1",0,1);;
  TGraph *graph;
  TGraph *g_sigma;


  // int hybridID = 28;


  float saturation[16] =	{1.088270E-08	,	1.054790E-08	,	1.101300E-08	,	1.028710E-08	,	1.072330E-08	,	1.100130E-08	,	1.074750E-08	,	1.075220E-08	,	1.036340E-08	,	1.138910E-08	,	1.077850E-08	,	1.163180E-08	,	1.129270E-08,1.120150E-08	,	1.084750E-08	,	1.090670E-08}; //saturation parameters with glass

  float chargeBinning = 320e-15;
  //
  //
  //( -sat * TMath::Log(1.0 - (  ch_dig   /sat) ) )


  // digitizer     = 16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
  // timing = 16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
  // mppc          = C2,D2   ,D1,C1,C4,D4,D3,C3,B1,A1,A2,B2,B3,A3,A4,B4

  // plotPositions = 6   ,2   ,1   ,5   ,8   ,4   ,3   ,7   ,9   ,13  ,14  ,10  ,11  ,15  ,16  ,12
  // xPositions    = -1.6,-1.6,-4.8,-4.8,4.8 ,4.8 ,1.6 ,1.6 ,-4.8,-4.8,-1.6,-1.6,1.6 ,1.6 ,4.8 ,4.8
  // yPositions    = 1.6, 4.8, 4.8, 1.6,1.6 ,4.8 ,4.8 ,1.6 ,-1.6,-4.8,-4.8,-1.6,-1.6,-4.8,-4.8,-1.6


  TH2F *scatter = new TH2F("Flood Map","Flood Map",500,-7,7,500,-7,7);
  TH2F *scatterReal = new TH2F("Flood Map Real","Flood Map Real",500,-7,7,500,-7,7);
  TH2F *Rz_w = new TH2F("Rz_w","Rz_w",500,0,1,100,-8,8);
  TH1F *spectrum = new TH1F("spectrum","spectrum",250,0,100000);
  TH1F *H_realZ = new TH1F("H_realZ","H_realZ",100,-8,8);
  TH1F *H_w = new TH1F("H_w","H_w",100,0,1);

  TH3I *scatter3D = new TH3I("Flood Map 3D","Flood Map 3D",100,-7,7,100,-7,7,100,0,1);

  TH1F *spectrum1Cry = new TH1F("spectrum1Cry","spectrum 1Cry",250,0,100000);
  TH1F *spectrum1Cry_a = new TH1F("spectrum1Cry_a","spectrum 1Cry_a",250,0,100000);
  TH1F *spectrum1Cry_b = new TH1F("spectrum1Cry_b","spectrum 1Cry_b",250,0,100000);
  TH2F *scatter1Cry = new TH2F("Flood Map1Cry","Flood Map 1Cry",500,-7,7,500,-7,7);
  TH2F *scatterReal1Cry = new TH2F("Flood Map Real 1Cry","Flood Map Real 1Cry",500,-7,7,500,-7,7);
  TH3I *scatter3D1Cry = new TH3I("Flood Map 3D 1Cry","Flood Map 3D 1Cry",100,-7,7,100,-7,7,100,0,1);

  TH1F *refCry = new TH1F("refCry","refCry",1000,0,60000);

  TH1F *ctr = new TH1F("ctr","ctr",100,-2e-9,2e-9);
  TH1F *ctr_9 = new TH1F("ctr_9","ctr_9",100,-2e-9,2e-9);
  TH1F *ctr_corrected = new TH1F("ctr_corrected","ctr_corrected",100,-2e-9,2e-9);

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
  // LOOP for tag
  for (long long int i=0;i<neventAnalysis;i++)
  {

    tree->GetEvent(i);              //read complete accepted event in memory
    // float totalCharge = 0;
    refCry->Fill(charge[62]);

    //LOOP COUNTER
    counterAnalysis++;
    int perc = ((100*counterAnalysis)/neventAnalysis);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }

  }
  std::cout << std::endl;

  float refMin = 41000;
  float refMax = 45000;


  // neighbour map

  // std::vector<std::vector<int>> map_n;
  //
  // // for ch 0
  // std::vector<int> el0;
  // el0.push_back(1);
  // el0.push_back(2);
  // el0.push_back(6);
  // el0.push_back(3);
  // el0.push_back(7);
  // el0.push_back(8);
  // el0.push_back(11);
  // el0.push_back(12);
  // map_n.push_back(el0);
  //
  // // for ch 1
  // std::vector<int> el1;
  // el1.push_back(2);
  // el1.push_back(6);
  // el1.push_back(3);
  // el1.push_back(0);
  // el1.push_back(7);
  // map_n.push_back(el1);
  //
  // // for ch 2
  // std::vector<int> el2;
  // el2.push_back(1);
  // el2.push_back(0);
  // el2.push_back(3);
  // map_n.push_back(el2);
  //
  // // for ch 3
  // std::vector<int> el3;
  // el3.push_back(2);
  // el3.push_back(1);
  // el3.push_back(0);
  // el3.push_back(8);
  // el3.push_back(11);
  // map_n.push_back(el3);
  //
  // // for ch 4
  // std::vector<int> el4;
  // el4.push_back(6);
  // el4.push_back(5);
  // el4.push_back(7);
  // el4.push_back(12);
  // el4.push_back(15);
  // map_n.push_back(el4);
  //
  // // for ch 5
  // std::vector<int> el5;
  // el5.push_back(6);
  // el5.push_back(7);
  // el5.push_back(4);
  // map_n.push_back(el5);
  //
  // // for ch 6
  // std::vector<int> el6;
  // el6.push_back(1);
  // el6.push_back(5);
  // el6.push_back(0);
  // el6.push_back(7);
  // el6.push_back(4);
  // map_n.push_back(el6);
  //
  // // for ch 7
  // std::vector<int> el7;
  // el7.push_back(1);
  // el7.push_back(6);
  // el7.push_back(5);
  // el7.push_back(0);
  // el7.push_back(4);
  // el7.push_back(11);
  // el7.push_back(12);
  // el7.push_back(15);
  // map_n.push_back(el7);
  //
  // // for ch 8
  // std::vector<int> el8;
  // el8.push_back(3);
  // el8.push_back(0);
  // el8.push_back(11);
  // el8.push_back(9);
  // el8.push_back(10);
  // map_n.push_back(el8);
  //
  // // for ch 9
  // std::vector<int> el9;
  // el9.push_back(8);
  // el9.push_back(11);
  // el9.push_back(10);
  // map_n.push_back(el9);
  //
  // // for ch 10
  // std::vector<int> el10;
  // el10.push_back(8);
  // el10.push_back(11);
  // el10.push_back(12);
  // el10.push_back(9);
  // el10.push_back(13);
  // map_n.push_back(el10);
  //
  // // for ch 11
  // std::vector<int> el11;
  // el11.push_back(3);
  // el11.push_back(0);
  // el11.push_back(7);
  // el11.push_back(8);
  // el11.push_back(12);
  // el11.push_back(9);
  // el11.push_back(10);
  // el11.push_back(13);
  // map_n.push_back(el11);
  //
  // // for ch 12
  // std::vector<int> el12;
  // el12.push_back(0);
  // el12.push_back(7);
  // el12.push_back(4);
  // el12.push_back(11);
  // el12.push_back(15);
  // el12.push_back(10);
  // el12.push_back(13);
  // el12.push_back(14);
  // map_n.push_back(el12);
  //
  // // for ch 13
  // std::vector<int> el13;
  // el13.push_back(10);
  // el13.push_back(11);
  // el13.push_back(12);
  // el13.push_back(15);
  // el13.push_back(14);
  // map_n.push_back(el13);
  //
  // // for ch 10
  // std::vector<int> el14;
  // el14.push_back(12);
  // el14.push_back(13);
  // el14.push_back(15);
  // map_n.push_back(el14);
  //
  // // for ch 10
  // std::vector<int> el15;
  // el15.push_back(7);
  // el15.push_back(4);
  // el15.push_back(12);
  // el15.push_back(13);
  // el15.push_back(14);
  // map_n.push_back(el15);


  // for(unsigned int iCh = 0 ; iCh < map_n.size(); iCh++)
  // {
  //   std::cout << iCh << ": ";
  //   for(unsigned int iN = 0 ; iN < map_n[iCh].size(); iN++)
  //   {
  //     std::cout << map_n[iCh][iN] << " ";
  //   }
  //   std::cout << std::endl;
  // }



  // MAIN LOOP
  goodEventsAnalysis = 0;
  counterAnalysis = 0;
  for (long long int i=0;i<neventAnalysis;i++)
  {

    tree->GetEvent(i);              //read complete accepted event in memory


    //only photopeak of tag
    if(charge[62] > refMin && charge[62] < refMax)
    {

      // correct sipm by saturation

      for(int iCh = 0 ; iCh < 16 ; iCh++)
      {
        //( -sat * TMath::Log(1.0 - (  charge   /sat) ) )
        int ch = channels[iCh];
        charge[ch] = ( -(saturation[iCh]/chargeBinning) * TMath::Log(1.0 - ( charge[ch]/(saturation[iCh]/chargeBinning) ) )  );
      }


      float totalCharge = 0;
      // float totalChargeFlood = 0;
      //find max charge
      // and calc sum
      float maxCharge = 0;
      int IDch = -1;
      // int relevant_id[9] = {16,23,20,27,28,31,26,29,30};

      // int relevant_id[9] = {0,7,,};
      // float rel_xmppc[9] = {-1.6,1.6,4.8,-1.6,1.6,4.8,-1.6,1.6,4.8};
      // float rel_ymppc[9] = {1.6,1.6,1.6,-1.6,-1.6,-1.6,-4.8,-4.8,-4.8};

      // find max charge id
      for(int iCh = 0 ; iCh < 16; iCh++)
      {
        int ch = channels[iCh];
        if(charge[ch] > maxCharge)
        {
          maxCharge = charge[ch];
          IDch = ch;
        }
        // totalCharge += charge[channels[iCh]];
      }


      float FloodX = 0;
      float FloodY = 0;
      float FloodZ = 0;

      for(unsigned int iCh = 0 ; iCh <16; iCh++)
      {
        // if(iCh == IDch) // find the "trigger"
        // {
        //   // maxCharge = charge[ channels[iCh] ];
        // }
        int ch = channels[iCh];
        totalCharge += charge[ ch ];

        FloodX += charge[ch]*xmppc[iCh];
        FloodY += charge[ch]*ymppc[iCh];
      }

      // std::cout << totalCharge << std::endl;

      if(totalCharge > 2000)
      {

        spectrum->Fill(totalCharge);

        FloodX = FloodX/totalCharge;
        FloodY = FloodY/totalCharge;
        FloodZ = maxCharge/totalCharge;

        scatter->Fill(FloodX,FloodY);
        scatter3D->Fill(FloodX,FloodY,FloodZ);
        if(IDch == hybridID)
          // if(true)
        {
          scatter1Cry->Fill(FloodX,FloodY);
          scatter3D1Cry->Fill(FloodX,FloodY,FloodZ);
          spectrum1Cry->Fill(totalCharge);
        }


        if( (FloodX > fxMin ) && (FloodX < fxMax ) && (FloodY > fyMin ) && (FloodY < fyMax )  )
        {
          if(IDch == hybridID)
          // if(true)
          {

            if(FloodZ > hardCut)
            {
              spectrum1Cry_b->Fill(totalCharge);
            }
            else
            {
              spectrum1Cry_a->Fill(totalCharge);
            }

            // select photopeak
            if( (totalCharge > enLow ) && (totalCharge < enUp ) )
            {
              //FIXME hardcut of "weird" w

              if(FloodZ < hardCut)
              {
                goodEventsAnalysis++;
                // H_realZ->Fill(RealZ);
                H_w->Fill(FloodZ);
                // Rz_w->Fill(FloodZ,RealZ);
                if(timeStamp[62] != 0 && timeStamp[hybridID] != 0)
                {
                  if(timeStamp[62] > -80e-9 && timeStamp[62] < -70e-9  )
                  {
                    if(timeStamp[hybridID] > -80e-9 && timeStamp[hybridID] < -70e-9  )
                    {
                      ctr->Fill(timeStamp[62] - timeStamp[hybridID]);
                      ctrVSw->Fill(FloodZ,timeStamp[62] - timeStamp[hybridID]);
                      // float t2 = 0;
                      // t2 += timeStamp[hybridID];
                      for(int iT = 0; iT < 8 ; iT++)
                      {
                        // t2 +=timeStamp[neighboursID[iT]];
                        delay_n[iT]->Fill(timeStamp[neighboursID[iT]] - timeStamp[hybridID]);
                        delayVSw_n[iT]->Fill(FloodZ,timeStamp[neighboursID[iT]] - timeStamp[hybridID]);
                      }
                      // t2 = t2 / 9.0;
                      // ctr_all->Fill(timeStamp[62] - t2);

                    }
                  }
                }
              }

            }
          }
        }
      }
    }


    //LOOP COUNTER
    counterAnalysis++;
    int perc = ((100*counterAnalysis)/neventAnalysis);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }



  std::cout << std::endl;

  pfx = ctrVSw->ProfileX("pfx");
  pfx->Fit(lineReg);


  std::vector <float>  cx_point;
  std::vector <float>  cy_point;
  std::vector <float>  cx_sigma;
  std::vector <float>  cy_sigma;

  int cpoints = 0;

  for(int bin = 0 ; bin < pfx->GetNbinsX() ; bin++)
  {
    if(pfx->GetBinContent(bin) != 0)
    {
      cpoints++;
      cx_point.push_back(pfx->GetBinCenter(bin));
      cy_point.push_back(pfx->GetBinContent(bin));
      cx_sigma.push_back(pfx->GetBinCenter(bin));
      cy_sigma.push_back(pfx->GetBinError(bin));
    }
  }
  std::stringstream csname;
  graph = new TGraph(cpoints,&cx_point[0],&cy_point[0]);
  csname.str("");
  csname << "Graph";
  graph->SetName(csname.str().c_str());
  graph->SetTitle(csname.str().c_str());
  graph->GetXaxis()->SetTitle("w");
  graph->GetYaxis()->SetTitle("delta T [s]");
  csname.str("");

  g_sigma = new TGraph(cpoints,&cx_sigma[0],&cy_sigma[0]);
  csname.str("");
  csname << "Graph_Sigma";
  g_sigma->SetName(csname.str().c_str());
  g_sigma->SetTitle(csname.str().c_str());
  g_sigma->GetXaxis()->SetTitle("w");
  g_sigma->GetYaxis()->SetTitle("sigma delta T [s]");
  csname.str("");



  for(int iT = 0; iT < 8 ; iT++)
  {
    Double_t res[4];
    // t2 +=timeStamp[neighboursID[iT]];
    fitWithEMG(delay_n[iT],res);
    //res[0] = mean
    //res[0] = sigma
    //res[0] = mean_err
    //res[0] = sigma_err
    std::cout << res[0] << " "
    << res[1] << " "
    << res[2] << " "
    << res[3] << " ";
    std::stringstream sname;
    sname << "Profile" << iT;
    pfx_n[iT] = delayVSw_n[iT]->ProfileX(sname.str().c_str());

    std::vector <float> x_point;
    std::vector <float> y_point;
    // std::vector <float> ex_point;
    // std::vector <float> ey_point;

    std::vector <float>  x_sigma;
    std::vector <float>  y_sigma;
    // std::vector <float> ex_sigma;
    // std::vector <float> ey_sigma;

    int points = 0;
    for(int bin = 0 ; bin < pfx_n[iT]->GetNbinsX() ; bin++)
    {

      if(pfx_n[iT]->GetBinContent(bin) != 0)
      {
        points++;
        x_point.push_back(pfx_n[iT]->GetBinCenter(bin));
        y_point.push_back(pfx_n[iT]->GetBinContent(bin));
        x_sigma.push_back(pfx_n[iT]->GetBinCenter(bin));
        y_sigma.push_back(pfx_n[iT]->GetBinError(bin));
      }
    }

    graph_n[iT] = new TGraph(points,&x_point[0],&y_point[0]);
    sname.str("");
    sname << "Graph_n" << iT;
    graph_n[iT]->SetName(sname.str().c_str());
    graph_n[iT]->SetTitle(sname.str().c_str());
    graph_n[iT]->GetXaxis()->SetTitle("w");
    graph_n[iT]->GetYaxis()->SetTitle("delay [s]");
    sname.str("");

    g_sigma_n[iT] = new TGraph(points,&x_sigma[0],&y_sigma[0]);
    sname.str("");
    sname << "Graph_Sigma_n" << iT;
    g_sigma_n[iT]->SetName(sname.str().c_str());
    g_sigma_n[iT]->SetTitle(sname.str().c_str());
    g_sigma_n[iT]->GetXaxis()->SetTitle("w");
    g_sigma_n[iT]->GetYaxis()->SetTitle("sigma [s]");
    sname.str("");




    pfx_n[iT]->Fit(lineReg_n[iT]);

  }
  std::cout << std::endl;


  goodEventsAnalysis = 0;
  counterAnalysis = 0;
  for (long long int i=0;i<neventAnalysis;i++)
  {

    tree->GetEvent(i);              //read complete accepted event in memory


    //only photopeak of tag
    if(charge[62] > refMin && charge[62] < refMax)
    {

      // correct sipm by saturation

      for(int iCh = 0 ; iCh < 16 ; iCh++)
      {
        //( -sat * TMath::Log(1.0 - (  charge   /sat) ) )
        int ch = channels[iCh];
        charge[ch] = ( -(saturation[iCh]/chargeBinning) * TMath::Log(1.0 - ( charge[ch]/(saturation[iCh]/chargeBinning) ) )  );
      }


      float totalCharge = 0;
      // float totalChargeFlood = 0;
      //find max charge
      // and calc sum
      float maxCharge = 0;
      int IDch = -1;
      // int relevant_id[9] = {16,23,20,27,28,31,26,29,30};

      // int relevant_id[9] = {0,7,,};
      // float rel_xmppc[9] = {-1.6,1.6,4.8,-1.6,1.6,4.8,-1.6,1.6,4.8};
      // float rel_ymppc[9] = {1.6,1.6,1.6,-1.6,-1.6,-1.6,-4.8,-4.8,-4.8};

      // find max charge id
      for(int iCh = 0 ; iCh < 16; iCh++)
      {
        int ch = channels[iCh];
        if(charge[ch] > maxCharge)
        {
          maxCharge = charge[ch];
          IDch = ch;
        }
        // totalCharge += charge[channels[iCh]];
      }


      float FloodX = 0;
      float FloodY = 0;
      float FloodZ = 0;

      for(unsigned int iCh = 0 ; iCh <16; iCh++)
      {
        // if(iCh == IDch) // find the "trigger"
        // {
        //   // maxCharge = charge[ channels[iCh] ];
        // }
        int ch = channels[iCh];
        totalCharge += charge[ ch ];

        FloodX += charge[ch]*xmppc[iCh];
        FloodY += charge[ch]*ymppc[iCh];
      }

      // std::cout << totalCharge << std::endl;

      if(totalCharge > 2000)
      {

        // spectrum->Fill(totalCharge);

        FloodX = FloodX/totalCharge;
        FloodY = FloodY/totalCharge;
        FloodZ = maxCharge/totalCharge;

        // scatter->Fill(FloodX,FloodY);
        // scatter3D->Fill(FloodX,FloodY,FloodZ);


        if( (FloodX > fxMin ) && (FloodX < fxMax ) && (FloodY > fyMin ) && (FloodY < fyMax )  )
        {
          if(IDch == hybridID)
          // if(true)
          {
            // scatter1Cry->Fill(FloodX,FloodY);
            // scatter3D1Cry->Fill(FloodX,FloodY,FloodZ);
            // spectrum1Cry->Fill(totalCharge);
            // if(FloodZ > hardCut)
            // {
            //   spectrum1Cry_b->Fill(totalCharge);
            // }
            // else
            // {
            //   spectrum1Cry_a->Fill(totalCharge);
            // }

            // select photopeak
            if( (totalCharge > enLow ) && (totalCharge < enUp ) )
            {
              //FIXME hardcut of "weird" w

              if(FloodZ < hardCut)
              {
                goodEventsAnalysis++;
                // H_realZ->Fill(RealZ);
                // H_w->Fill(FloodZ);
                // Rz_w->Fill(FloodZ,RealZ);

                if(timeStamp[62] != 0 && timeStamp[hybridID] != 0)
                {
                  if(timeStamp[62] > -80e-9 && timeStamp[62] < -70e-9  )
                  {
                    if(timeStamp[hybridID] > -80e-9 && timeStamp[hybridID] < -70e-9  )
                    {
                      // ctr->Fill(timeStamp[62] - timeStamp[hybridID]);
                      float t2 = 0.0;
                      float TOTweight = 0.0;

                      float weight0;
                      if(useLine)
                      {
                        weight0 = 1;
                      }
                      else
                      {
                        weight0 = (1.0/pow(g_sigma->Eval(FloodZ),2));
                      }

                      // float weight0 = (1.0/pow(g_sigma->Eval(FloodZ),2));
                      // float weight0 = 1;

                      t2 += timeStamp[hybridID] * weight0 ;
                      TOTweight += weight0;

                      for(int iT = 0; iT < 8 ; iT++)
                      {
                        float delay;
                        float weight;
                        if(useLine)
                        {
                          delay = lineReg_n[iT]->Eval(FloodZ);
                          weight = 1;
                        }
                        else
                        {
                          delay = graph_n[iT]->Eval(FloodZ);
                          weight = (1.0/pow(g_sigma_n[iT]->Eval(FloodZ),2));
                        }
                        // float delay = graph_n[iT]->Eval(FloodZ);
                        // float weight = (1.0/pow(g_sigma_n[iT]->Eval(FloodZ),2));

                        // float delay = lineReg_n[iT]->Eval(FloodZ);
                        // float weight = 1;


                        float t_corrected = timeStamp[neighboursID[iT]] - timeStamp[hybridID] - delay;

                        delayVSw_corr_n[iT]->Fill(FloodZ,t_corrected);

                        t2 += t_corrected * weight ;
                        TOTweight += weight;
                        // delay_n[iT]->Fill(timeStamp[neighboursID[iT]] - timeStamp[hybridID]);
                        // delayVSw_n[iT]->Fill(FloodZ,timeStamp[neighboursID[iT]] - timeStamp[hybridID]);
                      }
                      t2 = t2 / TOTweight;
                      ctr_9->Fill(timeStamp[62] - t2);
                      // calc a reference point
                      // to w = 0.4 in this example
                      // float ref_t = graph->Eval(0.4);
                      // float distance = ref_t - graph->Eval(FloodZ);
                      float ref_t = lineReg->Eval(0.4);
                      float distance = ref_t - lineReg->Eval(FloodZ);


                      ctr_corrected->Fill(timeStamp[62] - (t2 - distance));


                    }
                  }
                }
              }
            }
          }
        }
      }
    }


    //LOOP COUNTER
    counterAnalysis++;
    int perc = ((100*counterAnalysis)/neventAnalysis);
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }


  // timing correction




  std::cout << "Good events = " << goodEventsAnalysis << std::endl;









      //
      // // sort crystals struct (can be useful)
      // std::sort(crystal.begin(), crystal.end(), compareByNumber);

      TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
      outputFile->cd();
      refCry->Write();
      // write whatever you want to save
      spectrum->Write();
      scatter->Write();
      scatter3D->Write();
      scatterReal->Write();


      spectrum1Cry   ->Write();
      spectrum1Cry_a   ->Write();
      spectrum1Cry_b   ->Write();
      scatter1Cry    ->Write();
      scatter3D1Cry  ->Write();
      scatterReal1Cry->Write();

      H_realZ->Write();
      H_w->Write();

      ctr->Write();
      ctr_9->Write();
      ctr_corrected->Write();

      ctrVSw ->Write();
      pfx ->Write();
      graph ->Write();
      g_sigma->Write();

      for(int iT = 0; iT < 8 ; iT++)
      {
        // t2 +=timeStamp[neighboursID[iT]];
        delay_n[iT]->Write();
      }
      for(int iT = 0; iT < 8 ; iT++)
      {
        // t2 +=timeStamp[neighboursID[iT]];
        delayVSw_n[iT]->Write();
      }
      for(int iT = 0; iT < 8 ; iT++)
      {
        // t2 +=timeStamp[neighboursID[iT]];
        pfx_n[iT]->Write();
        graph_n[iT]->Write();
        g_sigma_n[iT]->Write();
      }
      for(int iT = 0; iT < 8 ; iT++)
      {
        delayVSw_corr_n[iT]->Write();
      }

      Rz_w->Write();

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
