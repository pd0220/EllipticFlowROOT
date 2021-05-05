// ROOT C++ code for ellipTic flow (v2) measurement

// including used libraries
#define particle_tree_cxx
#include "particle_tree.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TGraphErrors.h>

// ------------------------------------------------------------------------------------------------------------------------------

// main function
int main(int argc, const char **argv)
{
  // centrality bounds
  int NC = 2;
  std::pair<int, int> centrality1 = {0, 30};
  std::pair<int, int> centrality2 = {40, 70};
  std::vector<std::pair<int, int>> centralities = {centrality1, centrality2};
  // transverse momentum splitting
  int NpT = (int)(2. / 0.1) - 1;
  std::vector<double> pTSplitting(static_cast<size_t>(NpT));
  double pTInit = 0.1;
  double pTDelta = 0.1;
  std::generate(pTSplitting.begin(), pTSplitting.end(), [&]() { return pTInit += pTDelta; });

  // ------------------------------------------------------------------------------------------------------------------------------

  // contiainers for elliptic flow (v2) results
  std::vector<TGraph *> v2Result(NC);
  std::vector<TGraphErrors *> v2Error(NC);

  // ------------------------------------------------------------------------------------------------------------------------------

  // reading reaction plane resolution values from outer file
  // start reading
  std::ifstream fileToRead;
  std::string fileName = "reactionplane_resolution_rxnp.txt";
  fileToRead.open(fileName);

  // determine number of columns
  std::string firstLine;
  std::getline(fileToRead, firstLine);
  std::stringstream firstLineStream(firstLine);

  // number of columns in given file
  int numOfCols = 0;
  std::string tmpString;
  // count number of writes to a temporary string container
  while (firstLineStream >> tmpString)
    numOfCols++;

  fileToRead.close();

  // string for all the lines
  std::string line;

  // data structures
  std::vector<double> RPCentrality;
  std::vector<double> RPResolution;

  // reopen file
  fileToRead.open(fileName);
  // check if open
  if (fileToRead.is_open())
  {
    // read line by line
    int i = 0;
    while (std::getline(fileToRead, line))
    {
      // using stringstream to write data
      std::stringstream dataStream(line);
      double centralityTemp = 0.;
      double resolutionTemp = 0.;
      dataStream >> centralityTemp;
      dataStream >> resolutionTemp;
      RPCentrality.push_back(centralityTemp);
      RPResolution.push_back(resolutionTemp);
      i++;
    }
    // close file
    fileToRead.close();
  }
  // error check
  else
  {
    std::cout << "ERROR\nProblem occured while reading given file." << std::endl;
    std::exit(-1);
  }

  // ------------------------------------------------------------------------------------------------------------------------------

  // reaction plane resolution means
  // split resolution containers
  std::vector<double> RPResolution1, RPResolution2;
  for (int iRes = 0; iRes < static_cast<int>(RPResolution.size()); iRes++)
  {
    if (RPCentrality[iRes] >= centrality1.first && RPCentrality[iRes] <= centrality1.second)
      RPResolution1.push_back(RPResolution[iRes]);
    else if (RPCentrality[iRes] >= centrality2.first && RPCentrality[iRes] <= centrality2.second)
      RPResolution2.push_back(RPResolution[iRes]);
  }
  // means
  double RPMean1 = std::accumulate(RPResolution1.begin(), RPResolution1.end(), 0.) / static_cast<int>(RPResolution1.size());
  double RPMean2 = std::accumulate(RPResolution2.begin(), RPResolution2.end(), 0.) / static_cast<int>(RPResolution2.size());
  // container of means
  std::vector<double> RPMeans = {RPMean1, RPMean2};

  // ------------------------------------------------------------------------------------------------------------------------------

  // checking number of arguments
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " input file name> <output file name> <max events=-1>" << std::endl;
    std::exit(-1);
  }
  // reading argument ~ input/output file names
  std::string inFileName(argv[1]);
  std::string outFileName(argv[2]);
  std::cout << "Writing to " << outFileName << std::endl;
  // max number of events to process
  int NMaxEvent = -1;
  if (argc >= 4)
    NMaxEvent = atoi(argv[3]);
  if (NMaxEvent < 1)
    NMaxEvent = -1;

  // ------------------------------------------------------------------------------------------------------------------------------

  // CREATE HISTOGRAMS TO BE FILLED BY LOOPING THROUGH ALL EVENTS
  // transverse momentum ditsribution
  TH1 *pTDistribution = new TH1D("pTDistribution", "pT distribution", 100, 0, 2);
  // centrality distribution
  TH1 *centralityDistribution = new TH1D("centralityDistribution", "Centrality distribution", 100, 0, 100);
  // azimuthal distribution(s) for given pT splittings and centrality range
  std::vector<std::vector<TH1D *>> azimuthDistribution(2, std::vector<TH1D *>(NpT));
  // vector of vector ~ possible update to row/column-major vector (better memory handling)
  for (int iCentr = 0; iCentr < NC; iCentr++)
    for (int ipT = 0; ipT < NpT; ipT++)
      azimuthDistribution[iCentr][ipT] = new TH1D(Form("azimuthDistribution_Centrality%i_pT%i", iCentr, ipT),
                                                  "Azimuthal distribution", 100, -M_PI_2, M_PI_2);

  // ------------------------------------------------------------------------------------------------------------------------------

  // INITIALIZE particle_tree OBJECT
  particle_tree p(inFileName.c_str());
  if (p.fChain)
    std::cout << "Tree initialized" << std::endl;
  else
  {
    std::cout << "No tree found." << std::endl;
    std::exit(-1);
  }

  // ------------------------------------------------------------------------------------------------------------------------------

  // DETERMINE HOW MANY EVENTS TO RUN ON
  long unsigned int NEvents = p.fChain->GetEntries();
  if (NMaxEvent > 0 && NMaxEvent < (int)NEvents)
    NEvents = NMaxEvent;
  std::cout << "Will run on " << NEvents << " events (out of " << p.fChain->GetEntries() << ")." << std::endl;

  // ------------------------------------------------------------------------------------------------------------------------------

  // LOOP THROUGH EVENTS IN THE GIVEN DATASET
  for (long unsigned int iEvent = 0; iEvent < NEvents; iEvent++)
  {
    // MONITOR PROGRESS THROUGH STDERR OUTPUT
    if (iEvent > 0 && iEvent % 1000 == 0)
      std::cout << ".";
    if (iEvent > 0 && iEvent % 10000 == 0)
      std::cout << "Analyzing event #" << iEvent << std::endl;

    // ------------------------------------------------------------------------------------------------------------------------------

    // LOAD EVENT DATA INTO particle_tree OBJECT
    p.GetEntry(iEvent);

    // ------------------------------------------------------------------------------------------------------------------------------

    // reaction plane
    double reactionPlane = p.ReactionPlane;
    // centrality
    double centrality = p.Centrality;
    centralityDistribution->Fill(centrality);

    // ------------------------------------------------------------------------------------------------------------------------------

    // choose centrality range (initialize to some meaningless value) ~ [0, 30] vs. [40, 70]
    int centralityRange = -1;
    if (centrality >= centrality1.first && centrality <= centrality1.second)
      centralityRange = 0;
    else if (centrality >= centrality2.first && centrality <= centrality2.second)
      centralityRange = 1;
    else
      continue;

    // ------------------------------------------------------------------------------------------------------------------------------

    // LOOP THROUGH ALL PARTICLES OF THE GIVEN EVENT
    for (int iPart = 0; iPart < p.Ntracks; iPart++)
    {
      // CALCULATE SOME VARIABLES AND FILL THEM INTO A HISTOGRAMS
      // transverse momentum
      double pT = sqrt(p.px[iPart] * p.px[iPart] + p.py[iPart] * p.py[iPart]);
      pTDistribution->Fill(pT);
      if (pT > 2.)
        continue;

      // ------------------------------------------------------------------------------------------------------------------------------

      // choose transverse momentum range (initialize to some meaningless value)
      int pTRange = -1;
      double pTTemp = 0.1;
      while (pTTemp <= pT)
        pTTemp += 0.1, pTRange++;

      // ------------------------------------------------------------------------------------------------------------------------------

      // calculate azimuthal angle for given particle (in lab)
      double phi = std::atan2(p.py[iPart], p.px[iPart]);
      // calculate azimuthal angle for given particle (in reaction plane)
      double phiRP = phi - reactionPlane;
      while(phiRP > M_PI_2) phiRP -= M_PI;
      while(phiRP < -M_PI_2) phiRP += M_PI;

      // add to specific histogram
      azimuthDistribution[centralityRange][pTRange]->Fill(phiRP);
    }
  }
  // line break
  std::cout << std::endl;

  // ------------------------------------------------------------------------------------------------------------------------------

  // MAKE FIT(S) to determine elliptic flow (v2)
  // define ansatz ~ 1D Fourier expansion in the azimuthal angle [-pi, pi]
  TF1 *FourierFitFunc = new TF1("FourierFit", "[0] + [1] * 2 * cos(2 * x)", -M_PI_2, M_PI_2);
  // loop through centralities and pT splittings
  for (int iCentr = 0; iCentr < NC; iCentr++)
  {
    // save results to...
    v2Result[iCentr] = new TGraph(NpT);
    v2Result[iCentr]->SetName(Form("v_{2} with centrality %i-%i [%]", centralities[iCentr].first, centralities[iCentr].second));
    v2Error[iCentr] = new TGraphErrors(NpT);
    v2Error[iCentr]->SetName(Form("v_{2} errors with centrality %i-%i [%]", centralities[iCentr].first, centralities[iCentr].second));

    for (int ipT = 0; ipT < NpT; ipT++)
    {
      // make fit
      TFitResultPtr r = azimuthDistribution[iCentr][ipT]->Fit(FourierFitFunc, "S");

      // constant term
      double A = r->Parameter(0);
      // coefficient for 2 * cos(2 * x)
      double B = r->Parameter(1);
      // get errors
      double AErr = r->ParError(0);
      double BErr = r->ParError(1);
      // get covariance
      double ABCov = r->CovMatrix(1, 0);

      // elliptic flow (v2) and correction via reaction plane resolution
      double v2 = B / A / RPMeans[iCentr];
      // error estimation through error propagation
      double v2Err = std::abs(B / A) * std::sqrt((AErr * AErr) / (A * A) + (BErr * BErr) / (B * B) - 2 * ABCov / (A * B)) / RPMeans[iCentr];

      std::cout << v2 << " +/-" << v2Err << std::endl;

      v2Result[iCentr]->SetPoint(ipT, pTSplitting[ipT], v2);
      v2Error[iCentr]->SetPoint(ipT, pTSplitting[ipT], v2);
      v2Error[iCentr]->SetPointError(ipT, 0., v2Err);
    }
  }
  // ------------------------------------------------------------------------------------------------------------------------------

  // WRITE ALL HISTOGRAMS TO THE OUTPUT ROOT FILE
  TFile *f = new TFile(outFileName.c_str(), "RECREATE");
  if (!f->IsWritable())
    std::cout << "File " << outFileName << " was not opened!" << std::endl;
  else
    std::cout << "Analysis done, writing histograms to " << outFileName << std::endl;
  f->cd();
  pTDistribution->Write();
  centralityDistribution->Write();
  for (int i = 0; i < NC; i++)
  {
    v2Result[i]->Write();
    v2Error[i]->Write();
  }
  f->Write();
  f->Close();
}
