// plot function
void Plot_analyzetree(const char *filename = "analyzetree.root", const char *figdir = "figs")
{
  // file
  TFile *f = new TFile(filename);
  // canvas
  TCanvas *c1 = new TCanvas();
  // plot
  TH1 *pTDistribution = (TH1F *)f->Get("pTDistribution");
  TH1 *centralityDistribution = (TH1F *)f->Get("centralityDistribution");

  // transverse momenta
  pTDistribution->SetTitle("N(p_{T})");
  pTDistribution->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  pTDistribution->GetYaxis()->SetTitle("N(p_{T}) [c/GeV]");
  pTDistribution->Draw("e");
  c1->Print(Form("%s/pTDistribution.png", figdir));
  c1->Clear();

  // centrality
  centralityDistribution->SetTitle("Centrality distribution");
  centralityDistribution->GetXaxis()->SetTitle("Centrality [%]");
  centralityDistribution->GetYaxis()->SetTitle("N(centrality)");
  centralityDistribution->Draw("e");
  c1->Print(Form("%s/centrality.png", figdir));
  c1->Clear();

  
}
