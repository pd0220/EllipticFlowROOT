// plot function
void Plot_analyzetree(const char *filename = "analyzetree.root", const char *figdir = "figs")
{
  // file
  TFile *f = new TFile(filename);
  // canvas
  TCanvas *c1 = new TCanvas("canvas", "", 1200, 600);
  c1->SetGrid();
  // plots
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

  // v2
  TGraphErrors *v2_C1 = (TGraphErrors *)f->Get("v_{2} errors with centrality 0-30 [%]");
  v2_C1->SetName("v2_C1");
  v2_C1->GetXaxis()->SetTitle("p_{T} [GeV]");
  v2_C1->GetXaxis()->SetTitleOffset(1);
  v2_C1->GetYaxis()->SetTitle("v_{2}");
  v2_C1->GetXaxis()->SetTitleSize(0.045);
  v2_C1->GetYaxis()->SetTitleSize(0.055);
  v2_C1->SetTitle("Elliptic flow measurement results.");
  v2_C1->GetXaxis()->SetRangeUser(0.1, 2.1);
  v2_C1->GetYaxis()->SetRangeUser(0., 0.4);
  v2_C1->SetMarkerStyle(20);
  v2_C1->SetMarkerColor(2);
  v2_C1->SetMarkerSize(1.5);
  v2_C1->Draw("apl");

  TGraphErrors *v2_C2 = (TGraphErrors *)f->Get("v_{2} errors with centrality 40-70 [%]");
  v2_C2->SetName("v2_C2");
  v2_C2->SetMarkerStyle(21);
  v2_C2->SetMarkerColor(3);
  v2_C2->SetMarkerSize(1.5);
  v2_C2->Draw("pl same");

  TGraphErrors *v2_C1_NonCorr = (TGraphErrors *)f->Get("v_{2} errors with centrality 0-30 [%] (without correction)");
  v2_C1_NonCorr->SetName("v2_C1_NonCorr");
  v2_C1_NonCorr->SetMarkerStyle(22);
  v2_C1_NonCorr->SetMarkerColor(4);
  v2_C1_NonCorr->SetMarkerSize(1.5);
  v2_C1_NonCorr->Draw("pl same");

  TGraphErrors *v2_C2_NonCorr = (TGraphErrors *)f->Get("v_{2} errors with centrality 40-70 [%] (without correction)");
  v2_C2_NonCorr->SetName("v2_C2_NonCorr");
  v2_C2_NonCorr->SetMarkerStyle(23);
  v2_C2_NonCorr->SetMarkerColor(6);
  v2_C2_NonCorr->SetMarkerSize(1.5);
  v2_C2_NonCorr->Draw("pl same");

  auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  legend->AddEntry("v2_C1", "Centrality: 0-30%, RP correction", "lep");
  legend->AddEntry("v2_C2", "Centrality: 40-70%, RP correction", "lep");
  legend->AddEntry("v2_C1_NonCorr", "Centrality: 0-30%", "lep");
  legend->AddEntry("v2_C2_NonCorr", "Centrality: 40-70%", "lep");
  legend->Draw("pl same");

  c1->Print(Form("%s/v_2_test.pdf", figdir));
  c1->Clear();
}
