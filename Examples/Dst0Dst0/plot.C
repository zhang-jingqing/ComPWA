#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>

void plot() {
  TFile *f_data = new TFile("output/plot_data.root", "READ");
  TTree *t_data = (TTree *)f_data->Get("data");
  double event_weight_data;
  t_data->SetBranchAddress("weight", &event_weight_data);
  double mSq_1_2_data, theta_1_2_data, phi_1_2_data;
  t_data->SetBranchAddress("mSq_1_2", &mSq_1_2_data);
  t_data->SetBranchAddress("theta_1_2", &theta_1_2_data);
  t_data->SetBranchAddress("phi_1_2", &phi_1_2_data);

  TFile *f_fit = new TFile("output/plot_fit.root", "READ");
  TTree *t_fit = (TTree *)f_fit->Get("intensity_weighted_phspdata");
  double event_weight_fit, intensity_fit;
  t_fit->SetBranchAddress("weight", &event_weight_fit);
  t_fit->SetBranchAddress("intensity", &intensity_fit);
  double mSq_1_2_fit, theta_1_2_fit, phi_1_2_fit;
  t_fit->SetBranchAddress("mSq_1_2", &mSq_1_2_fit);
  t_fit->SetBranchAddress("theta_1_2", &theta_1_2_fit);
  t_fit->SetBranchAddress("phi_1_2", &phi_1_2_fit);

  int N_data = t_data->GetEntries();
  int N_fit = t_fit->GetEntries();
  double scale = 1.0; //global scale
  //N_data_net = scale * \sum_i intensity_fit (if no weight for phsp)
  //N_data_net = scale * \sum_i intensity_fit * event_fit (if has weight)
  double N_data_net = 0;
  for (int ievt = 0; ievt < N_data; ++ievt) {
    t_data->GetEntry(ievt);
    N_data_net += event_weight_data; 
  }
  double N_fit_net = 0;
  for (int ievt = 0; ievt < N_fit; ++ievt) {
    t_fit->GetEntry(ievt);
    N_fit_net += event_weight_data * intensity_fit;
  }
  scale = N_data_net/N_fit_net;
  std::cout << " N_data_net = " << N_data_net << " N_fit_net = "
      << N_fit_net << " scale = N_data_fit/N_fit_net = "
      << scale << std::endl;

  int N_cos_theta = 50;
  double low_cos_theta = -1, up_cos_theta = 1;

  TH1F *h_cos_theta_1_2_data = new TH1F("h_cos_theta_1_2_data",
      "h_cos_theta_1_2_data", N_cos_theta, low_cos_theta, up_cos_theta);
  h_cos_theta_1_2_data->SetXTitle("cos(#theta)");
  h_cos_theta_1_2_data->SetYTitle("Events");
  h_cos_theta_1_2_data->SetLineColor(1);

  for (int ievt = 0; ievt < N_data; ++ievt) {
    t_data->GetEntry(ievt);
    h_cos_theta_1_2_data->Fill(cos(theta_1_2_data), event_weight_data);
  }

  TH1F *h_cos_theta_1_2_fit = new TH1F("h_cos_theta_1_2_fit",
      "h_cos_theta_1_2_fit", N_cos_theta, low_cos_theta, up_cos_theta);
  h_cos_theta_1_2_fit->SetXTitle("cos(#theta)");
  h_cos_theta_1_2_fit->SetYTitle("Events");
  h_cos_theta_1_2_fit->SetLineColor(2);
  h_cos_theta_1_2_fit->SetLineStyle(2);

  for (int ievt = 0; ievt < N_fit; ++ievt) {
    t_fit->GetEntry(ievt);
    h_cos_theta_1_2_fit->Fill(cos(theta_1_2_fit),
        event_weight_fit * intensity_fit * scale);
  }
  
  std::cout << " ndata = " << h_cos_theta_1_2_data->Integral()
      << " nfit = " << h_cos_theta_1_2_fit->Integral() << std::endl;

  TLegend leg(0.3, 0.3, 0.5, 0.5);
  leg.AddEntry(h_cos_theta_1_2_data, "data", "ple");
  leg.AddEntry(h_cos_theta_1_2_fit, "Fit", "l");

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();
  h_cos_theta_1_2_data->Draw("pe");
  h_cos_theta_1_2_fit->Draw("HISTSAME");
  leg.Draw();
  c1->Print("figures/cos_theta_1_2.eps");

}
