#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMinuit.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include <fstream>
#include "THnSparse.h"
#include <TStyle.h>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>
#include <iterator>
#include <numeric>
#include <map>
#include "THStack.h"
#include <list>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <string_view>
#include "Fit/FitResult.h"
#include "TEllipse.h"
#include <tuple>
#include <TLegend.h>
#include "TProfile.h"
#include <TArrayD.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

using namespace std;

void LoadStyle();
void SetLegend(TLegend *);
void PlotSystErrors(const vector<double>&, const int, const vector<double>&, const vector<double>&, int);


void plot_results() {
    LoadStyle();

    const Double_t BrJpsiToMuMu = 0.05961;
    const Double_t BrPsi2sToMuMu = 0.008;
    const Double_t BrPsi2sToJpsiRatio = BrPsi2sToMuMu / BrJpsiToMuMu;

    TFile *fInShm = new TFile("psi2S2Jpsi_SHMc.root", "READ");
    TGraphErrors *graShmCor1Psi2sToJpsiRatioVsNpart = (TGraphErrors*) fInShm -> Get("psi2S2Jpsi_cor1");
    graShmCor1Psi2sToJpsiRatioVsNpart -> Scale(BrPsi2sToJpsiRatio);

    TGraphErrors *graShmCor2Psi2sToJpsiRatioVsNpart = (TGraphErrors*) fInShm -> Get("psi2S2Jpsi_cor2");
    graShmCor2Psi2sToJpsiRatioVsNpart -> Scale(BrPsi2sToJpsiRatio);

    const int n_TAMU = 29;

    Double_t Npart_TAMU_low[n_TAMU]={7.135433,11.518662,17.575771,25.432085,35.100117,46.516098,59.574142,74.144707,
    90.087494,107.259377,125.506401,144.669342,164.589844,185.09343,205.998627,227.12294,248.25795,269.18103,
    289.661224,309.43277,328.19632,345.635101,361.405121,375.141785,386.523132,395.330017,401.485535,405.069427,406.237549};

    Double_t TAMU_low[n_TAMU]={0.020186925320993143,0.01828644214016006,0.014835554809622198,0.011862822632827667,0.009857535890000912,
    0.008890723354884285,0.009020847076614425,0.009231138291748441,0.009431897227547788,0.009591002204358212,0.009731286872581868,
    0.00980457858747916,0.009831534254048223,0.009853643888792587,0.009828721510548168,0.009868223283023728,0.009828358164944163,
    0.009759878775327281,0.009694168434728694,0.009575054206910055,0.009509427703169494,0.009454868766140193,0.009408864580117202,
    0.009376827335054473,0.009336667294088214,0.009321536754976092,0.00931448337761218,0.009280747315049432,0.00927319336005245};

    Double_t Npart_TAMU_high[n_TAMU]={7.135433,11.518662,17.575771,25.432085,35.100117,46.516098,59.574142,74.144707,
    90.087494,107.259377,125.506401,144.669342,164.589844,185.09343,205.998627,227.12294,248.25795,269.18103,
    289.661224,309.43277,328.19632,345.635101,361.405121,375.141785,386.523132,395.330017,401.485535,405.069427,406.237549};

    Double_t TAMU_high[n_TAMU]={0.020186925320993143,0.01828644214016006,0.014960223337434148,0.012382070761537292,0.010848796482291508,
    0.00969531370374259,0.010127883036265283,0.010624948776393625,0.01101407831947193,0.011300018990594528,0.011523892349110106,
    0.011639636608564803,0.011688241309726564,0.011718819329702201,0.011689572325013164,0.011729212042654329,0.01166697117526832,
    0.011569679216458802,0.011467614911642207,0.011305963858968137,0.011207209581627929,0.01111908627917895,0.011043066624848077,
    0.010988740351986203,0.010928467161579906,0.01089988624690345,0.010872606142224293,0.01083338430602772,0.010824784475263497}; 

    Double_t Npart_TAMU_centr[29], Npart_TAMU_err_low[29], Npart_TAMU_err_high[29];
    Double_t TAMU_centr[29], TAMU_err_low[29], TAMU_err_high[29];


    for (int iBin = 0;iBin < n_TAMU;iBin++) {
        Npart_TAMU_centr[iBin] = (Npart_TAMU_low[iBin] + Npart_TAMU_high[iBin]) / 2.;
        Npart_TAMU_err_low[iBin] = Npart_TAMU_centr[iBin] - Npart_TAMU_low[iBin];
        Npart_TAMU_err_high[iBin] = Npart_TAMU_high[iBin] - Npart_TAMU_centr[iBin];

        TAMU_centr[iBin] = (TAMU_low[iBin] + TAMU_high[iBin]) / 2.;
        TAMU_err_low[iBin] = TAMU_centr[iBin] - TAMU_low[iBin];
        TAMU_err_high[iBin] = TAMU_high[iBin] - TAMU_centr[iBin];
    }


    TGraphAsymmErrors *graTamuPsi2sToJpsiRatioVsNpart = new TGraphAsymmErrors(n_TAMU, Npart_TAMU_centr, TAMU_centr, Npart_TAMU_err_low, Npart_TAMU_err_high, TAMU_err_low, TAMU_err_high); 
    graTamuPsi2sToJpsiRatioVsNpart -> SetFillStyle(0);
    graTamuPsi2sToJpsiRatioVsNpart -> SetFillColorAlpha(kRed+1, 0.7);
    graTamuPsi2sToJpsiRatioVsNpart -> SetLineColorAlpha(kRed+1, 0.7);
    graTamuPsi2sToJpsiRatioVsNpart -> SetMarkerStyle(20);

    TGraphAsymmErrors *graTamuLowPsi2sToJpsiRatioVsNpart = new TGraphAsymmErrors(n_TAMU, Npart_TAMU_centr, TAMU_low, Npart_TAMU_err_low, Npart_TAMU_err_high, 0, 0); 
    graTamuLowPsi2sToJpsiRatioVsNpart -> SetLineColorAlpha(kRed+1, 0.9);

    TGraphAsymmErrors *graTamuHighPsi2sToJpsiRatioVsNpart = new TGraphAsymmErrors(n_TAMU, Npart_TAMU_centr, TAMU_high, Npart_TAMU_err_low, Npart_TAMU_err_high, 0, 0); 
    graTamuHighPsi2sToJpsiRatioVsNpart -> SetLineColorAlpha(kRed+1, 0.9);

    vector<double> ratioRun2 = {0.0117, 0.0097, 0.0107, 0.0109};
    vector<double> statratioRun2 = {0.0020, 0.0022, 0.0026, 0.0028};
    vector<double> systratioRun2 = {0.0027, 0.0021, 0.0017, 0.0024};

    //Medium MC syst
    vector<double> ratioRun3 = {0.0073, 0.0090, 0.0081, 0.0112};
    vector<double> statratioRun3 = {0.0018, 0.0022, 0.0026, 0.0028};
    vector<double> systratioRun3 = {0.001014, 0.000730, 0.001412, 0.001227};

    vector<double> ratiopp = {0.0200823};
    double ratiopp_double = 0.0200823;
    vector<double> statratiopp = {0.000792887};
    double statratiopp_double = 0.000792887;
    vector<double> systratiopp = {0.00177286};
    double systratiopp_double = 0.00177286;
    Double_t Npart_data_pp[1] = {2};
    
    Double_t Npart_data_Run2[4] = {316.7, 167.4, 77.74, 25.0};
    Double_t Npart_data[4] = {309.7, 159.4, 70.74, 18.0};
    Double_t Npart_data_error[4] = { 0, 0, 0, 0};

    vector<double> Npart_data_vector = {309.7, 159.4, 70.74, 18.0};
    vector<double> Npart_data_vector_Run2 = {316.7, 167.4, 77.74, 25.0};
    vector<double> Npart_data_vector_error = {0, 0, 0, 0};
    vector<double> Npart_data_vector_error_Run2 = {0, 0, 0, 0};

    TGraphErrors *graRun2 = new TGraphErrors(4, Npart_data_Run2, &ratioRun2[0], 0, &statratioRun2[0]);
    graRun2->SetMarkerStyle(20);
    graRun2->SetMarkerColor(kBlue);
    graRun2->SetLineColor(kBlue);

    TGraphErrors *graRun3 = new TGraphErrors(4, Npart_data, &ratioRun3[0], 0, &statratioRun3[0]);
    graRun3->SetMarkerStyle(21);
    graRun3->SetMarkerColor(kRed);
    graRun3->SetLineColor(kRed);

    TGraphErrors *graPP = new TGraphErrors(1, Npart_data_pp, &ratiopp[0], 0, &statratiopp[0]);
    graPP->SetMarkerStyle(20);
    graPP->SetMarkerColor(kGreen+3);
    graPP->SetLineColor(kGreen+3);


    // Creazione Canvas
    TCanvas *canvasPsi2sToJpsiRatioVsNpart = new TCanvas("canvasPsi2sToJpsiRatioVsNpart", "", 800, 600);
    canvasPsi2sToJpsiRatioVsNpart -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtRun2VsRun3 = new TH2D("histGridV2JpsiVsPtRun2VsRun3", "", 100, -10, 400, 100, 0, 0.025);
    histGridV2JpsiVsPtRun2VsRun3 -> GetXaxis() -> SetTitle("<#it{N}_{part}>");
    histGridV2JpsiVsPtRun2VsRun3 -> GetYaxis() -> SetTitle("BR_{#psi(2S)#rightarrow#mu^{+}#mu^{-}} #sigma_{#psi(2S)} / BR_{J/#psi#rightarrow#mu^{+}#mu^{-}} #sigma_{J/#psi}");
    histGridV2JpsiVsPtRun2VsRun3 -> GetYaxis() -> SetTitleOffset(1.6);
    histGridV2JpsiVsPtRun2VsRun3 -> Draw();
    graTamuPsi2sToJpsiRatioVsNpart -> Draw("E3 SAME");
    graTamuLowPsi2sToJpsiRatioVsNpart -> Draw("L SAME");
    graTamuHighPsi2sToJpsiRatioVsNpart -> Draw("L SAME");
    graShmCor1Psi2sToJpsiRatioVsNpart -> Draw("L SAME");
    graShmCor2Psi2sToJpsiRatioVsNpart -> Draw("L SAME");
    PlotSystErrors(Npart_data_vector_Run2, 4, ratioRun2, systratioRun2, kBlue);
    PlotSystErrors(Npart_data_vector, 4, ratioRun3, systratioRun3, kRed);
    double ylow = ratiopp_double - systratiopp_double;
    double yhigh = ratiopp_double + systratiopp_double;
    TBox *box = new TBox(-2, ylow, 6, yhigh);
    box->SetFillStyle(0);
    box->SetLineColor(kGreen+3);
    box->SetLineWidth(2);
    box->Draw("same");
    graRun2->Draw("P SAME");
    graRun3->Draw("P SAME");
    graPP->Draw("P SAME");
    graShmCor1Psi2sToJpsiRatioVsNpart -> Draw("L SAME");
    graShmCor2Psi2sToJpsiRatioVsNpart -> Draw("L SAME");

    TLegend *legendPsi2sToJpsiRatioVsNpart = new TLegend(0.55, 0.65, 0.90, 0.90, " ", "brNDC");
    SetLegend(legendPsi2sToJpsiRatioVsNpart);
    legendPsi2sToJpsiRatioVsNpart->SetTextSize(0.03);
    legendPsi2sToJpsiRatioVsNpart->SetEntrySeparation(0.03); 
    legendPsi2sToJpsiRatioVsNpart -> AddEntry(graTamuPsi2sToJpsiRatioVsNpart, "TAMU, #sqrt{#it{s}_{NN}} = 5.02 TeV", "F");
    legendPsi2sToJpsiRatioVsNpart -> AddEntry(graShmCor1Psi2sToJpsiRatioVsNpart, "SHMc, #sqrt{#it{s}_{NN}} = 5.02 TeV", "L");
    legendPsi2sToJpsiRatioVsNpart -> AddEntry(graRun2, "Run 2, #sqrt{#it{s}_{NN}} = 5.02 TeV", "P");
    legendPsi2sToJpsiRatioVsNpart -> AddEntry(graRun3, "Run 3, #sqrt{#it{s}_{NN}} = 5.36 TeV", "P");
    legendPsi2sToJpsiRatioVsNpart -> AddEntry(graPP, "Extr. pp, #sqrt{#it{s}} = 5.36 TeV", "P");
    legendPsi2sToJpsiRatioVsNpart -> Draw();


    canvasPsi2sToJpsiRatioVsNpart -> SaveAs("Psi2S_to_Jpsi_ratio_vs_npart_Fin.pdf");
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LoadStyle(){
    int font = 42;
    gStyle -> SetFrameBorderMode(0);
    gStyle -> SetFrameFillColor(0);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(0);
    gStyle -> SetPadColor(10);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetTitleFillColor(10);
    gStyle -> SetTitleBorderSize(1);
    gStyle -> SetStatColor(10);
    gStyle -> SetStatBorderSize(1);
    gStyle -> SetLegendBorderSize(1);
    gStyle -> SetDrawBorder(0);
    gStyle -> SetTextFont(font);
    gStyle -> SetStatFontSize(0.05);
    gStyle -> SetStatX(0.97);
    gStyle -> SetStatY(0.98);
    gStyle -> SetStatH(0.03);
    gStyle -> SetStatW(0.3);
    gStyle -> SetTickLength(0.02,"y");
    gStyle -> SetEndErrorSize(3);
    gStyle -> SetLabelSize(0.05,"xyz");
    gStyle -> SetLabelFont(font,"xyz");
    gStyle -> SetLabelOffset(0.01,"xyz");
    gStyle -> SetTitleFont(font,"xyz");
    gStyle -> SetTitleOffset(0.9,"x");
    gStyle -> SetTitleOffset(1.02,"y");
    gStyle -> SetTitleSize(0.05,"xyz");
    gStyle -> SetMarkerSize(1.3);
    gStyle -> SetOptStat(0);
    gStyle -> SetEndErrorSize(0);
    gStyle -> SetCanvasPreferGL(kTRUE);
    gStyle -> SetHatchesSpacing(0.5);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetPadBottomMargin(0.15);
    gStyle -> SetPadTopMargin(0.05);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetEndErrorSize(0.0);
    gStyle -> SetTitleSize(0.05,"X");
    gStyle -> SetTitleSize(0.045,"Y");
    gStyle -> SetLabelSize(0.045,"X");
    gStyle -> SetLabelSize(0.045,"Y");
    gStyle -> SetTitleOffset(1.2,"X");
    gStyle -> SetTitleOffset(1.35,"Y");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.04);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PlotSystErrors(const vector<double>& x, const int xErr, const vector<double>& y, const vector<double>& yErr, int color) {
    for (size_t i = 0; i < x.size(); i++) {
        double xlow = x[i] - xErr;
        double xhigh = x[i] + xErr;
        double ylow = y[i] - yErr[i];
        double yhigh = y[i] + yErr[i];

        TBox *box = new TBox(xlow, ylow, xhigh, yhigh);
        box->SetFillStyle(0);
        box->SetLineColor(color);
        box->SetLineWidth(2);
        box->Draw("same");
    }
}