#include <TFile.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TLegend.h>
#include <TString.h>
#include <RooCurve.h>
#include <TPaveText.h>
#include <iostream>
#include <RooHist.h> 
#include <TPaveText.h>
#include <TPave.h>
#include "TAxis.h"

void plot_rescaled() {
    // Apri il file ROOT
    TFile *file = TFile::Open("/Users/saragaretti/dq_fitter_tails/analysis/muonLowPt210SigmaPDCA/output_chi2/Centr/0_90/Centr_0.0_90.0/multi_trial_NA60_VWG_1.05_width_2_8_MC_tails_MC_tails.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: impossible to open ROOT file!" << std::endl;
        return;
    }

    // Prendi il Canvas
    TCanvas *c = (TCanvas*)file->Get("fit_plot_NA60_NA60_VWG__2.0_5.0_1.05_2_8_data_tails_Mass_BkgSub_Int_ME");
    if (!c) {
        std::cerr << "Error: TCanvas not found!" << std::endl;
        return;
    }

    RooCurve *psi2S = nullptr;
    RooCurve *jpsi = nullptr;
    RooCurve *bkg = nullptr;
    RooHist *dataHist = nullptr;

    TIter next(c->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next())) {
        if (obj->InheritsFrom(RooCurve::Class())) {
            RooCurve *curve = (RooCurve*)obj;
            TString name = curve->GetName();
            if (name.Contains("Psi2sPdf")) {
                psi2S = curve;
            }
            if (name.Contains("JpsiPdf")) {
                jpsi = curve;
            }
            if (name.Contains("BkgPdf")) {
                bkg = curve;
            }
            if (obj->InheritsFrom(RooHist::Class())) { 
                dataHist = (RooHist*)obj;
            }
        }
    }

    if (!psi2S) {
        std::cerr << "Error: Ψ(2S) RooCurve not found!" << std::endl;
        return;
    }
    if (!jpsi) {
        std::cerr << "Error: J/Ψ RooCurve not found!" << std::endl;
        return;
    }
    if (!bkg) {
        std::cerr << "Error: background RooCurve not found!" << std::endl;
        return;
    }
    RooCurve *totalFit_notRescaled = new RooCurve(*bkg);

    for (int i = 0; i < bkg->GetN(); i++) {
        double x_bkg, y_bkg;
        bkg->GetPoint(i, x_bkg, y_bkg);

        double y_jpsi = 0;
        double y_psi2S = 0;

        for (int j = 0; j < jpsi->GetN(); j++) {
            double x_jpsi, y_jpsi_point;
            jpsi->GetPoint(j, x_jpsi, y_jpsi_point);
            if (x_bkg == x_jpsi) {
                y_jpsi = y_jpsi_point;
                break;
            }
        }
        for (int k = 0; k < psi2S->GetN(); k++) {
            double x_psi2S, y_psi2S_point;
            psi2S->GetPoint(k, x_psi2S, y_psi2S_point);
            if (x_bkg == x_psi2S) {
                y_psi2S = y_psi2S_point;
                break;
            }
        }

        double y_total = y_bkg + y_jpsi + y_psi2S;
        totalFit_notRescaled->SetPoint(i, x_bkg, y_total);
    }

    RooCurve *psi2S_scaled = new RooCurve(*psi2S);
    psi2S_scaled->SetName("Psi2S_scaled");
    psi2S_scaled->SetTitle("Riscalata #Psi(2S)");

    for (int i = 0; i < psi2S_scaled->GetN(); i++) {
        double x, y;
        psi2S_scaled->GetPoint(i, x, y);
        psi2S_scaled->SetPoint(i, x, y * 0.0200823 / 0.00922245);
    }

    RooCurve *totalFit_modified = new RooCurve(*bkg); 

    for (int i = 0; i < bkg->GetN(); i++) {
        double x_bkg, y_bkg;
        bkg->GetPoint(i, x_bkg, y_bkg);

        double y_jpsi = 0;
        double y_psi2S = 0;

        for (int j = 0; j < jpsi->GetN(); j++) {
            double x_jpsi, y_jpsi_point;
            jpsi->GetPoint(j, x_jpsi, y_jpsi_point);
            if (x_bkg == x_jpsi) {
                y_jpsi = y_jpsi_point;
                break;
            }
        }

        for (int k = 0; k < psi2S_scaled->GetN(); k++) {
            double x_psi2S, y_psi2S_point;
            psi2S_scaled->GetPoint(k, x_psi2S, y_psi2S_point);
            if (x_bkg == x_psi2S) {
                y_psi2S = y_psi2S_point;
                break;
            }
        }

        double y_total = y_bkg + y_jpsi + y_psi2S;
        totalFit_modified->SetPoint(i, x_bkg, y_total);
    }

    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.15);

    totalFit_modified->GetXaxis()->SetTitle("M [#mu^{+} #mu^{-}] (GeV)");
    totalFit_modified->GetYaxis()->SetTitle("Counts");
    totalFit_modified->GetXaxis()->SetTitleSize(0.045);
    totalFit_modified->GetYaxis()->SetTitleSize(0.045);
    totalFit_modified->GetXaxis()->SetLabelSize(0.035);
    totalFit_modified->GetYaxis()->SetLabelSize(0.035);

    c->cd();
    gPad -> SetLogy(1);
    c->Draw();

    totalFit_modified->SetLineColor(kBlue);
    totalFit_modified->SetLineStyle(9);  
    totalFit_modified->SetLineWidth(6);
    totalFit_modified->Draw("L SAME");

    totalFit_notRescaled->SetLineColor(kRed);  
    totalFit_notRescaled->SetLineStyle(1);     
    totalFit_notRescaled->SetLineWidth(6);

    totalFit_notRescaled->Draw("L SAME");

    

    jpsi->SetLineColor(kAzure+7);
    jpsi->SetLineStyle(1);
    jpsi->SetLineWidth(4);
    jpsi->Draw("L SAME");

    psi2S->SetLineColor(kGreen-3);
    psi2S->SetLineStyle(1);
    psi2S->SetLineWidth(4);
    psi2S->Draw("L SAME");

    TLegend *leg = new TLegend(0.58, 0.60, 0.75, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(dataHist, "Data", "lep");
    leg->AddEntry(totalFit_notRescaled, "Pb-Pb total fit", "l");
    leg->AddEntry(jpsi, "J/#psi signal", "l");
    leg->AddEntry(psi2S, "#psi(2S) signal", "l");
    leg->AddEntry(totalFit_modified, "pp overlaid", "l");
    leg->Draw("SAME");

    TPaveText *pt = new TPaveText(0.20, 0.91, 0.40, 0.99, "NDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.025);
    pt->AddText("Centr: (0-90)%");
    pt->AddText("p_{T} < 20 GeV/c");
    pt->AddText("p_{T}^{#mu} > 0.7 GeV/c");
    pt->AddText("2.5 < y < 4");
    pt->Draw();

    gPad->SetTicks(1, 1);

    c->Update();
}



