#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <bitset>

#include "TSystemDirectory.h"
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TKey.h"
#include "THashList.h"
#include "TProfile.h"
#include "TTreeReader.h"
#include "TLine.h"
#include "TGaxis.h"
#include <ROOT/RDataFrame.hxx>
#include "TH2F.h"
#include "TMinuit.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "THnSparse.h"
#include <TStyle.h>
#include "TSystem.h"
#include "TAxis.h"
#include "TBranch.h"
#include <vector>
#include <stdlib.h> 
#include <iterator>
#include <numeric>
#include <map>
#include "THStack.h"
#include <list>
#include <algorithm>
#include <cctype>
#include <string_view>
#include "Fit/FitResult.h"
#include "TEllipse.h"
#include <tuple>
#include <TLegend.h>
#include <TArrayD.h>

#endif

void LoadStyle();
void SetLegend(TLegend *);

const int kNsel = 51;

// Taken from O2Physics/Common/CCDB/EventSelectionParams.cxx
const char* selectionLabels[kNsel] = {
    "kIsBBV0A", "kIsBBV0C", "kIsBBFDA", "kIsBBFDC", "kIsBBT0A", "kIsBBT0C",
    "kNoBGV0A", "kNoBGV0C", "kNoBGFDA", "kNoBGFDC", "kNoBGT0A", "kNoBGT0C",
    "kIsBBZNA", "kIsBBZNC", "kIsBBZAC", "kNoBGZNA", "kNoBGZNC", "kNoV0MOnVsOfPileup",
    "kNoSPDOnVsOfPileup", "kNoV0Casymmetry", "kIsGoodTimeRange", "kNoIncompleteDAQ",
    "kNoTPCLaserWarmUp", "kNoTPCHVdip", "kNoPileupFromSPD", "kNoV0PFPileup", "kNoSPDClsVsTklBG",
    "kNoV0C012vsTklBG", "kNoInconsistentVtx", "kNoPileupInMultBins", "kNoPileupMV",
    "kNoPileupTPC", "kIsTriggerTVX", "kIsINT1", "kNoITSROFrameBorder", "kNoTimeFrameBorder",
    "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kIsVertexITSTPC", "kIsVertexTOFmatched",
    "kIsVertexTRDmatched", "kNoCollInTimeRangeNarrow", "kNoCollInTimeRangeStrict",
    "kNoCollInTimeRangeStandard", "kNoCollInTimeRangeVzDependent", "kNoCollInRofStrict",
    "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "kIsGoodITSLayer3",
    "kIsGoodITSLayer0123", "kIsGoodITSLayersAll"
};

constexpr uint64_t kTriggerMask =
    (1ULL << 37) |  // kIsGoodZvtxFT0vsPV
    (1ULL << 32) |  // kIsTriggerTVX
    (1ULL << 36) |  // kNoSameBunchPileup
    (1ULL << 34) |  // kNoITSROFrameBorder
    (1ULL << 35);   // kNoTimeFrameBorder

inline bool eventSelection(uint64_t selection) {
    return (selection & kTriggerMask) == kTriggerMask;
}

TH1D* ProjectTHnSparse(THnSparseD *, double , double , double , double , double , double );
THnSparseD* CutTHnSparse(THnSparseD *, double , double , double , double , double , double );
TH2D* CutTH2Y(TH2D *, string , double , double);
TH2D* CutTH2X(TH2D *, string , double , double);

TF1* PtJPsiPbPb5TeV_Func() {
    return new TF1("PtPsiPbPb5TeV", [](double* x, double* p) {
        Double_t px = x[0];
        Double_t p0 = p[0];
        Double_t p1 = p[1];
        Double_t p2 = p[2];
        Double_t p3 = p[3];

        return p0 * px / TMath::Power(1. + TMath::Power(px / p1, p2), p3);
    }, 0, 20, 4); // Range [0,20] con 4 parametri

}
static Double_t PtJPsiPbPb5TeV_tuned(const Double_t* px, const Double_t*){
    // jpsi pT in PbPb, tuned on data (2015) -> Castillo embedding https://alice.its.cern.ch/jira/browse/ALIROOT-8174?jql=text%20~%20%22LHC19a2%22
    Double_t x = *px;
    Float_t p0, p1, p2, p3;
    p0 = 1.00715e6;
    p1 = 3.50274;
    p2 = 1.93403;
    p3 = 3.96363;
    return p0 * x / TMath::Power(1. + TMath::Power(x / p1, p2), p3);
}

TF1* RapPsiPbPb5TeV_Func() {
    return new TF1("PtPsiPbPb5TeV", [](double* x, double* p) {
        Double_t px = x[0];
        Double_t p0 = p[0];
        Double_t p1 = p[1];
        Double_t p2 = p[2];
        Double_t p3 = p[3];

        return p0 * px / TMath::Power(1. + TMath::Power(px / p1, p2), p3);
    }, 0, 20, 4); // Range [0,20] con 4 parametri

}

inline void SetHistogram(TH1D *hist, Color_t color) {
    hist -> SetLineColor(color); 
    hist -> SetMarkerColor(color);
    hist -> SetMarkerStyle(20);
}

void AO2D_iterative() {
    //LoadStyle();
    bool isPbPb = false;
    double ptCut = 0.0;
    TDatabasePDG *database = TDatabasePDG::Instance();
    int muPdgCode = 13;
    int jpsiPdgCode = 443;
    double massMu = database -> GetParticle(muPdgCode) -> Mass();
    double massJpsi = database -> GetParticle(jpsiPdgCode) -> Mass();

    double weight = 1.0;
    TF1* weightFuncPt;
    TF1* weightFuncRap;
    std::vector<TF1*> weightFuncArrayRap_0_20;
    std::vector<TF1*> weightFuncArrayPt_0_20;
    std::vector<TF1*> weightFuncArrayRap_20_40;
    std::vector<TF1*> weightFuncArrayPt_20_40;
    std::vector<TF1*> weightFuncArrayRap_40_60;
    std::vector<TF1*> weightFuncArrayPt_40_60;
    std::vector<TF1*> weightFuncArrayRap_60_90;
    std::vector<TF1*> weightFuncArrayPt_60_90;
    std::vector<TF1*> weightFuncArrayRapGen_0_20;
    std::vector<TF1*> weightFuncArrayPtGen_0_20;
    std::vector<TF1*> weightFuncArrayRapGen_20_40;
    std::vector<TF1*> weightFuncArrayPtGen_20_40;
    std::vector<TF1*> weightFuncArrayRapGen_40_60;
    std::vector<TF1*> weightFuncArrayPtGen_40_60;
    std::vector<TF1*> weightFuncArrayRapGen_60_90;
    std::vector<TF1*> weightFuncArrayPtGen_60_90;
    std::map<std::tuple<double, double, double>, double> weightMap;
    std::map<int, double> impactParameterMap;

    uint64_t fSelection;
    float fMass, fPt, fEta, fTauz, fTauxy, fU2Q2, fCos2DeltaPhi, fR2EP, fR2SP, fCentFT0C, fImpactParameter = -99999;
    float fChi2pca, fSVertex, fEMC1, fEMC2, fPt1, fPt2, fPhi1, fPhi2, fEta1, fEta2, fPtMC1, fPtMC2, fPhiMC1, fPhiMC2, fPhi, fEtaMC1, fEtaMC2 = -99999;;
    int fSign, fSign1, fSign2 = -99999;
    UInt_t fMcDecision;
    Int_t fIsAmbig1, fIsAmbig2;

    std::vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+2};
    double ptBinsRun3[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0};
    double centrBins[] = {0.0, 20.0, 40.0, 60.0, 90.0};
    const int nBinsCentr = 4;
    //const int nBinsCentrRec = 10;
    const int nBinsPt = 8;
    const int nBinsRap = 6;
    double rapBins[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};

    TCanvas *canvasAxePt = new TCanvas("canvas_Axe_Pt", " Ax#it{#epsilon} vs #it{p}_{T}", 800, 600);
    TCanvas *canvasAxeRap = new TCanvas("canvas_Axe_Rap", " Ax#it{#epsilon} vs #it{y}", 800, 600);
    TCanvas *canvasAxeCentr = new TCanvas("canvas_Axe_Centr", " Ax#it{#epsilon} vs #it{CentrFT0}", 800, 600);
    std::vector<TH3F*> histCentrPtRapGenList;


    string pathToFile = "/Users/saragaretti/dq_fitter_tails/inputShapes/AO2D.root";
    string pathToFileOut = "/Users/saragaretti/dq_fitter_tails/inputShapes";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
    TIter next(fIn -> GetListOfKeys()); 
    TKey *key; 
    double centrBinImpactParam = 0.0;
    for (int iter = 0; iter < 3; iter ++) {
        cout << "Entro nell'iter" << endl;
        TH1D *histMassJpsiGen = new TH1D("histMassJpsiGen", " ; #it{m}_{#mu#mu} (GeV/c^{2})", 240, 2, 5); SetHistogram(histMassJpsiGen, kRed+1);
        TH1D *histMassJpsiGenFromRec = new TH1D("histMassJpsiGenFromRec", " ; #it{m}_{#mu#mu} (GeV/c^{2})", 240, 2, 5); SetHistogram(histMassJpsiGenFromRec, kMagenta);
        TH1D *histMassJpsiRec = new TH1D("histMassJpsiRec", " ; #it{m}_{#mu#mu} (GeV/c^{2})", 240, 2, 5); SetHistogram(histMassJpsiRec, kBlue+1);

        TH1D *histPtJpsiGen = new TH1D("histPtJpsiGen", " ; #it{p}_{T} (GeV/c)", 300, 0, 30); SetHistogram(histPtJpsiGen, kRed+1);
        TH1D *histPtJpsiRec = new TH1D("histPtJpsiRec", " ; #it{p}_{T} (GeV/c)", 300, 0, 30); SetHistogram(histPtJpsiRec, kBlue+1);

        TH1D *histRapJpsiGen = new TH1D("histRapJpsiGen", " ; #it{y}", 150, 2.5, 4); SetHistogram(histRapJpsiGen, kRed+1);
        TH1D *histRapJpsiRec = new TH1D("histRapJpsiRec", " ; #it{y}", 150, 2.5, 4); SetHistogram(histRapJpsiRec, kBlue+1);

        TH1D *histCentrJpsiGen = new TH1D("histCentrJpsiGen", " ; Centr (%)", 20, 0, 20); SetHistogram(histCentrJpsiGen, kRed+1);
        TH1D *histCentrJpsiRec = new TH1D("histCentrJpsiRec", " ; Centr (%)", 10, 0, 100); SetHistogram(histCentrJpsiRec, kBlue+1);

        TH2F *histPtRapRec = new TH2F("histPtRapRec", "p_T vs rapidity", 300, 0, 30, 150, 2.5, 4);
        TH2F *histPtRapGen = new TH2F("histPtRapGen", "p_T vs rapidity", 300, 0, 30, 150, 2.5, 4);
        TH3F *histCentrPtRapGen = new TH3F("histCentrPtRapGen", "centr vs p_T vs rapidity", 93, 0, 100, 300, 0, 30, 150, 2.5, 4);
        TH1F *hPtDimuonRec = new TH1F("hPtDimuonRec", "Pt Dimuon", 300, 0, 30);
        TH1F *hPtDimuonGen = new TH1F("hPtDimuonGen", "Pt Dimuon", 300, 0, 30);
        TH1F *hRapDimuonRec = new TH1F("hRapDimuonRec", "Rap Dimuon", 150, 2.5, 4);
        TH1F *hRapDimuonGen = new TH1F("hRapDimuonGen", "Rap Dimuon", 150, 2.5, 4);

        TH1D *histNormResoPx1 = new TH1D("histNormResoPx1", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); histNormResoPx1 -> SetLineColor(kBlack);
        TH1D *histNormResoPy1 = new TH1D("histNormResoPy1", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); histNormResoPy1 -> SetLineColor(kBlack);
        TH1D *histNormResoPz1 = new TH1D("histNormResoPz1", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); histNormResoPz1 -> SetLineColor(kBlack);
        TH1D *histNormResoPx2 = new TH1D("histNormResoPx2", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); histNormResoPx2 -> SetLineColor(kBlack);
        TH1D *histNormResoPy2 = new TH1D("histNormResoPy2", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); histNormResoPy2 -> SetLineColor(kBlack);
        TH1D *histNormResoPz2 = new TH1D("histNormResoPz2", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); histNormResoPz2 -> SetLineColor(kBlack);

        TH2D *histRapJpsiAll = new TH2D("histRapJpsiAll", " ; #it{y}_{rec} ; #it{y}_{gen}", 150,  2.5, 4, 150,  2.5, 4);
        TH1D *histCentrJpsiGenRebin = (TH1D*) histCentrJpsiGen -> Rebin(nBinsCentr, "histCentrJpsiGenRebin", centrBins);
        TH2D *histCentrJpsiAll = new TH2D("histCentrJpsiAll", " ; #it{Centr}_{rec} (%); #it{b}_{gen} (fm)", nBinsCentr,  0, 100, 20,  0, 20);
        TH3D *histCentrPtRapJpsiGen = new TH3D("histCentrPtRapJpsiGen", " ; #it{Centr} (%); #it{p}_{T} (GeV/c); #it{y}", nBinsPt,  0, 12, nBinsCentr,  0, 90, nBinsRap, 2.5, 4);
        TH1F* histCentrJpsiGenFin = new TH1F("histCentrJpsiGenFin", "Centrality of the generated", nBinsCentr, 0, 90);
        TH1D *histCentrJpsiRecRebin = new TH1D("histCentrJpsiRecRebin", " ; Centr (%)", 4, 0, 90); SetHistogram(histCentrJpsiRecRebin, kBlue+1);

        TH1D* hCentrality = new TH1D("hCentrality", "Generated centrality", 4, 0, 4);
        hCentrality->GetXaxis()->SetBinLabel(1, "0-20%");
        hCentrality->GetXaxis()->SetBinLabel(2, "20-40%");
        hCentrality->GetXaxis()->SetBinLabel(3, "40-60%");
        hCentrality->GetXaxis()->SetBinLabel(4, "60-90%");

        while ((key = (TKey*) next())) { 
            cout << "Sono nel while" << endl;
            TString dirName = key -> GetName();
            if (!dirName.Contains("DF_")) {
                continue;
            }

            TTree *treeGen = (TTree*) fIn -> Get(Form("%s/O2rtdilmtreegen", dirName.Data()));
            treeGen -> SetBranchAddress("fMcDecision", &fMcDecision);
            treeGen -> SetBranchAddress("fImpactParameter", &fImpactParameter);
            treeGen -> SetBranchAddress("fPtMC1", &fPtMC1);
            treeGen -> SetBranchAddress("fEtaMC1", &fEtaMC1);
            treeGen -> SetBranchAddress("fPhiMC1", &fPhiMC1);
            treeGen -> SetBranchAddress("fPtMC2", &fPtMC2);
            treeGen -> SetBranchAddress("fEtaMC2", &fEtaMC2);
            treeGen -> SetBranchAddress("fPhiMC2", &fPhiMC2);

            cout << "TTree Gen trovato" << endl;

            std::map<int, double> mapRapidityGen;
            std::map<int, double> mapEtaGen;

            float sumB = 0.0;

            for (int iEntry = 0;iEntry < treeGen -> GetEntries();iEntry++) {
                treeGen -> GetEntry(iEntry);
                if (fMcDecision < 1) continue;

                if (fPtMC2 == -999 && fPhiMC2 == -999 && fPhiMC2 == -999) {
                    ROOT::Math::PtEtaPhiMVector vecJpsiGen(fPtMC1, fEtaMC1, fPhiMC1, massJpsi);
                    if (TMath::Abs(vecJpsiGen.Rapidity()) > 4 || TMath::Abs(vecJpsiGen.Rapidity()) < 2.5) continue;
                    if (vecJpsiGen.Pt() > 30) continue;

                    double weightTotRap = 1;
                    double weightTotPt = 1;
                    double weightTot = 1;
                    int centralityBin = -1;

                    //Double_t centrImpParam = (fImpactParameter/20.0)*100;
                    /* if(fImpactParameter >= 0 && fImpactParameter <= 6.98){
                        sumB += fImpactParameter;
                        //cout << "sumB: " << sumB << " fImpactParameter: " << fImpactParameter << endl;
                    } */
                    //Double_t checkB = (fImpactParameter/20.0)*100;
                    impactParameterMap[iEntry] = fImpactParameter;
                    /*if(centrImpParam < 20.001){
                        std::string histCentrPtRapGenName = "histCentrPtRapGen_" + centrBins[0] + "_" + centrBins[1];
                        std::string histCentrPtRapGenTitle = "Centrality " + centrBins[0] + " - " + centrBins[1] + " vs p_T vs rapidity";
                    }
                    if(centrImpParam > 20.0 && centrImpParam < 40.001){
                        std::string histCentrPtRapGenName = "histCentrPtRapGen_" + centrBins[1] + "_" + centrBins[2];
                        std::string histCentrPtRapGenTitle = "Centrality " + centrBins[1] + " - " + centrBins[2] + " vs p_T vs rapidity";
                    }
                    if(centrImpParam > 40.0 && centrImpParam < 60.001){
                        std::string histCentrPtRapGenName = "histCentrPtRapGen_" + centrBins[2] + "_" + centrBins[3];
                        std::string histCentrPtRapGenTitle = "Centrality " + centrBins[2] + " - " + centrBins[3] + " vs p_T vs rapidity";
                    }
                    else{
                        std::string histCentrPtRapGenName = "histCentrPtRapGen_" + centrBins[3] + "_" + centrBins[4];
                        std::string histCentrPtRapGenTitle = "Centrality " + centrBins[3] + " - " + centrBins[4] + " vs p_T vs rapidity";
                    }
                    TH3F* histCentrPtRapGen = new TH3F(histCentrPtRapGenName.c_str(), histCentrPtRapGenTitle.c_str(), 93, 0, 100,300, 0, 30, 150, 2.5, 4); */
                    if(iter == 0){
                        histPtRapGen->Fill(vecJpsiGen.Pt(), -vecJpsiGen.Rapidity());
                        histCentrPtRapGen->Fill(fImpactParameter, vecJpsiGen.Pt(), -vecJpsiGen.Rapidity());
                        hPtDimuonGen->Fill(vecJpsiGen.Pt());
                        hRapDimuonGen->Fill(-vecJpsiGen.Rapidity());
                        //histCentrPtRapGenList.push_back(histCentrPtRapGen);
                    }
                    if(iter == 1){
                        if(fImpactParameter < 5.5){
                            weightTotRap = weightFuncArrayRapGen_0_20[0]->Eval(-vecJpsiGen.Rapidity());
                            weightTotPt = weightFuncArrayPtGen_0_20[0]->Eval(vecJpsiGen.Pt());
                            centralityBin = 1; 
                        }
                        if(fImpactParameter > 5.4 && fImpactParameter < 8.5){
                            weightTotRap = weightFuncArrayRapGen_20_40[0]->Eval(-vecJpsiGen.Rapidity());
                            weightTotPt = weightFuncArrayPtGen_20_40[0]->Eval(vecJpsiGen.Pt());
                            centralityBin = 2; 
                        }
                        if(fImpactParameter > 8.4 && fImpactParameter < 10.5){
                            weightTotRap = weightFuncArrayRapGen_40_60[0]->Eval(-vecJpsiGen.Rapidity());
                            weightTotPt = weightFuncArrayPtGen_40_60[0]->Eval(vecJpsiGen.Pt());
                            centralityBin = 3; 
                        }
                        if(fImpactParameter > 10.4 && fImpactParameter <= 13.5){
                            weightTotRap = weightFuncArrayRapGen_60_90[0]->Eval(-vecJpsiGen.Rapidity());
                            weightTotPt = weightFuncArrayPtGen_60_90[0]->Eval(vecJpsiGen.Pt());
                            centralityBin = 4; 
                        }
                        if(fImpactParameter > 13.5) continue;
                        weightTot = weightTotPt * weightTotRap;
                        hRapDimuonGen->Fill(-vecJpsiGen.Rapidity(), weightTot);
                        histPtRapGen->Fill(vecJpsiGen.Pt(), -vecJpsiGen.Rapidity(), weightTot);
                        hPtDimuonGen->Fill(vecJpsiGen.Pt(), weightTot);
                        histCentrPtRapGen->Fill(fImpactParameter, vecJpsiGen.Pt(), -vecJpsiGen.Rapidity(), weightTot);
                        //histCentrPtRapGenList.push_back(histCentrPtRapGen);
                    }
                    else{
                        if(fImpactParameter < 5.5){
                            for (int iArr = 0; iArr < weightFuncArrayRapGen_0_20.size(); iArr++) {
                                TF1* func = weightFuncArrayRapGen_0_20[iArr];
                                weightTotRap *= func->Eval(-vecJpsiGen.Rapidity());  
                                weightTotPt *= func->Eval(vecJpsiGen.Pt());  
                                weightTot *= weightTotPt * weightTotRap;
                            }
                            centralityBin = 1; 
                        }
                        if(fImpactParameter > 5.4 && fImpactParameter < 8.5){
                            for (int iArr = 0; iArr < weightFuncArrayRapGen_20_40.size(); iArr++) {
                                TF1* func = weightFuncArrayRapGen_20_40[iArr];
                                weightTotRap *= func->Eval(-vecJpsiGen.Rapidity());  
                                weightTotPt *= func->Eval(vecJpsiGen.Pt());  
                                weightTot *= weightTotPt * weightTotRap;
                            }
                            centralityBin = 2; 
                        }
                        if(fImpactParameter > 8.4 && fImpactParameter < 10.5){
                            for (int iArr = 0; iArr < weightFuncArrayRapGen_40_60.size(); iArr++) {
                                TF1* func = weightFuncArrayRapGen_40_60[iArr];
                                weightTotRap *= func->Eval(-vecJpsiGen.Rapidity());  
                                weightTotPt *= func->Eval(vecJpsiGen.Pt());  
                                weightTot *= weightTotPt * weightTotRap;
                            }
                            centralityBin = 3; 
                        }
                        if(fImpactParameter > 10.4 && fImpactParameter <= 13.5){
                            for (int iArr = 0; iArr < weightFuncArrayRapGen_60_90.size(); iArr++) {
                                TF1* func = weightFuncArrayRapGen_60_90[iArr];
                                weightTotRap *= func->Eval(-vecJpsiGen.Rapidity());  
                                weightTotPt *= func->Eval(vecJpsiGen.Pt());  
                                weightTot *= weightTotPt * weightTotRap;
                            }
                            centralityBin = 4; 
                        }
                        if(fImpactParameter > 13.5) continue;
                        hRapDimuonGen->Fill(-vecJpsiGen.Rapidity(), weightTot);
                        histPtRapGen->Fill(vecJpsiGen.Pt(), -vecJpsiGen.Rapidity(), weightTot);
                        hPtDimuonGen->Fill(vecJpsiGen.Pt(), weightTot);
                        histCentrPtRapGen->Fill(fImpactParameter, vecJpsiGen.Pt(), -vecJpsiGen.Rapidity(), weightTot);
                        //histCentrPtRapGenList.push_back(histCentrPtRapGen);
                        if (centralityBin > 0) {
                            hCentrality->Fill(centralityBin - 0.5, weightTot);
                        }
                    }

                    std::tuple<double, double, double> key = {fImpactParameter, vecJpsiGen.Pt(), -vecJpsiGen.Rapidity()};
                    weightMap[key] = weightTot;

                    histMassJpsiGen -> Fill(vecJpsiGen.M());
                    mapRapidityGen[iEntry] = -vecJpsiGen.Rapidity();
                    mapEtaGen[iEntry] = -vecJpsiGen.Eta();
                    histCentrJpsiGen -> Fill(fImpactParameter);
                } else {
                    ROOT::Math::PtEtaPhiMVector vecMuGen1(fPtMC1, fEtaMC1, fPhiMC1, massMu);
                    ROOT::Math::PtEtaPhiMVector vecMuGen2(fPtMC2, fEtaMC2, fPhiMC2, massMu);
                    auto vecJpsiGen = vecMuGen1 + vecMuGen2;
                    if (TMath::Abs(vecJpsiGen.Rapidity()) > 4 || TMath::Abs(vecJpsiGen.Rapidity()) < 2.5) continue;
                    histMassJpsiGen -> Fill(vecJpsiGen.M());
                }
            }
            /* sumB = histCentrJpsiGen -> Integral(0.0,6.98);
            float checkC = sumB/3232444.0;
            cout << "sumB: " << sumB << endl;
            cout << "0-20 CHECK B: " << checkC << endl;
            sumB = 0; */

            TTree *treeRec = (TTree*) fIn -> Get(Form("%s/O2rtdilmtreerec", dirName.Data()));
            treeRec -> SetBranchAddress("fMcDecision", &fMcDecision);
            treeRec -> SetBranchAddress("fMass", &fMass);
            treeRec -> SetBranchAddress("fPt", &fPt);
            treeRec -> SetBranchAddress("fEta", &fEta);
            treeRec -> SetBranchAddress("fPhi", &fPhi);
            treeRec -> SetBranchAddress("fCentFT0C", &fCentFT0C);
            treeRec -> SetBranchAddress("fPtMC1", &fPtMC1);
            treeRec -> SetBranchAddress("fEtaMC1", &fEtaMC1);
            treeRec -> SetBranchAddress("fPhiMC1", &fPhiMC1);
            treeRec -> SetBranchAddress("fPtMC2", &fPtMC2);
            treeRec -> SetBranchAddress("fEtaMC2", &fEtaMC2);
            treeRec -> SetBranchAddress("fPhiMC2", &fPhiMC2);
            treeRec -> SetBranchAddress("fPt1", &fPt1);
            treeRec -> SetBranchAddress("fEta1", &fEta1);
            treeRec -> SetBranchAddress("fPhi1", &fPhi1);
            treeRec -> SetBranchAddress("fPt2", &fPt2);
            treeRec -> SetBranchAddress("fEta2", &fEta2);
            treeRec -> SetBranchAddress("fPhi2", &fPhi2);

            cout << "TTree Rec trovato" << endl;
            const double bMax = 20.0;
            TH1D *histCentralitaCalcolata = new TH1D("histCentralitaCalcolata", "; Centralità ricostruita (%); Centralità calcolata (%)", nBinsCentr, 0, 90);


            for (int iEntry = 0;iEntry < treeRec -> GetEntries();iEntry++) {
                treeRec -> GetEntry(iEntry);
                if (fMcDecision < 1) continue;
                ROOT::Math::PtEtaPhiMVector vecMuGen1(fPtMC1, fEtaMC1, fPhiMC1, massMu);
                ROOT::Math::PtEtaPhiMVector vecMuGen2(fPtMC2, fEtaMC2, fPhiMC2, massMu);
                ROOT::Math::PtEtaPhiMVector vecMuRec1(fPt1, fEta1, fPhi1, massMu);
                ROOT::Math::PtEtaPhiMVector vecMuRec2(fPt2, fEta2, fPhi2, massMu);

                auto vecJpsiGen = vecMuGen1 + vecMuGen2;
                ROOT::Math::PtEtaPhiMVector vecJpsiRec(fPt, fEta, fPhi, fMass);

                if (TMath::Abs(vecJpsiRec.Rapidity()) > 4 || TMath::Abs(vecJpsiRec.Rapidity()) < 2.5) continue;
                if (fPt1 < ptCut || fPt2 < ptCut) continue;
                if (fPt > 30) continue;

                double pxGen1 = vecMuGen1.Px();
                double pyGen1 = vecMuGen1.Py();
                double pzGen1 = vecMuGen1.Pz();
                double pxGen2 = vecMuGen2.Px();
                double pyGen2 = vecMuGen2.Py();
                double pzGen2 = vecMuGen2.Pz();

                double pxRec1 = vecMuRec1.Px();
                double pyRec1 = vecMuRec1.Py();
                double pzRec1 = vecMuRec1.Pz();
                double pxRec2 = vecMuRec2.Px();
                double pyRec2 = vecMuRec2.Py();
                double pzRec2 = vecMuRec2.Pz();

                histNormResoPx1 -> Fill((pxRec1 - pxGen1) / pxGen1);
                histNormResoPy1 -> Fill((pyRec1 - pyGen1) / pyGen1);
                histNormResoPz1 -> Fill((pzRec1 - pzGen1) / pzGen1);

                histNormResoPx2 -> Fill((pxRec2 - pxGen2) / pxGen2);
                histNormResoPy2 -> Fill((pyRec2 - pyGen2) / pyGen2);
                histNormResoPz2 -> Fill((pzRec2 - pzGen2) / pzGen2);

                histMassJpsiGenFromRec -> Fill(vecJpsiGen.M());
                histMassJpsiRec -> Fill(vecJpsiRec.M());
                //histPtJpsiRec -> Fill(vecJpsiRec.Pt());
                //histRapJpsiRec -> Fill(-vecJpsiRec.Rapidity());
                histCentrJpsiRec -> Fill(fCentFT0C);

                double weightTot = 1;
                std::tuple<double, double, double> key = {fCentFT0C, vecJpsiRec.Pt(), -vecJpsiRec.Rapidity()};
                if (weightMap.find(key) != weightMap.end()) {
                    weightTot = weightMap[key];
                }
                if(iter == 0){
                    histPtRapRec->Fill(vecJpsiRec.Pt(), -vecJpsiRec.Rapidity());
                    hPtDimuonRec->Fill(vecJpsiRec.Pt());
                    hRapDimuonRec->Fill(-vecJpsiRec.Rapidity());
                }
                else{
                    hRapDimuonRec->Fill(-vecJpsiRec.Rapidity(), weightTot);
                    histPtRapRec->Fill(vecJpsiRec.Pt(), -vecJpsiRec.Rapidity(), weightTot);
                    hPtDimuonRec->Fill(vecJpsiRec.Pt(), weightTot);
                }

                if (mapRapidityGen.find(iEntry) != mapRapidityGen.end()) {
                    histRapJpsiAll->Fill(-vecJpsiRec.Rapidity(), mapRapidityGen[iEntry]);
                }
                if (mapEtaGen.find(iEntry) != mapEtaGen.end()) {
                    histRapJpsiAll->Fill(-vecJpsiRec.Eta(), mapEtaGen[iEntry]);
                }
                double impactParamGen = -1; 
                if (impactParameterMap.find(iEntry) != impactParameterMap.end()) {
                    impactParamGen = impactParameterMap[iEntry];
                }
                histCentrJpsiAll->Fill(fCentFT0C, impactParamGen);
            }
            double minVarBin = histCentrJpsiRec -> GetXaxis()-> FindBin(0.0);
            double maxVarBin = histCentrJpsiRec -> GetXaxis()-> FindBin(90 - 0.01);
            histCentrJpsiRec -> GetXaxis() -> SetRange(minVarBin,maxVarBin);
            histCentrJpsiRecRebin = (TH1D*) histCentrJpsiRec -> Rebin(nBinsCentr, "histCentrJpsiRecRebin", centrBins);
            
            std::vector<int> binsImpactParam;
            
            double minVarBin2 = histCentrJpsiAll -> GetXaxis()-> FindBin(0.0);
            double maxVarBin2 = histCentrJpsiAll -> GetXaxis()-> FindBin(90 - 0.01);
            histCentrJpsiAll -> GetXaxis() -> SetRange(minVarBin2,maxVarBin2);
            for (int iCentr = 1; iCentr <= nBinsCentr; iCentr++) {
                int bin = histCentrJpsiAll->GetXaxis()->FindBin(centrBins[iCentr]);
                binsImpactParam.push_back(bin);
                std::cout << "Centralità " << centrBins[iCentr] << "% -> Bin: " << bin << std::endl;
            }

            int nBinsGen = histCentrJpsiGen->GetNbinsX();
            int indexCentr = 0;
            for (int binGen = 1; binGen <= nBinsGen; binGen++) {
                Double_t impParam = histCentrJpsiGen->GetXaxis()->GetBinCenter(binGen); 
                Double_t centrImpParam = (impParam/20.0)*(impParam/20.0)*100;
                int centrBin = 0;
                if (centrImpParam < centrBins[indexCentr + 1]) {
                    centrBin = indexCentr + 1; 
                    double entriesCentr = histCentrJpsiGen -> GetBinContent(centrBin);
                    histCentrJpsiGenFin->SetBinContent(centrBin, entriesCentr);
                }
                else {
                    indexCentr++;
                }
            }

            /* TCanvas* cProva = new TCanvas("cProva", "Generated Centrality Distribution", 800, 600);
            histCentrJpsiGenFin->Draw(); */
            
            TFile *fOut = new TFile(Form("%s/centrRec.root", pathToFileOut.c_str()), "RECREATE");
            histCentrJpsiRec -> Write("histCentrJpsiRec");
            histCentrJpsiRecRebin -> Write("histCentrJpsiRecRebin");
            hCentrality -> Write("histCentrJpsiGenFin");
            histCentrPtRapGen -> Write("histCentrPtRapGen");
            histCentrJpsiAll -> Write("histCentrJpsiAll");
            fOut->Close();
        }
    
        //fIn -> Close();

        //--------------------------------------------------- Save all plots in allPlots.root ---------------------------------------------------//

        TFile *fOut = new TFile(Form("%s/allPlots_newTree.root", pathToFile.c_str()), "RECREATE");
        histMassJpsiGenFromRec -> Write("hMassGen");
        histMassJpsiRec -> Write("hMassRec");
        hPtDimuonGen-> Write("hPtGen");
        hPtDimuonRec -> Write("hPtRec");
        hRapDimuonGen-> Write("hRapDimuonGen");
        hRapDimuonRec -> Write("hRapDimuonRec");
        histPtRapGen->Write("histPtRapGen");
        histPtRapRec->Write("histPtRapRec");
        fOut->Close();

        cout << "Scrivo nel file fOut" << endl;

        //fIn->Close();
        
        //------------------------------------------------- Generated MC events -------------------------------------------------//
        //string pathInAxe = "/Users/saragaretti/dq_fitter_tails/inputShapes/plots_newTree";
        string pathOutAxe = "/Users/saragaretti/dq_fitter_tails/inputShapes/plotsAxe_newTree";

        //TFile *fInAxe = new TFile(Form("%s/allPlots_newTree.root",pathInAxe.c_str()),"READ");

        //----------------------------------------------- Generated J/Psi events ------------------------------------------------//
        //TH2D *histGenJPsi = (TH2D*) fInAxe -> Get("histPtRapGen");
        TH2D *histGenJPsi = (TH2D*) histPtRapGen -> Clone();
        //TH2D *histGenJPsi = (TH2D*) fInAxe -> Get("histPtRapGen");
        TH2D *histGenJpsiPtRapClone = (TH2D*) CutTH2X(histGenJPsi, "histGenJPsi",0., 20.);
        TH2D *histGenJpsiPtRapClone2 = (TH2D*) CutTH2Y(histGenJPsi, "histGenJPsi", 2.5, 4);
        TH1D *histGenJpsiPt = (TH1D*) histGenJpsiPtRapClone2 -> ProjectionX("histGenJpsiPt");
        TH1D *histGenJpsiRap = (TH1D*) histGenJpsiPtRapClone2 -> ProjectionY("histGenJpsiRap");
        TH1D *histGenJpsiCentr = (TH1D*) histCentrJpsiGenFin -> Clone("histGenJpsiCentr");
        
        //-------------------------------------------- Reconstructed J/Psi events ----------------------------------------------//

        /* TH2D *histRecJPsi = (TH2D*) fInAxe -> Get("histPtRapRec"); */
        TH2D *histRecJPsi = (TH2D*) histPtRapRec -> Clone();
        TH2D *histRecJpsiPtRapClone = (TH2D*) CutTH2X(histRecJPsi, "histRecJPsi",0., 20.);
        TH2D *histRecJpsiPtRapClone2 = (TH2D*) CutTH2Y(histRecJPsi, "histRecJPsi", 2.5, 4);
        TH1D *histRecJpsiPt = (TH1D*) histRecJpsiPtRapClone2 -> ProjectionX("histRecJpsiPt");
        TH1D *histRecJpsiRap = (TH1D*) histRecJpsiPtRapClone2 -> ProjectionY("histRecJpsiRap");
        TH1D *histRecJpsiCentr = (TH1D*) histCentrJpsiRecRebin -> Clone("histRecJpsiCentr");

        //------------------------------------------------ J/Psi Axe pT and rap -------------------------------------------------//

        TH1D *histGenJpsiPtRebin = (TH1D*) histGenJpsiPt -> Rebin(nBinsPt, "histGenJpsiPtRebin", ptBinsRun3); 
        TH1D *histRecJpsiPtRebin = (TH1D*) histRecJpsiPt -> Rebin(nBinsPt, "histRecJpsiPtRebin", ptBinsRun3); 

        TH1D *histGenJpsiRapRebin = (TH1D*) histGenJpsiRap -> Rebin(nBinsRap, "histGenJpsiRapRebin", rapBins); 
        TH1D *histRecJpsiRapRebin = (TH1D*) histRecJpsiRap -> Rebin(nBinsRap, "histRecJpsiRapRebin", rapBins);

        //histCentrPtRapJpsiGen -> Fill(histGenJpsiCentr, histGenJpsiPtRebin, histGenJpsiRapRebin);

        for (int icentr = 1; icentr <= histGenJpsiCentr->GetNbinsX(); icentr++) {
            double centr = histGenJpsiCentr->GetXaxis()->GetBinCenter(icentr);
            double weightCentr = histGenJpsiCentr->GetBinContent(icentr);

            for (int ipt = 1; ipt <= histGenJpsiPtRebin->GetNbinsX(); ipt++) {
                double pt = histGenJpsiPtRebin->GetXaxis()->GetBinCenter(ipt);
                double weightPt = histGenJpsiPtRebin->GetBinContent(ipt);

                for (int irap = 1; irap <= histGenJpsiRapRebin->GetNbinsX(); irap++) {
                    double rap = histGenJpsiRapRebin->GetXaxis()->GetBinCenter(irap);
                    double weightRap = histGenJpsiRapRebin->GetBinContent(irap);

                    double weightTotal = weightCentr * weightPt * weightRap;

                    if (weightTotal > 0) {
                        histCentrPtRapJpsiGen->Fill(centr, pt, rap, weightTotal);
                    }
                }
            }
        }

        double errRecJpsi, errGenJpsi;

        double intRecJpsi = histRecJpsiPtRapClone2->IntegralAndError(1, histRecJpsiPtRapClone2->GetNbinsX(), 1, histRecJpsiPtRapClone2->GetNbinsY(), errRecJpsi);
        double intGenJpsi = histGenJpsiPtRapClone2->IntegralAndError(1, histGenJpsiPtRapClone2->GetNbinsX(), 1, histGenJpsiPtRapClone2->GetNbinsY(), errGenJpsi);

        double AxeJpsi = intRecJpsi / intGenJpsi;
        double errAxeJpsi = AxeJpsi * TMath::Sqrt(TMath::Power(errRecJpsi / intRecJpsi, 2) + TMath::Power(errGenJpsi / intGenJpsi, 2));

        std::cout << "AxeJpsi: " << AxeJpsi << " ± " << errAxeJpsi << std::endl;

        std::ofstream outFile(Form("%s/acceptance_efficiency.txt", pathOutAxe.c_str()));
        if (outFile.is_open()) {
            outFile << AxeJpsi << "  " << errAxeJpsi << std::endl;
            outFile.close();
            std::cout << "Risultato salvato in 'acceptance_efficiency_result.txt'." << std::endl;
        } else {
            std::cerr << "Errore nell'aprire il file per la scrittura." << std::endl;
        }


        TH2D *histAxeJpsiPtRap = (TH2D*) histRecJpsiPtRapClone2 -> Clone("histAxeJpsiPtRap");
        histAxeJpsiPtRap -> Divide(histGenJpsiPtRapClone2);
        //cout << "bins: " << histGenJpsiPt -> GetNbinsX() << endl;

        //------------------------------------------------ J/Psi Axe centr bins -------------------------------------------------//

        TH1D *histAxeJpsiCentr = (TH1D*) histRecJpsiCentr -> Clone("histAxeJpsiCentr");
        histAxeJpsiCentr -> Divide(histGenJpsiCentr);

        canvasAxeCentr -> cd();
        if(iter == 0){
            histAxeJpsiCentr -> Draw();
        }
        else{
            histAxeJpsiCentr -> Draw("SAME");
        }


        //------------------------------------------------ J/Psi Axe pT cuts -------------------------------------------------//

        TH1D *histAxeJpsiPt = new TH1D("histAxeJpsiPt", "", nBinsPt, ptBinsRun3);
        histAxeJpsiPt -> Divide(histRecJpsiPtRebin, histGenJpsiPtRebin, 1, 1, "B");
        canvasAxePt -> cd();
        if(iter == 0){
            histAxeJpsiPt -> Draw();
        }
        else{
            histAxeJpsiPt -> Draw("SAME");
        }  

        //------------------------------------------------- J/Psi Axe y cuts -------------------------------------------------//

        TH1D *histAxeJpsiRap = new TH1D("histAxeJpsiRap", "", nBinsRap, rapBins);
        histAxeJpsiRap -> Divide(histRecJpsiRapRebin, histGenJpsiRapRebin, 1, 1, "B");

        canvasAxeRap -> cd();
        if(iter == 0){
            histAxeJpsiRap -> Draw();
        }
        else{
            histAxeJpsiRap -> Draw("SAME");
        }

        //--------------------------------Take the results from the fits ----------------------------------------------------//
        string pathInData = "/Users/saragaretti/dq_fitter_tails/analysis/muonLowPt210SigmaPDCA/output_chi2/inputShapes/";

        std::vector<TH1F*> hist_Centr_Arr_Pt;
        std::vector<TH1F*> hist_Centr_Arr_Rap;
        int index_centr = 0;
        for (size_t iCentr = 0; iCentr < nBinsCentr; ++iCentr) {
            std::ostringstream folderCentr;
            folderCentr << std::fixed << std::setprecision(1) << centrBins[iCentr] << "_" << centrBins[iCentr+1];
            std::string pathCentr = pathInData + "Centr_" + folderCentr.str();
            //std::string pathCentr = pathInData + "Centr_" + std::to_string(centrBins[iCentr]) + "_" + std::to_string(centrBins[iCentr+1]);
            TH1F* hist_Centr_pT = new TH1F(Form("hist_Centr_%.0f_%.0f_Pt", centrBins[iCentr], centrBins[iCentr+1]),Form("Results for Centrality %.0f-%.0f%%", centrBins[iCentr], centrBins[iCentr+1]),nBinsPt, ptBinsRun3);
            int index_hist_centr_pT = 0;
            for (size_t iPt = 0; iPt < nBinsPt; ++iPt) {
                TH1F* hist_res_pT = nullptr;
                std::ostringstream folderPt;
                folderPt << std::fixed << std::setprecision(1) << ptBinsRun3[iPt] << "_" << ptBinsRun3[iPt+1];
                std::string pathPt = pathCentr + "/Pt_" + folderPt.str();
                //std::string pathPt = pathCentr + "/Pt_" + std::to_string(ptBinsRun3[iPt]) + "_" + std::to_string(ptBinsRun3[iPt+1]);

                TFile *fInData = new TFile(Form("%s/output__Mass_BkgSub_CentrFT0C_%.0f_%.0f_Pt_%.0f_%.0f_ME_CB2_VWG__2.2_4.5.root",pathPt.c_str(), centrBins[iCentr], centrBins[iCentr+1], ptBinsRun3[iPt], ptBinsRun3[iPt+1]),"READ");
                //fInData -> ls();
                if (!fInData || fInData->IsZombie()) {
                    cout << "Error: impossible to open the data file" << endl;
                    return;
                }
                TH1F* hist_res_pT_Temp = (TH1F*) fInData->Get(Form("fit_results_CB2_VWG__2.2_4.5_Mass_BkgSub_CentrFT0C_%.0f_%.0f_Pt_%.0f_%.0f_ME",centrBins[iCentr], centrBins[iCentr+1], ptBinsRun3[iPt], ptBinsRun3[iPt+1]));
                if (!hist_res_pT_Temp) {
                    std::cout << "Error: " << Form("%s/output__Mass_BkgSub_CentrFT0C_%.0f_%.0f_Pt_%.0f_%.0f_ME_CB2_VWG__2.2_4.5.root",pathPt.c_str(), centrBins[iCentr], centrBins[iCentr+1], ptBinsRun3[iPt], ptBinsRun3[iPt+1]) << std::endl;
                    std::cout << "Error: Histogram " << Form("fit_results_CB2_VWG__2.2_4.5_Mass_BkgSub_CentrFT0C_%.0f_%.0f_Pt_%.0f_%.0f_ME",centrBins[iCentr], centrBins[iCentr+1], ptBinsRun3[iPt], ptBinsRun3[iPt+1]) << " not found for " << pathPt << std::endl;
                    fInData->Close();
                    continue;
                }
                hist_res_pT = (TH1F*)hist_res_pT_Temp -> Clone();
                hist_res_pT -> SetDirectory(0);
                int binToExtract = 10;
                double value = hist_res_pT->GetBinContent(binToExtract);
                cout << "Centr: " << iCentr << " Pt: " << iPt << " yield: " << value << endl;
                double error = hist_res_pT->GetBinError(binToExtract);

                double AxePt = histAxeJpsiPt->GetBinContent(index_hist_centr_pT + 1);
                cout << "Centr: " << iCentr << " Pt: " << iPt << " AxePt: " << AxePt << endl;

                double value_corr = value/AxePt;
                double error_corr = error/AxePt;

                double yield_div_pT = value_corr/(ptBinsRun3[iPt+1] - ptBinsRun3[iPt]);
                double error_yield_div_pT = error_corr/(ptBinsRun3[iPt+1] - ptBinsRun3[iPt]);

                hist_Centr_pT->SetBinContent(index_hist_centr_pT + 1, yield_div_pT);
                cout << "Centr: " << iCentr << " Pt: " << iPt << " yield_corr: " << yield_div_pT << endl;
                hist_Centr_pT->SetBinError(index_hist_centr_pT + 1, error_yield_div_pT);
                
                fInData->Close();
                index_hist_centr_pT ++;
            }
            hist_Centr_Arr_Pt.push_back(hist_Centr_pT);

            //---------------------------------------------- Rapidity --------------------------------------//

            TH1F* hist_Centr_Rap = new TH1F(Form("hist_Centr_%.0f_%.0f_Rap", centrBins[iCentr], centrBins[iCentr+1]),Form("Results for Centrality %.0f-%.0f%%", centrBins[iCentr], centrBins[iCentr+1]),nBinsRap, rapBins);
            int index_hist_centr_Rap = 0;
            for (size_t iRap = 0; iRap < nBinsRap; ++iRap) {
                TH1F* hist_res_pT = nullptr;
                std::ostringstream folderRap;
                folderRap << std::fixed << std::setprecision(2) << rapBins[iRap] << "_" << rapBins[iRap+1];
                std::string pathPt = pathCentr + "/Rap_" + folderRap.str();
                //std::string pathPt = pathCentr + "/Rap_" + std::to_string(rapBins[iRap]) + "_" + std::to_string(rapBins[iRap+1]);

                TFile *fInData = new TFile(Form("%s/output__Mass_BkgSub_CentrFT0C_%.0f_%.0f_Rap_%.2f_%.2f_ME_CB2_VWG__2.2_4.5.root",pathPt.c_str(), centrBins[iCentr], centrBins[iCentr+1], rapBins[iRap], rapBins[iRap+1]),"READ");
                //fInData -> ls();
                if (!fInData || fInData->IsZombie()) {
                    cout << "Error: impossible to open the data file" << endl;
                    return;
                }
                TH1F* hist_res_pT_Temp = (TH1F*) fInData->Get(Form("fit_results_CB2_VWG__2.2_4.5_Mass_BkgSub_CentrFT0C_%.0f_%.0f_Rap_%.2f_%.2f_ME",centrBins[iCentr], centrBins[iCentr+1], rapBins[iRap], rapBins[iRap+1]));
                if (!hist_res_pT_Temp) {
                    std::cout << "Error: " << Form("%s/output__Mass_BkgSub_CentrFT0C_%.0f_%.0f_Rap_%.2f_%.2f_ME_CB2_VWG__2.2_4.5.root",pathPt.c_str(), centrBins[iCentr], centrBins[iCentr+1], rapBins[iRap], rapBins[iRap+1]) << std::endl;
                    std::cout << "Error: Histogram " << Form("fit_results_CB2_VWG__2.2_4.5_Mass_BkgSub_CentrFT0C_%.0f_%.0f_Rap_%.2f_%.2f_ME",centrBins[iCentr], centrBins[iCentr+1], rapBins[iRap], rapBins[iRap+1]) << " not found for " << pathPt << std::endl;
                    fInData->Close();
                    continue;
                }
                hist_res_pT = (TH1F*)hist_res_pT_Temp -> Clone();
                hist_res_pT -> SetDirectory(0);
                int binToExtract = 10;
                double value = hist_res_pT->GetBinContent(binToExtract);
                cout << "Centr: " << iCentr << " Rap: " << iRap << " yield: " << value << endl;
                double error = hist_res_pT->GetBinError(binToExtract);

                double AxeRap = histAxeJpsiRap->GetBinContent(index_hist_centr_Rap + 1);
                cout << "Centr: " << iCentr << " Rap: " << iRap << " AxeRap: " << AxeRap << endl;

                double value_corr = value/AxeRap;
                double error_corr = error/AxeRap;

                double yield_div_rap = value_corr/(rapBins[iRap+1] - rapBins[iRap]);
                double error_yield_div_rap = error_corr/(rapBins[iRap+1] - rapBins[iRap]);

                hist_Centr_Rap->SetBinContent(index_hist_centr_Rap + 1, yield_div_rap);
                cout << "Centr: " << iCentr << " Rap: " << iRap << " yield_corr: " << yield_div_rap << endl;
                hist_Centr_Rap->SetBinError(index_hist_centr_Rap + 1, error_yield_div_rap);
                
                fInData->Close();
                index_hist_centr_Rap ++;
            }
            hist_Centr_Arr_Rap.push_back(hist_Centr_Rap);
            //index_centr ++;
        }
        TFile *fOutTest = new TFile(Form("%s/outputHistograms_newTree.root", pathInData.c_str()), "UPDATE");
        for (size_t centralityBin = 0; centralityBin < nBinsCentr; ++centralityBin) {
            //std::string histName = "hist_Centrality_" + std::to_string(centralityBin) + "_" + std::to_string(iter) + "_Pt";
            hist_Centr_Arr_Pt[centralityBin]->SetName(Form("hist_Centrality_%.0f_%.0f_Pt_iter_%i", centrBins[centralityBin], centrBins[centralityBin+1], iter));
            hist_Centr_Arr_Pt[centralityBin]->Write();
            TCanvas *canvasPt = new TCanvas(Form("canvas_Centrality_%.0f_%.0f_Pt", centrBins[centralityBin], centrBins[centralityBin+1]), "", 800, 600);
            hist_Centr_Arr_Pt[centralityBin]->Draw();

            TF1 *fitFunctionPt = PtJPsiPbPb5TeV_Func();
            //fitFunctionPt->SetParameters(1.00715e6, 3.50274, 1.93403, 3.96363);
            fitFunctionPt->SetLineColor(kRed); 
            hist_Centr_Arr_Pt[centralityBin]->Fit(fitFunctionPt, "R"); 
            fitFunctionPt->Draw("SAME");

            canvasPt->Write();
            /* std::string histNameRap = "hist_Centrality_" + std::to_string(centralityBin) + "_" + std::to_string(iter) + "_Rap";
            hist_Centr_Arr_Rap[centralityBin]->SetName(histNameRap.c_str()); */
            hist_Centr_Arr_Rap[centralityBin]->SetName(Form("hist_Centrality_%.0f_%.0f_Rap_iter_%i", centrBins[centralityBin], centrBins[centralityBin+1], iter));
            hist_Centr_Arr_Rap[centralityBin]->Write();

            TCanvas *canvasRap = new TCanvas(Form("canvas_Centrality_%.0f_%.0f_Rap", centrBins[centralityBin], centrBins[centralityBin+1]), "", 800, 600);
            hist_Centr_Arr_Rap[centralityBin]->Draw();

            TF1 *fitFunctionRap = PtJPsiPbPb5TeV_Func();
            //fitFunctionRap->SetParameters(1.00715e6, 3.50274, 1.93403, 3.96363);
            fitFunctionRap->SetLineColor(kRed);
            hist_Centr_Arr_Rap[centralityBin]->Fit(fitFunctionRap, "R");
            fitFunctionRap->Draw("SAME");

            canvasRap->Write();
        }

        fOutTest->Close();


        //------------------------------------------------------------- Fit the Pt generated distribution -------------------------------------------------------//

        TF1 *fitFunctionGenPt = PtJPsiPbPb5TeV_Func();

        histGenJpsiPtRebin->Fit(fitFunctionGenPt, "R");

        //---------------------------------------------------------- Fit the Rapidity generated distribution ----------------------------------------------------//
        std::vector<TH1F*> histListPtGen;
        std::vector<TH1F*> histListRapGen;
        for (int j = 0; j < nBinsCentr; j++) {
            int binMin = histCentrPtRapJpsiGen->GetXaxis()->FindBin(centrBins[j]);   
            int binMax = histCentrPtRapJpsiGen->GetXaxis()->FindBin(centrBins[j+1]) - 1; 

            // Proiezione dell'istogramma 3D su pT e rapidità per la centralità fissata
            TH1F* projPt = (TH1F*)histCentrPtRapJpsiGen->ProjectionY(Form("projPt_%d", j), binMin, binMax);
            TH1F* projRap = (TH1F*)histCentrPtRapJpsiGen->ProjectionZ(Form("projRap_%d", j), binMin, binMax);

            histListPtGen[j]->Add(projPt);
            histListRapGen[j]->Add(projRap);
        }
        TF1 *fitFunctionGenRap = PtJPsiPbPb5TeV_Func();
        histGenJpsiRapRebin->Fit(fitFunctionGenRap, "R");

        TF1 *fitFunctionGenPt_0_20 = PtJPsiPbPb5TeV_Func();
        histListPtGen[0]->Fit(fitFunctionGenPt_0_20, "R");
        weightFuncArrayPtGen_0_20.push_back(fitFunctionGenPt_0_20);

        TF1 *fitFunctionGenPt_20_40 = PtJPsiPbPb5TeV_Func();
        histListPtGen[1]->Fit(fitFunctionGenPt_20_40, "R");
        weightFuncArrayPtGen_20_40.push_back(fitFunctionGenPt_20_40);

        TF1 *fitFunctionGenPt_40_60 = PtJPsiPbPb5TeV_Func();
        histListPtGen[2]->Fit(fitFunctionGenPt_40_60, "R");
        weightFuncArrayPtGen_40_60.push_back(fitFunctionGenPt_40_60);

        TF1 *fitFunctionGenPt_60_90 = PtJPsiPbPb5TeV_Func();
        histListPtGen[3]->Fit(fitFunctionGenPt_60_90, "R");
        weightFuncArrayPtGen_60_90.push_back(fitFunctionGenPt_60_90);

        TF1 *fitFunctionGenRap_0_20 = PtJPsiPbPb5TeV_Func();
        histListRapGen[0]->Fit(fitFunctionGenRap_0_20, "R");
        weightFuncArrayRapGen_0_20.push_back(fitFunctionGenRap_0_20);

        TF1 *fitFunctionGenRap_20_40 = PtJPsiPbPb5TeV_Func();
        histListRapGen[1]->Fit(fitFunctionGenRap_20_40, "R");
        weightFuncArrayRapGen_20_40.push_back(fitFunctionGenRap_20_40);

        TF1 *fitFunctionGenRap_40_60 = PtJPsiPbPb5TeV_Func();
        histListRapGen[2]->Fit(fitFunctionGenRap_40_60, "R");
        weightFuncArrayRapGen_40_60.push_back(fitFunctionGenRap_40_60);

        TF1 *fitFunctionGenRap_60_90 = PtJPsiPbPb5TeV_Func();
        histListRapGen[3]->Fit(fitFunctionGenRap_60_90, "R");
        weightFuncArrayRapGen_60_90.push_back(fitFunctionGenRap_60_90);

        //---------------------------------------------------------- Fit the Pt corrected data distribution -----------------------------------------------------//
        TF1 *fitFunctionDataPt_0_20 = PtJPsiPbPb5TeV_Func();
        hist_Centr_Arr_Pt[0]->Fit(fitFunctionDataPt_0_20, "R");
        weightFuncArrayPt_0_20.push_back(fitFunctionDataPt_0_20);

        TF1 *fitFunctionDataPt_20_40 = PtJPsiPbPb5TeV_Func();
        hist_Centr_Arr_Pt[1]->Fit(fitFunctionDataPt_20_40, "R");
        weightFuncArrayPt_20_40.push_back(fitFunctionDataPt_20_40);

        TF1 *fitFunctionDataPt_40_60 = PtJPsiPbPb5TeV_Func();
        hist_Centr_Arr_Pt[2]->Fit(fitFunctionDataPt_40_60, "R");
        weightFuncArrayPt_40_60.push_back(fitFunctionDataPt_40_60);

        TF1 *fitFunctionDataPt_60_90 = PtJPsiPbPb5TeV_Func();
        hist_Centr_Arr_Pt[3]->Fit(fitFunctionDataPt_60_90, "R");
        weightFuncArrayPt_60_90.push_back(fitFunctionDataPt_60_90);

        TF1 *fitFunctionDataRap_0_20 = PtJPsiPbPb5TeV_Func();
        hist_Centr_Arr_Rap[0]->Fit(fitFunctionDataRap_0_20, "R");
        weightFuncArrayRap_0_20.push_back(fitFunctionDataRap_0_20);

        TF1 *fitFunctionDataRap_20_40 = PtJPsiPbPb5TeV_Func();
        hist_Centr_Arr_Rap[1]->Fit(fitFunctionDataRap_20_40, "R");
        weightFuncArrayRap_20_40.push_back(fitFunctionDataRap_20_40);

        TF1 *fitFunctionDataRap_40_60 = PtJPsiPbPb5TeV_Func();
        hist_Centr_Arr_Rap[2]->Fit(fitFunctionDataRap_40_60, "R");
        weightFuncArrayRap_40_60.push_back(fitFunctionDataRap_40_60);

        TF1 *fitFunctionDataRap_60_90 = PtJPsiPbPb5TeV_Func();
        hist_Centr_Arr_Rap[3]->Fit(fitFunctionDataRap_60_90, "R");
        weightFuncArrayRap_60_90.push_back(fitFunctionDataRap_60_90);


        if(iter == 0){
            //--------------------------------------------------- Create Canvases  ---------------------------------------------------//

            TCanvas *canvasNormReso = new TCanvas("canvasNormReso", "", 1800, 1200);
            canvasNormReso -> Divide(3, 2);

            canvasNormReso -> cd(1);
            gPad -> SetLogy(true);
            histNormResoPx1 -> Draw("EP");

            canvasNormReso -> cd(2);
            gPad -> SetLogy(true);
            histNormResoPy1 -> Draw("EP");

            canvasNormReso -> cd(3);
            gPad -> SetLogy(true);
            histNormResoPz1 -> Draw("EP");

            canvasNormReso -> cd(4);
            gPad -> SetLogy(true);
            histNormResoPx2 -> Draw("EP");

            canvasNormReso -> cd(5);
            gPad -> SetLogy(true);
            histNormResoPy2 -> Draw("EP");

            canvasNormReso -> cd(6);
            gPad -> SetLogy(true);
            histNormResoPz2 -> Draw("EP");

            TCanvas *canvasMass = new TCanvas("canvasMass", "", 800, 600);
            gPad -> SetLogy(true);
            histMassJpsiGen -> Draw("L");
            histMassJpsiGenFromRec -> Draw("EP SAME");
            histMassJpsiRec -> Draw("EP SAME");

            hPtDimuonGen-> Rebin(4);
            hPtDimuonRec -> Rebin(4);

            TH1D *histAxePtJpsi = new TH1D("histAxePtJpsi", " ; #it{p}_{T} (GeV/c)", hPtDimuonGen -> GetNbinsX(), 0, 30);
            histAxePtJpsi -> Divide(hPtDimuonRec, hPtDimuonGen, 1, 1, "B");
            SetHistogram(histAxePtJpsi, kBlack);

            TH1D *histAxeRapJpsi = new TH1D("histAxeRapJpsi", " ; #it{y}", hRapDimuonGen -> GetNbinsX(), 2.5, 4);
            histAxeRapJpsi -> Divide(hRapDimuonRec, hRapDimuonGen, 1, 1, "B");
            SetHistogram(histAxeRapJpsi, kBlack);

            TH1D *histAxeCentrJpsi = new TH1D("histAxeCentrJpsi", " ; #it{y}", histCentrJpsiGen -> GetNbinsX(), 0.0, 90.0);
            histAxeCentrJpsi -> Divide(histCentrJpsiRec, histCentrJpsiGen, 1, 1, "B");
            SetHistogram(histAxeCentrJpsi, kBlack);

            TCanvas *canvasPtRap = new TCanvas("canvasPtRap", "", 1200, 1200);
            canvasPtRap -> Divide(3, 2);

            canvasPtRap -> cd(1);
            gPad -> SetLogy(true);
            hPtDimuonGen-> GetYaxis() -> SetRangeUser(0.1, 1e5);
            hPtDimuonGen -> Draw("EP");
            hPtDimuonRec -> Draw("EP SAME");

            canvasPtRap -> cd(2);
            gPad -> SetLogy(true);
            hRapDimuonGen -> GetYaxis() -> SetRangeUser(0.1, 1e5);
            hRapDimuonGen -> Draw("EP");
            hRapDimuonRec -> Draw("EP SAME");

            canvasPtRap -> cd(3);
            histAxePtJpsi -> GetYaxis() -> SetRangeUser(0, 1.2);
            histAxePtJpsi -> Draw("EP");

            canvasPtRap -> cd(4);
            histAxeRapJpsi -> GetYaxis() -> SetRangeUser(0, 1.2);
            histAxeRapJpsi -> Draw("EP");

            canvasPtRap -> cd(5);
            histAxeJpsiCentr -> GetYaxis() -> SetRangeUser(0, 1.2);
            histAxeJpsiCentr -> Draw("EP");

            TCanvas *canvasRapTot = new TCanvas("canvasRapTot", "", 800, 600);
            histRapJpsiAll -> Draw("COLZ");

            TCanvas *canvasPtRapRec = new TCanvas("canvasPtRapRec", "PtRapRec", 800, 600);
            histPtRapRec->Draw("COLZ");
            canvasPtRapRec->SaveAs(Form("%s/PtRapRec.pdf", pathToFile.c_str()));

            TCanvas *canvasPtRapGen = new TCanvas("canvasPtRapGen", "PtRapGen", 800, 600);
            histPtRapGen->Draw("COLZ");
            canvasPtRapGen->SaveAs(Form("%s/histPtRapGen.pdf", pathToFile.c_str()));
        }
    }
    fIn -> Close();
    TFile *fOutAxe = new TFile("Axe_newTree.root", "RECREATE");
    canvasAxeCentr -> Write();
    canvasAxePt -> Write();
    canvasAxeRap -> Write();
}
////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTHnSparse(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2, double minVar3, double maxVar3) {
    // Pt bins
    double minVar1Bin = histSparse -> GetAxis(0) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(0) -> FindBin(maxVar1 - 0.01);
    // Rapidity bins
    double minVar2Bin = histSparse -> GetAxis(1) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(1) -> FindBin(maxVar2 - 0.01);
    // Centrality bins
    double minVar3Bin = histSparse -> GetAxis(2) -> FindBin(minVar3);
    double maxVar3Bin = histSparse -> GetAxis(2) -> FindBin(maxVar3 - 0.01);

    histSparse -> GetAxis(0) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(1) -> SetRange(minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar3Bin, maxVar3Bin);

    Printf("%3.2f - %3.2f §§ %3.2f - %3.2f §§ %3.2f - %3.2f \n", minVar1Bin, maxVar1Bin, minVar2Bin, maxVar2Bin, minVar3Bin, maxVar3Bin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%3.2f_%3.2f__%3.2f_%3.2f__%3.2f_%3.2f", minVar1, maxVar1, minVar2, maxVar2, minVar3, maxVar3));
    return histProj;
}
////////////////////////////////////////////////////////////////////////////////
THnSparseD* CutTHnSparse(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2, double minVar3, double maxVar3) {
    // Pt bins
    double minVar1Bin = histSparse -> GetAxis(0) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(0) -> FindBin(maxVar1 - 0.01);
    // Rapidity bins
    double minVar2Bin = histSparse -> GetAxis(1) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(1) -> FindBin(maxVar2 - 0.01);
    // Centrality bins
    double minVar3Bin = histSparse -> GetAxis(2) -> FindBin(minVar3);
    double maxVar3Bin = histSparse -> GetAxis(2) -> FindBin(maxVar3 - 0.01);

    histSparse -> GetAxis(0) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(1) -> SetRange(minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar3Bin, maxVar3Bin);

    Printf("%3.2f - %3.2f §§ %3.2f - %3.2f §§ %3.2f - %3.2f \n", minVar1Bin, maxVar1Bin, minVar2Bin, maxVar2Bin, minVar3Bin, maxVar3Bin);

    return histSparse;
}
////////////////////////////////////////////////////////////////////////////////
TH2D* CutTH2Y(TH2D *hist2D, string title, double minCentr, double maxCentr) {
    double minCentrBin = hist2D -> GetYaxis() -> FindBin(minCentr);
    double maxCentrBin = hist2D -> GetYaxis() -> FindBin(maxCentr - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f", minCentrBin, maxCentrBin);
    hist2D -> GetYaxis() -> SetRange(minCentrBin, maxCentrBin);

    TH2D *histProj = (TH2D*) hist2D -> Clone();
    return histProj;
}
////////////////////////////////////////////////////////////////////////////////
TH2D* CutTH2X(TH2D *hist2D, string title, double minCentr, double maxCentr) {
    double minCentrBin = hist2D -> GetXaxis() -> FindBin(minCentr);
    double maxCentrBin = hist2D -> GetXaxis() -> FindBin(maxCentr - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f", minCentrBin, maxCentrBin);
    hist2D -> GetXaxis() -> SetRange(minCentrBin, maxCentrBin);

    TH2D *histProj = (TH2D*) hist2D -> Clone();
    return histProj;
}
////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}