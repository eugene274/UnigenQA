
#include <cassert>

#include <TLorentzVector.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TString.h>
#include "UnigenQA.h"

using namespace TMath;
using namespace std;
using namespace qa;

UnigenQA::UnigenQA() {
}
UnigenQA::~UnigenQA() {
}

void UnigenQA::Init(TString filePath, TString treeName) {
  fChain = new TChain(treeName);
  fChain->Add(filePath);
  event_ = new UEvent;
  fChain->SetBranchAddress("event", &event_);

  // Get run description
  TString fileName;

  if (!filePath.Contains("*")) fileName = filePath;
  else {
    TString dirName = filePath.Remove(filePath.Last('/') + 1, 10000);
    TSystemDirectory dir(dirName, dirName);
    TList *files = dir.GetListOfFiles();
    if (files) {
      TSystemFile *sysFile;
      TIter next(files);
      while ((sysFile = (TSystemFile *) next()) && !fileName.EndsWith(".root")) {
        fileName = dirName + sysFile->GetName();
      }
    }
  }

  cout << "fileName = " << fileName << endl;

  TFile *file = new TFile(fileName, "open");
  URun *run = (URun *) file->Get("run");
  if (NULL == run) {
    cout << "No run description in input file." << endl;
  }
  cout << "Target Momentum: " << run->GetPTarg() << endl;
  if (Abs(run->GetPTarg() + run->GetPProj()) < 1e-5) cout << "Input data is in CM frame" << endl;
  double mProton = 0.938272029;
  fA = run->GetAProj() > run->GetATarg() ? run->GetAProj() : run->GetATarg();
  fZ = run->GetZProj() > run->GetZTarg() ? run->GetZProj() : run->GetZTarg();
  fSnn = run->GetNNSqrtS();
  fPcm = 0.5 * sqrt(fSnn * fSnn - 4. * mProton * mProton);
  fElab = 2 * fPcm * fPcm / mProton + mProton;
  fPlab = sqrt(fElab * fElab - mProton * mProton);
  fEkin = fSnn * fSnn / 1.87 - 1.87;
  fBeta = - run->GetPTarg()/(run->GetTargetEnergy()/run->GetATarg());
  cout << "Snn = " << fSnn << " AGeV" << endl;
  cout << "Pcm = " << fPcm << " AGeV" << endl;
  cout << "Plab = " << fPlab << " AGeV" << endl;
  cout << "Elab = " << fElab << " AGeV" << endl;
  cout << "Ekin = " << fEkin << " AGeV" << endl;
  cout << "beta = " << fBeta << endl;
}

void UnigenQA::Run(Int_t nEvents) {
  Long64_t nevents = nEvents < fChain->GetEntries() ? nEvents : fChain->GetEntries();
  Long64_t outputStep = nevents / 10;
  std::cout << "Entries = " << nevents << std::endl;

  for (int i = 0; i < nevents; i++) {
//        if ( (i + 1) % outputStep == 0) std::cout << i + 1 << "/" << nevents << "\r" << std::flush;
    fChain->GetEntry(i);
    FillTracks(); // !!! particle loop goes before the event loop (energy summ is calculated in the former)
    FillEventInfo(); // !!! particle loop goes before the event loop (energy summ is calculated in the former)
  }
}

void UnigenQA::Init_Histograms() {
  gMomentumAxes[kEcm].max = fSnn * fA * 0.5;
  gMomentumAxes[kPcm].max = fSnn * fA * 0.5;
  gMomentumAxes[kMcm].max = fSnn * fA * 0.5;
  gMomentumAxes[kElab].max = fElab * fA * 0.5;
  gMomentumAxes[kPlab].max = fElab * fA * 0.5;
  gMomentumAxes[kMlab].max = fSnn * fA * 0.5;
  gMomentumAxes[kA].max = fA + 2;
  gMomentumAxes[kA].nBins = fA + 2;
  gMomentumAxes[kZ].max = fZ + 2;
  gMomentumAxes[kZ].nBins = fZ + 2;

  if (fReferenceChain == nullptr) fReferenceChain = fChain;
  else {
    cout << "Using reference chain..." << endl;
  }
  Double_t fPSDMax = fElab * 200.;
  Double_t fMmax = fReferenceChain->GetMaximum("fNpa") + 10;

  TString name, title;

  hM = new TH1D("hM", "Track multiplicity;Multiplicity;nEvents", Nint(fMmax), 0, fMmax);
  hB = new TH1D("hB", "Impact parameter; B (fm);nEvents", 170, 0, 17);
  hPsi = new TH1D("hPsi", "#Psi_{RP};#Psi_{RP} (rad);Nevents", 100, -3.15, 3.15);

  hMBcorr = new TH2D("hMBcorr", "M : B;Multiplicity;B (fm)", Nint(fMmax), 0, fMmax, 170, 0, 17);
  h2BAmax = new TH2D("h2BAmax", "B : A_{max};B (fm);A_{max}", 170, 0, 17, 200, 0, 200);
  h2BZmax = new TH2D("h2BZmax", "B : Z_{max};B (fm);Z_{max}", 170, 0, 17, 82, 0, 82);

  /* PSD histogram initialization */
  Int_t nbins = 500;
  for (uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++) {
    for (uint psdCorr = 0; psdCorr < gPSDCorr.size(); psdCorr++) {
      auto group1 = gPSDGroups[gPSDCorr[psdCorr][0]];
      auto group2 = gPSDGroups[gPSDCorr[psdCorr][1]];

      name = "hPSDGroupsCorr_" + group1.name + "_" + group2.name + "_" + fPidGroupNames[pidGroup];
      hPSDGroupsCorr[psdCorr][pidGroup] = new TH2D(name, name, nbins, 0, fPSDMax, nbins, 0, fPSDMax);
      hPSDGroupsCorr[psdCorr][pidGroup]->SetXTitle(group1.displayName.c_str());
      hPSDGroupsCorr[psdCorr][pidGroup]->SetYTitle(group2.displayName.c_str());
    }
    for (auto group : gPSDGroups) {
      name = "hEnergy_" + group.name + "_" + fPidGroupNames[pidGroup];
      title = name + ";E (GeV);nEvents";
      hPSDGroupEnergy[group.id][pidGroup] = new TH1D(name, title, nbins, 0, fPSDMax);
      name = "hPSDMultCorr_" + group.name + "_" + fPidGroupNames[pidGroup];
      title = name + ";Multiplicity;E (GeV)";
      hPSDMultCorr[group.id][pidGroup] = new TH2D(name, title, Nint(fMmax), 0, Nint(fMmax), nbins, 0, fPSDMax);
      name = "hBPSDCorr_" + group.name + "_" + fPidGroupNames[pidGroup];
      title = name + ";E (GeV);B (fm)";
      hBPSDCorr[group.id][pidGroup] = new TH2D(name, title, nbins, 0, fPSDMax, 170, 0, 17);
    }
  }

  for (Int_t iPart = 0; iPart < kParticles; ++iPart) {
    auto particle = gParticles[iPart];

    for (Int_t iMom = 0; iMom < kAxes; ++iMom) {
      auto axis = gMomentumAxes[iMom];
      name = "hTrack" + axis.name + particle.name;
      hTrackMomentum[iMom][iPart] = new TH1D(name, name, axis.nBins, axis.min, axis.max);
      hTrackMomentum[iMom][iPart]->SetXTitle(axis.displayName.c_str());
      hTrackMomentum[iMom][iPart]->SetYTitle("Counts");
    }

    for (uint iCorr = 0; iCorr < gTH2Axes.size(); ++iCorr) {
      vector<int> axes = gTH2Axes[iCorr];

      auto xAxis = gMomentumAxes[axes[0]];
      auto yAxis = gMomentumAxes[axes[1]];

      TString name = "hTrack" + xAxis.name + yAxis.name + particle.name;
      TString title = yAxis.name + " : " + xAxis.name + " (" + particle.name + ")";
      TString xTitle = xAxis.displayName;
      TString yTitle = yAxis.displayName;
      hTrackMomentumCorr.push_back(new TH2D *[kParticles]);
      hTrackMomentumCorr[iCorr][iPart] = new TH2D(name, title,
                                                  xAxis.nBins, xAxis.min, xAxis.max,
                                                  yAxis.nBins, yAxis.min, yAxis.max);
      hTrackMomentumCorr[iCorr][iPart]->SetXTitle(xTitle);
      hTrackMomentumCorr[iCorr][iPart]->SetYTitle(yTitle);
    }

    name = "hYield" + particle.name;
    title = "hYield" + particle.name + ";B (fm);Nparticles";
    hYields[iPart] = new TProfile(name, title, 100, 0, 20);

    for (Int_t iHarm = 0; iHarm < 2; ++iHarm) {
      name = Form("hv%d%s_pT", iHarm + 1, particle.name.c_str());
      title = Form("hv%d%s_pT;p_{T} (GeV/#it{c});v%d", iHarm + 1, particle.name.c_str(), iHarm + 1);
      pVn_pT[iHarm][iPart] = new TProfile(name, name, 50, gMomentumAxes[kPT].min, gMomentumAxes[kPT].max);
      name = Form("hv%d%s_Y", iHarm + 1, particle.name.c_str());
      title = Form("hv%d%s_Y;#it{y};v%d", iHarm + 1, particle.name.c_str(), iHarm + 1);
      pVn_Y[iHarm][iPart] = new TProfile(name, name, 50, gMomentumAxes[kYM].min, gMomentumAxes[kYM].max);
    }
  }

  hPdg = new TH1D("hPdg", "PDG code;PDG code;nCounts", 1000, -500, 3500);

  cout << "Initialization finished" << endl;
}

void UnigenQA::FillEventInfo() {
  double M = event_->GetNpa();
  double B = event_->GetB();
  double psiRP = event_->GetPhi();

  for (uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++) {
    for (auto group : gPSDGroups) {
      hPSDGroupEnergy[group.id][pidGroup]->Fill(fPSDGroupEnergy[group.id][pidGroup]);
      hPSDMultCorr[group.id][pidGroup]->Fill(M, fPSDGroupEnergy[group.id][pidGroup]);
      hBPSDCorr[group.id][pidGroup]->Fill(fPSDGroupEnergy[group.id][pidGroup], B);
    }

    for (uint psdCorr = 0; psdCorr < gPSDCorr.size(); ++psdCorr) {
      auto group1 = gPSDCorr[psdCorr][0];
      auto group2 = gPSDCorr[psdCorr][1];
      hPSDGroupsCorr[psdCorr][pidGroup]->Fill(fPSDGroupEnergy[group1][pidGroup], fPSDGroupEnergy[group2][pidGroup]);
    }
  }

  hM->Fill(M);
  hB->Fill(B);
  hPsi->Fill(psiRP);
  hMBcorr->Fill(M, B);
  h2BAmax->Fill(B, fAmax);
  h2BZmax->Fill(B, fZmax);
}

void UnigenQA::FillTracks() {
  UParticle *track;
  double psiRP = event_->GetPhi();
  Int_t nTracks = event_->GetNpa();
  TLorentzVector momentum;
  Int_t yield[kParticles] = {0};
  Int_t pdg, A, Z;
  double y, theta, Elab, Ecm, Pcm, Plab, Mcm, Mlab;

  fAmax = 0;
  fZmax = 0;
  for (auto group : gPSDGroups) {
    for (uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++) {
      fPSDGroupEnergy[group.id][pidGroup] = 0.;
    }
  }

  for (int i = 0; i < nTracks; i++) {
    track = event_->GetParticle(i);
    pdg = track->GetPdg();
    if (pdg == 22) continue; // patch
    if (pdg / 1000000000 != 0) {
      A = abs(pdg % 10000 / 10);
      Z = abs(pdg % 10000000 / 10000);
    } else {
      A = 1;
      Z = 1;
    }
    if (A > fAmax) fAmax = A;
    if (Z > fZmax) fZmax = Z;
    momentum = track->GetMomentum();
    y = momentum.Rapidity();
    Ecm = momentum.E();
    Pcm = momentum.P();

//				Mcm = sqrt (fabs (Ecm * Ecm - Pcm * Pcm));
    Mcm = sqrt(Ecm * Ecm - Pcm * Pcm);
    if (Ecm - Pcm < 0.) {
      cout << "m^2 = " << Ecm * Ecm - Pcm * Pcm << "\tMcm = " << Mcm << "\tA = " << A << "\tZ = " << Z << "\ti = " << i
           << "\tpdg = " << pdg << endl;
//					momentum.SetE (momentum.P ());
//					Ecm = Pcm;
//					Mcm = 0.;
    }
    momentum.Boost(0., 0., fBeta);
    theta = momentum.Theta();
    Elab = momentum.E();
    Plab = momentum.P();
    Mlab = sqrt(Elab * Elab - Plab * Plab);

//				Elab = track -> E();


    if (abs(pdg) < 3500) hPdg->Fill(pdg);

    double mom[kAxes];
    mom[kPT] = momentum.Pt();
    mom[kETA] = momentum.PseudoRapidity();
    mom[kPHI] = momentum.Phi();
    mom[kYM] = y;
    mom[kEcm] = Ecm;
    mom[kElab] = Elab;
    mom[kPcm] = Pcm;
    mom[kPlab] = Plab;
    mom[kMcm] = Mcm;
    mom[kMlab] = Mlab;
    mom[kA] = A;
    mom[kZ] = Z;
    mom[kMcm_Ecm] = Mcm / Ecm;
    mom[kMlab_Elab] = Mlab / Elab;

    for (Int_t iMom = 0; iMom < kAxes; ++iMom) {
      hTrackMomentum[iMom][0]->Fill(mom[iMom]);
    }

    for (uint j = 0; j < gTH2Axes.size(); ++j) {
      vector<int> axes = gTH2Axes[j];
      hTrackMomentumCorr[j][0]->Fill(mom[axes[0]], mom[axes[1]]);
    }

    for (Int_t iHarm = 0; iHarm < 2; ++iHarm) {
      pVn_pT[iHarm][0]->Fill(mom[kPT], Cos((iHarm + 1) * (mom[kPHI] - psiRP)));
      pVn_Y[iHarm][0]->Fill(mom[kYM], Cos((iHarm + 1) * (mom[kPHI] - psiRP)));
    }
    for (auto group : gPSDGroups) {
      for (uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++) {
        if (theta > group.theta[0] && theta < group.theta[1])
          if (abs(pdg) > fPidGroups[pidGroup][0] && abs(pdg) < fPidGroups[pidGroup][1]) {
            fPSDGroupEnergy[group.id][pidGroup] += Elab;
          }
      }
    }

    for (Int_t iPart = 1; iPart < kParticles; ++iPart) {
      if (pdg / 1000000000 != 0) pdg = 999999;
      if (abs(pdg) < 38) pdg = 99999;
      if (gParticles[iPart].pdg == pdg) {
        yield[iPart] += 1;

        for (Int_t iMom = 0; iMom < kAxes; ++iMom) {
          hTrackMomentum[iMom][iPart]->Fill(mom[iMom]);
        }

        for (uint iCorr = 0; iCorr < gTH2Axes.size(); ++iCorr) {
          vector<int> axes = gTH2Axes[iCorr];
          hTrackMomentumCorr[iCorr][iPart]->Fill(mom[axes[0]], mom[axes[1]]);
        }

        for (Int_t iHarm = 0; iHarm < 2; ++iHarm) {
          pVn_pT[iHarm][iPart]->Fill(mom[kPT], Cos((iHarm + 1) * (mom[kPHI] - psiRP)));
          pVn_Y[iHarm][iPart]->Fill(mom[kYM], Cos((iHarm + 1) * (mom[kPHI] - psiRP)));
        }
      }
    }
  }

  for (Int_t i = 0; i < kParticles; i++) {
    hYields[i]->Fill(event_->GetB(), yield[i]);
  }
}

void UnigenQA::Write_Histograms(const TString filename) {
  auto outputFile = new TFile(filename, "RECREATE");
  if (outputFile->IsOpen()) std::cout << "File '" << filename << "' is opened successfully" << std::endl;
  TDirectory *outputDir;

  hM->Write();
  hB->Write();
  hPsi->Write();
  hMBcorr->Write();
  h2BAmax->Write();
  h2BZmax->Write();
  hPdg->Write();

  for (auto psdGroup : gPSDGroups) {
    outputDir = outputFile->mkdir(psdGroup.name.c_str());
    outputDir->cd();
    for (auto hist : hPSDGroupEnergy[psdGroup.id]) hist->Write();
  }
  for (auto psdGroup : gPSDGroups) {
    outputDir = outputFile->mkdir(Form("%s_Mult", psdGroup.name.c_str()));
    outputDir->cd();
    for (auto hist : hPSDMultCorr[psdGroup.id]) hist->Write();
  }
  for (auto psdGroup : gPSDGroups) {
    outputDir = outputFile->mkdir(Form("%s_B", psdGroup.name.c_str()));
    outputDir->cd();
    for (auto hist : hBPSDCorr[psdGroup.id]) hist->Write();
  }

  for (uint psdCorr = 0; psdCorr < gPSDCorr.size(); psdCorr++) {
    auto group1 = gPSDGroups[gPSDCorr[psdCorr][0]];
    auto group2 = gPSDGroups[gPSDCorr[psdCorr][1]];
    outputDir = outputFile->mkdir((group1.name + group2.name).c_str());
    outputDir->cd();
    for (auto hist : hPSDGroupsCorr[psdCorr]) hist->Write();
  }

  for (auto axis : gMomentumAxes) {
    outputDir = outputFile->mkdir(axis.name.c_str());
    outputDir->cd();
    for (auto hist : hTrackMomentum[axis.id]) hist->Write();
  }

  for (uint iCorr = 0; iCorr < gTH2Axes.size(); ++iCorr) {
    vector<int> axes = gTH2Axes[iCorr];
    auto xAxis = gMomentumAxes[axes[0]];
    auto yAxis = gMomentumAxes[axes[1]];
    outputDir = outputFile->mkdir((xAxis.name + yAxis.name).c_str());
    outputDir->cd();
    for (uint iPart = 0; iPart < kParticles; iPart++) hTrackMomentumCorr[iCorr][iPart]->Write();
  }

  for (Int_t iHarm = 0; iHarm < 2; ++iHarm) {
    outputDir = outputFile->mkdir(Form("V%i_pT", iHarm + 1));
    outputDir->cd();
    for (auto hist : pVn_pT[iHarm]) hist->Write();
    outputDir = outputFile->mkdir(Form("V%i_Y", iHarm + 1));
    outputDir->cd();
    for (auto hist : pVn_Y[iHarm]) hist->Write();
  }

  outputDir = outputFile->mkdir("Yields");
  outputDir->cd();
  for (auto hist : hYields) hist->Write();

  outputFile->Close();
}
