#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <fstream>
#include <iostream>
#include <string>

#include "URun.h"
#include "UEvent.h"

using namespace std;

void dat2root(TString inFile = "./dcmqgsm_12.dat",
              TString outFile = "./dcmqgsm_12.root",
              Bool_t isHyperNucleus = false,
              Double_t pBeam = 13.) {
  enum EParticles {
    kPion = 0,
    kKaon,
    kNucleon,
    kLambda,
    kDeutron,
    kSigma,
    kXi,
    kDeutron2,
    kDeutronLambda,
    kTriton,
    kTritonLambda,
    kLambdaLambda,
    kHe4Lambda,
    kParticles
  };
  Double_t masses[kParticles];
  masses[kPion] = 0.140;
  masses[kKaon] = 0.48;
  masses[kNucleon] = 0.938;
  masses[kLambda] = 1.115;
  masses[kDeutron] = 1.876;
  masses[kSigma] = 1.192;
  masses[kXi] = 1.31;
  masses[kDeutron2] = 2.054;
  masses[kDeutronLambda] = 2.809;
  masses[kTriton] = 2.231;
  masses[kTritonLambda] = 2.992;
  masses[kLambdaLambda] = 3.728;
  masses[kHe4Lambda] = 3.925;
  //   pi    K     p,n   Lambda   d    sigma   Xi d(Lambda)t(Lambda) (2Lambda) 4He(Lambda)

  double mProton = 0.938272029;
  double eBeam = TMath::Sqrt(pBeam * pBeam + mProton * mProton);
  double pCM = TMath::Sqrt(0.5 * mProton * (eBeam - mProton));

  const char *generator = "DCM_QGSM++";
  const char *comment = "";
  Int_t aProj = 197;
  Int_t zProj = 79;
  Double_t pProj = pCM;
  Int_t aTarg = 197;
  Int_t zTarg = 79;
  Double_t pTarg = -pCM;
  Double_t bMin = 0;
  Double_t bMax = 20.0;
  Int_t bWeight = 0;
  Double_t phiMin = 0;
  Double_t phiMax = 2 * TMath::Pi();
  Double_t sigma = 0;
  Int_t nEvents = 1000;

  URun *header = new URun("DCM_QGSM",
                          "",
                          aProj,
                          zProj,
                          pProj,
                          aTarg,
                          zTarg,
                          pTarg,
                          bMin,
                          bMax,
                          bWeight,
                          phiMin,
                          phiMax,
                          sigma,
                          nEvents);
  UEvent *event = new UEvent;

  TFile *f = new TFile(outFile, "recreate");
  TTree *tree = new TTree("events", "Input data, DCM_QGSM");

  header->Write();
  tree->Branch("event", "UEvent", event);

  ifstream *InputFile = new ifstream(inFile);
  FILE *fp = fopen(inFile.Data(), "r");

  Int_t eventId, nTracks, iMass, iCharge, iStrangeness;
  Float_t bx, by, px, py, pz, e;
  Float_t b;
  Float_t mass = 0.;
  Float_t ExEnergy;
  Int_t pdgType = 0;

  Int_t ii = 0;

  while (!InputFile->eof()) {
    ii++;

    if (ii != 3) continue; // patch

    std::cout << "Event # " << ii << "... \r" << std::flush;
    *InputFile >> eventId >> b >> bx >> by;
    // std::cout << eventId << "\t" << b << "\t" << bx << "\t" << by << std::endl;

    Int_t eventNr = eventId;

    Double_t phi = TMath::ATan2(by, bx);

    Int_t nes = 0;
    Int_t stepNr = 0;
    Double_t stepT = 0;

    if (InputFile->eof()) break;

    event->SetParameters(eventNr, b, phi, nes, stepNr, stepT);

    Double_t eTotal = .0;
    for (Int_t i = 0; i < 2; i++) {
      *InputFile >> nTracks;

      std::cout << "n_Spec = " << nTracks << std::endl;

      for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
        pdgType = -1;
        *InputFile >> iMass >> iCharge >> iStrangeness >> ExEnergy >> px >> py >> pz;

        if (iMass > 1 && iCharge > 0) pdgType = iMass * 10 + iCharge * 1e4 + 1e9;
        else if (iMass == 1 && iCharge == 1) pdgType = 2212;
        else if (iMass == 1 && iCharge == -1) pdgType = -2212;
        else if (iMass == 1 && iCharge == 0) pdgType = 2112;
        else if (iMass == 0 && iCharge == 0) pdgType = 22;
        else
          std::cout << "Undef spectator " << iMass << "  " << iCharge << "  " << iStrangeness << "   " << px << "  "
                    << py << "  " << pz << "  " << std::endl;
        if (pdgType == -1) continue;

        Int_t index = iTrack;
        Int_t pdg = pdgType;
        Int_t status = 0;
        Int_t parent = 0;
        Int_t parentDecay = 0;
        Int_t mate = 0;
        Int_t decay = 0;
        Int_t child[2] = {0, 0};
        Double_t x = 0;
        Double_t y = 0;
        Double_t z = 0;
        Double_t t = 2 * i - 1; // -1 for target spectators, +1 for projectile
        Double_t weight = 1;

        mass = iMass * 0.938271998;
        e = TMath::Sqrt(px * px + py * py + pz * pz + mass * mass);
        // cout << iMass << "  " << iCharge << "  "  << iStrangeness << "   " << ExEnergy << "   " << px << "  " << py << "  " << pz << "  " << e << "  " << endl;

        event->AddParticle(index,
                           pdg,
                           status,
                           parent,
                           parentDecay,
                           mate,
                           decay,
                           child,
                           px,
                           py,
                           pz,
                           e,
                           x,
                           y,
                           z,
                           t,
                           weight);
        eTotal += e;
      }
    }

    *InputFile >> nTracks;
    std::cout << "n_Prod = " << nTracks << std::endl;

    Int_t iTrack = 0;

    while (iTrack < nTracks) {
      *InputFile >> iMass >> iCharge >> iStrangeness >> px >> py >> pz >> mass;
      // cout << iMass << "  " << iCharge << "  "  << iStrangeness << "   " << px << "  " << py << "  " << pz << "  " << mass << "  " << endl;
      pdgType = -1;

      if (TMath::Abs(mass - masses[0]) < 0.05 && iCharge == 1) pdgType = 211;
      else if (TMath::Abs(mass - masses[0]) < 0.05 && iCharge == 0) pdgType = 111;
      else if (TMath::Abs(mass - masses[0]) < 0.05 && iCharge == -1) pdgType = -211;

      else if (TMath::Abs(mass - masses[1]) < 0.05 && iCharge == 1) pdgType = 321;
      else if (TMath::Abs(mass - masses[1]) < 0.05 && iCharge == 0) pdgType = 310;
      else if (TMath::Abs(mass - masses[1]) < 0.05 && iCharge == -1) pdgType = -321;

      else if (TMath::Abs(mass - masses[2]) < 0.05 && iCharge == 1) pdgType = 2212;
      else if (TMath::Abs(mass - masses[2]) < 0.05 && iCharge == 0) pdgType = 2112;
      else if (TMath::Abs(mass - masses[2]) < 0.05 && iCharge == -1) pdgType = -2212;

      else if (TMath::Abs(mass - masses[3]) < 0.02 && iCharge == 0) pdgType = 3122;

      else if (TMath::Abs(mass - masses[4]) < 0.1 && iCharge == 1) pdgType = 1000010020;

      else if (TMath::Abs(mass - masses[5]) < 0.03 && iCharge == 1) pdgType = 3222;
      else if (TMath::Abs(mass - masses[5]) < 0.03 && iCharge == 0) pdgType = 3212;
      else if (TMath::Abs(mass - masses[5]) < 0.03 && iCharge == -1) pdgType = 3112;

      else if (TMath::Abs(mass - masses[6]) < 0.05 && iCharge == -1) pdgType = 3322;
      else if (TMath::Abs(mass - masses[6]) < 0.05 && iCharge == 0) pdgType = 3312;

      else if (TMath::Abs(mass - masses[7]) < 0.05) pdgType = 1000010020 + 10000000 * isHyperNucleus;
      else if (TMath::Abs(mass - masses[8]) < 0.05) pdgType = 1000010030 + 10000000 * isHyperNucleus;
      else if (TMath::Abs(mass - masses[9]) < 0.05) pdgType = 1000000020 + 20000000 * isHyperNucleus;
      else if (TMath::Abs(mass - masses[10]) < 0.05) pdgType = 1000010030 + 10000000 * isHyperNucleus;
      else if (TMath::Abs(mass - masses[11]) < 0.05) pdgType = 1000020040 + 10000000 * isHyperNucleus;
      else if (TMath::Abs(mass - masses[12]) < 0.05) pdgType = 1000010040 + 10000000 * isHyperNucleus;

      else if (TMath::Abs(mass) < 0.01 && iCharge == -1) pdgType = 11;
      else if (TMath::Abs(mass) < 0.01 && iCharge == 1) pdgType = -11;
      else if (TMath::Abs(mass) < 0.01 && iCharge == 0) pdgType = 22;

      if (pdgType == -1) {
        std::cout << "Undef produced " << iMass << "  " << iCharge << "  " << iStrangeness << "   " << mass
                  << std::endl;
        iTrack++;
        continue;
      }

      Int_t index = iTrack;
      Int_t pdg = pdgType;
      Int_t status = 0;//iStrangeness;
      Int_t parent = 0;//iCharge;
      Int_t parentDecay = 0;//iMass;
      Int_t mate = 0;
      Int_t decay = 0;
      Int_t child[2] = {0, 0};
      Double_t x = 0;
      Double_t y = 0;
      Double_t z = 0;
      Double_t t = 0;//mass;
      Double_t weight = 1;

      e = TMath::Sqrt(px * px + py * py + pz * pz + mass * mass);
      event->AddParticle(index,
                         pdg,
                         status,
                         parent,
                         parentDecay,
                         mate,
                         decay,
                         child,
                         px,
                         py,
                         pz,
                         e,
                         x,
                         y,
                         z,
                         t,
                         weight);
      eTotal += e;

      iTrack++;
    }

    cout << "E_Total = " << eTotal << endl;

    tree->Fill();
    event->Clear();
  }

  InputFile->close();
  tree->Write();
  f->Close();
}
