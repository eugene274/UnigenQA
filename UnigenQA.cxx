#include <TLorentzVector.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include "UnigenQA.h"

using namespace TMath;
using namespace std;
using namespace qa;

UnigenQA::UnigenQA()
{
}
UnigenQA::~UnigenQA()
{
}

void UnigenQA::Init(TString filePath, TString treeName)
{
	fChain = new TChain (treeName);
	fChain -> Add (filePath);
	event_ = new UEvent;
	fChain -> SetBranchAddress("event",&event_);
	
	// Get run description
	TString fileName;
	
	if (!filePath.Contains("*")) fileName = filePath;
	else 
	{
		TString dirName = filePath.Remove (filePath.Last('/') + 1, 10000);
		TSystemDirectory dir (dirName, dirName);
		TList *files = dir.GetListOfFiles ();
		if (files) 
		{
      TSystemFile *sysFile;
      TIter next(files);
      while ((sysFile=(TSystemFile*)next()) && !fileName.EndsWith(".root")) 
			{
         fileName = dirName + sysFile -> GetName();
			}
		}
	}
	
	cout << "fileName = " << fileName << endl;

	TFile *file = new TFile (fileName, "open");
  URun* run = (URun*) file -> Get("run");
  if(NULL==run) {
    cout << "No run description in input file." << endl;
  }
  cout << "Target Momentum: " << run->GetPTarg() << endl;
  if(TMath::Abs(run->GetPTarg()) > 0.001) 
	{ 
    cout << "Input data is in CM frame" << endl;
    fElab = (TMath::Power(run->GetNNSqrtS(),2)-
                     2*TMath::Power(0.938271998,2))/(2*0.938271998);
    fPlab = TMath::Sqrt(fElab*fElab - TMath::Power(0.938271998,2));
    fBeta = fPlab / (fElab + 0.938271998);
    fGamma = 1. / TMath::Sqrt(1. - fBeta*fBeta);
		SetSnn (run->GetNNSqrtS());
  } 
	else 
	{
		cout << "Input data is in LAB frame" << endl;
		SetPlab (run->GetPProj());
	}
}

void UnigenQA::SetPlab (double fPlab) 
{ 
	UnigenQA::fPlab = fPlab;
	fElab = sqrt (fPlab * fPlab + 0.938271998 * 0.938271998); 
	fBeta = sqrt (1.- pow (0.938271998 / (fPlab + 0.938271998), 2));
	fGamma = 1. / sqrt (1. - fBeta * fBeta);
	fPc = sqrt( ( 0.938271998 * fElab - 0.938271998 * 0.938271998 )/2 );
	fSnn = sqrt ( 4 * (fPc*fPc + 0.938271998 * 0.938271998) );
	fEkin = fSnn * fSnn / 1.87 - 1.87; 
}


void UnigenQA::Run(Int_t nEvents)
{
    cout << "Snn = " << fSnn << " AGeV" << endl;
    cout << "Elab = " << fElab << " AGeV" << endl;
    cout << "Ekin = " << fEkin << " AGeV" << endl;
    cout << "Plab = " << fPlab << " AGeV" << endl;
    cout << "Pc = " << fPc << " AGeV" << endl;
		cout << "beta = " << fBeta << endl;
		cout << "gamma = " << fGamma << endl << endl;
    Long64_t nevents =  nEvents < fChain->GetEntries() ? nEvents : fChain->GetEntries();
    Long64_t outputStep = nevents / 10;
    std::cout << "Entries = " << nevents << std::endl;

    for (int i=0;i<nevents;i++)
    {
        if ( (i+1) % outputStep == 0) std::cout << i+1 << "/" << nevents << "\r" << std::flush;
        fChain->GetEntry(i);
				FillTracks(); // !!! particle loop goes before the event loop (energy summ is calculated in the former)
        FillEventInfo(); // !!! particle loop goes before the event loop (energy summ is calculated in the former)
    }
}

void UnigenQA::Init_Histograms()
{

    Int_t nbins = 500;

    if (fReferenceChain == nullptr) fReferenceChain = fChain;
    else {
        cout << "Using reference chain..." << endl;
    }
    Double_t fPSDMax = fElab * 200. * 5.;
    Double_t fMmax = fReferenceChain->GetMaximum("fNpa") + 10;

    TString name, title;

    hM = new TH1D("hM","Track multiplicity;Multiplicity;nEvents", Nint(fMmax), 0, fMmax);
    hB = new TH1D("hB","Impact parameter; B (fm);nEvents", 170, 0, 17);
    hPsi = new TH1D("hPsi", "#Psi_{RP};#Psi_{RP} (rad);Nevents", 100, -3.15, 3.15);

    hMBcorr = new TH2D("hMBcorr", "M : B;Multiplicity;B (fm)", Nint(fMmax), 0, fMmax, 170, 0, 17);
    h2BAmax = new TH2D("h2BAmax", "B : A_{max};B (fm);A_{max}", 170, 0, 17, 200, 0, 200);
    h2BZmax = new TH2D("h2BZmax", "B : Z_{max};B (fm);Z_{max}", 170, 0, 17, 82, 0, 82);

    /* PSD histogram initialization */
		for (uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++)
		{
			for (uint psdCorr = 0; psdCorr < gPSDCorr.size (); psdCorr++) {
					auto group1 = gPSDGroups[gPSDCorr [psdCorr][0]];
					auto group2 = gPSDGroups[gPSDCorr [psdCorr][1]];

					name = "hPSDGroupsCorr_" + group1.name + "_" + group2.name + "_" + fPidGroupNames [pidGroup];
					hPSDGroupsCorr [psdCorr][pidGroup] = new TH2D(name, name, nbins, 0, fPSDMax, nbins, 0, fPSDMax);
					hPSDGroupsCorr [psdCorr][pidGroup] -> SetXTitle (group1.displayName.c_str());
					hPSDGroupsCorr [psdCorr][pidGroup] -> SetYTitle (group2.displayName.c_str());
			}
			for (auto group : gPSDGroups) {
					name = "hEnergy_" + group.name + "_" + fPidGroupNames [pidGroup];
					title = name + ";E (GeV);nEvents";
					hPSDGroupEnergy [group.id][pidGroup] = new TH1D(name, title, nbins, 0, fPSDMax);
					name = "hPSDMultCorr_" + group.name + "_" + fPidGroupNames [pidGroup];
					title = name + ";Multiplicity;E (GeV)";
					hPSDMultCorr [group.id][pidGroup] = new TH2D(name, title, Nint(fMmax), 0, Nint(fMmax), nbins, 0, fPSDMax);
					name = "hBPSDCorr_" + group.name + "_" + fPidGroupNames [pidGroup];
					title = name + ";E (GeV);B (fm)";
					hBPSDCorr [group.id][pidGroup] = new TH2D(name, title, nbins, 0, fPSDMax, 170, 0, 17);
			}
		}
		
		for (Int_t iPart=0; iPart<kParticles; ++iPart) 
		{
			auto particle = gParticles[iPart];
			
			for (Int_t iMom=0; iMom<kAxes; ++iMom) 
			{
					auto axis = gMomentumAxes[iMom];
					name = "hTrack" + axis.name + particle.name;
					hTrackMomentum[iMom][iPart] = new TH1D(name, name, axis.nBins, axis.min, axis.max);
					hTrackMomentum[iMom][iPart]->SetXTitle(axis.displayName.c_str());
					hTrackMomentum[iMom][iPart]->SetYTitle("Counts");
			}
			
			for (uint iCorr = 0; iCorr < gTH2Axes.size(); ++iCorr) 
			{
					vector <int> axes = gTH2Axes[iCorr];

					auto xAxis = gMomentumAxes[axes[0]];
					auto yAxis = gMomentumAxes[axes[1]];

					TString name = "hTrack" + xAxis.name + yAxis.name + particle.name;
					TString xTitle = xAxis.displayName;
					TString yTitle = yAxis.displayName;
					hTrackMomentumCorr[iCorr][iPart] = new TH2D(name, "",
																					 xAxis.nBins, xAxis.min, xAxis.max,
																					 yAxis.nBins, yAxis.min, yAxis.max);
					hTrackMomentumCorr[iCorr][iPart]->SetXTitle(xTitle);
					hTrackMomentumCorr[iCorr][iPart]->SetYTitle(yTitle);
			}
					
			name = "hYield" + particle.name;
			title = "hYield" + particle.name + ";B (fm);Nparticles";
			hYields[iPart] = new TProfile(name, title, 100, 0, 20);

			for (Int_t iHarm=0; iHarm<2; ++iHarm) {
//				auto& particle = gParticles[iPart];
				name = Form("hv%d%s_pT", iHarm + 1, particle.name.c_str());
				title = Form("hv%d%s_pT;p_{T} (GeV/#it{c});v%d", iHarm + 1, particle.name.c_str(), iHarm + 1);
				pVn_pT[iHarm][iPart] = new TProfile(name, name, 50, gMomentumAxes[kPT].min, gMomentumAxes[kPT].max);
				name = Form("hv%d%s_Y", iHarm + 1, particle.name.c_str());
				title = Form("hv%d%s_Y;#it{y};v%d", iHarm + 1, particle.name.c_str(), iHarm + 1);
				pVn_Y[iHarm][iPart] = new TProfile(name, name, 50, gMomentumAxes[kYM].min, gMomentumAxes[kYM].max);
			}
		}
		
		hPdg = new TH1D("hPdg","PDG code;PDG code;nCounts", 1000, -500, 3500);
		hA = new TH1D("hA","Mass number; A", 200, 0, 200);
		hZ = new TH1D("hZ","Charge number; Z", 82, 0, 82);
		h2ZA = new TH2D("h2ZA","A : Z; A; Z", 200, 0, 200, 82, 0, 82);
		h2EcmA = new TH2D ("h2EcmA", "Ecm:A;A;Ecm (GeV)", 200, 0, 200, gMomentumAxes[kEcm].nBins, gMomentumAxes[kEcm].min, gMomentumAxes[kEcm].max);
		h2ElabA = new TH2D ("h2ElabA", "Elab :A;A;Elab (GeV)", 200, 0, 200, gMomentumAxes[kElab].nBins, gMomentumAxes[kElab].min, gMomentumAxes[kElab].max);
	

    cout << "Initialization finished" << endl;
}

void UnigenQA::FillEventInfo()
{
    double M = event_ -> GetNpa();
    double B = event_ -> GetB();
    double psiRP = event_ -> GetPhi();
		
		for (uint pidGroup = 0; pidGroup < fPidGroups.size(); pidGroup++)
		{
			for (auto group : gPSDGroups) {
					hPSDGroupEnergy [group.id][pidGroup] -> Fill (fPSDGroupEnergy [group.id][pidGroup]);
					hPSDMultCorr [group.id][pidGroup] -> Fill(M, fPSDGroupEnergy [group.id][pidGroup]);
					hBPSDCorr [group.id][pidGroup] -> Fill(fPSDGroupEnergy [group.id][pidGroup], B);
			}
			
			for (uint psdCorr = 0; psdCorr < gPSDCorr.size (); ++psdCorr) {
					auto group1 = gPSDCorr [psdCorr][0];
					auto group2 = gPSDCorr [psdCorr][1];
					hPSDGroupsCorr [psdCorr][pidGroup] -> Fill(fPSDGroupEnergy [group1][pidGroup], fPSDGroupEnergy [group2][pidGroup]);
			}
		}
		
    hM -> Fill(M);
    hB -> Fill(B);
    hPsi -> Fill(psiRP);
    hMBcorr -> Fill(M, B);
		h2BAmax -> Fill (B, fAmax);
		h2BZmax -> Fill (B, fZmax);
}

void UnigenQA::FillTracks()
{
    UParticle* track;
    double psiRP = event_ -> GetPhi();
    Int_t nTracks = event_->GetNpa();
		TLorentzVector momentum;
    Int_t yield[kParticles] = {0};
		Int_t pdg, mass, charge;
		double y, theta, Elab, Ecm, Pcm, Plab;
		
		fAmax = 0;
		fZmax = 0;
		for (auto group : gPSDGroups) 
		{
			for (uint pidGroup = 0; pidGroup < fPidGroups.size (); pidGroup++)
			{
				fPSDGroupEnergy [group.id][pidGroup] = 0.;
			}
		}

    for (int i=0;i<nTracks; i++)
    {
        track = event_ -> GetParticle(i);
				pdg = track -> GetPdg();
				if (pdg / 1000000000 != 0) 
				{
					mass = abs (pdg % 10000 / 10);
					charge = abs (pdg % 10000000 / 10000); 
				}
				else 
				{
					mass = 1;
					charge = 1;
				}
				if (mass > fAmax) fAmax = mass;
				if (charge > fZmax) fZmax = charge;
				momentum = track -> GetMomentum();
				y = momentum.Rapidity();
				Ecm = momentum.E ();
				Pcm = momentum.P ();
				momentum.Boost (0., 0., fBeta);
				theta = momentum.Theta ();
				Elab = momentum.E ();
				Plab = momentum.P ();
//				Elab = track -> E();

				if (abs (pdg) < 3500) hPdg -> Fill (pdg);
				hA -> Fill (mass);
				hZ -> Fill (charge);
				h2ZA -> Fill (mass, charge);
				h2EcmA -> Fill (mass, Ecm);
				h2ElabA -> Fill (mass, Elab);
				
        double mom[8] = {momentum.Pt(), momentum.PseudoRapidity(), momentum.Phi(), y, Ecm, Elab, Pcm, Plab};

        for (Int_t iMom=0; iMom<kAxes; ++iMom) {
            hTrackMomentum[iMom][0] -> Fill( mom[iMom] );
        }
				
        for (uint j = 0; j < gTH2Axes.size (); ++j) {
            vector <int> axes = gTH2Axes[j];
            hTrackMomentumCorr[j][0] -> Fill(mom[axes[0]], mom[axes[1]]);
        }

        for (Int_t iHarm=0; iHarm<2; ++iHarm) {
						pVn_pT[iHarm][0] -> Fill ( mom[kPT], Cos( (iHarm+1)*(mom[kPHI] - psiRP) ) );
						pVn_Y[iHarm][0] ->  Fill ( mom[kYM] , Cos( (iHarm+1)*(mom[kPHI] - psiRP) ) );
				}
				
				for (auto group : gPSDGroups) 
				{
					for (uint pidGroup = 0; pidGroup < fPidGroups.size (); pidGroup++)
					{
						if (theta > group.theta [0] && theta < group.theta [1])
							if (abs (pdg) > fPidGroups [pidGroup][0] && abs (pdg) < fPidGroups [pidGroup][1])
								fPSDGroupEnergy [group.id][pidGroup] += Elab;
					}
				}
				
        for(Int_t iPart=0; iPart<kParticles; ++iPart) {
            if (gParticles[iPart].pdg == pdg) 
						{
                yield[iPart] += 1;
								
								for (Int_t iMom=0; iMom<kAxes; ++iMom) {
										hTrackMomentum[iMom][iPart] -> Fill( mom[iMom] );
								}
								
								for (uint iCorr = 0; iCorr < gTH2Axes.size (); ++iCorr) {
										vector <int> axes = gTH2Axes[iCorr];
										hTrackMomentumCorr[iCorr][iPart] -> Fill(mom[axes[0]], mom[axes[1]]);
								}

                for (Int_t iHarm=0; iHarm<2; ++iHarm) {
                    pVn_pT[iHarm][iPart] -> Fill ( mom[kPT], Cos( (iHarm+1)*(mom[kPHI] - psiRP) ) );
                    pVn_Y[iHarm][iPart] ->  Fill ( mom[kYM] , Cos( (iHarm+1)*(mom[kPHI] - psiRP) ) );
                }
            }
        }
    }

    for (Int_t i = 0; i < kParticles; i++) {
        hYields[i]->Fill(event_->GetB (), yield[i]);
    }
}


void UnigenQA::Write_Histograms(const TString filename)
{
    auto outputFile = new TFile(filename, "RECREATE");
    if (outputFile->IsOpen() ) std::cout << "File '" << filename << "' is opened successfully" << std::endl;
    TDirectory *outputDir;

    hM -> Write();
    hB -> Write();
    hPsi -> Write();
    hMBcorr -> Write();
    h2BAmax -> Write();
    h2BZmax -> Write();
		hPdg -> Write ();
		hA -> Write ();
		hZ -> Write ();
		h2ZA -> Write ();
		h2EcmA -> Write ();
		h2ElabA -> Write ();
		
		for (auto psdGroup : gPSDGroups)
		{
			outputDir = outputFile -> mkdir (psdGroup.name.c_str ());
			outputDir -> cd ();
			for (auto hist : hPSDGroupEnergy [psdGroup.id]) hist -> Write();
		}
		for (auto psdGroup : gPSDGroups)
		{
			outputDir = outputFile -> mkdir (Form ("%s_Mult", psdGroup.name.c_str () ));
			outputDir -> cd ();
			for (auto hist : hPSDMultCorr [psdGroup.id]) hist -> Write();
		}
		for (auto psdGroup : gPSDGroups)
		{
			outputDir = outputFile -> mkdir (Form ("%s_B", psdGroup.name.c_str () ));
			outputDir -> cd ();
			for (auto hist : hBPSDCorr [psdGroup.id]) hist -> Write();
		}		
		
		for (uint psdCorr = 0; psdCorr < gPSDCorr.size (); psdCorr++)
		{
			auto group1 = gPSDGroups[gPSDCorr [psdCorr][0]];
			auto group2 = gPSDGroups[gPSDCorr [psdCorr][1]];
			outputDir = outputFile -> mkdir ((group1.name + group2.name).c_str ());
			outputDir -> cd ();
			for (auto hist : hPSDGroupsCorr [psdCorr]) hist -> Write ();
		}
	
//		for (uint pidGroup = 0; pidGroup < fPidGroups.size (); pidGroup++)
//		{
//			for (auto hist : hPSDGroupEnergy [pidGroup]) hist -> Write();
//			for (auto hist : hPSDGroupsCorr [pidGroup]) hist -> Write();
//			for (auto hist : hPSDMultCorr [pidGroup]) hist -> Write();
//			for (auto hist : hBPSDCorr [pidGroup]) hist -> Write();
//		}		
			
		for (auto axis : gMomentumAxes) {
			outputDir = outputFile -> mkdir (axis.name.c_str ());
			outputDir -> cd ();
			for (auto hist : hTrackMomentum [axis.id]) hist -> Write();
		}
		
		for (uint iCorr = 0; iCorr < gTH2Axes.size(); ++iCorr) 
		{
				vector <int> axes = gTH2Axes[iCorr];
				auto xAxis = gMomentumAxes[axes[0]];
				auto yAxis = gMomentumAxes[axes[1]];
				outputDir = outputFile -> mkdir ((xAxis.name + yAxis.name).c_str ());
				outputDir -> cd ();
				for (auto hist : hTrackMomentumCorr [iCorr]) hist -> Write();
		}

		for (Int_t iHarm=0; iHarm<2; ++iHarm)
		{
			outputDir = outputFile -> mkdir (Form ("V%i_pT", iHarm + 1));
			outputDir -> cd ();
			for (auto hist : pVn_pT[iHarm]) hist -> Write();
			outputDir = outputFile -> mkdir (Form ("V%i_Y", iHarm + 1));
			outputDir -> cd ();
			for (auto hist : pVn_Y[iHarm]) hist -> Write();
		}
		
		outputDir = outputFile -> mkdir ("Yields");
		outputDir -> cd ();
		for (auto hist : hYields) hist -> Write ();
	
		outputFile->Close();
}