// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief QA task for V0 analysis using derived data
///
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)

#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/v0qaanalysis.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-fill", VariantType::Int, 1, {"Add histogram filling"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using DauTracksMC = soa::Join<DauTracks, aod::McTrackLabels>;
using V0Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::CentFV0As>;

struct LfV0qaanalysis {

  // Produces
  Produces<aod::MyV0Candidates> myv0s;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    int nProc = 0;
    if (doprocessData) {
      nProc += 1;
    }
    if (doprocessMCReco) {
      nProc += 1;
    }
    if (doprocessMCGen) {
      nProc += 1;
    }
    if (nProc == 0) {
      LOG(fatal) << "Enable at least one process function";
    }
    LOG(info) << "Number of process functions enabled: " << nProc;

    registry.add("hNEvents", "hNEvents", {HistType::kTH1I, {{10, 0.f, 10.f}}});
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "all");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel8");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "TVX");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "zvertex");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "TFBorder");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(6, "ITSROFBorder");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(7, "isTOFVertexMatched");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(8, "isGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(9, "Applied selection");

    registry.add("hCentFT0M", "hCentFT0M", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    registry.add("hCentFV0A", "hCentFV0A", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    if (isMC) {
      registry.add("hCentFT0M_RecoColl_MC", "hCentFT0M_RecoColl_MC", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
      registry.add("hCentFT0M_RecoColl_MC_INELgt0", "hCentFT0M_RecoColl_MC_INELgt0", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
      registry.add("hCentFT0M_GenRecoColl_MC", "hCentFT0M_GenRecoColl_MC", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
      registry.add("hCentFT0M_GenRecoColl_MC_INELgt0", "hCentFT0M_GenRecoColl_MC_INELgt0", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
      registry.add("hCentFT0M_GenColl_MC", "hCentFT0M_GenColl_MC", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
      registry.add("hCentFT0M_GenColl_MC_INELgt0", "hCentFT0M_GenColl_MC_INELgt0", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
      registry.add("hNEventsMCGen", "hNEventsMCGen", {HistType::kTH1I, {{4, 0.f, 4.f}}});
      registry.get<TH1>(HIST("hNEventsMCGen"))->GetXaxis()->SetBinLabel(1, "all");
      registry.get<TH1>(HIST("hNEventsMCGen"))->GetXaxis()->SetBinLabel(2, "zvertex_true");
      registry.get<TH1>(HIST("hNEventsMCGen"))->GetXaxis()->SetBinLabel(3, "INELgt0_true");
      registry.add("hNEventsMCGenReco", "hNEventsMCGenReco", {HistType::kTH1I, {{2, 0.f, 2.f}}});
      registry.get<TH1>(HIST("hNEventsMCGenReco"))->GetXaxis()->SetBinLabel(1, "INEL");
      registry.get<TH1>(HIST("hNEventsMCGenReco"))->GetXaxis()->SetBinLabel(2, "INELgt0");
      registry.add("hNEventsMCReco", "hNEventsMCReco", {HistType::kTH1I, {{4, 0.f, 4.f}}});
      registry.get<TH1>(HIST("hNEventsMCReco"))->GetXaxis()->SetBinLabel(1, "all");
      registry.get<TH1>(HIST("hNEventsMCReco"))->GetXaxis()->SetBinLabel(2, "pass ev sel");
      registry.get<TH1>(HIST("hNEventsMCReco"))->GetXaxis()->SetBinLabel(3, "INELgt0");
      registry.get<TH1>(HIST("hNEventsMCReco"))->GetXaxis()->SetBinLabel(4, "check");
      registry.add("Reconstructed_MCRecoColl_INEL_K0Short", "Reconstructed_MCRecoColl_INEL_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Reconstructed_MCRecoColl_INEL_Lambda", "Reconstructed_MCRecoColl_INEL_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Reconstructed_MCRecoColl_INEL_AntiLambda", "Reconstructed_MCRecoColl_INEL_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Reconstructed_MCRecoColl_INELgt0_K0Short", "Reconstructed_MCRecoColl_INELgt0_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Reconstructed_MCRecoColl_INELgt0_Lambda", "Reconstructed_MCRecoColl_INELgt0_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Reconstructed_MCRecoColl_INELgt0_AntiLambda", "Reconstructed_MCRecoColl_INELgt0_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenRecoColl_INEL_K0Short", "Generated_MCGenRecoColl_INEL_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenRecoColl_INEL_Lambda", "Generated_MCGenRecoColl_INEL_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenRecoColl_INEL_AntiLambda", "Generated_MCGenRecoColl_INEL_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoColl_INEL_K0Short", "Generated_MCRecoColl_INEL_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoColl_INEL_Lambda", "Generated_MCRecoColl_INEL_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoColl_INEL_AntiLambda", "Generated_MCRecoColl_INEL_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoCollCheck_INEL_K0Short", "Generated_MCRecoCollCheck_INEL_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoCollCheck_INEL_Lambda", "Generated_MCRecoCollCheck_INEL_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoCollCheck_INEL_AntiLambda", "Generated_MCRecoCollCheck_INEL_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenColl_INEL_K0Short", "Generated_MCGenColl_INEL_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenColl_INEL_Lambda", "Generated_MCGenColl_INEL_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenColl_INEL_AntiLambda", "Generated_MCGenColl_INEL_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenRecoColl_INELgt0_K0Short", "Generated_MCGenRecoColl_INELgt0_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenRecoColl_INELgt0_Lambda", "Generated_MCGenRecoColl_INELgt0_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenRecoColl_INELgt0_AntiLambda", "Generated_MCGenRecoColl_INELgt0_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoColl_INELgt0_K0Short", "Generated_MCRecoColl_INELgt0_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoColl_INELgt0_Lambda", "Generated_MCRecoColl_INELgt0_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoColl_INELgt0_AntiLambda", "Generated_MCRecoColl_INELgt0_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoCollCheck_INELgt0_K0Short", "Generated_MCRecoCollCheck_INELgt0_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoCollCheck_INELgt0_Lambda", "Generated_MCRecoCollCheck_INELgt0_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCRecoCollCheck_INELgt0_AntiLambda", "Generated_MCRecoCollCheck_INELgt0_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenColl_INELgt0_K0Short", "Generated_MCGenColl_INELgt0_K0Short", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenColl_INELgt0_Lambda", "Generated_MCGenColl_INELgt0_Lambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
      registry.add("Generated_MCGenColl_INELgt0_AntiLambda", "Generated_MCGenColl_INELgt0_AntiLambda", {HistType::kTH2F, {{250, 0.f, 25.f}, {1000, 0.f, 100.f}}});
    }
    registry.print();
  }

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 15.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 0, "Apply sel8 event selection"};
  Configurable<bool> isMC{"isMC", 0, "Is MC"};
  Configurable<bool> isTriggerTVX{"isTriggerTVX", 1, "Is Trigger TVX"};
  Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", 1, "Is No Time Frame Border"};
  Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", 1, "Is No ITS Readout Frame Border"};
  Configurable<bool> isVertexTOFmatched{"isVertexTOFmatched", 0, "Is Vertex TOF matched"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", 0, "isGoodZvtxFT0vsPV"};

  // V0 selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 10, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.0, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.0, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 0.0, "Radius"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};

  // Event selection
  template <typename TCollision>
  bool AcceptEvent(TCollision const& collision)
  {
    registry.fill(HIST("hNEvents"), 0.5);
    if (sel8 && !collision.sel8()) {
      return false;
    }
    registry.fill(HIST("hNEvents"), 1.5);
    if (isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    registry.fill(HIST("hNEvents"), 2.5);
    if (TMath::Abs(collision.posZ()) > cutzvertex) {
      return false;
    }
    registry.fill(HIST("hNEvents"), 3.5);
    if (isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    registry.fill(HIST("hNEvents"), 4.5);
    if (!isMC && isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    registry.fill(HIST("hNEvents"), 5.5);
    if (isVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    registry.fill(HIST("hNEvents"), 6.5);
    if (isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    registry.fill(HIST("hNEvents"), 7.5);

    return true;
  }

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&&
                                                         nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  void processData(V0Collisions::iterator const& collision,
                   soa::Filtered<aod::V0Datas> const& V0s,
                   DauTracks const& /*tracks*/)
  {

    // Apply event selection
    if (!AcceptEvent(collision)) {
      return;
    }
    registry.fill(HIST("hNEvents"), 8.5);
    registry.fill(HIST("hCentFT0M"), collision.centFT0M());
    registry.fill(HIST("hCentFV0A"), collision.centFV0A());

    for (auto& v0 : V0s) { // loop over V0s

      // c tau
      float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
      float ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

      // ITS clusters
      int posITSNhits = 0, negITSNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (v0.posTrack_as<DauTracks>().itsClusterMap() & (1 << i)) {
          posITSNhits++;
        }
        if (v0.negTrack_as<DauTracks>().itsClusterMap() & (1 << i)) {
          negITSNhits++;
        }
      }

      // Event flags
      int evFlag = 0;
      if (collision.isInelGt0()) {
        evFlag = 1;
      }

      int lPDG = 0;
      bool isPhysicalPrimary = isMC;

      if (v0.v0radius() > v0radius &&
          v0.v0cosPA() > v0cospa &&
          TMath::Abs(v0.posTrack_as<DauTracks>().eta()) < etadau &&
          TMath::Abs(v0.negTrack_as<DauTracks>().eta()) < etadau) {

        // Fill table
        myv0s(v0.globalIndex(), v0.pt(), v0.yLambda(), v0.yK0Short(),
              v0.mLambda(), v0.mAntiLambda(), v0.mK0Short(),
              v0.v0radius(), v0.v0cosPA(),
              v0.dcapostopv(), v0.dcanegtopv(), v0.dcaV0daughters(),
              v0.posTrack_as<DauTracks>().eta(), v0.negTrack_as<DauTracks>().eta(),
              v0.posTrack_as<DauTracks>().phi(), v0.negTrack_as<DauTracks>().phi(),
              posITSNhits, negITSNhits, ctauLambda, ctauAntiLambda, ctauK0s,
              v0.negTrack_as<DauTracks>().tpcNSigmaPr(), v0.posTrack_as<DauTracks>().tpcNSigmaPr(),
              v0.negTrack_as<DauTracks>().tpcNSigmaPi(), v0.posTrack_as<DauTracks>().tpcNSigmaPi(),
              v0.negTrack_as<DauTracks>().tofNSigmaPr(), v0.posTrack_as<DauTracks>().tofNSigmaPr(),
              v0.negTrack_as<DauTracks>().tofNSigmaPi(), v0.posTrack_as<DauTracks>().tofNSigmaPi(),
              v0.posTrack_as<DauTracks>().hasTOF(), v0.negTrack_as<DauTracks>().hasTOF(), lPDG, isPhysicalPrimary,
              collision.centFT0M(), collision.centFV0A(), evFlag, v0.alpha(), v0.qtarm());
      }
    }
  }
  PROCESS_SWITCH(LfV0qaanalysis, processData, "Process data", true);

  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  SliceCache cache1;

  Service<o2::framework::O2DatabasePDG> pdgDB;

  void processMCReco(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::PVMults> const& collisions,
                     soa::Join<aod::McCollisions, aod::McCentFT0Ms> const& /*mcCollisions*/,
                     soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                     aod::McParticles const& mcParticles, DauTracksMC const& /*tracks*/)
  {
    for (const auto& collision : collisions) {
      // Apply event selection

      if (!AcceptEvent(collision)) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& mcCollision = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();

      registry.fill(HIST("hNEventsMCReco"), 3.5);
      const float cent = 0.f;

      // Event flags
      int evFlag = 0;
      if (collision.isInelGt0()) {
        evFlag = 1;
      }

      auto v0sThisCollision = V0s.sliceBy(perCol, collision.globalIndex());
      for (auto& v0 : v0sThisCollision) { // loop over V0s

        if (!v0.has_mcParticle()) {
          continue;
        }
        auto v0mcparticle = v0.mcParticle();

        if (abs(v0mcparticle.y()) > 0.5f) {
          continue;
        }

        // Highest numerator of efficiency
        if (v0mcparticle.isPhysicalPrimary()) {
          if (v0mcparticle.pdgCode() == 310) {
            registry.fill(HIST("Reconstructed_MCRecoColl_INEL_K0Short"), v0mcparticle.pt(), mcCollision.centFT0M()); // K0s
            if (evFlag == 1) {
              registry.fill(HIST("Reconstructed_MCRecoColl_INELgt0_K0Short"), v0mcparticle.pt(), mcCollision.centFT0M()); // K0s
            }
          }
          if (v0mcparticle.pdgCode() == 3122) {
            registry.fill(HIST("Reconstructed_MCRecoColl_INEL_Lambda"), v0mcparticle.pt(), mcCollision.centFT0M()); // Lambda
            if (evFlag == 1) {
              registry.fill(HIST("Reconstructed_MCRecoColl_INELgt0_Lambda"), v0mcparticle.pt(), mcCollision.centFT0M()); // Lambda
            }
          }
          if (v0mcparticle.pdgCode() == -3122) {
            registry.fill(HIST("Reconstructed_MCRecoColl_INEL_AntiLambda"), v0mcparticle.pt(), mcCollision.centFT0M()); // AntiLambda
            if (evFlag == 1) {
              registry.fill(HIST("Reconstructed_MCRecoColl_INELgt0_AntiLambda"), v0mcparticle.pt(), mcCollision.centFT0M()); // AntiLambda
            }
          }
        }

        int lPDG = 0;
        bool isprimary = false;
        if (TMath::Abs(v0mcparticle.pdgCode()) == 310 || TMath::Abs(v0mcparticle.pdgCode()) == 3122) {
          lPDG = v0mcparticle.pdgCode();
          isprimary = v0mcparticle.isPhysicalPrimary();
        }

        int posITSNhits = 0, negITSNhits = 0;
        for (unsigned int i = 0; i < 7; i++) {
          if (v0.posTrack_as<DauTracksMC>().itsClusterMap() & (1 << i)) {
            posITSNhits++;
          }
          if (v0.negTrack_as<DauTracksMC>().itsClusterMap() & (1 << i)) {
            negITSNhits++;
          }
        }

        float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
        float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0Bar;
        float ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

        if (v0.v0radius() > v0radius &&
            v0.v0cosPA() > v0cospa &&
            TMath::Abs(v0.posTrack_as<DauTracksMC>().eta()) < etadau &&
            TMath::Abs(v0.negTrack_as<DauTracksMC>().eta()) < etadau // &&
        ) {

          // Fill table
          myv0s(v0.globalIndex(), v0.pt(), v0.yLambda(), v0.yK0Short(),
                v0.mLambda(), v0.mAntiLambda(), v0.mK0Short(),
                v0.v0radius(), v0.v0cosPA(),
                v0.dcapostopv(), v0.dcanegtopv(), v0.dcaV0daughters(),
                v0.posTrack_as<DauTracksMC>().eta(), v0.negTrack_as<DauTracksMC>().eta(),
                v0.posTrack_as<DauTracksMC>().phi(), v0.negTrack_as<DauTracksMC>().phi(),
                posITSNhits, negITSNhits, ctauLambda, ctauAntiLambda, ctauK0s,
                v0.negTrack_as<DauTracksMC>().tpcNSigmaPr(), v0.posTrack_as<DauTracksMC>().tpcNSigmaPr(),
                v0.negTrack_as<DauTracksMC>().tpcNSigmaPi(), v0.posTrack_as<DauTracksMC>().tpcNSigmaPi(),
                v0.negTrack_as<DauTracksMC>().tofNSigmaPr(), v0.posTrack_as<DauTracksMC>().tofNSigmaPr(),
                v0.negTrack_as<DauTracksMC>().tofNSigmaPi(), v0.posTrack_as<DauTracksMC>().tofNSigmaPi(),
                v0.posTrack_as<DauTracksMC>().hasTOF(), v0.negTrack_as<DauTracksMC>().hasTOF(), lPDG, isprimary,
                mcCollision.centFT0M(), cent, evFlag, v0.alpha(), v0.qtarm());
        }
      }

      // Generated particles
      const auto particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache1);

      for (auto& mcParticle : particlesInCollision) {
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }

        if (abs(mcParticle.y()) > 0.5f) {
          continue;
        }

        if (mcParticle.pdgCode() == 310) {
          registry.fill(HIST("Generated_MCRecoCollCheck_INEL_K0Short"), mcParticle.pt(), mcCollision.centFT0M()); // K0s
          if (evFlag == 1) {
            registry.fill(HIST("Generated_MCRecoCollCheck_INELgt0_K0Short"), mcParticle.pt(), mcCollision.centFT0M()); // K0s
          }
        }
        if (mcParticle.pdgCode() == 3122) {
          registry.fill(HIST("Generated_MCRecoCollCheck_INEL_Lambda"), mcParticle.pt(), mcCollision.centFT0M()); // Lambda
          if (evFlag == 1) {
            registry.fill(HIST("Generated_MCRecoCollCheck_INELgt0_Lambda"), mcParticle.pt(), mcCollision.centFT0M()); // Lambda
          }
        }
        if (mcParticle.pdgCode() == -3122) {
          registry.fill(HIST("Generated_MCRecoCollCheck_INEL_AntiLambda"), mcParticle.pt(), mcCollision.centFT0M()); // AntiLambda
          if (evFlag == 1) {
            registry.fill(HIST("Generated_MCRecoCollCheck_INELgt0_AntiLambda"), mcParticle.pt(), mcCollision.centFT0M()); // AntiLambda
          }
        }
      }
    }
  }
  PROCESS_SWITCH(LfV0qaanalysis, processMCReco, "Process MC Reco", false);

  void processMCGen(soa::Join<aod::McCollisions, aod::McCentFT0Ms>::iterator const& mcCollision,
                    aod::McParticles const& mcParticles,
                    soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::PVMults>> const& collisions)
  {
    //====================================
    //===== Event Loss Denominator =======
    //====================================

    registry.fill(HIST("hNEventsMCGen"), 0.5);

    if (TMath::Abs(mcCollision.posZ()) > cutzvertex) {
      return;
    }
    registry.fill(HIST("hNEventsMCGen"), 1.5);
    registry.fill(HIST("hCentFT0M_GenColl_MC"), mcCollision.centFT0M());

    bool isINELgt0true = false;

    if (pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      isINELgt0true = true;
      registry.fill(HIST("hNEventsMCGen"), 2.5);
      registry.fill(HIST("hCentFT0M_GenColl_MC_INELgt0"), mcCollision.centFT0M());
    }

    //=====================================
    //===== Signal Loss Denominator =======
    //=====================================

    for (auto& mcParticle : mcParticles) {

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      if (abs(mcParticle.y()) > 0.5f) {
        continue;
      }

      if (mcParticle.pdgCode() == 310) {
        registry.fill(HIST("Generated_MCGenColl_INEL_K0Short"), mcParticle.pt(), mcCollision.centFT0M()); // K0s
        if (isINELgt0true) {
          registry.fill(HIST("Generated_MCGenColl_INELgt0_K0Short"), mcParticle.pt(), mcCollision.centFT0M()); // K0s
        }
      }
      if (mcParticle.pdgCode() == 3122) {
        registry.fill(HIST("Generated_MCGenColl_INEL_Lambda"), mcParticle.pt(), mcCollision.centFT0M()); // Lambda
        if (isINELgt0true) {
          registry.fill(HIST("Generated_MCGenColl_INELgt0_Lambda"), mcParticle.pt(), mcCollision.centFT0M()); // Lambda
        }
      }
      if (mcParticle.pdgCode() == -3122) {
        registry.fill(HIST("Generated_MCGenColl_INEL_AntiLambda"), mcParticle.pt(), mcCollision.centFT0M()); // AntiLambda
        if (isINELgt0true) {
          registry.fill(HIST("Generated_MCGenColl_INELgt0_AntiLambda"), mcParticle.pt(), mcCollision.centFT0M()); // AntiLambda
        }
      }
    }

    int recoCollIndex_INEL = 0;
    int recoCollIndex_INELgt0 = 0;
    for (auto& collision : collisions) { // loop on reconstructed collisions

      //=====================================
      //====== Event Split Numerator ========
      //=====================================

      registry.fill(HIST("hNEventsMCReco"), 0.5);
      if (!AcceptEvent(collision)) {
        continue;
      }
      registry.fill(HIST("hNEvents"), 8.5);
      registry.fill(HIST("hNEventsMCReco"), 1.5);
      registry.fill(HIST("hCentFT0M_RecoColl_MC"), mcCollision.centFT0M());

      recoCollIndex_INEL++;

      if (collision.isInelGt0() && isINELgt0true) {
        registry.fill(HIST("hNEventsMCReco"), 2.5);
        registry.fill(HIST("hCentFT0M_RecoColl_MC_INELgt0"), mcCollision.centFT0M());

        recoCollIndex_INELgt0++;
      }

      //=====================================
      //======== Sgn Split Numerator ========
      //=====================================

      for (auto& mcParticle : mcParticles) {

        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }

        if (abs(mcParticle.y()) > 0.5f) {
          continue;
        }

        if (mcParticle.pdgCode() == 310) {
          registry.fill(HIST("Generated_MCRecoColl_INEL_K0Short"), mcParticle.pt(), mcCollision.centFT0M()); // K0s
          if (recoCollIndex_INELgt0 > 0) {
            registry.fill(HIST("Generated_MCRecoColl_INELgt0_K0Short"), mcParticle.pt(), mcCollision.centFT0M()); // K0s
          }
        }
        if (mcParticle.pdgCode() == 3122) {
          registry.fill(HIST("Generated_MCRecoColl_INEL_Lambda"), mcParticle.pt(), mcCollision.centFT0M()); // Lambda
          if (recoCollIndex_INELgt0 > 0) {
            registry.fill(HIST("Generated_MCRecoColl_INELgt0_Lambda"), mcParticle.pt(), mcCollision.centFT0M()); // Lambda
          }
        }
        if (mcParticle.pdgCode() == -3122) {
          registry.fill(HIST("Generated_MCRecoColl_INEL_AntiLambda"), mcParticle.pt(), mcCollision.centFT0M()); // AntiLambda
          if (recoCollIndex_INELgt0 > 0) {
            registry.fill(HIST("Generated_MCRecoColl_INELgt0_AntiLambda"), mcParticle.pt(), mcCollision.centFT0M()); // AntiLambda
          }
        }
      }
    }

    // From now on keep only mc collisions with at least one reconstructed collision (INEL)
    if (recoCollIndex_INEL < 1) {
      return;
    }

    //=====================================
    //====== Event Loss Numerator =========
    //=====================================

    registry.fill(HIST("hNEventsMCGenReco"), 0.5);
    registry.fill(HIST("hCentFT0M_GenRecoColl_MC"), mcCollision.centFT0M());

    if (recoCollIndex_INELgt0 > 0) {
      registry.fill(HIST("hNEventsMCGenReco"), 1.5);
      registry.fill(HIST("hCentFT0M_GenRecoColl_MC_INELgt0"), mcCollision.centFT0M());
    }

    //=====================================
    //===== Signal Loss Numerator =========
    //=====================================

    for (auto& mcParticle : mcParticles) {

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }

      if (abs(mcParticle.y()) > 0.5f) {
        continue;
      }

      if (mcParticle.pdgCode() == 310) {
        registry.fill(HIST("Generated_MCGenRecoColl_INEL_K0Short"), mcParticle.pt(), mcCollision.centFT0M()); // K0s
        if (recoCollIndex_INELgt0 > 0) {
          registry.fill(HIST("Generated_MCGenRecoColl_INELgt0_K0Short"), mcParticle.pt(), mcCollision.centFT0M()); // K0s
        }
      }
      if (mcParticle.pdgCode() == 3122) {
        registry.fill(HIST("Generated_MCGenRecoColl_INEL_Lambda"), mcParticle.pt(), mcCollision.centFT0M()); // Lambda
        if (recoCollIndex_INELgt0 > 0) {
          registry.fill(HIST("Generated_MCGenRecoColl_INELgt0_Lambda"), mcParticle.pt(), mcCollision.centFT0M()); // Lambda
        }
      }
      if (mcParticle.pdgCode() == -3122) {
        registry.fill(HIST("Generated_MCGenRecoColl_INEL_AntiLambda"), mcParticle.pt(), mcCollision.centFT0M()); // AntiLambda
        if (recoCollIndex_INELgt0 > 0) {
          registry.fill(HIST("Generated_MCGenRecoColl_INELgt0_AntiLambda"), mcParticle.pt(), mcCollision.centFT0M()); // AntiLambda
        }
      }
    }
  }
  PROCESS_SWITCH(LfV0qaanalysis, processMCGen, "Process MC", false);
};

struct LfMyV0s {

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hPt", "hPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}});
    registry.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.4f, 0.6f}}});
    registry.add("hMassVsPtK0Short", "hMassVsPtK0Short", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 0.4f, 0.6f}}});
    registry.add("V0Radius", "V0Radius", {HistType::kTH1D, {{100, 0.0f, 20.0f}}});
    registry.add("CosPA", "CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("V0DCANegToPV", "V0DCANegToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("V0DCAPosToPV", "V0DCAPosToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("V0DCAV0Daughters", "V0DCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("CtauK0s", "CtauK0s", {HistType::kTH1F, {{150, 0.0f, 30.0f}}});
    registry.add("CtauLambda", "CtauLambda", {HistType::kTH1F, {{200, 0.0f, 40.0f}}});
    registry.add("CtauAntiLambda", "CtauAntiLambda", {HistType::kTH1F, {{200, 0.0f, 40.0f}}});
    registry.add("TPCNSigmaPosPi", "TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaNegPi", "TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaPosPr", "TPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaNegPr", "TPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("PosITSHits", "PosITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
    registry.add("NegITSHits", "NegITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
    registry.add("multft0m", "multft0m", {HistType::kTH1F, {{100, 0.f, 100.f}}});
  }

  void process(aod::MyV0Candidates const& myv0s)
  {
    for (auto& candidate : myv0s) {

      registry.fill(HIST("hMassLambda"), candidate.masslambda());
      registry.fill(HIST("hPt"), candidate.v0pt());
      registry.fill(HIST("hMassVsPtLambda"), candidate.v0pt(), candidate.masslambda());
      registry.fill(HIST("hMassAntiLambda"), candidate.massantilambda());
      registry.fill(HIST("hMassVsPtAntiLambda"), candidate.v0pt(), candidate.massantilambda());
      registry.fill(HIST("hMassK0Short"), candidate.massk0short());
      registry.fill(HIST("hMassVsPtK0Short"), candidate.v0pt(), candidate.massk0short());
      registry.fill(HIST("V0Radius"), candidate.v0radius());
      registry.fill(HIST("CosPA"), candidate.v0cospa());
      registry.fill(HIST("V0DCANegToPV"), candidate.v0dcanegtopv());
      registry.fill(HIST("V0DCAPosToPV"), candidate.v0dcapostopv());
      registry.fill(HIST("V0DCAV0Daughters"), candidate.v0dcav0daughters());
      registry.fill(HIST("CtauK0s"), candidate.ctauk0short());
      registry.fill(HIST("CtauLambda"), candidate.ctaulambda());
      registry.fill(HIST("CtauAntiLambda"), candidate.ctauantilambda());
      registry.fill(HIST("TPCNSigmaPosPi"), candidate.ntpcsigmapospi());
      registry.fill(HIST("TPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
      registry.fill(HIST("TPCNSigmaPosPr"), candidate.ntpcsigmapospr());
      registry.fill(HIST("TPCNSigmaNegPr"), candidate.ntpcsigmanegpr());
      registry.fill(HIST("PosITSHits"), candidate.v0positshits());
      registry.fill(HIST("NegITSHits"), candidate.v0negitshits());
      registry.fill(HIST("multft0m"), candidate.multft0m());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto w = WorkflowSpec{adaptAnalysisTask<LfV0qaanalysis>(cfgc)};
  if (cfgc.options().get<int>("add-fill")) {
    w.push_back(adaptAnalysisTask<LfMyV0s>(cfgc));
  }
  return w;
}
