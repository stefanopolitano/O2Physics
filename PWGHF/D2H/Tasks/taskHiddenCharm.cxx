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

/// \file taskHiddenCharm.cxx
/// \brief Analysis task for hidden charm hadrons at midrapidity
///
/// \author A. Palasciano, <antonio.palasciano@cern.ch>, INFN Bari
/// \author S. Politanò <stefano.politano@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"

#include "Common/Core/trackUtilities.h"
#include "DetectorsBase/Propagator.h"
#include "DCAFitter/DCAFitterN.h"
#include <DetectorsVertexing/PVertexer.h>
#include <Framework/AnalysisDataModel.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <array>
#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;

enum TrackType : uint8_t {
  Pion = 0,
  Kaon,
  Proton
};

struct HfTaskHiddenCharm {
  o2::vertexing::DCAFitterN<2> df2; // 2-prong vertex fitter
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::base::MatLayerCylSet* lut{};
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  Configurable<int> centEstimator{"centEstimator", 0, "Centrality estimation (None: 0, FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4, NTracksPV: 5, FT0CVariant2: 6)"};
  Configurable<float> centralityMin{"centralityMin", 0.f, "Minimum accepted centrality"};
  Configurable<float> centralityMax{"centralityMax", 100.f, "Maximum accepted centrality"};
  Configurable<bool> fillOnlyUnlikeSign{"fillOnlyUnlikeSign", true, "Fill only unlike-sign proton pairs"};
  Configurable<bool> fillOnlyLikeSign{"fillOnlyLikeSign", true, "Fill only like-sign proton pairs"};

  SliceCache cache;

  using SelectedCollisionsPP = aod::HfRedCollisions;
  using SelectedCollisionsPbPb = soa::Join<aod::Collisions,
                                           aod::HfRedCollisions,
                                           aod::CentFT0As,
                                           aod::CentFT0Cs,
                                           aod::CentFT0Ms,
                                           aod::CentFV0As,
                                           aod::CentNTPVs,
                                           aod::CentFT0CVariant2s>;
  Partition<aod::HcSelTracks> selectedProtons = aod::hf_track_vars_reduced::trackType == static_cast<uint8_t>(TrackType::Proton);
  Partition<aod::HcSelTracks> selectedPions = aod::hf_track_vars_reduced::trackType == static_cast<uint8_t>(TrackType::Pion);
  Partition<aod::HcSelTracks> selectedKaons = aod::hf_track_vars_reduced::trackType == static_cast<uint8_t>(TrackType::Kaon);

struct : ConfigurableGroup {
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {1400, 2.8, 4.2}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 10.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {100, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisSign{"thnConfigAxisSign", {2, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisDCA{"thnConfigAxisDCA", {400, -0.2, 0.2}, "DCA of eta_c candidate to primary vertex (cm)"};
  
  // Vertexing configuration
  Configurable<bool> storeDCA{"storeDCA", false, "Store DCA of the track to the primary vertex in the output table"};
  Configurable<bool> propagateToPCA{"propagateToPCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", false, "Minimise abs. distance rather than chi2"};
  Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<double> etacRadiusMax{"etacRadiusMax", 1., "Maximum radius of the secondary vertex for eta_c candidates"};
  Configurable<double> etaCMinCosPa{"etaCMinCosPa", 0.9, "Minimum cosine of the pointing angle for eta_c candidates"};
  } config;

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    const AxisSpec axisInvMass{config.thnConfigAxisInvMass, "M_{p#bar{p}} (GeV/#it{c}^{2})"};
    const AxisSpec axisPt{config.thnConfigAxisPt, "#it{p}_{T}^{p#bar{p}} (GeV/#it{c})"};
    const AxisSpec axisCent{config.thnConfigAxisCent, "Centrality"};
    const AxisSpec axisSign{config.thnConfigAxisSign, "q_{1} #times q_{2}"};
    const AxisSpec axisDCA{config.thnConfigAxisDCA, "DCA (cm)"};
    
    std::vector<AxisSpec> axes = {axisInvMass, axisPt, axisSign, axisCent};
    if (config.storeDCA) {
      axes.push_back(axisDCA);
    }
    registry.add("hSparseHiddenCharm", "Hidden-charm proton-pair candidates", HistType::kTHnSparseF, axes);
    registry.add("hPtVsInvMassLikeSign", "Hidden-charm LS M_{inv}", HistType::kTH2D, {axisInvMass, axisPt});
    registry.add("hPtVsInvMassUnlikeSign", "Hidden-charm proton-pair ULS", HistType::kTH2D, {axisInvMass, axisPt});
    registry.add("hPtVsInvMassAllSign", "Hidden-charm proton-pair LS+ULS M_{inv}", HistType::kTH2D, {axisInvMass, axisPt});
    registry.add("hDCAxyEtacToPVUnlikeSign", "DCA of eta_c candidate to primary vertex;DCA (cm);entries", HistType::kTH2D, {axisDCA, axisPt});
    registry.add("hDCAxyEtacToPVLikeSign", "DCA of eta_c candidate to primary vertex;DCA (cm);entries", HistType::kTH2D, {axisDCA, axisPt});

  }

  void configure()
  {
    // Vertexing
    df2.setPropagateToPCA(config.propagateToPCA);
    df2.setMaxR(config.maxR);
    df2.setMaxDZIni(config.maxDZIni);
    df2.setMinParamChange(config.minParamChange);
    df2.setMinRelChi2Change(config.minRelChi2Change);
    df2.setUseAbsDCA(config.useAbsDCA);
    df2.setWeightedFinalPCA(config.useWeightedFinalPCA);
  }

  float calculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
  {
    return std::sqrt((std::pow((pvY - Y) * Pz - (pvZ - Z) * Py, 2) + std::pow((pvX - X) * Pz - (pvZ - Z) * Px, 2) + std::pow((pvX - X) * Py - (pvY - Y) * Px, 2)) / (Px * Px + Py * Py + Pz * Pz));
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename TCollisions, typename TProtonIds>
  void fillEtac(TCollisions const& collision,
                aod::BCsWithTimestamps const&,
                TProtonIds const& protonIds)
  {
    float cent{-1.f};
    if constexpr (CentEstimator != o2::hf_centrality::CentralityEstimator::None) {
      cent = o2::hf_centrality::getCentralityColl(collision, centEstimator);
      if (cent < centralityMin || cent >= centralityMax) {
        return; // skip events outside the centrality range
      }
    }

    // Reduced collisions already store the magnetic field.
    df2.setBz(collision.bz());

    for (const auto& proton1 : protonIds) {
      for (const auto& proton2 : protonIds) {
        if (proton1.trackId() >= proton2.trackId()) {
          continue; // avoid double counting and self-pairs
        }
        const int sign = (proton1.sign() * proton2.sign() > 0) ? 1 : -1;
        if ((sign == 1 && !fillOnlyLikeSign) || (sign == -1 && !fillOnlyUnlikeSign)) {
          continue;
        }
        std::array<float, 3> pVec1{proton1.px(), proton1.py(), proton1.pz()};
        std::array<float, 3> pVec2{proton2.px(), proton2.py(), proton2.pz()};

        // Vertexing
        float dcaEtaCToPv = -1.f;
        auto trackParVar1 = getTrackParCov(proton1);
        auto trackParVar2 = getTrackParCov(proton2);
        int nVtxFrom2ProngFitter = df2.process(trackParVar1, trackParVar2);
        if (nVtxFrom2ProngFitter > 0) { // should it be this or > 0 or are they equivalent
          // get primary vertex
          std::array<float, 3> primVtx = {collision.posX(), collision.posY(), collision.posZ()};
          // get secondary vertex
          const auto& secVtx2 = df2.getPCACandidate();
          std::array<float, 3> momEtaC = RecoDecay::sumOfVec(pVec1, pVec2);

          float radiusEtac = std::hypot(secVtx2[0], secVtx2[1]);
          if (radiusEtac > config.etacRadiusMax.value) {
            continue;
          }
          double cosPA = RecoDecay::cpa(primVtx, secVtx2, momEtaC);
          if (cosPA < config.etaCMinCosPa.value) {
            continue;
          }

          dcaEtaCToPv = calculateDCAStraightToPV(secVtx2[0], secVtx2[1], secVtx2[2], momEtaC[0], momEtaC[1], momEtaC[2], primVtx[0], primVtx[1], primVtx[2]);
        }
        float invMass = RecoDecay::m(std::array{pVec1, pVec2}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassProton});
        float ptEtac = RecoDecay::pt(RecoDecay::sumOfVec(pVec1, pVec2));
        if (config.storeDCA) {
          registry.fill(HIST("hSparseHiddenCharm"), invMass, ptEtac, sign, cent, dcaEtaCToPv);
        } else {
          registry.fill(HIST("hSparseHiddenCharm"), invMass, ptEtac, sign, cent);
        }
        registry.fill(HIST("hPtVsInvMassAllSign"), invMass, ptEtac);
        if (sign == 1) {
          registry.fill(HIST("hPtVsInvMassLikeSign"), invMass, ptEtac);
          registry.fill(HIST("hDCAxyEtacToPVLikeSign"), dcaEtaCToPv, ptEtac);
        } else if (sign == -1) {
          registry.fill(HIST("hDCAxyEtacToPVUnlikeSign"), dcaEtaCToPv, ptEtac);
          registry.fill(HIST("hPtVsInvMassUnlikeSign"), invMass, ptEtac);
        }
      }
    }
  }

  void processEtacPP(SelectedCollisionsPP::iterator const& collision,
                     aod::BCsWithTimestamps const& bcs,
                     aod::HcSelTracks const& /*tracks*/)
  {
    auto candProtons = selectedProtons->sliceByCached(aod::hf_track_index_reduced::hfRedCollisionId, collision.globalIndex(), cache);
    fillEtac<CentralityEstimator::None>(collision, bcs, candProtons);
  }
  PROCESS_SWITCH(HfTaskHiddenCharm, processEtacPP, "Process Etac candidates for pp", true);

  // void processEtacPbPb(SelectedCollisionsPbPb::iterator const& collisions,
  //                      aod::HcSelTracks const& protonIds)
  //{
  //   fillEtac(collisions, protonIds, true);
  // }
  // PROCESS_SWITCH(HfTaskHiddenCharm, processEtacPbPb, "Process Etac candidates for PbPb", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskHiddenCharm>(cfgc)};
}
