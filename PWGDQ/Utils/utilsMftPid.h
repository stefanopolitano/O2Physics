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

/// \file utilsMftPid.h
/// \brief Utilities for MFT PID in DQ analyses

#ifndef PWGDQ_UTILS_UTILSMFTPID_H_
#define PWGDQ_UTILS_UTILSMFTPID_H_

#include <algorithm> // std::upper_bound
#include <iterator>  // std::distance
#include <string>    //std::string

namespace o2::dq_mftpid
{

/// @brief Function to apply track-quality selections on mft and global tracks for pid studies
/// \param Track is the track candidate
/// \param maxDCA is the maximum acceptable dca in x and y
template <bool isMuontrack, typename Track>
bool isPidTrackSelected(Track const& track, float maxDCA) {
  if constexpr (isMuontrack) {
    return (track.has_matchMFTTrack() && track.trackType() == 0 && TMath::Abs(track.fwdDcaX()) < maxDCA && TMath::Abs(track.fwdDcaY()) < maxDCA);
  } else {
    return (TMath::Abs(track.fwdDcaX()) < maxDCA && TMath::Abs(track.fwdDcaY()) < maxDCA);
  }
}

/// @brief Function to apply track-quality selections on mft and global tracks for pid studies
/// \param muon1 is the first muon candidate
/// \param muon2 is the second muon candidate
/// \param maxDCA is the maximum acceptable dca in x and y
template <typename MuonP1>
bool isMuonPairsSelected(MuonP1 const& muon1, MuonP1 const& muon2, float maxDCA) {
  return (muon1.sign() == muon2.sign() || isPidTrackSelected<true>(muon1, maxDCA) || isPidTrackSelected<true>(muon2, maxDCA));
}

} // namespace o2::dq_mftpid

#endif // PWGDQ_UTILS_UTILSMFTPID_H_
