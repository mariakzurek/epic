// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef EPIC_DRPROTOTYPEFASTSIMMODELOPFIBER_H
#define EPIC_DRPROTOTYPEFASTSIMMODELOPFIBER_H

#include <DDG4/Factories.h>
#include <DDG4/Geant4FastSimShowerModel.h>
#include <DDG4/Geant4FastSimShowerModel.inl.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

class G4FastStep;
class G4FastTrack;
class G4ParticleDefinition;

struct DRPrototypeFastSimCounters {
  void attempted() { ++m_attempted; }
  void accepted() { ++m_accepted; }
  void fallback() { ++m_fallback; }
  void absorbed() { ++m_absorbed; }
  void wlsFallback() { ++m_wlsFallback; }
  void foreignTrackInfoConflict() { ++m_foreignTrackInfoConflict; }
  void axialCosineFailure() { ++m_axialCosineFailure; }
  void oppositeDirectionReversal() { ++m_oppositeDirectionReversal; }
  void physicalVolumeChange() { ++m_physicalVolumeChange; }
  void endMarginRejection() { ++m_endMarginRejection; }
  void wlsRejection() { ++m_wlsRejection; }
  void nonCoreVolumeRejection() { ++m_nonCoreVolumeRejection; }
  void nonCoreVolumeSample(const std::string& name) {
    if (name.find("fiber_core") != std::string::npos) m_nearMissCoreSubstringSeen = true;
    for (auto& sample : m_nonCoreVolumeSamples) {
      if (sample.name == name) {
        ++sample.count;
        return;
      }
    }
    if (m_nonCoreVolumeSampleCount < m_nonCoreVolumeSamples.size()) {
      auto& sample = m_nonCoreVolumeSamples[m_nonCoreVolumeSampleCount++];
      sample.name = name;
      sample.count = 1;
    } else {
      ++m_nonCoreVolumeOther;
    }
  }
  std::size_t nonCoreVolumeSampleCount() const { return m_nonCoreVolumeSampleCount; }
  const std::string& nonCoreVolumeSampleName(std::size_t index) const {
    return m_nonCoreVolumeSamples[index].name;
  }
  std::uint64_t nonCoreVolumeSampleHits(std::size_t index) const {
    return m_nonCoreVolumeSamples[index].count;
  }
  std::uint64_t nonCoreVolumeOtherCount() const { return m_nonCoreVolumeOther; }
  bool nearMissCoreSubstringSeenInNonCoreSamples() const {
    return m_nearMissCoreSubstringSeen;
  }
  void maxCountReached(int count, bool accepted) {
    const auto bucket = count < 5 ? count : (accepted ? 6 : 5);
    ++m_maxCountReached[static_cast<std::size_t>(bucket)];
  }
  std::uint64_t attemptedCount() const { return m_attempted; }
  std::uint64_t acceptedCount() const { return m_accepted; }
  std::uint64_t fallbackCount() const { return m_fallback; }
  std::uint64_t absorbedCount() const { return m_absorbed; }
  std::uint64_t wlsFallbackCount() const { return m_wlsFallback; }
  std::uint64_t foreignTrackInfoConflictCount() const { return m_foreignTrackInfoConflict; }
  std::uint64_t axialCosineFailureCount() const { return m_axialCosineFailure; }
  std::uint64_t oppositeDirectionReversalCount() const { return m_oppositeDirectionReversal; }
  std::uint64_t physicalVolumeChangeCount() const { return m_physicalVolumeChange; }
  std::uint64_t endMarginRejectionCount() const { return m_endMarginRejection; }
  std::uint64_t wlsRejectionCount() const { return m_wlsRejection; }
  std::uint64_t nonCoreVolumeRejectionCount() const { return m_nonCoreVolumeRejection; }
  std::uint64_t maxCountReachedCount(std::size_t bucket) const { return m_maxCountReached[bucket]; }
 private:
  struct NonCoreVolumeSample {
    std::string name;
    std::uint64_t count{0};
  };
  static constexpr std::size_t MAX_NON_CORE_VOLUME_SAMPLES = 8;

  std::uint64_t m_attempted{0};
  std::uint64_t m_accepted{0};
  std::uint64_t m_fallback{0};
  std::uint64_t m_absorbed{0};
  std::uint64_t m_wlsFallback{0};
  std::uint64_t m_foreignTrackInfoConflict{0};
  std::uint64_t m_axialCosineFailure{0};
  std::uint64_t m_oppositeDirectionReversal{0};
  std::uint64_t m_physicalVolumeChange{0};
  std::uint64_t m_endMarginRejection{0};
  std::uint64_t m_wlsRejection{0};
  std::uint64_t m_nonCoreVolumeRejection{0};
  std::array<NonCoreVolumeSample, MAX_NON_CORE_VOLUME_SAMPLES> m_nonCoreVolumeSamples{};
  std::size_t m_nonCoreVolumeSampleCount{0};
  std::uint64_t m_nonCoreVolumeOther{0};
  bool m_nearMissCoreSubstringSeen{false};
  std::uint64_t m_maxCountReached[7]{};
};

namespace dd4hep::sim {
class DRPrototypeFastSimModelOpFiber
  : public Geant4FSShowerModel<DRPrototypeFastSimCounters> {
public:
  DRPrototypeFastSimModelOpFiber(Geant4Context* context, const std::string& name);
  ~DRPrototypeFastSimModelOpFiber() override;

  bool check_applicability(const G4ParticleDefinition& particle) override;
  bool check_trigger(const G4FastTrack& track) override;
  void modelShower(const G4FastTrack& track, G4FastStep& step) override;

 private:
  bool m_trace{false};
};
}

DECLARE_GEANT4ACTION(DRPrototypeFastSimModelOpFiber)

#endif
