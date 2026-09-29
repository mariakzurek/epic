// SPDX-License-Identifier: LGPL-3.0-or-later
#include "DRPrototypeFastSimModelOpFiber.h"

#include <G4FastStep.hh>
#include <G4FastTrack.hh>
#include <G4LogicalVolume.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4OpticalPhoton.hh>
#include <G4ParticleDefinition.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4TouchableHandle.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserTrackInformation.hh>
#include <Randomize.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>

namespace {
// Tunable Phase-2b starting values.  The 0.97 cut is stricter than the rejected
// 0.95 proposal: it is deliberately well inside the core's guided cone (core
// radii 0.485 mm, n(core)/n(clad) = 1.59/1.49 for S and 1.49/1.42 for C).
static constexpr G4double AXIAL_CONFIDENCE_THRESHOLD = 0.97;
// Five independent fast-sim encounters provide conservative evidence of
// established guided propagation rather than trusting the first sample.
static constexpr int MIN_AXIAL_SAMPLE_COUNT = 5;
// Geometry has 1.0 mm OD and 0.485 mm core radius; 5 core diameters (4.85 mm,
// rounded up) exceeds several transverse pitches (1.22/1.35 mm) and protects
// the normal end-window boundary treatment.  Tunable after Phase-3 validation.
static constexpr G4double END_SAFETY_MARGIN_MM = 5.0 * mm;
// Two absorption lengths retain exp(-2) survival while keeping every jump short.
// Tunable after Phase-3 attenuation/timing validation.
static constexpr G4double N_ABSORPTION_LENGTHS = 2.0;
// DRPrototype fiber half-length is module_dz/2 - window thickness = 499 mm.
static constexpr G4double FIBER_HALF_LENGTH_MM = 499.0 * mm;

class DRFiberTrackInfo final : public G4VUserTrackInformation {
 public:
  // This is an axial-consistency proxy for guided propagation, not a verified
  // TIR count. True per-boundary TIR detection is deferred pending framework study.
  int consecutiveAxialSampleCount{0};
  int maxAxialSampleCount{0};
  const G4VPhysicalVolume* lastPhysicalVolume{nullptr};
  G4ThreeVector lastAxialDirection{};
  bool everEnteredFiber{false};
};

const G4VPhysicalVolume* currentPhysicalVolume(const G4Track& track) {
  const auto& touchable = track.GetTouchableHandle();
  return touchable() == nullptr ? nullptr : touchable->GetVolume();
}

bool isCoreVolume(const G4VPhysicalVolume* volume) {
  if (volume == nullptr || volume->GetLogicalVolume() == nullptr) return false;
  const auto& name = volume->GetLogicalVolume()->GetName();
  return name == "DRPrototype_fiber_core_S" ||
         name == "DRPrototype_fiber_core_C";
}

DRFiberTrackInfo* fiberInfo(G4Track& track, DRPrototypeFastSimCounters& counters) {
  auto* userInfo = track.GetUserInformation();
  if (userInfo != nullptr) {
    auto* info = dynamic_cast<DRFiberTrackInfo*>(userInfo);
    if (info == nullptr) {
      counters.foreignTrackInfoConflict();
      return nullptr;
    }
    info->everEnteredFiber = true;
    return info;
  }

  auto* info = new DRFiberTrackInfo;
  // G4Track owns this object and deletes it with the track (G4 11.3.2).
  track.SetUserInformation(info);
  info->everEnteredFiber = true;
  return info;
}

G4double propertyValue(const G4MaterialPropertiesTable* table, const char* name,
                      G4double energy, G4double fallback) {
  if (table == nullptr) return fallback;
  auto* property = table->GetProperty(name);
  return property == nullptr ? fallback : property->Value(energy);
}
}  // namespace

dd4hep::sim::DRPrototypeFastSimModelOpFiber::DRPrototypeFastSimModelOpFiber(
    dd4hep::sim::Geant4Context* context, const std::string& name)
    : dd4hep::sim::Geant4FSShowerModel<DRPrototypeFastSimCounters>(context, name) {
  m_trace = std::getenv("DR_FASTSIM_TRACE") != nullptr;
  m_regionName = "DRPrototypeFastFiberRegion";
  m_applicablePartNames = {"opticalphoton"};
  m_enable = true;
  declareProperty("Trace", m_trace);
}

dd4hep::sim::DRPrototypeFastSimModelOpFiber::~DRPrototypeFastSimModelOpFiber() {
  if (m_trace) {
    info("DR fast-sim counters: attempted=%llu accepted=%llu fallback=%llu absorbed=%llu WLS-fallback=%llu foreign_track_info_conflict=%llu axial-cosine-failure=%llu opposite-direction-reversal=%llu PV-pointer-change=%llu end-margin-rejection=%llu WLS-rejection=%llu non-core-volume-rejection=%llu max-count-reached={0:%llu,1:%llu,2:%llu,3:%llu,4:%llu,>=5-rejected:%llu,accepted:%llu}",
         static_cast<unsigned long long>(locals.attemptedCount()),
         static_cast<unsigned long long>(locals.acceptedCount()),
         static_cast<unsigned long long>(locals.fallbackCount()),
          static_cast<unsigned long long>(locals.absorbedCount()),
          static_cast<unsigned long long>(locals.wlsFallbackCount()),
          static_cast<unsigned long long>(locals.foreignTrackInfoConflictCount()),
          static_cast<unsigned long long>(locals.axialCosineFailureCount()),
          static_cast<unsigned long long>(locals.oppositeDirectionReversalCount()),
          static_cast<unsigned long long>(locals.physicalVolumeChangeCount()),
          static_cast<unsigned long long>(locals.endMarginRejectionCount()),
          static_cast<unsigned long long>(locals.wlsRejectionCount()),
          static_cast<unsigned long long>(locals.nonCoreVolumeRejectionCount()),
          static_cast<unsigned long long>(locals.maxCountReachedCount(0)),
          static_cast<unsigned long long>(locals.maxCountReachedCount(1)),
          static_cast<unsigned long long>(locals.maxCountReachedCount(2)),
          static_cast<unsigned long long>(locals.maxCountReachedCount(3)),
          static_cast<unsigned long long>(locals.maxCountReachedCount(4)),
          static_cast<unsigned long long>(locals.maxCountReachedCount(5)),
           static_cast<unsigned long long>(locals.maxCountReachedCount(6)));
     const auto nonCoreCount = locals.nonCoreVolumeRejectionCount();
     const auto overflow = locals.nonCoreVolumeOtherCount();
      info("DR fast-sim non-core volume samples: first-8-distinct-names-encountered (sorted-by-count-among-retained) overflow=%llu overflow-fraction=%.6f near-miss-core-substring-seen=%s",
           static_cast<unsigned long long>(overflow),
           nonCoreCount == 0 ? 0.0 : static_cast<double>(overflow) / static_cast<double>(nonCoreCount),
           locals.nearMissCoreSubstringSeenInNonCoreSamples() ? "yes" : "no");
      std::array<bool, 8> emitted{};
      // The fixed table retains the first 8 distinct names encountered; output is sorted by count among those retained, not globally.
     for (std::size_t rank = 0; rank < locals.nonCoreVolumeSampleCount(); ++rank) {
       std::size_t best = locals.nonCoreVolumeSampleCount();
       for (std::size_t i = 0; i < locals.nonCoreVolumeSampleCount(); ++i) {
         if (!emitted[i] && (best == locals.nonCoreVolumeSampleCount() ||
                             locals.nonCoreVolumeSampleHits(i) > locals.nonCoreVolumeSampleHits(best))) {
           best = i;
         }
       }
       emitted[best] = true;
       const auto count = locals.nonCoreVolumeSampleHits(best);
       const auto fraction = nonCoreCount == 0
                                 ? 0.0
                                 : static_cast<double>(count) / static_cast<double>(nonCoreCount);
       info("DR fast-sim non-core volume sample: rank=%llu name=%s count=%llu fraction=%.6f",
            static_cast<unsigned long long>(rank + 1),
            locals.nonCoreVolumeSampleName(best).c_str(),
            static_cast<unsigned long long>(count), fraction);
     }
   }
 }

bool dd4hep::sim::DRPrototypeFastSimModelOpFiber::check_applicability(
    const G4ParticleDefinition& particle) {
  if (&particle != G4OpticalPhoton::OpticalPhotonDefinition() || !m_enable) {
    locals.fallback();
    return false;
  }
  return true;
}

bool dd4hep::sim::DRPrototypeFastSimModelOpFiber::check_trigger(const G4FastTrack& fastTrack) {
  locals.attempted();
  const G4Track* track = fastTrack.GetPrimaryTrack();
  if (track == nullptr || track->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
    locals.fallback();
    return false;
  }

   const G4VPhysicalVolume* volume = currentPhysicalVolume(*track);
   if (!isCoreVolume(volume)) {
     locals.nonCoreVolumeRejection();
     // Phase-3 diagnostic for zero accepted jumps; trace remains default-off.
     if (m_trace) {
       locals.nonCoreVolumeSample(volume == nullptr ? "<null>" : volume->GetName());
     }
     locals.fallback();
     return false;
   }
  auto* info = fiberInfo(*const_cast<G4Track*>(track), locals);
  if (info == nullptr) {
    locals.fallback();
    return false;
  }
  if (info->lastPhysicalVolume != volume) {
    // A track-death finalization hook is unavailable; reset events are the
    // bounded high-water-mark observation points captured here.
    locals.physicalVolumeChange();
    locals.maxCountReached(info->maxAxialSampleCount, false);
    // Samples are only comparable within one physical core volume. This also
    // prevents stale state across core/cladding transitions and re-entry.
    info->consecutiveAxialSampleCount = 0;
    info->maxAxialSampleCount = 0;
    info->lastAxialDirection = G4ThreeVector();
    info->lastPhysicalVolume = volume;
  }
  // FastTrack's local accessors intentionally preserve the established local
  // coordinate-frame handling; this model only advances the local fiber axis.
  const G4ThreeVector localDirection = fastTrack.GetPrimaryTrackLocalDirection();
  const G4double axialCosine = localDirection.z();
  const bool axialCosineFailure = std::abs(axialCosine) < AXIAL_CONFIDENCE_THRESHOLD;
  const bool oppositeDirectionReversal =
      info->everEnteredFiber && info->lastAxialDirection.mag2() > 0.0 &&
      localDirection.dot(info->lastAxialDirection) <= 0.0;
  if (axialCosineFailure || oppositeDirectionReversal) {
    if (axialCosineFailure) locals.axialCosineFailure();
    if (oppositeDirectionReversal) locals.oppositeDirectionReversal();
    locals.maxCountReached(info->maxAxialSampleCount, false);
    info->consecutiveAxialSampleCount = 0;
    info->maxAxialSampleCount = 0;
    info->lastAxialDirection = localDirection;
    locals.fallback();
    return false;
  }
  info->lastAxialDirection = localDirection;
  ++info->consecutiveAxialSampleCount;
  info->maxAxialSampleCount =
      std::max(info->maxAxialSampleCount, info->consecutiveAxialSampleCount);
  if (info->consecutiveAxialSampleCount < MIN_AXIAL_SAMPLE_COUNT) {
    locals.fallback();
    return false;
  }

  const G4double z = fastTrack.GetPrimaryTrackLocalPosition().z();
  if (z <= -FIBER_HALF_LENGTH_MM + END_SAFETY_MARGIN_MM ||
      z >= FIBER_HALF_LENGTH_MM - END_SAFETY_MARGIN_MM) {
    locals.endMarginRejection();
    locals.maxCountReached(info->maxAxialSampleCount, false);
    locals.fallback();
    return false;
  }

  const auto* material = track->GetMaterial();
  const auto* table = material == nullptr ? nullptr : material->GetMaterialPropertiesTable();
  if (table == nullptr || table->GetProperty("WLSABSLENGTH") != nullptr ||
      table->GetProperty("WLSCOMPONENT") != nullptr) {
    locals.wlsRejection();
    locals.maxCountReached(info->maxAxialSampleCount, false);
    locals.wlsFallback();
    locals.fallback();
    return false;
  }
  locals.maxCountReached(info->maxAxialSampleCount, true);
  return true;
}

void dd4hep::sim::DRPrototypeFastSimModelOpFiber::modelShower(
    const G4FastTrack& fastTrack, G4FastStep& step) {
  const G4Track* track = fastTrack.GetPrimaryTrack();
  const G4double energy = track->GetKineticEnergy();
  const G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  const G4ThreeVector localDirection = fastTrack.GetPrimaryTrackLocalDirection();
  const G4double end = localDirection.z() >= 0.0 ? FIBER_HALF_LENGTH_MM : -FIBER_HALF_LENGTH_MM;
  const G4double safeDistance = std::abs(end - localPosition.z()) - END_SAFETY_MARGIN_MM;
  const auto* table = track->GetMaterial()->GetMaterialPropertiesTable();
  const G4double absLength = propertyValue(table, "ABSLENGTH", energy,
                                            std::numeric_limits<G4double>::infinity());
  if (!(safeDistance > 0.0) || !(absLength > 0.0)) {
    // This should be unreachable after check_trigger; define a harmless state.
    step.ProposePrimaryTrackFinalPosition(track->GetPosition(), false);
    step.ProposePrimaryTrackFinalTime(track->GetGlobalTime());
    step.ProposePrimaryTrackFinalProperTime(track->GetProperTime());
    step.ProposePrimaryTrackFinalMomentumDirection(track->GetMomentumDirection(), false);
    step.ProposePrimaryTrackFinalPolarization(track->GetPolarization(), false);
    step.ProposePrimaryTrackPathLength(0.0);
    locals.fallback();
    return;
  }

  const G4double axialCosine = localDirection.z();
  const G4double maxOpticalPath = safeDistance / std::abs(axialCosine);
  const G4double jump = std::min(maxOpticalPath, N_ABSORPTION_LENGTHS * absLength);
  const G4double survival = std::exp(-jump / absLength);
  const G4double u = G4UniformRand();
  const G4double opticalPathLength = u < survival ? jump :
      -absLength * std::log(1.0 - G4UniformRand() * (1.0 - survival));
  // Conditional absorption inverse CDF: F(d|d<=L)=(1-exp(-d/a))/(1-exp(-L/a));
  // solving F=U gives d=-a ln(1-U(1-exp(-L/a))).
  const G4double axialDisplacement = opticalPathLength * std::abs(axialCosine);
  const G4double n = propertyValue(table, "RINDEX", energy, 1.0);
  const G4double groupVelocity = propertyValue(table, "GROUPVEL", energy, c_light / n);
  // Geant4 may auto-calculate GROUPVEL from RINDEX; c/n is only the fallback
  // for materials where GROUPVEL is genuinely unavailable.
  const G4double finalTime = track->GetGlobalTime() + opticalPathLength / groupVelocity;
  // The idealized fast transport advances only on the fiber's local axis;
  // transverse drift and the original direction/polarization are intentionally
  // ignored/preserved, respectively, as required by the reference model.
  const G4ThreeVector finalLocalPosition =
      localPosition + G4ThreeVector(0.0, 0.0,
                                    localDirection.z() >= 0.0 ? axialDisplacement :
                                                               -axialDisplacement);

  step.ProposePrimaryTrackFinalPosition(finalLocalPosition, true);
  step.ProposePrimaryTrackFinalTime(finalTime);
  // Optical photons have degenerate proper time; leave the framework's current
  // value unchanged rather than advancing it with the laboratory travel time.
  step.ProposePrimaryTrackFinalProperTime(track->GetProperTime());
  step.ProposePrimaryTrackFinalMomentumDirection(localDirection, true);
  step.ProposePrimaryTrackFinalPolarization(fastTrack.GetPrimaryTrackLocalPolarization(), true);
  step.ProposePrimaryTrackFinalKineticEnergy(energy);
  step.ProposePrimaryTrackPathLength(opticalPathLength);
  locals.accepted();
  if (opticalPathLength != jump) {
    locals.absorbed();
    // Match G4OpAbsorption::PostStepDoIt(): absorbed photon energy is deposited.
    step.ProposeTotalEnergyDeposited(energy);
    step.KillPrimaryTrack();
  }
}
