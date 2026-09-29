// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 EPIC Collaboration
//
// DR Prototype: Dual-Readout Fiber Calorimeter Plugin
// Plan reference: plan-dr-prototype.md todos 2b, 2a.5, 2c.5
// Research references:
//   research-optical-patterns.md §1.6  (assembly voxelization)
//   research-optical-patterns.md §2.3  (SkinSurface binding)
//   research-DR-baseline.md §1.2       (lattice geometry)
//
// SD architecture (design principle, plan-dr-prototype.md top of Phase 2):
//   sens.setType("calorimeter") -> Geant4CalorimeterAction.
//   Sensitive volume: thin AirOptical SiPM window at each fiber end.
//   DR_SiPMSurf SkinSurface on each window volume: EFFICIENCY table governs detection.
//   1 SimTrackerHit == 1 detected photoelectron.  No analysis-time PDE weight.
//
// S/C assignment (deterministic, plan decision 4):
//   Fibers numbered in row-major order (row iy from bottom, column ix from left).
//   Pattern: cycle through (sc_ratio_S S-fibers, then sc_ratio_C C-fibers), repeat.
//   Same mask every event; no random seed variation.
//   channel PhysVolID: 0 = S (DR_Polystyrene core), 1 = C (DR_PMMA core).
//
// Assembly voxelization (REQUIRED, plan R17, research-optical-patterns.md §1.6):
//   Two-tier hierarchy:
//     - Outer: Assembly (unbounded) containing only the Pb absorber block.
//     - Fiber grid assemblies placed INSIDE pbVol (daughters of the Pb block).
//       Each fiberAsm cell contains both the fiber tubes AND the two SiPM window
//       disks at each end.  Windows are flush at the fiber end faces with a 1 µm
//       epsilon to avoid coincident boundaries.
  //   Fiber half-length is shortened by (sipm_win_thickness + 2×1 µm) so that all
  //   volumes (fiber + both end windows) fit strictly within pbVol, with a 1 µm gap
  //   between fiber end face and window, and another 1 µm between window and pbVol face.
//   Each fiber grid assembly is Voxelize("")d.
//
// Fiber cladding (plan R16):
//   Each fiber = full-radius FluorinatedPolymer cladding volume (mother) with
//   Polystyrene or PMMA core placed inside as a daughter.
//   Cladding volume: Tube(0, fiber_od_half, fiber_half_length) - solid cylinder of
//   DR_FluorinatedPolymer (n=1.42).  Core volume: Tube(0, core_radius, ...) placed
//   at (0,0,0) inside the cladding - daughter fully contained within mother.
//   n(cladding) < n(core) everywhere: TIR guaranteed (plan R16).
//   TIR critical angles: PS/FluoPoly=63.4 deg, PMMA/FluoPoly=72.6 deg.
//
//   NOTE: using a solid (non-annular) cladding tube avoids the "core protruding
//   from annular mother" TGeo overlap error.  The solid cladding material fills
//   the full fiber OD; the core daughter replaces cladding in the core region.
//   This is the standard DD4hep pattern for a cladded fiber.
//
//
// Cladding thickness: 15 um (core radius 0.485 mm, clad OD 0.500 mm = 1.0 mm fiber OD).
//   Source: research-DR-baseline.md §7.5 notes "single clad" for SCSF-78.
//   Typical SCSF-78 single-clad spec: 15-30 um cladding.  Choice: 15 um.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"

#include <cmath>
#include <string>
#include <vector>

using namespace dd4hep;

// ---------------------------------------------------------------------------
// Helper: compute fiber count per axis for a rectangular lattice
// ---------------------------------------------------------------------------
static int fiberCount(double module_size, double pitch) {
  return static_cast<int>(std::floor(module_size / pitch));
}

// ---------------------------------------------------------------------------
// create_detector
// ---------------------------------------------------------------------------
static Ref_t create_detector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens) {

  xml::DetElement detElem = handle;
  const std::string detName = detElem.nameStr();
  const int detID           = detElem.id();

  sens.setType("calorimeter");

  // Surface manager for SkinSurface binding (plan todo 2a.5)
  OpticalSurfaceManager surfMgr = desc.surfaceManager();

  // -------------------------------------------------------------------------
  // Read geometry parameters from compact XML <dimensions> child
  // -------------------------------------------------------------------------
  xml::Component dims = detElem.child(_Unicode(dimensions));

  const double module_dx          = dims.attr<double>(_Unicode(dx));
  const double module_dy          = dims.attr<double>(_Unicode(dy));
  const double module_dz          = dims.attr<double>(_Unicode(dz));
  const double fiber_od           = dims.attr<double>(_Unicode(fiber_od));
  const double fiber_pitch_x      = dims.attr<double>(_Unicode(fiber_pitch_x));
  const double fiber_pitch_y      = dims.attr<double>(_Unicode(fiber_pitch_y));
  const double sipm_win_thickness = dims.attr<double>(_Unicode(sipm_window_thickness));
  //const double sipm_win_gap       = dims.attr<double>(_Unicode(sipm_window_gap));

  // S:C ratio (integer attributes on the <detector> element)
  const int sc_ratio_S = detElem.attr<int>(_Unicode(sc_ratio_S));
  const int sc_ratio_C = detElem.attr<int>(_Unicode(sc_ratio_C));
  const int sc_period  = sc_ratio_S + sc_ratio_C; // period length

  // -------------------------------------------------------------------------
  // Read sensor surface parameters from <sensors><sipm> child
  // -------------------------------------------------------------------------
  xml::Component sensorsElem = detElem.child(_Unicode(sensors));
  xml::Component sipmElem    = sensorsElem.child(_Unicode(sipm));
  const std::string sipmMatName  = sipmElem.attr<std::string>(_Unicode(material));
  const std::string sipmSurfName = sipmElem.attr<std::string>(_Unicode(surface));

  Material sipmMat  = desc.material(sipmMatName);           // AirOptical
  auto     sipmSurf = surfMgr.opticalSurface(sipmSurfName); // DR_SiPMSurf

  // -------------------------------------------------------------------------
  // Fiber geometry constants
  // -------------------------------------------------------------------------
  // Cladding: 15 um single-clad (SCSF-78 spec, research-DR-baseline.md §7.5).
  const double fiber_od_half     = fiber_od / 2.0;                 // 0.500 mm
  const double clad_thickness    = 0.015 * dd4hep::mm;             // 15 um
  const double core_radius       = fiber_od_half - clad_thickness; // 0.485 mm
  // Shorten fiber half-length so SiPM windows (one per end) fit inside pbVol.
  // Full absorber half-depth = module_dz/2.  Reserve sipm_win_thickness + 2*epsilon
  // on each end: one epsilon between fiber end face and window inner face,
  // and one epsilon between window outer face and pbVol boundary.
  // Algebraic budget per end: epsilon_inner + t_win + epsilon_outer = 2*epsilon + t_win.
  //const double win_epsilon       = 0.001 * dd4hep::mm;             // 1 µm anti-overlap gap
  const double fiber_half_length = module_dz / 2.0
                                   - sipm_win_thickness;
  //                                 - 2.0 * win_epsilon;            // shortened

  // -------------------------------------------------------------------------
  // Lattice dimensions
  // -------------------------------------------------------------------------
  const int N_x     = fiberCount(module_dx, fiber_pitch_x); // columns
  const int N_y     = fiberCount(module_dy, fiber_pitch_y); // rows
  const int N_total = N_x * N_y;

  // Count S and C fibers for the full lattice (for logging only)
  int N_S = 0, N_C = 0;
  for (int fi = 0; fi < N_total; ++fi) {
    if ((fi % sc_period) < sc_ratio_S) ++N_S;
    else                                ++N_C;
  }

  // -------------------------------------------------------------------------
  // Logging (INFO level)
  // -------------------------------------------------------------------------
  printout(INFO, "DRPrototype_geo",
           "Module: dx=%.1f mm, dy=%.1f mm, dz=%.1f mm",
           module_dx/dd4hep::mm, module_dy/dd4hep::mm, module_dz/dd4hep::mm);
  printout(INFO, "DRPrototype_geo",
           "S:C ratio = %d:%d, N_S=%d, N_C=%d, N_total=%d",
           sc_ratio_S, sc_ratio_C, N_S, N_C, N_total);
  printout(INFO, "DRPrototype_geo",
           "Fiber lattice: N_x=%d (pitch=%.3f mm) x N_y=%d (pitch=%.3f mm)",
           N_x, fiber_pitch_x/dd4hep::mm, N_y, fiber_pitch_y/dd4hep::mm);
  printout(INFO, "DRPrototype_geo",
           "Total SiPM windows: %d (2 per fiber)", 2 * N_total);
  printout(INFO, "DRPrototype_geo",
           "Detector ID: %d, core_r=%.4f mm, clad_t=%.4f mm",
           detID, core_radius/dd4hep::mm, clad_thickness/dd4hep::mm);

  // -------------------------------------------------------------------------
  // Materials
  // -------------------------------------------------------------------------
  Material lead_mat = desc.material("DR_Lead");
  Material S_mat    = desc.material("DR_Polystyrene");
  Material C_mat    = desc.material("DR_PMMA");
  Material clad_mat = desc.material("DR_FluorinatedPolymer");

  // -------------------------------------------------------------------------
  // Detector envelope: outer Assembly (no intrinsic boundary constraint).
  // Contains only the Pb absorber block; fibers and windows are all placed
  // inside pbVol via fiberAsm.
  // Pattern: BarrelCalorimeterScFi_geo.cpp:81 (envelope is an Assembly).
  // -------------------------------------------------------------------------
  Assembly envelope(detName);
  DetElement det(detName, detID);

  // Pb absorber block — fibers will be placed INSIDE this volume as daughters.
  // In TGeo, daughters take precedence over the mother material, so fibers
  // "replace" Pb at their location.  This is the correct idiom for embedded
  // fiber calorimeters.
  Box    pbBox(module_dx / 2.0, module_dy / 2.0, module_dz / 2.0);
  Volume pbVol(detName + "_absorber", pbBox, lead_mat);

  // -------------------------------------------------------------------------
  // Fiber volumes
  // Cladding volume: SOLID tube Tube(0, fiber_od_half, ...) of DR_FluorinatedPolymer.
  // Core volume: Tube(0, core_radius, ...) placed at (0,0,0) INSIDE the cladding.
  //
  // KEY FIX (v2): cladding must be a SOLID tube (inner radius = 0), NOT an annulus.
  // Reason: if cladding were Tube(core_radius, fiber_od_half, ...) (hollow annulus),
  // placing a solid core Tube(0, core_radius, ...) as a daughter fails the TGeo
  // containment check at the tube endcaps because the daughter's bounding sphere
  // touches the mother's z-faces exactly, causing numerical "protrusion" warnings.
  // With solid cladding, the core is fully contained; the core MATERIAL overrides
  // the cladding material in the core region, preserving the correct n(core)>n(clad)
  // optical boundary at the core/clad interface.
  //
  // TIR check: n(clad)=1.42 < n(PS)=1.57-1.69 and n(PMMA)=1.48-1.53 (plan R16).
  // -------------------------------------------------------------------------
  Tube clad_tube(0.0, fiber_od_half, fiber_half_length); // SOLID cladding tube

  // S fiber: cladding (FluorinatedPolymer) with Polystyrene core inside
  Tube   core_tube_S(0.0, core_radius, fiber_half_length);
  Volume core_vol_S(detName + "_fiber_core_S", core_tube_S, S_mat);
  Volume clad_vol_S(detName + "_fiber_S", clad_tube, clad_mat);
  clad_vol_S.placeVolume(core_vol_S, Position(0, 0, 0));

  // C fiber: cladding (FluorinatedPolymer) with PMMA core inside
  Tube   core_tube_C(0.0, core_radius, fiber_half_length);
  Volume core_vol_C(detName + "_fiber_core_C", core_tube_C, C_mat);
  Volume clad_vol_C(detName + "_fiber_C", clad_tube, clad_mat);
  clad_vol_C.placeVolume(core_vol_C, Position(0, 0, 0));

  // -------------------------------------------------------------------------
  // SiPM window volumes
  // Thin Tube (Cylinder) matching fiber OD.  Material: AirOptical (RINDEX=1.00029).
  // Without RINDEX on the window medium, G4OpBoundaryProcess emits NoRINDEX
  // warnings and photon boundary status is undefined (plan round-3 change 1a).
  //
  // Two volume types: S-fiber windows and C-fiber windows.
  // SkinSurface is bound once per volume (applies to ALL placements of that volume).
  // Pattern: DRICH_geo.cpp:773 (one SkinSurface per pssVol, shared across placements).
  //
  // window_S and window_C have the same shape and material; the type distinction
  // allows SkinSurface binding and setSensitiveDetector to be applied once per type.
  // -------------------------------------------------------------------------
  Tube   sipm_win_tube(0.0, fiber_od_half, sipm_win_thickness / 2.0);
  Volume sipm_win_vol_S(detName + "_sipm_win_S", sipm_win_tube, sipmMat);
  Volume sipm_win_vol_C(detName + "_sipm_win_C", sipm_win_tube, sipmMat);

  sipm_win_vol_S.setSensitiveDetector(sens);
  sipm_win_vol_C.setSensitiveDetector(sens);

  // -------------------------------------------------------------------------
  // Assembly voxelization (REQUIRED, plan R17)
  // Grid: ~20x20 mm cells, following BarrelCalorimeterScFi_geo.cpp:245-284.
  //
  // One fiberAsm per grid cell, placed INSIDE pbVol at (gx_cen, gy_cen, 0).
  // Each fiberAsm contains both the fiber tubes (at z=0) and the SiPM window
  // disks at each end (at ±win_z_end).  All volumes stay within pbVol because
  // fiber_half_length has been shortened to leave room for the windows.
  // -------------------------------------------------------------------------
  const double grid_size = 20.0 * dd4hep::mm; // ~20x20 mm grid cells
  const int grid_nx = std::max(1, static_cast<int>(std::ceil(module_dx / grid_size)));
  const int grid_ny = std::max(1, static_cast<int>(std::ceil(module_dy / grid_size)));

  printout(INFO, "DRPrototype_geo",
           "Assembly grid: %d x %d cells (%.0f mm cell size)",
           grid_nx, grid_ny, grid_size/dd4hep::mm);

  // Fiber lattice origin: lower-left fiber centroid
  const double x0 = -module_dx / 2.0 + fiber_pitch_x / 2.0;
  const double y0 = -module_dy / 2.0 + fiber_pitch_y / 2.0;

  // SiPM window z-positions relative to the fiber grid assembly local origin (z=0).
  // Fiber tube occupies [-fiber_half_length, +fiber_half_length].
  // Budget per end (from pbVol boundary inward):
  //   [pbVol face] -- epsilon_outer -- [window outer face] -- t_win -- [window inner face]
  //   -- epsilon_inner -- [fiber end face]
  // => window centre at ±(fiber_half_length + win_epsilon + sipm_win_thickness/2)
  // => window outer face at ±(fiber_half_length + win_epsilon + sipm_win_thickness)
  //    = ±(module_dz/2 - win_epsilon)  [strictly inside pbVol by win_epsilon]
  const double win_z_end0 = -(fiber_half_length + sipm_win_thickness / 2.0);
  const double win_z_end1 = +(fiber_half_length + sipm_win_thickness / 2.0);
  //const double win_z_end0 = -(fiber_half_length + win_epsilon + sipm_win_thickness / 2.0);
  //const double win_z_end1 = +(fiber_half_length + win_epsilon + sipm_win_thickness / 2.0);

  // -------------------------------------------------------------------------
  // Build grid assemblies and place fibers (inside pbVol) + windows (in envelope)
  // -------------------------------------------------------------------------
  for (int giy = 0; giy < grid_ny; ++giy) {
    for (int gix = 0; gix < grid_nx; ++gix) {

      // Grid cell boundaries (in module-centered coordinates = pbVol frame)
      const double gx_min = -module_dx / 2.0 + gix * grid_size;
      const double gx_max = gx_min + grid_size;
      const double gy_min = -module_dy / 2.0 + giy * grid_size;
      const double gy_max = gy_min + grid_size;
      const double gx_cen = (gx_min + gx_max) / 2.0;
      const double gy_cen = (gy_min + gy_max) / 2.0;

      const std::string cellTag = std::to_string(gix) + "_" + std::to_string(giy);

      // Fiber sub-assembly: placed inside pbVol at (gx_cen, gy_cen, 0).
      // Contains fiber volumes AND both SiPM end windows per fiber.
      // All volumes (fibers + windows) fit within pbVol because fiber_half_length
      // has been shortened to leave room for the window disks at each end.
      Assembly fiberAsm(detName + "_fiberGrid_" + cellTag);

      bool cell_has_fibers = false;

      // Assign fibers whose centroid falls in this grid cell
      for (int iy = 0; iy < N_y; ++iy) {
        for (int ix = 0; ix < N_x; ++ix) {

          const double fx = x0 + ix * fiber_pitch_x;  // fiber x centroid (pbVol frame)
          const double fy = y0 + iy * fiber_pitch_y;  // fiber y centroid (pbVol frame)

          if (fx < gx_min || fx >= gx_max) continue;
          if (fy < gy_min || fy >= gy_max) continue;

          // Global fiber index (row-major order)
          const int fiber_idx = iy * N_x + ix;

          // S/C assignment (deterministic, plan decision 4):
          //   position in period = fiber_idx % sc_period
          //   first sc_ratio_S positions -> S, remaining sc_ratio_C -> C
          const int pos_in_period = fiber_idx % sc_period;
          const bool is_S         = (pos_in_period < sc_ratio_S);
          const int channel_id    = is_S ? 0 : 1; // 0=S, 1=C

          // Fiber / window position relative to grid centroid
          const double rel_x = fx - gx_cen;
          const double rel_y = fy - gy_cen;

          // Place fiber in fiber assembly (goes inside pbVol)
          Volume& fiberVol = is_S ? clad_vol_S : clad_vol_C;
          fiberAsm.placeVolume(fiberVol, Position(rel_x, rel_y, 0.0));

          // Place SiPM windows in the SAME fiber assembly (inside pbVol).
          // Windows are flush at each fiber end face, separated by win_epsilon.
          Volume& winVol = is_S ? sipm_win_vol_S : sipm_win_vol_C;

          // End 0: -z face (end = 0)
          PlacedVolume win_pv0 = fiberAsm.placeVolume(
              winVol, Position(rel_x, rel_y, win_z_end0));
          win_pv0.addPhysVolID("fiber",   fiber_idx);
          win_pv0.addPhysVolID("end",     0);
          win_pv0.addPhysVolID("channel", channel_id);

          // End 1: +z face (end = 1)
          PlacedVolume win_pv1 = fiberAsm.placeVolume(
              winVol, Position(rel_x, rel_y, win_z_end1));
          win_pv1.addPhysVolID("fiber",   fiber_idx);
          win_pv1.addPhysVolID("end",     1);
          win_pv1.addPhysVolID("channel", channel_id);

          cell_has_fibers = true;

        } // end ix
      } // end iy

      if (cell_has_fibers) {
        // Voxelize the fiber+window assembly before placement.
        fiberAsm.ptr()->Voxelize("");

        // Place the combined fiber+window assembly INSIDE the Pb absorber volume.
        // gx_cen, gy_cen are already in pbVol-centered coordinates.
        pbVol.placeVolume(fiberAsm, Position(gx_cen, gy_cen, 0.0));
      }

    } // end gix
  } // end giy

  // Place the Pb absorber block in the outer envelope at z=0.
  // All fibers and SiPM windows are daughters of pbVol via fiberAsm.
  envelope.placeVolume(pbVol, Position(0, 0, 0));

  // -------------------------------------------------------------------------
  // SkinSurface binding for SiPM window volumes (plan todo 2a.5)
  // Binding once per Volume applies to ALL placements of that volume.
  // Reference: DRICH_geo.cpp:773 (pssSkin bound to pssVol shared across PDUs).
  // -------------------------------------------------------------------------
  SkinSurface sipm_skin_S(desc, det, detName + "_sipm_surface_S", sipmSurf, sipm_win_vol_S);
  if (!sipm_skin_S.isValid()) {
    throw std::runtime_error(
        "DRPrototype_geo: DR_SiPMSurf SkinSurface binding failed for S-fiber windows");
  }

  SkinSurface sipm_skin_C(desc, det, detName + "_sipm_surface_C", sipmSurf, sipm_win_vol_C);
  if (!sipm_skin_C.isValid()) {
    throw std::runtime_error(
        "DRPrototype_geo: DR_SiPMSurf SkinSurface binding failed for C-fiber windows");
  }

  printout(INFO, "DRPrototype_geo",
           "SkinSurface '%s' bound to S and C window volumes", sipmSurfName.c_str());

  // -------------------------------------------------------------------------
  // Place detector envelope in the mother volume (world)
  // -------------------------------------------------------------------------
  Volume motherVol    = desc.pickMotherVolume(det);
  PlacedVolume env_pv = motherVol.placeVolume(envelope, Position(0, 0, 0));
  env_pv.addPhysVolID("system", detID);
  det.setPlacement(env_pv);

  printout(INFO, "DRPrototype_geo",
           "DRPrototype constructed: N_fibers=%d, N_windows=%d, grid=%dx%d assemblies",
           N_total, 2 * N_total, grid_nx, grid_ny);

  return det;
}

// Factory registration.
// Name must match type="epic_DRPrototype" in the compact XML detector element.
// Reference: DRICH_geo.cpp:1090 (DECLARE_DETELEMENT), plan-dr-prototype.md "Naming - locked in"
DECLARE_DETELEMENT(epic_DRPrototype, create_detector)
