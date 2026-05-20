# CLAUDE.md — Cameca Raw Data Importers

## Purpose

This folder contains reverse-engineered importers for Cameca's proprietary raw
detector data formats (RHIT, STR, HITS). These are NOT part of the public
Atom Probe Toolbox — the folder is gitignored.

## Format Knowledge

### RHIT (IVAS ≤ 3.6)

- **Container**: CERN ROOT file (magic bytes: `root`)
- **Readable via**: Python `uproot` library
- **TTree `nth`**: 28 branches of per-event data (x, y as int16 LSB; tof as
  float32 ns; v/Vref as float32 V; plus thresholds, walk, pulse, etc.)
- **Histograms**: pMass (calibrated mass spectrum), pTofRaw, phV, phE, pXyHist
- **Custom classes** (C-prefixed): wrap standard ROOT objects in proprietary
  serialization. CRunHeader is partially decoded (instrument params at known
  byte offsets). CVoltage=TProfile, CBowl=2D lookup, CCalibMass=TGraphErrors.
- **Key CRunHeader offsets** (after 6-byte header + 10-byte TObject):
  - 448: detector_halfsize_mm (float32 BE)
  - 452: t0 (float32 BE, in TDC counts not ns)
  - 456: flight_path_mm (float32 BE)
  - 496: kf voltage correction factor (float32 BE, typically 1.03)
  - 504: MCP gain voltage (float32 BE)
- **Detector coords**: int16 LSB, convert to mm via detector_halfsize/750
- **mc formula**: `mc = 2e/(u·L²) · VDC · (tof−t0)² / sqrt(kf) · 1e-18`

### STR (IVAS ≤ 3.6, companion to RHIT)

- **Container**: Custom binary, fixed 4-byte TLV: `[int24_LE value][uint8 tag]`
- **Version byte**: 0x02 (first byte of file)
- **File is 4-byte aligned** (filesize mod 4 = 0)
- **Header** (~128 bytes): tags 0xa0-0xa3 (metadata), 0x27-0x2c (detector
  channel labels as ASCII), 0x2d-0x2f (thresholds), 0x33-0x38 (walk), 0x1b
  (voltage), 0x49 (timestamp)
- **Bulk data**: interleaved pulse counters (0x05/0x0b) and hit events
- **Per-event channel tags**: 0x01 (detxt1), 0x02 (detxt2), 0x03 (detyt1),
  0x04 (detyt2), 0x21 (detwt1), 0x22 (detwt2)
- **Event delimiter**: tag 0x18 (quality/chi2). Also used as segment marker
  (value=5000) between pulse groups — these are filtered out.
- **Detector geometry**: 3 delay lines at 0° (x), 90° (y), 45° (w).
  Constraint: `detwRaw ≈ 0.65·detxRaw + 0.65·detyRaw + offset`
  (coefficients fitted per-file from sum-consistent full hits).

### HITS (IVAS 3.8+ / AP Suite)

- **Container**: Same 4-byte TLV header as STR
- **Version byte**: 0x03
- **File is NOT 4-byte aligned** (filesize mod 4 = 2; 2-byte prefix before TLV)
- **Header**: identical tag structure to STR
- **Bulk data**: DIFFERENT encoding from STR — NOT standard TLV. The encoding
  has not been reverse-engineered. Only the initial ~10k events in the header
  region are extractable.
- **Known but unsolved**: the bulk uses some form of delta/varint encoding;
  multiples-of-16 tags (0x10,0x20,...,0xf0) sum to ~7.4M matching EPOS count.

## Code Architecture

- `rhitLoad.m` → calls `rhitExtract.py` → temp HDF5 → MATLAB table
- `strLoad.m` → pure MATLAB binary reader → table with raw DLD timings
- `strCalculatePositions.m` → computes detxRaw/detyRaw/tof from DLD timings,
  recovers partial hits using 45° geometric constraint

### Voltage + 2D bowl calibration pipeline

Mirrors the original Python prototype in
`010113_data-pool/APT_data/calibration/scripts/`.

Driver: `rhitCalibrateFromEpos(hits, epos)` returns a calibration
struct with fields:
- `tOffsetNs`         — TOF zero offset (median of tof_rhit − tof_epos
                        across matched events; std ~0.01 ns when matching is solid)
- `bowl`              — 2D polynomial: `kind='xy2d'`, `degree=4`,
                        `coeffs(15)`, `terms(15×2)` of (i, j) exponents.
                        Evaluated by `rhitEvaluateBowl(bowl, x, y)`.
- `ICF`               — image-compression factor `epos.detx / hits.detx`
- `nMatched`, `fitResidualPctC0`, etc.

Pipeline stages (each is a separate function so they can be reused):
1. `rhitTOffsetXCorr(rhitTof, eposTof)` — robust TOF-offset seed via
   FFT cross-correlation of TOF histograms (no toolbox deps).
2. `rhitMatchEvents(hits, epos, ...)` — drift-tracking matcher.
   Pass 1 builds a piecewise-linear `r_offset(ei)` by anchoring at
   ~200 EPOS indices with a wide TOF+position search; Pass 2 does a
   tight tolerance match with the model predicting each `r_predict`.
3. `rhitFitBowl2D(hits, epos, matchedE, matchedR)` — least-squares fit
   of `mc_e = C(x,y) · VDC · tof_e²` for a 2D polynomial of x, y.
   Falls back to `kind='radial'` (polynomial in r²) on request.
4. `rhitApplyCalibration(hits, calib)` — populates `hits.mc`. Accepts
   both new (`calib.bowl` struct) and legacy (`calib.C_poly` radial)
   formats so existing scripts keep working.

Pooling across runs: `rhitPoolBowl({calib1, calib2, ...})` samples
each per-run bowl on a common (x, y) grid, takes the median per cell,
and refits a single polynomial. On a single LEAP, the 2D bowl is
stable to ~0.1% in `C0` across runs; the per-run free parameter is
just `t_offset`.

Cross-correlation `t_offset` refinement against a reference spectrum:
`rhitRefineTOffset(hits, refSpectrum, bowl, tOffsetGuess)` — used to
apply a pooled global bowl to RHIT files that lack a paired EPOS or
where the EPOS-anchored matcher fails (e.g., R56_01700 in the data
pool).

Walkthrough: `Workflow_RHIT_Bowl_Calibration.m`.

### Relation to Cameca's IVAS workflow

The Cameca *2015 LEAP Training Class* slides
(`30129G_4_Advanced Reconstruction and Analysis_JULY2015_Erlangen.pdf`)
formalise the standard interactive workflow:

- Cameca's mass formula: `Mass = k · (T − dT₀)²` with
  voltage-corrected TOF `T = a₀ / sqrt(a₁ + V_DC + a₂·V_DC²)`.
- Cameca's bowl: "radial + N angular sectors" (`Bowl Subdivide` UI,
  default `4` = radial + 96 sectors).
- Cameca iterates voltage→bowl→voltage→... until mass-resolving power
  stops improving, with a chosen reference peak.

Our `mc = C(x,y) · VDC · (tof − t_offset)²` is the same physics
re-parameterised: Cameca's `dT₀` → our `t_offset`; Cameca's
`k / (a₁ + V_DC + a₂·V_DC²)` → our `C(x, y) · VDC`. We model `a₂ = 0`
(linear in `V_DC` inside the sqrt), which is correct for steady-state
runs but not for runs with large compositional or tip-shape changes.
Our 2D xy polynomial (15 terms, degree 4) is a different basis than
Cameca's "radial + sectors" but captures the same asymmetric bowl.

Cross-check: Cameca's `K = 0.193 amu·mm²/(kV·ns²)` divided by `L_eff²`
gives our pooled `C0 = 1.6432e-9` exactly (`L_eff = 342.7 mm` vs nominal
flight path 382 mm — the standard reflectron folding).

Things the standard IVAS workflow does that this pipeline currently
skips (none have been needed for the 17 LEAP runs tested so far):
- Z-range selection from the V(z) curve (drop cap / oxidation / etc.).
- XY FOV ROI (drop low-efficiency detector edges; IVAS default 95% of
  ions inside fiducial ellipse).
- Voltage↔bowl iteration loop.
- The single-peak anchor variant (largest peak in the spectrum used as
  the V-correction reference, instead of EPOS or pMass).

## Testing

Reference data pairs for validation:
- R56_03063: has both .RHIT and .EPOS (matched events confirm tof offset,
  bowl correction, and detector coordinate scale)
- R5121_12823: has .HITS and .EPOS (EPOS only has mc, other fields zero)
- STR files in `exampleFiles_backup/external_repos/17136940/`: Fig6a_P3.STR,
  Fig6b_P7.STR, FigS2a.STR, FigS2b.STR (all decode successfully)

## Known Issues / Future Work

- RHIT bowl correction: superseded by `rhitFitBowl2D` which fits a 2D
  xy polynomial directly from matched EPOS data — reproduces Cameca's
  bowl to ~0.06% of `C0` (radial-only stalls at ~2.6%).  The
  `bowl_correction_coefficients` from CRunHeader are still surfaced
  but not used.
- RHIT CCalibMass: TGraphErrors identified but data arrays not extracted
  (complex ROOT reference/tag serialization).  Not needed now that the
  bowl can be fit directly.
- STR/HITS: all values in raw TDC counts; no instrument-specific calibration
  for mm/ns conversion (detector size, TDC clock period not in header)
- HITS v3 bulk decoding: main unsolved problem
- strCalculatePositions: hitType=0 events (single DL ends) could potentially be
  filtered at import stage
