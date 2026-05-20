# Cameca Raw Data Import Tools

Import raw (pre-reconstruction) detector data from Cameca LEAP atom probes
into MATLAB. Supports three proprietary formats: RHIT, STR, and HITS.

## Supported Formats

| Format | Versions | Loader | Requirements | Status |
|--------|----------|--------|-------------|--------|
| **RHIT** | IVAS ≤ 3.6 | `rhitLoad` | Python 3 + uproot | Full decode |
| **STR** | IVAS ≤ 3.6 and v3-header STR companions | `strLoad` + `strCalculatePositions` | Pure MATLAB | Full decode |
| **HITS** | IVAS 3.8+ / AP Suite | `hitsLoad` (or `strLoad`) | Pure MATLAB | Header + event index + raw 8-byte ion payloads; TOF/position bit-packing **not decoded** |

## Quick Start

### RHIT files (CERN ROOT based)

```matlab
addpath('RHIT_imports');
[hits, histograms, metadata] = rhitLoad('mydata.RHIT');

% Voltage-corrected mass spectrum
massSpecPlot(hits, 0.01, 'normalised');

% Cameca's fully calibrated mass spectrum (stored in file)
ms = histograms.massSpectrum;
centers = (ms.edges(1:end-1) + ms.edges(2:end)) / 2;
figure; semilogy(centers, ms.values);
```

### STR files (delay-line timing data)

```matlab
addpath('RHIT_imports');
[hits, meta] = strLoad('mydata.STR');
hits = strCalculatePositions(hits);

% TOF spectrum (all events with at least 1 complete delay line)
good = hits.hitType >= 1;
figure; histogram(hits.tof(good), 1000);
xlabel('TOF (TDC counts)'); set(gca, 'YScale', 'log');

% Detector hit map (events with position)
hasPos = hits.hitType >= 2;
figure; histogram2(hits.detxRaw(hasPos), hits.detyRaw(hasPos), 200);
```

### HITS files (partial decode)

```matlab
addpath('RHIT_imports');
[hits, header, events] = hitsLoad('mydata.HITS');
fprintf('%d ion payloads, header: v%d, %s\n', ...
    height(hits), header.version, header.detectorLabels);

% hits has columns: ionIdx, eventIdx, byteOffset, markerOffset,
% markerVal16, parentType, hitInEvent, nHitsInEvent, b0..b7,
% eventSize, hasEventTrailer, isDelta, isMulti. The b0..b7 bytes are the raw
% hit-data block; mapping into TOF / detx / dety has not been
% reverse-engineered.
%
% Default payloadMode='recovered': SINGLE/STATUS/OTHER contribute one
% leading 8-byte ion block when present; MULTI events are split into
% complete leading 8-byte pairs. For forensic reverse engineering only,
% use hitsLoad(file,'payloadMode','allBlocks') to extract every complete
% post-marker 8-byte pair, including likely TLV trailers.

% Inter-event timing (= per-pulse interval, ~ns at ~117 kHz)
figure; histogram(double(hits.markerVal16), 200);
xlabel('inter-event time'); ylabel('count');

% events struct contains the full event index for STATUS / SINGLE / MULTI
% / OTHER events too (use it to walk status records or TLV trailers).

% If a paired EPOS exists, first check whether it is usable as ground truth.
report = hitsAuditEposAlignment('mydata.HITS', 'mydata.EPOS');
```

## Files

| File | Description |
|------|-------------|
| `rhitLoad.m` | Load RHIT files (requires Python 3 + uproot/h5py/numpy) |
| `rhitExtract.py` | Python helper — extracts RHIT to HDF5 |
| `rhitCalibrateFromStoredMass.m` | Voltage + 2D bowl calibration **from the embedded `pMass` only** — no EPOS needed; works on any RHIT |
| `rhitCalibrateFromEpos.m` | Voltage + 2D bowl calibration from a paired RHIT/EPOS file (advanced; higher accuracy when EPOS is available) |
| `rhitApplyCalibration.m` | Apply a saved calibration (new 2D or legacy radial) to an RHIT table |
| `rhitMatchEvents.m` | Drift-tracking RHIT↔EPOS event matcher (two-pass, TOF-anchored) |
| `rhitFitBowl2D.m` | Fit voltage offset + 2D xy polynomial bowl from matched events |
| `rhitTOffsetXCorr.m` | Robust TOF-offset estimator via cross-correlation of TOF histograms |
| `rhitRefineTOffset.m` | Refine `t_offset` against a reference spectrum (no EPOS needed) |
| `rhitPoolBowl.m` | Pool per-run 2D bowls into one global instrument bowl |
| `rhitEvaluateBowl.m` / `rhitDesignXY.m` | Polynomial bowl primitives |
| `Workflow_RHIT_Bowl_Calibration.m` | End-to-end example: per-run cal → pooled bowl → apply to new RHIT |
| `strLoad.m` | Load STR (v2) files; for HITS (v3) delegates to `hitsLoad` |
| `strCalculatePositions.m` | Compute detector positions, TOF, recover partial STR hits |
| `hitsLoad.m` | Load HITS (v3) files — partial decode (raw bytes per ion) |
| `hitsScanEvents.m` | Walk a HITS file: parse header + index every event by type |
| `hitsExtractPayloads.m` | Per-event payload byte ranges (offset+length) |
| `hitsAuditEposAlignment.m` | Audit HITS payload-stream alignment against paired EPOS raw columns |
| `hitsAuditStrAlignment.m` | Audit same-run HITS payloads against STR event/timing streams |
| `hitsStrEventMap.m` | Shared HITS/STR event-axis map; identifies missing payload classes |
| `hitsAnalyzeStrOnlyEvents.m` | Summarize STR-channel events without complete HITS payload blocks |
| `hitsSearchStrBitfields.m` | Exploratory HITS bitfield search against decoded STR timing targets |
| `hitsSearchStrDeltas.m` | Test HITS byte/word predictors against same-run STR previous-ion deltas |
| `hitsSearchStrMultiPairs.m` | Test HITS MULTI hit2-hit1 byte differences against nearby STR hit-pair differences |
| `hitsVerifyStrMultiCandidate.m` | Robust/shuffled verification for the leading multi-pair candidate |
| `HITS_format_analysis.md` | Reverse-engineering notes for the HITS v3 format |
| `README.md` | This file |

---

## RHIT Format

### Requirements

```
pip3 install uproot h5py numpy
```

### Output columns (rhitLoad)

| Column | Unit | Description |
|--------|------|-------------|
| `ionIdx` | — | Sequential event index (1-based) |
| `detx` | mm | Detector X position (converted from LSB) |
| `dety` | mm | Detector Y position (converted from LSB) |
| `mc` | Da | Voltage-corrected mass-to-charge ratio |
| `tof` | ns | Raw time-of-flight |
| `VDC` | V | Specimen DC voltage |
| `Vref` | V | Reference voltage (= VDC × kf) |
| `detxRaw` | LSB | Raw detector X (digitiser units) |
| `detyRaw` | LSB | Raw detector Y (digitiser units) |

Plus per-event instrument state: `pulse`, `freq`, `chi2`, `erate`,
`tElapsed`, `Temp`, `Pres`, `VMcpGain`, `VAnodeAccel`, etc.

### Mass-to-charge calculation

The provisional formula `rhitLoad` populates into `hits.mc`:
```
mc = (2·e / u·L²) · VDC · (tof − t0)² / sqrt(kf)
```
- `L` (flight path) and `t0` from `CRunHeader`
- `kf` = voltage correction factor (Vref/VDC, typically 1.03)
- Dominant peak accuracy ~0.3%; secondary peaks 2–5% (no bowl correction)
- `histograms.massSpectrum` has Cameca's spectrum *as it was when the
  RHIT was last saved* (sometimes uncorrected, sometimes bowl-corrected;
  do not rely on it as a reference).

For a full calibration of *any* RHIT file (no EPOS required), use
`rhitCalibrateFromStoredMass`:

```matlab
[hits, hist, meta] = rhitLoad('myRun.RHIT');
[hits, calib] = rhitCalibrateFromStoredMass(hits, hist, meta);
% hits.mc is now calibrated.  calib has fields tOffsetNs, bowl, ...
massSpecPlot(hits, 0.01, 'normalised');
```

This fits `mc = C(x, y) · VDC · (tof − t_offset)²` with a 2D xy
polynomial bowl `C(x, y)` and a single TOF offset by aligning local
detector-cell spectra to the embedded `pMass` histogram. Works on any
LEAP RHIT regardless of sample composition.

If a paired EPOS file *is* available, `rhitCalibrateFromEpos` gives a
slightly more accurate fit because it can match individual ions
across the two streams instead of going through histogrammed spectra.
Same calibration struct format; same `rhitApplyCalibration` consumer.

Equivalent to Cameca's standard workflow `Mass = k · (T − dT₀)²` with
voltage-corrected TOF and 2D bowl, just re-parameterised: our
`t_offset` ↔ Cameca's `dT₀`; our `C(x, y) · VDC` ↔ Cameca's
`k / (a₁ + V_DC + a₂·V_DC²)`. We model `a₂ = 0` (sufficient for
steady-state runs).

See `Workflow_RHIT_Calibration_LiveScript.m` for the full RHIT-only
walkthrough as a MATLAB live script (open in MATLAB and right-click →
*Run*; renders as a live-editor document with embedded plots).

### Instrument parameters (metadata.instrumentParams)

| Parameter | Example | Description |
|-----------|---------|-------------|
| `flight_path_mm` | 382.0 | Reflectron flight path |
| `t0_ns` | 45.0 | TOF zero offset |
| `detector_halfsize_mm` | 18.5 | Detector active half-width |
| `lsb_to_mm` | 0.0247 | LSB-to-mm conversion |
| `mcp_gain_voltage_V` | 1990 | MCP gain voltage |
| `ivas_version` | 215.41.342j | IVAS software version |
| `bowl_correction_coefficients` | [20×1] | 4th-order polynomial |

### RHIT file structure

CERN ROOT file (magic: `root`) containing a TTree with 28 per-event branches,
standard ROOT histograms (mass spectrum, TOF, voltage/erate history, detector
map), and custom Cameca classes (CVoltage, CBowl, CCalibMass, CRunHeader,
CHits, CDetectorCal, CEfficiency2D, etc.). Standard ROOT objects are fully
readable via `uproot`; custom classes are partially decoded.

---

## STR / HITS Format

### Binary encoding

Fixed 4-byte TLV records throughout: `[int24_LE value][uint8 tag]`

File header (first ~128 bytes) followed by interleaved pulse counters and
event data. Tag 0x18 delimits events; segment markers (0x18 with value=5000
between pulse counters) are filtered out automatically.

Some AP Suite-era `.STR` files carry a v3 header but still use STR-style
TLV channel records (`0x01..0x04`, `0x21`, `0x22`) in bulk data. `strLoad`
now uses the file extension to distinguish `.HITS` from `.STR`: v3 `.HITS`
delegates to `hitsLoad`, while v3 `.STR` continues through the TLV channel
parser.

### Detector geometry

Three delay lines at 0°, 90°, and 45°:

| Delay line | Angle | Tags (STR v2) | Header labels |
|------------|-------|---------------|---------------|
| **x** | 0° | 0x01, 0x02 | varies per instrument |
| **y** | 90° | 0x03, 0x04 | varies per instrument |
| **w** | 45° | 0x21, 0x22 | varies per instrument |

### strLoad output columns

| Column | Unit | Description |
|--------|------|-------------|
| `ionIdx` | — | Sequential event index (1-based) |
| `detxt1`, `detxt2` | TDC | DL-x end timings |
| `detyt1`, `detyt2` | TDC | DL-y end timings |
| `detwt1`, `detwt2` | TDC | DL-w end timings (45° diagonal) |
| `quality` | — | Hit finding quality / chi² |

Events with missing channels have NaN for those fields. Only events with at
least one channel value are imported (empty segment markers are discarded).

### strCalculatePositions output (added columns)

| Column | Unit | Description |
|--------|------|-------------|
| `detxRaw` | TDC | DL-x position (= detxt1 − detxt2) |
| `detyRaw` | TDC | DL-y position (= detyt1 − detyt2) |
| `detwRaw` | TDC | DL-w position (= detwt1 − detwt2) |
| `tof` | TDC | Time-of-flight (mean of offset-corrected DL sums) |
| `hitType` | — | 0–3: number of complete delay lines |

### Hit types

| hitType | Meaning | Has position? | Has TOF? |
|---------|---------|---------------|----------|
| 3 | Full hit (all 3 DLs) | Yes | Yes |
| 2 | Partial hit (2 DLs, 3rd recovered) | Yes | Yes |
| 1 | Single DL only | No | Yes |
| 0 | Incomplete (single ends only) | No | No |

### Partial hit recovery

`strCalculatePositions` fits the 45° geometric constraint from full hits:

```
detwRaw ≈ a · detxRaw + b · detyRaw + c
```

For events with only 2 complete delay lines, the missing position is recovered
using this relationship. The fit is calibrated per-file from sum-consistent
full hits (residual std typically 40–100 TDC counts).

### TOF calculation

Each complete delay line provides a TOF estimate: `tof_DL = (t1 + t2) / 2`.
Per-DL offsets (from different delay-line lengths) are removed using medians
from full hits. The final TOF is the mean of all available DL estimates.

### HITS v3 — current decode state

HITS shares the 4-byte TLV header with STR but uses an entirely different
bulk encoding (events delimited by tag `0xa4` rather than per-channel TLV).
What `hitsLoad` returns:

| Layer | Decoded? | Notes |
|---|---|---|
| File header | yes | version, channel labels, voltage, thresholds, walk corrections, 5 timing channels (with hi/lo extension), and ~16 unknown header fields exposed verbatim |
| Event index | yes | every `0xa4` marker classified as STATUS / SINGLE / MULTI / OTHER, with byte offset and record gap |
| Multi-hit grouping | yes | gap=3/5/7/9 → 1/2/3/4 ions per pulse |
| TLV trailers on SINGLE events | yes | standard TLV state-update records appended after some hits (pulse counter, voltage, timing channel update); tag values returned, semantic meaning per HITS_format_analysis.md |
| 8-byte hit-data block (d1, d2) | **partial** | default `payloadMode='recovered'` returns one `b0..b7` row from the leading complete block of SINGLE/STATUS/OTHER events and splits MULTI events by complete leading 8-byte pairs. Leftover records are flagged as trailers. Multi-hit byte-role assay shows bytes 0,1,4,5 are independent per ion (position-like) and bytes 2,6 are shared/delta-encoded within a multi (TOF-like). `d1.b3` bit 0 marks delta-encoded continuation hits. |
| Bit-level mapping → TOF / detx / dety | **no** | unknown bit-packing scheme; both direct decoding attempts and cumsum-residual / delta models failed. See `HITS_format_analysis.md` section 7. |

### HITS v3 — bitfield search ruled out (2026-05-02)

Direct bit-field extraction from the 8-byte ion payload was tested
exhaustively against R5124_00368, the same-run HITS+STR pair (18.3M
aligned ions, 0.003% event-count diff). All apparent correlations
collapse to **r ≈ 0** under detrending — they were artifacts of slow
experiment-time drift, not real bit-fields. See `HITS_format_analysis.md`
section 7 for the decisive experiment and detrending evidence.

Conclusion: the 8-byte payload is NOT a fixed bit-packing of detector
quantities. The encoding requires more involved reverse engineering
(predictive / residual coding, dictionary-based, or arithmetic coded).

### Open avenues for the bit-decode (ranked)

1. **Multi-hit pair-wise prediction test**. First-pass search found an
   apparent `h2.b2` ↔ next-STR `detyt2` signal, but robust verification
   showed it is driven by abnormal continuation-byte / extreme STR-delta
   outliers. With `h2.b2 <= 64` and central STR deltas, the raw correlation
   collapses to ~0.
2. **Delta correlation**. First-pass byte/word test against
   `STR_value(i) − STR_value(i−1)` did not fire: top detrended |r| was
   ~0.026 for SINGLE and ~0.042 for STATUS on R5124. A simple direct
   previous-ion residual model is therefore unlikely.
3. **Stratified bitfield search by `d1.b3` mode class** (~20% odds).
   Run `hitsSearchStrBitfields` separately per value of
   `bitand(b3, 0x3c)`.
4. **AP Suite SDK / DLL decompilation**. The CAMECA reader is a binary;
   static analysis of its bit-extraction routines reveals the algorithm.
5. **Reach out to FAIRmat / pynxtools-apm** community.
6. **CAMECA research collaboration**.

For users today: export to `.epos` or `.apt` from AP Suite for ranged
analysis. `hitsLoad` is useful for header inspection, ion accounting,
multi-hit grouping, file QC, and as the validated substrate for further
decode work.

---

## Limitations

- **RHIT mc**: no bowl correction applied; use `histograms.massSpectrum` for
  Cameca's calibrated spectrum
- **STR**: all values in raw TDC counts (no mm or ns conversion);
  requires instrument-specific calibration for physical units
- **HITS v3**: bit-level packing of the 8-byte hit-data block is unknown
  — `hitsLoad` returns the raw bytes only, not TOF / detx / dety. See
  `HITS_format_analysis.md` for the full decode state.
- **RHIT**: requires Python 3 accessible via `python3`
