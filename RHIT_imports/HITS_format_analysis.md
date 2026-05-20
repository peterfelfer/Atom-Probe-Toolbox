# HITS Format (v3) Reverse-Engineering Report

Status: **partially decoded** — file structure, event indexing, ion-payload extraction, and same-run STR alignment fully working. Bit-level decoding of the 8-byte ion payload into TOF / detx / dety has been definitively ruled out as a direct bit-extraction problem (see section 7); the encoding requires a more sophisticated decoder.

Analysis performed primarily on R5121_12823 (HITS + decimated EPOS) and R5124_00368 (HITS + same-run STR — the strongest test pair).

---

## 0. Quick state for whoever picks this up

### What works today (`hitsLoad.m`)
- Header parser (version, channel labels, voltage, thresholds, walk corrections, 5 timing channels with hi/lo extension, ~16 unknown header tags exposed verbatim).
- Vectorised event walker — 100% match against reference event counts on R5121.
- Per-ion payload extraction via `payloadMode = 'recovered'`: SINGLE + STATUS + OTHER first-block payloads + all complete pairs in MULTI events. Yields ~17.99M ions for R5124, ~7.54M for R5121.
- Multi-hit splitting per ion with `hitInEvent`, `nHitsInEvent`, `parentType`, `isMulti`, `isDelta`, `hasEventTrailer`.
- TLV trailer detection (size > 1 + 2·payloadCount).
- Same-run STR alignment infrastructure: `hitsScanEvents`, `strScanEvents`, `hitsStrEventMap`, `hitsAuditStrAlignment`, `hitsAnalyzeStrOnlyEvents`, `hitsSearchStrBitfields`, `hitsAuditEposAlignment`.
- Validated structural alignment to same-run STR for R5124: 18,307,963 HITS events vs 18,307,332 STR events (0.003% diff), 98.2% co-presence rate at offset 0.

### What does NOT work
- TOF in nanoseconds. detx / dety in mm or TDC counts. mc in Da. VDC per pulse. None of these can be recovered from the 8 raw payload bytes by direct bit extraction. Tested extensively — see section 7.

### Most useful next experiment to try
**Conditioned / stateful residual decoding.** A first-pass previous-ion
delta test does not support a direct `byte_i ≈ STR(i)-STR(i-1)` model. The
remaining plausible residual models are conditional on mode class
(`d1.b3`) or depend on hidden running state rather than the immediately
previous decoded STR ion.

### Other paths if predictive decoding fails
1. **AP Suite SDK / DLL inspection**. CAMECA provides a C# / C++ SDK that reads HITS. Decompile or trace it.
2. **FAIRmat / pynxtools-apm community**. Ask if HITS bit-packing has been documented elsewhere.
3. **CAMECA support directly**. The format is proprietary but a research collaboration may unlock it.

### Hard data we have
| File | Size | Pair | Notes |
|---|---|---|---|
| `data-pool/APT_data/R5124_00368/R5124_00368.HITS` | 226.6 MB | + same-run STR | Strongest test pair. PPfdUD detector labels. |
| `data-pool/APT_data/R5124_00368/R5124_00368.STR` | 522.8 MB | + HITS above | v3 header but STR-style TLV bulk; `strLoad` parses it via file-extension dispatch. |
| `data-pool/APT_data/R5121_12823.HITS` | 112.0 MB | + decimated EPOS | EPOS has only `x, y, z, mc, atomNum`; raw fields zeroed. RPilbR detector labels. |
| `data-pool/APT_data/R5121_12823.EPOS` | 326.6 MB | + HITS above | 7,423,312 reconstructed ions. |

---

---

## 1. File Structure

The HITS format is version 3 of the STR/HITS binary format used by CAMECA LEAP instruments.

### Basic encoding

The entire file uses **4-byte TLV records**: `[24-bit signed value LE][8-bit tag]`.

```
byte0 + byte1*256 + byte2*65536 = 24-bit value (sign-extended from bit 23)
byte3                            = tag
```

The file is NOT always a perfect multiple of 4 bytes (R5121_12823.HITS has 2 trailing bytes).

### Magic number

All HITS files start with `03 10 80 a0` (tag `0xa0`, value encodes version 3).

---

## 2. Header (fields 0 to first `0xa4` tag)

The header uses standard TLV tags. Confirmed across 4 different instruments (R5121, R5124, R6001, R6006):

| Tag    | Meaning                    | Example (R5121)      |
|--------|----------------------------|----------------------|
| `0xa0` | Version (v3 = HITS)        | 3                    |
| `0xa2` | Unknown (always 0)         | 0                    |
| `0xa3` | Unknown (varies)           | 6                    |
| `0xa1` | Unknown                    | -7389528             |
| `0x27`-`0x2c` | 6 detector channel labels (ASCII) | 'R','P','i','l','b','R' |
| `0x2d`-`0x2f` | 3 detection thresholds   | 2500, 2500, 2500     |
| `0x33`-`0x38` | 6 walk corrections       | 980 each             |
| `0x3b` | Timing channel 1 (lo 24-bit) | 1029              |
| `0x4c` | Timing channel 1 (hi 24-bit) | 0                 |
| `0x3c`/`0x4d` | Channel 2 lo/hi        | 1281 / 0             |
| `0x3d`/`0x4e` | Channel 3 lo/hi        | -1000 / 0            |
| `0x3e`/`0x4f` | Channel 4 lo/hi        | 85 / 0               |
| `0x3f`/`0x50` | Channel 5 lo/hi        | 653 / 0              |
| `0x09` | Unknown (always 0)         | 0                    |
| `0x1b` | Initial voltage (V)        | 500                  |
| `0x49` | Unknown counter            | 5925282              |
| `0x40` | Unknown                    | 0                    |
| `0x13` | Unknown (instrument param) | -641079              |
| `0x14` | Unknown                    | -522951              |
| `0x15` | Unknown                    | -2644710             |
| `0x20` | Voltage/pulse info         | 1441805              |
| `0x47` | Unknown param              | 5998                 |
| `0x4b` | Unknown param              | 2004                 |
| `0x1e` | Timing reference           | 499743               |
| `0x23`-`0x25` | Timing constants (=5000000 in all files) | 5000000 |
| `0xa8` | Multi-hit indicator        | 3                    |

The 5 timing channels (tags `0x3b`-`0x3f`) with high-extension tags (`0x4c`-`0x50`) provide initial/reference TDC values. Each full channel value = `hi * 2^24 + lo` (48 bits total).

### Channel labels across instruments

| File        | Labels (tags 0x27-0x2c) | Voltage | Param 0x47 |
|-------------|-------------------------|---------|------------|
| R5121_12823 | R P i l b R             | 500 V   | 5998       |
| R5124_00368 | P P f d U D             | 1000 V  | 4999       |
| R6001_70985 | P K x x ? ?             | 2000 V  | 5000       |
| R6006_34092 | z g < < D <             | 500 V   | 4001       |

---

## 3. Event Delimiting

After the header, events are delimited by **tag `0xa4`**. The byte2 of the a4 marker is a *flavor* flag, not a hit/no-hit discriminator:

| byte2  | Flavor      | Count (R5121) | Carries ion payload? |
|--------|-------------|---------------|----------------------|
| `0x00` | STATUS      | 1,225,266     | **Yes (~90%)** — first 8 bytes after marker |
| `0x20` | SINGLE      | 6,343,221     | Yes (~98%) |
| `0x40` | MULTI       | 352,538       | Yes (~46%) |
| other  | OTHER       | 6,200         | Yes (~98%) |
| **Total** |          | **7,927,225** |                                      |

### STATUS payloads are real ions, not bookkeeping

The original reading "STATUS = pulse-counter event with no ion" was wrong.
A byte-fingerprint comparison of payloads pulled from SINGLE-flavor vs
STATUS-flavor events shows they are statistically indistinguishable:

| Byte | SINGLE mean / entropy | STATUS mean / entropy |
|---|---|---|
| b0 | 127.5 / 8.00 | 127.0 / 8.00 |
| b1 | 128.2 / 8.00 | 128.1 / 8.00 |
| b2 | 42.4 / 5.58 | 32.2 / 4.71 |
| b3 | 4.8 / 2.17 | 6.2 / 2.37 |
| b4 | 126.8 / 8.00 | 126.5 / 8.00 |
| b5 | 127.6 / 5.24 | 128.0 / 5.26 |
| b6 | 127.1 / 7.18 | 126.1 / 7.17 |
| b7 | 60.5 / 6.00 | 60.3 / 6.00 |

The d1.b3 distribution clinches it: **98.73% of STATUS payloads carry an
ion-like b3** (multiples of 4 — `0x00, 0x04, 0x08, 0x0c, 0x10, ...`),
versus 99.82% for SINGLE. The two are the same kind of event with a
different flavor flag, and both embed an ordinary 8-byte ion payload in
the first 2 records after the marker.

### Ion count / alignment status

`hitsLoad` with `payloadMode='recovered'` (default) emits one ion-payload
row from the complete leading block of SINGLE/STATUS/OTHER events and all
complete leading 8-byte pairs from MULTI events. Counts on R5121:

| Source | Ion payloads |
|---|---|
| SINGLE flavor | 6,237,078 |
| STATUS flavor | 1,098,901 |
| MULTI flavor (per ion) | 194,524 |
| OTHER flavor | 6,098 |
| **Total** | **7,536,601** |
| EPOS rows | 7,423,312 |
| Diff (HITS − EPOS) | +113,289 (1.5%) |

The 1.5% overshoot likely reflects EPOS dropping low-quality ions during
reconstruction. Marker-only SINGLE (103,629), empty MULTI (171,179), and
incomplete MULTI (2,291) events are excluded — they carry no payload.

Same-run R5124 HITS+STR validation after the recovered-mode correction:

| Quantity | Count |
|---|---:|
| HITS events | 18,307,963 |
| STR events | 18,307,332 |
| STR channel-bearing events | 18,302,600 |
| STR complete six-channel events | 17,884,884 |
| HITS recovered ion payload rows | 17,987,943 |
| STATUS first-block payload rows | 377,094 |
| SINGLE payload rows | 17,398,749 |
| MULTI payload rows | 210,487 |
| OTHER first-block payload rows | 1,613 |
| STR channel events - HITS payload rows | +314,657 |

Forensic `payloadMode='allBlocks'` on R5124 returns 19,084,248 rows, but
this over-extracts complete TLV trailer pairs as if they were ions. It is
useful for byte-role analysis, not as a default ion table.

The remaining unrecovered R5124 STR-channel events are 322,455 event-axis
positions, dominated by marker-only classes with no complete post-marker
8-byte block: ~204k `MULTI_nf1_dl3` and ~104k `SINGLE_nf1_dl3`. This is now
a marker/event-axis semantics problem, not a hidden post-marker payload
problem.

---

## 4. STATUS Events (byte2 = 0x00)

STATUS events use traditional TLV encoding within the a4-delimited segment. Confirmed tags in bulk STATUS events:

| Tag    | Meaning          | Behavior                              |
|--------|------------------|---------------------------------------|
| `0x05` | Pulse counter A  | Monotonically increasing, paired with 0x0b |
| `0x0b` | Pulse counter B  | Same value as 0x05 partner            |
| `0x49` | Event counter    | Appears every ~500 pulse pairs        |
| `0x1e` | Timing reference | ~499,700 (nearly constant)            |
| `0x20` | Voltage/pulse    | Appears rarely, slowly changing       |

Pulse counters increment by ~20-28 per pair, with ~5 pairs per STATUS event. Tags `0x3b`-`0x50` (timing channels) appear ONLY in the header region, NOT in bulk STATUS events.

---

## 5. SINGLE HIT Events (byte2 = 0x20)

Most common event type. Variable size — the structure is:

```
[a4 marker: byte0, byte1, 0x20, 0xa4]   -- 4 bytes
[hit data:  d1, d2]                     -- 8 bytes (fixed)
[optional TLV trailers]                 -- 0..N x 4-byte TLV records
```

The dominant size is 12 bytes total (3 records: marker + d1 + d2, no trailer).

### SINGLE size distribution (R5121_12823, n = 6,343,221)

| Size (records) | Count | Trailer bytes | What's in the trailer |
|---|---|---|---|
| 3 | 5,728,946 (90.3%) | 0 | bare hit |
| 5 | 393,683 (6.2%) | 8 | 2 TLV records — usually pulse counters `0x05/0x0b` |
| 1 | 103,629 (1.6%) | -4 | marker only, NO hit data (no d1/d2) |
| 4 | 68,272 (1.1%) | 4 | 1 TLV — typically timing reference `0x1e` |
| 6, 7, 8+ | 36,691 (0.6%) | 12+ | multiple TLVs, can include timing channel updates `0x3b..0x4f` |

### Trailer contents are STANDARD TLV records, not hit data

Every 4th byte of the trailer is dominated by tags from the well-known
header/STATUS dictionary (top-5 mass ≥ 0.95):

| Trailer slot | Top tags observed |
|---|---|
| Trailer record 1 (byte 11) | `0x05`, `0x09`, `0x1e`, `0x3b`, `0x4c` |
| Trailer record 2 (byte 15) | `0x0b`, `0x4c`, `0x3c`, `0x20`, `0x49` |
| Trailer record 3 (byte 19) | `0x3c`, `0x4d`, `0x05`, `0x09` |
| ...further slots | timing-channel tags `0x3b..0x4f` and counters |

This means: **the variable-precision residual encoding hypothesis was wrong.**
SINGLE events embed standard TLV state-update records (pulse counters, voltage,
timing channel updates) that the instrument needed to record alongside that
hit. The actual hit data is ALWAYS exactly 8 bytes (d1 + d2). The remaining
bit-packing problem is therefore a **fixed 64-bit hit-data block**, not a
variable-length encoding.

### Hit data fingerprint (d1, d2) — same across all sizes

Per-byte stats over 50,000 sampled SINGLE events of size=3:

| Byte | mean | std | nUniq | entropy | interpretation |
|---|---|---|---|---|---|
| d1.b0 | 127 | 74 | 256 | 8.00 | high-entropy data |
| d1.b1 | 128 | 74 | 256 | 8.00 | high-entropy data |
| d1.b2 | 42 | 30 | 153 | 5.57 | constrained data |
| d1.b3 | 5 | 7 | 34 | 2.17 | **multiples of 4 (top-5 mass 0.96)** — precision/length code |
| d2.b0 | 127 | 74 | 256 | 8.00 | high-entropy data |
| d2.b1 | 128 | 74 | 93 | 5.24 | constrained, ~ 6-bit |
| d2.b2 | 127 | 99 | 173 | 7.18 | bimodal — two-mode structure |
| d2.b3 | 60 | 36 | 65 | 6.00 | small biased set — sub-channel index? |

### "Marker only" SINGLE (size = 1, 103,629 events)

These have NO d1/d2 — just the `0xa4` marker. They likely correspond to
non-hit triggers or events the instrument flagged but did not capture
position/timing data for. Should be excluded from any HITS→EPOS ion alignment.

### a4 marker data (bytes 0-1)

- 16-bit unsigned value
- Range: 8193-11613 (R5121), 8193-8679 (R6001), 8193-9597 (R6006)
- Consistent with **inter-event timing** at ~100-200 kHz pulse rate with typical detection efficiency
- Cumulative sum grows monotonically -- confirmed timing interpretation

### d1 data word

- byte3 distribution: **strongly biased towards multiples of 4** (0x00: 42%, 0x04: 29%, 0x08: 14%, 0x0c: 7%, 0x10: 3.5%, ...) with geometrically decreasing frequency
- The std of the 24-bit value INCREASES with byte3 value (0x00: std=1.49M, 0x20: std=2.12M), confirming byte3 encodes **data magnitude/precision**, NOT a channel selector
- d1 full 32-bit values are uniformly distributed mod 4

### d2 data word

- byte3 distribution: **uniformly distributed** across 0x00-0xFF
- This means byte3 is pure data (part of the encoded value), not a tag
- d1 and d2 are essentially uncorrelated (r=-0.06)

### 16-bit sub-field analysis

Splitting each 4-byte word into two 16-bit halves (lo=bytes 0-1, hi=bytes 2-3):

| Sub-field        | Signed mean | Signed range     | Interpretation           |
|------------------|-------------|------------------|--------------------------|
| a4_16 (a4 lo16)  | +8458      | 8193 to 11613    | Inter-event time (ns?)   |
| d1_lo (d1 lo16)  | -972       | -32763 to +32765 | Signed, fluctuating      |
| d1_hi (d1 hi16)  | +1138      | 6 to 64549       | Mostly positive          |
| d2_lo (d2 lo16)  | +285       | -29234 to +29280 | Signed, fluctuating      |
| d2_hi (d2 hi16)  | +15538     | 431 to 32326     | Always positive, large   |

### Cumulative sum correlations (strong evidence for delta encoding)

Computing cumulative sums of the signed 16-bit sub-fields over 10,000 events:

| Pair          | Correlation |
|---------------|-------------|
| a4_16 vs d1_hi | r = 1.0000 |
| a4_16 vs d2_hi | r = 1.0000 |
| d1_hi vs d2_hi | r = 1.0000 |
| a4_16 vs d1_lo | r = -0.9882 |
| a4_16 vs d2_lo | r = 0.6193  |

Three sub-fields (a4_16, d1_hi, d2_hi) are perfectly correlated in their cumulative sums, confirming they track **timing** (all increase proportionally to experiment time). d1_lo is anti-correlated (tof or position residual). d2_lo is partially correlated.

### Cumulative sum growth rates

| Sub-field | Cumsum after 10K events | Growth rate |
|-----------|------------------------|-------------|
| a4_16     | 84.6M                  | ~8458/event |
| d1_hi     | 11.4M                  | ~1138/event |
| d2_hi     | 155.4M                 | ~15538/event|
| d1_lo     | -9.7M                  | ~-972/event |
| d2_lo     | 2.8M                   | ~285/event  |

The ratio d2_hi/a4_16 ~ 1.84, close to 2.0, consistent with d2_hi encoding a **delay-line sum** (which advances at ~2x the trigger rate).

---

## 6. MULTI HIT Events (byte2 = 0x40)

Gap distribution for multi-hit events:

| Gap | Data words | Hits | Count   |
|-----|-----------|------|---------|
| 1   | 0         | 0    | 171,179 |
| 3   | 2         | 1    | 163,353 |
| 5   | 4         | 2    | 11,540  |
| 7   | 6         | 3    | 533     |

Multi-hit events with gap=1 are empty markers. Gap=2 contains only one
payload record and is incomplete. For gaps >=3, the data word structure
appears identical to SINGLE events (pairs of 4-byte words per ion); even
gaps are treated as complete ion pairs followed by one trailer record.

---

## 6b. Multi-hit byte roles (structural assay)

Within a gap=5 MULTI event, hit1 and hit2 share the same TOF reference
(same trigger pulse) but differ in detector position. Comparing the two
8-byte hit-data blocks partitions the bytes by physical role:

| Byte | hit1 mean ± std | hit2 mean ± std | Role |
|---|---|---|---|
| 0 | 126 ± 74 | 124 ± 76 | full-range, **independent per ion** → position-like |
| 1 | 127 ± 74 | 122 ± 82 | full-range, **independent per ion** → position-like |
| 2 | 47 ± 35 | **7 ± 12** | absolute in hit1, **delta in hit2** → shared / TOF-related |
| 3 | 15 ± 18 | always 0x05 | **mode/flag byte** (see below) |
| 4 | 126 ± 74 | 122 ± 71 | full-range, **independent per ion** → position-like |
| 5 | 129 ± 74 | 142 ± 73 | full-range, **independent per ion** → position-like |
| 6 | 128 ± 95 | **10 ± 30** | absolute in hit1, **delta in hit2** → shared / TOF-related |
| 7 | 56 ± 38 | 15 ± 14 | partial delta in hit2 → shared |

For gap=7 (3 ions per pulse) the cross-difference table confirms hit1
is the **anchor** and hit2/hit3 are both small residuals from it:
`|h3-h2| ≪ |h2-h1| ≈ |h3-h1|` for bytes 2, 6, 7.

### d1.b3 = 0x05 is a "delta-encoded continuation" marker

In gap=5 MULTI events, hit2's d1.b3 is **always** `0x05` (96.7% transition
to it from hit1's value). The 5,750 SINGLE events with d1.b3 = 0x05 in
section 5 share the same byte-distribution fingerprint as multi-hit hit2's
(byte 2 range [0,7], byte 1 entropy < 8). They are likely
delta-encoded continuation hits that the walker classifies as SINGLE
because their `0xa4` marker happens to carry byte2 = 0x20.

### d1.b3 bit decomposition (over 6.24M SINGLE hits)

| Bit | % set | Meaning |
|---|---|---|
| 0 | 0.11% | **delta-encoded flag** (set when hit is a residual from a prior anchor) |
| 1 | 0.00% | unused |
| 2–5 | 38.7 / 22.8 / 6.75 / 0.82% | **4-bit field** (raw values 0..15 × 4); geometric distribution — likely a counter, sub-channel index, or precision code |
| 6–7 | ~0% | unused |

### Tentative byte→role mapping

- **Bytes 0, 1, 4, 5** carry independent per-ion **position residuals**
  (consistent with DL-x and DL-y timing pair differences).
- **Bytes 2 and 6** carry **shared / TOF-related** values that
  delta-encode in multi-hit.
- **Byte 3** is a mode flag + 4-bit field; byte 7 is partial.

## 7. Unsolved: 8-byte hit-data decode

**Definitive negative result (2026-05-02): direct bit-field extraction
from the 8-byte ion payload does NOT recover any STR-decoded timing or
position quantity, even with a perfectly aligned same-run HITS+STR pair.**

### The decisive experiment

Using R5124_00368 (HITS + same-run STR, 18.3M aligned ions, offset 0):

1. `hitsSearchStrBitfields` brute-forced every contiguous bitfield
   extraction (4–24 bits, signed/unsigned, little-/big-endian) against
   STR-derived targets (`sumX`, `sumY`, `sumW`, `tofRaw`, `detxRaw`,
   `detyRaw`, `detwRaw`, `quality`).
2. Top raw correlations: |r| ≈ 0.38 at startBit=18 LE (= bits 2–5 of
   byte 2) against `sumW`, `sumX`, `sumY`, `tofRaw`. Looked promising.
3. **But: identical correlations appear at almost every other startBit**
   (0–17, 19+) for the same targets. A real bit-field at position 18
   would correlate at 18 and *not* at 0, 5, 10, etc.
4. **And: the same field's correlation with `quality` is 0.018.** Quality
   has no time drift; if the field encoded any real ion-physical timing
   it would show up here too.
5. **Detrending kills the signal completely.** Subtracting a `movmean`
   over 50 events from both the candidate field and the STR target
   collapses the correlation from 0.38 to **0.000**.

The 0.38 raw correlation was driven entirely by experiment-time drift
that's reflected in *all* HITS bytes (encoder state) and *all* STR
TOF-related sums (real TOF reference drifts slowly across the run).
There is zero residual ion-physical signal in any contiguous bitfield.

### What this means for the encoding

The 8-byte payload is **not** a fixed bit-packing of detector quantities.
Cameca's encoding is something more involved. Plausible mechanisms:

1. **Predictive / residual encoding.** Each event stores the difference
   from a prediction the decoder maintains in running state (per-channel
   previous values, running average, etc.). The raw bytes look random
   without subtracting the prediction. *This is the most likely model*
   because it matches every observation: bytes 0/1/4/5 look like random
   residuals; bytes 2/6 look smaller in multi-hit hit2/hit3 because the
   prediction is closer; cumulative sums of d1_hi / d2_hi grow linearly
   with experiment time (consistent with a clock-driven predictor).
2. **Lookup / dictionary encoding.** Bytes are indices into a per-run
   table built from header context. No way to find the table without SDK.
3. **Bitwise scrambling.** Bits XORed with a per-event key derived from
   running state — equivalent in effect to (1).
4. **Arithmetic / range coding.** Bytes do not decode independently;
   would need full-file decompression context.

### Evidence catalog (what we tried, what showed nothing)

- **Direct 16-bit channel interpretation** (initial 2026-04 work): all
  pairwise byte correlations against EPOS `mc` < 0.03.
- **Cumulative delta model**: cumsum(d1_hi) / cumsum(d2_hi) grow linearly
  with cumsum(a4_lo16) at r = 1.0000 — confirms clock-driven drift, not a
  channel value. TOF magnitudes wrong by ~10⁶.
- **`d2_u32` direct TOF interpretation**: spectrum correlation with EPOS
  mc r = -0.01.
- **All 120 channel-tag permutations**: best mapping gives 3 of 5
  monotonic channels but no TOF.
- **Mass-spectrum match against EPOS mc histogram peaks** (R5121, 2026-05):
  no candidate produced a recognisable mass spectrum.
- **Cumsum-residual / linear detrend** at W ∈ {100, 1000, 10000, 100000}
  events: residual histograms unimodal, ≤1 peak per channel.
- **Bitfield search vs same-run STR** (R5124, 2026-05): see above.
  Detrended correlation = 0.000 against every STR target. Definitive.

### What still might work

Things NOT yet tried that have a real chance:

1. **Conditional decoding stratified by `d1.b3` mode class**. The d1.b3
   bit decomposition shows bit 0 is a delta flag and bits 2–5 are a 4-bit
   counter / sub-channel index. Maybe the bit-packing of bytes 0,1,2,4,5,6
   depends on the value of bits 2–5 of b3. Run the bitfield search
   separately within each of the 16 mode classes.
2. **Pair-wise prediction inside multi-hits**: first pass complete,
   negative after robust verification. `hitsSearchStrMultiPairs` found an
   apparent `h2.b2` vs next-STR `detyt2` signal (detrended r ≈ 0.35), but
   `hitsVerifyStrMultiCandidate` showed it is driven by abnormal
   continuation-byte / extreme STR-delta outliers. Restricting to
   `h2.b2 <= 64` and central STR deltas collapses raw r to ~0.
3. **Hidden-state residual model**: test predictors based on previous N-ion
   mean, exponentially weighted per-channel state, pulse-counter phase, and
   marker interval rather than the immediately previous STR ion.
4. **Decompile AP Suite**. The C# / C++ DLLs that read HITS contain the
   exact algorithm.

---

## 8. Comparison with STR v2 Format

The existing `strLoad.m` successfully decodes STR v2. Key differences:

| Feature        | STR v2                    | HITS v3                          |
|----------------|---------------------------|----------------------------------|
| Event marker   | Tag `0x18` (quality)      | Tag `0xa4` with type byte        |
| Channel data   | Per-channel TLV tags      | Packed 8-byte records            |
| Channel tags   | `0x01-0x04`, `0x21-0x22`  | Not individual tags in bulk      |
| Compression    | None (raw TLV)            | Variable-precision residual      |
| File size      | ~40 bytes/event           | ~12 bytes/event                  |

---

## 9. Recommended next steps (ranked by yield × cheapness)

1. **Multi-hit pair-wise prediction test**. **First pass complete,
   negative after robust verification.** The nominal best candidate was
   `h2.b2` against next-STR `detyt2` (detrended r ≈ 0.35), but this was
   driven by one extreme STR delta and rare high `h2.b2` values. In the
   physically normal continuation-byte range (`h2.b2 <= 64`) and central
   STR-delta range, raw r collapses to ~0.

2. **Delta correlation against STR previous-ion values** (~2 hours, ~50%
   odds). **First pass complete, negative.** `hitsSearchStrDeltas` tested
   HITS bytes and adjacent 16-/24-bit words against
   `STR_value(i) − STR_value(i−1)` on R5124. Top detrended |r| was ~0.026
   for SINGLE and ~0.042 for STATUS. A simple direct previous-ion residual
   is unlikely. Running-mean / hidden-state residuals remain untested.

3. **Stratified bitfield search by `d1.b3` 4-bit mode class** (~3 hours,
   ~20% odds). Run `hitsSearchStrBitfields` 16 times, once per value of
   `bitand(b3, 0x3c)`. If the encoding is mode-conditional, restricting
   to one mode could expose a real bitfield.

4. **Reach out to FAIRmat / pynxtools-apm developers** (passive, async).
   `pynxtools-apm` and `ifes_apt_tc_data_modeling` are open-source APT
   tooling. Ask if a HITS parser exists or is in progress.

5. **AP Suite SDK / DLL decompilation** (high effort, high yield if
   accessible). The CAMECA reader for HITS exists as a binary; static
   analysis of its bit-extraction routines would directly reveal the
   algorithm.

6. **CAMECA research collaboration**. The format is proprietary but a
   research request through the user's lab may unlock direct
   documentation.

Steps 1–3 can be done with the existing `R5124_00368` pair and the
existing infrastructure. Steps 4–6 require external coordination.

5. **Brute-force bit-packing search**: With a paired HITS+STR or HITS+EPOS/APT
   dataset where raw detector/time targets are populated and row alignment
   passes audit, systematically test conditional bit extractions from the
   64 data bits against known detector/time values. Search by `d1.b3` mode
   class and `isDelta` rather than treating the payload as one uniform
   encoding.

---

## 10. Reference: Test Files

Authoritative path: `/Users/peterfelfer/Dropbox/01_research/0101_shared/010113_data-pool/APT_data/`

| File | Size | Pair | Detector labels | Notes |
|---|---|---|---|---|
| `R5124_00368/R5124_00368.HITS` | 226.6 MB | + same-run STR ✓ | PPfdUD | **Strongest test pair.** 18.3M events. |
| `R5124_00368/R5124_00368.STR` | 522.8 MB | + HITS above | — | v3 header, STR-style TLV bulk. `strLoad` parses via extension. |
| `R5121_12823.HITS` | 112.0 MB | + decimated EPOS | RPilbR | 7.93M events. |
| `R5121_12823.EPOS` | 326.6 MB | + HITS above | — | 7,423,312 ions; only `x, y, z, mc, atomNum` populated; raw fields zeroed. |
| `R5124_00368/...` | — | — | — | Symlinked into project at `exampleFiles/`. |

Other HITS files in the data pool (no STR/raw-EPOS pair found):
R5124_00368, R6001_70985, R6001_70994, R6006_* (12 files).

For STR-only files (older RHIT-era machines, encoding incompatible with HITS), see `data-pool/APT_data/R56_*.STR`.
