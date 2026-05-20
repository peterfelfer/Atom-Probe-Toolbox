#!/usr/bin/env python3
"""Extract data from Cameca RHIT files (CERN ROOT-based format).

Usage:
    python3 rhitExtract.py input.RHIT output.h5

Requires: pip install uproot h5py numpy

The RHIT format is a proprietary extension of CERN ROOT used by Cameca/AMETEK
LEAP atom probes. This script extracts the TTree hit data, histograms, and
decoded instrument parameters using the uproot library and saves them to an
HDF5 file readable by MATLAB.

(c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg
"""

import sys
import struct
import zlib
import numpy as np


def decompress_tkey(raw, offset):
    """Read and decompress a ROOT TKey at the given file offset."""
    raw.seek(offset)
    nbytes = struct.unpack('>i', raw.read(4))[0]
    raw.seek(offset + 6)
    objlen = struct.unpack('>I', raw.read(4))[0]
    raw.seek(offset + 14)
    keylen = struct.unpack('>h', raw.read(2))[0]

    raw.seek(offset + keylen)
    compressed = raw.read(nbytes - keylen)

    chunks = []
    pos = 0
    while pos < len(compressed):
        if compressed[pos:pos + 2] == b'ZL':
            csz = int.from_bytes(compressed[pos + 3:pos + 6], 'little')
            chunk = zlib.decompress(compressed[pos + 9:pos + 9 + csz])
            chunks.append(chunk)
            pos += 9 + csz
        else:
            chunks.append(compressed[pos:])
            break
    return b''.join(chunks)


def decode_run_header(raw, fBEGIN, fEND):
    """Decode CRunHeader to extract instrument parameters."""
    # Scan all TKeys to find CRunHeader with highest cycle number
    offset = fBEGIN
    best_offset = None
    best_cycle = -1

    while offset < fEND:
        raw.seek(offset)
        nbytes_raw = raw.read(4)
        if len(nbytes_raw) < 4:
            break
        nbytes = struct.unpack('>i', nbytes_raw)[0]
        if nbytes == 0:
            break
        if nbytes < 0:
            offset += abs(nbytes)
            continue
        raw.seek(offset)
        header = raw.read(min(nbytes, 200))
        if len(header) < 27:
            offset += nbytes
            continue
        cycle = struct.unpack('>h', header[16:18])[0]
        clen = header[26]
        if 27 + clen <= len(header):
            cn = header[27:27 + clen].decode('ascii', errors='replace')
            if cn == 'CRunHeader' and cycle > best_cycle:
                best_cycle = cycle
                best_offset = offset
        offset += nbytes

    if best_offset is None:
        return {}

    data = decompress_tkey(raw, best_offset)
    rest = data[6:]  # skip bytecount + version

    params = {}

    # Extract IVAS version string
    try:
        for i in range(30, 70):
            if rest[i:i + 1].isdigit():
                end = rest.index(b'\x00', i) if b'\x00' in rest[i:i + 20] \
                    else i + 15
                candidate = rest[i:end].decode('ascii', errors='replace')
                if '.' in candidate and len(candidate) > 5:
                    params['ivas_version'] = candidate.rstrip('\x00')
                    break
    except (ValueError, IndexError):
        pass

    # Extract date string
    for month in [b'Jan', b'Feb', b'Mar', b'Apr', b'May', b'Jun',
                  b'Jul', b'Aug', b'Sep', b'Oct', b'Nov', b'Dec']:
        idx = rest.find(month)
        if idx > 0:
            date_bytes = rest[idx:idx + 20]
            date_str = date_bytes.split(b'\x00')[0].decode('ascii',
                                                           errors='replace')
            params['run_date'] = date_str.strip()
            break

    # Instrument parameters stored as big-endian floats at known offsets
    # (empirically determined from CRunHeader v26 structure)
    float_fields = {
        448: 'detector_halfsize_mm',
        452: 't0_ns',
        456: 'flight_path_mm',
        432: 'max_voltage_V',
        504: 'mcp_gain_voltage_V',
        508: 'anode_accel_voltage_V',
        384: 'detector_param1',
        388: 'detector_param2',
        444: 'detector_param3',
    }

    for offset_in_rest, name in float_fields.items():
        if offset_in_rest + 4 <= len(rest):
            val = struct.unpack('>f',
                                rest[offset_in_rest:offset_in_rest + 4])[0]
            if abs(val) < 1e8 and val == val:
                params[name] = float(val)

    # Bowl correction polynomial coefficients (big-endian doubles)
    coefficients = []
    for i in range(1952, min(2112, len(rest)), 8):
        d = struct.unpack('>d', rest[i:i + 8])[0]
        if d == d and abs(d) < 1e15:
            coefficients.append(d)
    if coefficients:
        params['bowl_correction_coefficients'] = coefficients

    return params


def extract_rhit(input_path, output_path):
    import uproot
    import h5py

    f = uproot.open(input_path)

    # --- TTree 'nth' (main hit data) ---
    tree = None
    for cycle in sorted([int(k.split(';')[1]) for k in f.keys()
                         if k.startswith('nth;')], reverse=True):
        tree = f[f'nth;{cycle}']
        break

    if tree is None:
        print("ERROR: No TTree 'nth' found in RHIT file.", file=sys.stderr)
        sys.exit(1)

    # Identify scalar vs struct branches
    scalar_branches = []
    struct_branches = []
    for bname in tree.keys():
        b = tree[bname]
        if b.typename.startswith('struct'):
            struct_branches.append(bname)
        else:
            scalar_branches.append(bname)

    # --- Decode instrument parameters from CRunHeader ---
    with open(input_path, 'rb') as raw:
        raw.seek(4)
        struct.unpack('>I', raw.read(4))[0]  # fVersion
        fBEGIN = struct.unpack('>I', raw.read(4))[0]
        fEND = struct.unpack('>I', raw.read(4))[0]
        run_params = decode_run_header(raw, fBEGIN, fEND)

    with h5py.File(output_path, 'w') as h5:
        # --- Write scalar hit data ---
        hits_grp = h5.create_group('hits')
        for bname in scalar_branches:
            try:
                data = tree[bname].array(library='np')
                hits_grp.create_dataset(bname, data=data, compression='gzip',
                                        compression_opts=4)
            except Exception as e:
                print(f"Warning: could not read branch '{bname}': {e}",
                      file=sys.stderr)

        # Write struct branches as sub-datasets
        for bname in struct_branches:
            try:
                data = tree[bname].array(library='np')
                grp = hits_grp.create_group(bname)
                for field in data.dtype.names:
                    grp.create_dataset(field, data=data[field],
                                       compression='gzip',
                                       compression_opts=4)
            except Exception as e:
                print(f"Warning: could not read struct branch '{bname}': {e}",
                      file=sys.stderr)

        hits_grp.attrs['num_entries'] = tree.num_entries
        hits_grp.attrs['branch_names'] = ','.join(tree.keys())

        # --- Write histograms ---
        hist_names = {
            'phV;1': 'voltageHistory',
            'phE;1': 'erateHistory',
            'pMass;1': 'massSpectrum',
            'pTofRaw;1': 'tofRaw',
        }
        hist_grp = h5.create_group('histograms')
        for rkey, hname in hist_names.items():
            try:
                hist = f[rkey]
                vals = hist.values()
                edges = hist.axis().edges()
                g = hist_grp.create_group(hname)
                g.create_dataset('values', data=vals)
                g.create_dataset('edges', data=edges)
                g.attrs['title'] = hist.title
                g.attrs['nbins'] = len(vals)
            except Exception as e:
                print(f"Warning: could not read histogram '{rkey}': {e}",
                      file=sys.stderr)

        # 2D detector histogram
        try:
            h2 = f['pXyHist;1']
            g = hist_grp.create_group('detectorXY')
            g.create_dataset('values', data=h2.values(), compression='gzip')
            g.create_dataset('xedges', data=h2.axis(0).edges())
            g.create_dataset('yedges', data=h2.axis(1).edges())
            g.attrs['title'] = h2.title
        except Exception as e:
            print(f"Warning: could not read pXyHist: {e}", file=sys.stderr)

        # --- TNtuple pElf ---
        try:
            elf = f['pElf;1']
            elf_grp = h5.create_group('pElf')
            elf_grp.attrs['num_entries'] = elf.num_entries
            for bname in elf.keys():
                data = elf[bname].array(library='np')
                elf_grp.create_dataset(bname, data=data)
        except Exception as e:
            print(f"Warning: could not read pElf: {e}", file=sys.stderr)

        # --- Instrument parameters ---
        params_grp = h5.create_group('instrumentParams')
        for key, val in run_params.items():
            if key == 'bowl_correction_coefficients':
                params_grp.create_dataset(key, data=np.array(val))
            elif isinstance(val, str):
                params_grp.attrs[key] = val
            else:
                params_grp.attrs[key] = val

        # Detector coordinate conversion factor
        det_half = run_params.get('detector_halfsize_mm', 18.5)
        params_grp.attrs['lsb_to_mm'] = det_half / 750.0

        # --- File-level metadata ---
        h5.attrs['source_file'] = input_path
        h5.attrs['format'] = 'RHIT (Cameca ROOT)'
        h5.attrs['root_keys'] = ','.join(f.keys())

    print(f"Extracted {tree.num_entries} hits to {output_path}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.RHIT output.h5", file=sys.stderr)
        sys.exit(1)
    extract_rhit(sys.argv[1], sys.argv[2])
