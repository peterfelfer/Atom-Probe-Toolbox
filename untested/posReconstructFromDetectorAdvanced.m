function [posOut, info, hullRec] = posReconstructFromDetectorAdvanced(pos, options)
% posReconstructFromDetectorAdvanced reconstructs APT data from detector hits.
%
% [posOut, info] = posReconstructFromDetectorAdvanced(pos)
% [posOut, info, hullRec] = posReconstructFromDetectorAdvanced(pos, options)
%
% OPTIONS (name-value via options struct fields):
%   evolutionModel    'voltage' | 'shank' | 'custom' (default: 'voltage')
%   voltageField      'VDC' | 'VP' (default: auto)
%   kf                field factor (default: 10)
%   ICF               image compression factor (default: 1.65)
%   Fevap             evaporation field [V/nm] (default: 65)
%   flightLength      flight length [mm] (default: 110)
%   detEff            detector efficiency (default: 0.82)
%   shankAngle        shank angle [deg] for shank model (default: 6)
%   radius0           apex radius [nm] for shank model (default: from voltage or 50)
%   customRadiusFcn   function handle mapping ion index -> radius [nm]
%   ionVolumeTable    table with columns ion and ionVolume (optional)
%   isotopeTable      isotope table for ionVolumesFromAtomVolumes (optional)
%   defaultIonVolume  fallback volume [nm^3] (default: 1/60.2)
%   reconstructHull   true/false (default: false)
%   hullOptions       struct for scanHull (optional)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    options.evolutionModel (1,1) string = "voltage"
    options.voltageField (1,1) string = ""
    options.kf (1,1) double = 10
    options.ICF (1,1) double = 1.65
    options.Fevap (1,1) double = 65
    options.flightLength (1,1) double = 110
    options.detEff (1,1) double = 0.82
    options.shankAngle (1,1) double = 6
    options.radius0 (1,1) double = NaN
    options.customRadiusFcn = []
    options.ionVolumeTable table = table()
    options.isotopeTable = []
    options.defaultIonVolume (1,1) double = 1/60.2
    options.reconstructHull (1,1) logical = false
    options.hullOptions = struct()
end

validatePosColumns(pos);

detx = pos.detx;
dety = pos.dety;
[ang, rad] = cart2pol(detx, dety);

Adet = pi * max(rad)^2;
flightLength = options.flightLength;
ICF = options.ICF;

omega = resolveIonVolumes(pos, options);

Rspec = resolveRadiusEvolution(pos, omega, Adet, options);

thetaP = atan(rad / flightLength);
theta = thetaP + asin((ICF - 1) .* sin(thetaP));
[zP, d] = pol2cart(theta, Rspec);
[x, y] = pol2cart(ang, d);
zP = Rspec - zP;

dz = omega .* flightLength^2 ./ (options.detEff * Adet * ICF^2 .* Rspec.^2);
cumZ = cumsum(double(dz));
z = cumZ + zP;

posOut = pos;
posOut.x = x;
posOut.y = y;
posOut.z = z;

info = struct();
info.Rspec = Rspec;
info.dz = dz;
info.cumZ = cumZ;
info.omega = omega;
info.Adet = Adet;
info.options = options;

hullRec = [];
if options.reconstructHull
    hullOptions = resolveHullOptions(options.hullOptions);
    nv = structToNameValue(hullOptions);
    hull = scanHull(pos, hullOptions.alpha, hullOptions.numSeg, ...
        hullOptions.nGon, hullOptions.capLoops, nv{:});
    hullRec = reconstructObjects(hull, Rspec, cumZ, flightLength, ICF);
end
end

function validatePosColumns(pos)
    required = {'detx','dety'};
    missing = setdiff(required, pos.Properties.VariableNames);
    if ~isempty(missing)
        error('posReconstructFromDetectorAdvanced:missingColumns', ...
            'pos table missing columns: %s', strjoin(missing, ', '));
    end
end

function omega = resolveIonVolumes(pos, options)
    if ~isempty(options.ionVolumeTable) && istable(options.ionVolumeTable)
        omega = mapIonVolumes(pos, options.ionVolumeTable, options.defaultIonVolume);
        return;
    end

    if ismember('ion', pos.Properties.VariableNames)
        isotopeTable = options.isotopeTable;
        if isempty(isotopeTable)
            isotopeTable = loadIsotopeTable();
        end
        [~, posWithVolumes] = ionVolumesFromAtomVolumes(pos, isotopeTable);
        if ismember('ionVolume', posWithVolumes.Properties.VariableNames)
            omega = posWithVolumes.ionVolume;
            omega(~isfinite(omega)) = options.defaultIonVolume;
            return;
        end
    end

    omega = repmat(options.defaultIonVolume, height(pos), 1);
end

function hullOptions = resolveHullOptions(hullOptions)
    if isempty(hullOptions)
        hullOptions = struct();
    end
    if ~isfield(hullOptions, 'alpha')
        hullOptions.alpha = 3;
    end
    if ~isfield(hullOptions, 'numSeg')
        hullOptions.numSeg = 10;
    end
    if ~isfield(hullOptions, 'nGon')
        hullOptions.nGon = 64;
    end
    if ~isfield(hullOptions, 'capLoops')
        hullOptions.capLoops = 4;
    end
end

function nv = structToNameValue(s)
    nv = {};
    if isempty(s)
        return;
    end
    fields = fieldnames(s);
    for i = 1:numel(fields)
        f = fields{i};
        if ismember(f, {'alpha','numSeg','nGon','capLoops'})
            continue;
        end
        nv{end+1} = f; %#ok<AGROW>
        nv{end+1} = s.(f); %#ok<AGROW>
    end
end

function omega = mapIonVolumes(pos, ionVolumeTable, defaultVolume)
    if ~ismember('ion', pos.Properties.VariableNames)
        error('posReconstructFromDetectorAdvanced:missingIon', ...
            'pos table must contain ion column for ionVolumeTable mapping.');
    end
    if ~all(ismember({'ion','ionVolume'}, ionVolumeTable.Properties.VariableNames))
        error('posReconstructFromDetectorAdvanced:invalidIonVolumeTable', ...
            'ionVolumeTable must contain ion and ionVolume columns.');
    end
    ion = pos.ion;
    ionStr = string(ion);
    tableStr = string(ionVolumeTable.ion);
    [isMatch, idx] = ismember(ionStr, tableStr);
    omega = repmat(defaultVolume, numel(ionStr), 1);
    omega(isMatch) = ionVolumeTable.ionVolume(idx(isMatch));
end

function isotopeTable = loadIsotopeTable()
    fileName = 'isotopeTable_naturalAbundances.mat';
    if exist(fileName, 'file')
        data = load(fileName);
    else
        data = load(fullfile(fileparts(mfilename('fullpath')), '..', fileName));
    end
    fields = fieldnames(data);
    isotopeTable = [];
    for k = 1:numel(fields)
        val = data.(fields{k});
        if istable(val)
            isotopeTable = val;
            break;
        end
    end
    if isempty(isotopeTable)
        error('posReconstructFromDetectorAdvanced:missingIsotopeTable', ...
            'Could not find isotope table in file %s', fileName);
    end
end

function Rspec = resolveRadiusEvolution(pos, omega, Adet, options)
    n = height(pos);
    model = lower(string(options.evolutionModel));
    if model == "voltage"
        V = resolveVoltage(pos, options);
        Rspec = V ./ (options.kf * options.Fevap);
        Rspec = max(Rspec, eps);
        return;
    end

    if model == "custom"
        fcn = options.customRadiusFcn;
        if isempty(fcn)
            error('posReconstructFromDetectorAdvanced:missingCustomFcn', ...
                'customRadiusFcn required for custom model.');
        end
        idx = (1:n)';
        Rspec = fcn(idx);
        Rspec = max(Rspec, eps);
        return;
    end

    if model ~= "shank"
        error('posReconstructFromDetectorAdvanced:invalidModel', ...
            'evolutionModel must be voltage, shank, or custom.');
    end

    radius0 = options.radius0;
    if ~isfinite(radius0)
        if any(ismember({'VDC','VP'}, pos.Properties.VariableNames))
            try
                V = resolveVoltage(pos, options);
                radius0 = max(V(1) / (options.kf * options.Fevap), 1);
            catch
                radius0 = 50;
            end
        else
            radius0 = 50;
        end
    end

    shankAngle = deg2rad(options.shankAngle);
    Rspec = zeros(n, 1);
    zAccum = 0;
    for i = 1:n
        if i == 1
            Rspec(i) = radius0;
        else
            Rspec(i) = radius0 + zAccum * tan(shankAngle);
        end
        dz = omega(i) * options.flightLength^2 / (options.detEff * Adet * options.ICF^2 * Rspec(i)^2);
        zAccum = zAccum + dz;
    end
    Rspec = max(Rspec, eps);
end

function V = resolveVoltage(pos, options)
    if strlength(options.voltageField) > 0
        field = char(options.voltageField);
        if ~ismember(field, pos.Properties.VariableNames)
            error('posReconstructFromDetectorAdvanced:missingVoltageField', ...
                'Voltage field %s not found in pos.', field);
        end
        V = pos.(field);
        return;
    end

    if ismember('VDC', pos.Properties.VariableNames)
        V = pos.VDC;
    elseif ismember('VP', pos.Properties.VariableNames)
        V = pos.VP;
    else
        error('posReconstructFromDetectorAdvanced:missingVoltage', ...
            'No voltage field found in pos (VDC/VP).');
    end
end

function objectsR = reconstructObjects(objects, Rspec, cumZ, flightLength, ICF)
    objectsR = repmat(struct('vertices', [], 'faces', []), numel(objects), 1);
    n = numel(Rspec);
    for o = 1:numel(objects)
        verts = objects(o).vertices;
        idx = round(verts(:,3));
        idx = min(max(idx, 1), n);
        [ang, rad] = cart2pol(verts(:,1), verts(:,2));
        thetaP = atan(rad / flightLength);
        theta = thetaP + asin((ICF - 1) * sin(thetaP));
        [zP, d] = pol2cart(theta, Rspec(idx));
        [objectsR(o).vertices(:,1), objectsR(o).vertices(:,2)] = pol2cart(ang, d);
        zP = Rspec(idx) - zP;
        objectsR(o).vertices(:,3) = cumZ(idx) + zP;
        objectsR(o).faces = objects(o).faces;
    end
end
