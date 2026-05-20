%[text] # **Data Quality Checks**
%[text] **Assessing APT data quality before analysis**
%[text] Before investing time in a full analysis, it is good practice to check the quality of the data. This workflow covers the most common quality checks performed on atom probe datasets: basic statistics, voltage curve inspection, detector hit maps, mass spectrum quality, and coordinate validation.
%[text] This workflow requires an .epos file (for detector coordinates and voltage). A .pos file can be used for a subset of the checks.
%%
%[text] ## 1. Load data and display basic statistics
%[text] Load the data and display summary information: total atom count, coordinate ranges, and mass-to-charge range.
pos = posLoad; % load .epos file for full QC; .pos for basic checks
  %[control:button:0201]{"position":[1,2]}
%%
numAtoms = height(pos);
fprintf('Total atoms: %d\n', numAtoms);
fprintf('X range: %.1f to %.1f nm (span: %.1f nm)\n', min(pos.x), max(pos.x), max(pos.x)-min(pos.x));
fprintf('Y range: %.1f to %.1f nm (span: %.1f nm)\n', min(pos.y), max(pos.y), max(pos.y)-min(pos.y));
fprintf('Z range: %.1f to %.1f nm (span: %.1f nm)\n', min(pos.z), max(pos.z), max(pos.z)-min(pos.z));
fprintf('Mass-to-charge range: %.2f to %.2f Da\n', min(pos.mc), max(pos.mc)); %[output:02020000]
  %[control:button:0202]{"position":[1,2]}
%%
%[text] ## 2. Voltage curve
%[text] The DC voltage (VDC) over the hit sequence reveals the measurement stability. A smooth, monotonically increasing curve indicates stable evaporation. Jumps or drops can indicate specimen fracture, laser drift, or other issues.
%[text] This section requires .epos data (the VDC column).
if ismember('VDC', pos.Properties.VariableNames)
    figure;
    plot(pos.VDC, 'LineWidth', 0.5);
    xlabel('ion hit index');
    ylabel('V_{DC} [V]');
    title('Voltage curve');
    grid on;
else
    warning('VDC column not found. Load an .epos file for voltage curve analysis.');
end %[output:02030000]
  %[control:button:0203]{"position":[1,2]}
%%
%[text] ## 3. Detector hit map (field desorption map)
%[text] The detector hit map reveals pole figures, dead spots, and detector artefacts. Uniform coverage with crystallographic poles indicates good data. Large dead regions or asymmetric patterns may indicate detector problems.
%[text] This section requires .epos data (detx, dety columns).
if ismember('detx', pos.Properties.VariableNames)
    figure;
    FDM = hist3([pos.detx, pos.dety], [256 256]);
    imagesc(FDM'); axis equal; axis tight;
    colorbar;
    title('Detector hit map (FDM)');
    xlabel('det x'); ylabel('det y');
else
    warning('Detector coordinates not found. Load an .epos file for FDM analysis.');
end %[output:02040000]
  %[control:button:0204]{"position":[1,2]}
%%
%[text] ## 4. Mass spectrum overview
%[text] A quick mass spectrum at two resolutions: coarse (0.1 Da) for overall shape and fine (0.01 Da) for peak detail. Check for: well-defined peaks, reasonable background level, and no unusual artefacts.
figure;
subplot(2,1,1);
histogram(pos.mc, 'BinWidth', 0.1, 'EdgeColor', 'none');
set(gca, 'YScale', 'log');
xlabel('mass-to-charge [Da]'); ylabel('counts');
title('Mass spectrum — coarse (0.1 Da bins)');
xlim([0 min(200, max(pos.mc))]);

subplot(2,1,2);
histogram(pos.mc, 'BinWidth', 0.01, 'EdgeColor', 'none');
set(gca, 'YScale', 'log');
xlabel('mass-to-charge [Da]'); ylabel('counts');
title('Mass spectrum — fine (0.01 Da bins)');
xlim([0 min(100, max(pos.mc))]); %[output:02050000]
  %[control:button:0205]{"position":[1,2]}
%%
%[text] ## 5. Coordinate outlier detection
%[text] Check for atoms with coordinates far outside the main point cloud. These outliers can arise from reconstruction artefacts or detector noise and may need to be excluded before analysis.
figure;
subplot(1,3,1);
histogram(pos.x, 100); xlabel('x [nm]'); title('X distribution');
subplot(1,3,2);
histogram(pos.y, 100); xlabel('y [nm]'); title('Y distribution');
subplot(1,3,3);
histogram(pos.z, 100); xlabel('z [nm]'); title('Z distribution'); %[output:02060000]
  %[control:button:0206]{"position":[1,2]}
%%
%[text] ## 6. Multi-hit analysis
%[text] In .epos data, the *multi* column indicates whether a detection event had multiple hits. A high multi-hit fraction can indicate pile-up effects that affect composition accuracy.
%[text] This section requires .epos data (the multi column).
if ismember('multi', pos.Properties.VariableNames)
    multiHits = sum(pos.multi > 1);
    multiHitFraction = multiHits / height(pos) * 100;
    fprintf('Multi-hit events: %d (%.1f%% of total)\n', multiHits, multiHitFraction);

    figure;
    histogram(pos.multi, 'BinMethod', 'integers');
    xlabel('multiplicity'); ylabel('counts');
    title(sprintf('Multi-hit distribution (%.1f%% multi-hits)', multiHitFraction)); %[output:02070000]
else
    warning('Multi column not found. Load an .epos file for multi-hit analysis.');
end
  %[control:button:0207]{"position":[1,2]}
%%
%[text] ## 7. Automated quality metrics
%[text] The function *dataQualityMetrics* computes spatial resolution estimates, density variations, and other quality indicators in one call.
metrics = dataQualityMetrics(pos, 'showPlots', true); %[output:02080000]
disp(metrics.summary);
  %[control:button:0208]{"position":[1,2]}
%%
%[text] ##

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":35}
%---
%[control:button:0201]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0202]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0203]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0204]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0205]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0206]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0207]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0208]
%   data: {"label":"Run","run":"Section"}
%---
