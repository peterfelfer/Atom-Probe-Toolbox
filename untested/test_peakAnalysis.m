function results = test_peakAnalysis()
% TEST_PEAKANALYSIS Test suite for peakAnalysis function
%
% results = test_peakAnalysis()
%
% Tests the background-corrected peak analysis function with
% synthetic mass spectrum data.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

fprintf('Running peakAnalysis tests...\n\n');

results = struct();
results.passed = 0;
results.failed = 0;
results.tests = {};

% Test 1: Single peak with known counts
results = runTest(results, @test_singlePeakCounts, 'Single peak count calculation');

% Test 2: Background fitting methods
results = runTest(results, @test_backgroundMethods, 'Background fitting methods');

% Test 3: Multiple peaks
results = runTest(results, @test_multiplePeaks, 'Multiple peak analysis');

% Test 4: Normalized spectrum
results = runTest(results, @test_normalizedSpectrum, 'Normalized spectrum handling');

% Test 5: Range table integration
results = runTest(results, @test_rangeTableFraction, 'Range table fraction calculation');

% Test 6: Struct input
results = runTest(results, @test_structInput, 'Struct input format');

% Test 7: Table input
results = runTest(results, @test_tableInput, 'Table input format');

% Summary
fprintf('\n========================================\n');
fprintf('Test Summary: %d passed, %d failed\n', results.passed, results.failed);
fprintf('========================================\n');

end

function results = runTest(results, testFunc, testName)
    try
        success = testFunc();
        if success
            fprintf('[PASS] %s\n', testName);
            results.passed = results.passed + 1;
        else
            fprintf('[FAIL] %s\n', testName);
            results.failed = results.failed + 1;
        end
        results.tests{end+1} = struct('name', testName, 'passed', success, 'error', '');
    catch ME
        fprintf('[ERROR] %s: %s\n', testName, ME.message);
        results.failed = results.failed + 1;
        results.tests{end+1} = struct('name', testName, 'passed', false, 'error', ME.message);
    end
end

%% Test Functions

function success = test_singlePeakCounts()
    % Test single peak with known count values

    % Create synthetic mass spectrum with a Gaussian peak
    mc = 25:0.01:30;
    bgLevel = 10;  % Background counts per bin
    peakCenter = 27;
    peakSigma = 0.15;
    peakHeight = 500;

    % Generate peak + background
    counts = bgLevel + peakHeight * exp(-((mc - peakCenter).^2) / (2 * peakSigma^2));
    counts = round(counts);

    % Create struct input
    massSpec = struct('mc', mc, 'counts', counts);

    % Define peak regions
    peakRegions = [25.5, 26.5, 27.5, 28.5];  % [bgStart, peakStart, peakEnd, bgEnd]

    % Run analysis
    res = peakAnalysis(massSpec, 'peaks', peakRegions, 'showPlot', false);

    % Check results
    success = true;

    % Peak location should be near 27
    success = success && abs(res.peakLocation - peakCenter) < 0.1;

    % Corrected counts should be close to the integrated Gaussian
    % (approximate, as we're using discrete bins)
    expectedPeakCounts = sum(peakHeight * exp(-((mc(mc >= 26.5 & mc <= 27.5) - peakCenter).^2) / (2 * peakSigma^2)));
    success = success && abs(res.correctedCounts - expectedPeakCounts) / expectedPeakCounts < 0.1;

    % Background counts should be reasonable
    success = success && res.bgCounts > 0;
    success = success && res.bgCounts < res.rawCounts;

    if ~success
        fprintf('  Peak location: expected ~%.1f, got %.2f\n', peakCenter, res.peakLocation);
        fprintf('  Corrected counts: expected ~%.0f, got %.0f\n', expectedPeakCounts, res.correctedCounts);
    end
end

function success = test_backgroundMethods()
    % Test different background fitting methods

    mc = 25:0.01:30;
    bgLevel = 10;
    peakCenter = 27;
    counts = bgLevel + 400 * exp(-((mc - peakCenter).^2) / (2 * 0.15^2));
    counts = round(counts);

    massSpec = struct('mc', mc, 'counts', counts);
    peakRegions = [25.5, 26.5, 27.5, 28.5];

    success = true;

    % Test linear method
    res1 = peakAnalysis(massSpec, 'peaks', peakRegions, 'bgMethod', 'linear', 'showPlot', false);
    success = success && strcmp(res1.bgFit.type, 'linear');

    % Test constant method
    res2 = peakAnalysis(massSpec, 'peaks', peakRegions, 'bgMethod', 'constant', 'showPlot', false);
    success = success && strcmp(res2.bgFit.type, 'constant');

    % Test exponential method
    res3 = peakAnalysis(massSpec, 'peaks', peakRegions, 'bgMethod', 'exponential', 'showPlot', false);
    success = success && (strcmp(res3.bgFit.type, 'exponential') || strcmp(res3.bgFit.type, 'linear'));

    % All methods should give similar results for flat background
    success = success && abs(res1.correctedCounts - res2.correctedCounts) / res1.correctedCounts < 0.2;
end

function success = test_multiplePeaks()
    % Test analysis of multiple peaks

    mc = 20:0.01:40;
    bgLevel = 5;

    % Create two peaks
    peak1Center = 25;
    peak2Center = 35;
    peakHeight = 300;
    peakSigma = 0.15;

    counts = bgLevel + ...
        peakHeight * exp(-((mc - peak1Center).^2) / (2 * peakSigma^2)) + ...
        peakHeight * 0.5 * exp(-((mc - peak2Center).^2) / (2 * peakSigma^2));
    counts = round(counts);

    massSpec = struct('mc', mc, 'counts', counts);

    % Define two peak regions
    peakRegions = [
        23.5, 24.5, 25.5, 26.5;  % First peak
        33.5, 34.5, 35.5, 36.5   % Second peak
    ];

    res = peakAnalysis(massSpec, 'peaks', peakRegions, 'showPlot', false);

    success = true;

    % Should have two results
    success = success && length(res) == 2;

    % First peak should have more counts
    success = success && res(1).correctedCounts > res(2).correctedCounts;

    % Peak locations should be correct
    success = success && abs(res(1).peakLocation - peak1Center) < 0.1;
    success = success && abs(res(2).peakLocation - peak2Center) < 0.1;

    if ~success
        fprintf('  Number of peaks: expected 2, got %d\n', length(res));
        if length(res) >= 2
            fprintf('  Peak1 location: expected %.1f, got %.2f\n', peak1Center, res(1).peakLocation);
            fprintf('  Peak2 location: expected %.1f, got %.2f\n', peak2Center, res(2).peakLocation);
        end
    end
end

function success = test_normalizedSpectrum()
    % Test handling of normalized spectra

    mc = 25:0.01:30;
    totalAtoms = 1e6;
    binWidth = 0.01;

    % Create normalized spectrum (counts/Da/totalCounts)
    bgLevel = 10 / binWidth / totalAtoms;
    peakHeight = 500 / binWidth / totalAtoms;
    peakCenter = 27;

    normalized = bgLevel + peakHeight * exp(-((mc - peakCenter).^2) / (2 * 0.15^2));

    massSpec = struct('mc', mc, 'normalized', normalized, 'totalAtoms', totalAtoms);
    peakRegions = [25.5, 26.5, 27.5, 28.5];

    res = peakAnalysis(massSpec, 'peaks', peakRegions, 'normalized', true, 'showPlot', false);

    success = true;

    % Should still find the peak
    success = success && abs(res.peakLocation - peakCenter) < 0.1;

    % Should have positive counts
    success = success && res.correctedCounts > 0;

    if ~success
        fprintf('  Peak location: expected %.1f, got %.2f\n', peakCenter, res.peakLocation);
        fprintf('  Corrected counts: %.0f\n', res.correctedCounts);
    end
end

function success = test_rangeTableFraction()
    % Test fraction calculation with range table

    mc = 20:0.01:40;
    bgLevel = 5;

    % Create spectrum with multiple peaks
    counts = bgLevel + ...
        1000 * exp(-((mc - 25).^2) / (2 * 0.15^2)) + ...  % Main peak
        500 * exp(-((mc - 30).^2) / (2 * 0.15^2)) + ...   % Second peak
        200 * exp(-((mc - 35).^2) / (2 * 0.15^2));        % Third peak
    counts = round(counts);

    massSpec = struct('mc', mc, 'counts', counts);

    % Create range table covering all peaks
    rangeTable = table();
    rangeTable.mcbegin = [24.5; 29.5; 34.5];
    rangeTable.mcend = [25.5; 30.5; 35.5];
    rangeTable.ion = {'Peak1'; 'Peak2'; 'Peak3'};

    % Analyze first peak
    peakRegions = [23.5, 24.5, 25.5, 26.5];

    res = peakAnalysis(massSpec, 'peaks', peakRegions, 'rangeTable', rangeTable, 'showPlot', false);

    success = true;

    % Should have both fractionTotal and fractionRanged
    success = success && ~isnan(res.fractionRanged);

    % fractionRanged should be larger than fractionTotal
    % (since ranged atoms are a subset of total atoms)
    success = success && res.fractionRanged > res.fractionTotal;

    % fractionRanged should be roughly 1000/(1000+500+200) ≈ 0.59
    expectedFraction = 1000 / (1000 + 500 + 200);
    success = success && abs(res.fractionRanged - expectedFraction) < 0.1;

    if ~success
        fprintf('  fractionRanged: expected ~%.2f, got %.2f\n', expectedFraction, res.fractionRanged);
        fprintf('  fractionTotal: %.4f\n', res.fractionTotal);
    end
end

function success = test_structInput()
    % Test struct input format

    mc = 25:0.01:30;
    counts = 10 + 300 * exp(-((mc - 27).^2) / (2 * 0.15^2));
    counts = round(counts);

    % Test with 'counts' field
    massSpec1 = struct('mc', mc, 'counts', counts);
    res1 = peakAnalysis(massSpec1, 'peaks', [25.5, 26.5, 27.5, 28.5], 'showPlot', false);

    % Test with 'cts' field (alternative naming)
    massSpec2 = struct('mc', mc, 'cts', counts);
    res2 = peakAnalysis(massSpec2, 'peaks', [25.5, 26.5, 27.5, 28.5], 'showPlot', false);

    success = true;

    % Both should work and give same results
    success = success && abs(res1.correctedCounts - res2.correctedCounts) < 1;
    success = success && abs(res1.peakLocation - res2.peakLocation) < 0.01;
end

function success = test_tableInput()
    % Test table input format

    mc = (25:0.01:30)';
    counts = round(10 + 300 * exp(-((mc - 27).^2) / (2 * 0.15^2)));

    massSpec = table(mc, counts, 'VariableNames', {'mc', 'counts'});

    res = peakAnalysis(massSpec, 'peaks', [25.5, 26.5, 27.5, 28.5], 'showPlot', false);

    success = true;

    % Should find the peak
    success = success && abs(res.peakLocation - 27) < 0.1;
    success = success && res.correctedCounts > 0;
end
