# Contributing to Atom Probe Toolbox

Thank you for your interest in contributing to the Atom Probe Toolbox! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Development Setup](#development-setup)
- [Coding Standards](#coding-standards)
- [Submitting Changes](#submitting-changes)
- [Documentation](#documentation)

---

## Code of Conduct

This project follows standard academic collaboration principles. Please be respectful and constructive in all interactions.

---

## How Can I Contribute?

### Reporting Bugs

Before submitting a bug report:
1. Check the [troubleshooting guide](doc/TROUBLESHOOTING.md)
2. Search existing [GitHub issues](https://github.com/peterfelfer/Atom-Probe-Toolbox/issues)
3. Verify the bug with the latest version

When reporting bugs, include:
- MATLAB version and operating system
- Toolbox version
- Minimal code to reproduce the issue
- Expected vs. actual behavior
- Error messages (full stack trace)

### Suggesting Features

Feature requests are welcome! Please:
1. Check if the feature already exists
2. Search existing issues for similar requests
3. Describe the use case and benefit
4. Provide examples if possible

### Contributing Code

We welcome:
- Bug fixes
- New analysis functions
- Performance improvements
- Documentation improvements
- Test cases

### Contributing Workflows

Share your analysis workflows:
- Create a `.mlx` Live Script
- Include sample data or use existing test data
- Add explanatory comments
- Follow the naming convention: `Workflow_YourTopicName.mlx`

---

## Development Setup

### Prerequisites

```matlab
% Required
- MATLAB R2019b or later
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox

% Recommended for development
- Curve Fitting Toolbox
- Parallel Computing Toolbox
```

### Getting Started

```bash
# Clone the repository
git clone https://github.com/peterfelfer/Atom-Probe-Toolbox.git
cd Atom-Probe-Toolbox

# Create a feature branch
git checkout -b feature/my-new-feature
```

```matlab
% In MATLAB, initialize the toolbox
setupToolbox();

% Verify installation
checkDependencies();

% Run tests
runTests();
```

---

## Coding Standards

### File Organization

```
function_name.m          % Main functions in root or appropriate subfolder
/analysis/              % Analysis functions
/utilities/             % Utility functions
/utilities_Cluster/     % Cluster analysis
/crystallography/       % Crystallography functions
/correlative/           % Correlative microscopy
/utilities_IO/          % I/O functions
```

### Function Header Template

Every function should include a standardized header:

```matlab
function [output1, output2] = functionName(input1, input2, options)
% FUNCTIONNAME Brief one-line description
%
% [output1, output2] = functionName(input1, input2)
% [output1, output2] = functionName(input1, input2, 'option', value)
%
% Detailed description of what the function does. Include information
% about the algorithm, assumptions, or limitations.
%
% INPUT:
%   input1 - Description of first input
%            Additional details if needed
%   input2 - Description of second input
%
% OPTIONS:
%   'optionName' - Description (default: value)
%
% OUTPUT:
%   output1 - Description of first output
%   output2 - Description of second output
%
% EXAMPLES:
%   % Basic usage
%   result = functionName(data, params);
%
%   % With options
%   result = functionName(data, params, 'verbose', true);
%
% SEE ALSO:
%   relatedFunction1, relatedFunction2
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    input1
    input2
    options.optionName = defaultValue
end

% Function implementation...

end
```

### Naming Conventions

| Element | Convention | Example |
|---------|------------|---------|
| Functions | camelCase with prefix | `posCalculateConcentration` |
| Variables | camelCase | `atomPositions` |
| Constants | UPPER_CASE | `MAX_ITERATIONS` |
| Classes | PascalCase | `CrystalOrientation` |

**Function Prefixes:**
- `pos*` - Position/point cloud operations
- `ion*` - Ion-related operations
- `range*` - Range operations
- `patch*` - Mesh/patch operations
- `roi*` - Region of interest
- `cluster*` - Cluster analysis
- `ipf*` - Inverse pole figure
- `hdf5*` - HDF5 operations
- `mass*` - Mass spectrum

### Code Style

```matlab
% Use meaningful variable names
atomPositions = pos(:, 1:3);  % Good
ap = pos(:, 1:3);             % Avoid

% Preallocate arrays
results = zeros(nIterations, 3);
for i = 1:nIterations
    results(i, :) = calculate(data(i));
end

% Use vectorization when possible
distances = sqrt(sum((pos1 - pos2).^2, 2));  % Good
for i = 1:n                                    % Avoid if possible
    distances(i) = norm(pos1(i,:) - pos2(i,:));
end

% Add comments for complex logic
% Calculate the misorientation considering crystal symmetry
% by testing all symmetry-equivalent orientations
for symOp = symmetryOperators
    ...
end

% Use argument validation (R2019b+)
arguments
    pos (:,3) double
    epsilon (1,1) double {mustBePositive}
    options.method char {mustBeMember(options.method, {'fast','accurate'})} = 'fast'
end
```

### Input Validation

Use the `validateInputs` utility or MATLAB's `arguments` block:

```matlab
% Using arguments block (preferred for R2019b+)
arguments
    pos (:,3) double
    ranges table
    options.binWidth (1,1) double {mustBePositive} = 0.1
end

% Using validateInputs utility
validateInputs('pos', pos, 'posTable', 'ranges', ranges, 'ranges');
```

### Error Handling

```matlab
% Use meaningful error identifiers
if isempty(pos)
    error('functionName:emptyInput', 'Position data cannot be empty.');
end

% Provide helpful error messages
if ~ismember(method, validMethods)
    error('functionName:invalidMethod', ...
        'Method ''%s'' not recognized. Valid options: %s', ...
        method, strjoin(validMethods, ', '));
end
```

---

## Submitting Changes

### Pull Request Process

1. **Fork and Branch**
   ```bash
   git checkout -b feature/descriptive-name
   ```

2. **Make Changes**
   - Follow coding standards
   - Add/update tests
   - Update documentation

3. **Test**
   ```matlab
   runTests();
   ```

4. **Commit**
   ```bash
   git add .
   git commit -m "Add: descriptive message of changes"
   ```

   Commit message prefixes:
   - `Add:` - New feature
   - `Fix:` - Bug fix
   - `Update:` - Enhancement to existing feature
   - `Doc:` - Documentation only
   - `Refactor:` - Code restructuring

5. **Push and Create PR**
   ```bash
   git push origin feature/descriptive-name
   ```
   Then create a Pull Request on GitHub.

### Pull Request Checklist

- [ ] Code follows the style guidelines
- [ ] Function has proper documentation header
- [ ] Tests pass (`runTests()`)
- [ ] New functionality has tests
- [ ] Documentation updated if needed
- [ ] CHANGELOG.md updated

---

## Documentation

### Function Documentation

Every public function needs:
1. One-line description (H1 line)
2. Syntax examples
3. INPUT/OUTPUT sections
4. Usage examples
5. See also references

### Workflow Documentation

Workflows should include:
1. Purpose and use case
2. Step-by-step instructions
3. Expected outputs
4. Troubleshooting tips

### Updating Documentation

- Update `README.md` for major features
- Update `doc/API_Reference.md` for new functions
- Update `CHANGELOG.md` for all changes
- Update `doc/helptoc.xml` for new function categories

---

## Testing

### Writing Tests

Create test files with `test_` prefix:

```matlab
% test_myFunction.m
function success = test_basic()
    % Test basic functionality
    input = [1, 2, 3];
    expected = [2, 4, 6];
    result = myFunction(input);
    success = isequal(result, expected);
end

function success = test_empty_input()
    % Test empty input handling
    try
        myFunction([]);
        success = false;  % Should have thrown error
    catch
        success = true;
    end
end
```

### Running Tests

```matlab
% Run all tests
runTests();

% Run specific tests
runTests('pattern', 'test_cluster*.m');

% Verbose output
runTests('verbose', true);
```

---

## Questions?

- Open an issue on GitHub
- Contact the maintainers
- Check the documentation

---

## License

By contributing, you agree that your contributions will be licensed under the GNU General Public License v3.0.

---

**(c) Prof. Peter Felfer Group @ FAU Erlangen-Nürnberg**
