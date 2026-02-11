# Visualization API Contract

This document defines the stable, script-facing API for APT scatter visualisation and profile automation.

## Scope

The contract applies to the following public functions:

- `scatterPlotPosWidget`
- `scatterPlotPosWidgetGetState`
- `scatterPlotPosWidgetApplyState`
- `scatterPlotPosWidgetLoadState`
- `visualisationProfileFromWidget`
- `visualisationProfileResolve`
- `visualisationProfileApply`
- `visualisationProfileExportImages`
- `visualisationProfileExportTurntable`
- `configExport`
- `configImport`

## Stable Behavior

### 1) Widget creation and control handles

- `scatterPlotPosWidget` returns:
  - `scatterHandles`
  - `ax`
  - `controlFig`
  - `info`
- The widget state is persisted on target axes (`ax.UserData.scatterPlotPosWidgetState`) when `persistStateToAxis` is enabled.
- The widget can restore persisted axis state when `restoreStateFromAxis` is enabled and the `pos` signature matches.

### 2) State capture and application

- `scatterPlotPosWidgetGetState` returns a scalar struct that can be saved and reapplied.
- `scatterPlotPosWidgetApplyState` accepts the saved state struct and applies it to an active widget.
- `scatterPlotPosWidgetLoadState` accepts either:
  - a workspace variable name containing a state/profile struct, or
  - a struct payload directly.

### 3) Profile workflow

- `visualisationProfileFromWidget` converts current widget state to canonical profile schema.
- `visualisationProfileApply` applies a canonical profile to a dataset and returns a resolved profile.
- Profile resolution behavior for unmatched species is controlled via `visualisationProfileResolve` policy options.

### 4) Config serialization

- `configExport` and `configImport` support MAT/JSON/YAML for profile/config structs.
- YAML support is internal and dependency-free (`configYamlExport`, `configYamlImport`).

## Canonical Profile Schema

Canonical profile structs include top-level fields:

- `version`
- `schema` (must be `visualisationProfile`)
- `settings`
- `species`
- `policies`
- `export`
- `meta` (optional metadata; included by migration defaults)

The `species` block is aligned by species name, with fields:

- `name`
- `visible`
- `fraction`
- `markerSize`
- `markerSizeLinked`
- `color` (Nx3 double)

## Compatibility Rules

### Forward compatibility

- New fields may be added to profile/state structs.
- Existing fields in this contract are not removed or repurposed without a documented update.

### Backward compatibility

- Legacy widget state structs are accepted and migrated using `visualisationProfileMigrate` and `visualisationProfileToWidgetState`.
- Missing optional fields are defaulted by migration.

### Species mismatch behavior

- Profiles applied to datasets with different species sets are resolved through `visualisationProfileResolve`.
- Default behavior is to preserve matched species settings and expand to current dataset species.

## API Evolution Policy

- Interface changes are allowed, but they must include:
  1. updates to this contract document,
  2. updates to inline function help for affected functions,
  3. updates to example workflow calls where relevant.

- If a behavior changes semantically (not just internally), the corresponding docs must be updated in the same change.

