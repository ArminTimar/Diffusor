# Validation extract data audit

Reviewed `examples/validation/extract.py`, `manifest.json`, all 316 emitted profile CSVs, and the corresponding downloaded source workbook rows. The manifest covers 316 profiles / 11,758 points across eight studies: Lynn 19, Mutch 20, Mourey 56, Sundermeyer Réunion 100, Sundermeyer Eifel 61, Gordeychik 27, Ruprecht 11, and Ruth 22.

## Findings

- **Selected values match the sources.** Recomputed every mapped output field from its recorded `source_row` and source sheet, including unit/scale transforms; all 11,758 rows matched. No duplicate manifest keys, missing source rows, or duplicate `source_row` provenance within/across profiles were found.
- **Mourey’s varying layouts are handled correctly.** The six `KE62-*` sheets use two header layouts: five sheets put Fo and relative distance in columns H/I (0-based 7/8); `KE62-3311F` puts them in M/N (12/13). Header-based selection resolves both layouts; mapped oxide columns also follow the headers. All 56 output profile IDs are preserved as listed in the source sheets, with no cross-sheet key collisions.
- **Mutch Fo conversion is correct.** Source `XFo`, its standard deviation, and `Inital_XFo` are mole fractions. The extract multiplies these by 100; resulting Fo values are on a 0–100 mol% scale. No column offset or scaling mismatch found.
- **Gordeychik grouping is appropriately specific.** Grouping uses the source Sample + Thin section + Grain/profile fields, so distinct profiles/grains do not collapse together. Distances in source mm are converted to µm. No duplicate provenance across groups.
- **Ruprecht mapping is consistent with workbook headers.** The extract selects normalized distance (V), Ni (AK), and Fo (AN) in the normalized data block. Selected values match those source cells.
- **Ruth’s legacy `.xls` mapping and row numbering check out.** The selected distance/Fo/oxide columns match source rows in `3 Data`; the 1-based `source_row` values correctly point to Excel row numbers.

## Review notes for model consumers

- **Duplicate distances are retained with distinct source-row provenance** in 31 Réunion profiles (some examples: `sundermeyer2020_reunion_150915_1_2`, `..._150915_1_4`). These appear as separate source measurements at the same distance; keep them as distinct observations unless a replication protocol explicitly handles coincident positions. They are not accidental duplicate CSV rows.
- **Source order is not normalized.** Distances decrease along 34 profiles (Lynn and Gordeychik examples), matching their archived traverse orientation. Model loaders should use the recorded distance coordinate/orientation deliberately rather than assume every CSV is ascending.
- The 61 Eifel profiles group the workbook’s sequential `LineN` point labels after stripping the sequence prefix, retaining the crystal/traverse suffix (for example `Line1LS-A-2r` → `LS-A-2r`). The resulting profiles preserve increasing point order and do not show a merge of distinct named crystal suffixes.

No extraction-mapping or profile-splitting correction is indicated by this audit. Preserve duplicate-distance rows and source-row provenance in downstream examples.
