# Georgia data provenance

Canonical source: `pysal/libpysal`, example dataset `georgia`.

Upstream files used by this research repository:

- `GData_utm.csv`
- `G_utm.*` polygon shapefile components

Upstream description: socio-economic variables for **159 counties in Georgia, USA (1990)**. The libpysal dataset documentation cites Fotheringham, Brunsdon & Charlton (2002) and the MGWR software work by Oshan et al.

The Oshan et al. (2019) `mgwr` paper uses this dataset to demonstrate GWR/MGWR. Its Georgia example uses:

- response: `PctBach`;
- predictors: `PctFB`, `PctBlack`, `PctRural`;
- projected coordinates: `X`, `Y`.

Do not silently replace this file with the separate `GeorgiaEduc` variants used elsewhere in pyGWRx. This directory is intended to preserve the canonical PySAL/libpysal Georgia benchmark used by the MGWR paper/example.

The actual upstream data files are copied automatically by the reproducibility workflow so that the exact package-provided copy and its shapefile sidecars are preserved.
