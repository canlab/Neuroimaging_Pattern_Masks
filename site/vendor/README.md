# Viewer provenance

Copied from CanlabCore revision `82c7085959ad7d2a95fe96c4f8035396464a73a4`:
`CanlabCore/Visualization_functions/canlab_niivue/assets/`.

The CANlab wrapper and CSS are reused unchanged. The vendored NiiVue bundle
is CANlab's pinned 0.57.0 build, including its atlas-outline patch.
CANlab license: `CANLAB_LICENSE`. NiiVue license: `NIIVUE_LICENSE`.

Shared anatomy, atlas labels and HCP surface provenance/checksums are recorded
in `../data/shared/provenance.json`. Refresh shared data explicitly with
`prepare_shared.py /path/to/CanlabCore`, then rebuild and validate.
