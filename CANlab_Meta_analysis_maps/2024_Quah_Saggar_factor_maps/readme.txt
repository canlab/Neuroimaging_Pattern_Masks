2025_Quah_Saggar_factor_maps

These files come from Quah et al. 2025, from Manish Saggar's lab, and contain data-driven latent factors in task-fMRI activation maps.

Full citation:
Quah, S. K. L., Jo, B., Geniesse, C., Uddin, L. Q., Mumford, J. A., Barch, D. M., Fair, D. A., Gotlib, I. H., Poldrack, R. A., & Saggar, M. (2025). A data-driven latent variable approach to validating the research domain criteria framework. Nature Communications, 16, 830. https://doi.org/10.1038/s41467-025-55831-z

The visualizations and Neurosynth topic labels were created using CANlab Core Tools, a toolbox for visualization and analysis of neuroimaging data:
https://github.com/canlab/canlabcore

Contents:

- PA*.nii.gz: Data-driven factor maps from the Quah et al. bifactor analysis. These are the individual factor maps for PA1-PA8.

- quah_factor_obj_combined.mat: MATLAB file containing factor_obj_combined, which stores all of the factor maps together in a single CANlab fmri_data object.

- factor_obj_combined.metadata_table: Metadata table attached to the combined fmri_data object. This stores the map/image names together with Neurosynth topic labels for the factors, as well as the original names from the paper.

- factor_names.txt: Human-readable factor labels for F1-F8 and the general factor G, from the original paper.

- figures/: Saved montage and surface rendering images for the PA maps.

- quah_factor_map_montages.m: MATLAB script that loads the PA*.nii.gz maps, renders montages and surfaces, combines the maps into factor_obj_combined, computes spatial correlations, adds metadata_table entries, and saves quah_factor_obj_combined.mat.

- publish_quah_factor_map_report.m: MATLAB script that publishes the montage script into the HTML reports in published_output/.
- published_output/: Published HTML and image outputs from the MATLAB report-generation step.

- quah_factor_map_montages.html: An HTML report showing the figures and Neurosynth topic and term associations

- Quah_2025_Saggar_factor_maps_natcomms.pdf: Local copy of the paper associated with these maps.

Notes:

- The variable name inside quah_factor_obj_combined.mat is factor_obj_combined.
- In the montage script, factor_obj_combined is created by concatenating the individual PA maps into one fmri_data object, making it easier to inspect the full set together.
- The metadata table is populated from the image names, the original factor labels in factor_names.txt, and topic labels derived in the script.
