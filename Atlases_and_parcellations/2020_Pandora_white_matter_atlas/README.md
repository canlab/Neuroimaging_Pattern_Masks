# Pandora-WhiteMatterAtlas

## Pandora white matter atlas
Pandora is a new population-based collection of white matter atlases, represented in both volumetric and surface coordinates in a standard space. These atlases are based on 2443 subjects, and include 216 white matter fascicles derived from 6 different state-of-the-art tractography techniques: AFQ, AFQclipped, Recobundles, TractSeg, Tracula, and Xtract. Because these pathways may overlap, these atlases are represented as 4D volumes where each probabalistic segmentation is seperated across the 4th dimension.

These subjects scans are of all healthy adults and were drawn from three datasets: Human Connectome Project (HCP), Baltimore Longitudinal Study of Aging (BLSA), and Vanderbilt University (VU). 

<p align="center">
    <img src="https://github.com/MASILab/Pandora-WhiteMatterAtlas/blob/master/figures/pipeline.png?raw=true")
</p>
    
## Organization
A directory exists for each of the 6 methods. Within each are the corresponding 4D probablistic atlas created using all datasets, a csv file which describes the bundles within the atlas, and a supplementary folder which contains three similar atlas which were created using data from only one of the datasets (HCP, BLSA, VU). Additionally, a similar directory exists for template T1 images which are created from the same subject data as the corresponding atlases. 

    .
    ├── ...
    ├── <Method>                              # One of the 6 tractography methods
    |   ├── <Method>.nii.gz                   # 4D Atlas using all datasets
    |   ├── rh.<Method>.vtk.gz                # Surface reconstruction of the 4D Atlas (right hemisphere)
    |   ├── lh.<Method>.vtk.gz                # Surface reconstruction of the 4D Atlas (left hemisphere)
    │   ├── <Method>_info.csv                 # Description of each bundle in the atlas
    │   ├── supplementary
    |   |   ├── <Method>_<Dataset>.nii.gz     # 4D Atlas using the HCP/BLSA/VU dataset
    |   |   ├── rh.<Method>_<Dataset>.vtk.gz  # Surface reconstruction of the 4D Atlas (right hemisphere) using the HCP/BLSA/VU dataset
    |   |   ├── lh.<Method>_<Dataset>.vtk.gz  # Surface reconstruction of the 4D Atlas (left hemisphere) using the HCP/BLSA/VU dataset
    |   |   └── ...   
    |   └── ...   
    └── ...   

## Download
Please find the latest release using the release tab or at this link https://github.com/MASILab/Pandora-WhiteMatterAtlas/releases.

## Reference
Our manuscript is currently under review: 

Colin B Hansen*, Qi Yang*, Ilwoo Lyu, Francois Rheault, Cailey Kerley, Bramsh Qamar Chandio,Shreyas Fadnavis, Owen Williams, Andrea T. Shafer, Susan M. Resnick, David H. Zald, Laurie Cutting, Warren D Taylor, Brian Boyd, Eleftherios Garyfallidis, Adam W Anderson, Maxime Descoteaux,Bennett A Landman, Kurt G Schilling. “Pandora: 4-D white matter bundle population-based atlasesderived from diffusion MRI fiber tractography”. Scientific Data. Submitted June 2020.

A preprint is available at https://www.biorxiv.org/content/10.1101/2020.06.12.148999v1.full

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
