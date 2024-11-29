These are the HCP groupICAs released as part of the PTN1200 dataset. This release also included subject specific 
maps and timeseries but only the group ICAs are included here. Documentation on this release can be found here:

https://www.humanconnectome.org/storage/app/media/documentation/s1200/HCP1200-DenseConnectome+PTN+Appendix-July2017.pdf

ICAs are available at different counts. Only 15, 25, and 50 ICA decompositions are included here. 100, 200 and 300
dimensional ICAs are also available via the HCP aws s3 bucket. For higher d ICAs or subject specific data please
download this file:

s3://hcp-openaccess/HCP_Resources/GroupAvg/HCP_PTN1200/groupICA_3T_HCP1200_MSMAll.tar.gz

Note, these individualized ICAs are not the same as those used for functional alignment in Glasser. These are simple
dual regression ICAs (one spatial regression of groupICAs against subject resting state volumes to obtaine the
subject specific timeseries followed by temporal regression of said timeseries back on to the resting state timeseries).
Glasser et al. use an iterative weighted dual regression procedure to obtain more individualized maps.

Bogdan
11/14/2024
