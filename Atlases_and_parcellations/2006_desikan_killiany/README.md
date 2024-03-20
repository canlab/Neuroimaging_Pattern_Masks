## Overview

This is a cortical surface atlas projected to volumetric space using data from 
SpaceTop (N=88), PainGen (N=241) and BMRK5 (N=88) using registration fusion.
The sample is identical to the one used for the Glasser reconstruction and
the Destrieux atlas volumetric projection.

The Desikan Killiany atlas is a coarse anatomical labeling of cortical gyri
and sulci. The "gyri" include the sulcal banks. For a finer grained anatomical
atlas see the Destrieux atlas.

The source anatomical segmentations are individualized meaning probabilities can
be used to infer the location of the named structure in new participants and 
new studies. Unlike in the Glasser projection, where the cortical folding is
individualized by interregional boundaries are not, here the interregional boundaries
are also individualized based on participants' cortical folding based features.

## Methods

For methodological details refer to (METHODS.md)[METHODS.md]

Additionally, for an example script used for the registration fusion, refer 
src/single_subject_registration_fusion.sh

