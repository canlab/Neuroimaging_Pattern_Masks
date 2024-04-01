## Overview

This is a cortical surface atlas projected to volumetric space using data from 
SpaceTop (N=88), PainGen (N=241) and BMRK5 (N=88) using registration fusion.
The sample is identical to the one used for the Glasser reconstruction and
the Desikan Killiany atlas volumetric projection.

The Destrieux atlas is a medium grained anatomical labeling of cortical gyri
and sulci. The "gyri" and "sulci" are divided based on sulcal curvature. For 
a coarser anatomical atlas see the Desikan-Killiany.

The source anatomical segmentations are individualized meaning probabilities can
be used	to infer the location of the named structure in	new participants and 
new studies. Unlike in the Glasser projection, where the cortical folding is
individualized by interregional boundaries are not, here the interregional boundaries
are also individualized based on participants' cortical folding based features.
That said, the training	sample here is probably	not representative of the
general	population due to selection biases that	result in overrepresentation
of healthy 20-40 year olds, so take it with a grain of salt.

## Methods

For methodological details refer to (METHODS.md)[METHODS.md]

Additionally, for an example script used for the registration fusion, refer 
src/single_subject_registration_fusion.sh

