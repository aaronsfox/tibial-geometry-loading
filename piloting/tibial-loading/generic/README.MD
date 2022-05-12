***TODO: update based on changes to loads code***

This folder contains some pilot testing on applying loads to a sample tibia. The tibia (and accompanying fibula) comes from an old Musculoskeletal Atlas Project workshop, whereby these surface meshes were generated from sample data. They are anatomically accurate, but 'smoothed' in a sense given the way they were created (i.e. from a statistical shape model). The data can be found at this Github repository: https://github.com/juzhang/mapclient-fai-workshop.

The points marked on these surfaces as outlined below are to ensure that a tibial coordinate system consistent with Wu et al. (2002) can be created.

The following process was undertaken in Materialise 3matic to prepare the MAP surfaces for processing in Matlab:

- The tibia-fibula MAP output STL was loaded as a single part
- The triangles on the fibula were marked and extracted to create their own part
- The medial ('MM') and lateral ('LM') malleoli were marked
- Uniformly remeshed both bodies with a target triangle edge length of 2.0mm
- The medial ('MC') and lateral ('LC') points on the borders of the tibial condyles were marked.
- The joint surface ('ankleJointEstimate') of the tibia was extracted by marking the triangles on the underside of the tibia and copying these to a new part. These will serve as the nodes to apply the forces to on the tibia.
- The 'tibia', 'fibula', and 'jointSurface' parts were exported as ASCII STL files.
- The identified points for landmarks were exported as curves and primitives in separate .xml files.

Following this the data were processed using the 'applyTibiaLoads.m' code. This code contains relevant commentary on the processes involved. It leverages two Matlab packages for the bulk of the analysis - those being the GIBBON toolbox (https://www.gibboncode.org/) and geom3d (https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d). In particular, the code mostly adapts one of the relevant GIBBON examples ([DEMO_febio_0062_femur_load](https://www.gibboncode.org/html/DEMO_febio_0062_femur_load_01.html)).

## References

Wu et al. (2002). 'ISB recoomendation on definitions of joint coordinate system of various joints for reporting of human joint motion - part I: ankle, hip and spine.' J Biomech, 35: 543-548.

