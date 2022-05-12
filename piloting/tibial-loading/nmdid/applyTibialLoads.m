%% Import surfaces

%Tibia
%Cortical
[tibiaCorticalSTLstruct] = import_STL('102480-tibia-cortical-remesh.stl');
tibiaCorticalF = tibiaCorticalSTLstruct.solidFaces{1}; %Faces
tibiaCorticalV = tibiaCorticalSTLstruct.solidVertices{1}; %Vertices
[tibiaCorticalF,tibiaCorticalV] = mergeVertices(tibiaCorticalF,tibiaCorticalV);
%Trabecular
[tibiaTrabecularSTLstruct] = import_STL('102480-tibia-trabecular.stl');
tibiaTrabecularF = tibiaTrabecularSTLstruct.solidFaces{1}; %Faces
tibiaTrabecularV = tibiaTrabecularSTLstruct.solidVertices{1}; %Vertices
[tibiaTrabecularF,tibiaTrabecularV] = mergeVertices(tibiaTrabecularF,tibiaTrabecularV);

%Fibula
[fibulaSTLstruct] = import_STL('102480-fibula.stl');
fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
[fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);

%Joint surface (nodes only)
[jointSurfaceSTLstruct] = import_STL('102480-AJ.stl');
surfaceFaces = jointSurfaceSTLstruct.solidFaces{1}; %Faces
surfaceNodes = jointSurfaceSTLstruct.solidVertices{1}; %Vertices
[surfaceFaces,surfaceNodes] = mergeVertices(surfaceFaces,surfaceNodes);

%Tibial plateau (nodes only)
[tibialPlateauSTLstruct] = import_STL('102480-TP.stl');
plateauFaces = tibialPlateauSTLstruct.solidFaces{1}; %Faces
plateauNodes = tibialPlateauSTLstruct.solidVertices{1}; %Vertices
[plateauFaces,plateauNodes] = mergeVertices(plateauFaces,plateauNodes);

%% Import and create landmarks

%Import tibial landmarks
tibPoints = [{'LC'},{'LM'},{'MC'},{'MM'}];
for pp = 1:length(tibPoints)
    tree = xml_read([tibPoints{pp},'.txt']);
    landmarks.(char(tree.Point.Name)) = tree.Point.Coordinate;% / 1000;
    clear tree
end
clear pp

%Create new landmarks
landmarks.IM = midPoint3d(landmarks.LM,landmarks.MM);
landmarks.IC = midPoint3d(landmarks.LC,landmarks.MC);

%% Reduce mesh size for speed

%%%% Simplified mesh here during testing to speed up simulations...

%Set remesh options
% % % remeshOpt.pointSpacing = 2;
remeshOpt.pointSpacing = 4;
remeshOpt.disp_on = 0;

%Remesh surfaces
[tibiaCorticalF,tibiaCorticalV] = ggremesh(tibiaCorticalF,tibiaCorticalV,remeshOpt);
[tibiaTrabecularF,tibiaTrabecularV] = ggremesh(tibiaTrabecularF,tibiaTrabecularV,remeshOpt);
[fibulaF,fibulaV] = ggremesh(fibulaF,fibulaV,remeshOpt);

%% Create and align tibial coordinate system

%Create the tibial planes
%Frontal
planes.frontal = createPlane(landmarks.IM,...
    landmarks.LC,landmarks.MC);
%Torsional
planes.torsional = createPlane(landmarks.IC,...
    landmarks.MM,landmarks.LM);

%Create transform to get tibia aligned to the global plane
globalTransform = createBasisTransform3d('global',planes.torsional);

%Transform surfaces
for pp = 1:length(tibiaCorticalV)
    tibiaCorticalV(pp,:) = transformPoint3d(tibiaCorticalV(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(tibiaTrabecularV)
    tibiaTrabecularV(pp,:) = transformPoint3d(tibiaTrabecularV(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(plateauNodes)
    plateauNodes(pp,:) = transformPoint3d(plateauNodes(pp,:),globalTransform);    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),globalTransform);
end
clear ff

%Do the secondary rotation around the X-axis to make the tibia vertical
%along the Y-axis

%Identify distance between IC and IM, along XY plane
pt1 = [landmarks.IC(1),landmarks.IC(2)];
pt2 = [landmarks.IM(1),landmarks.IM(2)];
dIC_IM = sqrt((pt2(2) - pt1(2))^2 + (pt2(1) - pt1(1))^2);

%Identify distance of IC from IM along the x-axis
pt3 = [landmarks.IC(1),landmarks.IM(2)];
dIC_IMx = sqrt((pt3(2) - pt1(2))^2 + (pt3(1) - pt1(1))^2);

%Calculate angle to rotate about x-axis
rotAng = asin(dIC_IMx/dIC_IM);

%Create rotation matrix around Z-axis by specified angle (in radians)
rotZ = createRotationOz(rotAng*-1); %-ve for anti-clockwise

%Transform surfaces
for pp = 1:length(tibiaCorticalV)
    tibiaCorticalV(pp,:) = transformPoint3d(tibiaCorticalV(pp,:),rotZ);    
end
clear pp
for pp = 1:length(tibiaTrabecularV)
    tibiaTrabecularV(pp,:) = transformPoint3d(tibiaTrabecularV(pp,:),rotZ);    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),rotZ);    
end
clear pp
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),rotZ);    
end
clear pp
for pp = 1:length(plateauNodes)
    plateauNodes(pp,:) = transformPoint3d(plateauNodes(pp,:),rotZ);    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),rotZ);
end
clear ff

%Create the transform to make the IM landmark the global origin
transMatrix = createTranslation3d([0,0,0] - landmarks.IM);

%Transform surfaces
for pp = 1:length(tibiaCorticalV)
    tibiaCorticalV(pp,:) = transformPoint3d(tibiaCorticalV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(tibiaTrabecularV)
    tibiaTrabecularV(pp,:) = transformPoint3d(tibiaTrabecularV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(plateauNodes)
    plateauNodes(pp,:) = transformPoint3d(plateauNodes(pp,:),transMatrix);    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),transMatrix);
end
clear ff

%The current bodies are aligned so that X is vertical and Y is lateral
%This needs to be shifted to align with the ISB recommendations so that Y
%is vertical and Z is lateral. This can easily be done by a few rotations
%about specific axes

%First, rotate about the z-axis by -90 degrees

%Transform surfaces
for pp = 1:length(tibiaCorticalV)
    tibiaCorticalV(pp,:) = transformPoint3d(tibiaCorticalV(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp
for pp = 1:length(tibiaTrabecularV)
    tibiaTrabecularV(pp,:) = transformPoint3d(tibiaTrabecularV(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp
for pp = 1:length(plateauNodes)
    plateauNodes(pp,:) = transformPoint3d(plateauNodes(pp,:),createRotationOz(deg2rad(-90)));    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOz(deg2rad(-90)));
end
clear ff

%Second, rotate about the y-axis by -90 degrees

%Transform surfaces
for pp = 1:length(tibiaCorticalV)
    tibiaCorticalV(pp,:) = transformPoint3d(tibiaCorticalV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(tibiaTrabecularV)
    tibiaTrabecularV(pp,:) = transformPoint3d(tibiaTrabecularV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(surfaceNodes)
    surfaceNodes(pp,:) = transformPoint3d(surfaceNodes(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(plateauNodes)
    plateauNodes(pp,:) = transformPoint3d(plateauNodes(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOy(deg2rad(-90)));
end
clear ff

%% Visualise imported and rotated surfaces

%Visualise
cFigure; hold on
gpatch(tibiaCorticalF,tibiaCorticalV,'gw','k');
gpatch(fibulaF,fibulaV,'bw','k');
axisGeom; camlight headlight
title('Imported and Aligned Surfaces');

%% Join and mesh surfaces

%%%%% IMPORTANT: for multi-region meshing, the order in joining the element
%%%%% sets and regions must be from outside-inside (i.e.
%%%%% cortical-trabecular) --- otherwise it won't work...

%Join surfaces
[F,V,C] = joinElementSets({tibiaCorticalF tibiaTrabecularF},{tibiaCorticalV tibiaTrabecularV});

%Mesh surfaces
[V_region1] = getInnerPoint(F(C==2,:),V); %First interior point
[V_region2] = getInnerPoint({tibiaCorticalF tibiaTrabecularF},{tibiaCorticalV tibiaTrabecularV}); %Second interior point
V_regions = [V_region1; V_region2]; %Collect region points
V_holes = []; %Define hole points
[regionTetVolume1] = tetVolMeanEst(tibiaCorticalF,tibiaCorticalV); %Volume estimate for regular tets
[regionTetVolume2] = tetVolMeanEst(tibiaTrabecularF,tibiaTrabecularV); %Volume estimate for regular tets
regionTetVolumes = [regionTetVolume1 regionTetVolume2];
stringOpt = '-pq1.2AaY'; %Tetgen options

%Create tetgen input structure
inputStruct.stringOpt = stringOpt; %Tetgen options
inputStruct.Faces = F; %Boundary faces
inputStruct.Nodes = V; %Nodes of boundary
inputStruct.faceBoundaryMarker = C;
inputStruct.regionPoints = V_regions; %Interior points for regions
inputStruct.holePoints = V_holes; %Interior points for holes
inputStruct.regionA = regionTetVolumes; %Desired tetrahedral volume for each region

% Mesh model using tetrahedral elements using tetGen
[meshOutput] = runTetGen(inputStruct); %Run tetGen

%Access mesh output structure
E = meshOutput.elements; %The elements
V = meshOutput.nodes; %The vertices or nodes
CE = meshOutput.elementMaterialID; %Element material or region id
Fb = meshOutput.facesBoundary; %The boundary faces
Cb = meshOutput.boundaryMarker; %The boundary markers

%Define material regions in bone elements
%For some reason CE works in -3 for cortical and -2 for trabecular
corticalE = E(CE == -3,:);
trabecularE = E(CE == -2,:);
E = [corticalE; trabecularE];
elementMaterialID = [ones(size(corticalE,1),1); 2*ones(size(trabecularE,1),1);];

% % % %Visualise
% % % cFigure; hold on;
% % % title('Input boundaries','FontSize',15);
% % % hp(1) = gpatch(Fb,V,Cb,'none',0.3);
% % % hp(2) = plotV(V_regions,'r.','MarkerSize',25);
% % % legend(hp,{'Input mesh','Interior point(s)'},'Location','NorthWestOutside');
% % % axisGeom(gca,15); camlight headlight;
% % % colormap(gjet(4)); icolorbar;
% % % 
% % % % Visualizing using |meshView|
% % % meshView(meshOutput);

%% Identify ankle joint contact surface nodes for prescribed force

%Note this process could be performed better, as the joint contact surface
%was manually brushed/labelled in 3matic

%Find the closest nodes on the tibia mesh to the surface nodes
for pp = 1:length(surfaceNodes)
    checkPt = surfaceNodes(pp,:);
    ptDist = distancePoints3d(V(:,:),checkPt);
    jointNode(pp,:) = V(find(ptDist == min(ptDist)),:);
end
clear pp

%Find index on mesh that match joint node points and corresponding faces
logicJointNodes = ismember(V,jointNode,'rows');
logicJointFaces = all(logicJointNodes(Fb),2);
bcForcePrescribeList = unique(Fb(logicJointFaces,:));

%Visualise to confirm
cFigure; hold on;
gpatch(Fb,V,'w','none',1);
plotV(V(bcForcePrescribeList,:),'r.','markerSize',15)
axisGeom; camlight headlight;

%% Work out the force distribution on the ankle joint surface nodes

% SEE: https://www.gibboncode.org/html/DEMO_febio_0062_femur_load_01.html

%Get nodal normal directions
[~,~,N] = patchNormal(fliplr(Fb),V); %Nodal normal directions

%Set total force vectors
%These are bases on some estimates of tibial contact force from Edwards et al.
axialContactF = 13.80;
apContactF = -0.66;
mlContactF = 0.68;
mass = 70.1;
forceTotal = [apContactF * (mass*9.80665), ...
    axialContactF * (mass*9.80665), ...
    mlContactF * (mass*9.80665)];
FX = [forceTotal(1) 0 0]; %X force vector
FY = [0 forceTotal(2) 0]; %Y force vector
FZ = [0 0 forceTotal(3)]; %Z force vector

wx = dot(N(bcForcePrescribeList,:),FX(ones(numel(bcForcePrescribeList),1),:),2);
wy = dot(N(bcForcePrescribeList,:),FY(ones(numel(bcForcePrescribeList),1),:),2);
wz = dot(N(bcForcePrescribeList,:),FZ(ones(numel(bcForcePrescribeList),1),:),2);

%Force zero
wx(wx>0) = 0; wy(wy>0) = 0; wz(wz>0) = 0;

force_X = forceTotal(1).*ones(numel(bcForcePrescribeList),1).*wx;
force_Y = forceTotal(2).*ones(numel(bcForcePrescribeList),1).*wy;
force_Z = forceTotal(3).*ones(numel(bcForcePrescribeList),1).*wz;

force_X = force_X./sum(force_X(:)); %sum now equal to 1
force_X = force_X.*forceTotal(1); %sum now equal to desired

force_Y = force_Y./sum(force_Y(:)); %sum now equal to 1
force_Y = force_Y.*forceTotal(2); %sum now equal to desired

force_Z = force_Z./sum(force_Z(:)); %sum now equal to 1
force_Z = force_Z.*forceTotal(3); %sum now equal to desired

force_joint = [force_X(:) force_Y(:) force_Z(:)];

%% Visualise force distribution

cFigure;
hold on;
title('F_x');
gpatch(Fb,V,'w','none',0.5);
quiverVec([0 0 0],FX,100,'k');
quiverVec(V(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_X);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar;

cFigure;
hold on;
title('F_y');
gpatch(Fb,V,'w','none',0.5);
quiverVec([0 0 0],FY,100,'k');
quiverVec(V(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_Y);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar;

cFigure;
hold on;
title('F_z');
gpatch(Fb,V,'w','none',0.5);
quiverVec([0 0 0],FZ,100,'k');
quiverVec(V(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_Z);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar;

%% Identify a tibial plateau nodes

%Note this process could be performed better, as the joint contact surface
%was manually brushed/labelled in 3matic

%Find the closest nodes on the tibia mesh to the surface nodes
for pp = 1:length(plateauNodes)
    checkPt = plateauNodes(pp,:);
    ptDist = distancePoints3d(V(:,:),checkPt);
    plateauNode(pp,:) = V(find(ptDist == min(ptDist)),:);
end
clear pp

%Find index on mesh that match joint node points and corresponding faces
logicPlateauNodes = ismember(V,plateauNode,'rows');
logicPlateauFaces = all(logicPlateauNodes(Fb),2);
bcPlateauPrescribeList = unique(Fb(logicPlateauFaces,:));

%Visualise to confirm
cFigure; hold on;
gpatch(Fb,V,'w','none',1);
plotV(V(bcPlateauPrescribeList,:),'r.','markerSize',15)
axisGeom; camlight headlight;

%% Identify medial malleolus node to constrain

%Find the closest node on the tibia mesh to the malleoli
%MM
mmDist = distancePoints3d(V, landmarks.MM);
mmPoint = V(find(mmDist == min(mmDist)),:);
% % % %LM
% % % lmDist = distancePoints3d(tibiaV, landmarks.LM);
% % % lmPoint = tibiaV(find(lmDist == min(lmDist)),:);

%Find index in tibia surface that corresponds to points
%MM
logicMedMalNodes = ismember(V,mmPoint,'rows');
logicMedMalFaces = all(logicMedMalNodes(Fb),2);
mmPointID = find(V(:,1) == mmPoint(:,1) & ...
    V(:,2) == mmPoint(:,2) & V(:,3) == mmPoint(:,3));
% % % %LM
% % % logicLatMalNodes = ismember(tibiaMeshV,lmPoint,'rows');
% % % logicLatMalFaces = all(logicLatMalNodes(tibiaMeshFb),2);
% % % lmPointID = find(tibiaMeshV(:,1) == lmPoint(:,1) & ...
% % %     tibiaMeshV(:,2) == lmPoint(:,2) & tibiaMeshV(:,3) == lmPoint(:,3));

%Visualise to confirm
cFigure; hold on;
gpatch(Fb,V,'w','none',1);
plotV(mmPoint,'r.','markerSize',15)
% % % plotV(lmPoint,'b.','markerSize',15)
axisGeom; camlight headlight;

%% Define FEBio input structure

%Get a template with default settings
[febio_spec] = febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version = '3.0';

%Module section
febio_spec.Module.ATTR.type = 'solid';

%FEA control settings
numTimeSteps = 10; %Number of time steps desired
max_refs = 25; %Max reforms
max_ups = 0; %Set to zero to use full-Newton iterations
opt_iter = 6; %Optimum number of iterations
max_retries = 5; %Maximum number of retires
dtmin = (1/numTimeSteps)/100; %Minimum time step size
dtmax = 1/numTimeSteps; %Maximum time step size
runMode = 'external'; %'external' or 'internal'

%File details
febioFebFileNamePart = 'tibiaModel';
febioFebFileName = fullfile(pwd,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName = fullfile(pwd,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp = [febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force = [febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress = [febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stresses
febioLogFileName_strainEnergy = [febioFebFileNamePart,'_energy_out.txt']; %Log file name for exporting strain energy density

%Control section
febio_spec.Control.analysis = 'STATIC';
febio_spec.Control.time_steps = numTimeSteps;
febio_spec.Control.step_size = 1/numTimeSteps;
febio_spec.Control.solver.max_refs = max_refs;
febio_spec.Control.solver.max_ups = max_ups;
febio_spec.Control.time_stepper.dtmin = dtmin;
febio_spec.Control.time_stepper.dtmax = dtmax;
febio_spec.Control.time_stepper.max_retries = max_retries;
febio_spec.Control.time_stepper.opt_iter = opt_iter;

%Material section
%%%% Material parameters adapted from early Edwards et al. work

%Cortical
matName1 = 'corticalMat';
febio_spec.Material.material{1}.ATTR.name = matName1;
febio_spec.Material.material{1}.ATTR.type = 'neo-Hookean';
febio_spec.Material.material{1}.ATTR.id = 1;
febio_spec.Material.material{1}.E = 18600;
febio_spec.Material.material{1}.v = 0.3;

matName2 = 'trabecularMat';
febio_spec.Material.material{2}.ATTR.name = matName2;
febio_spec.Material.material{2}.ATTR.type = 'neo-Hookean';
febio_spec.Material.material{2}.ATTR.id = 2;
febio_spec.Material.material{2}.E = 10400;  %%%% this is a lot higher than GIBBON femur example?
febio_spec.Material.material{2}.v = 0.3;

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name = 'tibiaNodes'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id = (1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL = V; %The nodel coordinates

% -> Elements
partName1 = 'CorticalBone';
febio_spec.Mesh.Elements{1}.ATTR.name = partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type = 'tet4'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id = (1:1:size(corticalE,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL = corticalE; %The element matrix

partName2 = 'TrabecularBone';
febio_spec.Mesh.Elements{2}.ATTR.name = partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type = 'tet4'; %Element type
febio_spec.Mesh.Elements{2}.elem.ATTR.id = size(corticalE,1)+(1:1:size(trabecularE,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL = trabecularE; %The element matrix

% -> NodeSets
forceNodes = 'bcForcePrescribeList';
fixedNodes = 'bcPlateauPrescribeList';
medMalNodes = 'medMalNodes';
% % nodeSetName4 = 'indicesVastusLateralis';
% % nodeSetName5 = 'indicesVastusMedialis';

febio_spec.Mesh.NodeSet{1}.ATTR.name = forceNodes;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id = bcForcePrescribeList(:);

febio_spec.Mesh.NodeSet{2}.ATTR.name = fixedNodes;
febio_spec.Mesh.NodeSet{2}.node.ATTR.id = bcPlateauPrescribeList(:);

febio_spec.Mesh.NodeSet{3}.ATTR.name = medMalNodes;
febio_spec.Mesh.NodeSet{3}.node.ATTR.id = mmPointID(:);

% % % febio_spec.Mesh.NodeSet{4}.ATTR.name = nodeSetName4;
% % % febio_spec.Mesh.NodeSet{4}.node.ATTR.id = bcPrescibeList_VastusLateralis(:);
% % % 
% % % febio_spec.Mesh.NodeSet{5}.ATTR.name = nodeSetName5;
% % % febio_spec.Mesh.NodeSet{5}.node.ATTR.id = bcPrescibeList_VastusMedialis(:);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name = partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat = matName1;

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name = partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat = matName2;

%Boundary condition section
% -> Fix boundary conditions
%Plateau
%Fixed translation
febio_spec.Boundary.bc{1}.ATTR.type = 'fix';
febio_spec.Boundary.bc{1}.ATTR.node_set = fixedNodes;
febio_spec.Boundary.bc{1}.dofs = 'x,y,z';
%Fixed rotation
febio_spec.Boundary.bc{2}.ATTR.type = 'fix';
febio_spec.Boundary.bc{2}.ATTR.node_set = fixedNodes;
febio_spec.Boundary.bc{2}.dofs = 'u,v,w';
%Medial malleolus (anterior-posterior, medial-lateral directions)
febio_spec.Boundary.bc{3}.ATTR.type = 'fix';
febio_spec.Boundary.bc{3}.ATTR.node_set = medMalNodes;
febio_spec.Boundary.bc{3}.dofs = 'x,z';



%MeshData secion
%-> Node data
loadDataName1 = 'force_joint';
febio_spec.MeshData.NodeData{1}.ATTR.name = loadDataName1;
febio_spec.MeshData.NodeData{1}.ATTR.node_set = forceNodes;
febio_spec.MeshData.NodeData{1}.ATTR.datatype = 'vec3';
febio_spec.MeshData.NodeData{1}.node.ATTR.lid = (1:1:numel(bcForcePrescribeList))';
febio_spec.MeshData.NodeData{1}.node.VAL = force_joint;

% % % loadDataName2 = 'force_abductor';
% % % febio_spec.MeshData.NodeData{2}.ATTR.name = loadDataName2;
% % % febio_spec.MeshData.NodeData{2}.ATTR.node_set = nodeSetName3;
% % % febio_spec.MeshData.NodeData{2}.ATTR.datatype = 'vec3';
% % % febio_spec.MeshData.NodeData{2}.node.ATTR.lid = (1:1:numel(bcPrescibeList_abductor))';
% % % febio_spec.MeshData.NodeData{2}.node.VAL = forceAbductor_distributed;
% % % 
% % % loadDataName3 = 'force_VL';
% % % febio_spec.MeshData.NodeData{3}.ATTR.name = loadDataName3;
% % % febio_spec.MeshData.NodeData{3}.ATTR.node_set = nodeSetName4;
% % % febio_spec.MeshData.NodeData{3}.ATTR.datatype = 'vec3';
% % % febio_spec.MeshData.NodeData{3}.node.ATTR.lid = (1:1:numel(bcPrescibeList_VastusLateralis))';
% % % febio_spec.MeshData.NodeData{3}.node.VAL = forceVastusLateralis_distributed;
% % % 
% % % loadDataName4 = 'force_VM';
% % % febio_spec.MeshData.NodeData{4}.ATTR.name = loadDataName4;
% % % febio_spec.MeshData.NodeData{4}.ATTR.node_set = nodeSetName5;
% % % febio_spec.MeshData.NodeData{4}.ATTR.datatype = 'vec3';
% % % febio_spec.MeshData.NodeData{4}.node.ATTR.lid = (1:1:numel(bcPrescibeList_VastusMedialis))';
% % % febio_spec.MeshData.NodeData{4}.node.VAL = forceVastusMedialis_distributed;

%Loads section
% -> Prescribed nodal forces
febio_spec.Loads.nodal_load{1}.ATTR.name = 'jointForce';
febio_spec.Loads.nodal_load{1}.ATTR.type = 'nodal_force';
febio_spec.Loads.nodal_load{1}.ATTR.node_set = forceNodes;
febio_spec.Loads.nodal_load{1}.value.ATTR.lc = 1;
febio_spec.Loads.nodal_load{1}.value.ATTR.type = 'map';
febio_spec.Loads.nodal_load{1}.value.VAL = loadDataName1;

% % % febio_spec.Loads.nodal_load{2}.ATTR.name = 'PrescribedForce2';
% % % febio_spec.Loads.nodal_load{2}.ATTR.type = 'nodal_force';
% % % febio_spec.Loads.nodal_load{2}.ATTR.node_set = nodeSetName3;
% % % febio_spec.Loads.nodal_load{2}.value.ATTR.lc = 1;
% % % febio_spec.Loads.nodal_load{2}.value.ATTR.type = 'map';
% % % febio_spec.Loads.nodal_load{2}.value.VAL = loadDataName2;
% % % 
% % % febio_spec.Loads.nodal_load{3}.ATTR.name = 'PrescribedForce3';
% % % febio_spec.Loads.nodal_load{3}.ATTR.type = 'nodal_force';
% % % febio_spec.Loads.nodal_load{3}.ATTR.node_set = nodeSetName4;
% % % febio_spec.Loads.nodal_load{3}.value.ATTR.lc = 1;
% % % febio_spec.Loads.nodal_load{3}.value.ATTR.type = 'map';
% % % febio_spec.Loads.nodal_load{3}.value.VAL = loadDataName3;
% % % 
% % % febio_spec.Loads.nodal_load{4}.ATTR.name = 'PrescribedForce4';
% % % febio_spec.Loads.nodal_load{4}.ATTR.type = 'nodal_force';
% % % febio_spec.Loads.nodal_load{4}.ATTR.node_set = nodeSetName5;
% % % febio_spec.Loads.nodal_load{4}.value.ATTR.lc = 1;
% % % febio_spec.Loads.nodal_load{4}.value.ATTR.type = 'map';
% % % febio_spec.Loads.nodal_load{4}.value.VAL = loadDataName2;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id = 1;
febio_spec.LoadData.load_controller{1}.ATTR.type = 'loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate = 'LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL = [0 0; 1 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file = febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file = febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data = 'ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim = ',';

febio_spec.Output.logfile.element_data{1}.ATTR.file = febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data = 's1;s2;s3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim = ',';

febio_spec.Output.logfile.element_data{2}.ATTR.file = febioLogFileName_strainEnergy;
febio_spec.Output.logfile.element_data{2}.ATTR.data = 'sed';
febio_spec.Output.logfile.element_data{2}.ATTR.delim = ',';

%% Export the FEBio file

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Run the FEBio analysis

%Settings
febioAnalysis.run_filename = febioFebFileName; %The input file name
febioAnalysis.run_logname = febioLogFileName; %The name for the log file
febioAnalysis.disp_on = 1; %Display information on the command window
febioAnalysis.runMode = runMode;

%Run
[runFlag] = runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!


%%%% The above approach seems to do OK! -- still could reduce some of the
%%%% boundary constraint on the tibial end (singular node perhaps?) and
%%%% apply the ankle force to a better area -- but the stress pattern seems
%%%% to be more appropriate in the .xplt file!

%% Import FEBio results

%Note this just comes from the example and could be done in a somewhat more
%efficient manner if desired -- a differnt variable could be plotted too

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacements from log file
    [time_mat, N_disp_mat,~] = importFEBio_logfile(fullfile(pwd,febioLogFileName_disp)); %Nodal displacement
    time_mat = [0; time_mat(:)]; %Time
    
    %Get nodal displacements
    N_disp_mat = N_disp_mat(:,2:end,:);
    sizImport = size(N_disp_mat);
    sizImport(3) = sizImport(3)+1;
    N_disp_mat_n = zeros(sizImport);
    N_disp_mat_n(:,:,2:end) = N_disp_mat;
    N_disp_mat = N_disp_mat_n;
    DN = N_disp_mat(:,:,end);
    DN_magnitude = sqrt(sum(DN(:,3).^2,2));
    V_DEF = N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    
    %Import strain energy density from log file
    [~,E_energy,~] = importFEBio_logfile(fullfile(pwd,febioLogFileName_strainEnergy)); %Element strain energy
% % %     [~,E_energy,~] =
% importFEBio_logfile(fullfile(pwd,febioLogFileName_strainEnergy)); %stress

    %Remove nodal index column
    E_energy = E_energy(:,2:end,:);

    %Add initial state i.e. zero energy
    sizImport = size(E_energy);
    sizImport(3) = sizImport(3)+1;
    E_energy_mat_n = zeros(sizImport);
    E_energy_mat_n(:,:,2:end) = E_energy;
    E_energy = E_energy_mat_n;
    
    %Convert to appropriate data for visualisation
    [FE_face,C_energy_face] = element2patch(E,E_energy(:,:,end),'tet4');
    [CV] = faceToVertexMeasure(FE_face,V,C_energy_face);
    [indBoundary] = tesBoundary(FE_face,V);
    tibiaMeshFb = FE_face(indBoundary,:);
    
    %Plot the simulated results using anim8
    axLim = [min(min(V_DEF,[],3),[],1); max(max(V_DEF,[],3),[],1)];

    % Create basic view and store graphics handle to initiate animation
    hf = cFigure; %Open figure
    title('Strain energy density')
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1 = gpatch(tibiaMeshFb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    hp1.FaceColor = 'Interp';
    
    %Set colour bar
    axisGeom(gca);
    colormap(gjet(250)); colorbar;
    caxis([0 max(E_energy(:))/25]);
    axis(axLim(:)'); %Set axis limits statically
    camlight headlight;

    % Set up animation features
    animStruct.Time = time_mat; %The time vector
    for qt = 1:1:size(N_disp_mat,3) %Loop over time increments
        
        DN = N_disp_mat(:,:,qt); %Current displacement

        [FE_face,C_energy_face] = element2patch(E,E_energy(:,:,qt),'tet4');
        [CV] = faceToVertexMeasure(FE_face,V,C_energy_face);

        %Set entries in animation structure
        animStruct.Handles{qt} = [hp1 hp1]; %Handles of objects to animate
        animStruct.Props{qt} = {'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt} = {V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
        
    end
    
    %Initiate animation
    anim8(hf,animStruct);
    
end
    
%%

















