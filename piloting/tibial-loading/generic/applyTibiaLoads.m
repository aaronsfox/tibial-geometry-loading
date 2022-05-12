%% Add code notes

%   - Run from current directory

%% Import surfaces

%Tibia
[tibiaSTLstruct] = import_STL('tibia.stl');
tibiaF = tibiaSTLstruct.solidFaces{1}; %Faces
tibiaV = tibiaSTLstruct.solidVertices{1}; %Vertices
[tibiaF,tibiaV] = mergeVertices(tibiaF,tibiaV);

%Fibula
[fibulaSTLstruct] = import_STL('fibula.stl');
fibulaF = fibulaSTLstruct.solidFaces{1}; %Faces
fibulaV = fibulaSTLstruct.solidVertices{1}; %Vertices
[fibulaF,fibulaV] = mergeVertices(fibulaF,fibulaV);

%Joint surface (nodes only)
[jointSurfaceSTLstruct] = import_STL('ankleSurfaceEstimate.stl');
jointFaces = jointSurfaceSTLstruct.solidFaces{1}; %Faces
jointNodes = jointSurfaceSTLstruct.solidVertices{1}; %Vertices
[jointFaces,jointNodes] = mergeVertices(jointFaces,jointNodes);

%Tibial plateau (nodes only)
[tibialPlateauSTLstruct] = import_STL('tibialPlateauEstimate.stl');
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

%% Remesh to edge length of 2mm

%%%% Simplified mesh here during testing to speed up simulations...

%Set remesh options
% % % remeshOpt.pointSpacing = 2;
remeshOpt.pointSpacing = 4;
remeshOpt.disp_on = 0;

%Remesh surfaces
[tibiaF,tibiaV] = ggremesh(tibiaF,tibiaV,remeshOpt);
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
for pp = 1:length(tibiaV)
    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),globalTransform);    
end
clear pp
for pp = 1:length(jointNodes)
    jointNodes(pp,:) = transformPoint3d(jointNodes(pp,:),globalTransform);    
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
for pp = 1:length(tibiaV)
    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),rotZ);    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),rotZ);    
end
clear pp
for pp = 1:length(jointNodes)
    jointNodes(pp,:) = transformPoint3d(jointNodes(pp,:),rotZ);    
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

%Create the transform to make the IC landmark the global origin
transMatrix = createTranslation3d([0,0,0] - landmarks.IC);

%Transform surfaces
for pp = 1:length(tibiaV)
    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(jointNodes)
    jointNodes(pp,:) = transformPoint3d(jointNodes(pp,:),transMatrix);    
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
%This needs to be shifted to align with the coordinate system that we want
%to match Haider et al. (2020), whereby X is medial-lateral, Y is
%anterior-posterior, and Z is vertical. This can easily be done by a few rotations
%about specific axes

%First, rotate about the y-axis by -90 degrees

%Transform surfaces
for pp = 1:length(tibiaV)
    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOy(deg2rad(-90)));    
end
clear pp
for pp = 1:length(jointNodes)
    jointNodes(pp,:) = transformPoint3d(jointNodes(pp,:),createRotationOy(deg2rad(-90)));    
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

%Second, rotate about the z-axis by 90 degrees

%Transform surfaces
for pp = 1:length(tibiaV)
    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),createRotationOz(deg2rad(90)));    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOz(deg2rad(90)));    
end
clear pp
for pp = 1:length(jointNodes)
    jointNodes(pp,:) = transformPoint3d(jointNodes(pp,:),createRotationOz(deg2rad(90)));    
end
clear pp
for pp = 1:length(plateauNodes)
    plateauNodes(pp,:) = transformPoint3d(plateauNodes(pp,:),createRotationOz(deg2rad(90)));    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOz(deg2rad(90)));
end
clear ff

%Third, rotate about the y-axis by 180 degrees

%Transform surfaces
for pp = 1:length(tibiaV)
    tibiaV(pp,:) = transformPoint3d(tibiaV(pp,:),createRotationOy(deg2rad(180)));    
end
clear pp
for pp = 1:length(fibulaV)
    fibulaV(pp,:) = transformPoint3d(fibulaV(pp,:),createRotationOy(deg2rad(180)));    
end
clear pp
for pp = 1:length(jointNodes)
    jointNodes(pp,:) = transformPoint3d(jointNodes(pp,:),createRotationOy(deg2rad(180)));    
end
clear pp
for pp = 1:length(plateauNodes)
    plateauNodes(pp,:) = transformPoint3d(plateauNodes(pp,:),createRotationOy(deg2rad(180)));    
end
clear pp

%Transform landmarks
currLandmarks = fieldnames(landmarks);
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),createRotationOy(deg2rad(180)));
end
clear ff

%% Visualise imported and rotated surfaces

%Visualise
cFigure; hold on
gpatch(tibiaF,tibiaV,'gw','k');
gpatch(fibulaF,fibulaV,'bw','k');
%Axes lines
plot3([0,0,0;50,0,0], [0,0,0;0,0,0], [0,0,0;0,0,0], 'r', 'LineWidth', 2);
plot3([0,0,0;0,0,0], [0,0,0;0,50,0], [0,0,0;0,0,0], 'g', 'LineWidth', 2);
plot3([0,0,0;0,0,0], [0,0,0;0,0,0], [0,0,0;0,0,50], 'b', 'LineWidth', 2);
% % % plotV(tibiaV(1234,:),'r.','MarkerSize',25)
axisGeom;
% % % camlight headlight
title('Imported and Aligned Surfaces');

%% Export rotated and remeshed bodies

%Doing this to test in FEBio

%Tibia
stlStruct.solidNames = {'tibia'}; %names of parts
stlStruct.solidVertices = {tibiaV}; %Vertices
stlStruct.solidFaces = {tibiaF}; %Faces
stlStruct.solidNormals = {[]}; %Face normals (optional)
export_STL_txt('alignedTibia.stl', stlStruct);

%Fibula
stlStruct.solidNames = {'fibula'}; %names of parts
stlStruct.solidVertices = {fibulaV}; %Vertices
stlStruct.solidFaces = {fibulaF}; %Faces
stlStruct.solidNormals = {[]}; %Face normals (optional)
export_STL_txt('alignedFibula.stl', stlStruct);

%% Cut top of tibia surface for fixed boundary conditions

% % % %Set slicing parameters
% % % snapTolerance = mean(patchEdgeLengths(tibiaF,tibiaV))/100;
% % % 
% % % %Set normal direction of cutting plane as along the vertical y-axis
% % % n = [0 1 0];
% % % 
% % % %Set point on the plane as 75% along the length of the tibia
% % % tibiaLength = distancePoints3d(landmarks.IM,landmarks.IC);
% % % p = [0,landmarks.IM(2)+(tibiaLength*0.75),0];
% % % 
% % % %Slice surface (note 3rd color data output is supressed)
% % % [tibiaFc,tibiaVc,~,logicSide,tibiaEb] = triSurfSlice(tibiaF,tibiaV,[],...
% % %     p,n,snapTolerance);
% % % 
% % % %Visualise sliced mesh
% % % cFigure;
% % % subplot(1,2,1); hold on;
% % % hp1=gpatch(tibiaFc(~logicSide,:),tibiaVc,'bw','none',1);
% % % hp2=gpatch(tibiaFc(logicSide,:),tibiaVc,'rw','none',1);
% % % legend([hp1 hp2],{'Surface above plane','Surface below plane'})
% % % axisGeom; axis manual; camlight headlight;
% % % colormap gjet;
% % % 
% % % subplot(1,2,2); hold on;
% % % gpatch(tibiaFc(logicSide,:),tibiaVc,'w','none',1);
% % % gpatch(tibiaFc(~logicSide,:),tibiaVc,'w','none',0.25);
% % % hp1=gpatch(tibiaEb,tibiaVc,'none','b',1,3);
% % % hp2=quiverVec(p,n,50,'k');
% % % legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
% % % axisGeom; axis manual; camlight headlight;
% % % 
% % % %Extract tibia
% % % tibiaModF = tibiaFc(logicSide==1,:);
% % % [tibiaModF,tibiaModV] = patchCleanUnused(tibiaModF,tibiaVc);
% % % 
% % % %Fill the hole mad by the cut
% % % 
% % % %Get the surface boundary
% % % Eb = patchBoundary(tibiaModF,tibiaModV);
% % % indCurve = edgeListToCurve(Eb);
% % % indCurve = indCurve(1:end-1);
% % % 
% % % %Smooth edge boundary
% % % cparSmooth.n = 5;
% % % cparSmooth.Method = 'HC';
% % % [Eb_smooth] = patchSmooth(Eb,tibiaModV(:,[1 2]),[],cparSmooth);
% % % tibiaModV(indCurve,[1 2]) = Eb_smooth(indCurve,:);
% % % cparSmooth.RigidConstraints = indCurve;
% % % [tibiaModV] = patchSmooth(tibiaModF,tibiaModV,[],cparSmooth);
% % % 
% % % %Mesh the hole region
% % % pointSpacing = mean(patchEdgeLengths(tibiaModF,tibiaModV));
% % % [fillF,fillV] = regionTriMesh3D({tibiaModV(indCurve,:)},pointSpacing,0,'linear');
% % % %Flip normals if inverted (i.e. -y)
% % % if dot(mean(patchNormal(fillF,fillV)),[0 -1 0]) > 0
% % %     fillF = fliplr(fillF);
% % % end
% % % 
% % % %Join element sets
% % % [tibiaModF,tibiaModV,tibiaModC] = joinElementSets({tibiaModF,fillF},...
% % %     {tibiaModV,fillV});
% % % [tibiaModF,tibiaModV] = mergeVertices(tibiaModF,tibiaModV);

%% Mesh the tibia and fibula using tetgen

%Find interior point
innerPointT = getInnerPoint(tibiaF,tibiaV);
innerPointF = getInnerPoint(fibulaF,fibulaV);

%Visualise interior point to confirm
% % % cFigure; hold on;
% % % gpatch(tibiaF,tibiaV,'w','none',0.5);
% % % plotV(innerPointT,'r.','MarkerSize',25)
% % % gpatch(fibulaF,fibulaV,'w','none',0.5);
% % % plotV(innerPointF,'b.','MarkerSize',25)
% % % axisGeom; camlight headlight;

%Set mesh parameters
%Tibia
tetVolumeT = tetVolMeanEst(tibiaF,tibiaV); %Volume for regular tets
tetGenStructT.stringOpt = '-pq1.2AaY';
tetGenStructT.Faces = tibiaF;
tetGenStructT.Nodes = tibiaV;
tetGenStructT.holePoints = [];
tetGenStructT.regionPoints = innerPointT; %region points
tetGenStructT.regionA = tetVolumeT;
%Fibula
tetVolumeF = tetVolMeanEst(fibulaF,fibulaV); %Volume for regular tets
tetGenStructF.stringOpt = '-pq1.2AaY';
tetGenStructF.Faces = fibulaF;
tetGenStructF.Nodes = fibulaV;
tetGenStructF.holePoints = [];
tetGenStructF.regionPoints = innerPointF; %region points
tetGenStructF.regionA = tetVolumeF;

%Run tetgen
[meshOutputT] = runTetGen(tetGenStructT);
[meshOutputF] = runTetGen(tetGenStructF);

%Access elements, nodes, and boundary faces
%Tibia
tibiaMeshE = meshOutputT.elements;
tibiaMeshV = meshOutputT.nodes;
tibiaMeshFb = meshOutputT.facesBoundary;
tibiaMeshCb = meshOutputT.boundaryMarker;
tibiaMeshCE = meshOutputT.elementMaterialID;
%Fibula
fibulaMeshE = meshOutputF.elements;
fibulaMeshV = meshOutputF.nodes;
fibulaMeshFb = meshOutputF.facesBoundary;
fibulaMeshCb = meshOutputF.boundaryMarker;
fibulaMeshCE = meshOutputF.elementMaterialID;

% % % %Visualise solid mesh
% % % hFig = cFigure; hold on;
% % % optionStruct.hFig = hFig;
% % % meshView(meshOutputT,optionStruct);
% % % % % % meshView(meshOutputF,optionStruct);
% % % axisGeom;

%% Identify ankle joint contact surface nodes for prescribed force

% % % %Find the closest nodes on the tibia mesh to the surface nodes
% % % for pp = 1:length(jointNodes)
% % %     checkPt = jointNodes(pp,:);
% % %     ptDist = distancePoints3d(tibiaModV(:,:),checkPt);
% % %     jointNode(pp,:) = tibiaModV(find(ptDist == min(ptDist)),:);
% % % end
% % % clear pp
% % % 
% % % %Find index on mesh that match joint node points and corresponding faces
% % % logicJointNodes = ismember(tibiaMeshV,jointNode,'rows');
% % % logicJointFaces = all(logicJointNodes(tibiaMeshFb),2);
% % % bcForcePrescribeList = unique(tibiaMeshFb(logicJointFaces,:));
% % % 
% % % %Visualise to confirm
% % % cFigure; hold on;
% % % gpatch(tibiaMeshFb,tibiaMeshV,'w','none',1);
% % % plotV(tibiaMeshV(bcForcePrescribeList,:),'r.','markerSize',15)
% % % axisGeom; camlight headlight;
% % % 
% % % %Note this process could be performed better, as the joint contact surface
% % % %was manually brushed/labelled in 3matic

%% Identify a tibial plateau nodes

%Note this process could be performed better, as the joint contact surface
%was manually brushed/labelled in 3matic

%Find the closest nodes on the tibia mesh to the surface nodes
for pp = 1:length(plateauNodes)
    checkPt = plateauNodes(pp,:);
    ptDist = distancePoints3d(tibiaV(:,:),checkPt);
    plateauNode(pp,:) = tibiaV(find(ptDist == min(ptDist)),:);
end
clear pp

%Find index on mesh that match joint node points and corresponding faces
logicPlateauNodes = ismember(tibiaMeshV,plateauNode,'rows');
logicPlateauFaces = all(logicPlateauNodes(tibiaMeshFb),2);
bcPlateauPrescribeList = unique(tibiaMeshFb(logicPlateauFaces,:));

% % % %Visualise to confirm
% % % cFigure; hold on;
% % % gpatch(tibiaMeshFb,tibiaMeshV,'w','none',1);
% % % plotV(tibiaMeshV(bcPlateauPrescribeList,:),'r.','markerSize',15)
% % % axisGeom; camlight headlight;

%%%%%

% % % %This point could be found more accurately. Refering to Haider et al.
% % % %(2020) the point for the concentrated force is based on centre of pressure
% % % %identified in a previous study
% % % 
% % % %Note this process could be performed better, as the joint contact surface
% % % %was manually brushed/labelled in 3matic
% % % 
% % % %Calculate the mean of the tibial plateau nodes
% % % meanPlateau = mean(plateauNodes);
% % % 
% % % %Find the closest node on the tibia mesh to the mean of the plateau
% % % plateauDist = distancePoints3d(tibiaV, meanPlateau);
% % % plateauPoint = tibiaV(find(plateauDist == min(plateauDist)),:);
% % % 
% % % %Find index in tibia surface that corresponds to plateau point
% % % logicPlateauNodes = ismember(tibiaMeshV,plateauPoint,'rows');
% % % logicPlateauFaces = all(logicPlateauNodes(tibiaMeshFb),2);
% % % plateauPointID = find(tibiaMeshV(:,1) == plateauPoint(:,1) & ...
% % %     tibiaMeshV(:,2) == plateauPoint(:,2) & tibiaMeshV(:,3) == plateauPoint(:,3));

% % % %Visualise to confirm
% % % cFigure; hold on;
% % % gpatch(tibiaMeshFb,tibiaMeshV,'w','none',1);
% % % plotV(plateauPoint,'r.','markerSize',15)
% % % axisGeom; camlight headlight;

%% Identify central point on ankle joint surface

%This is based on the closest point on the tibial surface to the
%mid-malleoli point

%Find the closest node on the tibia mesh to the mean of the plateau
jointDist = distancePoints3d(tibiaV, landmarks.IM);
jointPoint = tibiaV(find(jointDist == min(jointDist)),:);

%Find index in tibia surface that corresponds to joint point
logicJointNodes = ismember(tibiaMeshV,jointPoint,'rows');
logicJointFaces = all(logicJointNodes(tibiaMeshFb),2);
jointPointID = find(tibiaMeshV(:,1) == jointPoint(:,1) & ...
    tibiaMeshV(:,2) == jointPoint(:,2) & tibiaMeshV(:,3) == jointPoint(:,3));

% % % %Visualise to confirm
% % % cFigure; hold on;
% % % gpatch(tibiaMeshFb,tibiaMeshV,'w','none',1);
% % % plotV(jointPoint,'r.','markerSize',15)
% % % axisGeom; camlight headlight;

%% Identify ankle joint surface nodes to constrain as BC

% % % %This is a basic completely constrained application. This was trialled
% % % %alongside a pinned centre of rotation at the ankle in Haider et al. (2020)
% % % 
% % % %Note this process could be performed better, as the joint contact surface
% % % %was manually brushed/labelled in 3matic
% % % 
% % % %Find the closest nodes on the tibia mesh to the surface nodes
% % % for pp = 1:length(jointNodes)
% % %     checkPt = jointNodes(pp,:);
% % %     ptDist = distancePoints3d(tibiaV(:,:),checkPt);
% % %     jointNode(pp,:) = tibiaV(find(ptDist == min(ptDist)),:);
% % % end
% % % clear pp
% % % 
% % % %Find index on mesh that match joint node points and corresponding faces
% % % logicJointNodes = ismember(tibiaMeshV,jointNode,'rows');
% % % logicJointFaces = all(logicJointNodes(tibiaMeshFb),2);
% % % bcJointConstraintPrescribeList = unique(tibiaMeshFb(logicJointFaces,:));
% % % 
% % % % % % %Visualise to confirm
% % % % % % cFigure; hold on;
% % % % % % gpatch(tibiaMeshFb,tibiaMeshV,'w','none',1);
% % % % % % plotV(tibiaMeshV(bcJointConstraintPrescribeList,:),'r.','markerSize',15)
% % % % % % axisGeom; camlight headlight;

%% Work out force distribution on joint surface

% % % %Note that this is based on surface normal directions. Forces are assumed
% % % %to only be able to act in a compressive sense on the bone, as per the
% % % %GIBBON example.
% % % 
% % % %Get nodal normal directions
% % % [~,~,N] = patchNormal(fliplr(tibiaMeshFb),tibiaMeshV);
% % % 
% % % %Define applied force
% % % %Relatively arbitrary here, but tried to take the axial (Y),
% % % %anterior-posterior (X), and medial-lateral (Z) mean body weight forces
% % % %from Edwards et al. (2010) and convert back to newtons using the mean body
% % % %weight from the same study (70.1kg), using the 4.5m/s running
% % % axialContactF = 13.80;
% % % apContactF = -0.66;
% % % mlContactF = 0.68;
% % % mass = 70.1;
% % % forceTotal = [apContactF * (mass*9.80665), ...
% % %     axialContactF * (mass*9.80665), ...
% % %     mlContactF * (mass*9.80665)];
% % % %%%% Estimate seems somewhat large????
% % % 
% % % %Set force vectors
% % % FX = [forceTotal(1) 0 0]; %X force vector
% % % FY = [0 forceTotal(2) 0]; %Y force vector
% % % FZ = [0 0 forceTotal(3)]; %Z force vector
% % % 
% % % %Set force on prescribe nodes
% % % wx = dot(N(bcForcePrescribeList,:),FX(ones(numel(bcForcePrescribeList),1),:),2);
% % % wy = dot(N(bcForcePrescribeList,:),FY(ones(numel(bcForcePrescribeList),1),:),2);
% % % wz = dot(N(bcForcePrescribeList,:),FZ(ones(numel(bcForcePrescribeList),1),:),2);
% % % 
% % % %Force zero
% % % wx(wx>0) = 0; wy(wy>0) = 0; wz(wz>0) = 0;
% % % 
% % % %Calculate forces
% % % force_X = forceTotal(1).*ones(numel(bcForcePrescribeList),1).*wx;
% % % force_Y = forceTotal(2).*ones(numel(bcForcePrescribeList),1).*wy;
% % % force_Z = forceTotal(3).*ones(numel(bcForcePrescribeList),1).*wz;
% % % 
% % % force_X = force_X./sum(force_X(:)); %sum now equal to 1
% % % force_X = force_X.*forceTotal(1); %sum now equal to desired
% % % 
% % % force_Y = force_Y./sum(force_Y(:)); %sum now equal to 1
% % % force_Y = force_Y.*forceTotal(2); %sum now equal to desired
% % % 
% % % force_Z = force_Z./sum(force_Z(:)); %sum now equal to 1
% % % force_Z = force_Z.*forceTotal(3); %sum now equal to desired
% % % 
% % % %Visualise force distributions
% % % %X axis
% % % cFigure;
% % % subplot(1,3,1);hold on;
% % % title('F_x');
% % % gpatch(tibiaMeshFb,tibiaMeshV,'w','none',0.5);
% % % quiverVec([0 0 0],FX,100,'k');
% % % quiverVec(tibiaMeshV(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_X);
% % % axisGeom; camlight headlight;
% % % colormap(gca,gjet(250)); colorbar;
% % % %Y axis
% % % subplot(1,3,2);hold on;
% % % title('F_y');
% % % gpatch(tibiaMeshFb,tibiaMeshV,'w','none',0.5);
% % % quiverVec([0 0 0],FY,100,'k');
% % % quiverVec(tibiaMeshV(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_Y);
% % % axisGeom; camlight headlight;
% % % colormap(gca,gjet(250)); colorbar;
% % % %Z axis
% % % subplot(1,3,3);hold on;
% % % title('F_z');
% % % gpatch(tibiaMeshFb,tibiaMeshV,'w','none',0.5);
% % % quiverVec([0 0 0],FZ,100,'k');
% % % quiverVec(tibiaMeshV(bcForcePrescribeList,:),N(bcForcePrescribeList,:),10,force_Z);
% % % axisGeom; camlight headlight;
% % % colormap(gca,gjet(250)); colorbar;
% % % 
% % % %%%%% Unsure if forces are being applied in exactly the right directions,
% % % %%%%% so this is something to confirm/check

%% Visualise boundary conditions

% % % %Set support list based on filled upper surface of the tibia
% % % topSupport = tibiaMeshFb(tibiaMeshCb == 2,:);
% % % bcSupportList = unique(topSupport(:));
% % % 
% % % %Visualise
% % % hFig = cFigure; hold on;
% % % gpatch(tibiaMeshFb,tibiaMeshV,'kw','none',0.25);
% % % hl(1) = plotV(tibiaMeshV(bcSupportList,:),'k.','MarkerSize',25);
% % % hl(2) = plotV(tibiaMeshV(bcForcePrescribeList,:),'r.','MarkerSize',25);
% % % legend(hl,{'BC support','BC force prescribe'});
% % % axisGeom; camlight headlight;

%% Define FEBio input structure

%Settings

%File details
febioFebFileNamePart = 'tibiaModel';
febioFebFileName = fullfile(pwd,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName = fullfile(pwd,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp = [febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force = [febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress = [febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stresses
febioLogFileName_strainEnergy = [febioFebFileNamePart,'_energy_out.txt']; %Log file name for exporting strain energy density

%Material parameters (MPa if spatial units are mm)
%In an ideal setting we'd estimate these from images as per Haider et al.
%(2020). Here we take an estimate of axial elastic bone modulus as 17.7GPa,
%converted to 17700MPa --- and apply this value to the directional
%calculations outlined in Haider et al. (2020). Note that their coordinate
%system is aligned to ours, so the values should match.
%These parameters are related to the FEBio "orthotropic elastic" material,
%whereby E1, E2, and E3 are the Young's modulus in the x, y and z
%direction; G12, G23, and G31 are the shear modulus in the xy, yz, and xz
%plane; and v12, v23, v31 are the Poisson's ratio between the x-y, y-z and
%z-x directions. Our data is set-up that the Z-axis is the
%longitudinal/axial direction, Y is the anterior-posterior direction, and X
%is the medial-lateral direction. See section 4.1.3.13 of FEBio User Manual v3.0

%%%% Units vs. force/material properties issue?

%Youngs modulus
E3 = 17700; E1 = 0.574*E3; E2 = 0.577*E3;
%Shear modulus
G12 = 0.195*E3; G23 = 0.265*E3; G31 = 0.216*E3;
%Poisson's ratio
v12 = 0.427; v23 = 0.234; v31 = 0.405;

%Get a template with default settings
[febio_spec] = febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version = '2.5';

%Module section
febio_spec.Module.ATTR.type = 'solid';

%Material section
febio_spec.Material.material{1}.ATTR.type = 'orthotropic elastic';
febio_spec.Material.material{1}.ATTR.id = 1;
%%%% TODO: need mat_axis or is their a default for this? Default may work
%%%% given our materials are in the XYZ axis already?
febio_spec.Material.material{1}.E1 = E1;
febio_spec.Material.material{1}.E2 = E2;
febio_spec.Material.material{1}.E3 = E3;
febio_spec.Material.material{1}.v12 = v12;
febio_spec.Material.material{1}.v23 = v23;
febio_spec.Material.material{1}.v31 = v31;
febio_spec.Material.material{1}.G12 = G12;
febio_spec.Material.material{1}.G23 = G23;
febio_spec.Material.material{1}.G31 = G31;

% % % febio_spec.Material.material{1}.ATTR.type = 'neo-Hookean';
% % % febio_spec.Material.material{1}.ATTR.id = 1;
% % % febio_spec.Material.material{1}.E = 15000;
% % % febio_spec.Material.material{1}.v = 0.3;

%Geometry section

%%%%% TODO: just tibia for now...

% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name = 'nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id = (1:size(tibiaMeshV,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL = tibiaMeshV; %The node coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat = 1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name = 'tibia'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id = (1:1:size(tibiaMeshE,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL = tibiaMeshE;

% % % febio_spec.Geometry.Elements{2}.ATTR.type='tet4'; %Element type of this set
% % % febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set
% % % febio_spec.Geometry.Elements{2}.ATTR.name='CancellousBone'; %Name of the element set
% % % febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
% % % febio_spec.Geometry.Elements{2}.elem.VAL=E2;

% -> NodeSets
%Ankle joint surface constraint
febio_spec.Geometry.NodeSet{1}.ATTR.name = 'ankleJointPoint';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id = jointPointID(:);
%Force application node
febio_spec.Geometry.NodeSet{2}.ATTR.name = 'plateauSurface';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id = bcPlateauPrescribeList(:);
%Create a nodeset for concentrated force at a random node to simulate
%muscle force application
% % % muscPointID = 1234; %this is a random point on the lateral aspect of tibia
% % % febio_spec.Geometry.NodeSet{3}.ATTR.name = 'simMuscleForce';
% % % febio_spec.Geometry.NodeSet{3}.node.ATTR.id = muscPointID(:);

%Boundary condition section

% -> Fix boundary conditions

%Fix the tibial plateaue surface nodes set in the XYZ directions
febio_spec.Boundary.fix{1}.ATTR.bc = 'x';
febio_spec.Boundary.fix{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc = 'y';
febio_spec.Boundary.fix{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc = 'z';
febio_spec.Boundary.fix{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
% % % %Fix the plateau point node set in the XYZ directions
% % % febio_spec.Boundary.fix{4}.ATTR.bc = 'x';
% % % febio_spec.Boundary.fix{4}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
% % % febio_spec.Boundary.fix{5}.ATTR.bc = 'y';
% % % febio_spec.Boundary.fix{5}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
% % % febio_spec.Boundary.fix{6}.ATTR.bc = 'z';
% % % febio_spec.Boundary.fix{6}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;

% -> Force boundary conditions

%Values extracted from Brent's e-mail of JRF's from their speed study
%Took the BW forces for 4.5 m/s from the document and converted to a random
%weight of 75kg
%Note these are total JRFs including muscle force
mlAnkleForce = -0.68 * 75 * 9.80665; apAnkleForce = -0.66 * 75 * 9.80665; axialAnkleForce = 13.80 * 75 * 9.80665;

%Apply XYZ force directions to the tibial plateau point
%Estimates here for the knee contact force from Haider et al. (2020) from
%the highest force conditions. Note again that the coordinate system
%definitions are different
%The forces in the paper are reported in BW so are converted based on a
%radnom weight of 75kg
% % % axialKneeForce = -8.5 * 75 * 9.80665; mlKneeForce = 0.8 * 75 * 9.80665; apKneeForce = 2.5 * 75 * 9.80665;
% % % %Force X
% % % febio_spec.MeshData.NodeData{1}.ATTR.name = 'apKneeForce';
% % % febio_spec.MeshData.NodeData{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
% % % febio_spec.MeshData.NodeData{1}.node.VAL = apKneeForce;
% % % febio_spec.MeshData.NodeData{1}.node.ATTR.lid = (1:1:numel(plateauPointID))'; %just 1 for one point
% % % %Force Y
% % % febio_spec.MeshData.NodeData{2}.ATTR.name = 'axialKneeForce';
% % % febio_spec.MeshData.NodeData{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
% % % febio_spec.MeshData.NodeData{2}.node.VAL = axialKneeForce;
% % % febio_spec.MeshData.NodeData{2}.node.ATTR.lid = (1:1:numel(plateauPointID))'; %just 1 for one point
% % % %Force Z
% % % febio_spec.MeshData.NodeData{3}.ATTR.name='mlKneeForce';
% % % febio_spec.MeshData.NodeData{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{2}.ATTR.name;
% % % febio_spec.MeshData.NodeData{3}.node.VAL = mlKneeForce;
% % % febio_spec.MeshData.NodeData{3}.node.ATTR.lid = (1:1:numel(plateauPointID))'; %just 1 for one point

%Simulated muscle force
%Make this all in the medial to lateral direction
% % % axialMuscForce = 0; mlMuscForce = 1500; apMuscForce = 0;
% % % %Force X
% % % febio_spec.MeshData.NodeData{1}.ATTR.name = 'apMuscForce';
% % % febio_spec.MeshData.NodeData{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{3}.ATTR.name;
% % % febio_spec.MeshData.NodeData{1}.node.VAL = apMuscForce;
% % % febio_spec.MeshData.NodeData{1}.node.ATTR.lid = (1:1:numel(muscPointID))'; %just 1 for one point
% % % %Force Y
% % % febio_spec.MeshData.NodeData{2}.ATTR.name = 'axialMuscForce';
% % % febio_spec.MeshData.NodeData{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{3}.ATTR.name;
% % % febio_spec.MeshData.NodeData{2}.node.VAL = axialMuscForce;
% % % febio_spec.MeshData.NodeData{2}.node.ATTR.lid = (1:1:numel(muscPointID))'; %just 1 for one point
% % % %Force Z
% % % febio_spec.MeshData.NodeData{3}.ATTR.name='mlMuscForce';
% % % febio_spec.MeshData.NodeData{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{3}.ATTR.name;
% % % febio_spec.MeshData.NodeData{3}.node.VAL = mlMuscForce;
% % % febio_spec.MeshData.NodeData{3}.node.ATTR.lid = (1:1:numel(muscPointID))'; %just 1 for one point

%Loads section
%Seems like forces are applied consistently with this approach (i.e. the
%boundary condition force values multiplied by load curve of 1)
% -> Prescribed nodal forces

%Knee contact force
%X-axis
febio_spec.Loads.nodal_load{1}.ATTR.bc = 'x';
febio_spec.Loads.nodal_load{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Loads.nodal_load{1}.scale.ATTR.lc = 1;
febio_spec.Loads.nodal_load{1}.scale.VAL = 1;
febio_spec.Loads.nodal_load{1}.value = mlAnkleForce;
% % % febio_spec.Loads.nodal_load{1}.value.ATTR.node_data = febio_spec.MeshData.NodeData{1}.ATTR.name;
%Y-axis
febio_spec.Loads.nodal_load{2}.ATTR.bc = 'y';
febio_spec.Loads.nodal_load{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Loads.nodal_load{2}.scale.ATTR.lc = 1;
febio_spec.Loads.nodal_load{2}.scale.VAL = 1;
febio_spec.Loads.nodal_load{2}.value = apAnkleForce;
% % % febio_spec.Loads.nodal_load{2}.value.ATTR.node_data = febio_spec.MeshData.NodeData{2}.ATTR.name;
%Z-axis
febio_spec.Loads.nodal_load{3}.ATTR.bc = 'z';
febio_spec.Loads.nodal_load{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Loads.nodal_load{3}.scale.ATTR.lc = 1;
febio_spec.Loads.nodal_load{3}.scale.VAL = 1;
febio_spec.Loads.nodal_load{3}.value = axialAnkleForce;
% % % febio_spec.Loads.nodal_load{3}.value.ATTR.node_data = febio_spec.MeshData.NodeData{3}.ATTR.name;

% % % %Simulated muscle force
% % % %X-axis
% % % febio_spec.Loads.nodal_load{4}.ATTR.bc = 'x';
% % % febio_spec.Loads.nodal_load{4}.ATTR.node_set = febio_spec.Geometry.NodeSet{3}.ATTR.name;
% % % febio_spec.Loads.nodal_load{4}.scale.ATTR.lc = 1;
% % % febio_spec.Loads.nodal_load{4}.scale.VAL = 1;
% % % febio_spec.Loads.nodal_load{4}.value = apMuscForce;
% % % %Y-axis
% % % febio_spec.Loads.nodal_load{5}.ATTR.bc = 'y';
% % % febio_spec.Loads.nodal_load{5}.ATTR.node_set = febio_spec.Geometry.NodeSet{3}.ATTR.name;
% % % febio_spec.Loads.nodal_load{5}.scale.ATTR.lc = 1;
% % % febio_spec.Loads.nodal_load{5}.scale.VAL = 1;
% % % febio_spec.Loads.nodal_load{5}.value = axialMuscForce;
% % % %Z-axis
% % % febio_spec.Loads.nodal_load{6}.ATTR.bc = 'z';
% % % febio_spec.Loads.nodal_load{6}.ATTR.node_set = febio_spec.Geometry.NodeSet{3}.ATTR.name;
% % % febio_spec.Loads.nodal_load{6}.scale.ATTR.lc = 1;
% % % febio_spec.Loads.nodal_load{6}.scale.VAL = 1;
% % % febio_spec.Loads.nodal_load{6}.value = mlMuscForce;

%Load curve
%Set to just be at value of one from beginning
febio_spec.LoadData.loadcurve{1}.ATTR.id = 1;
febio_spec.LoadData.loadcurve{1}.ATTR.type = 'linear';
febio_spec.LoadData.loadcurve{1}.point.VAL = [0 1; 1 1];

%FEA control settings
%%%%% TODO: this is where some edits probably need to be made to reduce the
%%%%% steps and just apply the force in one instance
numTimeSteps = 10; %Number of time steps desired
max_refs = 25; %Max reforms
max_ups = 0; %Set to zero to use full-Newton iterations
opt_iter = 6; %Optimum number of iterations
max_retries = 5; %Maximum number of retires
dtmin = (1/numTimeSteps)/100; %Minimum time step size
dtmax = 1/numTimeSteps; %Maximum time step size
runMode = 'external';%'internal'; %'external' or 'internal'

%Control section
febio_spec.Control.analysis.ATTR.type = 'dynamic';%'static';
febio_spec.Control.time_steps = numTimeSteps;
febio_spec.Control.step_size = 1/numTimeSteps;
febio_spec.Control.time_stepper.dtmin = dtmin;
febio_spec.Control.time_stepper.dtmax = dtmax;
febio_spec.Control.time_stepper.max_retries = max_retries;
febio_spec.Control.time_stepper.opt_iter = opt_iter;
febio_spec.Control.max_refs = max_refs;
febio_spec.Control.max_ups = max_ups;

%Output section
% -> log file
%Displacement
febio_spec.Output.logfile.ATTR.file = febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file = febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data = 'ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim = ',';
febio_spec.Output.logfile.node_data{1}.VAL = 1:size(tibiaMeshV,1);
%Stress
febio_spec.Output.logfile.element_data{1}.ATTR.file = febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data = 's1;s2;s3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim = ',';
febio_spec.Output.logfile.element_data{1}.VAL = 1:1:size(tibiaMeshE,1);
%Strain energy
febio_spec.Output.logfile.element_data{2}.ATTR.file = febioLogFileName_strainEnergy;
febio_spec.Output.logfile.element_data{2}.ATTR.data = 'sed';
febio_spec.Output.logfile.element_data{2}.ATTR.delim = ',';
febio_spec.Output.logfile.element_data{2}.VAL = 1:1:size(tibiaMeshE,1);

%Export to XML file
%Takes a while and likely needs JavaHeap memory increase
%   - Worked pretty well with 2,000Mb allocated to heap (see: https://au.mathworks.com/help/matlab/matlab_external/java-heap-memory-preferences.html)
%Or conversely reduce the number of mesh elements to help here
%This has now been done above to reduce file size (i.e. parts remeshed to
%3.0 edge length). It seems like this is necessary for the simulation to
%not take forever as well.
febioStruct2xml(febio_spec,febioFebFileName);

%% Run the FEBio analysis

%Analysis structure
febioAnalysis.run_filename = febioFebFileName; %The input file name
febioAnalysis.run_logname = febioLogFileName; %The name for the log file
febioAnalysis.disp_on = 1; %Display information on the command window
febioAnalysis.disp_log_on = 1; %Display convergence information in the command window
febioAnalysis.runMode = runMode;%'internal';
febioAnalysis.t_check = 10; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi = 1e99; %Max analysis time
febioAnalysis.maxLogCheckTime = 10; %Max log file checking time

%Run
[runFlag] = runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%%%%% FEBio seems to need to be in dynamic mode for the above to work.
%%%%% Runs in a reasonable amount of time when using point spacing of 4,
%%%%% seemed to take a long time with denser mesh.
%%%%%
%%%%% Couldn't get any steps to converge when fixing the plateau point with
%%%%% a rigid constraint --- but everything seemed to run OK when removing
%%%%% this rigid constraint.
%%%%%
%%%%% Force doesn't propogate through the bone though right now --- it is
%%%%% just concentrated at the top where the force is applied --- material
%%%%% too dense?
%%%%%
%%%%% If you apply to the ankle in the reverse direction it propogates a
%%%%% little more, but still doesn't look right --- perhaps muscle forces
%%%%% are needed for better overall stress???
%%%%%
%%%%% Nodal forces just seem too concentrated???

%% Import FEBio results
%Note this just comes from the example and could be done in a somewhat more
%efficient manner if desired

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
    V_DEF = N_disp_mat+repmat(tibiaMeshV,[1 1 size(N_disp_mat,3)]);
    
    %Import strain energy density from log file
    [~,E_energy,~] = importFEBio_logfile(fullfile(pwd,febioLogFileName_strainEnergy)); %Element strain energy

    %Remove nodal index column
    E_energy = E_energy(:,2:end,:);

    %Add initial state i.e. zero energy
    sizImport = size(E_energy);
    sizImport(3) = sizImport(3)+1;
    E_energy_mat_n = zeros(sizImport);
    E_energy_mat_n(:,:,2:end) = E_energy;
    E_energy = E_energy_mat_n;
    
    %Convert to appropriate data for visualisation
    [FE_face,C_energy_face] = element2patch(tibiaMeshE,E_energy(:,:,end),'tet4');
    [CV] = faceToVertexMeasure(FE_face,tibiaMeshV,C_energy_face);
    [indBoundary] = tesBoundary(FE_face,tibiaMeshV);
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

        [FE_face,C_energy_face] = element2patch(tibiaMeshE,E_energy(:,:,qt),'tet4');
        [CV] = faceToVertexMeasure(FE_face,tibiaMeshV,C_energy_face);

        %Set entries in animation structure
        animStruct.Handles{qt} = [hp1 hp1]; %Handles of objects to animate
        animStruct.Props{qt} = {'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt} = {V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
        
    end
    
    %Initiate animation
    anim8(hf,animStruct);

else
    error('FEBio simulation failed...no results to extract...')
end

%% Notes

%%%%% The above simulation seems to work, but it's clear that the set-up,
%%%%% materials etc. aren't great, given that the tibia bends...

