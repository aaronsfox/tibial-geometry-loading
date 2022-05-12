%% This script runs the finite element simulations for the surfaces generated
%  from the statistical shape model data. See README.MD at top of repo and
%  comments throughout this script for details.
%
% Authors:
%     Meghan Keast
%     Centre for Sport Research
%     Deakin University
%     mfkeast@deakin.edu.au
%
%     Aaron Fox
%     Centre for Sport Research
%     Deakin University
%     aaron.f@deakin.edu.au

%%%%% TODO: this just uses mean model --- adapt to use variants of mean...
%%%%% e.g. mean points etc. need to be adaptable to a new model ---
%%%%% therefore need to identify the index matching/closest to that point
%%%%% to identify it in an altered surface...

%% Set-up

%Set home directory
homeDir = pwd;

%Add supplementary code path
addpath(genpath([pwd,'\Supplementary']));

%Set options for CPD algorithm - non rigid
optCPD.method = 'nonrigid'; % use nonrigid registration
optCPD.beta = 2;            % the width of Gaussian kernel (smoothness)
optCPD.lambda = 3;          % regularization weight
optCPD.viz = 0;             % don't visualise
optCPD.outliers = 0;        % don't account for outliers
optCPD.fgt = 0;             % do not use FGT (default)
optCPD.normalize = 1;       % normalize to unit variance and zero mean before registering (default)
optCPD.corresp = 0;         % compute correspondence vector at the end of registration (not being estimated by default)
optCPD.max_it = 100;        % max number of iterations
optCPD.tol = 1e-4;          % tolerance
optCPD.corresp = 0;         % estimate correspondence

%Set options for CPD algorithm - rigid
optCPD_rig.method = 'rigid';    % use rigid registration
optCPD_rig.viz = 0;             % don't visualise
optCPD_rig.outliers = 0;        % don't account for outliers
optCPD_rig.fgt = 0;             % do not use FGT (default)
optCPD_rig.normalize = 1;       % normalize to unit variance and zero mean before registering (default)
optCPD_rig.corresp = 0;         % compute correspondence vector at the end of registration (not being estimated by default)
optCPD_rig.max_it = 100;        % max number of iterations
optCPD_rig.tol = 1e-4;          % tolerance
optCPD_rig.corresp = 0;         % estimate correspondence

%Load the tibia and trabecular shape models required
load('..\ShapeModels\tibia\tibiaShapeModel.mat');
load('..\ShapeModels\trabecular-tibia\trabShapeModel.mat');
load('..\ShapeModels\tibia-fibula\tibiaFibulaShapeModel.mat');

%Load and organise the landmark data
meanLandmarks = readtable('..\ShapeModels\meanShapeModelLandmarks.csv',...
    'ReadRowNames', true);
landmarkNames = meanLandmarks.Properties.RowNames;
for landmark = 1:length(landmarkNames)
    landmarks.(landmarkNames{landmark}) = [meanLandmarks.X(landmark), ...
        meanLandmarks.Y(landmark), meanLandmarks.Z(landmark)];
end

%Load and organise the tib-fib interactions data
meanInteractions = readtable('..\ShapeModels\meanShapeModelTibFibInteractions.csv',...
    'ReadRowNames', true);
interactionNames = meanInteractions.Properties.RowNames;
for interaction = 1:length(interactionNames)
    interactions.(interactionNames{interaction}) = [meanInteractions.X(interaction), ...
        meanInteractions.Y(interaction), meanInteractions.Z(interaction)];
end

%% Separate out tibia and fibula from shape model
%  This is necessary to apply different properties to these surfaces within
%  the FE simulations

%Group the vertices and faces from the combined shape model
[groupIndexVertices,groupIndexFaces] = groupVertices(tibiaFibulaShapeModel.F, ...
    tibiaFibulaShapeModel.meanPoints,0);

% % % %Visualise
% % % cFigure; hold on
% % % title('Grouped Faces')
% % % gpatch(tibiaFibulaShapeModel.F, tibiaFibulaShapeModel.meanPoints, groupIndexFaces,'none');
% % % axisGeom; camlight headlight;
% % % colormap gjet; icolorbar;

%Identify which grouped section contains a higher volume
%This will be indicative of the tibia
if tetVolMeanEst(tibiaFibulaShapeModel.F(groupIndexFaces == 1,:),tibiaFibulaShapeModel.meanPoints) > ...
        tetVolMeanEst(tibiaFibulaShapeModel.F(groupIndexFaces == 2,:),tibiaFibulaShapeModel.meanPoints)
    %First index is tibia
    logicKeep = groupIndexFaces == 1;
else
    %Second index is tibia
    logicKeep = groupIndexFaces == 1;
end
%Separate the surfaces
[tibiaF, tibiaV] = patchCleanUnused(tibiaFibulaShapeModel.F(logicKeep,:), ...
    tibiaFibulaShapeModel.meanPoints);
[fibulaF, fibulaV] = patchCleanUnused(tibiaFibulaShapeModel.F(~logicKeep,:), ...
    tibiaFibulaShapeModel.meanPoints);

% % % %Visualise
% % % cFigure; hold on
% % % gpatch(tibiaF,tibiaV,'gw')
% % % gpatch(fibulaF,fibulaV,'bw')
% % % axisGeom;

% % % %Additionally, extract the trabecular mdeol
% % % trabV = trabShapeModel.meanPoints;
% % % trabF = trabShapeModel.F;

%% Generate trabeculer for the tibia surface based on developed model

%Create wait bar
wbar = waitbar(0,'Predicting trabecular surface. This can take a moment. Please be patient with me...');

%Create trabecular surface
[trabF, trabV] = generateTrabecularModel(tibiaShapeModel, trabShapeModel, tibiaF, tibiaV);

%Update waitbar
waitbar(1, wbar, 'Trabecular surface predicted with no issues (hopefully).');
pause(1); close(wbar);

%% Visualise models and landmarks
       
% % % %Visualise models and landmarks
% % % cFigure; hold on;
% % % %Surfaces
% % % gpatch(tibiaF, tibiaV, 'gw', 'none', 0.3)
% % % gpatch(fibulaF, fibulaV, 'bw', 'none')
% % % gpatch(trabF, trabV, 'rw', 'none', 1);
% % % %Landmarks
% % % for landmark = 1:length(landmarkNames)
% % %     plotV(landmarks.(landmarkNames{landmark}),'k.','MarkerSize',20);
% % % end
% % % %Interaction connections
% % % plotV([interactions.AntProxTib;interactions.AntProxFib], ...
% % %     '-o', 'Color', 'r', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'LineWidth', 2);
% % % plotV([interactions.PostProxTib;interactions.PostProxFib], ...
% % %     '-o', 'Color', 'r', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'LineWidth', 2);
% % % plotV([interactions.AntDistTib;interactions.AntDistFib], ...
% % %     '-o', 'Color', 'r', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'LineWidth', 2);
% % % plotV([interactions.PostDistTib;interactions.PostDistFib], ...
% % %     '-o', 'Color', 'r', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'LineWidth', 2);
% % % %Axes properties
% % % axisGeom; camlight headlight

%% Identify ankle joint centre for contact force application

%Take mean of medial and lateral malleoli landmarks
landmarks.AJC = mean([landmarks.MM; landmarks.LM]);

%Project point up longitudinal axis to find appropriate point on surface
%Create line for intersection
intersectLineAJC = createLine3d(landmarks.AJC, 0, 1, 0); %project along Y-axis

%Idnentify tibia triangles that are intersected by the lines
for triangleNo = 1:length(tibiaF)
    %Get the current triangle
    intersectTriangle = [tibiaV(tibiaF(triangleNo,1),:); ...
        tibiaV(tibiaF(triangleNo,2),:); ...
        tibiaV(tibiaF(triangleNo,3),:)];
    %Check intersection
    checkIntersect(triangleNo, :) = intersectLineTriangle3d(intersectLineAJC, intersectTriangle);
end

%Remove NaN's from intersection checks
intersectPointsAJC = rmmissing(checkIntersect);

%There will likely be a couple of intersections at the bottom and top of tibia
%Find the closes to the ankle joint centre to identify the appropriate one
intersectPointsIndAJC = find(distancePoints3d(intersectPointsAJC, landmarks.AJC) == min(distancePoints3d(intersectPointsAJC, landmarks.AJC)));

%Find the point on the tibia surface that is closest to the minimum intersection
contactPointAJC = tibiaV( ...
    find(distancePoints3d(intersectPointsAJC(intersectPointsIndAJC,:), tibiaV) == ...
    min(distancePoints3d(intersectPointsAJC(intersectPointsIndAJC,:), tibiaV))), :);

%%



%% Mesh the surfaces using tetgen
%  IMPORTANT NOTE: for multi-region meshing to work, the order in joining
%  the element sets and regions must be from outside-inside (i.e. cortical
%  then trabecular). It won't work the other way around

%%%% TODO: will need to ensure that the trabecular mapping process works
%%%% for altered tibia surfaces --- works on the basic mean but there may
%%%% be some intersecting faces when these get more dramatically tweaked...

%Join the tibia and trabecular surfaces
[fullTibiaF, fullTibiaV, fullTibiaC] = joinElementSets({tibiaF trabF},{tibiaV trabV});

%Mesh the two surfaces
%Set properties
[V_region1] = getInnerPoint(fullTibiaF(fullTibiaC==2,:),fullTibiaV); %First interior point
[V_region2] = getInnerPoint({tibiaF trabF},{tibiaV trabV}); %Second interior point
V_regions = [V_region1; V_region2]; %Collect region points
V_holes = []; %Define hole points
[regionTetVolume1] = tetVolMeanEst(tibiaF, tibiaV); %Volume estimate for regular tets
[regionTetVolume2] = tetVolMeanEst(trabF,trabV); %Volume estimate for regular tets
regionTetVolumes = [regionTetVolume1 regionTetVolume2];
stringOpt = '-pq1.2AaY'; %Tetgen options
%Create tetgen input structure
tetGenStructFullTibia.stringOpt = stringOpt; %Tetgen options
tetGenStructFullTibia.Faces = fullTibiaF; %Boundary faces
tetGenStructFullTibia.Nodes = fullTibiaV; %Nodes of boundary
tetGenStructFullTibia.faceBoundaryMarker = fullTibiaC;
tetGenStructFullTibia.regionPoints = V_regions; %Interior points for regions
tetGenStructFullTibia.holePoints = V_holes; %Interior points for holes
tetGenStructFullTibia.regionA = regionTetVolumes; %Desired tetrahedral volume for each region
%Mesh model using tetgen
[meshOutputFullTibia] = runTetGen(tetGenStructFullTibia);

%Access mesh output structure
fullTibiaE = meshOutputFullTibia.elements; %The elements
fullTibiaV = meshOutputFullTibia.nodes; %The vertices or nodes
fullTibiaCE = meshOutputFullTibia.elementMaterialID; %Element material or region id
fullTibiaFb = meshOutputFullTibia.facesBoundary; %The boundary faces
fullTibiaCb = meshOutputFullTibia.boundaryMarker; %The boundary markers

% % % %Visualise
% % % meshView(meshOutputFullTibia);

%Mesh the fibula
innerPointFibula = getInnerPoint(fibulaF, fibulaV);

%Set mesh parameters
tetVolumeFibula = tetVolMeanEst(fibulaF, fibulaV);
tetGenStructFibula.stringOpt = '-pq1.2AaY';
tetGenStructFibula.Faces = fibulaF;
tetGenStructFibula.Nodes = fibulaV;
tetGenStructFibula.holePoints = [];
tetGenStructFibula.regionPoints = innerPointFibula;
tetGenStructFibula.regionA = tetVolumeFibula;

%Run tetgen
[meshOutputFibula] = runTetGen(tetGenStructFibula);

%Access elements, nodes, and boundary faces
fibulaMeshE = meshOutputFibula.elements;
fibulaMeshV = meshOutputFibula.nodes;
fibulaMeshFb = meshOutputFibula.facesBoundary;
fibulaMeshCb = meshOutputFibula.boundaryMarker;
fibulaMeshCE = meshOutputFibula.elementMaterialID;

% % % %Visualise solid mesh
% % % hFig = cFigure; hold on;
% % % optionStruct.hFig = hFig;
% % % meshView(meshOutputFibula,optionStruct);
% % % axisGeom;





%% NOTES:

%%%%% Use a generic OpenSim model to register the models tibia to the
%%%%% simulation and map the muscle attachment points to the simulation
%%%%% tibia...

%% OpenSim tibia model...

%Read in VTP as xml tree
[vtpXML, ~, ~] = xml_read('..\GenericOpenSimModel\Geometry\r_tibia.vtp');

%Extract points
openSimTibiaV = vtpXML.PolyData.Piece.Points.DataArray.CONTENT * 1000; %convert to mm
openSimTibiaF = vtpXML.PolyData.Piece.Polys.DataArray(1).CONTENT + 1; %add 1 for zero indexing
% % % cFigure; hold on;
% % % % % % plotV(openSimTibiaV);
% % % gpatch(openSimTibiaF,openSimTibiaV);
% % % axisGeom;

%Remesh to match shape model points
%Set options for remeshing
optionStruct_tib.nb_pts = 3500; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
%Remesh
[openSimTibiaF,openSimTibiaV] = ggremesh(openSimTibiaF, openSimTibiaV, optionStruct_tib);

% % % %Apply procrustes to rigidly transform and align OpenSim to shape model
% % % [~, openSimTibiaV_procReg] = procrustes(tibiaShapeModel.meanPoints, ...
% % %     openSimTibiaV, 'scaling', false);

%Visualise and compare to shape model mean
cFigure; hold on;
% % % plotV(openSimTibiaV);
gpatch(openSimTibiaF,openSimTibiaV, 'rw', 'k', 0.5);
gpatch(tibiaShapeModel.F, tibiaShapeModel.meanPoints, 'gw', 'k', 0.5);
axisGeom;


%%%% Seems to work, need to figure out best way to align and identify
%%%% corresponding points so that muscle attachment sites can be translated
%%%% onto shape model surface...




%%