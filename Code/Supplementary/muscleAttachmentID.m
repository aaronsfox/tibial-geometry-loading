%% This script identifies the relevant muscle attachment point indices on the
%  tibia-fibula shape model by registering the mean shape model surface to
%  the OpenSim tibia model. These indices are used in subsequent FE
%  simulations by applying the muscle forces at these relevant points on
%  the experimental models.
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

%% Set-up
    
%Add supplementary code path (current folder)
addpath(genpath(pwd));

%Import OpenSim libraries
import org.opensim.modeling.*

%Set list of muscles to apply forces from (as per Haider et al. 2020)
%Quadriceps forces also included in place of patellar tendon in Haider et al.
muscleList = [{'semimem_r'};
    {'semiten_r'};
    {'bflh140_r'};
    {'bfsh140_r'};
    {'sart_r'};
    {'tfl_r'};
    {'grac_r'};
    {'soleus_r'};
    {'tibpost_r'};
    {'tibant_r'};
    {'fdl_r'};
    {'perbrev_r'};
    {'perlong_r'};
    {'edl_r'};
    {'recfem_r'};
    {'vasint_r'};
    {'vasmed_r'}];

%Set corresponding muscle force names (these are different for model labelling reasons)
muscleForceList = [{'semimem_r'};
    {'semiten_r'};
    {'bifemlh_r'};
    {'bifemsh_r'};
    {'sar_r'};
    {'tfl_r'};
    {'grac_r'};
    {'soleus_r'};
    {'tib_post_r'};
    {'tib_ant_r'};
    {'flex_dig_r'};
    {'per_brev_r'};
    {'per_long_r'};
    {'ext_dig_r'};
    {'rect_fem_r'};
    {'vas_int_r'};
    {'vas_med_r'}];

%Set colour order for muscles
muscleColours = [[0,31,63];
    [0,0,128];
    [0,116,217];
    [0,0,255];
    [127,219,255];
    [0,255,255];
    [57,204,204];
    [0,128,128];
    [61,153,112];
    [128,128,0];
    [46,204,64];
    [0,128,0];
    [1,255,112];
    [0,255,0];
    [255,220,0];
    [255,255,0];
    [255,65,54]] / 255;

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
optCPD.corresp = 1;         % estimate correspondence

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

%Set array to rotate views of surfaces later in reconstruction process
surfaceRot = [-90, 0, 90, 180];

%Set list to label reconstruction subplots with
viewLabel = [{'Anterior'}, {'Lateral'}, {'Posterior'}, {'Medial'}];

%Set options for tibial remeshing
optionStruct_tib.nb_pts = 3500; %Set desired number of points
optionStruct_tib.disp_on = 0; % Turn off command window text display
optionStruct_fib.nb_pts = 2000; %Set desired number of points
optionStruct_fib.disp_on = 0; % Turn off command window text display

%Load the tibia and trabecular shape models required
load('..\..\ShapeModels\tibia-fibula\tibiaFibulaShapeModel.mat');

%% Register generic OpenSim model and muscles to shape model

%Read in VTP as xml tree
[vtpXML, ~, ~] = xml_read('..\..\GenericOpenSimModel\Geometry\r_tibia.vtp');

%Extract points
openSimTibiaV = vtpXML.PolyData.Piece.Points.DataArray.CONTENT * 1000; %convert to mm
openSimTibiaF = vtpXML.PolyData.Piece.Polys.DataArray(1).CONTENT + 1; %add 1 for zero indexing

%Remesh to match shape model points
[openSimTibiaF, openSimTibiaV] = ggremesh(openSimTibiaF, openSimTibiaV, optionStruct_tib);

%Identify the tibia and fibula from the combined structure by grouping
%the vertices and faces
[~, groupIndexFaces] = groupVertices(tibiaFibulaShapeModel.F, tibiaFibulaShapeModel.meanPoints, 0);

%Identify which grouped section contains a higher volume
%This will be indicative of the tibia
if tetVolMeanEst(tibiaFibulaShapeModel.F(groupIndexFaces == 1,:), tibiaFibulaShapeModel.meanPoints) > ...
        tetVolMeanEst(tibiaFibulaShapeModel.F(groupIndexFaces == 2,:), tibiaFibulaShapeModel.meanPoints)
    %First index is tibia
    logicKeep = groupIndexFaces == 1;
else
    %Second index is tibia
    logicKeep = groupIndexFaces == 1;
end

%Separate the surfaces
[tibiaF, tibiaV] = patchCleanUnused(tibiaFibulaShapeModel.F(logicKeep,:), tibiaFibulaShapeModel.meanPoints);
[fibulaF, fibulaV] = patchCleanUnused(tibiaFibulaShapeModel.F(~logicKeep,:), tibiaFibulaShapeModel.meanPoints);

%Visualise and compare to shape model mean
cFigure; hold on;
gpatch(openSimTibiaF,openSimTibiaV, 'rw', 'k', 0.5);
gpatch(tibiaF, tibiaV, 'gw', 'k', 0.5);
axisGeom;

%Non-rigidly register opensim tibia to experimental tibia
[cpdTformOsim, cpdCorrespondOsim] = cpd_register(openSimTibiaV, tibiaV, optCPD);

%Sort the registered points so they align with the reference mesh
openSimTibiaV_shapeModel = openSimTibiaV(cpdCorrespondOsim,:);

%Identify muscle attachment points on mean shape model tibia

%Load generic model
genericModel = Model('..\..\GenericOpenSimModel\LaiArnold2017_refined.osim');
genericModelState = genericModel.initSystem();

%Get indices of muscle attachment sites against shape model

%Loop through muscles
for muscleInd = 1:length(muscleList)

    %Get the current muscles path point set
    muscPathPointSet = genericModel.updMuscles().get((muscleList{muscleInd})).getGeometryPath().getPathPointSet();

    %Identify the relevant attachment point to the tibia being origin or insertion

    %Get the body of the first and last attachment point
    originBodyName = char(muscPathPointSet.get(0).getBody().getName());
    insertionBodyName = char(muscPathPointSet.get(muscPathPointSet.getSize()-1).getBody().getName());

    %Check which is related to the tibia and extract location
    if strcmp(originBodyName, 'tibia_r')
        muscLoc = muscPathPointSet.get(0).getLocation(genericModelState);
    elseif strcmp(insertionBodyName, 'tibia_r')
        muscLoc = muscPathPointSet.get(muscPathPointSet.getSize()-1).getLocation(genericModelState);
    else
        error('Muscle doesn''t seem to attach to the tibia?')
    end

    %Get the 3D point in a useable form & convert to mm
    muscPt = [muscLoc.get(0), muscLoc.get(1), muscLoc.get(2)] * 1000;

    %Calculate distance between muscle point and new OpenSim tibia vertices
    muscPtDist = distancePoints3d(openSimTibiaV_shapeModel, muscPt);

    %Find and store indices in structure
    %First check for multiple matched distances
    minDistPts = find(muscPtDist == min(muscPtDist));
    muscleAttachmentInd.(muscleList{muscleInd}) = minDistPts(1);   

end

% % % %Visualise matching muscle points on experimental tibia
% % % %Use subplots to create different perspectives
% % % cFigure; hold on;    
% % % %Loop through four views to create subplot
% % % for viewNo = 1:4    
% % %     %Create subplot for current view
% % %     subplot(1,4,viewNo); hold on
% % %     %Add surface
% % %     hp = gpatch(tibiaF, tibiaV, 'kw', 'none', 0.7);
% % %     %Add muscle points
% % %     for muscleInd = 1:length(muscleList)
% % %         h1(muscleInd) = plotV(tibiaV(muscleAttachmentInd.(muscleList{muscleInd}),:), ...
% % %             'o', 'Color', muscleColours(muscleInd,:), 'MarkerFaceColor', muscleColours(muscleInd,:), 'MarkerSize', 5); 
% % %     end
% % %     
% % %     %Set axis view
% % %     axis equal; axis tight; view(0,90);
% % %     rotate(hp,[0 1 0], surfaceRot(viewNo));
% % %     rotate(h1,[0 1 0], surfaceRot(viewNo));
% % %     %Set axis parameters
% % %     camlight headlight; axis off
% % %     %Add colorbar on last view
% % %     if viewNo == 4
% % %         legend(h1,muscleList, 'Location', 'northeastoutside');
% % %     end
% % %     %Add title
% % %     title([viewLabel{viewNo},' View'], 'FontSize', 12);
% % % end

%% Write muscle indices to file

%Convert structure to table and write to file
writetable(struct2table(muscleAttachmentInd), ...
    '..\..\ShapeModels\meanMuscleAttachmentInds.csv');

%% ----- end of muscleAttachmentID.m ----- %%