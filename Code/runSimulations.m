function runSimulations(caseID)

    %% This script runs the finite element simulations for the surfaces segmented
    %  for the statistical shape model data. See README.MD at top of repo and
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
    %
    %  Inputs:
    %       caseID - number for case ID of the sample to use in simulations
    %
    % The primary unit set for FE simulations in FEBio was chosen as:
    %   [L]ength = mm
    %   [F]orce = N
    %   [T]ime = s
    % Therefore the set of secondary units can be derived as:
    %   [P]ressure/Stress = MPa
    %   [S]train = mm/mm
    
    % The case ID variable must correspond to one from the designated list 
    % included in the function below. The order of these matches to the
    % nodes within the shape model structure. It must be in integer format
    % otherwise an error will be through.

    %% Set-up
    
    %Turn off warnings
    warning off
    
    %Input checking
    
    %Check for case ID and integer value
    if nargin < 1
        error('A case ID number must be provided to analyse')
    else
        if ~isnumeric(caseID) || length(caseID) ~= 1
            error('Case ID must be a single number reflective of an actual case ID')
        end
    end
    %Check if case ID is in list
    %Provide list of cases for checking against. Note that this is ordered
    %on the data that is in the shape model
    availableCases = [147211; 102480; 102924; 103559; 103862; 107215; 107813;
        108421; 112802; 113033; 116070; 132433; 132592; 134434; 135418; 135881;
        136283; 140368; 140694; 141682; 144977; 145250; 146511; 149327; 170206;
        174615; 176087; 181140; 182401; 184171];
    %See if case can be found in list
    caseInd = find(caseID == availableCases);
    if isempty(caseInd)
        error('Case ID must be in the available cases list. Check the availableCases variable in function for a list of case IDs')
    end
    
    %Set the running trial labels
	runTrials = {'Run_3', 'Run_4', 'Run_5'};
    
    %Create directory to store results
    mkdir(['..\Simulations\case-',num2str(caseID)]);

    %Set home directory
    homeDir = pwd;

    %Add supplementary code path
    addpath(genpath([pwd,'\Supplementary']));
    
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
    load('..\ShapeModels\tibia-fibula\tibiaFibulaShapeModel.mat');
    load('..\ShapeModels\tibia-plus-trabecular\tibTrabShapeModel.mat');
    
    %Load and organise the landmark data
    meanLandmarks = readtable('..\ShapeModels\meanShapeModelLandmarks.csv',...
        'ReadRowNames', true);
    landmarkNames = meanLandmarks.Properties.RowNames;
    for landmark = 1:length(landmarkNames)
        landmarks.(landmarkNames{landmark}) = [meanLandmarks.X(landmark), ...
            meanLandmarks.Y(landmark), meanLandmarks.Z(landmark)];
    end

    %Convert landmarks coordinates to indices on shape model
    for landmark = 1:length(landmarkNames)
         landmarksInd.(landmarkNames{landmark}) = find(distancePoints3d(tibiaFibulaShapeModel.meanPoints,landmarks.(landmarkNames{landmark})) == ...
            min(distancePoints3d(tibiaFibulaShapeModel.meanPoints,landmarks.(landmarkNames{landmark}))));
    end

    %Load and organise the tib-fib interactions data
    meanInteractions = readtable('..\ShapeModels\meanShapeModelTibFibInteractions.csv',...
        'ReadRowNames', true);
    interactionNames = meanInteractions.Properties.RowNames;
    for interaction = 1:length(interactionNames)
        interactions.(interactionNames{interaction}) = [meanInteractions.X(interaction), ...
            meanInteractions.Y(interaction), meanInteractions.Z(interaction)];
    end

    %Convert interactions coordinates to indices on shape model
    for interaction = 1:length(interactionNames)
         interactionsInd.(interactionNames{interaction}) = find(distancePoints3d(tibiaFibulaShapeModel.meanPoints,interactions.(interactionNames{interaction})) == ...
            min(distancePoints3d(tibiaFibulaShapeModel.meanPoints,interactions.(interactionNames{interaction}))));
    end
    
    %Load muscle attachment index data
    muscleAttachmentInd = readtable('..\ShapeModels\meanMuscleAttachmentInds.csv',...
        'ReadRowNames', false);

    %Load HamnerDelp2013 generic running data
    genericRunningData = load('..\HamnerDelp2013-RunningData\Data\compiledHamnerDelp2013-RunningData.mat');

    %Set demographic data for current case
    demographicTable = readtable('..\ShapeModels\participantCharacteristics.csv');
    demographicInd = find(demographicTable.deidentified_record_number == caseID);
    massKG = demographicTable.Weight(demographicInd);
    heightM = demographicTable.Height(demographicInd) / 100;

    %Set the number of simulation time step points
    nSimPts = 1;
    
    %% Extract tibia, fibula and trabecular for current case
    
    %Extract and reshape the cases tibia-fibula from the shape model nodes
    tibFibV = reshape(tibiaFibulaShapeModel.nodes(caseInd,:), ...
        [3, length(tibiaFibulaShapeModel.mean)/3])';
    
    %Extract and reshape the cases tibia-trabecular from the shape model nodes
    tibTrabV = reshape(tibTrabShapeModel.nodes(caseInd,:), ...
        [3, length(tibTrabShapeModel.mean)/3])';
    
    %Rigidly register the tibial component of the tibia-trabecular model to
    %the tibial component of the tibial-fibula model
    [cpdRigTformTib] = cpd_register(tibFibV(1:optionStruct_tib.nb_pts,:), ...
        tibTrabV(1:optionStruct_tib.nb_pts,:), ...
        optCPD_rig);
    
    %Transform the trabecular component of the tibia-trabecular model using
    %the above rigid transformation. This will align the trabecular to the
    %tibia-fibula model. Extract the faces from the shape model while we're
    %doing this
    trabV = tibTrabV(optionStruct_tib.nb_pts+1:end,:) * cpdRigTformTib.R' + repmat(cpdRigTformTib.t,1,optionStruct_tib.nb_pts)';
    trabF = tibTrabShapeModel.F2;
    
    %Identify the tibia and fibula from the combined structure by grouping
    %the vertices and faces
    [~, groupIndexFaces] = groupVertices(tibiaFibulaShapeModel.F, tibFibV, 0);
    
% % %     %Visualise
% % %     cFigure; hold on
% % %     title('Grouped Faces')
% % %     gpatch(tibiaFibulaShapeModel.F, tibFibV, groupIndexFaces,'none');
% % %     axisGeom; camlight headlight;
% % %     colormap gjet; icolorbar;

    %Identify which grouped section contains a higher volume
    %This will be indicative of the tibia
    if tetVolMeanEst(tibiaFibulaShapeModel.F(groupIndexFaces == 1,:), tibFibV) > ...
            tetVolMeanEst(tibiaFibulaShapeModel.F(groupIndexFaces == 2,:), tibFibV)
        %First index is tibia
        logicKeep = groupIndexFaces == 1;
    else
        %Second index is tibia
        logicKeep = groupIndexFaces == 1;
    end

    %Separate the surfaces
    [tibiaF, tibiaV] = patchCleanUnused(tibiaFibulaShapeModel.F(logicKeep,:), tibFibV);
    [fibulaF, fibulaV] = patchCleanUnused(tibiaFibulaShapeModel.F(~logicKeep,:), tibFibV);

% % %     %Visualise all together
% % %     cFigure; hold on
% % %     gpatch(tibiaF, tibiaV, 'gw', 'none', 0.3)
% % %     gpatch(trabF, trabV, 'rw', 'none', 1.0)
% % %     gpatch(fibulaF, fibulaV, 'bw', 'none', 1.0)
% % %     axisGeom; camlight headlight
    
    %% Option to check OpenSim muscle attachments against experimental model
    
% %     %Visualise matching muscle points on experimental tibia
% %     %Use subplots to create different perspectives
% %     cFigure; hold on;
% %     %Loop through four views to create subplot
% %     for viewNo = 1:4
% %         %Create subplot for current view
% %         subplot(1,4,viewNo); hold on
% %         %Add surface
% %         hp = gpatch(tibiaF, tibiaV, 'kw', 'none', 0.7);
% %         %Add muscle points
% %         for muscleInd = 1:length(muscleList)
% %             h1(muscleInd) = plotV(tibiaV(muscleAttachmentIndExp.(muscleList{muscleInd}),:), ...
% %                 'o', 'Color', muscleColours(muscleInd,:), 'MarkerFaceColor', muscleColours(muscleInd,:), 'MarkerSize', 5); 
% %         end
% %         %Set axis view
% %         axis equal; axis tight; view(0,90);
% %         rotate(hp,[0 1 0], surfaceRot(viewNo));
% %         rotate(h1,[0 1 0], surfaceRot(viewNo));
% %         %Set axis parameters
% %         camlight headlight; axis off
% %         %Add colorbar on last view
% %         if viewNo == 4
% %             legend(h1,muscleList, 'Location', 'northeastoutside');
% %         end
% %         %Add title
% %         title([viewLabel{viewNo},' View'], 'FontSize', 12);
% %     end

    %% Identify landmarks on experimental surfaces
    
    %Create combined surface to look up landmarks
    tibiaFibulaV = [tibiaV; fibulaV];
    
    %Identify tibio-femoral articular surface for fixed boundary conditions
    %Use the simplified GIBOC algorithm within built in function to identify nodes
    [tibiaProxSurfacePts, tibiaDistSurfacePts, contactPointAJC] = extractTibioFemoralSurface(tibiaF, tibiaV);

    %% Work out the joint contact and muscle forces to apply in simulation
    
    %This section provides the force data for the JCFs and muscles for
    %simulating over the stance phase. The JCFs are applied to a single
    %point as per Haider et al.
    
    %Loop through run trials
    for runInd = 1:length(runTrials)
        
        %Calculate the resultant JRF for trial
        resultantJRF = sqrt((genericRunningData.meanDataNorm.(runTrials{runInd}).jrf.ankle_r_on_tibia_r_in_tibia_r_fx.^2 + ...
            genericRunningData.meanDataNorm.(runTrials{runInd}).jrf.ankle_r_on_tibia_r_in_tibia_r_fy.^2 + ...
            genericRunningData.meanDataNorm.(runTrials{runInd}).jrf.ankle_r_on_tibia_r_in_tibia_r_fz.^2));
        
        %Find peak resultant ankle joint reaction force as this is where we
        %will simulate the force data from
        trialSimInd = find(resultantJRF == max(resultantJRF), 1);
        
        %Extract the joint contact force at peak resultant JRF and normalise to body mass
        jointContactForcePeak = [genericRunningData.meanDataNorm.(runTrials{runInd}).jrf.ankle_r_on_tibia_r_in_tibia_r_fx(trialSimInd) * (massKG * 9.80665), ...
            genericRunningData.meanDataNorm.(runTrials{runInd}).jrf.ankle_r_on_tibia_r_in_tibia_r_fy(trialSimInd) * (massKG * 9.80665), ...
            genericRunningData.meanDataNorm.(runTrials{runInd}).jrf.ankle_r_on_tibia_r_in_tibia_r_fz(trialSimInd) * (massKG * 9.80665)];
        
        %Store the joint force data in a relevant structure
        simForces.(runTrials{runInd}).jointContactForcePeak = jointContactForcePeak;
        
        %Extract muscle forces
        
        %Loop through for each muscle
        for muscleInd = 1:length(muscleList)
            
            %Set variable names to extract data
            varX = [muscleForceList{muscleInd},'_X'];
            varY = [muscleForceList{muscleInd},'_Y'];
            varZ = [muscleForceList{muscleInd},'_Z'];
            
            %Extract the joint contact force at peak resultant JRF and normalise to body mass
            muscleForceContactPeak = [genericRunningData.meanDataNorm.(runTrials{runInd}).forcesVec.(char(varX))(trialSimInd) * (massKG * 9.80665), ...
                genericRunningData.meanDataNorm.(runTrials{runInd}).forcesVec.(char(varY))(trialSimInd) * (massKG * 9.80665), ...
                genericRunningData.meanDataNorm.(runTrials{runInd}).forcesVec.(char(varZ))(trialSimInd) * (massKG * 9.80665)];

            %Store the joint force data in a relevant structure
            simForces.(runTrials{runInd}).(muscleList{muscleInd}).muscleForceContactPeak = muscleForceContactPeak;
            
        end

    end

    %% Visualise final model, landmarks and muscle sites

% % %     %Visualise models and landmarks
% % %     cFigure; hold on;
% % % 
% % %     %Surfaces
% % %     gpatch(tibiaF, tibiaV, 'kw', 'none', 0.5)
% % %     gpatch(fibulaF, fibulaV, 'kw', 'none')
% % %     gpatch(trabF, trabV, 'bw', 'none', 0.5);
% % %     axis equal; axis tight; view(90,0);
% % %     camlight headlight
% % % 
% % %     %Landmarks
% % %     for landmark = 1:length(landmarkNames)
% % %         plotV(tibiaFibulaV(landmarksInd.(landmarkNames{landmark}),:),'b.','MarkerSize',20);
% % %     end
% % % 
% % %     %Interaction connections
% % %     plotV([tibiaFibulaV(interactionsInd.AntProxTib,:); tibiaFibulaV(interactionsInd.AntProxFib,:)], ...
% % %         '-o', 'Color', 'g', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'LineWidth', 2);
% % %     plotV([tibiaFibulaV(interactionsInd.PostProxTib,:);tibiaFibulaV(interactionsInd.PostProxFib,:)], ...
% % %         '-o', 'Color', 'g', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'LineWidth', 2);
% % %     plotV([tibiaFibulaV(interactionsInd.AntDistTib,:);tibiaFibulaV(interactionsInd.AntDistFib,:)], ...
% % %         '-o', 'Color', 'g', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'LineWidth', 2);
% % %     plotV([tibiaFibulaV(interactionsInd.PostDistTib,:);tibiaFibulaV(interactionsInd.PostDistFib,:)], ...
% % %         '-o', 'Color', 'g', 'MarkerSize', 5, 'MarkerFaceColor', 'g', 'LineWidth', 2);
% % % 
% % %     %Muscle attachments
% % %     for muscleInd = 1:length(muscleList)
% % %         plotV(tibiaFibulaV(muscleAttachmentInd.(muscleList{muscleInd}),:), ...
% % %             'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5);   
% % %     end
% % % 
% % %     %Ankle joint centre contact point
% % %     plotV(contactPointAJC, 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
% % % 
% % %     %Proximal and distal surface estimates
% % %     plotV(tibiaProxSurfacePts, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 2.5);
% % %     plotV(tibiaDistSurfacePts, 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 2.5);

    %% Mesh the surfaces using tetgen
    %  IMPORTANT NOTE: for multi-region meshing to work, the order in joining
    %  the element sets and regions must be from outside-inside (i.e. cortical
    %  then trabecular). It won't work the other way around

    %Use function to generate meshes
    [meshOutputFullTibia, meshOutputFibula] = generateVolumetricMeshes(tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV);
    
    %Save mesh outputs to file for later when analysing simulations
    save(['..\Simulations\case-',num2str(caseID),'\tetGenMeshOutputs.mat'], ...
        'meshOutputFullTibia', 'meshOutputFibula');

    %% Set-up & run FEBio simulation
    
    %Set FEA control settings
    settingsFEA.numTimeSteps = nSimPts; %Number of time steps desired
    settingsFEA.max_refs = 25; %Max reforms
    settingsFEA.max_ups = 0; %Set to zero to use full-Newton iterations
    settingsFEA.opt_iter = 6; %Optimum number of iterations
    settingsFEA.max_retries = 5; %Maximum number of retires
    settingsFEA.dtmin = (1/settingsFEA.numTimeSteps)/100; %Minimum time step size
    settingsFEA.dtmax = 1/settingsFEA.numTimeSteps; %Maximum time step size
    settingsFEA.runMode = 'internal'; %'external' or 'internal'

    %Set material parameters for FEA
    %Parameters currently from Edwards et al. (2010)
    matParameters.corticalYoungs = 18600;
    matParameters.trabYoungs = 10400;
    
    %Package up landmarks data to a structure to feed into function
    simLandmarks.contactPointAJC = contactPointAJC;
    simLandmarks.tibiaDistSurfacePts = tibiaDistSurfacePts;
    simLandmarks.tibiaProxSurfacePts = tibiaProxSurfacePts;
    simLandmarks.tibiaFibulaV = tibiaFibulaV;
    simLandmarks.muscleAttachmentInd = muscleAttachmentInd;
    simLandmarks.interactionsInd = interactionsInd;
    simLandmarks.landmarksInd = landmarksInd;
    
    %Package up original face and vertices data
    originalSurfaces.tibiaF = tibiaF;
    originalSurfaces.tibiaV = tibiaV;
    originalSurfaces.fibulaF = fibulaF;
    originalSurfaces.fibulaV = fibulaV;
    originalSurfaces.trabF = trabF;
    originalSurfaces.trabV = trabV;
    
    %Loop through run trials
    for runInd = 1:length(runTrials)
        
        %Set run trial variable
        runTrial = runTrials{runInd};
    
        %Use convenience function to create set-up file
        [febioFebFileName, febioLogFileName] = setUpFEBioSim(runTrial, caseID, ...
            meshOutputFullTibia, meshOutputFibula, originalSurfaces, ...
            simLandmarks, simForces, ...
            settingsFEA, matParameters);

        %Run FEBio simulation

        %Settings for the run
        febioAnalysis.run_filename = febioFebFileName; %The input file name
        febioAnalysis.run_logname = febioLogFileName; %The name for the log file
        febioAnalysis.disp_on = 1; %Display information on the command window
        febioAnalysis.runMode = 'internal';

        %Navigate to simulation directory for ease of running
        cd(fileparts(febioFebFileName));

        %Run the simulation
        [runFlag] = runMonitorFEBio(febioAnalysis);

        %Return to home directory
        cd(homeDir);
        
    end

end
    
%% ----- end of runSimulations.m ----- %%