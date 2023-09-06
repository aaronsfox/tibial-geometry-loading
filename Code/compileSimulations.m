function compileSimulations()

    %% This script compiles the results from the FE simulations run from the
    %  easyRun/runSimulations function. It brings the data from each case
    %  together to create a database for each case that will be used in
    %  subsequent analyses. See the README.MD at top of repo and comments
    %  throughout this script for details.
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

    %Turn off warnings
    warning off
    
    %Input checking
    if nargin > 0
        error('No inputs are required for this function.');
    end
    
    %Add supplementary code path
    addpath(genpath([pwd,'\Supplementary']));

    %Load the shape models required
    load('..\ShapeModels\tibia-fibula\tibiaFibulaShapeModel.mat');

    %Set options for tibia and fibula data
    %Note that these are higher than the original due to the finer meshes
    optionStruct_tib.nb_pts = 7000; %Set desired number of points
    optionStruct_tib.disp_on = 0; % Turn off command window text display
    optionStruct_fib.nb_pts = 3000; %Set desired number of points
    optionStruct_fib.disp_on = 0; % Turn off command window text display
    
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
    
    %Provide list of cases for checking against. Note that this is ordered
    %on the data that is in the shape model
    simulateCases = [147211; 102480; 102924; 103559; 103862; 107215; 107813;
        108421; 112802; 113033; 116070; 132433; 132592; 134434; 135418; 135881;
        136283; 140368; 140694; 141682; 144977; 145250; 146511; 149327; 170206;
        174615; 176087; 181140; 182401; 184171];

    %Set the running trial labels
	runTrials = {'Run_3', 'Run_4', 'Run_5'};
    
    %Create a variable for run names
    runNames = {'3 m/s', '4 m/s', '5 m/s'};
    
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
    
    %Read in muscle attachment indices
    muscleAttachmentInd = readtable('..\ShapeModels\meanMuscleAttachmentInds.csv',...
        'ReadRowNames', false);

    %Set array to rotate views of surfaces later in reconstruction process
    surfaceRot = [-90, 0, 90, 180];

    %Set list to label reconstruction subplots with
    viewLabel = [{'Anterior'}, {'Lateral'}, {'Posterior'}, {'Medial'}];
    
    %Create a variable of simulation data labels to extract and process
    simVars = [{'E1'}, {'E3'}, {'effectiveStress'}, {'effectiveStrain'}];
    
    %Extract the mean isolated tibia surface from the shape model for
    %registration against the individual cases. Note the use of 3500 points
    %to begin with here for tibia and 2000 for fibula.
    shapeModelTibiaLogicVertices = logical([ones(3500,1); zeros(2000,1)]);
    shapeModelTibiaLogicFaces = all(shapeModelTibiaLogicVertices(tibiaFibulaShapeModel.F),2);
    tibiaFibulaShapeModel.tibiaF = tibiaFibulaShapeModel.F(shapeModelTibiaLogicFaces,:);
    [tibiaFibulaShapeModel.tibiaF,tibiaFibulaShapeModel.tibiaV] = ...
        patchCleanUnused(tibiaFibulaShapeModel.tibiaF,tibiaFibulaShapeModel.meanPoints);
    %Remesh the tibia to the finer number of points used in simulations
    [tibiaFibulaShapeModel.tibiaF,tibiaFibulaShapeModel.tibiaV] = ...
        ggremesh(tibiaFibulaShapeModel.tibiaF,tibiaFibulaShapeModel.tibiaV, optionStruct_tib);
    [tibiaFibulaShapeModel.tibiaF,tibiaFibulaShapeModel.tibiaV] = ...
        mergeVertices(tibiaFibulaShapeModel.tibiaF,tibiaFibulaShapeModel.tibiaV);
    
    %Define the generic shaft proportion length
    %Bruce et al. (2022) defines the tibia shaft at 20%-80% of the tibia
    %length. We use a similar definition here, however the 15%-75% length
    %appears to fit better for the area of interest.
    lowerShaftProp = 0.15; upperShaftProp = 0.75;
    
% % %     %Calculate the nodal boundaries for the shaft length based on
% % %     %the generic proportions (along the y-axis as this is the
% % %     %longitudinal axis).
% % %     tibiaFibulaShapeModel.minShaftVal = ...
% % %         min(tibiaFibulaShapeModel.tibiaV(:,2)) + ...
% % %         ((max(tibiaFibulaShapeModel.tibiaV(:,2)) - ...
% % %         min(tibiaFibulaShapeModel.tibiaV(:,2))) * lowerShaftProp);
% % %     tibiaFibulaShapeModel.maxShaftVal = ...
% % %         min(tibiaFibulaShapeModel.tibiaV(:,2)) + ...
% % %         ((max(tibiaFibulaShapeModel.tibiaV(:,2)) - ...
% % %         min(tibiaFibulaShapeModel.tibiaV(:,2))) * upperShaftProp);
    
    %% Extract simulation data
    
    %Create waitbar
    caseWaitBar = waitbar(0, 'Loading simulation data...');
    
    %Loop through cases
    for caseInd = 1:length(simulateCases)
        
        %Update waitbar
        waitbar(caseInd/(length(simulateCases)+1), caseWaitBar, ...
            ['Loading data for case ',num2str(simulateCases(caseInd)),' ...']);

        %% Extract general case data

        %Create a variable for the caseID
        caseID = ['case-',num2str(simulateCases(caseInd))];
        caseID_label = ['case',num2str(simulateCases(caseInd))];

        %Load the mesh outputs for the case
        meshData.(char(caseID_label)) = load(['..\Simulations\',caseID,'\tetGenMeshOutputs.mat']);

        %Access mesh data from tibia output for cortical elements (element ID = -3)
        meshData.(char(caseID_label)).corticalTibiaE = ...
            meshData.(char(caseID_label)).meshOutputFullTibia.elements(meshData.(char(caseID_label)).meshOutputFullTibia.elementMaterialID == -3,:);

        %Extract out the mesh data for easier use throughout script
        meshData.(char(caseID_label)).E = meshData.(char(caseID_label)).meshOutputFullTibia.elements;
        meshData.(char(caseID_label)).V = meshData.(char(caseID_label)).meshOutputFullTibia.nodes;
        meshData.(char(caseID_label)).Fb = meshData.(char(caseID_label)).meshOutputFullTibia.facesBoundary;
        meshData.(char(caseID_label)).Cb = meshData.(char(caseID_label)).meshOutputFullTibia.boundaryMarker; %Cb == 1 is the outer cortical part; Cb == 2 is the inner cortical part
        meshData.(char(caseID_label)).CE = meshData.(char(caseID_label)).meshOutputFullTibia.elementMaterialID; %CE == -3 is the cortical part

        %Select out cortical elements
        meshData.(char(caseID_label)).corticalE = meshData.(char(caseID_label)).E(meshData.(char(caseID_label)).CE == -3,:);

        %Extract tibia surface faces and vertices for later use
        [meshData.(char(caseID_label)).tibiaF, meshData.(char(caseID_label)).tibiaV] = ...
            patchCleanUnused(meshData.(char(caseID_label)).Fb(meshData.(char(caseID_label)).Cb == 1,:), meshData.(char(caseID_label)).V);

        %% Extract isolated outer cortical surface as this will be used to
        %  map simulation data to for analyses

        %Extract faces for cortical parts of tibia
        [meshData.(char(caseID_label)).corticalF] = element2patch(meshData.(char(caseID_label)).corticalE);

        %Clean-up unused vertices with respect to cortical faces
        [meshData.(char(caseID_label)).corticalF, meshData.(char(caseID_label)).corticalV] = ...
            patchCleanUnused(meshData.(char(caseID_label)).corticalF, meshData.(char(caseID_label)).V);

        %Identify where cortical surface links to tibia outer surface
        corticalKNN = knnsearch(meshData.(char(caseID_label)).corticalV, ...
            meshData.(char(caseID_label)).tibiaV);

        %Set logical for cortical vertices of outer surface to keep
        checkInd = linspace(1, length(meshData.(char(caseID_label)).corticalV), length(meshData.(char(caseID_label)).corticalV))';
        meshData.(char(caseID_label)).outerCorticalLogicVertices = ismember(checkInd, corticalKNN, 'rows');

        %Register the cases outer cortical surface to align with shape
        %model mean surface. We'll work with these vertices registered to
        %the shape model from here on.
        [~, meshData.(char(caseID_label)).ptCorrespondence] = ...
            cpd_register(meshData.(char(caseID_label)).corticalV(meshData.(char(caseID_label)).outerCorticalLogicVertices,:), ...
            tibiaFibulaShapeModel.tibiaV, optCPD);

        %Convert the points and faces to the registered shape model surface
        meshData.(char(caseID_label)).outerCorticalV = meshData.(char(caseID_label)).corticalV(meshData.(char(caseID_label)).outerCorticalLogicVertices,:);
        meshData.(char(caseID_label)).outerCorticalV = meshData.(char(caseID_label)).outerCorticalV(meshData.(char(caseID_label)).ptCorrespondence,:);
        meshData.(char(caseID_label)).outerCorticalF = tibiaFibulaShapeModel.tibiaF;

% % %         %Visualise cases outer tibial surface
% % %         cFigure;
% % %         gpatch(meshData.(char(caseID_label)).outerCorticalF, meshData.(char(caseID_label)).outerCorticalV);
% % %         axisGeom;

        %% Extract the shaft of the outer cortical tibia for subsequent analyses

        %Calculate the nodal boundaries for the shaft length based on
        %the generic proportions (along the y-axis as this is the
        %longitudinal axis).
        meshData.(char(caseID_label)).minShaftVal = ...
            min(meshData.(char(caseID_label)).tibiaV(:,2)) + ...
            ((max(meshData.(char(caseID_label)).tibiaV(:,2)) - ...
            min(meshData.(char(caseID_label)).tibiaV(:,2))) * lowerShaftProp);
        meshData.(char(caseID_label)).maxShaftVal = ...
            min(meshData.(char(caseID_label)).tibiaV(:,2)) + ...
            ((max(meshData.(char(caseID_label)).tibiaV(:,2)) - ...
            min(meshData.(char(caseID_label)).tibiaV(:,2))) * upperShaftProp);

        %Extract shaft of tibial surface
        meshData.(char(caseID_label)) = extractTibialShaftSurface(meshData.(char(caseID_label)));

% % %         %Visualise final shaft in the context of the entire tibia
% % %         cFigure; hold on
% % %         gpatch(meshData.(char(caseID_label)).outerCorticalF, meshData.(char(caseID_label)).outerCorticalV, 'rw', 'none', 0.3);
% % %         gpatch(meshData.(char(caseID_label)).shaftF, meshData.(char(caseID_label)).shaftV, 'gw', 'k');
% % %         axisGeom; camlight headlight

        %% Loop through run trials
        for runInd = 1:length(runTrials)

            %% Load simulation data
            %For each variable the data is loaded, the final step taken,
            %and the indexing column removed.

            %Displacement
            [~, simData.(char(caseID_label)).(runTrials{runInd}).displacement] = ...
                importFEBio_logfile(['..\Simulations\',caseID,'\',runTrials{runInd},'\',caseID,'_',runTrials{runInd},'_disp_out.txt']);
            simData.(char(caseID_label)).(runTrials{runInd}).displacement = simData.(char(caseID_label)).(runTrials{runInd}).displacement(:,2:end,end);

            %Strain
            [~, simData.(char(caseID_label)).(runTrials{runInd}).strain] = ...
                importFEBio_logfile(['..\Simulations\',caseID,'\',runTrials{runInd},'\',caseID,'_',runTrials{runInd},'_strain_out.txt']);
            simData.(char(caseID_label)).(runTrials{runInd}).strain = simData.(char(caseID_label)).(runTrials{runInd}).strain(:,2:end,end);

            %Effective stress
            [~, simData.(char(caseID_label)).(runTrials{runInd}).effectiveStress] = ...
                importFEBio_logfile(['..\Simulations\',caseID,'\',runTrials{runInd},'\',caseID,'_',runTrials{runInd},'_stress_out.txt']);
            simData.(char(caseID_label)).(runTrials{runInd}).effectiveStress = simData.(char(caseID_label)).(runTrials{runInd}).effectiveStress(:,2:end,end); 

            %Effective strain
            [~, simData.(char(caseID_label)).(runTrials{runInd}).effectiveStrain] = ...
                importFEBio_logfile(['..\Simulations\',caseID,'\',runTrials{runInd},'\',caseID,'_',runTrials{runInd},'_effectiveStrain_out.txt']);
            simData.(char(caseID_label)).(runTrials{runInd}).effectiveStrain = simData.(char(caseID_label)).(runTrials{runInd}).effectiveStrain(:,2:end,end);

            %% Clean up simulation data to focus points for analysis

            %The first step is to map the cortical elements and simulation
            %data to the outer surface of the tibia, as this is where we
            %will focus our 3D surface analysis.

            %Extract the colour data out for variables for entire cortical surface
            [~, simData.(char(caseID_label)).(runTrials{runInd}).corticalC.E1] = element2patch(meshData.(char(caseID_label)).corticalE, ...
                simData.(char(caseID_label)).(runTrials{runInd}).strain(:,1));
            [~, simData.(char(caseID_label)).(runTrials{runInd}).corticalC.E3] = element2patch(meshData.(char(caseID_label)).corticalE, ...
                simData.(char(caseID_label)).(runTrials{runInd}).strain(:,3));
            [~, simData.(char(caseID_label)).(runTrials{runInd}).corticalC.effectiveStress] = element2patch(meshData.(char(caseID_label)).corticalE, ...
                simData.(char(caseID_label)).(runTrials{runInd}).effectiveStress);
            [~, simData.(char(caseID_label)).(runTrials{runInd}).corticalC.effectiveStrain] = element2patch(meshData.(char(caseID_label)).corticalE, ...
                simData.(char(caseID_label)).(runTrials{runInd}).effectiveStrain);
            %Special case to invert principal compression for ease of visualisation
            simData.(char(caseID_label)).(runTrials{runInd}).corticalC.E3 = simData.(char(caseID_label)).(runTrials{runInd}).corticalC.E3 * -1;

            %Convert face colouring data to vertices on cortical surface
            for varInd = 1:length(simVars)
                simData.(char(caseID_label)).(runTrials{runInd}).corticalCV.(simVars{varInd}) = ...
                    faceToVertexMeasure(meshData.(char(caseID_label)).corticalF, meshData.(char(caseID_label)).corticalV, ...
                    simData.(char(caseID_label)).(runTrials{runInd}).corticalC.(simVars{varInd}));
            end

            %Extract and reorder the vertices colour data so that it aligns
            %with just the outer cortical surface
            for varInd = 1:length(simVars)
                %Extract outer data
                simData.(char(caseID_label)).(runTrials{runInd}).outerCorticalCV.(simVars{varInd}) = ...
                    simData.(char(caseID_label)).(runTrials{runInd}).corticalCV.(simVars{varInd})(meshData.(char(caseID_label)).outerCorticalLogicVertices,:);
                %Reorder to match to registered surface
                simData.(char(caseID_label)).(runTrials{runInd}).outerCorticalCV.(simVars{varInd}) = ...
                    simData.(char(caseID_label)).(runTrials{runInd}).outerCorticalCV.(simVars{varInd})(meshData.(char(caseID_label)).ptCorrespondence,:);
                %Extract just the shaft colour data
                simData.(char(caseID_label)).(runTrials{runInd}).shaftCV.(simVars{varInd}) = ...
                    simData.(char(caseID_label)).(runTrials{runInd}).outerCorticalCV.(simVars{varInd})(meshData.(char(caseID_label)).shaftLogicVertices,:);
            end

% % %             %Take a similar approach to Bruce et al. (2022) and zero out elements
% % %             %within a specified radius of muscle attachments. Specifically, those
% % %             %elements within a 1cm radius of the soleus and 0.5cm radius of other
% % %             %muscles.
% % % 
% % %             %Loop through muscles
% % %             for muscleInd = 1:length(muscleList)
% % %                 %Set current radius on the basis of soleus vs. other muscles
% % %                 if contains(muscleList{muscleInd}, 'soleus')
% % %                     searchRadius = 1.0 * 10; %convert to mm
% % %                 else
% % %                     searchRadius = 0.5 * 10; %convert to mm
% % %                 end
% % %                 %Extract coordinate from tibial surface of current muscle
% % %                 muscCoord = meshData.(char(caseID_label)).outerCorticalV(muscleAttachmentInd.(muscleList{muscleInd}),:);
% % %                 %Loop through nodes on the shaft to zero out strain where appropriate
% % %                 for nodeInd = 1:length(meshData.(char(caseID_label)).shaftV)
% % %                     %Calculate distance from current shaft node to muscle coordinate 
% % %                     ptDist = distancePoints3d(meshData.(char(caseID_label)).shaftV(nodeInd,:), muscCoord);
% % %                     %Zero out simulation data at the node if less than the search radius
% % %                     if ptDist <= searchRadius
% % %                         for varInd = 1:length(simVars)
% % %                             simData.(char(caseID_label)).(runTrials{runInd}).shaftCV.(simVars{varInd})(nodeInd,:) = 0;
% % %                         end                        
% % %                     end
% % %                 end
% % %             end

            %% Create visualisation of simulation data on tibial shaft

            %Create visual of stress in shaft of variables
            %Add variable title names for each figure
            simVarTitles = [{['Principal Strain (Tension) with ',runNames{runInd},' run data for ',caseID]}, ...
                {['Principal Strain (Compression) with ',runNames{runInd},' run data for ',caseID]}, ...
                {['Effective Stress with ',runNames{runInd},' run data for ',caseID]}, ...
                {['Effective Strain with ',runNames{runInd},' run data for ',caseID]}];

            %Use subplots to create different perspectives
            for varInd = 1:length(simVars)

                %Create figure
                cFigure; hold on;

                %Loop through four views to create subplot
                for viewNo = 1:4                      

                    %Create subplot for current view
                    subplot(1,4,viewNo); hold on   

                    %Add generic whole bone view
                    hpPlain = gpatch(meshData.(char(caseID_label)).outerCorticalF, ...
                        meshData.(char(caseID_label)).outerCorticalV, [200/255 200/255 200/255], 'none', 0.2);
                    %Add variable heatmap
                    hpMap = gpatch(meshData.(char(caseID_label)).shaftF, ...
                        meshData.(char(caseID_label)).shaftV, ...
                        simData.(char(caseID_label)).(runTrials{runInd}).shaftCV.(simVars{varInd}), ...
                        'none', 1);
                    hpMap.FaceColor = 'Interp';
                    colormap bloodbone
                    %Scale axes colouring to better highlight stress areas
                    caxis([max(simData.(char(caseID_label)).(runTrials{runInd}).shaftCV.(simVars{varInd})) / 2 ...
                        max(simData.(char(caseID_label)).(runTrials{runInd}).shaftCV.(simVars{varInd}))])

                    %Set axis view
                    axis equal; axis tight; view(0,90);                    
                    %Rotate surfaces for current view
                    rotate(hpPlain,[0 1 0], surfaceRot(viewNo));
                    rotate(hpMap,[0 1 0], surfaceRot(viewNo));
                    %Set axis parameters
                    camlight headlight; axis off        
                    %Add colorbar on last view
                    if viewNo == 4
                        colorbar
                    end

                    %Add title
                    title([viewLabel{viewNo},' View'], 'FontSize', 12);

                end

                %Add overall figure title for variable
                figTitle = suptitle(simVarTitles{varInd});
                set(figTitle, 'FontSize', 14, 'FontWeight', 'bold');

            %Save figure
            export_fig(['..\Results\Figures\simulationOutputs\',caseID,'_',runTrials{runInd},'_',simVars{varInd},'.png'],'-m3');

            %Close figure
            close all

            end            

        end

        %% Save the processed results to file

        %Extract the current cases data to temporary space for saving
        meshDataCase = meshData.(char(caseID_label));
        simDataCase = simData.(char(caseID_label));

        %Save desired variables to file
        % % % save(['..\Results\processedSimulationData_',caseID,'.mat'], 'meshDataCase', 'simDataCase');
        save(['..\Results\processedSimulationData_',caseID,'_meshData.mat'], 'meshDataCase');
        save(['..\Results\processedSimulationData_',caseID,'_simData.mat'], 'simDataCase');

        %Cleanup variables that will take up substantial memory in the
        %loading of simulation data
        clear meshDataCase simDataCase meshData simData

    end
    
    %% Update the waitbar

    %Update after completing all cases
    waitbar(1, caseWaitBar, 'All data compiled & saved. Let''s go get a drink.');


%% ----- end of compileSimulations.m ----- %% 