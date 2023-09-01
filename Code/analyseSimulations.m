function analyseSimulations()

    %% This script pulls together the compiled data for the sets of simulations
    %  from all cases for the subsequent analyses. The script therefore
    %  needs the outputs from the compileSimulations function for it to
    %  work. See the README.MD at top of repo and comments throughout this
    %  script for details.
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
    
    %Add supplementary code path
    addpath(genpath([pwd,'\Supplementary']));
    
    %Provide list of cases for checking against. Note that this is ordered
    %on the data that is in the shape model
    simulateCases = [147211; 102480; 102924; 103559; 103862; 107215; 107813;
        108421; 112802; 113033; 116070; 132433; 132592; 134434; 135418; 135881;
        136283; 140368; 140694; 141682; 144977; 145250; 146511; 149327; 170206;
        174615; 176087; 181140; 182401; 184171];

    %Set the running trial labels
	runTrials = {'Run_3', 'Run_4', 'Run_5'};
    
    %Create a variable for run names
    runNames = {'3 m^{.}s^{-1}', '4 m^{.}s^{-1}', '5 m^{.}s^{-1}'};
    
    %Set a list of variables to extract from the mesh and simulation data
    meshDataVars = [{'outerCorticalV'}; {'outerCorticalF'};
        {'minShaftVal'}; {'maxShaftVal'};
        {'shaftLogicVerticesLower'}; {'shaftLogicFacesLower'}; {'indBoundaryTop'};
        {'shaftLogicVerticesUpper'}; {'shaftLogicFacesUpper'}; {'indBoundaryBottom'};
        {'shaftV'}; {'shaftF'};
        {'shaftLogicVertices'}; {'shaftLogicFaces'}];
    simDataVars = [{'outerCorticalCV'}; {'shaftCV'}];
    
    %Create a variable of simulation data labels to extract and process
    simVars = [{'E1'}, {'E3'}, {'effectiveStress'}, {'effectiveStrain'}];
    
    %Create a variable for variable names
    simNames = [{'Principal Strain (Tension)'}, {'Principal Strain (Compression)'}, ...
        {'Effective Stress'}, {'Effective Strain'}];
    
    %Load the shape models required
    load('..\ShapeModels\tibia\tibiaShapeModel.mat');
    load('..\ShapeModels\tibia-fibula\tibiaFibulaShapeModel.mat');
    load('..\ShapeModels\tibia-plus-trabecular\tibTrabShapeModel.mat');
    
    %Define the generic shaft proportion length
    %Bruce et al. (2022) defines the tibia shaft at 20%-80% of the tibia
    %length. We use a similar definition here, however the 15%-75% length
    %appears to fit better for the area of interest.
    lowerShaftProp = 0.15; upperShaftProp = 0.75;
    
    %Set list to label reconstruction subplots with
    viewLabel = [{'Anterior'}, {'Lateral'}, {'Posterior'}, {'Medial'}];
    
    %Set array to rotate views of surfaces later in reconstruction process
    surfaceRot = [-90, 0, 90, 180];
    
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

    %% Load the compiled simulation data
        
    %Loop through cases and extract data
    for caseInd = 1:length(simulateCases)

        %Create a variable for the caseID
        caseID = ['case-',num2str(simulateCases(caseInd))];
        caseID_label = ['case',num2str(simulateCases(caseInd))];

        %Load the data for the current case
        load(['..\Results\processedSimulationData_',caseID,'.mat']);
        
        %Create a variable that collates the element numbers for the meshes
        tibialElementN(caseInd) = length(meshDataCase.meshOutputFullTibia.elements);
        fibulaElementN(caseInd) = length(meshDataCase.meshOutputFibula.elements);
        
        %Extract the necessary variables for analysis
        
        %Mesh data
        for varInd = 1:length(meshDataVars)
            meshData.(char(caseID_label)).(meshDataVars{varInd}) = ...
                meshDataCase.(meshDataVars{varInd});
        end
        clear meshDataCase
        
        %Sim data
        %Loop through run variables
        for runInd = 1:length(runTrials)
            for varInd = 1:length(simDataVars)
                simData.(char(caseID_label)).(runTrials{runInd}).(simDataVars{varInd}) = ...
                    simDataCase.(runTrials{runInd}).(simDataVars{varInd});
            end
        end
        clear simDataCase

    end
    
    %Display minimum for element data
    disp(['Minimum tibia elements number: ',num2str(min(tibialElementN))]);
    disp(['Minimum fibula elements number: ',num2str(min(fibulaElementN))]);
    
    %% Calculate simulation variables means across shaft area from dataset
    
    %Calculate mean tibial shape
    %Loop through and extract points from individual axes to a variable
    %Set variable for each axes
    ptAxes = [{'Vx'}, {'Vy'}, {'Vz'}];
    for caseInd = 1:length(simulateCases)
        caseID_label = ['case',num2str(simulateCases(caseInd))];
        for axInd = 1:length(ptAxes)
            tibiaPtsGroup.(ptAxes{axInd})(:,caseInd) = ...
                meshData.(char(caseID_label)).outerCorticalV(:,axInd);
        end
    end
    
    %Calculate and extract mean tibia points and faces (using 1st case)
    tibiaV_m = [mean(tibiaPtsGroup.Vx, 2), ...
        mean(tibiaPtsGroup.Vy, 2), ...
        mean(tibiaPtsGroup.Vz, 2)];
    tibiaF_m = meshData.case147211.outerCorticalF;
    
% % %     %Visualise mean surface
% % %     cFigure;
% % %     gpatch(tibiaF_m, tibiaV_m);
% % %     axisGeom;

    %Determine shaft boundaries for mean surface so that data can be
    %extracted for the group across these points
    
    %Calculate the nodal boundaries for the shaft length based on
    %the generic proportions (along the y-axis as this is the
    %longitudinal axis).
    minShaftVal_m = min(tibiaV_m(:,2)) + ((max(tibiaV_m(:,2)) - ...
        min(tibiaV_m(:,2))) * lowerShaftProp);
    maxShaftVal_m = min(tibiaV_m(:,2)) + ((max(tibiaV_m(:,2)) - ...
        min(tibiaV_m(:,2))) * upperShaftProp);

    %Extract shaft of tibial surface
    shaftData_m = extractMeanTibialShaftSurface(tibiaF_m, tibiaV_m, minShaftVal_m, maxShaftVal_m);

% % %     %Visualise mean surface with shaft extracted
% % %     cFigure;
% % %     gpatch(tibiaF_m, tibiaV_m, 'rw', 'none', 0.3);
% % %     gpatch(shaftData_m.shaftF, shaftData_m.shaftV);
% % %     axisGeom;
    
    %Extract simulation variables over shaft area
    for caseInd = 1:length(simulateCases)
        
        %Create case label
        caseID_label = ['case',num2str(simulateCases(caseInd))];
        
        %Loop through run trials
        for runInd = 1:length(runTrials)
            
            %Loop through simulation variables
            for varInd = 1:length(simVars)
                
                %Extract current case, trial and variable data to structure
                simDataGroup.shaftCV.(runTrials{runInd}).(simVars{varInd})(:,caseInd) = ...
                    simData.(char(caseID_label)).(runTrials{runInd}).outerCorticalCV.(simVars{varInd})(shaftData_m.shaftLogicVertices,1);
                
            end
            
        end
        
    end
    
    %Take a similar approach to Bruce et al. (2022) and zero out values
    %within a specified radius of muscle attachments. Specifically, those
    %elements within a 1cm radius of the soleus and 0.5cm radius of other
    %muscles.

    %Loop through muscles
    for muscleInd = 1:length(muscleList)
        %Set current radius on the basis of soleus vs. other muscles
        if contains(muscleList{muscleInd}, 'soleus')
            searchRadius = 1.0 * 10; %convert to mm
        else
            searchRadius = 0.5 * 10; %convert to mm
        end
        %Extract coordinate from tibial surface of current muscle
        muscCoord = tibiaFibulaShapeModel.meanPoints(muscleAttachmentInd.(muscleList{muscleInd}),:);
        %Identify point distances on mean surface from the muscle coordinate
        ptDist = distancePoints3d(shaftData_m.shaftV, muscCoord);
        %Identify points within search radius
        withinRadius = ptDist <= searchRadius;
        %Zero out points in simulation variables within muscle search radius
        %Loop through run trials
        for runInd = 1:length(runTrials)
            %Loop through simulation variables
            for varInd = 1:length(simVars)
                simDataGroup.shaftCV.(runTrials{runInd}).(simVars{varInd})(withinRadius,:) = 0;
            end
        end
    end
    
    %Calculate the group mean for each run trial and simulation variable
    %Loop through run trials
    for runInd = 1:length(runTrials)
        %Loop through simulation variables
        for varInd = 1:length(simVars)
            %Calculate group mean
            simData_m.shaftCV.(runTrials{runInd}).(simVars{varInd}) = ...
                mean(simDataGroup.shaftCV.(runTrials{runInd}).(simVars{varInd}),2);
        end
    end
    
    %% Visualise group mean for each variable
    
    %Note that colouring of simulation variables (i.e. stress/strain) is
    %normalised to the maximum values for the variables across all running
    %trials
    
    %Loop through simulation variables
    for varInd = 1:length(simVars)
        
        %Get overall max for variable across all trials
        cMax = max([max(simData_m.shaftCV.Run_3.(simVars{varInd})); ...
             max(simData_m.shaftCV.Run_4.(simVars{varInd})); ...
             max(simData_m.shaftCV.Run_5.(simVars{varInd}))]);

        %Create figure
        fig = cFigure; hold on;
        
        %Adjust size of figure (A4 size)
        fig.Units = 'centimeters';
        fig.Position = [1, 1, 1+21, 1+29.7];

        %Loop through run trials
        for runInd = 1:length(runTrials)
            
            %Loop through four views to create subplot
            for viewNo = 1:4  

                %Create subplot for current view
                %Note the subplot size is relative to fitting all 4 views
                %across the 3 running trials
                subplot(3,4,viewNo+(4*(runInd-1))); hold on
                
                %Add speed title on first column view
                if viewNo == 1
                    ylabel(runNames{runInd});
                end
                
                %Add generic whole bone view
                hpPlain = gpatch(tibiaF_m, tibiaV_m, [200/255 200/255 200/255], 'none', 0.2);
                %Add variable heatmap
                hpMap = gpatch(shaftData_m.shaftF, shaftData_m.shaftV, ...
                   simData_m.shaftCV.(runTrials{runInd}).(simVars{varInd}), ...
                    'none', 1);
                hpMap.FaceColor = 'Interp';
                colormap bloodbone
                %Scale axes colouring to better highlight stress areas
                %Note that this is scaled to the relative to maximum of
                %variable across ALL run trials
                caxis([cMax / 2, cMax]);

                %Set axis view
                axis equal; axis tight; view(0,90);                    
                %Rotate surfaces for current view
                rotate(hpPlain,[0 1 0], surfaceRot(viewNo));
                rotate(hpMap,[0 1 0], surfaceRot(viewNo));
                %Set axis parameters
                camlight headlight;
                hAx = get(gca);
                %Turn off axis lines without removing labels
                set(gca,'XColor','none')
                set(gca,'YColor','none')
                hAx.YAxis.Label.Color=[0 0 0];
                hAx.YAxis.Label.Visible='on';
      
                %Add colorbar on last view
                if viewNo == 4
                    axPos = get(gca(), 'Position');
                    colorbar('eastoutside');
                    set(gca(), 'Position', axPos);
                end

                %Add title
                title(viewLabel{viewNo}, 'FontSize', 8);

            end

        end
    
        %Add overall figure title for variable
        figTitle = suptitle(simNames{varInd});
        set(figTitle, 'FontSize', 14, 'FontWeight', 'bold');

        %Save figure
        export_fig(['..\Results\Figures\meanOutputs\groupData_',simVars{varInd},'.png'],'-m3');

        %Close figure
        close all
        
    end
    
    %% Run regression analysis on desired components against effective strain
    
    %Note that certain components have been selected here on the basis of 
    %minimising the number of statistical tests and avoiding repeating the 
    %analysis on correlated components that represent the same shape
    %characteristics. The components selected for analysis were:
    %
    %   Tibia-fibula model PC1, PC2, PC3, PC4, PC5, PC6, PC7
    %   Tibia-trabecular model PC2, PC3, PC4
    
    %Set a variable for these two models on whether to run analysis
    tibiaFibulaShapeModel.runRegressionOnComponent = [true, true, true, true, true, true, true];
    tibTrabShapeModel.runRegressionOnComponent = [false, true, true, true, false];
    
    %Set parameters for statistical testing
    alpha = 0.05;
    two_tailed = true; %true;
    iterations = 1000; %10000;
    
    %Loop through run trials
	for runInd = 1:length(runTrials)
        
        %Extract variable data
        %Note that the variable of interest here is effective strain
        y = simDataGroup.shaftCV.(runTrials{runInd}).effectiveStrain';

        %Remove non-zero variance elements that have been zeroed out
        %Elements that have zero variance (i.e. all the same) across the
        %entire dataset need to be removed
        iY = std(y, [], 1) > eps;

        %Get y non-zero elements
        ynz       = y(:,iY);
        
        %%%%% Tibia-fibula shape model %%%%%
    
        %Loop through principal components
        for nPC = 1:tibiaFibulaShapeModel.retainPCs
            
            %Check whether to run analysis on this component
            if tibiaFibulaShapeModel.runRegressionOnComponent(nPC)
            
                %Extract current PC
                %Normalise to z-scores for ease of interpretation
                x = zscore(tibiaFibulaShapeModel.score(:,nPC));

                %Conduct non-parametric linear regression test
                rng(12345);
                snpm = spm1d.stats.nonparam.regress(ynz, x);
                snpmi = snpm.inference(alpha, 'two_tailed', two_tailed, 'iterations', iterations);
                
                %Display regression output and save to text file
                diary(['..\Results\Statistics\tibiaFibulaShapeModel_PC',num2str(nPC),'_',runTrials{runInd},'_effectiveStrain_regression.txt']);
                diary on
                disp(['Non-Parametric Regression Results: Tibia-Fibula PC', num2str(nPC)]);
                disp(snpmi);
                diary off

                %Organize results for plotting
                %This reinserts the zero variance nodes back in
                znz = snpmi.z;
                zstar = snpmi.zstar;
                z = zeros(1,size(simDataGroup.shaftCV.(runTrials{runInd}).effectiveStrain',2));
                z(iY) = znz;
                zi = z;
                zi(abs(z) < zstar) = 0;

                %Create plot of SPM t-values and inferential values < alpha
                cFigure;

                %Loop through views
                for viewNo = 1:4

                    %Raw SPM values
                    subplot(2,4,viewNo); hold on
                    %Add plain aspect
                    hpPlain = gpatch(tibiaF_m, tibiaV_m, [200/255 200/255 200/255], 'none', 0.2);
                    %Add z-statistic heatmap
                    hpZ = gpatch(shaftData_m.shaftF, shaftData_m.shaftV, z', 'none');
                    hpZ.FaceColor = 'Interp';
                    colormap(warmcold); 
                    caxis([zstar*-1 zstar]); %scale color axis so it's centred around zero & relative to zstar
                    %Set axis view
                    axis equal; axis tight; view(0,90);
                    %Rotate surfaces for current view
                    rotate(hpPlain,[0 1 0], surfaceRot(viewNo));
                    rotate(hpZ,[0 1 0], surfaceRot(viewNo));
                    camlight headlight;
                    %Add title
                    title(viewLabel{viewNo}, 'FontSize', 12);
                    %Add y-label
                    ylabel('SPM \{t\}');
                    %Turn off axis lines without removing labels
                    hAx = get(gca);
                    set(gca,'XColor','none')
                    set(gca,'YColor','none')
                    hAx.YAxis.Label.Color=[0 0 0];
                    hAx.YAxis.Label.Visible='on';

                    %Add colorbar if last view number
                    if viewNo == 4
                        axPos = get(gca(), 'Position');
                        colorbar('eastoutside');
                        set(gca(), 'Position', axPos);
                    end

                    %Add inferential statistic values
                    subplot(2,4,viewNo+4); hold on
                    %Add plain aspect
                    hpPlain = gpatch(tibiaF_m, tibiaV_m, [200/255 200/255 200/255], 'none', 0.2);
                    %Add inferential z-statistic heatmap
                    hpZ = gpatch(shaftData_m.shaftF, shaftData_m.shaftV, zi', 'none');
                    hpZ.FaceColor = 'Interp';
                    colormap(warmcold); 
                    caxis([zstar*-1 zstar]); %scale color axis so it's centred around zero & relative to zstar
                    %Set axis view
                    axis equal; axis tight; view(0,90);
                    %Rotate surfaces for current view
                    rotate(hpPlain,[0 1 0], surfaceRot(viewNo));
                    rotate(hpZ,[0 1 0], surfaceRot(viewNo));
                    camlight headlight;
                    %Add y-label
                    ylabel(['SPM \{t\} > Critical Threshold (\alpha < ',num2str(alpha),')']);
                    %Turn off axis lines without removing labels
                    hAx = get(gca);
                    set(gca,'XColor','none')
                    set(gca,'YColor','none')
                    hAx.YAxis.Label.Color=[0 0 0];
                    hAx.YAxis.Label.Visible='on';

                    %Add colorbar if last view number
                    if viewNo == 4
                        axPos = get(gca(), 'Position');
                        colorbar('eastoutside');
                        set(gca(), 'Position', axPos);
                    end

                end

                %Add figure title
                %Add some notation to denote statistical significance
                if sum(zi) == 0
                    %Non-statistically significant results
                    figTitle = suptitle(['Tibia-Fibula Shape Model PC',num2str(nPC),' Regression vs. Effective Strain at ',runNames{runInd},' ({\itn.s.})']);
                else
                    %Statistically significant result
                    figTitle = suptitle(['Tibia-Fibula Shape Model PC',num2str(nPC),' Regression vs. Effective Strain at ',runNames{runInd},' *']);
                end
                set(figTitle, 'FontSize', 14, 'FontWeight', 'bold');

                %Save figure
                export_fig(['..\Results\Figures\statisticalOutputs\tibiaFibulaShapeModel_PC',num2str(nPC),'_',runTrials{runInd},'_effectiveStrain_regression.png'],'-m3');
                
                %Close figure
                close all

                %Save statistical results
                save(['..\Results\Statistics\tibiaFibulaShapeModel_PC',num2str(nPC),'_',runTrials{runInd},'_effectiveStress_regression.mat'], ...
                    'snpm','snpmi');

                %Cleanup
                clear x snpm snpmi znz zstar z zi
                
            end
            
        end
        
        %%%%% Tibia-trabecular shape model %%%%%
    
        %Loop through principal components
        for nPC = 1:tibTrabShapeModel.retainPCs
            
            %Check whether to run analysis on this component
            if tibTrabShapeModel.runRegressionOnComponent(nPC)
            
                %Extract current PC
                %Normalise to z-scores for ease of interpretation
                x = zscore(tibTrabShapeModel.score(:,nPC));

                %Conduct non-parametric linear regression test
                rng(12345);
                snpm = spm1d.stats.nonparam.regress(ynz, x);
                snpmi = snpm.inference(alpha, 'two_tailed', two_tailed, 'iterations', iterations);
                
                %Display regression output and save to text file
                diary(['..\Results\Statistics\tibTrabShapeModel_PC',num2str(nPC),'_',runTrials{runInd},'_effectiveStrain_regression.txt']);
                diary on
                disp(['Non-Parametric Regression Results: Cortical-Trabecular PC', num2str(nPC)]);
                disp(snpmi);
                diary off

                %Organize results for plotting
                %This reinserts the zero variance nodes back in
                znz = snpmi.z;
                zstar = snpmi.zstar;
                z = zeros(1,size(simDataGroup.shaftCV.(runTrials{runInd}).effectiveStrain',2));
                z(iY) = znz;
                zi = z;
                zi(abs(z) < zstar) = 0;

                %Create plot of SPM t-values and inferential values < alpha
                cFigure;

                %Loop through views
                for viewNo = 1:4

                    %Raw SPM values
                    subplot(2,4,viewNo); hold on
                    %Add plain aspect
                    hpPlain = gpatch(tibiaF_m, tibiaV_m, [200/255 200/255 200/255], 'none', 0.2);
                    %Add z-statistic heatmap
                    hpZ = gpatch(shaftData_m.shaftF, shaftData_m.shaftV, z', 'none');
                    hpZ.FaceColor = 'Interp';
                    colormap(warmcold); 
                    caxis([zstar*-1 zstar]); %scale color axis so it's centred around zero & relative to zstar
                    %Set axis view
                    axis equal; axis tight; view(0,90);
                    %Rotate surfaces for current view
                    rotate(hpPlain,[0 1 0], surfaceRot(viewNo));
                    rotate(hpZ,[0 1 0], surfaceRot(viewNo));
                    camlight headlight;
                    %Add title
                    title(viewLabel{viewNo}, 'FontSize', 12);
                    %Add y-label
                    ylabel('SPM \{t\}');
                    %Turn off axis lines without removing labels
                    hAx = get(gca);
                    set(gca,'XColor','none')
                    set(gca,'YColor','none')
                    hAx.YAxis.Label.Color=[0 0 0];
                    hAx.YAxis.Label.Visible='on';

                    %Add colorbar if last view number
                    if viewNo == 4
                        axPos = get(gca(), 'Position');
                        colorbar('eastoutside');
                        set(gca(), 'Position', axPos);
                    end

                    %Add inferential statistic values
                    subplot(2,4,viewNo+4); hold on
                    %Add plain aspect
                    hpPlain = gpatch(tibiaF_m, tibiaV_m, [200/255 200/255 200/255], 'none', 0.2);
                    %Add inferential z-statistic heatmap
                    hpZ = gpatch(shaftData_m.shaftF, shaftData_m.shaftV, zi', 'none');
                    hpZ.FaceColor = 'Interp';
                    colormap(warmcold); 
                    caxis([zstar*-1 zstar]); %scale color axis so it's centred around zero & relative to zstar
                    %Set axis view
                    axis equal; axis tight; view(0,90);
                    %Rotate surfaces for current view
                    rotate(hpPlain,[0 1 0], surfaceRot(viewNo));
                    rotate(hpZ,[0 1 0], surfaceRot(viewNo));
                    camlight headlight;
                    %Add y-label
                    ylabel(['SPM \{t\} > Critical Threshold (\alpha < ',num2str(alpha),')']);
                    %Turn off axis lines without removing labels
                    hAx = get(gca);
                    set(gca,'XColor','none')
                    set(gca,'YColor','none')
                    hAx.YAxis.Label.Color=[0 0 0];
                    hAx.YAxis.Label.Visible='on';

                    %Add colorbar if last view number
                    if viewNo == 4
                        axPos = get(gca(), 'Position');
                        colorbar('eastoutside');
                        set(gca(), 'Position', axPos);
                    end

                end

                %Add figure title
                %Add some notation to denote statistical significance
                if sum(zi) == 0
                    %Non-statistically significant results
                    figTitle = suptitle(['Cortical-Trabecular Shape Model PC',num2str(nPC),' Regression vs. Effective Strain at ',runNames{runInd},' ({\itn.s.})']);
                else
                    %Statistically significant result
                    figTitle = suptitle(['Cortical-Trabecular Shape Model PC',num2str(nPC),' Regression vs. Effective Strain at ',runNames{runInd},' *']);
                end
                set(figTitle, 'FontSize', 14, 'FontWeight', 'bold');

                %Save figure
                export_fig(['..\Results\Figures\statisticalOutputs\tibTrabShapeModel_PC',num2str(nPC),'_',runTrials{runInd},'_effectiveStrain_regression.png'],'-m3');
                
                %Close figure
                close all

                %Save statistical results
                save(['..\Results\Statistics\tibiaFibulaShapeModel_PC',num2str(nPC),'_',runTrials{runInd},'_effectiveStress_regression.mat'], ...
                    'snpm','snpmi');

                %Cleanup
                clear x snpm snpmi znz zstar z zi
                
            end
            
        end
        
        %Cleanup
        clear y iY ynz
        
    end
    
    %% Create specific figures for paper
    
    %% PC6 of tibia-fibula model
    
    %Set parameters for statistical testing
    alpha = 0.05;
    two_tailed = true; %true;
    iterations = 1000; %10000;
    
    %Set the PC
    nPC = 6;
    
    %Create figure
    pcFig = cFigure; hold on
    
    %Adjust to A4 size
    pcFig.Units = 'centimeters';
    pcFig.Position = [1, 1, 22, 30];
    
    %Loop through run trials
	for runInd = 1:length(runTrials)
        
        %Extract variable data
        %Note that the variable of interest here is effective strain
        y = simDataGroup.shaftCV.(runTrials{runInd}).effectiveStrain';

        %Remove non-zero variance elements that have been zeroed out
        %Elements that have zero variance (i.e. all the same) across the
        %entire dataset need to be removed
        iY = std(y, [], 1) > eps;

        %Get y non-zero elements
        ynz       = y(:,iY);
        
        %Extract current PC
        %Normalise to z-scores for ease of interpretation
        x = zscore(tibiaFibulaShapeModel.score(:,nPC));
        
        %Conduct non-parametric linear regression test
        rng(12345);
        snpm = spm1d.stats.nonparam.regress(ynz, x);
        snpmi = snpm.inference(alpha, 'two_tailed', two_tailed, 'iterations', iterations);

        %Organize results for plotting
        %This reinserts the zero variance nodes back in
        znz = snpmi.z;
        zstar = snpmi.zstar;
        z = zeros(1,size(simDataGroup.shaftCV.(runTrials{runInd}).effectiveStrain',2));
        z(iY) = znz;
        zi = z;
        zi(abs(z) < zstar) = 0;

        %Loop through views
        for viewNo = 1:4

            %Add subplot
            subplot(3,4,((runInd-1)*4)+viewNo); hold on
 
            %Add plain aspect
            hpPlain = gpatch(tibiaF_m, tibiaV_m, [200/255 200/255 200/255], 'none', 0.2);
            %Add inferential z-statistic heatmap
            hpZ = gpatch(shaftData_m.shaftF, shaftData_m.shaftV, zi', 'none');
            hpZ.FaceColor = 'Interp';
            colormap(warmcold); 
            caxis([zstar*-1 zstar]); %scale color axis so it's centred around zero & relative to zstar
            %Set axis view
            axis equal; axis tight; view(0,90);
            %Rotate surfaces for current view
            rotate(hpPlain,[0 1 0], surfaceRot(viewNo));
            rotate(hpZ,[0 1 0], surfaceRot(viewNo));
            camlight headlight;
            %Add y-label
            yLab = ylabel(['SPM \{t\} (\alpha < ',num2str(alpha),')']);
            %Turn off axis lines without removing labels
            hAx = get(gca);
            set(gca,'XColor','none')
            set(gca,'YColor','none')
            hAx.YAxis.Label.Color=[0 0 0];
            hAx.YAxis.Label.Visible='on';
            
            %Add plot title
            title(viewLabel{viewNo}, 'FontSize', 8, 'FontWeight', 'bold');
            
            %Add running speed if first view
            if viewNo == 1
                text(yLab.Position(1) - 220, yLab.Position(2), runNames{runInd}, ...
                    'FontSize', 12, 'FontWeight', 'bold');
            end

            %Add colorbar if last view number
            if viewNo == 4
                axPos = get(gca(), 'Position');
                colorbar('eastoutside');
                set(gca(), 'Position', axPos);
            end

        end

        %Save figure
        export_fig(['..\Results\Figures\selectOutputs\tibiaFibulaShapeModel_PC',num2str(nPC),'_effectiveStrain_regression_allSpeeds.png'],'-m3');

        %Close figure
        close all

        %Cleanup
        clear x snpm snpmi znz zstar z zi
                
    end

    %% PC2 of tibia-trabecular model
    
    %Set parameters for statistical testing
    alpha = 0.05;
    two_tailed = true; %true;
    iterations = 1000; %10000;
    
    %Set the PC
    nPC = 2;
    
    %Create figure
    pcFig = cFigure; hold on
    
    %Adjust to A4 size
    pcFig.Units = 'centimeters';
    pcFig.Position = [1, 1, 5, 7];
    
    %Set specific run index
    runInd = 3;
        
    %Extract variable data
    %Note that the variable of interest here is effective strain
    y = simDataGroup.shaftCV.(runTrials{runInd}).effectiveStrain';

    %Remove non-zero variance elements that have been zeroed out
    %Elements that have zero variance (i.e. all the same) across the
    %entire dataset need to be removed
    iY = std(y, [], 1) > eps;

    %Get y non-zero elements
    ynz       = y(:,iY);

    %Extract current PC
    %Normalise to z-scores for ease of interpretation
    x = zscore(tibTrabShapeModel.score(:,nPC));

    %Conduct non-parametric linear regression test
    rng(12345);
    snpm = spm1d.stats.nonparam.regress(ynz, x);
    snpmi = snpm.inference(alpha, 'two_tailed', two_tailed, 'iterations', iterations);

    %Organize results for plotting
    %This reinserts the zero variance nodes back in
    znz = snpmi.z;
    zstar = snpmi.zstar;
    z = zeros(1,size(simDataGroup.shaftCV.(runTrials{runInd}).effectiveStrain',2));
    z(iY) = znz;
    zi = z;
    zi(abs(z) < zstar) = 0;

    %Set the specific view number
    viewNo = 4;

    %Add plain aspect
    hpPlain = gpatch(tibiaF_m, tibiaV_m, [200/255 200/255 200/255], 'none', 0.2);
    %Add inferential z-statistic heatmap
    hpZ = gpatch(shaftData_m.shaftF, shaftData_m.shaftV, zi', 'none');
    hpZ.FaceColor = 'Interp';
    colormap(warmcold); 
    caxis([zstar*-1 zstar]); %scale color axis so it's centred around zero & relative to zstar
    %Set axis view
    axis equal; axis tight; view(0,90);
    %Rotate surfaces for current view
    rotate(hpPlain,[0 1 0], surfaceRot(viewNo));
    rotate(hpZ,[0 1 0], surfaceRot(viewNo));
    camlight headlight;
    %Add y-label
    yLab = ylabel(['SPM \{t\} (\alpha < ',num2str(alpha),')']);
    %Turn off axis lines without removing labels
    hAx = get(gca);
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    hAx.YAxis.Label.Color=[0 0 0];
    hAx.YAxis.Label.Visible='on';

    %Add plot title
    title(viewLabel{viewNo}, 'FontSize', 8, 'FontWeight', 'bold');

    %Add the colourbar
    axPos = get(gca(), 'Position');
    colorbar('eastoutside');
    set(gca(), 'Position', axPos);

    %Save figure
    export_fig(['..\Results\Figures\selectOutputs\tibiaTrabShapeModel_PC',num2str(nPC),'_effectiveStrain_regression_5ms.png'],'-m3');

    %Close figure
    close all

    %Cleanup
    clear x snpm snpmi znz zstar z zi
    
%% ----- end of analyseSimulations.m ----- %%

end