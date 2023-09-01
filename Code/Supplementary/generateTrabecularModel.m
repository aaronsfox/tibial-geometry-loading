function [trabF, trabV] = generateTrabecularModel(tibiaShapeModel, trabShapeModel, tibiaF, tibiaV, smoothIterations)

    %% This script essentially replicates the process from our shape model paper
    %  that generates a trabecular model from a tibial surface. For details
    %  please see the below paper/data repository:
    %
    %  TODO: add paper and repo details...
    %
    %  Inputs:
    %       tibiaShapeModel - the structure with shape model data for the trabecular  
    %       trabshapeModel - the structure with shape model data for the trabecular 
    %       tibiaF - faces for the tibia surface being used with prediction
    %       tibiaV - vertices for the tibia surface being used with prediction
    %       smoothIterations - number of iterations for smoothing predicted trabecular surface
    %           defaults to zero as this will ensure all points lie within tibia
    %           recommend using zero if script generates error of points still outside as this is a result of smoothing
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
    
    %Create wait bar
    wbar = waitbar(0,'Predicting trabecular surface. This can take a moment. Please be patient with me...');
    
    %Set options for remeshing
    optionStruct_tib.nb_pts = 3500; %Set desired number of points
    optionStruct_tib.disp_on = 0; % Turn off command window text display

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
    optCPD_rig.scale = 0;           % turn off scaling
    
    %Set default for smoothing
    if nargin < 5
        smoothIterations = 0;
    end

    %% Create a regression model between tibia and trabecular

    %Here we create a regression model that use the retained tibia
    %shape model scores to predict each of retained trabecular shape model
    %scores. The theoretical basis for this is that the shape of the tibia can
    %be effective in predicting the shape of the trabecular.
    %
    %This process and it's accuracy was established in our original paper
    %of the shape model data.

    %Get the retained PC scores from the tibia to use as predictors
    %Leave out appropriate case
    tibiaPCs = tibiaShapeModel.score(:,1:tibiaShapeModel.retainPCs);

    %Generate regression models for each of the retained trabecular PCs
    for predictPC = 1:trabShapeModel.retainPCs

        %Get the current scores as training for the regression model
        %Leave out appropriate case
        trabPCs = trabShapeModel.score(:,predictPC);

        %Fit the linear model on the data
        linearModel{predictPC} = fitlm(tibiaPCs, trabPCs);

    end

    %% Apply trabecular prediction to the input surface
    
    %Update waitbar
    waitbar((1/3), wbar, 'Predicting the trabecular shape. This is hard work, trust me...');
    
% % %     %Remesh the surfaces to match the shape model points
% % %     %This also avoids any weirdness with the surface
% % %     [tibiaF, tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);
% % %     [tibiaF, tibiaV] = mergeVertices(tibiaF, tibiaV);

    %%% Given we're using shape model data for these tibia models we can be
    %%% fairly confident that they are relatively well aligned in the first
    %%% place
    
% % %     %Apply a simple procrustes transform to maximise chance of best alignment
% % %     [~, tibiaV] = procrustes(tibiaShapeModel.meanPoints, tibiaV, ...
% % %         'scaling', false, 'reflection', false);

    %%% Given we're using the shape model data here we don't then need to
    %%% align it to the original shape model

% % %     %To map the new surface to the shape model it's important to identify the
% % %     %corresponding points on the surface to that of the shape model.
% % % 
% % %     %Perform non-rigid registration using CPD to align surfaces
% % %     [cpdTformTib] = cpd_register(tibiaShapeModel.meanPoints, tibiaV, optCPD_rig);
% % % 
% % %     %Identify the matching points in the target mesh against the
% % %     %registration to identify the corresponding point indices
% % %     regSortIdxTib = knnsearch(cpdTformTib.Y, tibiaShapeModel.meanPoints);
% % % 
% % %     %Sort the registered points so they align with the reference mesh
% % %     tibiaV_shapeModel = tibiaV(regSortIdxTib,:);

    %Visualise to confirm that points now correspond. This can be done
    %by using the target mesh faces against the registered moving
    %points. The surface meshes should essentially come out looking the
    %same (i.e. overlapped) despite using different faces
% % %     cFigure; hold on;
% % %     gpatch(tibiaF, tibiaV, 'rw', 'none', 0.3);
% % %     gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'k');
% % %     title('Registered Tibia');
% % %     axisGeom;

    %Reshape the registered data points and remove the mean
% % %     newData = reshape(tibiaV_shapeModel', [1, length(tibiaV_shapeModel)*3]) - tibiaShapeModel.mean;
    newData = reshape(tibiaV', [1, length(tibiaV)*3]) - tibiaShapeModel.mean;

    %Project the new data against the loadings to find the estimated scores
    %Set the upper and lower bounds on the score variables. Here we set them as
    %-5 / +5 standard deviations about the PC score mean (i.e. zero)
    sdRange = 5;
    %Loop through retained PCs
    for nPC = 1:length(tibiaShapeModel.varExplained)
        %Calculate scores and set to an initial guess variable
        x0(nPC,1) = dot(tibiaShapeModel.loadings(:,nPC), newData);
        %Set the lower bound
        lb(nPC,1) = std(tibiaShapeModel.score(:,nPC)) * -sdRange;
        %Set the upper bound
        ub(nPC,1) = std(tibiaShapeModel.score(:,nPC)) * sdRange;
    end

    %Create the function handle for fmincon
% % %     optFunc = @(pcScores)calcReconstructionError(pcScores, tibiaShapeModel, tibiaV_shapeModel);
    optFunc = @(pcScores)calcReconstructionError(pcScores, tibiaShapeModel, tibiaV);

    %Set fmincon options
    fminconOpts = optimoptions('fmincon');
    fminconOpts.MaxIterations = 5000; %up potential iterations
    fminconOpts.MaxFunctionEvaluations = 20000; %up potential function evals
    fminconOpts.StepTolerance = 1e-4; %scale optimisation tolerance to problem
    fminconOpts.Display = 'iter'; %set to display function iterations

    %Run fmincon to optimise PC scores for trabecular reconstruction
    [x,fval] = fmincon(optFunc, x0, [], [], [], [], lb, ub, [], fminconOpts);

    %Reconstruct using the optimised PC scores
    optReconstructed = x(1:tibiaShapeModel.retainPCs)' * ...
        tibiaShapeModel.loadings(:,1:tibiaShapeModel.retainPCs)' + ...
        tibiaShapeModel.mean;

    %Reshape reconstructed points
    optReconstructedV = reshape(optReconstructed', [3, length(optReconstructed)/3])';

% % %     %Visualise original vs. reconstructed
% % %     cFigure; hold on;
% % % % % %     gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % %     gpatch(tibiaShapeModel.F, tibiaV, 'gw', 'none', 0.3);
% % %     gpatch(tibiaShapeModel.F, optReconstructedV, 'rw', 'k', 1);
% % %     axisGeom; camlight headlight

    %Reconstruct a trabecular surface from previously created linear model
    for predictPC = 1:trabShapeModel.retainPCs

        %Output predicted values for left out case
        newTrabPCs(1, predictPC) = predict(linearModel{predictPC}, x(1:tibiaShapeModel.retainPCs)');

    end

    %Reconstruct the new trabecular using the predicted PCs
    newTrabReconstructed = newTrabPCs(:,1:trabShapeModel.retainPCs) * ...
        trabShapeModel.loadings(:,1:trabShapeModel.retainPCs)';

    %Reshape the predicted points to 3D
    newTrabV = reshape(newTrabReconstructed(1,:) + trabShapeModel.mean, ...
            [3, length(trabShapeModel.mean)/3])';

% % %     %Visualise new predicted trabecular within reconstructed tibia
% % %     cFigure; hold on;
% % %     gpatch(tibiaShapeModel.F, optReconstructedV, 'gw', 'none', 0.3);
% % %     gpatch(trabShapeModel.F, newTrabV, 'bw', 'k', 1);
% % %     axisGeom; camlight headlight

    %Update waitbar
    waitbar((2/3), wbar, 'Got the shape right, now just have to align it. I''m starting to sweat...');

% % %     %Rigidly align the predicted trabecular to the original tibia surface
% % %     [cpdRigTformTrabToOrig] = cpd_register(tibiaV, newTrabV, optCPD_rig);

% % %     %Create new trabecular aligned points
% % %     alignedNewTrabV = cpdRigTformTrabToOrig.Y;

    %Rigidly align the reconstructed tibia to the original tibia surface
    %This transformation will be used to shift the predicted trabecular
    %to align with the original surface
    [cpdRigTformToOrig] = cpd_register(tibiaV, optReconstructedV, optCPD_rig);
    
    %Apply the rigid transformation to the predicted trabecular to align it
    %with the original tibia
    alignedNewTrabV = cpdRigTformToOrig.s * newTrabV * cpdRigTformToOrig.R' + repmat(cpdRigTformToOrig.t,1,length(newTrabV))';

% % %     %Visualise new predicted aligned trabecular within original tibia
% % %     cFigure; hold on;
% % %     gpatch(tibiaShapeModel.F, tibiaV, 'gw', 'none', 0.3);
% % %     gpatch(trabShapeModel.F, alignedNewTrabV, 'bw', 'k', 1);
% % %     axisGeom; camlight headlight

    %Identify trabecular points within (and hence also outside) of original tibia
    trabPtsIn = intriangulation(tibiaV, tibiaShapeModel.F, alignedNewTrabV, 0);

    %Rescale trabecular if not all points are inside tibia
    if sum(~trabPtsIn) > 0

        %Generate a generic trabecular to support an initial guess of rescaling 
        [genericTrabF, genericTrabV] = generateGenericTrabecular(...
            tibiaV, tibiaShapeModel, 1, length(alignedNewTrabV));

        %Non-rigidly align the predicted trabecular to the generic trabecular
        %to identify matching points where corrections are needed
        [cpdRigTformTrab] = cpd_register(alignedNewTrabV, genericTrabV, optCPD);

        %Identify the matching points in the target mesh against the
        %registration to identify the corresponding point indices
        regSortIdxTrab = knnsearch(cpdRigTformTrab.Y, alignedNewTrabV);

        %Identify displacement of each point from predicted trabecular
        %Only extract displacement if the point in the tibia boundaries
        for ptNo = 1:length(alignedNewTrabV)
            if ~trabPtsIn(ptNo)
                translationVals(ptNo,:) = alignedNewTrabV(ptNo,:) - genericTrabV(regSortIdxTrab(ptNo),:);
            else
                translationVals(ptNo,:) = [0 0 0];
            end
        end

        %Create new trabecular vertices based on necessary translations
        alignedNewTrabV = alignedNewTrabV - translationVals;

        %Smooth out new trabecular surface if required
        %Note that this is the main reason how trabecular points can end up
        %outside the tibial surface. If this occurs then the will
        %progressively decrease smoothing iterations to zero until all
        %points lie within the boundaries
        while smoothIterations > 0
            
            %Set predicted points variable
            toSmoothV = alignedNewTrabV;
            
            %Run smoothing
            cPar.n = smoothIterations; %Number of iterations
            cPar.Method = 'LAP'; %Smooth method
            [predictedV] = patchSmooth(trabShapeModel.F, toSmoothV, [], cPar);
            
            %Check points within tibia
            trabPtsIn = intriangulation(tibiaV_shapeModel, tibiaShapeModel.F, predictedV, 0);
            if sum(~trabPtsIn) == 0
                %Break the while loop
                break
            else
                %Reduce smoothing iterations
                smoothIterations = smoothIterations - 1;
            end
            
        end
            
        %If smoothing iterations have progressed to zero then need to use
        %the original points
        if smoothIterations == 0
            predictedV = alignedNewTrabV;    
        end

    else

        %Use original
        predictedV = alignedNewTrabV;

    end

% % %     %Review newly fixed trabecular within original tibia
% % %     cFigure; hold on;
% % %     gpatch(tibiaShapeModel.F, tibiaV_shapeModel, 'gw', 'none', 0.3);
% % %     gpatch(trabShapeModel.F, predictedV, 'bw', 'k', 1);
% % % % % %     gpatch(trabShapeModel.F, alignedNewTrabV, 'bw', 'k', 1);
% % %     axisGeom; camlight headlight

    %One last final error check for if trabecular points are inside, as
    %later meshing functions will fail if this isn't the case
    trabPtsIn = intriangulation(tibiaV, tibiaShapeModel.F, predictedV, 0);
    if sum(~trabPtsIn) > 0
        error('Something has gone wrong with containing trabecular points within tibia...whoops...')
    end
    
    %Check that any tibia points aren't outside the trabecular surface
    tibiaPtsIn = intriangulation(predictedV, trabShapeModel.F, tibiaV, 0);
    
    %Correct if necessary
    if sum(tibiaPtsIn) > 0
        
        %Get point indices
        tibiaPtsInInd = find(tibiaPtsIn);
        
        %Loop through and correct points around
        for ptInd = 1:length(tibiaPtsInInd)
            
            %Get tibia point index as variable
            currTibiaPtInd = tibiaPtsInInd(ptInd);
            
            %Get the distance of the trabecular points to the current tibia point
            currPtDist = distancePoints3d(tibiaV(currTibiaPtInd,:), predictedV);
            
            %Find points less than the current 3mm radius
            manipulatePts = currPtDist <= 5;
            
            %Identify displacement of each point from predicted trabecular
            %Only extract displacement if the point in the tibia boundaries
            for ptNo = 1:length(alignedNewTrabV)
                if manipulatePts(ptNo)
                    translationVals2(ptNo,:) = predictedV(ptNo,:) - cpdRigTformTrab.Y(regSortIdxTrab(ptNo),:);
                else
                    translationVals2(ptNo,:) = [0 0 0];
                end
            end
            
            %Replace predicted points by applying translations extracted
            predictedV = predictedV - translationVals2;
            
        end
  
    end
    
    %Check final points step
    tibiaPtsInNew = intriangulation(predictedV, trabShapeModel.F, tibiaV, 0);
    if sum(tibiaPtsInNew) > 0
        error('Something has gone wrong with tibia points staying outside of trabecular...whoops...')
    end
    
    %% Store data to output faces and vertices variables
    
    %Set outputs
    trabF = trabShapeModel.F;
    trabV = predictedV;
    
    %Merge vertices to remove any issues with respect to duplicate points
    [trabF, trabV] = mergeVertices(trabF,trabV);

    %Update waitbar
    waitbar(1, wbar, 'Trabecular surface predicted with no issues (hopefully).');
    pause(1); close(wbar);

end