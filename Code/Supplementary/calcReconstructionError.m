function [meanError] = calcReconstructionError(pcScores, shapeModel, origPoints)

    %% This script serves as the function input to the fmincon procedure
    %  where we aim to identify the PC scores that minimise the point error
    %  in the reconstruction of the tibia surface. The function
    %  reconstructs the surface and calculates the error
    
    %  Inputs:
    %       pcScores - n x 1 array of PC scores to use in reconstruction
    %       shapeModel - the structure with shape model data for reconstructing surface
    %       origPoints - n x 3 array of the XYZ original points to compare to  
    
    %% Set-up
    
    %Calculate the number of PCs to use based on the score input
    nPC = length(pcScores);
    
    %% Reconstruct surface
    
    %Use the provided scores and shape model loadings to reconstruct surface
    %Reconstruct the tibia points for new surface
    reconstructedPoints = pcScores' * ...
        shapeModel.loadings(:,1:nPC)' + ...
        shapeModel.mean;

    %Reshape points to visualise
    reconstructedV = transpose(reshape(reconstructedPoints,...
        [3, length(reconstructedPoints)/3]));
        
% % %     %Visualise original vs. reconstructed
% % %     cFigure; hold on;
% % %     gpatch(shapeModel.F, origPoints, 'gw', 'k', 0.3);
% % %     gpatch(shapeModel.F, reconstructedV, 'rw', 'k', 1);
% % %     axisGeom;

    %Calculate error between reconstructed and original surface
    reconError = distancePoints3d(origPoints, reconstructedV);
    
    %Convert error to colour map for visualisation
    errorColF = vertexToFaceMeasure(shapeModel.F, reconError);
    errorColV = faceToVertexMeasure(shapeModel.F, origPoints, errorColF);
    
% % %     %Visualise error
% % %     cFigure; hold on;
% % %     hpOrig = gpatch(shapeModel.F, origPoints, [200/255 200/255 200/255], 'none', 0.5);
% % %     hp = gpatch(shapeModel.F, reconstructedV, errorColV, 'none', 1);
% % %     hp.FaceColor = 'Interp'; colormap viridis
% % %     axis equal; axis tight; axis off
% % %     view(0,90);
% % %     colorbar; camlight headlight
    
    %Calculate mean error across points
    meanError = mean(reconError);
    
    %%%% TODO: consider addition of rotation component to get better match?
    %%%% Or consider a procrustes - does this help with PC scores though?
    
    
    
    
    
end