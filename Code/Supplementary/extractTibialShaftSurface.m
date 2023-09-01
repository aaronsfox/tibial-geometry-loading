function [meshData] = extractTibialShaftSurface(meshData)

    %% This script extracts the tibial shaft surface from the outer tibial
    %  surface and associated mesh inputs
    %
    %  Inputs:
    %       meshData - structure containing cases relevant mesh/surface data
    %
    %  Outputs:
    %       meshData - updated input structure with data related to shaft surface
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

    %% Extract the tibial shaft surface
    
    %Create logic for cutting away the ends of the bone
    %Note that longitudinal axis for these surfaces is the Y-axis

    %Start with under maximum shaft value
    meshData.shaftLogicVerticesLower = meshData.outerCorticalV(:,2) < meshData.maxShaftVal;
    meshData.shaftLogicFacesLower = all(meshData.shaftLogicVerticesLower(meshData.outerCorticalF), 2);
    meshData.shaftLogicFacesLower = triSurfLogicSharpFix(meshData.outerCorticalF, meshData.shaftLogicFacesLower, 3); %Altered logic so it is smoother

    %Cut away the faces for the lower part of the shaft
    %Use a temporary throughout here to later manipulate the actual surface
    shaftF = meshData.outerCorticalF(meshData.shaftLogicFacesLower,:);

    %Clean up unused points
    [shaftF, shaftV] = patchCleanUnused(shaftF, meshData.outerCorticalV);

    %Extract boundaries of surface
    shaftEbLower = patchBoundary(shaftF, shaftV);
    indBoundaryEbLower = edgeListToCurve(shaftEbLower); %Convert boundary edges to a curve list
    indBoundaryEbLower = indBoundaryEbLower(1:end-1); %Trim off last point since it is equal to first on a closed loop

    %Identify the indices of those around the top boundary in the original
    %shape model tibia surface
    meshData.indBoundaryTop = logical(zeros(length(meshData.outerCorticalV),1));
    for indPt = 1:length(indBoundaryEbLower)
        ind = find(distancePoints3d(meshData.outerCorticalV, shaftV(indBoundaryEbLower(indPt),:)) == 0);
        meshData.indBoundaryTop(ind,1) = true;
    end

    %Repeat the above process to get the points along the lower boundary
    meshData.shaftLogicVerticesUpper = meshData.outerCorticalV(:,2) > meshData.minShaftVal;
    meshData.shaftLogicFacesUpper = all(meshData.shaftLogicVerticesUpper(meshData.outerCorticalF), 2);
    meshData.shaftLogicFacesUpper = triSurfLogicSharpFix(meshData.outerCorticalF, meshData.shaftLogicFacesUpper, 3); %Altered logic so it is smoother

    %Use a temporary throughout here to later manipulate the actual surface
    shaftF = meshData.outerCorticalF(meshData.shaftLogicFacesUpper,:);

    %Clean up unused points
    [shaftF, shaftV] = patchCleanUnused(shaftF, meshData.outerCorticalV);

    %Extract boundaries of surface
    shaftEbUpper = patchBoundary(shaftF, shaftV);
    indBoundaryEbUpper = edgeListToCurve(shaftEbUpper); %Convert boundary edges to a curve list
    indBoundaryEbUpper = indBoundaryEbUpper(1:end-1); %Trim off last point since it is equal to first on a closed loop

    %Identify the indices of those around the bottom boundary in the original
    %shape model tibia surface
    meshData.indBoundaryBottom = logical(zeros(length(meshData.outerCorticalV),1));
    for indPt = 1:length(indBoundaryEbUpper)
        ind = find(distancePoints3d(meshData.outerCorticalV, shaftV(indBoundaryEbUpper(indPt),:)) == 0);
        meshData.indBoundaryBottom(ind,1) = true;
    end

    %Create a shaft vertices variable and shift the boundaries to the desired levels
    meshData.shaftV = meshData.outerCorticalV;
    meshData.shaftV(meshData.indBoundaryTop,2) = meshData.maxShaftVal;
    meshData.shaftV(meshData.indBoundaryBottom,2) = meshData.minShaftVal;

    %Extract the shaft while keeping the boundaries at the desired levels
    %Create the logicals for extracting vertices and faces
    meshData.shaftLogicVertices = meshData.shaftV(:,2) <= meshData.maxShaftVal & ...
        meshData.shaftV(:,2) >= meshData.minShaftVal;
    meshData.shaftLogicFaces = all(meshData.shaftLogicVertices(meshData.outerCorticalF), 2);

    %Cut away the faces for the shaft
    meshData.shaftF = meshData.outerCorticalF(meshData.shaftLogicFaces,:);

    %Clean up unused points
    [meshData.shaftF, meshData.shaftV] = ...
        patchCleanUnused(meshData.shaftF, meshData.shaftV);
    
    %Check if shaft logic vertices requires correction due to mismatch in
    %number of points between shaft vertices. This happens with some
    %specific cases and needs correcting
    if length(meshData.shaftV) ~= sum(meshData.shaftLogicVertices)
        %Recreate the altered vertices used to extract shaft vertices
        testOuterCorticalV = meshData.outerCorticalV;
        testOuterCorticalV(meshData.indBoundaryTop,2) = meshData.maxShaftVal;
        testOuterCorticalV(meshData.indBoundaryBottom,2) = meshData.minShaftVal;
        %Create new logical vertices variable
        newLogicVertices = logical(zeros(length(testOuterCorticalV),1));
        %Properly identify points that logically match between shaft
        %vertices and original outer cortical surface
        for ptInd = 1:length(newLogicVertices)
            ptDist = distancePoints3d(testOuterCorticalV(ptInd,:), meshData.shaftV);
            if min(ptDist) == 0
                newLogicVertices(ptInd,1) = true;
            else
                newLogicVertices(ptInd,1) = false;
            end
        end
        %Replace the shaft logic vertices variable
        meshData.shaftLogicVertices = newLogicVertices;        
    end
    
% % %     %Visualise final shaft in the context of the entire tibia
% % %     cFigure; hold on
% % %     gpatch(meshData.outerCorticalF, meshData.outerCorticalV, 'rw', 'none', 0.3);
% % %     gpatch(meshData.shaftF, meshData.shaftV, 'gw', 'k');
% % %     axisGeom; camlight headlight

end