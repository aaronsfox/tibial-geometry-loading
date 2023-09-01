function [shaftData] = extractMeanTibialShaftSurface(tibiaF, tibiaV, minShaftVal, maxShaftVal)

    %% This script extracts the tibial shaft surface from the outer tibial
    %  surface and associated mesh inputs
    %
    %  Inputs:
    %       tibiaF - faces of tibia surface
    %       tibiaV - 3D points of tibia surface
    %       minShaftVal - value relating to minimum boundary of shaft level
    %       maxShaftVal - value relating to maximum boundary of shaft level
    %
    %  Outputs:
    %       shaftData - structure with data related to shaft surface
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
    shaftData.shaftLogicVerticesLower = tibiaV(:,2) < maxShaftVal;
    shaftData.shaftLogicFacesLower = all(shaftData.shaftLogicVerticesLower(tibiaF), 2);
    shaftData.shaftLogicFacesLower = triSurfLogicSharpFix(tibiaF, shaftData.shaftLogicFacesLower, 3); %Altered logic so it is smoother

    %Cut away the faces for the lower part of the shaft
    %Use a temporary throughout here to later manipulate the actual surface
    shaftF = tibiaF(shaftData.shaftLogicFacesLower,:);

    %Clean up unused points
    [shaftF, shaftV] = patchCleanUnused(shaftF, tibiaV);

    %Extract boundaries of surface
    shaftEbLower = patchBoundary(shaftF, shaftV);
    indBoundaryEbLower = edgeListToCurve(shaftEbLower); %Convert boundary edges to a curve list
    indBoundaryEbLower = indBoundaryEbLower(1:end-1); %Trim off last point since it is equal to first on a closed loop

    %Identify the indices of those around the top boundary in the original
    %shape model tibia surface
    shaftData.indBoundaryTop = logical(zeros(length(tibiaV),1));
    for indPt = 1:length(indBoundaryEbLower)
        ind = find(distancePoints3d(tibiaV, shaftV(indBoundaryEbLower(indPt),:)) == 0);
        shaftData.indBoundaryTop(ind,1) = true;
    end

    %Repeat the above process to get the points along the lower boundary
    shaftData.shaftLogicVerticesUpper = tibiaV(:,2) > minShaftVal;
    shaftData.shaftLogicFacesUpper = all(shaftData.shaftLogicVerticesUpper(tibiaF), 2);
    shaftData.shaftLogicFacesUpper = triSurfLogicSharpFix(tibiaF, shaftData.shaftLogicFacesUpper, 3); %Altered logic so it is smoother

    %Use a temporary throughout here to later manipulate the actual surface
    shaftF = tibiaF(shaftData.shaftLogicFacesUpper,:);

    %Clean up unused points
    [shaftF, shaftV] = patchCleanUnused(shaftF, tibiaV);

    %Extract boundaries of surface
    shaftEbUpper = patchBoundary(shaftF, shaftV);
    indBoundaryEbUpper = edgeListToCurve(shaftEbUpper); %Convert boundary edges to a curve list
    indBoundaryEbUpper = indBoundaryEbUpper(1:end-1); %Trim off last point since it is equal to first on a closed loop

    %Identify the indices of those around the bottom boundary in the original
    %shape model tibia surface
    shaftData.indBoundaryBottom = logical(zeros(length(tibiaV),1));
    for indPt = 1:length(indBoundaryEbUpper)
        ind = find(distancePoints3d(tibiaV, shaftV(indBoundaryEbUpper(indPt),:)) == 0);
        shaftData.indBoundaryBottom(ind,1) = true;
    end

    %Create a shaft vertices variable and shift the boundaries to the desired levels
    shaftData.shaftV = tibiaV;
    shaftData.shaftV(shaftData.indBoundaryTop,2) = maxShaftVal;
    shaftData.shaftV(shaftData.indBoundaryBottom,2) = minShaftVal;

    %Extract the shaft while keeping the boundaries at the desired levels
    %Create the logicals for extracting vertices and faces
    shaftData.shaftLogicVertices = shaftData.shaftV(:,2) <= maxShaftVal & ...
        shaftData.shaftV(:,2) >= minShaftVal;
    shaftData.shaftLogicFaces = all(shaftData.shaftLogicVertices(tibiaF), 2);

    %Cut away the faces for the shaft
    shaftData.shaftF = tibiaF(shaftData.shaftLogicFaces,:);

    %Clean up unused points
    [shaftData.shaftF, shaftData.shaftV] = ...
        patchCleanUnused(shaftData.shaftF, shaftData.shaftV);
    
% % %     %Visualise final shaft in the context of the entire tibia
% % %     cFigure; hold on
% % %     gpatch(tibiaF, tibiaV, 'rw', 'none', 0.3);
% % %     gpatch(shaftData.shaftF, shaftData.shaftV, 'gw', 'k');
% % %     axisGeom; camlight headlight

end