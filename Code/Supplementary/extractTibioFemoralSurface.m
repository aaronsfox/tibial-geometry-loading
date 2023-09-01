function [tibiaProxSurfacePts, tibiaDistSurfacePts, ankleCentrePt] = extractTibioFemoralSurface(tibiaF, tibiaV)

    %% This script uses the GIBOC algorithm to extract the articulating surface
    %  at the proximal end of the tibia. The nodes are subsequently used as
    %  a fixed boundary condition in simulations to fix the proximal end of
    %  the bone.
    
    %  Inputs:
    %       tibiaV - n x 3 array of the XYZ points of the outer tibia surface 
    %       tibiaF - connectivity list of faces for the tibia
    
    %% Set-up
    
    %Extract proximal and distal ends of the tibia and convert to
    %triangulations as per GIBOC function needs
    
    %Identify y-level of upper and lower thirds/fifths
    upperThirdPt = max(tibiaV(:,2)) * (2/3);
    lowerFifthPt = max(tibiaV(:,2)) * (1/5);
    
    %% Extract proximal tibia
       
    %Create a logic for cutting away faces
    logicVerticesProx = tibiaV(:,2) > upperThirdPt;
    logicFacesProx = all(logicVerticesProx(tibiaF),2);
    logicFacesProx = triSurfLogicSharpFix(tibiaF, logicFacesProx, 3); %Altered logic so it is smoother
    
    %Extract the faces to keep
    tibiaProxF = tibiaF(logicFacesProx, :); %The faces to keep
    [tibiaProxF, tibiaProxV] = patchCleanUnused(tibiaProxF, tibiaV); %Remove unused points
    
    %Attempt to self triangulate potentially jagged edge
    proxEb = patchBoundary(tibiaProxF, tibiaProxV); %Get boundary edges
    indBoundaryProx = edgeListToCurve(proxEb); %Convert boundary edges to a curve list
    indBoundaryProx = indBoundaryProx(1:end-1); %Trim off last point since it is equal to first on a closed loop
    angleThreshold = pi*(120/180); %threshold for self triangulation
    [tibiaProxF, tibiaProxV, indBoundaryTop] = ...
        triSurfSelfTriangulateBoundary(tibiaProxF, tibiaProxV, indBoundaryProx, angleThreshold, 1);
    
    %Force boundary to have the Y level chosen
    tibiaProxV(indBoundaryTop, 2) = upperThirdPt;
    
    %Close over the end hole
    [proxCloseF, proxCloseV] = ...
        regionTriMesh2D({tibiaProxV(indBoundaryTop, [1 3])}, 0.5, 0, 0);
    proxCloseV = [proxCloseV(:,1), ...
        ones(length(proxCloseV),1) * mean(tibiaProxV(indBoundaryTop, 2)), ... %Add/set y-level
        proxCloseV(:,2)];
    
    %Join surface sets
    [tibiaProxF, tibiaProxV] = joinElementSets({tibiaProxF, proxCloseF}, ...
        {tibiaProxV, proxCloseV});
    
    %Remesh for consistent point spacing
    optionStruct.nb_pts = 20000; %Set desired number of points
    optionStruct.disp_on = 0; % Turn off command window text display
    [tibiaProxF, tibiaProxV] = ggremesh(tibiaProxF, tibiaProxV, optionStruct);
    
    %% Extract distal tibia

    %Create a logic for cutting away faces
    logicVerticesDist = tibiaV(:,2) < lowerFifthPt;
    logicFacesDist = all(logicVerticesDist(tibiaF),2);
    logicFacesDist = triSurfLogicSharpFix(tibiaF, logicFacesDist, 3); %Altered logic so it is smoother
    
    %Extract the faces to keep
    tibiaDistF = tibiaF(logicFacesDist, :); %The faces to keep
    [tibiaDistF, tibiaDistV] = patchCleanUnused(tibiaDistF, tibiaV); %Remove unused points
    
    %Attempt to self triangulate potentially jagged edge
    distEb = patchBoundary(tibiaDistF, tibiaDistV); %Get boundary edges
    indBoundaryDist = edgeListToCurve(distEb); %Convert boundary edges to a curve list
    indBoundaryDist = indBoundaryDist(1:end-1); %Trim off last point since it is equal to first on a closed loop
    angleThreshold = pi*(120/180); %threshold for self triangulation
    [tibiaDistF, tibiaDistV, indBoundaryTop] = ...
        triSurfSelfTriangulateBoundary(tibiaDistF, tibiaDistV, indBoundaryDist, angleThreshold, 1);
    
    %Force boundary to have the Y level chosen
    tibiaDistV(indBoundaryTop, 2) = lowerFifthPt;
    
    %Close over the end hole
    [distCloseF, distCloseV] = ...
        regionTriMesh2D({tibiaDistV(indBoundaryTop, [1 3])}, 0.5, 0, 0);
    distCloseV = [distCloseV(:,1), ...
        ones(length(distCloseV),1) * mean(tibiaDistV(indBoundaryTop, 2)), ... %Add/set y-level
        distCloseV(:,2)];
    
    %Join surface sets
    [tibiaDistF, tibiaDistV] = joinElementSets({tibiaDistF, distCloseF}, ...
        {tibiaDistV, distCloseV});
    
    %Remesh for consistent point spacing
    optionStruct.nb_pts = 20000; %Set desired number of points
    optionStruct.disp_on = 0; % Turn off command window text display
    [tibiaDistF, tibiaDistV] = ggremesh(tibiaDistF, tibiaDistV, optionStruct);
    
    %% Reorient surfaces to match algorithm
    
% % %     ours = algorithm
% % %     -X = +Y
% % %     +Z = -X
% % %     +Y = +Z
    
    %Proximal surface
    tibiaProxV_rot = [tibiaProxV(:,3)*-1, ...
        tibiaProxV(:,1)*-1, ...
        tibiaProxV(:,2)];
    
    %Distal surface
    tibiaDistV_rot = [tibiaDistV(:,3)*-1, ...
        tibiaDistV(:,1)*-1, ...
        tibiaDistV(:,2)];

    %% Identify articulating surfaces
    
    %Convert faces and vertices to triangulation objects
    proxTib = triangulation(tibiaProxF, tibiaProxV_rot);
    distTib = triangulation(tibiaDistF, tibiaDistV_rot);
    
    %Run function for joint surfaces
    [ankleAS, ankleCentre, tibAS] = RTibiaFun_modified(proxTib, distTib);
    
% % %     %Visualise
% % %     cFigure; hold on;
% % %     gpatch(tibiaProxF, tibiaProxV_rot, 'kw', 'none', 0.8);
% % %     gpatch(tibiaDistF, tibiaDistV_rot, 'kw', 'none', 0.8);
% % %     axisGeom; camlight headlight
% % %     plotV(ankleAS.Points, 'b.');
% % %     plotV(tibAS.Points, 'r.');

    %% Reduce point number from articulating surfaces to match originals
    
    %Rotate the articulating surface points back to original coordinate system
    proxSurfacePoints = [tibAS.Points(:,2)*-1, ...
        tibAS.Points(:,3), ...
        tibAS.Points(:,1)*-1];
    distSurfacePoints = [ankleAS.Points(:,2)*-1, ...
        ankleAS.Points(:,3), ...
        ankleAS.Points(:,1)*-1];
    ankleCentre = [ankleCentre(:,2)*-1, ...
        ankleCentre(:,3), ...
        ankleCentre(:,1)*-1];
    
    %Find the nearest point in the proximal and distal surfaces to the
    %identified articulating surfaces
    
    %Proximal tibia
    %Find nearest neighbours to surface points
    proxNeighbours = knnsearch(tibiaV, proxSurfacePoints);
    %Keep the unique neighbours identified
    proxNeighboursUnique = unique(proxNeighbours);
    %Extract the points
    tibiaProxSurfacePts = tibiaV(proxNeighboursUnique,:);
    
    %Distal tibia
    %Find nearest neighbours to surface points
    distNeighbours = knnsearch(tibiaV, distSurfacePoints);
    %Keep the unique neighbours identified
    distNeighboursUnique = unique(distNeighbours);
    %Extract the points
    tibiaDistSurfacePts = tibiaV(distNeighboursUnique,:);
    
    %Ankle centre
    %Find nearest neighbours to surface points
    distNeighbours = knnsearch(tibiaV, ankleCentre);
    %Keep the unique neighbours identified
    distNeighboursUnique = unique(distNeighbours);
    %Extract the points
    ankleCentrePt = tibiaV(distNeighboursUnique,:);
    
% % %     %Visualise
% % %     cFigure; hold on;
% % %     gpatch(tibiaF, tibiaV, 'kw', 'none', 0.8);
% % %     axisGeom; camlight headlight
% % %     plotV(tibiaProxSurfacePts, 'r.');
% % %     plotV(tibiaDistSurfacePts, 'r.');
% % %     plotV(ankleCentrePt, 'b.', 'MarkerSize', 20);

end