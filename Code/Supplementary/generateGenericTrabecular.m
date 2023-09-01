function [F, V] = generateGenericTrabecular(tibiaV, shapeModel, corticalThickness, nb_pts)

    %% This script generates a generic internal trabecular surface of the 
    %  tibia based on a given 'cortical' thickness.
    
    %  Inputs:
    %       tibiaV - n x 3 array of the XYZ points of the outer tibia surface 
    %       shapeModel - the structure with shape model data for the tibia
    %       corticalThickness - desired outer 'cortical' thickness (in mm, larger value = skinnier internal surface)
	%       nb_pts - desired number of points to remesh the created surface to
    
    %% Set-up
    
    %Set ggremesh options structure
    opts.nb_pts = nb_pts; %Set desired number of points
    opts.disp_on = 0; % Turn off command window text display
    
    %% Create internal surface
    
% % %     %Set volume factor for meshing
% % %     volumeFactor = 2;
    
    %Remesh surface
    [tibiaF, tibiaV] = ggremesh(shapeModel.F, tibiaV, opts);
    [tibiaF, tibiaV] = mergeVertices(tibiaF, tibiaV);
    
    %Find interior point
    innerPoint = getInnerPoint(tibiaF, tibiaV);

    %Create a volumetric mesh using tetGen
    tetVolume = tetVolMeanEst(tibiaF, tibiaV); %Volume for regular tets
    tetGenStruct.stringOpt = '-pq1.2AaY';
    tetGenStruct.Faces = tibiaF;
    tetGenStruct.Nodes = tibiaV;
    tetGenStruct.holePoints = [];
    tetGenStruct.faceBoundaryMarker = ones(size(tibiaF,1),1); %Face boundary markers
    tetGenStruct.regionPoints = innerPoint; %region points
    tetGenStruct.regionA = tetVolume;

    %Run tetgen
    [meshOutput] = runTetGen(tetGenStruct); %Run tetGen

    %Access elements, nodes, and boundary faces
    E = meshOutput.elements;
    V = meshOutput.nodes;
    Fb = meshOutput.facesBoundary;
    Cb = meshOutput.boundaryMarker;
    CE = meshOutput.elementMaterialID;

    %Define a generic trabecular area
    indBoundary = unique(Fb(Cb == 1,:));
    DE = minDist(V,V(indBoundary,:));
    logicCorticalNodes = DE <= corticalThickness;
    logicCorticalElements = any(logicCorticalNodes(E),2);
    logicTrabecularElements =~ logicCorticalElements;

    %Extract separate elements
    E1 = E(logicCorticalElements,:);
    E2 = E(logicTrabecularElements,:);

    %Get face boundary for generic trabecular
    [genTrabF, ~] = meshBoundary(E2,[]);
    
    %Remesh new trabecular surface
    [F, V] = ggremesh(genTrabF, V, opts);
    
end