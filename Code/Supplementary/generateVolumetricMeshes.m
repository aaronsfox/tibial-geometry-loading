function [meshOutputFullTibia, meshOutputFibula] = generateVolumetricMeshes(tibiaF, tibiaV, fibulaF, fibulaV, trabF, trabV)

    %% This script generates and outputs the volumetric meshes of the structures
    %  to be used in the FE simulations. Note that the remeshing of
    %  surfaces within this function is necessary to avoid errors in the
    %  tetgen function, and that this remeshing alters the node order of
    %  the surfaces so that they no longer match the original shape model
    %  data structure. This will get fixed later when analysing the FE
    %  simulation data anyway.
    %
    %  Inputs:
    %       tibiaF - faces for the tibia surface
    %       tibiaV - vertices for the tibia surface
    %       fibulaF - faces for the fibula surface
    %       fibulaV - vertices for the fibula surface
    %       trabF - faces for the trabecular surface
    %       trabV - vertices for the trabecular surface
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

    %% Generate volumetric meshes
    
    %Set options for tibial remeshing
    %Note that this increases the mesh density to include a larger number
    %of elements in the volumetric mesh to a similar scale to Bruce et al.
    %(2022)
    optionStruct_tib.nb_pts = 7000; %Set desired number of points
    optionStruct_tib.disp_on = 0; % Turn off command window text display
    optionStruct_fib.nb_pts = 3000; %Set desired number of points
    optionStruct_fib.disp_on = 0; % Turn off command window text display
    
    %Tibia
    [tibiaF, tibiaV] = ggremesh(tibiaF, tibiaV, optionStruct_tib);
    [tibiaF, tibiaV] = mergeVertices(tibiaF, tibiaV);
    %Trabecular
    [trabF, trabV] = ggremesh(trabF, trabV, optionStruct_tib);
    [trabF, trabV] = mergeVertices(trabF, trabV);
    %Fibula
    [fibulaF, fibulaV] = ggremesh(fibulaF, fibulaV, optionStruct_fib);
    [fibulaF, fibulaV] = mergeVertices(fibulaF, fibulaV);
    
    %Join the tibia and trabecular surfaces
    [fullTibiaF, fullTibiaV, fullTibiaC] = joinElementSets({tibiaF trabF},{tibiaV trabV});

    %Merge vertices to avoid meshing issues
    [fullTibiaF, fullTibiaV] = mergeVertices(fullTibiaF, fullTibiaV);

    %Mesh the two surfaces
    
    %Set properties
    [V_region1] = getInnerPoint(fullTibiaF(fullTibiaC==2,:),fullTibiaV); %First interior point
    [V_region2] = getInnerPoint({tibiaF trabF},{tibiaV trabV}); %Second interior point
    V_regions = [V_region1; V_region2]; %Collect region points
    V_holes = []; %Define hole points
    [regionTetVolume1] = tetVolMeanEst(tibiaF, tibiaV); %Volume estimate for regular tets
    [regionTetVolume2] = tetVolMeanEst(trabF,trabV); %Volume estimate for regular tets
    regionTetVolumes = [regionTetVolume1 regionTetVolume2];
    stringOpt = '-pq1.2AaY'; %Tetgen options
    
    %Create tetgen input structure
    tetGenStructFullTibia.stringOpt = stringOpt; %Tetgen options
    tetGenStructFullTibia.Faces = fullTibiaF; %Boundary faces
    tetGenStructFullTibia.Nodes = fullTibiaV; %Nodes of boundary
    tetGenStructFullTibia.faceBoundaryMarker = fullTibiaC;
    tetGenStructFullTibia.regionPoints = V_regions; %Interior points for regions
    tetGenStructFullTibia.holePoints = V_holes; %Interior points for holes
    tetGenStructFullTibia.regionA = regionTetVolumes; %Desired tetrahedral volume for each region
    
    %Mesh model using tetgen
    [meshOutputFullTibia] = runTetGen(tetGenStructFullTibia);

    % % % %Visualise
    % % % meshView(meshOutputFullTibia);

    %Drop error message if meshing is unsuccessful
    if isempty(meshOutputFullTibia.elements)
        error('Volumetric meshing of tibia unsuccessful. Check output. Intersecting faces is a common issue...');
    end
    
    %Mesh the fibula
    innerPointFibula = getInnerPoint(fibulaF, fibulaV);

    %Set mesh parameters
    tetVolumeFibula = tetVolMeanEst(fibulaF, fibulaV);
    tetGenStructFibula.stringOpt = '-pq1.2AaY';
    tetGenStructFibula.Faces = fibulaF;
    tetGenStructFibula.Nodes = fibulaV;
    tetGenStructFibula.holePoints = [];
    tetGenStructFibula.regionPoints = innerPointFibula;
    tetGenStructFibula.regionA = tetVolumeFibula;

    %Run tetgen
    [meshOutputFibula] = runTetGen(tetGenStructFibula);

    % % % %Visualise
    % % % meshView(meshOutputFibula);

    %Drop error message if meshing is unsuccessful
    if isempty(meshOutputFibula.elements)
        error('Volumetric meshing of fibula unsuccessful. Check output. Intersecting faces is a common issue...');
    end

end