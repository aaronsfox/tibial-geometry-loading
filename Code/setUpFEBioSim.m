function [febioFebFileName, febioLogFileName] = setUpFEBioSim(runTrial, caseID, meshOutputFullTibia, meshOutputFibula, originalSurfaces, simLandmarks, simForces, settingsFEA, matParameters)

    %% This script generates the FEBio model to run a simulation
    %
    %  Inputs:
    %       runTrial - label (Run_3, Run_4, Run_5) for running trial to simulate
    %       caseID - case number for the current participant
    %       meshOutputFullTibia - volumetric mesh of the tibia (inc. trabecular) generated by TetGen
    %       meshOutputFibula - volumetric mesh of the fibula generated by TetGen
    %       originalSurfaces - structure with original vertices and faces data [not really used]
    %       simLandmarks - a series of reference points used in creating node sets
    %       simForces - a structure containing the forces to apply in simulation
    %       settingsFEA - settings to input to the FEBio model
    %       matParameters - material parameters for the various parts in model 
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
    
    %% Set-up the FEBio template

    %Get FEBio template with default settings
    [febio_spec] = febioStructTemplate;

    %febio_spec version
    febio_spec.ATTR.version = '3.0';

    %Module section
    febio_spec.Module.ATTR.type = 'solid';

    %Control section
    febio_spec.Control.analysis = 'STATIC';
    febio_spec.Control.time_steps = settingsFEA.numTimeSteps;
    febio_spec.Control.step_size = 1/settingsFEA.numTimeSteps;
    febio_spec.Control.solver.max_refs = settingsFEA.max_refs;
    febio_spec.Control.solver.max_ups = settingsFEA.max_ups;
    febio_spec.Control.time_stepper.dtmin = settingsFEA.dtmin;
    febio_spec.Control.time_stepper.dtmax = settingsFEA.dtmax;
    febio_spec.Control.time_stepper.max_retries = settingsFEA.max_retries;
    febio_spec.Control.time_stepper.opt_iter = settingsFEA.opt_iter;

    %% Materials section
       
    %Material parameters (MPa if spatial units are mm)
    %In an ideal setting we'd estimate these from images as per Haider et al.
    %(2020). Here we take an estimate of axial elastic bone modulus for
    %cortical and trabecular as 18600MPa and 10400MPa, respectively, from
    %Edwards et al. (2010). Our coordinate system is in the format of X
    %being anterior-posterior, Y being longitudinal/axial, and Z being
    %medial-lateral. Therefore our definitions of the material parameters
    %are slightly different with respect to the axes compared to Haider et
    %al. (2020). Section 4.1.3.13 of the FEBio User Manual (v3.0) has some
    %details on the definition of orthotropic elastic materials.
    
    %Set the material parameters to enter into orthotropic elastic material
    %Cortical bone (inc. whole fibula)
    corticalE2 = matParameters.corticalYoungs; %longitudinal Young's modulus
    corticalE1 = 0.577 * corticalE2; %anterior-posterior Young's modulus
    corticalE3 = 0.574 * corticalE2; %medial-lateral Young's modulus modulus
    corticalG12 = 0.265 * corticalE2; %anterior-posterior & axial shear modulus
    corticalG23 = 0.216 * corticalE2; %axial & medial-lateral shear modulus
    corticalG31 = 0.195 * corticalE2; %anterior-posterior & medial-lateral shear modulus
    corticalV12 = 0.234; %Poisson's ratio anterior-posterior & axial
    corticalV23 = 0.405; %Poisson's ratio axial & medial-lateral
    corticalV31 = 0.427; %Poisson's ratio medial-lateral & anterior-posterior
    %Trabecular bone
    trabecularE2 = matParameters.trabYoungs; %longitudinal Young's modulus
    trabecularE1 = 0.577 * trabecularE2; %anterior-posterior Young's modulus
    trabecularE3 = 0.574 * trabecularE2; %medial-lateral Young's modulus modulus
    trabecularG12 = 0.265 * trabecularE2; %anterior-posterior & axial shear modulus
    trabecularG23 = 0.216 * trabecularE2; %axial & medial-lateral shear modulus
    trabecularG31 = 0.195 * trabecularE2; %anterior-posterior & medial-lateral shear modulus
    trabecularV12 = 0.234; %Poisson's ratio anterior-posterior & axial
    trabecularV23 = 0.405; %Poisson's ratio axial & medial-lateral
    trabecularV31 = 0.427; %Poisson's ratio medial-lateral & anterior-posterior
    
    %Cortical tibia properties
    febio_spec.Material.material{1}.ATTR.name = 'corticalTibia';
    febio_spec.Material.material{1}.ATTR.type = 'orthotropic elastic';
    febio_spec.Material.material{1}.ATTR.id = 1;
    febio_spec.Material.material{1}.mat_axis.ATTR.type = 'vector';
    febio_spec.Material.material{1}.mat_axis.a = [1,0,0];
    febio_spec.Material.material{1}.mat_axis.d = [0,1,0];
    febio_spec.Material.material{1}.E1 = corticalE1;
    febio_spec.Material.material{1}.E2 = corticalE2;
    febio_spec.Material.material{1}.E3 = corticalE3;
    febio_spec.Material.material{1}.v12 = corticalV12;
    febio_spec.Material.material{1}.v23 = corticalV23;
    febio_spec.Material.material{1}.v31 = corticalV31;
    febio_spec.Material.material{1}.G12 = corticalG12;
    febio_spec.Material.material{1}.G23 = corticalG23;
    febio_spec.Material.material{1}.G31 = corticalG31;
    
    %Trabecular tibia properties
    febio_spec.Material.material{2}.ATTR.name = 'trabTibia';
    febio_spec.Material.material{2}.ATTR.type = 'orthotropic elastic';
    febio_spec.Material.material{2}.ATTR.id = 2;
    febio_spec.Material.material{2}.mat_axis.ATTR.type = 'vector';
    febio_spec.Material.material{2}.mat_axis.a = [1,0,0];
    febio_spec.Material.material{2}.mat_axis.d = [0,1,0];
    febio_spec.Material.material{2}.E1 = trabecularE1;
    febio_spec.Material.material{2}.E2 = trabecularE2;
    febio_spec.Material.material{2}.E3 = trabecularE3;
    febio_spec.Material.material{2}.v12 = trabecularV12;
    febio_spec.Material.material{2}.v23 = trabecularV23;
    febio_spec.Material.material{2}.v31 = trabecularV31;
    febio_spec.Material.material{2}.G12 = trabecularG12;
    febio_spec.Material.material{2}.G23 = trabecularG23;
    febio_spec.Material.material{2}.G31 = trabecularG31;
    
    %Fibula properties (assumed cortical)
    febio_spec.Material.material{3}.ATTR.name = 'corticalFibula';
    febio_spec.Material.material{3}.ATTR.type = 'orthotropic elastic';
    febio_spec.Material.material{3}.ATTR.id = 3;
    febio_spec.Material.material{3}.mat_axis.ATTR.type = 'vector';
    febio_spec.Material.material{3}.mat_axis.a = [1,0,0];
    febio_spec.Material.material{3}.mat_axis.d = [0,1,0];
    febio_spec.Material.material{3}.E1 = corticalE1;
    febio_spec.Material.material{3}.E2 = corticalE2;
    febio_spec.Material.material{3}.E3 = corticalE3;
    febio_spec.Material.material{3}.v12 = corticalV12;
    febio_spec.Material.material{3}.v23 = corticalV23;
    febio_spec.Material.material{3}.v31 = corticalV31;
    febio_spec.Material.material{3}.G12 = corticalG12;
    febio_spec.Material.material{3}.G23 = corticalG23;
    febio_spec.Material.material{3}.G31 = corticalG31;

    %% Mesh section
    
	%Access mesh data from tibia output
    fullTibiaE = meshOutputFullTibia.elements; %The elements
    fullTibiaV = meshOutputFullTibia.nodes; %The vertices or nodes
    fullTibiaCE = meshOutputFullTibia.elementMaterialID; %Element material or region id
    fullTibiaFb = meshOutputFullTibia.facesBoundary; %The boundary faces
    fullTibiaCb = meshOutputFullTibia.boundaryMarker; %The boundary markers

    %Extract the element matrix components
    %Cortical elements = -3; Trabecular elements = -2;
    corticalTibiaE = fullTibiaE(find(fullTibiaCE == -3),:);
    trabTibiaE = fullTibiaE(find(fullTibiaCE == -2),:);
    
    %Extract just tibial surface nodes for later calculations
    [tibiaF, tibiaV] = patchCleanUnused(fullTibiaFb(fullTibiaCb==1,:), fullTibiaV);

    %Access elements, nodes, and boundary faces
    fibulaMeshE = meshOutputFibula.elements;
    fibulaMeshV = meshOutputFibula.nodes;
    fibulaMeshFb = meshOutputFibula.facesBoundary;
    fibulaMeshCb = meshOutputFibula.boundaryMarker;
    fibulaMeshCE = meshOutputFibula.elementMaterialID;

    %Extract fibula elements and add the tibial node length so it works
    %appropriately with the full node set
    fibulaE = size(fullTibiaV,1) + fibulaMeshE;

    %Nodes
    
    %Combine nodes of both bones into one structure
    V = [fullTibiaV; fibulaMeshV];
    febio_spec.Mesh.Nodes{1}.ATTR.name = 'allNodes'; %The node set name
    febio_spec.Mesh.Nodes{1}.node.ATTR.id = (1:size(V,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL = V; %The nodel coordinates
    
    %Surfaces
    
    %Tibia surface
    febio_spec.Mesh.Surface{1}.ATTR.name = 'tibiaSurface';
    febio_spec.Mesh.Surface{1}.tri3.ATTR.id = (1:1:size(fullTibiaFb,1))';
    febio_spec.Mesh.Surface{1}.tri3.VAL = fullTibiaFb;
    %Fibula surface
    febio_spec.Mesh.Surface{2}.ATTR.name = 'fibulaSurface';
    febio_spec.Mesh.Surface{2}.tri3.ATTR.id = (1:1:size(fibulaMeshFb,1))';
    febio_spec.Mesh.Surface{2}.tri3.VAL = size(fullTibiaV,1) + fibulaMeshFb;
    
    %Surface pairs
    febio_spec.Mesh.SurfacePair{1}.ATTR.name = 'tibiaFibulaSurfacePair';
    febio_spec.Mesh.SurfacePair{1}.primary = 'tibiaSurface';
    febio_spec.Mesh.SurfacePair{1}.secondary = 'fibulaSurface';
    
    %% Contact section
    
    %Tied contact between tibia and fibula
    febio_spec.Contact.contact{1}.ATTR.type = 'tied-elastic';
    febio_spec.Contact.contact{1}.ATTR.surface_pair = 'tibiaFibulaSurfacePair';

    %% Node Sets section

    %Tibia
    febio_spec.Mesh.NodeSet{1}.ATTR.name = 'tibiaNodes';
    febio_spec.Mesh.NodeSet{1}.node.ATTR.id = (1:size(fullTibiaV,1))';

    %Fibula
    febio_spec.Mesh.NodeSet{2}.ATTR.name = 'fibulaNodes';
    febio_spec.Mesh.NodeSet{2}.node.ATTR.id = size(fullTibiaV,1) + (1:size(fibulaMeshV,1))';

    %Ankle joint contact node
    febio_spec.Mesh.NodeSet{3}.ATTR.name = 'ankleJointCentreNode';
    febio_spec.Mesh.NodeSet{3}.node.ATTR.id = find(distancePoints3d(fullTibiaV, simLandmarks.contactPointAJC) == ...
        min(distancePoints3d(fullTibiaV, simLandmarks.contactPointAJC)), 1);

    %Ankle joint contact nodes
    %These are the surface points estimated across the whole ankle surface
    febio_spec.Mesh.NodeSet{4}.ATTR.name = 'ankleJointContactNode';
    for ptInd = 1:length(simLandmarks.tibiaDistSurfacePts)
    febio_spec.Mesh.NodeSet{4}.node.ATTR.id(ptInd,1) = find(distancePoints3d(fullTibiaV, simLandmarks.tibiaDistSurfacePts(ptInd,:)) == ...
        min(distancePoints3d(fullTibiaV, simLandmarks.tibiaDistSurfacePts(ptInd,:))), 1);
    end

    %Proximal tibial surface
    febio_spec.Mesh.NodeSet{5}.ATTR.name = 'proximalSurfaceNodes';
    for ptInd = 1:length(simLandmarks.tibiaProxSurfacePts)
    febio_spec.Mesh.NodeSet{5}.node.ATTR.id(ptInd,1) = find(distancePoints3d(fullTibiaV, simLandmarks.tibiaProxSurfacePts(ptInd,:)) == ...
        min(distancePoints3d(fullTibiaV, simLandmarks.tibiaProxSurfacePts(ptInd,:))), 1);
    end

    %Lateral malleolus node
    febio_spec.Mesh.NodeSet{6}.ATTR.name = 'lateralMalleolus';
    febio_spec.Mesh.NodeSet{6}.node.ATTR.id = find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.landmarksInd.LM,:)) == ...
        min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.landmarksInd.LM,:))), 1);

    %Medial malleolus node
    febio_spec.Mesh.NodeSet{7}.ATTR.name = 'medialMalleolus';
    febio_spec.Mesh.NodeSet{7}.node.ATTR.id = find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.landmarksInd.MM,:)) == ...
        min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.landmarksInd.MM,:))), 1);

    %Muscle force application points
    for muscleInd = 1:length(muscleList)
        %Get muscle attachment point coordinates
        muscleCoord = simLandmarks.tibiaFibulaV(simLandmarks.muscleAttachmentInd.(muscleList{muscleInd}),:);
        %Specify node set details
        nodeSetInd = size(febio_spec.Mesh.NodeSet,2) + 1;
        febio_spec.Mesh.NodeSet{nodeSetInd}.ATTR.name = [muscleList{muscleInd},'_Node'];
        febio_spec.Mesh.NodeSet{nodeSetInd}.node.ATTR.id = find(distancePoints3d(fullTibiaV, muscleCoord) == ...
            min(distancePoints3d(fullTibiaV, muscleCoord)), 1);    
    end
    
    %% Discrete Sets section
    %For connecting springs between tibia-fibula nodes

    %Anterior-proximal
    febio_spec.Mesh.DiscreteSet{1}.ATTR.name = 'AntProxSet';
    febio_spec.Mesh.DiscreteSet{1}.delem.VAL = [...
        find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.AntProxFib,:)) == min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.AntProxFib,:)))), ...
        find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.AntProxTib,:)) == min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.AntProxTib,:))))];

    %Posterior-proximal
    febio_spec.Mesh.DiscreteSet{2}.ATTR.name = 'PostProxSet';
    febio_spec.Mesh.DiscreteSet{2}.delem.VAL = [...
        find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.PostProxFib,:)) == min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.PostProxFib,:)))), ...
        find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.PostProxTib,:)) == min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.PostProxTib,:))))];

    %Anterior-distal
    febio_spec.Mesh.DiscreteSet{3}.ATTR.name = 'AntDistSet';
    febio_spec.Mesh.DiscreteSet{3}.delem.VAL = [...
        find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.AntDistFib,:)) == min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.AntDistFib,:)))), ...
        find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.AntDistTib,:)) == min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.AntDistTib,:))))];

    %Posterior-distal
    febio_spec.Mesh.DiscreteSet{4}.ATTR.name = 'PostDistSet';
    febio_spec.Mesh.DiscreteSet{4}.delem.VAL = [...
        find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.PostDistFib,:)) == min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.PostDistFib,:)))), ...
        find(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.PostDistTib,:)) == min(distancePoints3d(V, simLandmarks.tibiaFibulaV(simLandmarks.interactionsInd.PostDistTib,:))))];

    %% Elements section

    %Cortical bone
    febio_spec.Mesh.Elements{1}.ATTR.name = 'corticalTibiaPart'; %Name of this part
    febio_spec.Mesh.Elements{1}.ATTR.type = 'tet4'; %Element type
    febio_spec.Mesh.Elements{1}.elem.ATTR.id = (1:1:size(corticalTibiaE,1))'; %Element id's
    febio_spec.Mesh.Elements{1}.elem.VAL = corticalTibiaE; %The element matrix

    %Trabecular bone
    febio_spec.Mesh.Elements{2}.ATTR.name = 'trabTibiaPart'; %Name of this part
    febio_spec.Mesh.Elements{2}.ATTR.type = 'tet4'; %Element type
    febio_spec.Mesh.Elements{2}.elem.ATTR.id = size(corticalTibiaE,1) + (1:1:size(trabTibiaE,1))'; %Element id's
    febio_spec.Mesh.Elements{2}.elem.VAL = trabTibiaE; %The element matrix

    %Fibula
    febio_spec.Mesh.Elements{3}.ATTR.name = 'fibulaPart'; %Name of this part
    febio_spec.Mesh.Elements{3}.ATTR.type = 'tet4'; %Element type
    febio_spec.Mesh.Elements{3}.elem.ATTR.id = size(corticalTibiaE,1) + size(trabTibiaE,1) + (1:1:size(fibulaE,1))'; %Element id's
    febio_spec.Mesh.Elements{3}.elem.VAL = fibulaE; %The element matrix

    %% Mesh Domains section

    %Cortical tibia
    febio_spec.MeshDomains.SolidDomain{1}.ATTR.name = 'corticalTibiaPart';
    febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat = 'corticalTibia';

    %Trabecular tibia
    febio_spec.MeshDomains.SolidDomain{2}.ATTR.name = 'trabTibiaPart';
    febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat = 'trabTibia';

    %Fibula
    febio_spec.MeshDomains.SolidDomain{3}.ATTR.name = 'fibulaPart';
    febio_spec.MeshDomains.SolidDomain{3}.ATTR.mat = 'corticalFibula';

    %% Boundary Conditions section

    %Fixed proximal end
    febio_spec.Boundary.bc{1}.ATTR.type = 'fix';
    febio_spec.Boundary.bc{1}.ATTR.node_set = 'proximalSurfaceNodes';
    febio_spec.Boundary.bc{1}.dofs = 'x,y,z';

    % % % %Fixed distal end
    % % % febio_spec.Boundary.bc{2}.ATTR.type = 'fix';
    % % % febio_spec.Boundary.bc{2}.ATTR.node_set = 'ankleJointCentreNode';
    % % % febio_spec.Boundary.bc{2}.dofs = 'x,z';

    %Fixed malleoli
    %Medial
    febio_spec.Boundary.bc{2}.ATTR.type = 'fix';
    febio_spec.Boundary.bc{2}.ATTR.node_set = 'medialMalleolus';
    febio_spec.Boundary.bc{2}.dofs = 'x,y,z';
    
    %% Mesh Data section

    %Ankle joint contact force - isolated point
    febio_spec.MeshData.NodeData{1}.ATTR.name = 'ankleJointContactForce';
    febio_spec.MeshData.NodeData{1}.ATTR.node_set = 'ankleJointContactNode';
    febio_spec.MeshData.NodeData{1}.ATTR.datatype = 'vec3';
    febio_spec.MeshData.NodeData{1}.node.ATTR.lid = (1:1:numel(simLandmarks.contactPointAJC,1))';
    febio_spec.MeshData.NodeData{1}.node.VAL = simForces.(char(runTrial)).jointContactForcePeak;

    %Muscle forces
    %Loop through for each muscle
    for muscleInd = 1:length(muscleList)
        %Increment size in mesh data structure
        nodeDataInd = size(febio_spec.MeshData.NodeData,2) + 1;
        %Specify load details
        febio_spec.MeshData.NodeData{nodeDataInd}.ATTR.name = [muscleList{muscleInd},'_Force'];
        febio_spec.MeshData.NodeData{nodeDataInd}.ATTR.node_set = [muscleList{muscleInd},'_Node'];
        febio_spec.MeshData.NodeData{nodeDataInd}.ATTR.datatype = 'vec3';
        febio_spec.MeshData.NodeData{nodeDataInd}.node.ATTR.lid = (1:1:numel(originalSurfaces.tibiaV(simLandmarks.muscleAttachmentInd.(muscleList{muscleInd}),:),1))';
        febio_spec.MeshData.NodeData{nodeDataInd}.node.VAL = simForces.(char(runTrial)).(muscleList{muscleInd}).muscleForceContactPeak;
    end
    
    %% Loads section

    %Ankle joint contact forces
    febio_spec.Loads.nodal_load{1}.ATTR.name = 'ankleJointContactNodalForce';
    febio_spec.Loads.nodal_load{1}.ATTR.type = 'nodal_force';
    febio_spec.Loads.nodal_load{1}.ATTR.node_set = 'ankleJointContactNode';
    febio_spec.Loads.nodal_load{1}.value.ATTR.lc = 1;
    febio_spec.Loads.nodal_load{1}.value.ATTR.type = 'map';
    febio_spec.Loads.nodal_load{1}.value.VAL = 'ankleJointContactForce';

    %Muscle forces
    %Loop through for each muscle
    for muscleInd = 1:length(muscleList)
        %Increament size in mesh data structure
        nodeDataInd = size(febio_spec.Loads.nodal_load,2) + 1;
        %Specify load details
        febio_spec.Loads.nodal_load{nodeDataInd}.ATTR.name = [muscleList{muscleInd},'_NodalForce'];
        febio_spec.Loads.nodal_load{nodeDataInd}.ATTR.type = 'nodal_force';
        febio_spec.Loads.nodal_load{nodeDataInd}.ATTR.node_set = [muscleList{muscleInd},'_Node'];
        febio_spec.Loads.nodal_load{nodeDataInd}.value.ATTR.lc = 1;
        febio_spec.Loads.nodal_load{nodeDataInd}.value.ATTR.type = 'map';
        febio_spec.Loads.nodal_load{nodeDataInd}.value.VAL = [muscleList{muscleInd},'_Force'];
    end

    %% Load Data section

    %Generic load controller for all forces
    febio_spec.LoadData.load_controller{1}.ATTR.id = 1;
    febio_spec.LoadData.load_controller{1}.ATTR.type = 'loadcurve';
    febio_spec.LoadData.load_controller{1}.interpolate = 'LINEAR';
    febio_spec.LoadData.load_controller{1}.points.point.VAL = [0 0; 1 1];

    %% Discrete section
    
    %NOTE: spring parameters based on Haider et al. (2020)

    %Anterior-proximal spring
    febio_spec.Discrete.discrete_material{1}.ATTR.id = 1;
    febio_spec.Discrete.discrete_material{1}.ATTR.type = 'linear spring';
    febio_spec.Discrete.discrete_material{1}.ATTR.name = 'AntProxSpring';
    febio_spec.Discrete.discrete_material{1}.E = 133;
    febio_spec.Discrete.discrete{1}.ATTR.dmat = 1;
    febio_spec.Discrete.discrete{1}.ATTR.discrete_set = 'AntProxSet';

    %Posterior-proximal spring
    febio_spec.Discrete.discrete_material{2}.ATTR.id = 2;
    febio_spec.Discrete.discrete_material{2}.ATTR.type = 'linear spring';
    febio_spec.Discrete.discrete_material{2}.ATTR.name = 'PostProxSpring';
    febio_spec.Discrete.discrete_material{2}.E = 109;
    febio_spec.Discrete.discrete{2}.ATTR.dmat = 2;
    febio_spec.Discrete.discrete{2}.ATTR.discrete_set = 'PostProxSet';

    %Anterior-distal spring
    febio_spec.Discrete.discrete_material{3}.ATTR.id = 3;
    febio_spec.Discrete.discrete_material{3}.ATTR.type = 'linear spring';
    febio_spec.Discrete.discrete_material{3}.ATTR.name = 'AntDistSpring';
    febio_spec.Discrete.discrete_material{3}.E = 78;
    febio_spec.Discrete.discrete{3}.ATTR.dmat = 3;
    febio_spec.Discrete.discrete{3}.ATTR.discrete_set = 'AntDistSet';

    %Posterior-distal spring
    febio_spec.Discrete.discrete_material{4}.ATTR.id = 4;
    febio_spec.Discrete.discrete_material{4}.ATTR.type = 'linear spring';
    febio_spec.Discrete.discrete_material{4}.ATTR.name = 'PostDistSpring';
    febio_spec.Discrete.discrete_material{4}.E = 101;
    febio_spec.Discrete.discrete{4}.ATTR.dmat = 4;
    febio_spec.Discrete.discrete{4}.ATTR.discrete_set = 'PostDistSet';

    %% Output section
    
    %Get current directory to reference back to
    homeDir = pwd;
    
    %Create directory to store results
    mkdir(['..\Simulations\case-',num2str(caseID),'\',runTrial]);
    
    %Navigate to results directory
    cd(['..\Simulations\case-',num2str(caseID),'\',runTrial]);

    %Defining file name
    febioFebFileNamePart = ['case-',num2str(caseID),'_',runTrial];
    
    %FEBio model file
    febioFebFileName = fullfile(pwd,[febioFebFileNamePart,'.feb']); %FEB file name

    %Log and output files
    febioLogFileName = [febioFebFileNamePart,'.txt']; %FEBio log file name
    febioLogFileName_disp = [febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
    febioLogFileName_stress = [febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stresses
    febioLogFileName_strain = [febioFebFileNamePart,'_strain_out.txt']; %Log file name for exporting strains
    febioLogFileName_effectiveStrain = [febioFebFileNamePart,'_effectiveStrain_out.txt']; %Log file name for exporting effective strains

    %Settings for collecting displacements to log file
    febio_spec.Output.logfile.ATTR.file = febioLogFileName;
    febio_spec.Output.logfile.node_data{1}.ATTR.file = febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
    febio_spec.Output.logfile.node_data{1}.VAL = 1:size(V,1);

    %Settings for collecting stress to log file
    febio_spec.Output.logfile.element_data{1}.ATTR.file = febioLogFileName_stress;
% % %     febio_spec.Output.logfile.element_data{1}.ATTR.data = 's1;s2;s3';
    febio_spec.Output.logfile.element_data{1}.ATTR.data = 'effective stress';
    febio_spec.Output.logfile.element_data{1}.ATTR.delim = ',';
    febio_spec.Output.logfile.element_data{1}.VAL = 1:size(corticalTibiaE,1);

    %Settings for collecting strain to log file
    febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_strain;
    febio_spec.Output.logfile.element_data{2}.ATTR.data='E1;E2;E3';
    febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';
    febio_spec.Output.logfile.element_data{2}.VAL = 1:size(corticalTibiaE,1);
    
    %Settings for collecting strain to log file
    febio_spec.Output.logfile.element_data{3}.ATTR.file=febioLogFileName_effectiveStrain;
    febio_spec.Output.logfile.element_data{3}.ATTR.data='effective strain';
    febio_spec.Output.logfile.element_data{3}.ATTR.delim=',';
    febio_spec.Output.logfile.element_data{3}.VAL = 1:size(corticalTibiaE,1);

    %Export the FEB model file
    febioStruct2xml(febio_spec,febioFebFileName);

    %% Finish-up
    
    %Return to starting directory
    cd(homeDir);

%% ----- end of setUpFEBioSim.m ----- %%