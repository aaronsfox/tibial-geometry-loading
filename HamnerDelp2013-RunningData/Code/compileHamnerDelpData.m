%% This script compiles the data from the Hamner & Delp (2013) dataset.
%  Specifically it extracts the individual subject data to create a
%  'generic' mean of the calculated results.
%
%  See the README.MD in the above folder for more descriptive details on this
%  code.
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

%Set home directory
homeDir = pwd;

%Import OpenSim libraries
import org.opensim.modeling.*


%Set list of subjects to work through
subList = [{'subject01'};
    {'subject02'};
    {'subject03'};
    {'subject04'};
    {'subject08'};
    {'subject10'};
    {'subject11'};
    {'subject17'};
    {'subject19'};
    {'subject20'}];

%Set running speed labels to work through
runLabels = [{'Run_3'};
    {'Run_4'};
    {'Run_5'}];

%Set labels for gait cycles in dataset
cycleLabels = [{'cycle1'};
    {'cycle2'};
    {'cycle3'}];

% % % %Set muscles that contribute to ankle joint reaction force
% % % ankleMuscles_r = [{'ext_dig_r'};
% % %     {'ext_hal_r'};
% % %     {'flex_dig_r'};
% % %     {'flex_hal_r'};
% % %     {'lat_gas_r'};
% % %     {'med_gas_r'};
% % %     {'per_brev_r'};
% % %     {'per_long_r'};
% % %     {'per_tert_r'};
% % %     {'soleus_r'};
% % %     {'tib_ant_r'};
% % %     {'tib_post_r'}];

%% Compile subject data

%Navigate to data directory
cd('..\Data\');

%Loop through subjects
for subInd = 1:length(subList)

    %Navigate to subjects data directory
    cd(subList{subInd});
    
    %Loop through running trials
    for runInd = 1:length(runLabels)
        
        %Navigate to folder
        cd(runLabels{runInd});

        %Loop through cycles
        for cycleInd = 1:length(cycleLabels)
            
            %Navigate to folder
            cd(cycleLabels{cycleInd});
            
            %%%%% MODEL DATA %%%%%
            
            %Load the model for current cycle
            trialModel = Model([runLabels{runInd},'_',cycleLabels{cycleInd},'_model.osim']);
            
            %Get the total mass of the trial model
            trialMass = 0;
            for bodyInd = 0:trialModel.updBodySet().getSize()-1
                trialMass = trialMass + trialModel.updBodySet().get(bodyInd).getMass();
            end
            
            %%%%% FORCES DATA %%%%%
            
            %Load in forces file
            forcesData = TimeSeriesTable([runLabels{runInd},'_',cycleLabels{cycleInd},'_forces.sto']);

            %Convert data to structure
            forcesStruct = osimTableToStruct(forcesData);
            
            %Get fieldnames from structure
            forcesNames = fieldnames(forcesStruct);
            
            %Keep right hand side labels
            forcesNames = forcesNames(endsWith(forcesNames,'_r'));
            
            %Create a 101 point normalised version of the time array
            normTime = linspace(forcesStruct.time(1), forcesStruct.time(end), 101);
            
            %Loop through labels, extract and normalise data
            for forcesInd = 1:length(forcesNames)                
                %Extract and interpolate the current variable data
                subjectData.(subList{subInd}).(runLabels{runInd}).forces.(forcesNames{forcesInd})(cycleInd,:) = ...
                    interp1(forcesStruct.time, forcesStruct.(forcesNames{forcesInd}), normTime);
                %Normalise force to body weights
                subjectDataNorm.(subList{subInd}).(runLabels{runInd}).forces.(forcesNames{forcesInd})(cycleInd,:) = ...
                    subjectData.(subList{subInd}).(runLabels{runInd}).forces.(forcesNames{forcesInd})(cycleInd,:) / (trialMass * 9.80665); 
            end
            
            %%%%% GRF DATA %%%%%
            
            %Load in grf file
            grfData = TimeSeriesTable([runLabels{runInd},'_',cycleLabels{cycleInd},'_grf.mot']);

            %Convert data to structure
            grfStruct = osimTableToStruct(grfData);
            
            %Get fieldnames from structure
            grfNames = fieldnames(grfStruct);
            
            %Keep right hand side labels
            grfNames = grfNames(startsWith(grfNames,'R_'));
            
            %Identify the indices of where to extract based on muscle force
            %data - i.e. the GRF is from the whole trial
            grfStartInd = find(grfStruct.time > forcesStruct.time(1), 1) - 1;
            grfEndInd = find(grfStruct.time > forcesStruct.time(end), 1) - 1;
            grfTime = grfStruct.time(grfStartInd:grfEndInd);
            
            %Create a 101 point normalised version of the time array
            normTime = linspace(grfTime(1), grfTime(end), 101);
            
            %Loop through labels, extract and normalise data
            for grfInd = 1:length(grfNames)                
                %Extract and interpolate the current variable data
                subjectData.(subList{subInd}).(runLabels{runInd}).grf.(grfNames{grfInd})(cycleInd,:) = ...
                    interp1(grfTime, grfStruct.(grfNames{grfInd})(grfStartInd:grfEndInd), normTime);
                %Normalise force to body weights
                subjectDataNorm.(subList{subInd}).(runLabels{runInd}).grf.(grfNames{grfInd})(cycleInd,:) = ...
                    subjectData.(subList{subInd}).(runLabels{runInd}).grf.(grfNames{grfInd})(cycleInd,:) / (trialMass * 9.80665); 
            end
            
            %%%%% JOINT REACTIONS DATA %%%%%

            %Load in joint reactions file
            jrfData = TimeSeriesTable([runLabels{runInd},'_',cycleLabels{cycleInd},'_JointReactions_ReactionLoads.sto']);

            %Convert data to structure
            jrfStruct = osimTableToStruct(jrfData);
            
            %Get fieldnames from structure
            jrfNames = fieldnames(jrfStruct);
            
            %Keep right hand side labels
            jrfNames = jrfNames(startsWith(jrfNames,'ankle_r'));
            
            %Create a 101 point normalised version of the time array
            normTime = linspace(jrfStruct.time(1), jrfStruct.time(end), 101);
            
            %Loop through labels, extract and normalise data
            for jrfInd = 1:length(jrfNames)                
                %Extract and interpolate the current variable data
                subjectData.(subList{subInd}).(runLabels{runInd}).jrf.(jrfNames{jrfInd})(cycleInd,:) = ...
                    interp1(jrfStruct.time, jrfStruct.(jrfNames{jrfInd}), normTime);
                %Normalise force to body weights
                subjectDataNorm.(subList{subInd}).(runLabels{runInd}).jrf.(jrfNames{jrfInd})(cycleInd,:) = ...
                    subjectData.(subList{subInd}).(runLabels{runInd}).jrf.(jrfNames{jrfInd})(cycleInd,:) / (trialMass * 9.80665); 
            end
            
            %%%%% JOINT CONTACT DATA %%%%%

            %Load in joint reactions file
            jcfData = TimeSeriesTable([runLabels{runInd},'_',cycleLabels{cycleInd},'_JointContacts_ReactionLoads.sto']);

            %Convert data to structure
            jcfStruct = osimTableToStruct(jcfData);
            
            %Get fieldnames from structure
            jcfNames = fieldnames(jcfStruct);
            
            %Keep right hand side labels
            jcfNames = jcfNames(startsWith(jcfNames,'ankle_r'));
            
            %Create a 101 point normalised version of the time array
            normTime = linspace(jcfStruct.time(1), jcfStruct.time(end), 101);
            
            %Loop through labels, extract and normalise data
            for jcfInd = 1:length(jrfNames)                
                %Extract and interpolate the current variable data
                subjectData.(subList{subInd}).(runLabels{runInd}).jcf.(jrfNames{jcfInd})(cycleInd,:) = ...
                    interp1(jcfStruct.time, jcfStruct.(jcfNames{jcfInd}), normTime);
                %Normalise force to body weights
                subjectDataNorm.(subList{subInd}).(runLabels{runInd}).jcf.(jcfNames{jcfInd})(cycleInd,:) = ...
                    subjectData.(subList{subInd}).(runLabels{runInd}).jcf.(jcfNames{jcfInd})(cycleInd,:) / (trialMass * 9.80665); 
            end
            
            %%%%% MUSCLE FORCE DIRECTION DATA %%%%%
                        
            %Load in muscle force directions file
            mfdData = TimeSeriesTable([runLabels{runInd},'_',cycleLabels{cycleInd},'_MuscleForceDirection_vectors.sto']);

            %Convert data to structure
            mfdStruct = osimTableToStruct(mfdData);
            
            %Get fieldnames from structure
            mfdNames = fieldnames(mfdStruct);
            
            %Keep right hand side labels
            mfdNames = mfdNames(endsWith(mfdNames,'_on_tibia_r'));
            
            %Create a 101 point normalised version of the time array
            normTime = linspace(mfdStruct.time(1), mfdStruct.time(end), 101);
            
            %Loop through labels, extract and normalise data
            for mfdInd = 1:length(mfdNames)                
                %Extract and interpolate the current variable data
                subjectData.(subList{subInd}).(runLabels{runInd}).mfd.(mfdNames{mfdInd})(cycleInd,:) = ...
                    interp1(mfdStruct.time, mfdStruct.(mfdNames{mfdInd}), normTime);
            end
            
            %%%%% MUSCLE ATTACHMENT SITES DATA %%%%%
            
            %Load in muscle force directions file
            masData = TimeSeriesTable([runLabels{runInd},'_',cycleLabels{cycleInd},'_MuscleForceDirection_attachments.sto']);

            %Convert data to structure
            masStruct = osimTableToStruct(masData);
            
            %Get fieldnames from structure
            masNames = fieldnames(masStruct);
            
            %Keep right hand side labels
            masNames = masNames(endsWith(masNames,'_on_tibia_r'));
            
            %Create a 101 point normalised version of the time array
            normTime = linspace(masStruct.time(1), masStruct.time(end), 101);
            
            %Loop through labels, extract and normalise data
            for masInd = 1:length(masNames)                
                %Extract and interpolate the current variable data
                subjectData.(subList{subInd}).(runLabels{runInd}).mas.(masNames{masInd})(cycleInd,:) = ...
                    interp1(masStruct.time, masStruct.(masNames{masInd}), normTime);
            end
            
            %%%%% DECOMPOSE TIBIAL MUSCLE FORCE VECTORS %%%%%
                        
            %Loop through force names
            for forcesInd = 1:length(forcesNames)
                
                %Check if current force is in muscle directions
                if sum(startsWith(mfdNames, forcesNames{forcesInd})) > 0
                    
                    %Extract the three component name for the muscle
                    muscDirNames = mfdNames(startsWith(mfdNames, forcesNames{forcesInd}));
                    
                    %Extract the force for the current muscle
                    F = subjectData.(subList{subInd}).(runLabels{runInd}).forces.(forcesNames{forcesInd})(cycleInd,:);
                    
                    %Extract the directional vectors
                    vx = subjectData.(subList{subInd}).(runLabels{runInd}).mfd.(muscDirNames{1})(cycleInd,:);
                    vy = subjectData.(subList{subInd}).(runLabels{runInd}).mfd.(muscDirNames{2})(cycleInd,:);
                    vz = subjectData.(subList{subInd}).(runLabels{runInd}).mfd.(muscDirNames{3})(cycleInd,:);
                    
                    %Calculate vector components
                    Fx = vx .* F;
                    Fy = vy .* F;
                    Fz = vz .* F;
                    
                    %Store data in structure
                    %Create variable names
                    varX = [forcesNames{forcesInd},'_X'];
                    varY = [forcesNames{forcesInd},'_Y'];
                    varZ = [forcesNames{forcesInd},'_Z'];
                    %Store data
                    subjectData.(subList{subInd}).(runLabels{runInd}).forcesVec.(char(varX))(cycleInd,:) = Fx;
                    subjectData.(subList{subInd}).(runLabels{runInd}).forcesVec.(char(varY))(cycleInd,:) = Fy;
                    subjectData.(subList{subInd}).(runLabels{runInd}).forcesVec.(char(varZ))(cycleInd,:) = Fz;
                    %Normalise force to body weights
                    subjectDataNorm.(subList{subInd}).(runLabels{runInd}).forcesVec.(char(varX))(cycleInd,:) = Fx / (trialMass * 9.80665);
                    subjectDataNorm.(subList{subInd}).(runLabels{runInd}).forcesVec.(char(varY))(cycleInd,:) = Fy / (trialMass * 9.80665);
                    subjectDataNorm.(subList{subInd}).(runLabels{runInd}).forcesVec.(char(varZ))(cycleInd,:) = Fz / (trialMass * 9.80665);                    
                    
                end
                
            end
            
            %Jump back up a directory to cycle list
            cd('..');
            
        end
        
        %Calculate mean data for current run trial
        
        %Forces
        for forcesInd = 1:length(forcesNames)                
            groupData.(runLabels{runInd}).forces.(forcesNames{forcesInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).forces.(forcesNames{forcesInd}));
            groupDataNorm.(runLabels{runInd}).forces.(forcesNames{forcesInd})(subInd,:) = ...
                mean(subjectDataNorm.(subList{subInd}).(runLabels{runInd}).forces.(forcesNames{forcesInd}));
        end
        
        %Force vectors
        forcesVecNames = fieldnames(subjectDataNorm.(subList{subInd}).(runLabels{runInd}).forcesVec);
        for forcesVecInd = 1:length(forcesVecNames)                
            groupData.(runLabels{runInd}).forcesVec.(forcesVecNames{forcesVecInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).forcesVec.(forcesVecNames{forcesVecInd}));
            groupDataNorm.(runLabels{runInd}).forcesVec.(forcesVecNames{forcesVecInd})(subInd,:) = ...
                mean(subjectDataNorm.(subList{subInd}).(runLabels{runInd}).forcesVec.(forcesVecNames{forcesVecInd}));
        end
        
        %GRF
        for grfInd = 1:length(grfNames)                
            groupData.(runLabels{runInd}).grf.(grfNames{grfInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).grf.(grfNames{grfInd}));
            groupDataNorm.(runLabels{runInd}).grf.(grfNames{grfInd})(subInd,:) = ...
                mean(subjectDataNorm.(subList{subInd}).(runLabels{runInd}).grf.(grfNames{grfInd}));
        end
        
        %JRF
        for jrfInd = 1:length(jrfNames)                
            groupData.(runLabels{runInd}).jrf.(jrfNames{jrfInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).jrf.(jrfNames{jrfInd}));
            groupDataNorm.(runLabels{runInd}).jrf.(jrfNames{jrfInd})(subInd,:) = ...
                mean(subjectDataNorm.(subList{subInd}).(runLabels{runInd}).jrf.(jrfNames{jrfInd}));
        end
        
        %JCF
        for jcfInd = 1:length(jcfNames)                
            groupData.(runLabels{runInd}).jcf.(jrfNames{jcfInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).jcf.(jrfNames{jcfInd}));
            groupDataNorm.(runLabels{runInd}).jcf.(jrfNames{jcfInd})(subInd,:) = ...
                mean(subjectDataNorm.(subList{subInd}).(runLabels{runInd}).jcf.(jrfNames{jcfInd}));
        end
        
        %Muscle force directions
        for mfdInd = 1:length(mfdNames)                
            groupData.(runLabels{runInd}).mfd.(mfdNames{mfdInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).mfd.(mfdNames{mfdInd}));
        end
        
        %Jump back up a directory to run trial list
        cd('..');
        
    end
    
    %Jump back up a directory to subject list
    cd('..');
    
end

%% Compile group means

%Loop through running trials
for runInd = 1:length(runLabels)

    %Forces
    for forcesInd = 1:length(forcesNames)                
        meanData.(runLabels{runInd}).forces.(forcesNames{forcesInd}) = ...
            mean(groupData.(runLabels{runInd}).forces.(forcesNames{forcesInd}));
        meanDataNorm.(runLabels{runInd}).forces.(forcesNames{forcesInd}) = ...
            mean(groupDataNorm.(runLabels{runInd}).forces.(forcesNames{forcesInd}));
    end
    
    %Force vectors
    for forcesVecInd = 1:length(forcesVecNames)                
        meanData.(runLabels{runInd}).forcesVec.(forcesVecNames{forcesVecInd}) = ...
            mean(groupData.(runLabels{runInd}).forcesVec.(forcesVecNames{forcesVecInd}));
        meanDataNorm.(runLabels{runInd}).forcesVec.(forcesVecNames{forcesVecInd}) = ...
            mean(groupDataNorm.(runLabels{runInd}).forcesVec.(forcesVecNames{forcesVecInd}));
    end
    
    %GRF
    for grfInd = 1:length(grfNames)                
        meanData.(runLabels{runInd}).grf.(grfNames{grfInd}) = ...
            mean(groupData.(runLabels{runInd}).grf.(grfNames{grfInd}));
        meanDataNorm.(runLabels{runInd}).grf.(grfNames{grfInd}) = ...
            mean(groupDataNorm.(runLabels{runInd}).grf.(grfNames{grfInd}));
    end

    %JRF
    for jrfInd = 1:length(jrfNames)                
        meanData.(runLabels{runInd}).jrf.(jrfNames{jrfInd}) = ...
            mean(groupData.(runLabels{runInd}).jrf.(jrfNames{jrfInd}));
        meanDataNorm.(runLabels{runInd}).jrf.(jrfNames{jrfInd}) = ...
            mean(groupDataNorm.(runLabels{runInd}).jrf.(jrfNames{jrfInd}));
    end
    
    %JCF
    for jcfInd = 1:length(jcfNames)                
        meanData.(runLabels{runInd}).jcf.(jrfNames{jcfInd}) = ...
            mean(groupData.(runLabels{runInd}).jcf.(jrfNames{jcfInd}));
        meanDataNorm.(runLabels{runInd}).jcf.(jrfNames{jcfInd}) = ...
            mean(groupDataNorm.(runLabels{runInd}).jcf.(jrfNames{jcfInd}));
    end
    

    %Muscle force directions
    for mfdInd = 1:length(mfdNames)                
        meanData.(runLabels{runInd}).mfd.(mfdNames{mfdInd}) = ...
            mean(groupData.(runLabels{runInd}).mfd.(mfdNames{mfdInd}));
    end

end
    
%% Save compiled outputs

%Save desired variables
save('compiledHamnerDelp2013-RunningData.mat', ...
    'groupData', 'groupDataNorm', 'meanData', 'meanDataNorm', 'subjectData', 'subjectDataNorm');

%% ----- End of compileHamnerDelpData.m ----- %%
