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

%% Compile subject data

%Navigate to data directory
cd('..\Data\');

%%%%%% TODO: consider we eventually need to convert forces to BWs to match
%%%%%% up with shape model data mean bodyweight to force...

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
            
            %Jump back up a directory to cycle list
            cd('..');
            
        end
        
        %Calculate mean data for current run trial
        
        %Forces
        for forcesInd = 1:length(forcesNames)                
            groupData.(runLabels{runInd}).forces.(forcesNames{forcesInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).forces.(forcesNames{forcesInd}));
        end
        
        %GRF
        for grfInd = 1:length(grfNames)                
            groupData.(runLabels{runInd}).grf.(grfNames{grfInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).grf.(grfNames{grfInd}));
        end
        
        %JRF
        for jrfInd = 1:length(jrfNames)                
            groupData.(runLabels{runInd}).jrf.(jrfNames{jrfInd})(subInd,:) = ...
                mean(subjectData.(subList{subInd}).(runLabels{runInd}).jrf.(jrfNames{jrfInd}));
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
    end

    %GRF
    for grfInd = 1:length(grfNames)                
        meanData.(runLabels{runInd}).grf.(grfNames{grfInd}) = ...
            mean(groupData.(runLabels{runInd}).grf.(grfNames{grfInd}));
    end

    %JRF
    for jrfInd = 1:length(jrfNames)                
        meanData.(runLabels{runInd}).jrf.(jrfNames{jrfInd}) = ...
            mean(groupData.(runLabels{runInd}).jrf.(jrfNames{jrfInd}));
    end

    %Muscle force directions
    for mfdInd = 1:length(mfdNames)                
        meanData.(runLabels{runInd}).mfd.(mfdNames{mfdInd}) = ...
            mean(groupData.(runLabels{runInd}).mfd.(mfdNames{mfdInd}));
    end

end
    
%% Save compiled outputs

%%%%% TODO: consider this after normalising to BW...










