%% This script extracts the relevant data from the Hamner & Delp (2013) 
%  dataset. Note that the dataset needs to be downloaded and placed in
%  appropriate subject folders for this script to work.
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

%Import OpenSim libraries
import org.opensim.modeling.*

%Set home directory
homeDir = pwd;

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

%% Extract data

%Loop through subjects
for subInd = 1:length(subList)
    
    %Create folders to store extracted data in
    
    %Navigate to data directory
    cd('..\Data\');
    
    %Make subject directory
    mkdir(subList{subInd});
    cd(subList{subInd});
    
    %Create directories for each run trial and cycle
    %Loop through run labels
    for runInd = 1:length(runLabels)
        %Create trial level directory
        mkdir(runLabels{runInd});
        cd(runLabels{runInd});
        %Loop through cycle labels
        for cycleInd = 1:length(cycleLabels)
            %Create directory
            mkdir(cycleLabels{cycleInd});
            %Store full path label to variable to store data to later
            exportPaths.(runLabels{runInd}).(cycleLabels{cycleInd}) = ...
                [pwd,'\',cycleLabels{cycleInd}];
        end
        %Jump back up to subject path level
        cd('..');
    end
    
    %Return to home directory
    cd(homeDir);
    
    %Navigate to subject directory
    cd(['..\',subList{subInd}]);
    subDir = pwd;
    
    %Identify cmc folder for subject
    dirList = dir();
    for folderInd = 1:length(dirList)
        if startsWith(dirList(folderInd).name, 'cmc')
            cmcFolderName = dirList(folderInd).name;
            break
        end
    end
    
    %Jump into cmc folder
    cd(cmcFolderName);
    
    %Identify folders for each run label
    dirList = dir();
    for runInd = 1:length(runLabels)
        for folderInd = 1:length(dirList)
            if contains(dirList(folderInd).name, runLabels{runInd})
                cmcRunFolders{runInd} = dirList(folderInd).name;
                break
            end
        end
    end
    
    %Loop through run labels and extract data
    for runInd = 1:length(cmcRunFolders)
        
        %Navigate to folder
        cd(cmcRunFolders{runInd});
        runFolder = pwd;
        
        %Extract the directory list
        dirList = dir();
        
        %Get the cycle labels for the current trial
        currCycleLabels = {};
        for dirInd = 1:length(dirList)
            if contains(dirList(dirInd).name, '_cycle') && ...
                    contains(dirList(dirInd).name, 'CMC_Results')
                %Identify cycle labelling
                fileSplit = strsplit(dirList(dirInd).name, 'cycle');
                currCycleLabels{length(currCycleLabels)+1} = ['cycle',fileSplit{2}];
            end
        end
        
        %Loop through cycles
        for cycleInd = 1:length(currCycleLabels)
            
            %Identify the current results folder
            for dirInd = 1:length(dirList)
                if contains(dirList(dirInd).name, currCycleLabels{cycleInd}) && ...
                        contains(dirList(dirInd).name, 'CMC_Results')
                    cmcResultsFolder = dirList(dirInd).name;
                    break
                end
            end
            
            %Navigate to results folder
            cd(cmcResultsFolder);
            
            %Get file list
            fileList = dir();
            
            %Copy and rename relevant files across to pre-specified folder
            for fileInd = 1:length(fileList)
                %Kinematic data
                if contains(fileList(fileInd).name, '_Kinematics_q.sto')
                    copyfile(fileList(fileInd).name, ...
                        [exportPaths.(runLabels{runInd}).(cycleLabels{cycleInd}),'\',runLabels{runInd},'_',cycleLabels{cycleInd},'_kinematics.sto']);
                %Forces data
                elseif contains(fileList(fileInd).name, '_Actuation_force.sto')
                    copyfile(fileList(fileInd).name, ...
                        [exportPaths.(runLabels{runInd}).(cycleLabels{cycleInd}),'\',runLabels{runInd},'_',cycleLabels{cycleInd},'_forces.sto']);
                %States
                elseif contains(fileList(fileInd).name, '_states.sto')
                    copyfile(fileList(fileInd).name, ...
                        [exportPaths.(runLabels{runInd}).(cycleLabels{cycleInd}),'\',runLabels{runInd},'_',cycleLabels{cycleInd},'_states.sto']);                    
                end
            end
            
            %Jump back out of results directory
            cd('..');
            
            %Identify current cmc setup file
            for dirInd = 1:length(dirList)
                if contains(dirList(dirInd).name, currCycleLabels{cycleInd}) && ...
                        contains(dirList(dirInd).name, 'Setup_CMC')
                    cmcSetupFile = dirList(dirInd).name;
                    break
                end
            end
            
            %Read in cmc setup file
            fid = fopen(cmcSetupFile, 'rt');
            
            %Search for model file indicator
            while 1
                %Get text line
                line = fgetl(fid); 
                %Check for model file indicator
                if contains(line, '<model_file>')
                    break
                end   
            end
            
            %Extract model file path
            modelPath = extractBetween(line, '<model_file>', '</model_file>');
            
            %For some weird reason need to navigate to the path of the
            %model file to extract and sort it out
            %Split model file path into parts
            [modelFolder, modelFileName] = fileparts(modelPath{1});
            %Navigate to folder
            cd(modelFolder);
            %Copy model file to desired directory
            copyfile([modelFileName,'.osim'], ...
                [exportPaths.(runLabels{runInd}).(cycleLabels{cycleInd}),'\',runLabels{runInd},'_',cycleLabels{cycleInd},'_model.osim']);
            
            %Jump back to run folder
            cd(runFolder);
            
            %Search for external loads file indicator
            while 1
                %Get text line
                line = fgetl(fid); 
                %Check for model file indicator
                if contains(line, '<external_loads_file>')
                    break
                end   
            end
            
            %Extract external loads file name
            loadsFile = extractBetween(line, '<external_loads_file>', '</external_loads_file>');
            
            %Read in external loads file
            extLoads = ExternalLoads(loadsFile{1}, true);
            
            %Get the specified datafile path
            motDataFile = char(extLoads.getDataFileName());
            [~, motDataFileName] = fileparts(motDataFile);
            
            %Navigate to subjects experimental data folder
            cd([subDir,'\ExportedData']);
            
            %Copy grfs file over to the desired directory
            copyfile([motDataFileName,'.mot'], ...
                [exportPaths.(runLabels{runInd}).(cycleLabels{cycleInd}),'\',runLabels{runInd},'_',cycleLabels{cycleInd},'_grf.mot']);
            
            %Return to run folder
            cd(runFolder);
            
            %Update external loads file to use relative new filename
            extLoads.setDataFileName([runLabels{runInd},'_',cycleLabels{cycleInd},'_grf.mot']);
            
            %Save new external loads file to desired directory
            extLoads.print([exportPaths.(runLabels{runInd}).(cycleLabels{cycleInd}),'\',runLabels{runInd},'_',cycleLabels{cycleInd},'_grf.xml']);

        end
        
        %Return to top level cmc folder
        cd('..');
        
    end
    
    %Return to home directory to restart
    cd(homeDir);
    
    %Print confirmation
    disp(['Extracted data for ',subList{subInd},'...'])
    
end

%% ----- End of extractHamnerDelpData.m ----- %%