%% This script processes the data from the Hamner & Delp (2013) dataset.
%  Specifically it runs an analysis to extract the joint reaction forces
%  and muscle force directions.
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

%%%% TODO: I think force direction analysis can be set up better to also
%%%% extract the muscle force vectors contributing tok ankle joint reaction
%%%% forces - i.e. the distal attachments of those muscles that cross the
%%%% ankle joint...

%Turn off warnings
warning off

%Set home directory
homeDir = pwd;

%Import OpenSim libraries
import org.opensim.modeling.*

%Load in muscle force direction plugin
cd('..\..\MuscleForceDirectionPlugin');
opensimCommon.LoadOpenSimLibraryExact([pwd,'\MuscleForceDirection.dll']);

%Return to home directory
cd(homeDir);

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

%Create analyses to append to the tool

%Joint reactions
jointReactionAnalysis = JointReaction();
jointReactionAnalysis.setName('JointReactions');
%Set array string of parent to use
parentArrayStr = ArrayStr();
parentArrayStr.append('parent');
%Set in analysis
jointReactionAnalysis.setOnBody(parentArrayStr);
jointReactionAnalysis.setInFrame(parentArrayStr);
%Set to only record right ankle joint reaction forces
jointsArrayStr = ArrayStr();
jointsArrayStr.append('ankle_r');
jointReactionAnalysis.setJointNames(jointsArrayStr);

%Joint contacts
jointContactAnalysis = JointReaction();
jointContactAnalysis.setName('JointContacts');
%Set in analysis
jointContactAnalysis.setOnBody(parentArrayStr);
jointContactAnalysis.setInFrame(parentArrayStr);
%Set to only record right ankle joint contact forces
jointContactAnalysis.setJointNames(jointsArrayStr);

%Muscle force direction
%Given the custom plugin approach here we need to load it in from a blank
%setup file
mfdAnalyzeTool = AnalyzeTool('..\..\MuscleForceDirectionPlugin\muscleForceDirectionSetupFile.xml', false);

%Create blank analyze tool to iteratively change
analyzeTool = AnalyzeTool();

%Append analyses to tool
%Joint reactions
analyzeTool.updAnalysisSet().cloneAndAppend(jointReactionAnalysis);
%Muscle force direction
analyzeTool.updAnalysisSet().cloneAndAppend(mfdAnalyzeTool.updAnalysisSet().get(0));

%Create second analyses tool just for the joint contacts
analyzeTool2 = AnalyzeTool();

%Append analyses to tool
%Joint contacts
analyzeTool2.updAnalysisSet().cloneAndAppend(jointContactAnalysis);

%% Process data

%Loop through subjects
for subInd = 1:length(subList)

    %Navigate to subjects data directory
    cd(['..\Data\',subList{subInd}]);
    
    %Loop through running trials
    for runInd = 1:length(runLabels)
        
        %Navigate to folder
        cd(runLabels{runInd});

        %Loop through cycles
        for cycleInd = 1:length(cycleLabels)
            
            %Navigate to folder
            cd(cycleLabels{cycleInd});
            
            %Set logger to record analysis
            Logger.removeFileSink();
            Logger.addFileSink([runLabels{runInd},'_',cycleLabels{cycleInd},'_analyzeLog.log']);
            
            %Set tool name
            analyzeTool.setName([runLabels{runInd},'_',cycleLabels{cycleInd}]);
            
            %Set general inputs in tool
            %Model file
            analyzeTool.setModelFilename([runLabels{runInd},'_',cycleLabels{cycleInd},'_model.osim']);
            %States file
            analyzeTool.setStatesFileName([runLabels{runInd},'_',cycleLabels{cycleInd},'_states.sto']);
            %Initial and final time
            analyzeTool.setInitialTime(Storage([runLabels{runInd},'_',cycleLabels{cycleInd},'_states.sto']).getFirstTime());
            analyzeTool.setFinalTime(Storage([runLabels{runInd},'_',cycleLabels{cycleInd},'_states.sto']).getLastTime());
            %External loads
            analyzeTool.setExternalLoadsFileName([runLabels{runInd},'_',cycleLabels{cycleInd},'_grf.xml']);
            
            %For some reason the analyze tool needs to be saved and then
            %imported back in to run properly
            analyzeTool.print([runLabels{runInd},'_',cycleLabels{cycleInd},'_setupAnalyze.xml']);
            runAnalysis = AnalyzeTool([runLabels{runInd},'_',cycleLabels{cycleInd},'_setupAnalyze.xml']);
            
            %Run analyze tool
            runAnalysis.run();
            
            %Set logger to record second analysis
            Logger.removeFileSink();
            Logger.addFileSink([runLabels{runInd},'_',cycleLabels{cycleInd},'_analyzeLog2.log']);
            
            %Set second tool name
            analyzeTool2.setName([runLabels{runInd},'_',cycleLabels{cycleInd}]);
            
            %Read in and update model to disable all muscle forces to
            %extract joint contact instead of reaction forces
            trialModel = Model([runLabels{runInd},'_',cycleLabels{cycleInd},'_model.osim']);
            %Loop through muscles and disable
            for muscleInd = 0:trialModel.updMuscles().getSize()-1
                trialModel.updMuscles().get(muscleInd).set_appliesForce(false)                
            end
            %Finalise connections
            trialModel.finalizeConnections();
            %Print to file
            trialModel.print([runLabels{runInd},'_',cycleLabels{cycleInd},'_model_noMuscles.osim']);
            
            %Set general inputs in tool
            %Model file
            analyzeTool2.setModelFilename([runLabels{runInd},'_',cycleLabels{cycleInd},'_model_noMuscles.osim']);
            %States file
            analyzeTool2.setStatesFileName([runLabels{runInd},'_',cycleLabels{cycleInd},'_states.sto']);
            %Initial and final time
            analyzeTool2.setInitialTime(Storage([runLabels{runInd},'_',cycleLabels{cycleInd},'_states.sto']).getFirstTime());
            analyzeTool2.setFinalTime(Storage([runLabels{runInd},'_',cycleLabels{cycleInd},'_states.sto']).getLastTime());
            %External loads
            analyzeTool2.setExternalLoadsFileName([runLabels{runInd},'_',cycleLabels{cycleInd},'_grf.xml']);
            
            %For some reason the analyze tool needs to be saved and then
            %imported back in to run properly
            analyzeTool2.print([runLabels{runInd},'_',cycleLabels{cycleInd},'_setupAnalyze2.xml']);
            runAnalysis = AnalyzeTool([runLabels{runInd},'_',cycleLabels{cycleInd},'_setupAnalyze2.xml']);
            
            %Run analyze tool
            runAnalysis.run();            
            
            %Return to run trial directory
            cd('..');
            
        end
        
        %Return to subject directory
        cd('..');

    end
        
    %Return to home directory to restart
    cd(homeDir);
    
    %Print confirmation
    disp(['Processed data for ',subList{subInd},'...'])
        
end

%% ----- End of processHamnerDelpData.m ----- %%