function easyRun()

    %% This script is a convenient function for running the FE simulations
    %  for all cases in the dataset using the runSimulations.m function
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

    %Set the list of cases to simulate
    simulateCases = [147211; 102480; 102924; 103559; 103862; 107215; 107813;
        108421; 112802; 113033; 116070; 132433; 132592; 134434; 135418; 135881;
        136283; 140368; 140694; 141682; 144977; 145250; 146511; 149327; 170206;
        174615; 176087; 181140; 182401; 184171];

    %Create waitbar
    simWaitBar = waitbar(0, 'Simulating cases...');

    %% Run the simulations
    
    %Loop through cases
    for caseInd = 1:length(simulateCases)

        %Update waitbar
        waitbar(caseInd/(length(simulateCases)+1), simWaitBar, ...
            ['Simulating case ',num2str(simulateCases(caseInd)),' ...']);

        %Run simulations
        runSimulations(simulateCases(caseInd));

    end

    %% Update waitbar for end
    waitbar(1, simWaitBar, 'Simulations Complete. Time for a well earned rest.');

end