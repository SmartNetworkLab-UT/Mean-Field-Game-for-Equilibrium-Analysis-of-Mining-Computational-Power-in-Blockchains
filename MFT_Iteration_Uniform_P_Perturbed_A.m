% This figure shows the Mean Field Term (MFT) evolution for different
% iterations.
% Parameter p is chosen according to a uniform distribution.
% Also, Parameter A is perturbed mid simulation to analyze the convergence
% of the mean field term under dynamic changes of paramter.

% The parameters are chosen as follows:
% A_i = 100000 $
% p_i = uniform(2, 9) $ per EH
% x_min = 0
% x_max = 100
% x_initial = uniform(x_min, x_max)
% N = 1000
% Number of runs to calculate mean and standard deviation = NumRun 
% Maximum number of iterations = MaxIter

tic;

addpath('Core/')

% Determining the number of simulations to be run. The mean and standard
% deviation are calculated on these different runs.
NumRun = 5;

% Determining the maximum number of iteration of the simulation
MaxIter = 1200;

% Determining the perturbation parameters
IterChange = 500; % At this iteration, the A parameter increases by chosen percent.
A_PercentChange = 50; % The proportion of change. new_A = (1+ A_PercentChange/100)*old_A

% Determining the number of agents. Other parameters are set below, when
% the agents are initialized
N = 1000;

% Allocating the results buffer
AllRunResults = zeros(MaxIter, NumRun);


for simulationRun = 1:NumRun
%% Creating agents
ListOfAgents = [];
for i = 1:N
    x_min = 0;
    x_max = 100;
    A = 100000;
    p = 2 + 7*rand; %uniform(2, 9)
    gamma = 0.5;
    x = x_min + (x_max-x_min)*rand; %uniform(x_min, x_max)
    NewAgent = Agent_V2(x_min, x_max, A, p, gamma, x);
    ListOfAgents = [ListOfAgents; NewAgent]; %#ok<AGROW>
end

%% Simulating the system
MF_Iter = zeros(1, MaxIter);
MF_Iter(1) = 35000; 
for iterationIndex = 1:(MaxIter-1)
    for player = 1:N
        % Players update their strategy for the next round
        ListOfAgents(player).UpdateStrategyMF(MF_Iter(iterationIndex)); 
    end
    % MFT is estimated after the round is finished
    MF_Iter(iterationIndex+1) = SummerOfStrategies(ListOfAgents);
    % Applying the change in A in the designated iteration
    if iterationIndex == IterChange
        for player = 1:N
            ListOfAgents(player).A = (1 + A_PercentChange/100)*ListOfAgents(player).A; %#ok<SAGROW>
        end
    end
end
AllRunResults(:, simulationRun) = MF_Iter;
end


%% Calculating Mean and STD and Ploting the results
Mean_MF_Iter = mean(AllRunResults, 2)'; % Mean calculation
STD_MF_Iter = std(AllRunResults, 1, 2)'; % Standard deviation calculation
plot(1:MaxIter, Mean_MF_Iter-STD_MF_Iter, 'DisplayName','\mu-\sigma', 'LineWidth',1.5);
hold on
plot(1:MaxIter, Mean_MF_Iter, 'DisplayName','\mu', 'LineWidth',2.0);
hold on
plot(1:MaxIter, Mean_MF_Iter+STD_MF_Iter, 'DisplayName','\mu+\sigma', 'LineWidth',1.5)
title('MFT-iteration: mean \pm standard deviation')
legend('FontSize', 14);
toc;