% This figure shows the Mean Field Term (MFT) evolution for different
% iterations.
% Parameter A is chosen according to a uniform distribution.
% Also, Parameter p is chosen according to a uniform distribution.
% The learning rate, gamma, is changed in this simulation to see the effect
% of the learning rate on mean field term evolution.

% The parameters are chosen as follows:
% A_i = uniform(80000, 100000) $
% p_i = uniform(2, 9) $ per EH
% x_min = 0
% x_max = 100
% x_initial = uniform(x_min, x_max)
% N = 1000
% Maximum number of iterations = MaxIter
% To mitigate the effect of randomness, the simulation is run
% SimulationRunNum times and averaged on those runs.
tic;

addpath('Core/')

% Determining the N values to be simulated
gamma_Vector = [0.3, 0.5, 0.7];

% Determining the maximum number of iteration of the simulation
MaxIter = 200;
% Determining the number of simulation to be average on.
SimulationRunNum = 5;

% Determining the number of players in the network
N = 1000;

% Allocating the results buffer
AllRunResults = zeros(MaxIter, length(gamma_Vector));


for simulationRun = 1:length(gamma_Vector)
gamma = gamma_Vector(simulationRun);
%% Creating agents
for iSim = 1:SimulationRunNum
SumMF_Iter = 0;
ListOfAgents = [];
for i = 1:N
    x_min = 0;
    x_max = 100;
    A = 80000+ 20000*rand; %uniform(80000,100000)
    p = 2 + 7*rand; %uniform(2,9)
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
end
SumMF_Iter = SumMF_Iter + MF_Iter;
end
AllRunResults(:, simulationRun) = SumMF_Iter/SimulationRunNum;
end


%% Ploting the results
figure
for iPlot = 1:length(gamma_Vector)
plot(1:MaxIter, AllRunResults(:,iPlot), 'DisplayName',['gamma=',num2str(gamma_Vector(iPlot))], 'LineWidth',1.5);
hold on
end
hold off
title('MFT-iteration for different learning rates')
legend('FontSize', 14);
toc;