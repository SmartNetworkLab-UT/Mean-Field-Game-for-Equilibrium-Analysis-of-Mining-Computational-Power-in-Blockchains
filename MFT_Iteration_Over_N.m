% This figure shows the Mean Field Term (MFT) evolution for different
% iterations.
% Parameter A is chosen according to a uniform distribution.
% Also, Parameter p is chosen according to a uniform distribution.
% The number of players in the network is different to show the behavior of
% the algorithm

% The parameters are chosen as follows:
% A_i = uniform(80000, 100000) $
% p_i = uniform(2, 9) $ per EH
% x_min = 0
% x_max = 100
% x_initial = uniform(x_min, x_max)
% N = 100, 300, 1000
% Number of runs to calculate mean and standard deviation = NumRun 
% Maximum number of iterations = MaxIter

tic;

addpath('Core/')

% Determining the N values to be simulated
N_Vector = [100, 300, 1000];

% Determining the maximum number of iteration of the simulation
MaxIter = 50;


% Allocating the results buffer
AllRunResults = zeros(MaxIter, length(N_Vector));


for simulationRun = 1:length(N_Vector)
N = N_Vector(simulationRun);
%% Creating agents
ListOfAgents = [];
for i = 1:N
    x_min = 0;
    x_max = 100;
    A = 80000+ 20000*rand; %uniform(80000,100000)
    p = 2 + 7*rand; %uniform(2,9)
    gamma = 0.5;
    x = x_min + (x_max-x_min)*rand; %uniform(x_min, x_max)
    NewAgent = Agent_V2(x_min, x_max, A, p, gamma, x);
    ListOfAgents = [ListOfAgents; NewAgent]; %#ok<AGROW>
end

%% Simulating the system
MF_Iter = zeros(1, MaxIter);
MF_Iter(1) = 25000; 
for iterationIndex = 1:(MaxIter-1)
    for player = 1:N
        % Players update their strategy for the next round
        ListOfAgents(player).UpdateStrategyMF(MF_Iter(iterationIndex)); 
    end
    % MFT is estimated after the round is finished
    MF_Iter(iterationIndex+1) = SummerOfStrategies(ListOfAgents);
end
AllRunResults(:, simulationRun) = MF_Iter;
end


%% Ploting the results
for iPlot = 1:length(N_Vector)
plot(1:MaxIter, AllRunResults(:,iPlot), 'DisplayName',['gamma=',num2str(N_Vector(iPlot))], 'LineWidth',1.5);
hold on
end
hold off
title('MFT-iteration for different network sizes')
legend('FontSize', 14);
toc;