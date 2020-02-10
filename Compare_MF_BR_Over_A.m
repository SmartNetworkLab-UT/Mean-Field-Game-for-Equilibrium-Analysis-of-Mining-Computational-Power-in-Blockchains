% This figure compares the social welfare for Mean Field and Best Response
% methods 
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
A_Vector = 40000:20000:140000;
A_Range = 20000;
MF_SWF_Vector = zeros(size(A_Vector));
BR_SWF_Vector = zeros(size(A_Vector));

% Determining the maximum number of iteration of the simulation
MaxIter = 600;
% Determining the number of simulation to be average on.
SimulationRunNum = 5;


%% Running simulation for different number of agents for BR and MF Scenario
for simulationRun = 1:length(A_Vector)
A_base = A_Vector(simulationRun);
disp(['BR and MF, A=', num2str(A_base), ':', num2str(A_base+A_Range)]);
%% Creating agents
for iSim = 1:SimulationRunNum
SumMFSWF = 0;
SumBRSWF = 0;
ListOfAgents = [];
for i = 1:N
    x_min = 0;
    x_max = 100;
    A = A_base+ A_Range*rand; %80000 + 40*i; %80000+ 20000*rand; %uniform(80000,100000)
    p = 2 + 7*rand; %2 + 7*i/N; %2 + 7*rand; %uniform(2,9)
    x = x_min + (x_max-x_min)*rand; %uniform(x_min, x_max)
    gamma = 0.5;
    NewAgent = Agent_V2(x_min, x_max, A, p, gamma, x);
    ListOfAgents = [ListOfAgents; NewAgent]; %#ok<AGROW>
end

%% Simulating the system under Mean Field scenario
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
MeanFieldSocialWelfare = 0;
for player = 1:N
    % Calculating the Social Welfare for MF case
    MeanFieldSocialWelfare = MeanFieldSocialWelfare + ListOfAgents(player).UtilityCalc(MF_Iter(end)); 
end

SumMFSWF = SumMFSWF + MeanFieldSocialWelfare;

%% Simulating the system under Best Response scenario
TotalComputationalPower_Iter = zeros(1, MaxIter);
TotalComputationalPower_Iter(1) = 35000; 
for iterationIndex = 1:(MaxIter-1)
    for player = 1:N
        % Players update their strategy for the next round
        ListOfAgents(player).UpdateStrategyBR(TotalComputationalPower_Iter(iterationIndex)); 
    end
    % MFT is estimated after the round is finished
    TotalComputationalPower_Iter(iterationIndex+1) = SummerOfStrategies(ListOfAgents);
end
BestResponseSocialWelfare = 0;
for player = 1:N
    % Calculating the Social Welfare for MF case
    BestResponseSocialWelfare = BestResponseSocialWelfare + ListOfAgents(player).UtilityCalc(TotalComputationalPower_Iter(end)); 
end

SumBRSWF = SumBRSWF + BestResponseSocialWelfare;

end
MF_SWF_Vector(simulationRun) = SumMFSWF/SimulationRunNum;
BR_SWF_Vector(simulationRun) = SumBRSWF/SimulationRunNum;
end


%% Ploting the results
figure
plot(A_Vector, MF_SWF_Vector, 'o-', 'DisplayName', 'MF', 'LineWidth',2);
hold on
plot(A_Vector, BR_SWF_Vector, 'o-', 'DisplayName', 'BR', 'LineWidth',2);
hold off

xticks(A_Vector)
TickLabelCell = cell(1, length(A_Vector));
for i = 1:length(A_Vector)
    TickLabelCell{i} = ['A=', num2str(A_Vector(i)),':', num2str(A_Vector(i)+A_Range)];
end
xticklabels(TickLabelCell)
ylabel('Social Welfare');
title('Social Welfare for Different A parameter')
legend('FontSize', 14);
toc;