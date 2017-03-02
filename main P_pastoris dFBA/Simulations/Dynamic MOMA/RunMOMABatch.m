clear all

% Initial Conditions
x0 = [  0.5; ...    % Volume
        1; ...    % Biomass
        40; ...     % Glucose
        0.2; ...    % Etanol
        0.1; ...   % Piruvato
        0.3; ...    % Arabitol
        1; ...    % Citrato
        0.00001];     % HSA

simTime = 30;

% Directory where the batch model scripts are kept
oldFolder = cd('C:\Users\Francisco\Dropbox\Paper dFBA Pichia pastoris\Supplementary Material\main P_pastoris dFBA\Batch model');

% it.results file with parameter information
load('Parameter_Values_Example.mat')

k = it_results.k_SS;
kfixed = it_results.kfixed;
assignin('base','kfixed',kfixed)

[tWT,xWT] = run_PPdFBA(k,x0,simTime,[],[],[]);
figure(1)
subplot(2,2,1)
plot(tWT,xWT(:,2))
title('Biomass')
subplot(2,2,2)
plot(tWT,xWT(:,3))
title('Glucose')
subplot(2,2,3)
plot(tWT,xWT(:,6))
title('Arabitol')
subplot(2,2,4)
plot(tWT,xWT(:,8))
title('HSA')


%% Simulación MOMA - análisis de todos los genes demora 8 horas aprox

ngenes = 669;
DynamicMOMAdata = cell(1,ngenes+1); % Modelo contiene 669 genes

% WT phenotype
DynamicMOMAdata{1,end} = [tWT,xWT]; 

for i=1:ngenes
    [tKO,xKO] = run_PPdFBA(k,x0,simTime,i,[],1);
    DynamicMOMAdata{1,i} = [tKO,xKO];
end

cd(oldFolder)
save('DynamicMOMAdata.mat','DynamicMOMAdata')

