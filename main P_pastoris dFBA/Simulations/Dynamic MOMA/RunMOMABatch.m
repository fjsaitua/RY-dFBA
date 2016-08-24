clear all

% Cargar condiciones iniciales
x0 = [  0.5; ...    % Volumen
        1; ...    % Biomasa
        40; ...     % Glucosa
        0.2; ...    % Etanol
        0.1; ...   % Piruvato
        0.3; ...    % Arabitol
        1; ...    % Citrato
        0.00001];     % HSA

simTime = 40;

oldFolder = cd('C:\Users\Francisco\Dropbox\Tesis Magister\RPP_dFBA BATCH\RPP_dFBA - Modelo Completo BATCH_ Lab Biotec');
load('it_results_struct_31_dataset_9_post.mat')

k = it_results.k_SS;
kfixed = it_results.kfixed;
kfixed(7) = 1e-7;
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
title('Thaumatin')


%% Simulación MOMA - análisis de todos los genes demora 8 horas aprox

ngenes = 670;
DynamicMOMAdata = cell(1,ngenes+1); % Modelo contiene 669 genes

% WT phenotype
DynamicMOMAdata{1,end} = [tWT,xWT]; 

for i=351:ngenes
    [tKO,xKO] = run_PPdFBA(k,x0,simTime,i,[],1);
    DynamicMOMAdata{1,i} = [tKO,xKO];
end

cd(oldFolder)
save('DynamicMOMAdata_351_670_HSA_revised.mat','DynamicMOMAdata')

