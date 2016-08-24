clear all

% Cargar condiciones iniciales
x0 = [  0.5; ...    % Volumen
        0.2; ...    % Biomasa
        40; ...     % Glucosa
        0.2; ...    % Etanol
        0.1; ...   % Piruvato
        0.3; ...    % Arabitol
        1; ...    % Citrato
        0.01];     % Taumatina

simTime = 30;

k = [5.67 0.01 0.05 0.03 0.08 -.20 1e-3 0.000318153 2.18];

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

ngenes = 669;
DynamicMOMAdata = zeros(length(x0),ngenes+1); % Modelo contiene 669 genes

% WT phenotype
DynamicMOMAdata(:,1) = xWT(end,:)'; 

for i=1:ngenes
    [tKO,xKO] = run_PPdFBA(k,x0,simTime,i,[],1);
    DynamicMOMAdata(:,i+1) = xKO(end,:)';
end


