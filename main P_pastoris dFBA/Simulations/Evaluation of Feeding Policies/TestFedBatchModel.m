clear all % Sacar si quiero probar un k_best

k = [   2.3102; ...      % Vmax
        0.9083;...       % Kg
        0.5255;...     % pEtOH
        0.0703;...    % pPyr
        0.1174;...    % pArab
        1.7082e-5;...    % pCit
        1.6104e-6;...       % pAcet
        0;...       % pTau
        -.8056;...   % pEtOH FED
        -.4399;...   % pPyr FED
        -.1233;...   % pArab FED
        2.4255e-5;...       % pCit FED
        2.3177e-8;...       % pAcet FED
        0;...           % pTau FED
        1.1604e-5;...         % w F.O. Batch
        2.0217e-6];         % w F.O. FedBatch
        %2;         % m_atp
        %31.24 ]; % Tsubopt
    
% k = k_best;

%% Experimental data
filename = 'Datos Calibración RPP dFBA.xlsx';
dataset = 5;

expdata = xlsread(filename,dataset);
% Peso seco
DWRelation = 0.8;
expdata(:,3) = expdata(:,3)*DWRelation;


% Initial conditions
x0 = expdata(1,2:end)';
texp    = expdata(:,1);
ydata   = expdata(:,2:end);

simTime = texp(end)+10;

[tWT,xWT] = run_PPdFBA(k,x0,simTime,[],[],[]);

%% Plots
figure(1)
printResultsFedBatch(tWT,xWT,[texp ydata])

