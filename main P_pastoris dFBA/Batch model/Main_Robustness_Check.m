

% Experimental data
filename = 'Datos Calibración RPP dFBA.xlsx';
dataset = 2;

initCobraToolbox

% Change Cobra Solver
changeCobraSolver('glpk','LP');
changeCobraSolver('gurobi5','QP');


%Define initial variables and constraints

% model
model = readCbModel('PP_iFS618.xml');

% excMet
metNames ={ 'Volume' 'Biomass' 'glc-D[e]' 'etoh[e]'...
            'pyr[e]' 'abt_D[e]' 'cit[e]' 'ac[e]' 'TAU[e]'};
        
excMet = findMetIDs(model,metNames);

% excRxn
rxnNames = {'Volume' 'BIOMASS' 'EX_glc(e)' 'EX_etoh(e)'...
            'EX_pyr(e)' 'EX_abt_D(e)' 'EX_cit(e)' 'EX_ac(e)' 'EX_tau(e)'};

rxnIDs = findRxnIDs(model,rxnNames);

excRxn = [rxnIDs' zeros(size(rxnIDs'))];

% Molecular weight vector
PM = [  0       0       180.16  46.07 ...
        88.06   152.14  192.124 60.05   1000]./1000; % Pesos moleculares de los ácidos no ionizados

% feed
feed = [0 0 300 0 ...
        0 0 0 0   0];   % 300 g/L glucose feed
    
assignin('base','model',model);
assignin('base','excMet',excMet);
assignin('base','excRxn',excRxn);
assignin('base','PM',PM);
assignin('base','feed',feed);
% assignin('base','Vout',expdata(:,1:2));
% assignin('base','dataset',dataset); % Ojo a DATASET
assignin('base','skip_delays',true);

% Read experimental data
cd data
    expdata = xlsread(filename,dataset);
    % Peso seco
    DWRelation = 0.8;
    expdata(:,3) = expdata(:,3)*DWRelation;
cd ..

% Initial conditions
x0 = (expdata(1,2:end)+expdata(2,2:end))./2;
texp    = expdata(:,1);
ydata   = expdata(:,2:end);
weights = ones(9,1);
weights(2) = 0.1; % Privilegiar ajuste biomasa 
weights(3) = 0.2; % Privilegiar ajuste glucosa
weights(9) = 0.5; % Privilegiar ajuste taumatina

% Save in base workspace
assignin('base','x0',x0);
assignin('base','texp',texp);
assignin('base','ydata',ydata);
assignin('base','weights',weights);

% Detalles del experimento
experiment = 'CaC4';
PreRegResults = ['it_results_d2_pre_' experiment '.mat'];
RobustStruct = ['cmp_d2_' experiment '.mat'];
Robustness_Info = Test_Robust_Structures(PreRegResults,RobustStruct);

save(['Robustness_Info_' experiment '.mat'],'Robustness_Info')