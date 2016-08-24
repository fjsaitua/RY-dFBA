clear all
%% Inicialización del modelo y declaración de variables globales

initCobraToolbox

% Change Cobra Solver
changeCobraSolver('gurobi5','LP');
changeCobraSolver('gurobi5','QP');


%Define initial variables and constraints

% model
model = readCbModel('PP_iFS618.xml');

% excMet
metNames ={ 'Volume' 'Biomass' 'glc-D[e]' 'etoh[e]'...
            'pyr[e]' 'abt_D[e]' 'cit[e]' 'TAU[e]'};
        
excMet = findMetIDs(model,metNames);

% excRxn
rxnNames = {'Volume' 'BIOMASS' 'EX_glc(e)' 'EX_etoh(e)'...
            'EX_pyr(e)' 'EX_abt_D(e)' 'EX_cit(e)' 'EX_tau(e)'};

rxnIDs = findRxnIDs(model,rxnNames);

excRxn = [rxnIDs' zeros(size(rxnIDs'))];

% Molecular weight vector
PM = [  0       0       180.16  46.07 ...
        88.06   152.14  192.124 1000]./1000; % Pesos moleculares de los ácidos no ionizados

% feed
feed = [0 0 300 0 ...
        0 0 0 0];   % 300 g/L glucose feed

    
   
assignin('base','model',model);
assignin('base','excMet',excMet);
assignin('base','excRxn',excRxn);
assignin('base','PM',PM);
assignin('base','feed',feed);
% assignin('base','Vout',expdata(:,1:2));
% assignin('base','dataset',dataset); % Ojo a DATASET
assignin('base','skip_delays',true);


%% Inicializar datos experimentales
% Read experimental data
filename = 'Calibration data Picha pastoris dFBA.xlsx';
dataset = 3;
cd data
    expdata = xlsread(filename,9);
    % Peso seco
    DWRelation = 0.57;
    expdata(:,3) = expdata(:,3)*DWRelation;
cd ..

% Initial conditions

x0 = expdata(1,2:9);
weights = expdata(:,10:17);
texp    = expdata(:,1);
ydata   = expdata(:,2:9);


assignin('base','x0',x0);
assignin('base','texp',texp);
assignin('base','ydata',ydata);
assignin('base','weights',weights);

%Normalize data with maximum measures:
% for i = 1:8
%     ydata(:,i) = ydata(:,i)./(max(ydata(:,i))*weights(:,i));
% end

%% Preparar vectores con parámetros
% load(['it_results_d' num2str(dataset) '_pre.mat'])
% kfixed = NaN(1,9);
% assignin('base','kfixed',kfixed);
% k = it_results.k_SS;

% % Análisis post
%load(['it_results_d' num2str(dataset) '_post.mat'])
load('it_results_struct_4_dataset_9_post_intento 0.mat')
kfixed = it_results.kfixed;
assignin('base','kfixed',kfixed);
k = it_results.k_SS;
k(find(~isnan(kfixed))) = [];

kL = k.*0.9999999999;
kU = k.*1.0000000001 + ones(1,length(k))*1e-10;
for i = 1:length(k)
    if k(i) < 0
        kL(i) = k(i)*1.0000000001;
        kU(i) = k(i)*0.9999999999;
    end
end

%% Análisis de significancia
% Falta yexp
options = optimset('MaxIter',1,'MaxFunEvals',1,'TolFun',1,'TolX',1);

[k2,~,res,~,~,~,Jac] = lsqcurvefit(@minSquares2,k,texp,ydata,kL,kU,options);
diff = sum((k-k2).^2);
%Confidence intervals:
[CI,CC] = intconfianza(Jac,res,k2,0.05);
CC      = CC./100;
m       = length(k2);
disp('Unsignificant Parameters:');

    for i = 1:m
        if k2(i) == 0
            CC(i) = 0;
        elseif CC(i) >= 2
            disp(num2str(i));
        end
    end
    