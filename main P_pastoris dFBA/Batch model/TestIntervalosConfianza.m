    
clear all
filename = 'Datos Cepa 8 copias.xlsx';
dataset = 1;

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


assignin('base','model',model);
assignin('base','excMet',excMet);
assignin('base','excRxn',excRxn);
assignin('base','PM',PM);
% assignin('base','Vout',expdata(:,1:2));
% assignin('base','dataset',dataset); % Ojo a DATASET
assignin('base','skip_delays',true);

cd data
    expdata = xlsread(filename,dataset);
    % Peso seco
    DWRelation = 0.8;
    expdata(:,3) = expdata(:,3)*DWRelation;
cd ..

x0 = expdata(1,2:end)';
texp    = expdata(:,1);
ydata   = expdata(:,2:end);
weights = ones(9,1);
weights(2) = 0.1; % privilegiar ajuste biomasa
assignin('base','x0',x0);
assignin('base','texp',texp);
assignin('base','ydata',ydata);
assignin('base','weights',weights); 

experiment = 'CaC2';
load(['Robustness_Info_' experiment '.mat'])

% Select Robust structure to test
n = 1;

% NaN indexes
kfixed_n = Robustness_Info{n,1}.kfixed; % Variable auxiliar que guarda el valor de los parámetros fijos y NaNs
assignin('base','kfixed',kfixed_n)

NaNindex = find(isnan(Robustness_Info{n,1}.kfixed));
k = Robustness_Info{n,1}.k_SS;
CI_n = Robustness_Info{n,1}.CI;

% -95% CI
KSS_inf = k*0.99999999;
%KSS_inf(NaNindex) = kfixed_n(NaNindex)*.8;
kL = zeros(1,3);
% +95% CI
KSS_sup = k*1.00000001+ones(1,length(k))*1e-15;
%KSS_sup(NaNindex) = kfixed_n(NaNindex)*1.2;
kU = KSS_sup;

%options = optimset('MaxIter',1,'MaxFunEvals',1,'TolFun',1,'TolX',1);
%[k2,~,res,~,~,~,Jac] = lsqcurvefit(@minSquares2,k,texp,ydata,kL,kU,options); % Aquí está el problema de las tolerancias
[k2,~,res,~,~,~,Jac] = lsqcurvefit(@minSquares2,k,texp,ydata,kL,[]); % Aquí está el problema de las tolerancias


% kfixed_n(NaNindex) = Robustness_Info{n,1}.k_SS;
% [t_0,x_0] = run_dFBA(filename,dataset,kfixed_n);
% res = x_0-ydata;
%Confidence intervals:
[CI,CC] = intconfianza(Jac,res,k2,0.01);
CC      = CC./100;
m       = length(k2);

    