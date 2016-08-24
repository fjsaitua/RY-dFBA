%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [t,x] = run_dFBA(dataset,k)
% Simulates a batch or fed-batch run of the procedure, for a given set of
% parameters.
%
% INPUTS:
% dataset       Number indicating wich sheet will be analyzed
% k             Parameter values
%
% OUTPUTS:
% t             Time vector of the simulation
% x             Temporal evolution of variables: First column is volume
%               (in [L]), the rest is metabolites ordered as DATA.xls
%               indicates (in [g/L])
%
% Benjamín J. Sánchez
% Last Update: 2014-11-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,x] = run_dFBA(filename,dataset,k)

%Initialize COBRA
initCobraToolbox

% Change Cobra Solver
changeCobraSolver('glpk','LP');
changeCobraSolver('gurobi5','QP');

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


% Read experimental data

cd data
    expdata = xlsread(filename,dataset);
    % Peso seco
    DWRelation = 0.8;
    expdata(:,3) = expdata(:,3)*DWRelation;
cd ..

% Initial conditions
x0 = expdata(1,2:end)';

texp    = expdata(:,1);
ydata   = expdata(:,2:end);
assignin('base','texp',texp);
assignin('base','ydata',ydata);



[t,~]   = size(expdata);
simTime = expdata(t,1);
kfixed  = NaN(size(k));

assignin('base','model',model)
assignin('base','excMet',excMet)
assignin('base','excRxn',excRxn)
assignin('base','x0',x0)
assignin('base','feed',feed)
assignin('base','PM',PM)
%assignin('base','Vout',expdata(:,1:2))
%assignin('base','dataset',dataset)
%assignin('base','trans',trans)
assignin('base','skip_delays',false)
assignin('base','kfixed',kfixed)
%assignin('base','t_limit',100)
assignin('base','expdata',expdata)

%Integrate
odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',0.7,'NonNegative',1:length(x0));
total_time = tic;

% Simtime continuo
%[t,x]=ode113(@pseudoSteadyState,[0 simTime+2],x0,odeoptions,k);

% Simtime discreto
[t,x]=ode113(@pseudoSteadyState,texp,x0,odeoptions,k);
clear pseudoSteadyState

%Show results
disp(['Simulation time: ',num2str(toc(total_time)),' seconds'])
printResults(t,x,expdata)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%