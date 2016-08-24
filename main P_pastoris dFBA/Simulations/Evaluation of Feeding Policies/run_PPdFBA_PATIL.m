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

function [t,x] = run_PPdFBA_PATIL(k,x0,simTime,geneID,expdata,ifMOMA)

%Initialize COBRA
initCobraToolbox

% Change Cobra Solver
changeCobraSolver('gurobi5','LP');
changeCobraSolver('gurobi5','QP');

% model
model = readCbModel('PP_iFS618.xml');
%model = modelConstraints(model);

% Apply singleGeneDeletion
if isempty(geneID)==0 && isempty(ifMOMA)
    [~,results] = findRxnsFromGenes(model,model.genes(geneID),[],1);
    model = changeRxnBounds(model,results(:,1),zeros(size(results(:,1))),'b');
end

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
feed = [0 0 534 0 ...
        0 0 0 0];   % 300 g/L glucose feed

kfixed  = NaN(size(k));

assignin('base','model',model)
assignin('base','excMet',excMet)
assignin('base','excRxn',excRxn)
assignin('base','x0',x0)
assignin('base','feed',feed)
assignin('base','PM',PM)
assignin('base','skip_delays',false)
assignin('base','kfixed',kfixed)
assignin('base','geneID',geneID)

%Integrate
%odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',0.7,'NonNegative',1:length(x0));
odeoptions = odeset('Events',@StopEvents,'RelTol',1e-3,'AbsTol',1e-3,'MaxStep',0.7,'NonNegative',1:length(x0));
total_time = tic;

% Simtime continuo Simulación Parental o SGD
if isempty(ifMOMA)
    %[t,x]=ode15s(@pseudoSteadyState_PATIL,simTime,x0,odeoptions,k);
    [t,x,TE,YE,IE]=ode15s(@pseudoSteadyState_PATIL,simTime,x0,odeoptions,k);
else
    [t,x]=ode113(@pseudoSteadyState_simulationMOMA,[0 simTime],x0,odeoptions,k);
end
% Simtime discreto
%[t,x]=ode113(@pseudoSteadyState,texp,x0,odeoptions,k);
clear pseudoSteadyState

%Show results
disp(['Simulation time: ',num2str(toc(total_time)),' seconds'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%