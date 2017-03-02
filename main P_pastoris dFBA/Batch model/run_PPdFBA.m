%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [t,x] = run_dFBA(k,x0,simTime,geneID,ifMOMA)
% Simulates a batch or fed-batch run of the procedure, for a given set of
% parameters.
%
% INPUTS:
% k             Parameter values
% x0            Initial Conditions of the simulation
% simTime       Time of the simulation
% geneID        Index of the gene to be deleted
% ifMOMA        1 if MOMA is going to be executed
%
% OUTPUTS:
% t             Time vector of the simulation
% x             Temporal evolution of variables: First column is volume
%               (in [L]), the rest is metabolites ordered as DATA.xls
%               indicates (in [g/L])
%
% Benjamín J. Sánchez
% Last Update: 2016-12-22 Francisco Saitua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,x] = run_PPdFBA(k,x0,simTime,geneID,expdata,ifMOMA)
if isempty(ifMOMA)
    %Initialize COBRA
    initCobraToolbox
else
    if geneID==1
        initCobraToolbox
    end
end

% Change Cobra Solver
changeCobraSolver('gurobi5','LP');
changeCobraSolver('gurobi5','QP');

model = readCbModel('iFS670.xml');
model = changeRxnBounds(model,{'EX_tau(e)' 'EX_fab(e)'},[0 0],'b');

% Apply singleGeneDeletion
if isempty(geneID)==0 && isempty(ifMOMA)
    [~,results] = findRxnsFromGenes(model,model.genes(geneID),[],1);
    model = changeRxnBounds(model,results(:,1),zeros(size(results(:,1))),'b');
end

% excMet

metNames ={ 'Volume' 'Biomass' 'glc-D[e]' 'etoh[e]'...
            'pyr[e]' 'abt_D[e]' 'cit[e]'};

% metNames ={ 'Volume' 'Biomass' 'glc-D[e]' 'etoh[e]'...
%             'pyr[e]' 'abt_D[e]' 'cit[e]' 'HSA[e]'}; % Activate if runing dynamic MOMA
        
excMet = findMetIDs(model,metNames);

% excRxn
rxnNames = {'Volume' 'BIOMASS' 'EX_glc(e)' 'EX_etoh(e)'...
            'EX_pyr(e)' 'EX_abt_D(e)' 'EX_cit(e)'}; 


% rxnNames = {'Volume' 'BIOMASS' 'EX_glc(e)' 'EX_etoh(e)'...
%             'EX_pyr(e)' 'EX_abt_D(e)' 'EX_cit(e)' 'EX_HSA(e)'}; % Activate if runing dynamic MOMA

rxnIDs = findRxnIDs(model,rxnNames);

excRxn = [rxnIDs' zeros(size(rxnIDs'))];

% Molecular weight vector
PM = [  0       0       180.16  46.07 ...
        88.06   152.14  192.124 1000]./1000; % Pesos moleculares de los ácidos no ionizados


%assignin('base','kfixed',NaN(1,9))
assignin('base','model',model)
assignin('base','excMet',excMet)
assignin('base','excRxn',excRxn)
assignin('base','x0',x0)
assignin('base','PM',PM)
assignin('base','skip_delays',false)
assignin('base','geneID',geneID)

%Integrate
odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',0.7,'NonNegative',1:length(x0));
total_time = tic;

% Simtime continuo Simulación Parental o SGD
if isempty(ifMOMA)
%    [t,x]=ode113(@pseudoSteadyState_simulation,[0 simTime+2],x0,odeoptions,k);
    [t,x]=ode113(@pseudoSteadyState,[0 simTime+2],x0,odeoptions,k);
else
    [t,x]=ode113(@pseudoSteadyState_simulationMOMA,[0 simTime+2],x0,odeoptions,k);
end
% Simtime discreto
clear pseudoSteadyState

%Show results
disp(['Simulation time: ',num2str(toc(total_time)),' seconds'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%