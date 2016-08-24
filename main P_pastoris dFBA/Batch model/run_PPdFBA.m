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
changeCobraSolver('glpk','LP');
changeCobraSolver('gurobi5','QP');

% model
model = readCbModel('iFS670.xml');

% Apply singleGeneDeletion
if isempty(geneID)==0 && isempty(ifMOMA)
    [~,results] = findRxnsFromGenes(model,model.genes(geneID),[],1);
    model = changeRxnBounds(model,results(:,1),zeros(size(results(:,1))),'b');
end

% excMet
metNames ={ 'Volume' 'Biomass' 'glc-D[e]' 'etoh[e]'...
            'pyr[e]' 'abt_D[e]' 'cit[e]'};
        
excMet = findMetIDs(model,metNames);

% excRxn
rxnNames = {'Volume' 'BIOMASS' 'EX_glc(e)' 'EX_etoh(e)'...
            'EX_pyr(e)' 'EX_abt_D(e)' 'EX_cit(e)'};

rxnIDs = findRxnIDs(model,rxnNames);

excRxn = [rxnIDs' zeros(size(rxnIDs'))];

% Molecular weight vector
PM = [  0       0       180.16  46.07 ...
        88.06   152.14  192.124]./1000; % Pesos moleculares de los ácidos no ionizados

% feed
feed = [0 0 300 0 ...
        0 0 0];   % 300 g/L glucose feed

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
odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',0.7,'NonNegative',1:length(x0));
total_time = tic;

% Simtime continuo Simulación Parental o SGD
if isempty(ifMOMA)
    %[t,x]=ode113(@pseudoSteadyState,[0 simTime+2],x0,odeoptions,k); % MultiObj
    [t,x]=ode113(@pseudoSteadyState,[0 simTime+2],x0,odeoptions,k);
else
    [t,x]=ode113(@pseudoSteadyState_simulationMOMA,[0 simTime+2],x0,odeoptions,k);
end
% Simtime discreto
%[t,x]=ode113(@pseudoSteadyState,texp,x0,odeoptions,k);
clear pseudoSteadyState

%Show results
disp(['Simulation time: ',num2str(toc(total_time)),' seconds'])

%printResults(t,x,expdata)

% Gráficos variables de estado

% variables = {   'Biomass' 'Glucose' 'Ethanol' 'Pyruvate'...
%                 'Arabitol' 'Citrate' 'Acetate' 'Thaumatin'};

% figure(1)
% for i=1:8
%     if isempty(expdata)
%         subplot(3,3,i)
%         plot(t,x(:,i+1),'r-','LineWidth',2)
%         xlim([0 t(end)])
%         title(variables(i))
%         ylabel('[g/L]')
%         xlabel('Time [h]')
%     else
%         subplot(3,3,i)
%         texp = expdata(:,1);
%         plot(t,x(:,i+1),'r-',texp,expdata(:,2+i),'k+','LineWidth',2)
%         xlim([0 t(end)])
%         title(variables(i))
%         ylabel('[g/L]')
%         xlabel('Time [h]')
%     end
% end
% 
% % Gráfico crecimiento subóptimo
% for i=1:length(t)
%     w(i) = Obj_weight(t(i),k(9),k(10),k(11));
% end
% 
% subplot(3,3,9)
% [AX,H1,H2] = plotyy(t,x(:,2),t,w);
% %legend('Biomass','w','Location','Best')
% xlabel('Time [h]')
% ylabel(AX(1),'Biomass [g/L]')
% ylabel(AX(2),'W')
% tlim = t(end);
% ylimW = 0.002;
% xlim(AX(1),[0 tlim])
% xlim(AX(2),[0 tlim])
% ylim(AX(2),[0 ylimW])
% set(AX(1),'YColor','k')
% set(AX(2),'YColor','k','YTick',[0 0.001 0.002])
% set(H1,'Linewidth',2)
% set(H2,'Linewidth',2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%