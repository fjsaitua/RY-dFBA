%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it_results = iteration_complete(dataset,kfixed)
% Does a complete iteration of the procedure, including parameter
% estimation and all the pre/post-regression analysis metrics.
%
% INPUTS:
% dataset       Number indicating wich sheet will be analyzed
% kfixed        Vector indicating which parameters are fixed. Complete the
%               rest of the parameters with NaN
%
% OUTPUT:
% it_results    Structure containing all iteration results
%
% Benjamín J. Sánchez
% Last Update: 2014-11-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function it_results = iteration_complete_valid(filename,dataset,kfixed,n_struct)

%Initialize SSm
cd ssm
    ssm_startup
cd ..

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

%========================= PROBLEM SPECIFICATIONS=========================

problem.f = 'minSquares';  %Objective function (.m file)

ParamValues =  [1       6           8; ... % 1.- Vmax glucose uptake [mmol/gDCW*h]
                1e-5    2.7e-3      1e-2;  ... % 2.- Kg glucose saturation constant [g/L]
                0       0.4         2;  ... % 3.- Glucose-ethanol minimum yield BATCH
                0       0.117       1;  ... % 4.- Glucose-pyruvate minimum yield BATCH
                0       0.1364      1;  ... % 5.- Glucose-arabitol minimum yield BATCH
                0       0           1;  ... % 6.- Glucose-citrate minimum yield BATCH
                0       0           1e-3;  ... % 7.- Glucose-thaumatin minimum yield BATCH
                0       0           1e-3;  ... % 8.- Parámetro 'a' función objetivo
                0       2.18      10];    % 9.- ATP mantención

%Decide the parameters to be estimated depending on kfixed:
m = length(kfixed);
j = 1;
for i = 1:m
    if isnan(kfixed(i))
        problem.x_L(j) = ParamValues(i,1);
        problem.x_0(j) = ParamValues(i,2);
        problem.x_U(j) = ParamValues(i,3);
        j = j+1;
    end    
end

%Problem specifications:
opts.maxeval      = 4000;
opts.local.n1     = 1000;    
opts.maxtime      = 1.5e5;
opts.strategy     = 3;
opts.local.solver = 'n2fb';
opts.local.finish = 'lsqnonlin';

%========================== DATA EXPERIMENTAL ============================
% Read experimental data

cd data
    expdata = xlsread(filename,dataset);
    % Peso seco
    DWRelation = 0.57;
    expdata(:,3) = expdata(:,3)*DWRelation;
cd ..

% Initial conditions

x0 = expdata(1,2:9);
weights = expdata(:,10:17);
texp    = expdata(:,1);
ydata   = expdata(:,2:9);

weights(:,2) = 0.5*weights(:,2); % Privilegiar ajuste biomasa 
% weights(:,3) = 0.07*weights(:,3); % Privilegiar ajuste glucosa
% weights(:,4) = 0.15*weights(:,4); % Privilegiar ajuste etanol
% weights(:,6) = 0.15*weights(:,6); % Privilegiar ajuste arabitol
%weights(:,7) = 0.5*weights(:,7); % Privilegiar ajuste citrato
% weights(:,8) = 0.05*weights(:,8); % Privilegiar ajuste taumatina

assignin('base','x0',x0);
assignin('base','texp',texp);
assignin('base','ydata',ydata);
assignin('base','weights',weights);

%================================== OPTIMIZATION =========================

assignin('base','kfixed',kfixed);

time1   = tic;
results = ess_kernel(problem,opts,texp,ydata);
t1      = toc(time1);

%============================= RESULTS ===================================

k_SS       = results.xbest;
simTime    = texp(length(texp));
odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',0.7,'NonNegative',1:length(x0));

time2 = tic;
if feedFunction(20)==0
    %Batch fermentation, works faster with ode113
    [t,x] = ode113(@pseudoSteadyState,[0 simTime],x0,odeoptions,k_SS);
else
    %Fed-Batch fermentation, works faster with ode15s
    [t,x] = ode15s(@pseudoSteadyState,[0 simTime],x0,odeoptions,k_SS);
end
t2 = toc(time2);

%Save fitting results as a figure
printResults(t,x,expdata,k_SS)
saveas(gcf,['fitting_d' num2str(dataset) '_struct_' num2str(n_struct) '.fig'])

%======================== REGRESSION ANALYSIS ============================

close all
time3 = tic;
[AICc,CI_SS,CC_SS,Mc,Ms,diff] = reg_analysis_complete(k_SS);
t3    = toc(time3);

%Adjust the results to easier integration with Excel:
k  = zeros(1,m);
CC = zeros(1,m);
CI = zeros(1,m);
j  = 0;
for i = 1:m
    if isnan(kfixed(i))
        j     = j+1;
        k(i)  = k_SS(j);
        CC(i) = CC_SS(j);
        CI(i) = CI_SS(j);
    end    
end

%============================== SAVE RESULTS =============================

it_results.kfixed     = kfixed;
it_results.k_SS       = k;
it_results.J_SS       = results.fbest;
it_results.AICc       = AICc;
it_results.CI         = CI;
it_results.CC         = CC;
it_results.Mc         = Mc;
it_results.Ms         = Ms;
it_results.diff       = diff;
it_results.time       = [t1;t2;t3];
it_results.ess_report = load('ess_report.mat');
it_results.ktofix     = decision(Mc,Ms);

if j == length(kfixed)
    save(['it_results_struct_' num2str(n_struct) '_dataset_' num2str(dataset) '_pre.mat'],'it_results');
else
    save(['it_results_struct_' num2str(n_struct) '_dataset_' num2str(dataset) '_post.mat'],'it_results');
end

%Delete other files:
delete('Sensib*')
delete('Mc.txt')
delete('GraficoBarra.txt')
delete('checkpoint.mat')
delete('ess_report.mat')
close all
set(0,'ShowHiddenHandles','on');delete(get(0,'Children'))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%