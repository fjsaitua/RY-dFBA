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

function it_results = iteration_complete(filename,dataset,kfixed)

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
feed = [0 0 534 0 ...
        0 0 0 0 0];   % 300 g/L glucose feed

assignin('base','model',model);
assignin('base','excMet',excMet);
assignin('base','excRxn',excRxn);
assignin('base','PM',PM);
assignin('base','feed',feed);
assignin('base','dataset',dataset); % Ojo a DATASET
assignin('base','skip_delays',true);

%========================= PROBLEM SPECIFICATIONS=========================

problem.f = 'minSquares';  %Objective function (.m file)

ParamValues =  [0       2.33057   10; ... % 1.- Vmax glucose uptake [mmol/gDCW*h]
                0       .937      1;  ... % 2.- Kg glucose saturation constant [g/L]
                0       0.53802   0.7;  ... % 3.- Glucose-ethanol minimum yield BATCH
                0       0.0697    0.2;  ... % 4.- Glucose-pyruvate minimum yield BATCH
                0       0.1788    0.6;  ... % 5.- Glucose-arabitol minimum yield BATCH
                0       0         0.1;  ... % 6.- Glucose-citrate minimum yield BATCH
                0       0         0.01;  ... % 7.- Glucose-acetate minimum yield BATCH
                0       0         1e-3;  ... % 8.- Glucose-thaumatin minimum yield BATCH
                -1.5   -.672      0;  ... % 9.- Glucose-ethanol minimum yield FEDBATCH
                -1     -.444      0;  ... % 10.- Glucose-pyruvate minimum yield FEDBATCH
                -.5    -.08246    0;  ... % 11.- Glucose-arabitol minimum yield FEDBATCH
                0       0        2;  ... % 12.- Glucose-citrate minimum yield FEDBATCH
                0       0        0.01;  ... % 13.- Glucose-acetate minimum yield FEDBATCH
                0       0        1e-3;  ... % 14.- Glucose-thaumatin minimum yield FEDBATCH
                0       0       1;  ... % 15.- Parámetro 'a' función objetivo BATCH
                0       0       1];% 16.- Parámetro 'a' función objetivo FEDBATCH   
                
                

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
opts.maxeval      = 6000;
opts.local.n1     = 3000;    
opts.maxtime      = 1e5;
opts.strategy     = 1;
opts.local.solver = 'n2fb';
opts.local.finish = 'lsqnonlin';

%========================== DATA EXPERIMENTAL ============================
% Read experimental data

expdata = xlsread(filename,dataset);
% Peso seco
DWRelation = 0.8;
expdata(:,3) = expdata(:,3)*DWRelation;


% Initial conditions
x0 = expdata(1,2:end)';
texp    = expdata(:,1);
ydata   = expdata(:,2:end);
weights = ones(9,1);
weights(2) = 0.05; % privilegiar ajuste biomasa 
weights(3) = 0.08; % Privilegiar tambien ajuste a la glucosa
weights(4) = 0.15; % Etanol
weights(6) = 0.15; % Arabitol

assignin('base','x0',x0);
assignin('base','texp',texp);
assignin('base','ydata',ydata);
assignin('base','weights',weights);
% assignin('base','Vout',expdata(:,1:2));

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
if feedFunction(35,k_SS)==0
    %Batch fermentation, works faster with ode113
    [t,x] = ode113(@pseudoSteadyState,[0 simTime],x0,odeoptions,k_SS);
else
    %Fed-Batch fermentation, works faster with ode15s
    [t,x] = ode15s(@pseudoSteadyState,[0 simTime],x0,odeoptions,k_SS);
end
t2 = toc(time2);

%Save fitting results as a figure
printResults(t,x,expdata)
saveas(gcf,['fitting_d' num2str(dataset) '.fig'])

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
    save(['it_results_d' num2str(dataset) '_pre.mat'],'it_results');
else
    save(['it_results_d' num2str(dataset) '_post.mat'],'it_results');
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