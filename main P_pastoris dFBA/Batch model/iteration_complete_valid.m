%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it_results = iteration_complete_valid(filename,dataset,kfixed,n_struct)
% 
% performs the same routines as iteration_complete.m, but it also receives 
% an input of the number of the model structure being tested in the cross 
% – calibration stage.
%
% INPUTS:
% filename      Name of the file from which data will be extracted (.xlsx
%               file)
% dataset       Sheet of the .xlsx file with the data which will be used 
%               for calibration 
% dataset       Number indicating wich sheet will be analyzed
% kfixed        Vector indicating which parameters are fixed. Complete the
%               rest of the parameters with NaN
%
% OUTPUT:
% it_results    Structure containing all iteration results
%
% Francisco Saitua
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
model = readCbModel('iFS670.xml');

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

   
assignin('base','model',model);
assignin('base','excMet',excMet);
assignin('base','excRxn',excRxn);
assignin('base','PM',PM);
assignin('base','skip_delays',true);

%========================= PROBLEM SPECIFICATIONS=========================

problem.f = 'minSquares';  %Objective function (.m file)

ParamValues =  [1       6           8; ... % 1.- Vmax glucose uptake [mmol/gDCW*h]
                1e-5    2.7e-3      1e-2;  ... % 2.- Kg glucose saturation constant [g/L]
                0       0.4         2;  ... % 3.- Glucose-ethanol minimum yield BATCH
                0       0.117       1;  ... % 4.- Glucose-pyruvate minimum yield BATCH
                0       0.1364      1;  ... % 5.- Glucose-arabitol minimum yield BATCH
                0       0           1;  ... % 6.- Glucose-citrate minimum yield BATCH
                0       0           1e-3;  ... % 7.- Suboptimal grwoth coefficient alpha
                0       2.18      10];    % 8.- Maintenance ATP

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
cd ..

% Initial conditions

x0 = expdata(1,2:8);
texp    = expdata(:,1);
ydata   = expdata(:,2:8);

assignin('base','x0',x0);
assignin('base','texp',texp);
assignin('base','ydata',ydata);


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