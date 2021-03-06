%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cmp_group = recalculate_CC(dataset)
% Calculates for all solutions of the reparametrization the CCs, and
% discards the ones that have any CC > 2.
%
% INPUTS:
% dataset       Dataset to check (NOTE: the file "it_dXX.mat" must be
%               present in the folder)
%
% OUTPUTS:
% cmp_group     All reparametrizations with no sensitivity,
%               identifiability and significance problems.
%
% Benjam�n J. S�nchez
% Last update: 2016-12-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmp_group = recalculate_CC(filename,dataset)

load(['it_d' num2str(dataset) '.mat']);
assignin('base','it',it);
n = length(it.codes{1,1});

cmp_group.codes = find_best(dataset);
m               = length(cmp_group.codes);
cmp_group.sols  = zeros(m,n);
cmp_group.CC    = zeros(m,1);

%Initialize COBRA
initCobraToolbox

%Define initial variables and constraints

% Change Cobra Solver
changeCobraSolver('gurobi5','LP');
changeCobraSolver('gurobi5','QP');

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
        88.06   152.14  192.124]./1000; % Pesos moleculares de los �cidos no ionizados


cd data
    expdata = xlsread(filename,dataset);
cd ..

% Initial conditions
x0 = expdata(1,2:8);
texp    = expdata(:,1);
ydata   = expdata(:,2:8);

assignin('base','texp',texp);
assignin('base','ydata',ydata);
assignin('base','model',model);
assignin('base','excMet',excMet);
assignin('base','excRxn',excRxn);
assignin('base','x0',x0);
assignin('base','feed',feed);
assignin('base','PM',PM);
assignin('base','skip_delays',false);

for i = 1:m
    %Calculate CC:
    fixed_values = it.codes{1,2}.k_SS;
    kfixed       = it.codes{cmp_group.codes(i),2}.kfixed;
    CC           = reg_analysis_onlyCC(fixed_values,kfixed);
    
    %Discard analysis if any CC is NaN or larger than 2:
    discard = false;
    for j = 1:length(kfixed)
        if isnan(CC(j)) || abs(CC(j)) > 2
            discard = true;
        end
    end
    
    %Update cmp_group:
    cmp_group.sols(i,:) = kfixed;
    
    if discard
        cmp_group.CC(i) = NaN;
    else
        N = 0;
        for j = 1:length(kfixed)
            if isnan(kfixed(j))
                N = N+1;
            end
        end
        cmp_group.CC(i) = sum(CC)/N;
    end
    
    save(['cmp_d' num2str(dataset) '.mat'],'cmp_group');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%