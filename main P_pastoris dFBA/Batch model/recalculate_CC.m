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
% Benjamín J. Sánchez
% Last update: 2014-11-29
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
model = readCbModel('PP_iFS618.xml');

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

cd data
    expdata = xlsread(filename,dataset);
    % Peso seco
    DWRelation = 0.72;
    expdata(:,3) = expdata(:,3)*DWRelation;
cd ..

% Initial conditions
x0 = expdata(1,2:8);
weights = expdata(:,9:15);

texp    = expdata(:,1);
ydata   = expdata(:,2:8);

% weights(:,2) = 0.07*weights(:,2); % Privilegiar ajuste biomasa 
% weights(:,3) = 0.13*weights(:,3); % Privilegiar ajuste glucosa
% weights(:,4) = 0.15*weights(:,4); % Privilegiar ajuste etanol
% weights(:,6) = 0.15*weights(:,5); % Privilegiar ajuste arabitol
% weights(:,7) = 0.2*weights(:,6); % Privilegiar ajuste citrato
% weights(:,8) = 0.07*weights(:,7); % Privilegiar ajuste taumatina

assignin('base','texp',texp);
assignin('base','ydata',ydata);
assignin('base','weights',weights);
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