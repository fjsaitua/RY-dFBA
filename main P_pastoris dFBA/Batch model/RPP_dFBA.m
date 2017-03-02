%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RY_dFBA(dataset)
% Main function to use, performs a complete analysis of the model for a
% given dataset.
%
% INPUTS:
% dataset       Number indicating wich sheet will be analyzed 
%
% Benjamín J. Sánchez
% Last Update: 2016 - 12 - 22 Francisco Saitua
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RPP_dFBA(filename,dataset)

profile on % Start Profiler
%(1) Perform a parameter estimation and pre/post regression analysis with
%    no parameters fixed. If everything goes ok, a file called
%    "it_results_d[i]_pre.mat" should appear in the main folder, with all
%    results from the iteration.

% Parameter number
K = 8;

time1 = tic;
iteration_complete(filename,dataset,NaN(1,K));
t1 = toc(time1);

%(2) Perform the reparametrization procedure, generating all possible
%    solutions with no sensitivity/identifiability problems. If everything
%    goes ok, a file called "it_d[i].mat" should appear in the main folder,
%    with all results from the reparametrization.

time2 = tic;
reparam_dFBA(dataset,K);
t2 = toc(time2);

%(3) Calculate significance for all aforementioned combinations, and
%    discards the ones that have an insignificant parameter. If everything
%    goes ok, a file called "cmp_d[i].mat" should appear in the main
%    folder, with all the CC's of each solution.

time3 = tic;
cmp_group = recalculate_CC(filename,dataset);
m         = length(cmp_group.codes);
del_pos   = zeros(1,m);
for i = 1:m
    if isnan(cmp_group.CC(i))
        del_pos(i) = 1;
    end
end

cmp_group.CC(find(del_pos),:)    = [];
cmp_group.codes(find(del_pos),:) = [];
cmp_group.sols(find(del_pos),:)  = [];
t3 = toc(time3);

%(4) Repeat the parameter estimation and the pre/post regression analysis
%    on the best solution found. If everything goes ok, a file called
%    "it_results_d[i]_post.mat" should appear in the main folder, with all
%    results from the iteration.

time4 = tic;
best = 1;
for i = 2:length(cmp_group.codes)
    if cmp_group.CC(i) < cmp_group.CC(best)
        best = i;
    end
end

iteration_complete(filename,dataset,cmp_group.sols(best,:));

t4 = toc(time4);

mainTimes = [t1,t2,t3,t4];
assignin('base','mainTimes',mainTimes);

profile viewer 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%