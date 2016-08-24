%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Mc,Ms] = reg_analysis_noSS_noCC
% Updates Ms and Mc matrices, removing columns and rows of the original one
% that are fixed in this iteration.
%
% Benjamín J. Sánchez
% Last Update: 2014-11-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mc,Ms] = reg_analysis_noSS_noCC

dataset = evalin('base','dataset');
kfixed  = evalin('base','kfixed');

%Identifiability -Sensitivity analysis:
it_ori  = load(['it_results_d' num2str(dataset) '_pre.mat']);
it_ori  = it_ori.it_results;
Mc      = it_ori.Mc;
Ms      = it_ori.Ms;

fixed = find(~isnan(kfixed));
Mc(fixed,:) = [];
Mc(:,fixed) = [];
Ms(fixed,:) = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%