%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% best_group = find_best(dataset)
% Find all solutions with no identification/sensitivity problems.
%
% INPUTS:
% dataset       Dataset to check (NOTE: the file "it_d[i].mat" must be
%               present in the folder)
%
% OUTPUTS:
% best_group    Positions in "it_d[i].mat" of the solutions free of
%               problems.
%
% Benjamín J. Sánchez
% Last Update: 2014-11-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function best_group = find_best(dataset)

load(['it_d' num2str(dataset) '.mat'])

best_group = [];
for i = 1:length(it.codes)
    Mc     = it.codes{i,2}.Mc;
    Ms     = it.codes{i,2}.Ms;
    kfixed = it.codes{i,2}.kfixed;
    assignin('base','kfixed',kfixed);
    if sum(decision(Mc,Ms)) == 0
        best_group = [best_group;i];
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%