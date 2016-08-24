%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reparam_dFBA(dataset,size_k)
% Main function parfor reparametrization of the dFBA model.
%
% INPUTS:
% dataset   Dataset to reparametrize (NOTE: the file "it_results_dXX.mat"
%           must be present in the folder)
% size_k    Number of parameters to study
%
% Benjamín J. Sánchez
% Last update: 2014-11-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function reparam_dFBA(dataset,size_k)

id = ['_d' num2str(dataset)];
assignin('base','dataset',dataset);

name_it = ['it' id '.mat'];

%Initialize procedure with first iteration:
last_results  = load(['it_results' id '_pre.mat'],'it_results');
last_results  = last_results.it_results;
fixed_values  = last_results.k_SS;
it.last       = zeros(1,size_k);	%Last iteration performed
it.past       = zeros(1,size_k);    %All past iterations (including it.last)
it.pending    = [];                 %All pending iterations
it.tree       = zeros(1,1);         %Iteration tree
it.codes      = cell(1,2);          %Position in it_tree - iterations - results
it.codes{1,1} = it.last;            %First iteration has code #1
it.codes{1,2} = last_results;       %First iteration results
it.remaining  = true;               %True if an iteration still remains
    
while it.remaining

    %1. Create new iteration vectors from previous one:
    for j = 1:size_k
        if last_results.ktofix(j) == 1
            it.new    = it.last;
            it.new(j) = 1;
            it        = add_it(it,last_results);
        end
    end

    %2. Analize a remaining iteration:
    if isempty(it.pending)
        it.remaining = false;   %Finish HIPPO
    else
        %Get a pending iteration:
        it.next = it.pending(1,:);  %OBS: Tree is analized from a top down approach

        %Construct kfixed:
        kfixed = NaN(1,size_k);
        name   = 'kfixed = [ ';
        for j = 1:length(it.next)
            if it.next(j) == 1
                kfixed(j) = fixed_values(j);
            end
            name = [name num2str(kfixed(j)) ' '];
        end
        disp([name ']'])
        
        %Perform new iteration:
        assignin('base','it',it);
        last_results = iteration_noSS_noCC(fixed_values,kfixed);

        %Update it.last, it.past and it.pending
        it.last = it.next;
        it.past = [it.past;it.next];
        it.pending(1,:) = [];

        %Find iteration in it.codes and save the iteration results
        for j = 1:length(it.codes(:,1))
            R    = it.next - it.codes{j,1};
            diff = sum(R.^2);
            if diff ==0
                it.codes{j,2} = last_results;
            end
        end
    end   
    save(name_it,'it','-v7.3');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%