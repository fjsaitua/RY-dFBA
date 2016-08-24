%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it = add_it(it,last_results)
% Adds a new iteration to the pending iteration group, only if it is not
% present already in the past group or the pending group. Returns the
% modified pending iterations group, and the modified tree and codes.
%
% Benjamín J. Sánchez
% Last Update: 2014-11-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function it = add_it(it,last_results)

all_it      = [it.past;it.pending];  %Past + pending iterations
[N,~]       = size(all_it);
repeated_it = false;

%Check if iteration is repeated in past_it or pending_it:
for j = 1:N
    R    = it.new - all_it(j,:);
    diff = sum(R.^2);
    if diff == 0
        repeated_it = true;
    end
end

if permited_it(it,last_results)
    if ~repeated_it
        %Add iteration to pending group:
        it.pending = [it.pending;it.new];
        %New iteration position:
        new_pos = N+1;
        %Update iteration codes:
        it.codes{N+1,1} = it.new;
    else
        %Find repeated iteration:
        for i = 1:N
            R    = it.new - it.codes{i,1};
            diff = sum(R.^2);
            if diff ==0
                new_pos = i;
            end
        end
    end
    %Find last iteration and update iteration tree:
    for i = 1:N
        R    = it.last - it.codes{i,1};
        diff = sum(R.^2);
        if diff ==0
            it.tree(i,new_pos) = 1;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%