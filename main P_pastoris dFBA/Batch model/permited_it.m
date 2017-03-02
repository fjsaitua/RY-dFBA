%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% answer = permited_it(it,last_results)
% Returns a boolean if the proposed iteration is allowed or not.
% Customizable for each dynamic model.
%
% Francisco J. Saitua
% Last Update: 2016-12-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function answer = permited_it(it,last_results)

answer = true;

% Avoid iterations with all the parameters fixed. This could be used to
% reduce the number of modeling structures explored by HIPPO. For instance,
% one can set that the algorithm should not explore structures that fix all
% the parameters associated to one of the state variables.

if it.new(1) + it.new(2) + it.new(3) + it.new(4) + it.new(5) + ...
   it.new(6) + it.new(7) + it.new(8) == 8
    answer = false;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%