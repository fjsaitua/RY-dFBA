%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% answer = permited_it(it,last_results)
% Returns a boolean if the proposed iteration is allowed or not.
% Customizable for each dynamic model.
%
% Benjamín J. Sánchez
% Last Update: 2014-11-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function answer = permited_it(it,last_results)

answer = true;

%Avoid iterations with all parameters in the biomass dynamic equation
%fixed:
if it.new(1) + it.new(2) + it.new(3) + it.new(4) + it.new(5) + ...
   it.new(6) + it.new(7) + it.new(8) + it.new(9) + it.new(10) + ...
   it.new(11) + it.new(12) + it.new(13) + it.new(14) == 14
    answer = false;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%