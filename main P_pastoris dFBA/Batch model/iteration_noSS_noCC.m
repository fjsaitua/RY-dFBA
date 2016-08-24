%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it_results = iteration_noSS_noCC(fixed_values,kfixed)
% Does an iteration of the procedure, but skips the parameter estimation
% and pre/post regression diagnostics, only deciding which parameters to
% fix according to sensitivity and identifiability.
%
% INPUTS:
% fixed_values      Values from the first iteration
% kfixed            Vector indicating which parameters are fixed
%                   Complete the rest of the parameters with NaN
%
% OUTPUT:
% it_results        Structure containing all iteration results
%
% Benjamín J. Sánchez
% Last Update: 2014-11-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function it_results = iteration_noSS_noCC(fixed_values,kfixed)

%Decide the parameters to be estimated depending on kfixed:
assignin('base','kfixed',kfixed);
k = fixed_values;
for i = 1:length(kfixed)
    if ~isnan(kfixed(i))
        k(i) = 0;
    end
end

%Regression analysis:
[Mc,Ms] = reg_analysis_noSS_noCC;

%Decision:
ktofix = decision(Mc,Ms);

%Save results
it_results.kfixed     = kfixed;
it_results.k_SS       = k;
it_results.Mc         = Mc;
it_results.Ms         = Ms;
it_results.ktofix     = ktofix;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%