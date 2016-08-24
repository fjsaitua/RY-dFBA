%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [J,g,R] = minSquares(k,texp,yexp)
% Integrates the dynamic model and calculates afterwards the cuadratic 
% difference between the model predictions and the experimental data.
% Returns the cuadratic difference (objective function). To be used in the
% parameter estimation with SSm.
%
% Benjamín J. Sánchez
% Last Update: 2014-11-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J,g,R] = minSquares(k,texp,yexp)

%Integrate:
x0         = evalin('base','x0');
weights    = evalin('base','weights');
odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',3,'NonNegative',1:length(x0));


try
    post_starv = 0;
    assignin('base','post_starv',post_starv)
    [~,xmod] = ode15s(@pseudoSteadyState,texp,x0,odeoptions,k);
    ymod = xmod(:,1:end);
    
catch exception
    ymod = 1e3*ones(size(yexp));
    disp(exception.identifier);
    
end

clear pseudoSteadyState

%Define optimization function (difference between experimental and model data, normalized by maximum measure):
[mmod,nmod] = size(ymod);
[mexp,nexp] = size(yexp);

% Checks the dimensions of the experimental and model data, apart it checks
% that the model doesn't give a ¿? constant solution // Cambié nmod por
% mmod en el tercer término de la condición

if mmod == mexp && nmod == nexp && sum(abs(ymod(1,:)-ymod(mmod,:))) ~= 0
    R = ymod-yexp;
else
    R = 1e3*ones(size(yexp));
end


for i = 1:7 % Number of state variables (not counting thaumatin)
    R(:,i) = R(:,i)./(max(yexp(:,i))*weights(:,i));
end

[m,n] = size(R);
for i = 1:m
    for j = 1:n
        if isnan(R(i,j))
            R(i,j) = 0;
        end
    end
end

J = sum(sum(R.^2));      %For visualization purposes, J is without semi-colon
R = reshape(R,numel(R),1);
g = 0;

%Save in checkpoint.mat all the solutions (and the best so far):
if exist(fullfile(cd, 'checkpoint.mat'), 'file') == 2
    load('checkpoint.mat','k_tot','J_tot','k_best','J_best');
    k_tot = [k_tot;k];
    J_tot = [J_tot;J];
    if J < J_best
        k_best = k';
        J_best = J;
    end
else
    k_tot = k;
    J_tot = J;
    k_best = k';
    J_best = J;
end
save('checkpoint.mat','k_tot','J_tot','k_best','J_best');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%