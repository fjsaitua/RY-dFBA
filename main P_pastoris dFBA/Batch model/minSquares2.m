%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ymod = minSquares2(k,texp)
% Integrates the dynamic model and returns the simulation output, weighted
% appropiately. To be used in the lsqcurvefit function.
%
% Benjamín J. Sánchez
% Last Update: 2014-11-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ymod = minSquares2(k,texp)

%Integrate:
x0         = evalin('base','x0');
yexp       = evalin('base','ydata');
odeoptions = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',0.7,'NonNegative',1:length(x0));
try
    [~,xmod] = ode113(@pseudoSteadyState,texp,x0,odeoptions,k);
    ymod = xmod(:,1:end);
catch exception
    ymod = 1e3*ones(size(yexp));
    disp(exception.identifier);
end

clear pseudoSteadyState

[m,n] = size(ymod);
for i = 1:m
    for j = 1:n
        if isnan(ymod(i,j))
            ymod(i,j) = 0;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%