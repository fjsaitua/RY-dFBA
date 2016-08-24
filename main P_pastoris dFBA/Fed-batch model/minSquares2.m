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
weights    = evalin('base','weights');
odeoptions = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',0.7,'NonNegative',1:length(x0));
try
    Time1 = tic;
    post_starv = 0;
    assignin('base','post_starv',post_starv);
    [~,xmod] = ode15s(@pseudoSteadyState,texp,x0,odeoptions,k);
    ymod = xmod(:,1:end);
    display([num2str(toc(Time1)) 'seconds'])
catch exception
    ymod = 1e3*ones(size(yexp));
    disp(exception.identifier);
end

clear pseudoSteadyState

%Return simulation output, appropiately weighted:
% for i = 1:8
%     ymod(:,i) = ymod(:,i)./(max(yexp(:,i))*weights(:,i));
% end

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