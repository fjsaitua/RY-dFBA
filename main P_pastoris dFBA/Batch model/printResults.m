%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% printResults(t,x,expdata)
% Displays the model simulation & experimental results in a 2x2 graph.
%
% INPUTS:
% t             Time vector
% x             Variable vector
% expdata       Experimental data
%
% Benjamín J. Sánchez
% Last Update: 2014-11-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printResults(t,x,expdata,p)

variables = {   'Biomass' 'Glucose' 'Ethanol' 'Pyruvate'...
                'Arabitol' 'Citrate'};

% Gráficos variables de estado

texp    = expdata(:,1);
ydata   = expdata(:,2:end);


for i=1:6
    subplot(3,3,i)
    plot(t,x(:,i+1),'r-',texp(),ydata(:,i+1),'k+')
    title(variables(i))
    ylabel('[g/L]')
    xlabel('Time [h]')    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%