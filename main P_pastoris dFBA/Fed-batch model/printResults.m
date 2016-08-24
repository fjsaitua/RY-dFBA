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

function printResults(t,x,texp,ydata)

variables = {   'Volume' 'Biomass' 'Glucose' 'Ethanol' 'Pyruvate'...
                'Arabitol' 'Citrate'};

% Gráficos variables de estado

for i=1:7
    subplot(3,3,i)
    plot(t,x(:,i),'r-',texp,ydata(:,i),'k+')
    title(variables(i))
    ylabel('[g/L]')
    xlabel('Time [h]')    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%