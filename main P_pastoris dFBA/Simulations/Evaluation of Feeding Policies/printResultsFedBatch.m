%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% printResults(t,x,expdata)
% Displays the model simulation & experimental results in a 2x2 graph.
%
% INPUTS:
% t             Time vector
% x             Variable vector
% expdata       Experimental data
%
% Benjam�n J. S�nchez
% Last Update: 2014-11-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printResultsFedBatch(t,x,expdata)

variables = {   'Volume' 'Biomass' 'Glucose' 'Ethanol' 'Pyruvate'...
                'Arabitol' 'Citrate' 'Acetate' 'Thaumatin'};

% Gr�ficos variables de estado

texp    = expdata(:,1);
ydata   = expdata(:,2:end);

for i=1:9
    subplot(3,3,i)
    plot(t,x(:,i),'r-',texp(),ydata(:,i),'k+')
    title(variables(i))
    ylabel('[g/L]')
    xlabel('Time [h]')    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%