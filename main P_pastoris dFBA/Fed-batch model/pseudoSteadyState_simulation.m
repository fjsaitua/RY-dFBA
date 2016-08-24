%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dx = pseudoSteadyState(t,x,k)
% Optimizes biomass using FBA, under the pseudo-steady state assumption
%
% INPUTS:
% t             Time of simulation [h]
% x             Concentrations of each variable in the last iteration 
%               ([L] or [g/L])
% k             Kinetic parameters for glucose consumption
% 
% OUTPUT:
% dx            Derivatives of concentrations [g/Lh]
%
% Benjamín J. Sánchez
% Last Update: 2014-11-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = pseudoSteadyState_simulation(t,x,k) % cambié k por p ya que no estoy fijando nada 

model  = evalin('base','model');
excMet = evalin('base','excMet');
excRxn = evalin('base','excRxn');
PM     = evalin('base','PM');
feed   = evalin('base','feed');
kfixed = evalin('base','kfixed');

%Construct p from the fixed values of kfixed and the estimated values of k.
p = zeros(length(kfixed),1);
j = 1;
for i = 1:length(kfixed)
    if isnan(kfixed(i))
        p(i) = k(j); % ¿?
        j    = j+1;
    else
        p(i) = kfixed(i);
    end    
end

% Fixed modifications of the model
model = fixedModifications(model,p);

%Kinetic constraints
model = kineticConstraints(model,t,x,excRxn,p);

%FBA & MOMA
FBAsol = solveFBA_simulation(model,t,excRxn,p);

if numel(FBAsol.x) == 0
    FBAsol.x = zeros(length(model.rxns),1);
end


%Integration:
N  = length(excMet);
dx = zeros(N,1);

% Volume:
Fin   = feedFunction(t);
dx(1) = Fin;              %[L/h]

% Biomass:
mu    = FBAsol.x(excRxn(2,1));  %[1/h]
dx(2) = mu*x(2)-x(2)/x(1)*Fin;  %[gDW/Lh]

%Extracellular Metabolites:
for q = 3:N
    i     = excMet(q);
    j     = excRxn(q,1);
    v     = -FBAsol.x(j)*model.S(i,j);              %[mmol/gDWh]
    dx(q) = (feed(q)-x(q))*Fin/x(1) + v*PM(q)*x(2); %[g/Lh]
end


skip_delays = evalin('base','skip_delays');
persistent odetime
maxtime = 300;

if skip_delays
    if isempty(odetime)
        odetime = tic;
    elseif toc(odetime) > maxtime
        odetime = [];
        error('Slow integration');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%