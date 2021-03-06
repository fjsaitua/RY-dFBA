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
% Benjam�n J. S�nchez
% Last Update: 2014-11-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = pseudoSteadyState(t,x,k) % cambi� k por p ya que no estoy fijando nada 

model       = evalin('base','model');
excMet      = evalin('base','excMet');
excRxn      = evalin('base','excRxn');
PM          = evalin('base','PM');
feed        = evalin('base','feed');
kfixed      = evalin('base','kfixed');
post_starv  = evalin('base','post_starv');
t_feed      = evalin('base','t_feed');

%Construct p from the fixed values of kfixed and the estimated values of k.
p = zeros(length(kfixed),1);
j = 1;
for i = 1:length(kfixed)
    if isnan(kfixed(i))
        p(i) = k(j); % �?
        j    = j+1;
    else
        p(i) = kfixed(i);
    end    
end


% Fixed modifications of the model
model = fixedModifications(model,p);

%Kinetic constraints
model = kineticConstraints(model,t,x,excRxn,p);

%FBA
FBAsol = solveFBA(model,t,excRxn,p);

if numel(FBAsol.x) == 0
    FBAsol.x = zeros(length(model.rxns),1);
end

% Save the data from the end of the batch
if t >= (t_feed-1) && post_starv == 0
    assignin('base','X0_FB',x(2))
    assignin('base','V0_FB',x(1))
    post_starv = 1;
    assignin('base','post_starv',post_starv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oxygen Demand Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
O2 = FBAsol.x(findRxnIDs(model,'EX_o2(e)'));
O2_instant_demand = abs(O2*32*x(2));
assignin('base','O2_instant_demand',O2_instant_demand);


%Integration:
N  = length(excMet);
dx = zeros(N,1);

% Volume:
SampligRate = 0; % For simulation
Fin   = feedFunction(t);
dx(1) = Fin-SampligRate;              %[L/h]

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