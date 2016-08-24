%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = kineticConstraints(model,x,excRxn,p)
% Calculate relevant kinetic constraints, as LB or UB for the SBML model
%
% INPUTS:
% model     COBRA model used in the simulation
% x         Concentrations of each variable in the last iteration 
%           ([L] or [g/L])
% excRxn    Numeric matrix with the position in rxns of each
%           exchange reaction (the first row (volume) has zeros)
% p         Kinetic parameters
% 
% OUTPUT:
% model     COBRA model used in the simulation (changed)
%
% Benjamín J. Sánchez
% Last Update: 2014-11-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = kineticConstraints(model,t,x,excRxn,k)

t_feed  = evalin('base','t_feed');
ateEtOH = evalin('base','ateEtOH');
atePyr = evalin('base','atePyr');
ateArab = evalin('base','ateArab');

p = k;

%Parameters used in kineticConstraints:
vmax  = p(1);	% Maximum glucose uptake rate [mmol/gDWh]
Kg    = p(2);	% Glucose saturation constant [g/L]
pEtOH = p(3);   % Ethanol specific production rate [fraction of glucose consumption] - BATCH PHASE
pPyr  = p(4);   % Glucose pyruvate minimum yield BATCH
pArab = p(5);   % Glucose arabitol minimum yield BATCH
pCit  = p(6);   % Glucose citrate minimum yield BATCH
pTAU  = p(7);   % Glucose Thaumatin minimum yield BATCH

if feedFunction(t) ~= 0
    vEtOH   = p(8);  % Ethanol specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
    vPyr    = p(9);  % Pyruvate specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
    vArab   = p(10);  % Arabitol specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
    vCit    = p(11);  % Citrate specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
    vTAU    = p(12);  % Thaumatin specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
end


%Glucose uptake
G  = x(3);      % Glucose [g/L]
E  = x(4);      % Ethanol
P  = x(5);      % Pyruvate
A  = x(6);      % Arabitol

if G < 1e-3
    v = 0;
else
    % v = vmax*G/(Kg+G)/(1+E/Ke);   %[mmol/gDWh]
    v = vmax*G/(Kg+G);
end

if E < 1e-3
    vEtOH = 0;
    if ateEtOH == 0 && t > t_feed
        assignin('base','V0_FB',x(1))
        assignin('base','X0_FB',x(2))
        assignin('base','t_feed_carb',t)
        assignin('base','ateEtOH',1)
    end
end

if P < 1e-3 
    vPyr = 0;
    if atePyr == 0 && t > t_feed
        assignin('base','V0_FB',x(1))
        assignin('base','X0_FB',x(2))
        assignin('base','t_feed_carb',t)
        assignin('base','atePyr',1)
    end
end

if A < 1e-3
    vArab = 0;
    if ateArab == 0 && t > t_feed
        assignin('base','V0_FB',x(1))
        assignin('base','X0_FB',x(2))
        assignin('base','t_feed_carb',t)
        assignin('base','ateArab',1)
    end
end

%Assign flux uptake/production rates to model:
model = changeRxnBounds(model,model.rxns(excRxn(3,1)),-v,'l'); % Glucose

if feedFunction(t) == 0
    model = changeRxnBounds(model,model.rxns(excRxn(4,1)),v*pEtOH,'l'); % Proportional to glucose consumption
    model = changeRxnBounds(model,model.rxns(excRxn(5,1)),v*pPyr,'l');
    model = changeRxnBounds(model,model.rxns(excRxn(6,1)),v*pArab,'l'); % Proportional to glucose consumption
    model = changeRxnBounds(model,model.rxns(excRxn(7,1)),v*pCit,'l');
    model = changeRxnBounds(model,model.rxns(excRxn(8,1)),v*pTAU,'l');
else
    instant_mu = evalin('base','instant_mu_set');
    model = changeRxnBounds(model,model.rxns(excRxn(4,1)),-vEtOH,'l'); 
    model = changeRxnBounds(model,model.rxns(excRxn(5,1)),-vPyr,'l');
    model = changeRxnBounds(model,model.rxns(excRxn(6,1)),-vArab,'l'); 
    model = changeRxnBounds(model,model.rxns(excRxn(7,1)),-vCit,'l');
    model = changeRxnBounds(model,model.rxns(excRxn(8,1)),qp_serum(instant_mu),'l');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%