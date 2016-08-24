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

function model = kineticConstraints(model,t,x,excRxn,p)

%Parameters used in kineticConstraints:
vmax  = p(1);	% Maximum glucose uptake rate [mmol/gDWh]
Kg    = p(2);	% Glucose saturation constant [g/L]
pEtOH = p(3);   % Ethanol specific production rate [fraction of glucose consumption] - BATCH PHASE
pPyr  = p(4);   % Glucose pyruvate minimum yield BATCH
pArab = p(5);   % Glucose arabitol minimum yield BATCH
pCit  = p(6);   % Glucose citrate minimum yield BATCH

%Glucose uptake, with inhibition by ethanol
G  = x(3);      % Glucose [g/L]
Cit = x(7);

if G < 1e-3
    v = 0;
else
    % v = vmax*G/(Kg+G)/(1+E/Ke);   %[mmol/gDWh]
    v = vmax*G/(Kg+G);
end

%Assign flux uptake/production rates to model:
model = changeRxnBounds(model,model.rxns(excRxn(3,1)),-v,'l'); % Glucose
model = changeRxnBounds(model,model.rxns(excRxn(4,1)),pEtOH,'l'); % Proportional to glucose consumption
model = changeRxnBounds(model,model.rxns(excRxn(5,1)),pPyr,'l');
model = changeRxnBounds(model,model.rxns(excRxn(6,1)),pArab,'l'); % Proportional to glucose consumption

if Cit > 1e-3 % Use only when citrate is being consumed
    model = changeRxnBounds(model,model.rxns(excRxn(7,1)),-pCit,'l');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%