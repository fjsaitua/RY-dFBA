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
t_cons = p(14); % Time where secondary products start to be consumed

if t >= t_cons
    % FED BATCH   
    vEtOH   = p(7);  % Ethanol specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
    vPyr    = p(8);  % Pyruvate specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
    vArab   = p(9);  % Arabitol specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
    vCit    = p(10);  % Citrate specific production rate [fraction of glucose consumption] - FEDBATCH PHASE
end



%Glucose uptake
G  = x(3);      % Glucose [g/L]
E  = x(4);      % Ethanol
P  = x(5);      % Pyruvate
A  = x(6);      % Arabitol
C  = x(7);       % Citrate

if G < 1e-3
    v = 0;
else
    v = vmax*G/(Kg+G);
end

if E < 1e-3
    vEtOH = 0;
end

if P < 1e-3
    vPyr = 0;
end

if A < 1e-3
    vArab = 0;
end

if C < 1e-3
    pCit = 0;
    vCit = 0;
end
    

%Assign flux uptake/production rates to model:
model = changeRxnBounds(model,model.rxns(excRxn(3,1)),-v,'l'); % Glucose

if t < t_cons

    model = changeRxnBounds(model,model.rxns(excRxn(4,1)),pEtOH,'l'); % Proportional to glucose consumption
    model = changeRxnBounds(model,model.rxns(excRxn(5,1)),pPyr,'l');
    model = changeRxnBounds(model,model.rxns(excRxn(6,1)),pArab,'l'); % Proportional to glucose consumption
    model = changeRxnBounds(model,model.rxns(excRxn(7,1)),-pCit,'l');

else
    model = changeRxnBounds(model,model.rxns(excRxn(4,1)),-vEtOH,'l'); % Proportional to glucose consumption
    model = changeRxnBounds(model,model.rxns(excRxn(5,1)),-vPyr,'l');
    model = changeRxnBounds(model,model.rxns(excRxn(6,1)),-vArab,'l'); % Proportional to glucose consumption
    model = changeRxnBounds(model,model.rxns(excRxn(7,1)),vCit,'l');
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%