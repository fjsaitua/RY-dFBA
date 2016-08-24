%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = fixedModifications(model,p)
% Modifies model according to the parameters values. This modifications
% will be fixed during the whole integration.
%
% INPUTS:
% model     SBML model used in the simulation
% p         Kinetic parameters
% 
% OUTPUTS:
% model     SBML model used in the simulation (changed)
%
% Benjamín J. Sánchez
% Last Update: 2014-11-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = fixedModifications(model,p)

%Parameters used in fixedModifications:
vATP =  p(8);    % NGAM

%ATP maintenance requirements:
model = changeRxnBounds(model,'ATPM',vATP,'l');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%