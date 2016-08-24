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
%prot = p(5);    %Protein fraction [-]
%carb = p(6);    %Carbohidrate fraction [-]
%lipi = p(7);    %Carbohidrate fraction [-]
vATP = p(15);    %Non-Growth ATP maintenance [mmol/gDWh]
%t    = p(9);    %Threshold for expression [# standard deviations below average]
%perc = p(10);   %Arrays that should be below threshold in order to delete gen [%]

%Protein requirements:
%XrxnID   = findRxnIDs(model,'r_2110');  % Biomass Reaction
% protList = {'s_0955'                    % Essential aminoacids participating of the biomass reaction
%             's_0965'
%             's_0969'
%             's_0973'
%             's_0981'
%             's_0991'
%             's_0999'
%             's_1003'
%             's_1006'
%             's_1016'
%             's_1021'
%             's_1025'
%             's_1029'
%             's_1032'
%             's_1035'
%             's_1039'
%             's_1045'
%             's_1048'
%             's_1051'
%             's_1056'};
% protID   = findMetIDs(model,protList);
% Svalues  = zeros(size(protList));
% for i = 1:length(protID)
%     Svalues(i) = model.S(protID(i),XrxnID);     % Changes the S Matrix
% end
% model = changeRxnMets(model,protList,protList,'r_2110',prot*Svalues);
% 
% %Carbohidrate requirements:
% carbList = {'s_0002'
%             's_0773'
%             's_1107'
%             's_1520'};
% carbID   = findMetIDs(model,carbList);
% Svalues  = zeros(size(carbList));
% for i = 1:length(carbID)
%     Svalues(i) = model.S(carbID(i),XrxnID);
% end
% model = changeRxnMets(model,carbList,carbList,'r_2110',carb*Svalues);
% 
% %Lipid requirement:
% lipiList = 's_1096';
% lipiID   = findMetIDs(model,lipiList);
% Svalue   = model.S(lipiID,XrxnID);
% model    = changeRxnMets(model,lipiList,lipiList,'r_2110',lipi*Svalue);

%ATP maintenance requirements:
model = changeRxnBounds(model,'ATPM',vATP,'l');

% 
% %Transcriptomic data:
% tdata = evalin('base','trans');
% [M,N] = size(tdata.values);
% for i = 1:M
%     count = 0;
%     for j = 1:N
%         if tdata.values(i,j) <= tdata.tave(j) - t*tdata.tstd(j)
%             count = count + 1;
%         end
%     end
%     %Delete gen only if a perc% of the arrays have the expression value
%     %t standard deviations below the average:
%     if count/N*100 >= perc
%         model = deleteModelGenes(model,tdata.names(i,1));
%     end
% end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%