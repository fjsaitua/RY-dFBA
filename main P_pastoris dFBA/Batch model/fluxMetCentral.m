function [X_Central,X_O2,X_CO2,RxnList,O2List,CO2List] = fluxMetCentral(model,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function that recovers the central carbon metabolism fluxes of a flux
% vector (result of a cobra optimization). The script is only useful for
% the Pichia pastoris model developed by Chung et al. in 2009.
%
% In order to run this script, COBRA must be initialized.
%
% Input:
%   x = flux matrix (n x i, where n is the number of reactions of the model
%   and i the number of time points considered) 
% 
% Output:
%   X_Central = flux vector of the central metabolism. It only includes the
%               fluxes of the Glycolysis/Gluconeogenesis, PPP, TCA,
%               fermentative pathways, O2 and CO2 metabolism
%
% Last update: Francisco Saitua 30-09-15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of reaction indentificators

GlycList = {'EX_glc(e)','GLCt1','GLUK','G6PI','G6PI3','HEX1','PGI','PFK',...
            'FBA','TPI','GAPD','PGK','PGM','ENO','PYK','FBA3','PFK_3'};
        
PPPList  = {'G6PDH2','PGL','GND','RPE','RPI','TKT1','TKT2','TALA','RIBK',...
            'ABTDH','ABTt','ABTDt'}; 

PyrFermList  = {'PYRt2m','D-LACt2m','PDHm','PC','OAAt2m',...
                'PYRDC','ALCD2x','ALCD2if','ALCD2ir','ALDD2y',...
                'PYRt','PYRt2','ETOHt'};
            
TCAList = {'CSm','ACONTm','ICDHxm','ICDHym','AKGDam','AKGDbm','SUCOASm'...
           'SUCD2_u6m','FUMm','MDHm','MDH','MALtm'};

O2List = findRxnsFromMets(model,'o2[c]');
O2List = ['EX_o2(e)' O2List'];
O2mList = findRxnsFromMets(model,'o2[m]');
O2List = [O2List O2mList'];

CO2List = findRxnsFromMets(model,'co2[c]');
CO2List = CO2List';
       
% Get reaction indexes
GlycIDs     = findRxnIDs(model,GlycList);
PPPIDs      = findRxnIDs(model,PPPList);
PyrFermIDs  = findRxnIDs(model,PyrFermList);
TCAIDs      = findRxnIDs(model,TCAList);
O2IDs       = findRxnIDs(model,O2List);
CO2IDs      = findRxnIDs(model,CO2List);

% Get reaction Fluxes
Glyc    = x(GlycIDs,:);
PPP     = x(PPPIDs,:);
PyrFerm = x(PyrFermIDs,:);
TCA     = x(TCAIDs,:);
X_O2      = x(O2IDs,:);
X_CO2     = x(CO2IDs,:);

% Concatenate Fluxes
X_Central = [Glyc;PPP;PyrFerm;TCA];


RxnList = [GlycList PPPList PyrFermList TCAList];

end