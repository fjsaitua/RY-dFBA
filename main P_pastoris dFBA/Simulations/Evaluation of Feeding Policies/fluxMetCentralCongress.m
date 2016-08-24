function [X_Central,RxnList] = fluxMetCentralCongress(model,x)
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

GlycList = {'EX_glc(e)','GLCt1','HEX1','PGI','PFK',...
            'FBA','TPI','GAPD','ENO','PYK','FBA3','PFK_3','PGMT','TRE6PS'};
        
PPPList  = {'G6PDH2','RPE','RPI','TKT1','TKT2','TALA',...
            'ABTDH','ABTt','ABTDt','AMPN','YUMPS'}; 

PyrFermList  = {'PYRt2m','D-LACt2m','PDHm','PC','OAAt2m',...
                'PYRDC','ALCD2x','ALDD2y','PYRt','PYRt2','ETOHt'};
            
PyrAcCoAList = {'ACS', 'ACACT1', 'ACCOACr', 'ACOATA',...
            'CSNATr', 'FAS80COA_L', 'FAS80_L', 'HMGCOAS','HSERTA'}; % Includes lipid biosynthesis
            
TCAList = {'CSm','ICDHxm','SUCOASm','SUCD2_u6m','FUMm',...
            'MDHm','MDH','MALtm','ASPTAm'};

Ex_rxnList = { 'BIOMASS','EX_etoh(e)', 'EX_pyr(e)','EX_abt_D(e)',...
            'EX_cit(e)','EX_tau(e)','EX_co2(e)'};

% Get reaction indexes
GlycIDs     = findRxnIDs(model,GlycList);
PPPIDs      = findRxnIDs(model,PPPList);
PyrFermIDs  = findRxnIDs(model,PyrFermList);
TCAIDs      = findRxnIDs(model,TCAList);
PyrAcCoAIDs = findRxnIDs(model,PyrAcCoAList);
Ex_rxnIDs   = findRxnIDs(model,Ex_rxnList);


% Get reaction Fluxes
Glyc    = x(GlycIDs,:);
PPP     = x(PPPIDs,:);
PyrFerm = x(PyrFermIDs,:);
TCA     = x(TCAIDs,:);
AcCoA   = x(PyrAcCoAIDs,:);
EXC     = x(Ex_rxnIDs,:);

% Concatenate Fluxes
X_Central = [Glyc;PPP;PyrFerm;TCA;AcCoA;EXC];
RxnList = [GlycList PPPList PyrFermList TCAList PyrAcCoAList Ex_rxnList];

end