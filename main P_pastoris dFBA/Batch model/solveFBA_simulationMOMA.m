%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FBAsol = solveFBA(model,t,excRxn,p);
% Solves the LP problem of finding the fluxes via FBA. Maximizes biomass
% formation rate and minimizes absolute flux sum (Schuetz 2012).
%
% INPUTS:
% model     COBRA model used in the simulation
% t         Time of simulation [h]
% excRxn    Numeric matrix with the position in rxns of each
%           exchange reaction (the first row (volume) has zeros)
% p         Kinetic parameters
% 
% OUTPUT:
% FBAsol    FBA solution
%
% Benjamín J. Sánchez
% Last Update: 2014-11-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FBAsol = solveFBA_simulationMOMA(model,t,excRxn,p,geneID)

%Parameters used in solveFBA:
w = p(8);

% Minimize sum of absolut fluxes
QPproblem.A             = model.S;
QPproblem.b             = model.b;
N                       = length(model.rxns);
LinearCoeff             = zeros(N,1);
BiomCoeff               = findRxnIDs(model,'BIOMASS');
LinearCoeff(BiomCoeff)  = -(1-w);
QPproblem.F             = w*eye(N);                  %Objective quadratic term
QPproblem.c             = LinearCoeff;              %Objective linear term
QPproblem.ub            = model.ub;
QPproblem.lb            = model.lb;
QPproblem.osense        = +1;                      %minimization

for i=1:length(model.mets) % N° mets
    QPproblem.csense(i,1) = 'E';            %Equality constraints for gnpk
end

FBAsol_min       = solveCobraQP(QPproblem);
%FBAsol_min.full(238)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(geneID)
    try
    modelDel = model;
    [~,results] = findRxnsFromGenes(modelDel,modelDel.genes(geneID),[],1);
    modelDel = changeRxnBounds(modelDel,results(:,1),zeros(size(results(:,1))),'b');
    
    % Redefinimos el problema de balance de flujos
    QPproblem.F             = eye(N);
    QPproblem.ub            = modelDel.ub;
    QPproblem.lb            = modelDel.lb;
    QPproblem.c             = -FBAsol_min.full;
    MOMAsol_min = solveCobraQP(QPproblem);
    FBAsol.x = MOMAsol_min.full;
    
    catch
        FBAsol.x = zeros(size(model.rxns));
    end

else 
    FBAsol.x=FBAsol_min.full;
end

fluxDistrib = [FBAsol.x;t];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metabolic Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if t==0
            save(['metMovie_KO' num2str(geneID) '.mat'],'fluxDistrib')
        else
            if length(fluxDistrib)<1000
                fluxDistrib = [zeros(length(model.rxns),1);fluxDistrib];
            end
            S = load(['metMovie_KO' num2str(geneID) '.mat'],'fluxDistrib');
            fluxDistrib = [S.fluxDistrib,fluxDistrib];
            save(['metMovie_KO' num2str(geneID) '.mat'],'fluxDistrib');
        end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%