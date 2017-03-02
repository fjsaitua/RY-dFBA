%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FBAsol = solveFBA_simulationMOMA(model,t,excRxn,p);
% Solves two sequential QP problems to find the fluxes of a single gene
% deletion mutant using the minimization of metabolic adjustment as
% objective function.
%
% INPUTS:
% model     COBRA model used in the simulation
% t         Time of simulation [h]
% excRxn    Numeric matrix with the position in rxns of each
%           exchange reaction (the first row (volume) has zeros)
% p         Kinetic parameters
% geneID    Index of the gene to be deleted
% 
% OUTPUT:
% FBAsol    FBA solution
%
% Francisco Saitua
% Last Update: 2016-12-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FBAsol = solveFBA_simulationMOMA(model,t,excRxn,p,geneID)

%Parameters used in solveFBA:
w = p(8); % w=alpha

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

% Determine the parental Flux distributions
FBAsol_min       = solveCobraQP(QPproblem);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(geneID)
    try
    % Gene deletions    
    modelDel = model;
    [~,results] = findRxnsFromGenes(modelDel,modelDel.genes(geneID),[],1);
    modelDel = changeRxnBounds(modelDel,results(:,1),zeros(size(results(:,1))),'b');
    
    % We write the MOMA objective function in its matrix form with the
    % model with a single gene deletion.
    
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Metabolic Movie
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if t==0
%             save(['metMovie_KO' num2str(geneID) '.mat'],'fluxDistrib')
%         else
%             if length(fluxDistrib)<1000
%                 fluxDistrib = [zeros(length(model.rxns),1);fluxDistrib];
%             end
%             S = load(['metMovie_KO' num2str(geneID) '.mat'],'fluxDistrib');
%             fluxDistrib = [S.fluxDistrib,fluxDistrib];
%             save(['metMovie_KO' num2str(geneID) '.mat'],'fluxDistrib');
%         end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%