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

function FBAsol = solveFBA_simulation(model,t,excRxn,p)

%Parameters used in solveFBA:
a = p(9);
b = p(10);
Tsubopt = p(11);
w = Obj_weight(t,a,b,Tsubopt);
% c = p(4); % Coeficiente en la FO de proteína

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

% % Generates error if FBA does not find a feasible solution:
% if FBAsol_min.stat~=1
%    disp('ERROR! FBA is not converging')
% else
%    disp(['Iteration achieved at t = ',num2str(t),' h'])
% end

%Returns the QP solution:
% FBAsol   = FBAsol_X;
% FBAsol.x = FBAsol_min.full;
FBAsol.x = FBAsol_min.full;
fluxDistrib = [FBAsol.x;t];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metabolic Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if t==0
            save('metMovie.mat','fluxDistrib')
        else
            if length(fluxDistrib)<1000
                fluxDistrib = [zeros(length(model.rxns),1);fluxDistrib];
            end
            S = load('metMovie.mat','fluxDistrib');
            fluxDistrib = [S.fluxDistrib,fluxDistrib];
            save('metMovie.mat','fluxDistrib');
        end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%