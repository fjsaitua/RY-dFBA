function [ylimit,Fluxes,RxnNames] = StackedBarPlot(Results,n_TopFlux,model,ylimit)
% Selection of the reactions that carry most of the flux

[~,times] = size(Results);
RxnIdx = cell(1,times);
metFluxes = cell(1,times);
NeedRest = 0;

for i=1:times
    RxnIdxAux = Results{1,i}{1,2}; % Reactions with the top fluxes
    metFluxesAux = Results{1,i}{1,3};
    RxnIdx{1,i} = RxnIdxAux(find(RxnIdxAux));
    metFluxes{1,i} = metFluxesAux(find(metFluxesAux));
end


FluxToConsider = [];
% Make Group of fluxes that I'll include in the Graph

for i=1:times
    IdList = RxnIdx{1,i};
    if length(IdList) < n_TopFlux
        if i==1
            FluxToConsider = IdList;
        else
            newIdx = find(~ismember(IdList,FluxToConsider));
            FluxToConsider = [FluxToConsider;IdList(newIdx)];
        end
    else % 5 or more fluxes involves
        if length(IdList)>n_TopFlux
            NeedRest = 1;
        end
        if i==1
            FluxToConsider = IdList(1:n_TopFlux);
        else
            newIdx = find(~ismember(IdList(1:n_TopFlux),FluxToConsider));
            FluxToConsider = [FluxToConsider;IdList(newIdx)];
        end
    end
end

% Get the Fluxes
Fluxes = zeros(length(FluxToConsider),times);
Type3 = zeros(length(FluxToConsider),1);
Rest   = zeros(1,times); %%%%%% OJO AL RESTO %%%%%%

for j=1:times
    for i=1:length(FluxToConsider)
        if sum(ismember(RxnIdx{1,j},FluxToConsider(i)))>0
            
            idx = find(ismember(RxnIdx{1,j},FluxToConsider(i)));
            FluxVector = Results{1,j}{1,3};
            Fluxes(i,j)= abs(FluxVector(idx));
            
            if Fluxes(i,j)>100 %Manage Type 3 Pathways
                Fluxes(i,j) = 1000 - Fluxes(i,j); 
                Type3(i) = 1; 
            end
        end
    end 
end

% Get RxnIDs
RxnNames = model.rxns(FluxToConsider);
for i=1:length(Type3)
    if Type3(i) == 1;
        RxnNames{i} = [RxnNames{i} '*'];
    end
end

% Rest of the Fluxes
if NeedRest
    for i=1:times
        idx_notMember = find(~ismember(RxnIdx{1,i},FluxToConsider));
        Rest(i) = sum(abs(Results{1,i}{1,3}(idx_notMember)));
    end
    Fluxes = [Fluxes; Rest];
    RxnNames = [RxnNames; 'Rest'];    
end

totFluxes = sum(Fluxes);
ylimit = 1.1*max(totFluxes); % 10% more of the major total flux evaluated at different times
Fluxes = Fluxes';

end