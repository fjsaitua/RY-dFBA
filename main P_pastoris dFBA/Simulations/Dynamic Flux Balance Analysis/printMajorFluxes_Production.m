function [MetProdFlux,SortRxnIDs,SortProdFlux] = printMajorFluxes_Production(model,metName,flux,N)
% Function that prints the N fluxes that produce the most of a certain
% metabolite.

metID = findMetIDs(model,metName);
rxns = findRxnsFromMets(model,metName);
rxnIDs = findRxnIDs(model,rxns);

ProductionFluxes = zeros(length(rxns),2);
j=1;

for i=1:length(rxnIDs)
    % Conditions for production
    if model.S(metID,rxnIDs(i)) > 0 && flux(rxnIDs(i)) > 0 
        % Product and Positive flux 
        ProductionFluxes(j,1) = rxnIDs(i);
        ProductionFluxes(j,2) = flux(rxnIDs(i))*model.S(metID,rxnIDs(i));
        j=j+1;
    elseif model.S(metID,rxnIDs(i)) < 0 && flux(rxnIDs(i)) < 0
        % Substrate and Negative flux        
        ProductionFluxes(j,1) = rxnIDs(i);
        ProductionFluxes(j,2) = flux(rxnIDs(i))*model.S(metID,rxnIDs(i));
        j=j+1;
    end
end

% sort flux vectors
[SortProdFlux,I]=sort(abs(ProductionFluxes(:,2)),'descend');

SortRxnIDs = ProductionFluxes(I,1);
SortProdFlux = ProductionFluxes(I,2);

% for i=1:N
%     if SortRxnIDs(i)~=0
%     disp([model.rxns{SortRxnIDs(i)} ' ' num2str(SortProdFlux(i))]);
%     else
%         break
%     end
% end

MetProdFlux = sum(abs(SortProdFlux));
end