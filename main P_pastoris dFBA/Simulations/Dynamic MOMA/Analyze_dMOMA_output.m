clear all
initCobraToolbox
model = readCbModel('iFS670.xml');

%% Unite Simulation outputs

% Matrices containing the state variables at the end of simulated batch
% cultivations of single gene deletion mutantes

matName = {'KO_info.mat'};

dMOMA_DATA = load(matName);

final_Prot = zeros(670,1);
final_Biomass = zeros(670,1);

for i=1:671
    final_Biomass(i) = dMOMA_DATA{1,i}(end,3);
    final_Prot(i) = dMOMA_DATA{1,i}(end,9);
end

InitialProtProd = final_Prot(end);

% Isolating and sorting
improved = find(final_Prot > final_Prot(end));
Reduced_Prot = final_Prot(improved);
Reduced_Biomass = final_Biomass(improved);
[sort_Prot,I] = sort(Reduced_Prot,'descend');
improved = improved(I);
sort_Biomass = Reduced_Biomass(I);
genes = model.genes(improved);
sort_genes = genes;

[results,ListResults] = findRxnsFromGenes(model,sort_genes,[],1);

%% Plots
figure(1)
mat = [sort_Biomass, sort_Prot];
plot(mat(:,1),mat(:,2),'bo','MarkerFaceColor','b')
hold on
plot([0 18],[InitialProtProd InitialProtProd],'k--','LineWidth',2)

xlabel('Final Biomass Concentration [g/L]')
ylabel('Final Protein Concentration [g/L]')


%% XLS write
filename = 'MOMAoutput_curated_HSA.xlsx'; 
xlswrite(filename,ListResults,1,'A1');
xlswrite(filename,improved,2,'A1');
xlswrite(filename,sort_genes,2,'B1');
xlswrite(filename,sort_Biomass,2,'C1');
xlswrite(filename,sort_Protm,2,'D1');
