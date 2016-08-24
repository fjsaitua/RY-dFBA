clear all
initCobraToolbox
model = readCbModel('PP_iFS618.xml');

%% Unite Simulation outputs

matNames = {'DynamicMOMAdata_1_350_HSA_revised.mat'...
            'DynamicMOMAdata_351_670_HSA_revised.mat'};

dMOMA_DATA = cell(1,671);

for i=1:2
    load(matNames{i})
    if i==1
        dMOMA_DATA(1:350) = DynamicMOMAdata(1:350); 
    else
        dMOMA_DATA(351:end) = DynamicMOMAdata(351:end);
    end
end

save('DynamicMOMAdata_complete.mat','dMOMA_DATA')

final_Thau = zeros(670,1);
final_Biomass = zeros(670,1);

for i=1:671
    final_Biomass(i) = dMOMA_DATA{1,i}(end,3);
    final_Thau(i) = dMOMA_DATA{1,i}(end,9);
end

InitialThauProd = final_Thau(end);

% Isolating and sorting
improved = find(final_Thau > final_Thau(end));
Reduced_Thau = final_Thau(improved);
Reduced_Biomass = final_Biomass(improved);
[sort_Thau,I] = sort(Reduced_Thau,'descend');
improved = improved(I);
sort_Biomass = Reduced_Biomass(I);
genes = model.genes(improved);
sort_genes = genes;

[results,ListResults] = findRxnsFromGenes(model,sort_genes,[],1);

%% Plots
figure(1)
mat = [sort_Biomass, sort_Thau];
%mat([1 4 5],:)=[];
plot(mat(:,1),mat(:,2),'bo','MarkerFaceColor','b')
%semilogy(mat(:,1),mat(:,2),'b.')
hold on
% Nocon = mat([121 181 189 233 253],:);
% plot(Nocon(:,1),Nocon(:,2),'ro','MarkerFaceColor','r')
% hold on
plot([0 18],[InitialThauProd InitialThauProd],'k--','LineWidth',2)

%semilogy(Nocon(:,1),Nocon(:,2),'r^','MarkerFaceColor','r')
xlabel('Final Biomass Concentration [g/L]')
ylabel('Final Thaumatin Concentration [g/L]')


%% XLS write
filename = 'MOMAoutput_curated_HSA.xlsx'; 
xlswrite(filename,ListResults,1,'A1');
xlswrite(filename,improved,2,'A1');
xlswrite(filename,sort_genes,2,'B1');
xlswrite(filename,sort_Biomass,2,'C1');
xlswrite(filename,sort_Thau,2,'D1');
