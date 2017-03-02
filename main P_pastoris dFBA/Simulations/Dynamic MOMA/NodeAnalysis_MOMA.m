clear all
initCobraToolbox
model_PP = readCbModel('iFS670.xml');

KOgenes = [ 308 133 646 668 669 348 298 571 ...
            572 455 599 67 101 167 338 533 ...
            632 460 435 266 582 652 389 20 ...
            144 179 340 463 314 188 257 433];


cd metMoviesCurated_HSA

FluxDistributions = zeros(length(model_PP.rxns),length(KOgenes)+1);

for i=1:length(KOgenes)+1
    if i==1
        load('metMovie_parental_final.mat')
        FluxDistributions(:,i) = fluxDistrib(1:end-1,1);
    else
        load(['metMovie_KO_' num2str(KOgenes(i-1)) '.mat'])
        FluxDistributions(:,i) = fluxDistrib(1:end-1,1);
    end
end

cd ..

load('DynamicMOMAdata_complete.mat')
DynamicMOMAdata = dMOMA_DATA;

nKOs = length(KOgenes)+1; % Parental + (nKOs-1) knock outs
times = 1:nKOs;

metList = { 'ala-L[c]' 'arg-L[c]' 'asn-L[c]' 'asp-L[c]' 'cys-L[c]' 'glu-L[c]' ...
            'gln-L[c]' 'gly[c]' 'his-L[c]' 'ile-L[c]' 'leu-L[c]' 'lys-L[c]' 'met-L[c]' ...
            'phe-L[c]' 'pro-L[c]' 'ser-L[c]' 'trp-L[c]' 'tyr-L[c]' 'thr-L[c]' ...
            'val-L[c]'};
        

finalHSA = zeros(nKOs,1);
finalBiom = zeros(nKOs,1);


for i=1:nKOs
    if i==1
        finalBiom(i) = DynamicMOMAdata{1,end}(end,3);
        finalHSA(i) = DynamicMOMAdata{1,end}(end,9);
    else
        finalBiom(i) = DynamicMOMAdata{1,KOgenes(i-1)}(end,3);
        finalHSA(i) = DynamicMOMAdata{1,KOgenes(i-1)}(end,9);
    end
end


Results_Prod = cell(length(metList),length(times));
Results_Cons = cell(length(metList),length(times));

Fluxes_P = cell(length(metList),1);
Fluxes_C = cell(length(metList),1);
RxnIDs_P = cell(length(metList),1);
RxnIDs_C = cell(length(metList),1);

for i=1:length(times)
    
        N = 40; % los primeros 40 flujos
        FLUX = FluxDistributions(:,times(i));
        
    for j=1:length(metList)
        
        metName = metList(j);
        
        % Production
        [Totflux_P,SortRxnIDs_P,SortProdFlux_P] = printMajorFluxes_Production(model_PP,metName,FLUX,N);
        % Consumption
        [Totflux_C,SortRxnIDs_C,SortProdFlux_C] = printMajorFluxes_Consumption(model_PP,metName,FLUX,N);
               
        % Results
        Results_Prod{j,i} = {Totflux_P SortRxnIDs_P SortProdFlux_P};
        Results_Cons{j,i} = {Totflux_C SortRxnIDs_C SortProdFlux_C};
        
      
    end
end
%% Stacked Bar Plot
n_TopFlux = 7; % Consider only the 5 reactions that carry more flux in the network

% Make labels for the plot
labels = cell(1,nKOs+2);
labels(1) = {'Parent'};
labels(end) = {' '};

for i = 1:nKOs-1
    labels{i+1} = num2str(KOgenes(i));
end

Productions = zeros(length(KOgenes)+1,length(metList));

for i=1:length(metList)
    
    % Get Data For Plots
    [ylimit_P,Fluxes_P{i,1},RxnIDs_P{i,1}]    = StackedBarPlot(Results_Prod(i,:),n_TopFlux,model_PP); % Production
    [~,Fluxes_C{i,1},RxnIDs_C{i,1}]           = StackedBarPlot(Results_Cons(i,:),n_TopFlux,model_PP); % Consumption
    
    Productions(:,i) = sum(Fluxes_P{i,1},2);
    
    % Production Plot
    figure(i)
    subplot(1,2,1)
    bar(Fluxes_P{i,1},'stacked'), legend(RxnIDs_P{i,1})
    ylim([0 ylimit_P]), ylabel('mmol/g_D_C_Wh')
    xlim([0 length(times)+1])
    set(gca,'XTickLabel',labels)
    title (['Production of ' metList(i)])
    
    % Consumption Plot
    subplot(1,2,2)
    bar(Fluxes_C{i,1},'stacked'), legend(RxnIDs_C{i,1})
    ylim([0 ylimit_P])
    title (['Consumption of ' metList(i)])
    xlim([0 length(times)+1])
    set(gca,'XTickLabel',labels)
end


figure
subplot(1,2,1)
title('Final Biomass')
plot(1:nKOs,finalBiom,'ob')
ylabel('g/L')
set(gca,'XTickLabel',labels)
ylim([0 22])
subplot(1,2,2)
title('Final HSA')
plot(1:nKOs,finalHSA,'ok')
ylabel('g/L')
set(gca,'XTickLabel',labels)
ylim([0 1])

figure
for i=1:length(Productions(1,:))
    
    Productions(:,i) = Productions(:,i)/Productions(1,i);
 
end

Productions(1,:)=[];
plot([0 21],[1 1],'r-')
hold on
boxplot(Productions,metList,'PlotStyle','compact')
ylabel('KO production flux relative to parental strain')
xlim([0 21])
set(gca,'XTick',0:22)

% 
% % Write node information in excel
% for i=1:length(metList)
%     
%         filename = strcat(metIDs{i},' Balance MOMA.xlsx');
%         Data_P = Results_Prod(i,:);
%         Data_C = Results_Cons(i,:);
%         writePieChartData(filename,Data_P,times,model_PP,3)
%         writePieChartData(filename,Data_C,times,model_PP,length(times)+3) 
%         
%         % Write graph information in the first page two sheets
%         % Production
%         xlswrite(filename,RxnIDs_P{i,1},1,'A1')
%         xlswrite(filename,Fluxes_P{i,1}',1,'B1')
%         % Consumption
%         xlswrite(filename,RxnIDs_C{i,1},2,'A1')
%         xlswrite(filename,Fluxes_C{i,1}',2,'B1')
%                 
% end

