clear all
initCobraToolbox
model_PP = readCbModel('PP_iFS618.xml');
model_PP = modelConstraints(model_PP);
indexes = findRxnIDs(model_PP,{'EX_o2(e)' 'BIOMASS'});
load('metMovie_FB_valid.mat')

times = [20 28 35 45];

metList = {'g6p[c]' 'pyr[c]' 'atp[c]' 'nadh[c]' 'nadph[c]'};
metIDs  = {'glucose 6 phosphate cyt balance' 'pyruvate cit balance' 'atp cyt balance' 'nadh cyt balance' 'nadph cyt balance'};

Results_Prod = cell(length(metList),length(times));
Results_Cons = cell(length(metList),length(times));

Fluxes_P = cell(length(metList),1);
Fluxes_C = cell(length(metList),1);
RxnIDs_P = cell(length(metList),1);
RxnIDs_C = cell(length(metList),1);

for i=1:length(times)
    
        time_idx = find(fluxDistrib(end,:) >= times(i));
        time_idx = time_idx(1);
        N = 40; % los primeros 40 flujos
        FLUX = fluxDistrib(:,time_idx);
        
    for j=1:length(metList)
        
        metName = metList(j);
        
        % Production
        [Totflux_P,SortRxnIDs_P,SortProdFlux_P] = printMajorFluxes(model_PP,metName,FLUX,N);
        % Consumption
        [Totflux_C,SortRxnIDs_C,SortProdFlux_C] = printMajorFluxes_Consumption(model_PP,metName,FLUX,N);
               
        % Results
        Results_Prod{j,i} = {Totflux_P SortRxnIDs_P SortProdFlux_P};
        Results_Cons{j,i} = {Totflux_C SortRxnIDs_C SortProdFlux_C};
        
        display(['mu = ' num2str(FLUX(238)) 'h-1'])
    end
end
%% Stacked Bar Plot
n_TopFlux = 5; % Consider only the 5 reactions that carry more flux in the network

for i=1:length(metList)
    
    % Get Data For Plots
    [ylimit_P,Fluxes_P{i,1},RxnIDs_P{i,1}]    = StackedBarPlot(Results_Prod(i,:),n_TopFlux,model_PP); % Production
    [~,Fluxes_C{i,1},RxnIDs_C{i,1}]           = StackedBarPlot(Results_Cons(i,:),n_TopFlux,model_PP); % Consumption
    
    % Production Plot
    figure(i)
    subplot(1,2,1)
    bar(Fluxes_P{i,1},'stacked'), legend(RxnIDs_P{i,1})
    ylim([0 ylimit_P]), ylabel('mmol/g_D_C_Wh')
    xlim([0 length(times)+1])
    title (['Production of ' metList(i)])
    
    % Consumption Plot
    subplot(1,2,2)
    bar(Fluxes_C{i,1},'stacked'), legend(RxnIDs_C{i,1})
    ylim([0 ylimit_P])
    title (['Consumption of ' metList(i)])
    xlim([0 length(times)+1])
end

%% Write node information in excel
for i=1:length(metList)
    
        filename = strcat(metIDs{i},' Balance.xlsx');
        Data_P = Results_Prod(i,:);
        Data_C = Results_Cons(i,:);
        writePieChartData(filename,Data_P,times,model_PP,3)
        writePieChartData(filename,Data_C,times,model_PP,length(times)+3) 
        
        % Write graph information in the first page two sheets
        cd('C:\Users\Francisco\Dropbox\RPP_dFBA_fed-batch\Cofactor Pie Charts')
        % Production
        xlswrite(filename,RxnIDs_P{i,1},1,'A1')
        xlswrite(filename,Fluxes_P{i,1}',1,'B1')
        % Consumption
        xlswrite(filename,RxnIDs_C{i,1},2,'A1')
        xlswrite(filename,Fluxes_C{i,1}',2,'B1')
        cd ..
        
end

