function writePieChartData(filename,metArray,times,model,Start_sheet)


for i=1:length(times)
    TotFlux = metArray{1,i}{1,1};
    
    SortRxnIDs = metArray{1,i}{1,2};
    idx = find(SortRxnIDs);
    SortRxnIDs = SortRxnIDs(idx);
    
    SortProdFlux = metArray{1,i}{1,3};
    SortProdFlux = SortProdFlux(idx);
    SortRxnNames = model.rxnNames(SortRxnIDs);
    SortFormulas = printRxnFormula(model,model.rxns(SortRxnIDs));
    
    cd('C:\Users\Francisco\Dropbox\RPP_dFBA_fed-batch\Cofactor Pie Charts')
    xlswrite(filename,TotFlux,Start_sheet,'A1')
    xlswrite(filename,SortRxnNames,Start_sheet,'B1')
    xlswrite(filename,SortProdFlux,Start_sheet,'C1')
    xlswrite(filename,SortFormulas,Start_sheet,'D1')
    Start_sheet = Start_sheet+1;
    cd ..
end

end