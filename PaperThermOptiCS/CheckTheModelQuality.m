% This code check number of thermodynamically blocked reactions in the
% context-specific models built from fastcore and thermoptics
clear
nModels=64;
FnBlkd=[];FnTICs=[];FnTICRxns=[];FSize=[];
TnBlkd=[];TnTICs=[];TnTICRxns=[];TSize=[];
nComRxns=[];
for i=1:nModels
    f = load(['./FastCore_models/m',num2str(i)]);
    f = f.m;
    FSize(i) = numel(f.rxns);
    [a,TICs,~,~,TIC_Rxns] = ThermOptCC(f,1e-4);
    FnBlkd(i)=sum(ismember(a,'Blocked'));
    FnTICs(i)=numel(TICs);
    FnTICRxns(i)=numel(TIC_Rxns);
    
    if exist(['./ThermOptiCS_models/m',num2str(i),'.mat'],'file')
        t = load(['./ThermOptiCS_models/m',num2str(i)]);
        t = t.m;
        TSize(i) = numel(t.rxns);
        [a,TICs,~,~,TIC_Rxns] = ThermOptCC(t,1e-4);
        TnBlkd(i)=sum(ismember(a,'Blocked'));
        TnTICs(i)=numel(TICs);
        TnTICRxns(i)=numel(TIC_Rxns);
    end
    nComRxns(i) = numel(intersect(f.rxns,t.rxns));
end

tbl = table('Size',[64,9],'VariableTypes',{'string','double','double','double',...
    'double','double','double','double','double'},'VariableNames',{'Model',...
    'FnBlkd','FnTICRxns','FnTICs','FSize','TnBlkd','TnTICRxns','TnTICs','TSize'});
tbl.Model = arrayfun(@(x)['E',num2str(x)],[2199:2262]','UniformOutput',false);
tbl.FnBlkd=FnBlkd';tbl.FnTICRxns=FnTICRxns';tbl.FnTICs=FnTICs';tbl.FSize=FSize';
tbl.TnBlkd=TnBlkd';tbl.TnTICRxns=TnTICRxns';tbl.TnTICs=TnTICs';tbl.TSize=TSize';
tbl.nComRxns = nComRxns';

writetable(tbl,'FC_vs_TC.csv')