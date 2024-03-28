% This code counts the number of TICs obtained after applying the
% thermodynamically feasible direction constraints on all the BiGG models
clear
p=dir('./BiggBlkdResults/');
p={p(3:end).name}';
p1=p(~contains(p,'Incomp'));
nTICsBef=[];nTICsAft=[];nTICrxnsBef=[];nTICrxnsAft=[];
for z=1:numel(p1)
    load(['./BiggBlkdResults/',p1{z}])
    load(['./Bigg models/',p1{z}]) % directory of bigg models
    if exist('TIC_rxns', 'var') 
        TIC_Rxns = TIC_rxns;  
    end
    model = getModModel(model);
    model.lb(ismember(Tcc,'Forward')|ismember(Tcc,'Blocked'))=0;
    model.ub(ismember(Tcc,'Reverse')|ismember(Tcc,'Blocked'))=0;
   [TICs2,Direction2,TIC_Rxns2] = ThermOptEnumMILP(model);

    nTICsBef(z) = numel(TICs); nTICsAft(z) = numel(TICs2);
    nTICrxnsBef(z)=numel(TIC_Rxns);nTICrxnsAft(z)=numel(TIC_Rxns2);
    
    clear TIC_Rxns TIC_rxns TICs TIC_Rxns2 TICs2
end
tbl = table('Size',[99,5],'VariableTypes',{'string','double','double','double','double'},...
    'VariableNames',{'Model','nTICrxnsAft','nTICrxnsBef','nTICsBef','nTICsAft'});
tbl.Model =p1;
tbl.nTICrxnsAft=nTICrxnsAft';
tbl.nTICrxnsBef=nTICrxnsBef';
tbl.nTICsBef=nTICsBef';
tbl.nTICsAft=nTICsAft';

function model =getModModel(model)
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
% normalising the bounds to max of 1000
temp1 = max(abs(model.lb));
temp2 = max(abs(model.ub));
model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*1000;
model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*1000;

end