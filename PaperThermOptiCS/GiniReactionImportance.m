function [RxnImp,Contexts] = GiniReactionImportance(geneExpression,model,ut,lt,coreRxn)

%%INPUT
%       geneExpression: matlab structure with fields
%                       .value : mRNA expression matrix with dimension N_genes*N_samples
%                       .genes : cell array with geneIDs in the same format as model.genes
%                       .context : cell array with names of the samples
%
%       model: COBRA model structure

%%OPTIONAL INPUTS        
%       ut: Scalar value [0 100] (upper threshold percentile) (default:100)
%         
%       lt: Scalar value [0 100] (lower threshold percentile) (default:0)
%         
%       coreRxn: Integer vector denoting indices of reactions for which 
%                higher importance has to be given manually (default:none)

%%OUTPUT
%       RxnImp: Reaction importance matrix with dimension N_genes*N_samples

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and Control (BiSECt) lab, IIT Madras


% Setting the default values
if ~exist('ut','var') || isempty(ut)
    ut =100;
end

if ~exist('lt','var') || isempty(lt)
    lt =0;
end


if ~exist('coreRxn','var') || isempty(coreRxn)
    coreRxn =[];
end

Contexts = geneExpression.context;

% getting expression values of genes available in the model
[idxa,~] = ismember(geneExpression.genes,model.genes);
model_exp_value = geneExpression.value(idxa,:);
model_exp_gene = geneExpression.genes(idxa);
parsedGPR = GPRparser(model);


temp_var = model_exp_value;
temp_var(temp_var==0)=[];    
LT = prctile(temp_var,lt); % lower threshold value
UT = prctile(temp_var,ut); % upper threshold value   
clear temp_var
Loc_Gini = ginicoeff(model_exp_value); % gini coefficient based thresholding
Loc_Gini(Loc_Gini>=UT)=UT;Loc_Gini(Loc_Gini<=LT)=LT;
gene_exp = log2((model_exp_value./repmat(Loc_Gini,1,numel(Contexts)))+1); 
    
    
for i=1:numel(Contexts)
    rxn_exp(:,i) = selectGeneFromGPR(model, model_exp_gene, gene_exp(:,i), parsedGPR); % gene to reaction mapping
end
    
if ~isfield(model,'rxnGeneMat')
    model=buildRxnGeneMat(model);
end

RxnImp = rxn_exp;
RxnImp(find(sum(model.rxnGeneMat,2)==0),:)=0.001;    
RxnImp(coreRxn,:)=1;
end

function [gcp]=ginicoeff(data)
    data_sort=sort(data,2);
    n=size(data,2);
    
    for i=1:size(data,1)
        x=data_sort(i,:);
        G_num=sum(((2*[1:n])-n-1).*x);
        G_den=sum(x)*n;
        gc(i)=(G_num/G_den)*100;
    end
    
    for i=1:numel(gc)
        gcp(i)=prctile(data(i,:),gc(i));
    end
    gcp=gcp';
end 

