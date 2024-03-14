function [order,bins,binsizes] =getOrderOfRxns(model)
S = model.S~=0; S=logical(S'*S);
g = graph(S);
[bins,binsizes]=conncomp(g);
[binsizes,i]=sort(binsizes);
bins = arrayfun(@(x)find(i==x),bins);
order=[];
for i =1:numel(binsizes)
    order=[order;find(bins==i)'];
end
end