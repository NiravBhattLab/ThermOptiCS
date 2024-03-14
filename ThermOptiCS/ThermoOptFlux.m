function [flux] = ThermoOptFlux(model,flux,TICs,Dir)


for i =1:numel(TICs)
    t = TICs{i}; c = Dir{i};
    ids = find(ismember(model.rxns,t));
    f=flux(ids);
    if sum(f>0&c>0)+sum(f<0&c<0)==numel(t)
        [~,mID] = min(abs(f));
        c = (c/abs(c(mID)))*abs(f(mID));
        flux(ids) = flux(ids)-c;
    end
end

end