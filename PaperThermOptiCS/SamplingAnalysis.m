clear
nSamples = 10000;
load('./ThermOptiCS_models/m4.mat')
[TICs,Direction,TIC_Rxns,m,opt] = ThermOptEnumMILP(m);
chrr_samples= chrrSampler(m,100,nSamples);

chrr_samples_no_tic=zeros(size(chrr_samples));
for i=1:nSamples
    f = chrr_samples(:,i);
    [flux] = ThermOptFlux(m,f,TICs,Direction);
    chrr_samples_no_tic(:,i)=flux;
    clear flux
end

% visualizing the results
for i =1:numel(TIC_Rxns)
    id = find(ismember(m.rxns,TIC_Rxns(i)));
    figure()
    t = tiledlayout(1, 2, "TileSpacing", "tight");
    sgtitle(TIC_Rxns(i),'fontweight','bold','fontsize',25)

    nexttile
    histogram(chrr_samples(id,:))
    xlabel('Flux','fontweight','bold','fontsize',20)
    ylabel('Count','fontweight','bold','fontsize',20)
    title('Before implying ThermOptFlux','fontweight','bold','fontsize',20)
    set(get(gca, 'XAxis'), 'FontWeight', 'bold');
    set(get(gca, 'YAxis'), 'FontWeight', 'bold');
    nexttile
    histogram(chrr_samples_no_tic(id,:))
    xlabel('Flux','fontweight','bold','fontsize',20)
    ylabel('Count','fontweight','bold','fontsize',20)
    title('After implying ThermOptFlux','fontweight','bold','fontsize',20)
    set(get(gca, 'XAxis'), 'FontWeight', 'bold');
    set(get(gca, 'YAxis'), 'FontWeight', 'bold');
    set(gcf, 'Units', 'inches', 'Position', [0 0 15 5]);
    exportgraphics(gcf, ['./sampling_figs/samp_',TIC_Rxns{i},'.png'], 'Resolution', 300);
end