clear
% initCobraToolbox(0);
p = dir('../../Microbial_community_gap_filling/With_AGORA2/AGORA2/');
agoraModels = {p(3:end).name}';
p2 = dir('./AGORA_TIC_Results/');
p2 = {p2(3:end).name}';
agoraModels = setdiff(agoraModels,p2);
agoraModels = agoraModels(randperm(length(agoraModels)));
for i =1:numel(agoraModels)
    load(['../../Microbial_community_gap_filling/With_AGORA2/AGORA2/',agoraModels{i}])
    
    [TICs,Direction,TIC_Rxns,modModel,opt] = ThermoOptEnumMILP(model,300);
    if opt
        save(['./AGORA_TIC_Results/',agoraModels{i}],'TICs','Direction','TIC_Rxns')
    end
    clear model TICs Direction TIC_Rxns modModel opt
end