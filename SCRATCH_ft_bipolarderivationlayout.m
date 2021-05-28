%% bipolar derivation
 
distThreshold = 1.5;
databip = rawdata;
pos = layout.pos(1:end-2,:);
n = size(pos,1);
for i = 1:n
    d = vecnorm(pos(i,:) - pos,2,2); % distance to all other electrodes
    id = 1 ./ d;
    ind = d > 0 & d <= distThreshold;
    databip.trial{1}(i,:) = rawdata.trial{1}(i,:) - mean(id(ind).*rawdata.trial{1}(ind,:),1);
end



%%

layoutbip = layout;


layoutbip.label(1:32) = [];
layoutbip.label = [data.label; layoutbip.label];

layoutbip.width(32) = [];
layoutbip.height(32) = [];
layoutbip.pos(1:32,:) = [];

pos = layout.pos;
posbip = [];
for i = 1:31
    posbip(i,:) = mean(pos(i:i+1,:));
end

layoutbip.pos = [posbip; layoutbip.pos];

layout = layoutbip;

%%
cfg = [];
cfg.layout = layoutbip;
ft_layoutplot(cfg)
%%

save('G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\ECoG_channel_map_layoutbip.mat','layout')