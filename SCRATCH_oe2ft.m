%%
addpath c:\src\fieldtrip\  
addpath c:\src\analysis-tools\
addpath c:\src\oeCode\
addpath M:\EEG\NoiseTools

%% OE/Baphy to Fieldtrip

cfg = [];



cfg.path = 'G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\';



% cfg.experiment = 'Mangrove_2021-05-12_12-12-24_ToneClouds';
% cfg.trialdef.skiptriggers = 1:2;
% cfg.trialdef.prestim = 0.1;
% cfg.trialdef.poststim = 3.5;
% cfg.baphystimformat = ????????????????

% cfg.experiment = 'Mangrove_2021-05-06_12-54-34_TORCS'; % OpenEphys data folder
% cfg.trialdef.skiptriggers = 1:2;
% cfg.trialdef.prestim = 0;
% cfg.trialdef.poststim = 3.5;
% cfg.baphystimformat = '%*s%s%*[^\n]';

% cfg.experiment = 'Mangrove_2021-05-19_12-19-06_TONES'; % OpenEphys data folder
% cfg.trialdef.skiptriggers = 1:2;
% cfg.trialdef.prestim = 0;
% cfg.trialdef.poststim = 0.25;
% cfg.baphystimformat = '%*s%f%f%*[^\n]';
% cfg.node = '100';

cfg.experiment = 'Mangrove_2021-05-20_11-49-21_Clicks'; % OpenEphys data folder
cfg.trialdef.skiptriggers = 1:2;
cfg.trialdef.prestim = 1;
cfg.trialdef.poststim = 1;
cfg.baphystimformat = '%*s%f%*[^\n]';
cfg.node = '100';

% cfg.experiment = 'Mangrove_2021-05-26_12-30-31_Clicks'; % OpenEphys data folder
% cfg.trialdef.skiptriggers = 1:2;
% cfg.trialdef.prestim = 1;
% cfg.trialdef.poststim = 1;
% cfg.baphystimformat = '%*s%f%*[^\n]';
% cfg.node = '100';


% cfg.experiment = 'Mangrove_2021-05-27_11-32-07_ClicksSE'; % OpenEphys data folder
% cfg.trialdef.skiptriggers = 1:2;
% cfg.trialdef.prestim = 1;
% cfg.trialdef.poststim = 1;
% cfg.baphystimformat = '%*s%f%*[^\n]';
% cfg.node = '100';



cfg.channelmapfile = 'G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\ECoG_channel_map.mat';

cfg.trialfun = 'ft_oe_trialfun';

tcfg = ft_definetrial(cfg);
[rawdata,tcfg] = ft_read_oe_data(tcfg);


load('G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\ECoG_channel_map_layout.mat')



%% standard bipolar derivation


fprintf('Applying bipolar derivation with x-axis neighbors ...')

databip = rmfield(rawdata,'trial');

pos = layout.pos(1:end-2,:);
ux = unique(pos(:,1));
uy = unique(pos(:,2));

posNew = nan(size(pos));
for i = 1:length(ux)-1 % derive bipolar electrodes across x direction
    for j = 1:length(uy)
        indCur = pos(:,1) == ux(i)   & pos(:,2) == uy(j);
        indNbr = pos(:,1) == ux(i+1) & pos(:,2) == uy(j);
        
        databip.trial{1}(indCur,:) = rawdata.trial{1}(indCur,:) - rawdata.trial{1}(indNbr,:);
        
        posNew(indCur,:) = [mean(ux(i:i+1)) uy(j)];
    end
end

ind = pos(:,1) == ux(end); % remove lost channels

databip.trial{1}(ind,:) = [];
databip.label(ind) = [];

posNew(ind,:) = [];
layout.pos(1:size(posNew,1),:) = posNew;
layout.label(ind) = [];
layout.width(ind,:) = [];
layout.height(ind,:) = [];
layout.cfg.channel(ind) = [];
layout.cfg.layout = structfun(@(a) a(~ind,:),layout.cfg.layout,'uni',0);

rawdata = databip;
clear databip
fprintf(' done\n')



%% weighted mean neighbour channel derivation - **** UNVERIFIED ****
% reject common sources by subtracting out mean activity of neighboring
% channels inversly weighted by distance.
distThreshold = 1.5;
fprintf('Applying bipolar derivation with distance threshold <= %.2f mm ...',distThreshold)
databip = rawdata;
pos = layout.pos(1:end-2,:);
n = size(pos,1);
for i = 1:n
    d = vecnorm(pos(i,:) - pos,2,2); % distance to all other electrodes
    id = 1 ./ d;
    ind = d > 0 & d <= distThreshold;
    databip.trial{1}(i,:) = rawdata.trial{1}(i,:) - mean(id(ind).*rawdata.trial{1}(ind,:),1);
end

rawdata = databip;
clear databip
fprintf(' done\n')


% %% rawdata -> CSD rawdata
% cfg = [];
% cfg.method = 'distance';
% cfg.neighbourdist = 1.5;
% cfg.layout = layout;
% neighbours = ft_prepare_neighbours(cfg);
% % cfg = [];
% % cfg.neighbours = neighbours;
% % cfg.layout = layout;
% % ft_neighbourplot(cfg);
% 
% cfg = [];
% % cfg.method = 'hjorth';
% cfg.method = 'spline';
% cfg.conductivity = 0.33; % S/m (INCORRECT! VALUE FOR ECOG? ASSUMES HOMOGENEITY)
% cfg.neighbours = neighbours;
% cfg.elec = 'c:\src\oeCode\NN_32chECoG_electrodes.elc';
% rawdata = ft_scalpcurrentdensity(cfg,rawdata);
 

% %% resample the data
% cfg            = [];
% cfg.resamplefs = 400;
% cfg.detrend    = 'yes';
% rawdata = ft_resampledata(cfg, rawdata);
% rawdata.sampleinfo = [1 length(rawdata.time{1})];
%% Preprocess rawdata -> data
cfg = [];
cfg.lpfilter = 'yes';
% cfg.lpfreq = 50;
% cfg.lpfreq = 100;
cfg.lpfreq = 150;

cfg.hpfilter = 'yes';
cfg.hpfreq = 2;
% cfg.hpfreq = 30;
% cfg.hpfreq = 70;

cfg.demean = 'yes';


% cfg.hilbert = 'abs';



data = ft_preprocessing(cfg,rawdata);

%% Raw -> Trial-based  -  note that sampleinfo is updated automatically if resampled
data = ft_redefinetrial(tcfg,data);





%% Standardize to z-score

% data.trial = zscore(data.trial,0,2);
% fprintf('Standardized trial data to z score manually\n')






%%
cfg = [];
cfg.layout = layout;
ft_databrowser(cfg,data)



%% Timelock Average over all trials with baseline correction
cfg = [];
% cfg.latency = [0 1];
% cfg.keeptrials = 'yes';
% cfg.removemean = 'yes';
tldata = ft_timelockanalysis(cfg,data);

% bcfg = [];
% bcfg.baseline = [-0.1 0];
% cfg.baselinetype = 'avg';
% tldata = ft_timelockbaseline(bcfg,tldata);
%%
cfg = [];
cfg.layout = layout;
cfg.showscale = 'no';
cfg.showlabels = 'yes';
cfg.xlim = [-0.5 0.8];
ft_multiplotER(cfg,tldata);
sgtitle(data.cfg.experiment,'interpreter','none')

%% NSL implementation of DSS

NKEEP = 5;

x = permute(tldata.trial,[3,2,1]); % -> time * channels * trials
% x = nt_demean2(x); %to remove the mean from channels

nt_dss1(x,[],NKEEP);
[todss,pwr0,pwr1] = nt_dss1(x);
fromdss = pinv(todss);
z = nt_mmat(x,todss);
xx = nt_mmat(z,todss(:,1:NKEEP)*fromdss(1:NKEEP,:));

tldata.trial = permute(xx,[3,2,1]); % -> trials * channels * time

cfg = [];
cfg.layout = layout;
cfg.showscale = 'no';
cfg.showlabels = 'yes';
ft_multiplotER(cfg,tldata);

%% Timelock Average for individual stimuli
cfg = [];
% cfg.latency = [0 3.5];
% cfg.latency = [-.25 .5];
cfg.latency = [0 .25];
cfg.removemean = 'no';

bcfg = [];
bcfg.baseline = [-.25 0];
cfg.baselinetype = 'absolute';

clear tldata
uti = unique(data.trialinfo);
for i = 1:length(uti)
    ind = data.trialinfo == uti(i);
    cfg.trials = find(ind);
%     tldata(i) = ft_timelockanalysis(cfg,data);
    tldata = ft_timelockanalysis(cfg,data);
    tldata(i) = ft_timelockbaseline(bcfg,tldata);
end

%%
trialIdx = 1;
trialIdx = length(tldata);

cfg = [];
cfg.layout = layout;
cfg.showscale = 'no';
cfg.showlabels = 'yes';
cfg.xlim = [0 3.5];
ft_multiplotER(cfg,tldata(trialIdx));


%% Stimulus selectivity

% win = [0 0.1];
win = [0.01 0.03];


wind = cellfun(@(a) a >= win(1) & a < win(2),data.time,'uni',0);

M = cell2mat(cellfun(@(a,b) rms(a(:,b),2),data.trial,wind,'uni',0));

zM = zscore(M,0,2);


uti = unique(data.trialinfo);
clear R
for i = 1:length(uti)
    ind = data.trialinfo == uti(i);
    R(:,i) = mean(zM(:,ind),2);
end

[mR,mRi] = max(R,[],2);


f = figure;
ax = axes(f);

pos = layout.pos(1:end-2,:) + 1;

m = zeros(max(pos(:,2)),max(pos(:,1)));
for i = 1:length(mRi)
    m(pos(i,2),pos(i,1)) = mRi(i);
%     m(pos(i,2),pos(i,1)) = uti(mRi(i));
end

h = heatmap(m);
colormap(flipud(jet));
caxis([1 length(uti)]);
% caxis(uti([1 end]));


%%

val = 500;

tpcfg = [];

tpcfg.trials = find(data.trialinfo == val)';


tpcfg.layout = layout;
tpcfg.interplimits = 'electrodes';
tpcfg.shading = 'interp';
tpcfg.colorbar = 'no';
tpcfg.colorbartext = 'mV';
tpcfg.interactive = 'yes';
tpcfg.comment = 'no';
% tpcfg.zlim = [-2 2];
cfg.latency = [0 0.1];
ft_topoplotER(tpcfg,data);
set(gcf,'windowstyle','docked');
title(val)



%%

f = figure;
tiledlayout(f,'flow')

v = cellfun(@(a) rms(a,2),data.trial,'uni',0);

m = cell(size(uti));
for i = 1:length(uti)
    ind = data.trialinfo == uti(i)';
    
    trv = cell2mat(v(ind));
    
    pos = layout.pos(1:end-2,:) + 1;
    
    for j = 1:size(pos,1)
        m{i}(pos(j,2),pos(j,1)) = mean(trv(j,:));
    end
    
    ax = nexttile;
    imagesc(m{i})
    
    ax.YDir = 'normal';
    ax.XAxis.TickValues = [];
    ax.YAxis.TickValues = [];
    ax.Title.String = sprintf('%d Hz',uti(i));
end

C = cat(3,m{:});
[mCv,mCi]= max(C,[],3)

%%
cfg = [];
dur = 0.05;
lat = 0:0.01:0.45;
clear tldata
for j = 1:length(lat)
    for i = 1:length(uti)
        ind = data.trialinfo == uti(i);
        cfg.trials = find(ind);
        cfg.latency = [0 dur]+lat(j);
        tldata(i,j) = ft_timelockanalysis(cfg,data);
    end
end

%%
for j = 1:length(lat)
    for i = 1:length(uti)
        tpcfg = [];
        tpcfg.layout = layout;
        tpcfg.interplimits = 'electrodes';
        tpcfg.shading = 'interp';
        tpcfg.colorbar = 'no';
        tpcfg.colorbartext = 'mV';
        tpcfg.interactive = 'no';
        tpcfg.comment = 'no';
        % tpcfg.zlim = [-2 2];
        cfg.latency = [0 0.1];
        ft_topoplotER(tpcfg,tldata(i));
        set(gcf,'windowstyle','docked');
        title(num2str(uti(i),'%d Hz'))
    end
    pause
    close all
end
%% create topoplot animation

dur = 0.05;
lat = 0:0.01:0.45;

tpcfg = [];
tpcfg.layout = layout;
tpcfg.interplimits = 'electrodes';
tpcfg.shading = 'interp';
tpcfg.colorbar = 'no';
tpcfg.colorbartext = 'mV';
tpcfg.interactive = 'no';
tpcfg.comment = 'no';
tpcfg.zlim = [-7 7];

clear M
for i = 1:length(lat)
    cfg = [];
    cfg.latency = [0 dur]+lat(i);
    fprintf('plotting latency: %s s\n',mat2str(cfg.latency))
    tldata = ft_timelockanalysis(cfg,data);
    
    ft_topoplotER(tpcfg,tldata);
    set(gcf,'position',[2100 90 1200 870],'renderer','painters');
    title(sprintf('Latency: %s s',mat2str(cfg.latency)));
    drawnow
    M(i) = getframe(gcf);
    close
end
%%
figure;
set(gcf,'position',[2100 90 1200 870]);
axes('Position',[0 0 1 1])
movie(M,1,5);








%% frequency analysis

cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.pad = 'nextpow2';

dt = diff(data.time{1}([1 end]));
% cfg.foi       = round(logspace(log10(2),log10(150),50));
cfg.foi       = 2:150;
cfg.t_ftimwin = min(dt,7./cfg.foi);
cfg.tapsmofrq = min(30,max(2,cfg.foi*0.5));


% cfg.toi = 'all';
cfg.toi = -.5:0.05:0.75;

tfdata = ft_freqanalysis(cfg,data);

%%

cfg = [];
cfg.baseline = [-0.5 0];
cfg.baselinetype = 'relative';
cfg.zlim = [0 3];

cfg.layout = layout;

ft_multiplotTFR(cfg,tfdata);











%%

M = squeeze(mean(tldata.trial,1));

for i = 1:size(M,1)
    t = squeeze(tldata.trial(:,i,:));
%     r = corrcoef([M(i,:)' t']);
%     rm = r(:,1);
%     ind = rm > 0.1;

    t = t - M(i,:);
end







%% ICA
cfg                 = [];
cfg.method          = 'fastica';
cfg.runica.maxsteps = 1000;
comp = ft_componentanalysis(cfg, data);

%%

cfg = [];
cfg.layout = layout;
cfg.viewmode = 'component';
ft_databrowser(cfg,comp)

%% project the original data thru the components again
cfg           = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp = ft_componentanalysis(cfg, data);

%%
cfg = [];
cfg.layout = lcfg.layout;
ft_databrowser(cfg,comp)
%%
cfg           = [];
cfg.component = [17 20];
data_clean     = ft_rejectcomponent(cfg, comp, data);

%%
cfg = [];
cfg.layout = lcfg.layout;
ft_databrowser(cfg,data_clean)



