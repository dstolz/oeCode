%% BEHAVIOR/ECOG ANALYSIS
addpath c:\src\fieldtrip\  
addpath c:\src\analysis-tools\
addpath c:\src\oeCode\

addpath c:\src\DS_RefTar\ % Need object code for +bhvr

%% OE/Baphy to Fieldtrip

cfg = [];



cfg.path = 'G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\';




cfg.experiment = 'Mangrove_2021-05-14_10-26-09_Bhvr';
% cfg.experiment = 'Mangrove_2021-05-13_10-20-40_Bhvr';
% cfg.experiment = 'Mangrove_2021-05-12_11-39-03_Bhvr';
% cfg.experiment = 'Mangrove_2021-05-05_12-18-31_Bhvr';

cfg.trialdef.prestim = 1;
cfg.trialdef.poststim = 3;

cfg.channelmapfile = 'G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\ECoG_channel_map.mat';

cfg.trialfun = 'ft_oe_trialfun_bhvr';

tcfg = ft_definetrial(cfg);
[rawdata,tcfg] = ft_read_oe_data(tcfg);


load('G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\ECoG_channel_map_layout.mat')





%% Preprocess rawdata -> data
cfg = [];
cfg.lpfilter = 'yes';
% cfg.lpfreq = 60;
% cfg.lpfreq = 100;
cfg.lpfreq = 150;

cfg.hpfilter = 'yes';
cfg.hpfreq = 2;
% cfg.hpfreq = 30;
% cfg.hpfreq = 70;

% cfg.demean = 'yes';

cfg.reref = 'no';
% cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.refmethod = 'avg';

% cfg.hilbert = 'abs';

cfg.precision = 'single';

data = ft_preprocessing(cfg,rawdata);


% %%
% cfg = [];
% cfg.layout = layout;
% ft_databrowser(cfg,data)

%% Raw -> Trial-based  -  note that sampleinfo is updated automatically if resampled

data = ft_redefinetrial(tcfg,data);
data.trial = cellfun(@single,data.trial,'uni',0);
fprintf('\nYou can find trial data in data.cfg.event\n\n')




%% Identify trials

trialData = [data.cfg.event.trialData];

indTarg = [trialData.TrialType] == bhvr.TrialType.Target;
indRef  = [trialData.TrialType] == bhvr.TrialType.Reference;
indHit  = [trialData.Response]  == bhvr.Response.HIT;
indMiss = [trialData.Response]  == bhvr.Response.MISS;
indFA   = [trialData.Response]  == bhvr.Response.FALSE_ALARM;
indCR   = [trialData.Response]  == bhvr.Response.CORRECT_REJECT;


fprintf('#References = %d:\t#Correct Rejects = %d\t#False Alarms = %d\n',sum(indRef),sum(indCR),sum(indFA))
fprintf('#Targets = %d:\t\t#Hits = %d\t\t\t\t#Misses = %d\n',sum(indTarg),sum(indHit),sum(indMiss))




%% Timelock Average over all trials with baseline correction
cfg = [];
cfg.latency = [-.2 .8];
tldata = ft_timelockanalysis(cfg,data);

bcfg = [];
bcfg.baseline = 'yes';
tldata = ft_timelockbaseline(bcfg,tldata);


cfg = [];
cfg.layout = layout;
cfg.showscale = 'no';
cfg.showlabels = 'yes';
ft_multiplotER(cfg,tldata);






%%
ind = indTarg;
% ind = indRef;
% ind = indHit;
% ind = indMiss;
% ind = indCR;
% ind = indFA;



cfg = [];
% cfg.latency = [-.2 .8];
cfg.latency = [0 .8];
cfg.removemean = 'no';
cfg.trials = find(ind);
cfg.keeptrials = 'yes';
tldata = ft_timelockanalysis(cfg,data);


% bcfg = [];
% bcfg.baseline = 'yes';
% tldata = ft_timelockbaseline(bcfg,tldata);

%%
cfg = [];
cfg.layout = layout;
cfg.showscale = 'no';
cfg.showlabels = 'yes';
ft_multiplotER(cfg,tldata);


%% compute neurometic d' on time-locked data

cfg = [];
cfg.latency = [-0.4 0.8];
cfg.removemean = 'no';
cfg.keeptrials = 'yes';


tldata = ft_timelockanalysis(cfg,data);

bcfg = [];
bcfg.baseline = 'yes';
tldata = ft_timelockbaseline(bcfg,tldata);

% compute mean and std for Target and Reference trials seperately
mTarg = squeeze(mean(tldata.trial(indTarg,:,:),1));
mRef  = squeeze(mean(tldata.trial(indRef,:,:),1));
sTarg = squeeze(std(tldata.trial(indTarg,:,:),0,1));
sRef  = squeeze(std(tldata.trial(indRef,:,:),0,1));

% rectify
mTarg = sqrt(mTarg.^2);
mRef  = sqrt(mRef.^2);
sTarg = sqrt(sTarg.^2);
sRef  = sqrt(sRef.^2);

dprime = 2.*(mTarg - mRef) ./ (sTarg + sRef);

% dprime = dprime - mean(dprime(:));

tldata_dPrime = rmfield(tldata,{'sampleinfo','trial'});
tldata_dPrime.avg = dprime;
% tldata_dPrime.var = 
% tldata_dPrime.dof = 
tldata_dPrime.dimord = 'chan_time';
tldata_dPrime.cfg.keeptrials = 'no';


cfg = [];
cfg.layout = layout;
cfg.showscale = 'no';
cfg.showlabels = 'yes';
% cfg.xlim = [0 0.8];
% cfg.ylim = [0 5];
cfg.linewidth = 1;
cfg.linecolor = [1 0.2 0.2];

ft_multiplotER(cfg,tldata_dPrime);




%% Compute time-frequency response by condition

ind = indHit;
% ind = indMiss;
% ind = indCR;

cfg = [];
cfg.method = 'superlet';
cfg.superlet.basewidth = 7;
cfg.superlet.gwidth = 3;
cfg.superlet.combine = 'multiplicative';%'additive';
% cfg.superlet.order 

cfg.keeptrials = 'no';

cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.toi = 0:0.01:0.8;
cfg.foi = 20:2:150;
% cfg.foi = logspace(log10(2),log10(150),50);


cfg.trials = find(ind);
freq = ft_freqanalysis(cfg,data);

freq.powspctrm = 10.*log10(freq.powspctrm);

%
cfg = [];
cfg.layout = layout;

cfg.xlim = [0 0.8];
% cfg.zlim = [0 32]
cfg.colormap = jet;
ft_multiplotTFR(cfg,freq);


%% TFR d'

cfg = [];
cfg.method = 'superlet';
cfg.superlet.basewidth = 7;
cfg.superlet.gwidth = 3;
cfg.superlet.combine = 'multiplicative';%'additive';
% cfg.superlet.order 

cfg.keeptrials = 'yes';

cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.toi = 0:0.01:0.8;
cfg.foi = 20:2:150;
% cfg.foi = logspace(log10(2),log10(150),75);


freq = ft_freqanalysis(cfg,data);

%%

dfreq = freq;

dfreq.powspctrm = 10.*log10(dfreq.powspctrm);

mTarg = squeeze(mean(freq.powspctrm(indTarg,:,:,:),1));
mRef  = squeeze(mean(freq.powspctrm(indRef,:,:,:),1));
sTarg = squeeze(std(freq.powspctrm(indTarg,:,:,:),0,1));
sRef  = squeeze(std(freq.powspctrm(indRef,:,:,:),0,1));

dprime = 2.*(mTarg - mRef) ./ (sTarg + sRef);

freq.powspctrm = dprime;
freq.cfg.keeptrials = 'no';

%%
cfg = [];
cfg.layout = layout;
ft_multiplotTFR(cfg,freq);










