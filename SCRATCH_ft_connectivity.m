%% limit trials

%% freeup some memory if you need to
clear rawdata tldata

data.trial = cellfun(@single,data.trial,'uni',0);

%% optionally select subset of trials
ind = [];
% ind = indHit;
% ind = indMiss;
% ind = indFA;
% ind = indCR;


%% perform time-frequency analysis using Slepians
cfg            = [];
if ~isempty(ind)
    cfg.trials = find(ind);
end
cfg.output     = 'powandcsd';
cfg.method     = 'mtmconvol';
cfg.taper      = 'dpss';
% cfg.foi        = 2:50;
cfg.foi        = 2:2:150;

dt = diff(data.time{1}([1 end]));
cfg.t_ftimwin = min(dt,7./cfg.foi);
cfg.tapsmofrq = min(30,max(2,cfg.foi*0.2));
cfg.toi        = -0.5:0.01:0.6;
cfg.pad        = 'maxperlen';%'nextpow2';
% cfg.keeptrials = 'no';
cfg.keeptrials = 'yes';

tfr = ft_freqanalysis(cfg, data);

% save memory where possible
tfr.powspctrm = single(tfr.powspctrm);
tfr.crsspctrm = single(tfr.crsspctrm);

%% plot tfr
cfg         = [];
cfg.layout  = layout;
cfg.xlim    = [-0.2 0.6];
cfg.zlim    = [-2.0 2.0];

cfg.baseline = [-0.5 -0.1];
cfg.baselinetype = 'db';
cfg.hotkeys = 'yes';
cfg.colormap = jet;

cfg.showlabels = 'yes';
cfg.showcomment = 'no';
cfg.showscale = 'no';


ft_multiplotTFR(cfg,tfr);
drawnow

%% Computation and inspection of the connectivity measures


cfg = [];
% cfg.method = 'amplcorr'; parameter = 'amplcorrspctrm';
cfg.method = 'coh'; parameter = 'cohspctrm';
% cfg.method = 'granger';
% cfg.method = 'plv';
% cfg.method = 'pdc';
% cfg.method = 'ppc'; parameter = 'ppcspctrm';

coh = ft_connectivityanalysis(cfg, tfr);



%% Manual baseline correction for connectivity analysis

baseline = [-0.5 -0.1];

indb = coh.time >= baseline(1) & coh.time <= baseline(2);
mbl = mean(coh.(parameter)(:,:,indb),3);

coh.(parameter) = coh.(parameter)-mbl;
% coh.(parameter) = (coh.(parameter)-mbl)./mbl;

%% plot
cfg           = [];
cfg.parameter = parameter;

cfg.xlim = [0 0.8];
cfg.ylim = [0 60];

figure;
cfg.zlim = [-1 1]*quantile(abs(coh.(parameter)(:)),0.99);
ft_connectivityplot(cfg, coh);

colormap(jet);
%%
cfg         = [];
cfg.layout  = layout;
cfg.xlim    = [-0.2 0.8];
% cfg.zlim    = [-2.25 2.25];

cfg.parameter = parameter;

cfg.baseline = [-0.5 -0.1];
cfg.baselinetype = 'db';
cfg.hotkeys = 'yes';
cfg.colormap = jet;
cfg.showlabels = 'yes';

cfg.refchannel = 'CH20';

ft_multiplotTFR(cfg,coh);
%% compute the difference between mean coherence spectra for Hit and Miss trials

coh.cohspctrm = diff(cat(4,cohH.cohspctrm,cohM.cohspctrm),1,4);















