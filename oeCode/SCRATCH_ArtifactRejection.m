%% NSL implementation of DSS
addpath M:\EEG\NoiseTools

NKEEP = 12;

x = permute(tldata.trial,[3,2,1]); % -> time * channels * trials


nt_dss1(x);

[todss,pwr0,pwr1] = nt_dss1(x,[],[],NKEEP);
% z = nt_mmat(x,todss);

fromdss = pinv(todss); 
xx = nt_mmat(x ,todss(:,NKEEP)*fromdss(NKEEP,:));


tldata.trial = permute(xx,[3,2,1]); % -> trials * channels * time


%%


%% vvvvvvvv FieldTrip implementation of DSS vvvvvvvvv
data = ft_redefinetrial(tcfg,rawdata);



%%
cfg = [];
cfg.artfctdef.zvalue.channel         = 'all';
cfg.artfctdef.zvalue.cutoff          = 11;
cfg.artfctdef.zvalue.interactive     = 'yes';
% cfg.artfctdef.zvalue.interactive     = 'no';
cfg.artfctdef.zvalue.bpfilter        = 'yes';
cfg.artfctdef.zvalue.bpfreq          = [2 30];
cfg.artfctdef.zvalue.hilbert         = 'yes';
cfg.artfctdef.zvalue.artfctpeak      = 'yes';
cfg.artfctdef.zvalue.artfctpeakrange = [-.1 .25];
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.continuous = 'yes';
[cfg,art] = ft_artifact_zvalue(cfg,data);



%%
params.tr  = cellfun(@abs,cfg.artfctdef.zvalue.peaks_indx,'uni',0);
% params.tr  = cfg.artfctdef.zvalue.peaks;
params.pre = 0.1*data.fsample;
params.pst = 0.25*data.fsample;
params.demean = true;

%% DSS component rejection
cfg                   = [];
cfg.method            = 'dss';
cfg.dss.denf.function = 'denoise_avg2';
cfg.dss.denf.params   = params;
cfg.dss.wdim          = 16;
cfg.numcomponent      = 6;
cfg.channel           = 'all';
cfg.cellmode          = 'yes';
comp                  = ft_componentanalysis(cfg, data);

%% Visualize artifactual components ranked from largest offender
cfg          = [];
cfg.layout   = layout;
ft_databrowser(cfg, comp)




%%
cfg           = [];
cfg.component = [1 2];
data_clean     = ft_rejectcomponent(cfg, comp, data);

%%
data_clean = ft_redefinetrial(tcfg,data_clean);

%%
cfg = [];
cfg.layout = layout;
ft_databrowser(cfg,data_clean)
