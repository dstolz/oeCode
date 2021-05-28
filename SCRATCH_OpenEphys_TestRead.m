%% notes
% https://open-ephys.atlassian.net/wiki/spaces/OEW/pages/65667092/Open+Ephys+format
% from: https://github.com/open-ephys/analysis-tools 
addpath('c:\src\analysis-tools\')



%%
pth = 'G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\';


expt = 'Mangrove_2021-04-21_11-25-55_TONES';
mfile = 'Mangrove_2021_04_21_TONES.m'; % for Baphy events

rcrdn = 'Record Node 105';

% origFs = 2000;
origFs = 1000;
Fs = 500;

load('G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\ECoG_channel_map.mat');

%%
% efn = 'all_channels.events';
clear data info
n = round(origFs/Fs);
Fs = origFs/n;
for i = 1:32
    efn = sprintf('107_CH%d.continuous',ChannelData.Map(i));
    
    ffn = fullfile(pth,expt,rcrdn,efn);
    [d,t,info(i)] = load_open_ephys_data(ffn);
    
    d = downsample(d,n);
    
    data(:,i) = d;
    
    fprintf('\bdone\n')
end
timestamps = downsample(t,n);


%% Read OpenEphys event (TTL) data
% eventIDs seem to be TTL input channel # - 1
% All events are stored as on/off times synchronized with physiology data
%  event 0 : Lick IR TTL
%  event 1 : Shock TTL
%  event 3 : Trial TTL

evID = 3;


ffn = fullfile(pth,expt,rcrdn,'all_channels.events');

[eventData,eventTimestamps,eventInfo] = load_open_ephys_data(ffn);

ind = evID == eventData;

eventTimestamps = eventTimestamps(ind);
evOnTimestamp   = eventTimestamps(1:2:end);
evOffTimestamp  = eventTimestamps(2:2:end);

%********************************************************************
%FOR SOME REASON, OPENEPHYS IS DETECTING ONE MORE TTL EVENT (evID= 3) THAN
%BAPHY IS REPORTING IT SENDS.  DON'T KNOW WHO'S TO BLAME. MAY HAVE TO
%REMOVE FIRST OR LAST EVENT TIMESTAMP TO FIGURE THIS OUT
% 
% evOnTimestamp(1) = [];
evOnTimestamp(end) = [];
evOffTimestamp(end) = [];
%********************************************************************

evOnSamples = round(evOnTimestamp .* Fs);
evOffSamples = round(evOffTimestamp.*Fs);

% Read Baphy stimulus event values
% NOTE: VALUES FROM BAPHY ARE STORED IN AN M FILE
ffn = fullfile(pth,expt,mfile);
run(ffn);

clear globalparams

ind = startsWith({exptevents.Note},'Stim');
exptevents(~ind) = [];


trialNum = [exptevents.Trial];

c = cellfun(@(a) textscan(a,'%*s%s%*[^\n]','delimiter',','),{exptevents.Note},'uni',1);
stimFreqs = cell2mat(cellfun(@(a) str2double(a{1}(1:find(a{1}==' ')-1)),c,'uni',0))';
stimAtten = cell2mat(cellfun(@(a) str2double(a{1}(find(a{1}==' ')+1:end)),c,'uni',0))';
uStimFreqs = unique(stimFreqs);
uStimAtten = unique(stimAtten);


% since baphy can't present a set number of repetitions, make sure to
% remove extras
c = arrayfun(@(a) sum(a==stimFreqs),uStimFreqs,'uni',1);
mc = min(c);
idx = find(c > mc);
eidx = [];
for i = 1:length(idx)
    x = find(stimFreqs == uStimFreqs(idx(i)),c(idx(i))-mc,'last');
    eidx(end+1:end+length(x)) = x;
end
stimFreqs(eidx) = [];
stimAtten(eidx) = [];
evOnTimestamp(eidx) = [];
evOffTimestamp(eidx) = [];
evOnSamples(eidx)   = [];
evOffSamples(eidx)  = [];

fprintf('\b done\n')





%% Port to FieldTrip

addpath('c:\src\fieldtrip\')
ft_defaults;


cfg = [];
cfg.datasetpath = pth;
cfg.dataset = expt;
cfg.trialdef.eventtype = 'stimulus';
cfg.trialdef.eventvalue = 3;
cfg.trialdef.prestim  = 0.2;
cfg.trialdef.poststim = 0.4;
cfg.trialdef.fsample = Fs;
cfg.trialfun = 'ft_oe_trialfun';

trl = [evOnSamples evOffSamples zeros(size(evOnSamples)) stimFreqs];


cfg = ft_definetrial(cfg);
% cfg = ft_redefinetrial(cfg,data');


%% Remove artifacts with ICA
[wIC,A] = wICA(data','fastica',1,0,Fs);
artifacts = A*wIC;
data = data - artifacts';


%% setup bandpass filter
bpFilt = designfilt('bandpassiir','FilterOrder',30, ...
         'HalfPowerFrequency1',1,'HalfPowerFrequency2',150, ...
         'SampleRate',Fs);
% fvtool(bpFilt)



% filtfilt continuous data
fprintf('Filtering ')

[nSamples,nChannels] = size(data);
data = [zeros(size(data)); data; zeros(size(data))];

for i = 1:size(data,2)
    data(:,i) = filtfilt(bpFilt,data(:,i));
    fprintf('.')
end

data(1:nSamples,:) = [];
data(nSamples+1:end,:) = [];
fprintf(' done\n')


%% Quick plot of all data
figure
fprintf('Plotting data ...')
cla
drawnow expose

xOffset = 0.2;
fprintf('Note: Plotting events with xOffset = %g s\n',xOffset)
ax = gca;

nChannels = size(data,2);

hold(ax,'on');
mv = median(abs(data(:)));
ysteps = linspace(0,(nChannels-1)*mv,nChannels);
plot(ax,xOffset+timestamps(:)*ones(1,nChannels),data+ysteps);
axis(ax,'tight');
y = ylim(ax);


txt = cellfun(@(a,b) sprintf('%d Hz\n%d dB',a,b),num2cell(stimFreqs),num2cell(stimAtten),'uni',0)';
ht = text(ax,xOffset+evOnTimestamp,max(y).*ones(size(evOnTimestamp(:,1),1),1),txt);
set(ht,'clipping','on');

h = plot((xOffset+evOnTimestamp(:,1)*[1 1])',y'*ones(1,length(evOnTimestamp))*1.1,'-','color',[.6 .6 .6]);
uistack(h,'bottom');
hold(ax,'off');
axis(ax,'tight')
set(ax,'ytick',[])

axis(ax,'tight');


xlabel(ax,'time (s)')
drawnow expose
fprintf(' done\n')








%% Cut raw stream by event onsets -> evData: time x channel x rep x freq x attenuation

win = [0.0 0.5];
[evData,evAB,evDataDims] = oe_cut_on_events(data,[evOnTimestamp stimFreqs stimAtten],win,Fs);


% %% Denoised Source Separation (DSS) for ECoG Data
% 
% addpath('M:\fUS\Denoising\NoiseTools')
% 
% L = 28;
% [todss,pwr0,pwr1] = nt_dss1(evData,[],L);
% 
% NKEEP = 28;
% fromdss = pinv(todss); 
% evData = oe_nt_mmat(evData ,todss(:,NKEEP)*fromdss(NKEEP,:));

%% baseline
bwin = [0.0 0.2];
evBData = oe_cut_on_events(data,[evOnTimestamp stimFreqs stimAtten],bwin,Fs);

%%
% rEvData = mean(evData,evDataDims.Repetition);
rEvData = median(evData,evDataDims.Repetition);
rEvBData = median(evBData,evDataDims.Repetition);

% rEvData = max(abs(rEvData),[],evDataDims.Time); % largest peak
rEvData = sqrt(mean(rEvData.^2)); % RMS
rEvBData = sqrt(mean(rEvBData.^2)); % RMS baseline

% re baseline
% rEvData = (rEvData - rEvBData)./rEvBData;
% rEvData = rEvData - rEvBData;

rEvData = squeeze(rEvData);

rEvData = shiftdim(rEvData,1);

m = [];
for i = 1:ChannelData.N
    x = ChannelData.Position(i,1);
    y = ChannelData.Position(i,2);
    m(:,:,y,x) = rEvData(:,:,i);
end
m = reshape(m,size(m,1),size(m,2),[]);

figure;
for i = 1:ChannelData.N
    subplot(4,8,i)
    idx = sub2ind([4 8],ChannelData.Position(i,2),ChannelData.Position(i,1));
    imagesc(evAB{1},evAB{2},rEvData(:,:,idx)');
end
ax = findobj(gcf,'type','axes');
set(ax,'ydir','normal')



% q = quantile(rEvData(:),[0.05 0.95]);
% % montage(m,'Size',[4 8],'DisplayRange',q,'Interpolation','nearest');
% t = imtile(m,'GridSize',[4 8]);
% imagesc(t);
% colormap parula
% set(gca,'ydir','normal')


%%


% rEvData = squeeze(sqrt(mean(evData.^2,1)));

rEvData = squeeze(max(abs(evData),[],1));

% rM = shiftdim(rEvData,2);
% montage(rM,'DisplayRange',[0 200],'Size',[4 8],'BorderSize',[5 5],'Interpolation','bilinear');

figure
colormap(jet)

for i = 1:nChannels
    subplot(4,8,i)
    y = squeeze(rEvData(i,:,:));
    imagesc(uStimFreqs,uStimAtten,y');
%     surf(uStimFreqs,uStimAtten,y');
%     view(2)
%     axis tight
%     shading flat
end

ax = findobj(gcf,'type','axes');
set(ax,'ydir','normal');%,'xscale','log')
set(ax,'clim',[0 60])

%%




%%
d = dir(fullfile(pth,expt,rcrdn,'*.continuous'));

clear data info
for i = 1:length(d)
    n = d(i).name;
    ffn = fullfile(d(i).folder,n);
%     ch = str2double(n(find(n=='H',1,'first')+1:find(n=='.',1,'last')-1));
    ch = str2double(n(find(n=='_',1,'first')+1:find(n=='.',1,'last')-1));
    fprintf('Loading channel %02d, "%s" ...',ch,n)
    
    
    % note that you could also use 'load_open_ephys_data' and manually
    % index into any position because it will memmap the large files.
    [data(:,ch),timestamps,info(ch)] = load_open_ephys_data_faster(ffn);
    
    if i == 1 % preallocate now that we know how many samples
        data(:,length(d)) = 0;
    end
    
    fprintf(' done\n')
end


%% event-based data reading

d = dir(fullfile(pth,expt,rcrdn,'*.continuous'));

clear data info

% eventIDs seem to be TTL input channel # - 1
% All events are stored as on/off times synchronized with physiology data
%  event 0 : Lick IR
%  event 1 : Shock TTL
%  event 3 : In Trial

eventID = 3;
ind = eventData == eventID;

eidx = find(ind);

ets = eventTimestamps(ind);

clear event
for e = 1:2:length(eidx)
    k = 1;
    for i = 1:length(d)
        n = d(i).name;
        ffn = fullfile(d(i).folder,n);
        ch = str2double(n(find(n=='_',1,'first')+1:find(n=='.',1,'last')-1));
        
%         eventSamples = 
        
        [event(k).data(:,ch),event(k).timestamps,event(k).info(ch)] = load_open_ephys_data(ffn,'Indices',s);
    end
    k = k + 1;
end










