%% TORCs

% WARNING: DO NOT ADD FULL PATH TO BAPHY! MATLAB'S GET IS OVERLOADED! WTF!
addpath c:\src\baphy\RemoteAnalysis
addpath c:\src\baphy\Utilities
addpath c:\src\baphy\Meska

%% Setup
cfg = [];

cfg.path = 'G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG_Data\';

cfg.baphymfile = 'Mangrove_2021_05_06_TORCS.m'; % for Baphy events
cfg.experiment = 'Mangrove_2021-05-06_12-54-34_TORCS'; % OpenEphys data folder

baphymfile = fullfile(cfg.path,cfg.experiment,cfg.baphymfile);
% [ststims,stimparam] = loadtorc(baphymfile);
[ststims,stimparam,waveparams,TorcNames] = loadtorc(baphymfile);


% options = [];
% options.lfp = 1;
% [strfest,snr] = strf_online(baphymfile,1,gca,options)


LoadMFile(baphymfile);

TorcObject = exptparams.TrialObject.ReferenceHandle;



%

%% compute STRF quick & dirty

utrl = unique(data.trialinfo);


onsetSample  = round(data.fsample*(stimparam.stonset./1000));
offsetSample = round(data.fsample*(stimparam.stdur./1000))+onsetSample;

svec = onsetSample:offsetSample;

figure
pos = layout.pos(1:end-2,:) + 1;
pos(:,2) = max(pos(:,2)) - pos(:,2) + 1;
pos = fix(pos);
mlp = [7 4];

t = tiledlayout(mlp(2),mlp(1));
t.TileSpacing = 'compact';
t.Padding = 'tight';

smooth = [100 250]; % for plotting


strfest = [];
ax = [];
snr = [];
for ch = 1:length(data.label)
    
    spdata = [];
    for i = 1:length(utrl)
        
        ind = data.trialinfo == utrl(i);
        
        idx = find(ind,5);% limit to first 5 trials for each stimulus
        
        r = data.trial(idx);
        r = cat(3,r{:}); % -> channel x samples x reps

        spdata(:,:,i) = squeeze(r(ch,:,:));
    end
    spdata = spdata(svec,:,:); % limit to response window
    
%     spdata = sqrt(spdata .^ 2); % rectify
    
    [strfest(:,:,ch),indstrfs,resp]=pastspec(spdata,ststims,stimparam.basep, ...
        0,1000*TorcObject.Duration, ...
        stimparam.hfreq,stimparam.nrips,stimparam.mf,0);
    
    
    snr(ch) = get_snr(spdata,ststims,stimparam.basep,stimparam.mf,waveparams,stimparam.a1rv);
    
    
    % plot
    pidx = sub2ind(mlp,pos(ch,1),pos(ch,2));
    nexttile(pidx);
    lfreq = stimparam.lfreq;

    stdata = interpft(interpft(strfest(:,:,ch),smooth(1),1),smooth(2),2);
   
    fvec = stimparam.lfreq.*2.^linspace(0,stimparam.octaves,size(stdata,2));
    tvec = linspace(0,stimparam.basep,size(stdata,1));
    
    
%     timeRestrict = tvec >= 0 & tvec <= 100; % ms
%     imagesc(tvec(timeRestrict),fvec,stdata(:,timeRestrict)); axis xy


    imagesc(tvec,fvec,stdata); axis tight
    set(gca,'ydir','normal','clim',[-.9, .9]*max(abs(stdata(:))),'yscale','log');
    title(gca,sprintf('%s - %.2f',data.label{ch},snr(ch)))
    drawnow
end
%     strftemp = real(ifft(conj(fft(stim)).*fft(mr(ch,:)')));

%     strftemp = strftemp/stimT; % Normalization
%     strftemp = strftemp*(2*StimParam.nrips/mean(mean(stim.^2))/...
%         StimParam.numrecs);%Normalization

%     strfest = strfest + strftemp./realrepcount;



%% Continuous sample strf

skipFirstTorc = 1;
rectify = 0;

fsample = data.fsample;

utrl = unique(data.trialinfo);

onsetSample  = round(fsample*(stimparam.stonset./1000));
offsetSample = round(fsample*(stimparam.stdur./1000))+onsetSample-1;
stepSamples  = round(fsample.*stimparam.basep/1000);

% svec = onsetSample:stepSamples:offsetSample;
svec = onsetSample:offsetSample;


splitSamples = 1:stepSamples:offsetSample-stepSamples;

nSamples = length(svec);

[nFreqs,nDur,nStims] = size(ststims);

stimLengthSamples = round(fsample.*stimparam.basep./1000);
nTorcsInTrain = nSamples ./ stimLengthSamples;
%%%%%%
saf = fsample*nDur/stimparam.basep;  % stimulus temporal sampling frequency
maxlag = saf;

STRF = repmat({zeros(nFreqs,stimLengthSamples)},size(data.label));
for stimID = 1:nStims
    fprintf('Processing stim %d of %d ...',stimID,nStims)
    
    ind = data.trialinfo == utrl(stimID);
    
    idx = find(ind,5); % limit to 5 trials each
    
    r = data.trial(idx);
    r = cat(3,r{:}); % channel x samples x reps
    
    r = r(:,svec,:);
    
    if rectify
        r = sqrt(r.^2); % rectify
    end
    
%     r = sqrt(mean(r.^2,3)); % RMS
%     r = mean(r,3); % average time-locked resposne
    
    [nChan,nSamp,nReps] = size(r);
    
    % prepare stimulus parameters for interpolation
    stim = ststims(:,:,stimID);
    
    stim = cat(2,stim,stim(:,end));
    
    [m,n] = size(stim);
    
    % interpolate stimulus parameters to match sampling rate of the
    % recorded data
    istim = interp1(1:n,stim',linspace(1,n,stimLengthSamples),'makima')';
    
    % replicate the stimulus parameters to the number of torcs in the
    % stimulus train
    istim = repmat(istim,1,nTorcsInTrain,nReps);

    [m,n] = size(istim);
    
    parfor ch = 1:nChan
        strfTmp = zeros(m,nSamp,nReps);
        for i = 1:m
            strfTmp(i,:,:) = real(ifft(conj(fft(squeeze(istim(i,:,:)))).*squeeze(fft(r(ch,:,:)))));
        end
        strfTmp = strfTmp ./ nSamp;
        strfTmp = strfTmp*(2*stimparam.nrips/mean(stim(:).^2)/nStims); % ????? IS THIS CORRECT???
        
%         strfStd{ch} = std(srfTmp,0,3);
        
        strfTmp = mean(strfTmp,3); % average of strfs
        
        
        
        % breakup torc train
        y = zeros(nFreqs,stimLengthSamples);
        for i = splitSamples(1+skipFirstTorc:end-1) % skip first torc response?
            y = y + strfTmp(:,i:i+stimLengthSamples-1);
        end
        
        STRF{ch} = STRF{ch} + y;
    end
    
    fprintf(' done\n')

end

%%
strfPlotAmpThreshold = 0;
% strfPlotAmpThreshold = 30;

figure
pos = layout.pos(1:32,:) + 1;
pos(:,2) = max(pos(:,2)) - pos(:,2) + 1;
mlp = max(pos);

t = tiledlayout(mlp(2),mlp(1));
t.TileSpacing = 'compact';
t.Padding = 'tight';

smooth = [100 250]; % for plotting

for ch = 1:length(STRF)
    
    % plot
    pidx = sub2ind(mlp,pos(ch,1),pos(ch,2));
    nexttile(pidx);
    lfreq = stimparam.lfreq;

    fvec = stimparam.lfreq.*2.^linspace(0,stimparam.octaves,size(STRF{ch},1));
    tvec = linspace(0,stimparam.basep-1,size(STRF{ch},2));
    
    
%     timeRestrict = tvec >= 0 & tvec <= 100; % ms
%     imagesc(tvec(timeRestrict),fvec,stdata(:,timeRestrict)); axis xy

    stdata = interpft(interpft(STRF{ch},smooth(1),1),smooth(2),2);
    
    stdata(abs(stdata) < strfPlotAmpThreshold) = nan;

    imagesc(tvec,fvec,stdata); axis tight
    set(gca,'ydir','normal','yscale','log');
    set(gca,'clim',[-1 1]*min([max([abs(stdata(:)); 1]) inf]));
    colormap(gca,'parula');
    title(gca,data.label{ch})
    drawnow
end








