%% Tone Receptive Field Map

baseSoundLevel = 70; % dB
interpMultiplier = 3; 
maxNumberReps = 10; % ensure equal number of stimulus presentations because baphy doesn't automatically stop

% aTimeWindow = [0 0.05]; % seconds; restrict analysis to a time window
% aTimeWindow = [0.05 0.1];
% aTimeWindow = [0.15 3.2];
aTimeWindow = [.04 0.07];

metric = 'rms';
% metric = 'max';
% metric = 'min';
% metric = 'maxabs';


% recover the frequency and attenuation data from trialinfo
Freq = fix(data.trialinfo);
Attn = round((data.trialinfo - Freq) * 100);

uFreq = unique(Freq);
uAttn = unique(Attn);

% compute some metric over all trials
switch metric
    case 'rms'
        rTrial = cellfun(@(a,b) rms(a(:,b>=aTimeWindow(1)&b<aTimeWindow(2)),2),data.trial,data.time,'uni',0);
    case 'max'
        rTrial = cellfun(@(a,b) max(a(:,b>=aTimeWindow(1)&b<aTimeWindow(2)),[],2),data.trial,data.time,'uni',0);
    case 'min'
        rTrial = cellfun(@(a,b) min(a(:,b>=aTimeWindow(1)&b<aTimeWindow(2)),[],2),data.trial,data.time,'uni',0);
    case 'maxabs'
        rTrial = cellfun(@(a,b) max(abs(a(:,b>=aTimeWindow(1)&b<aTimeWindow(2))),[],2),data.trial,data.time,'uni',0);
end


rTrial = cat(2,rTrial{:}); % -> channel x metric

R = nan(length(uAttn),length(uFreq),size(rTrial,1),maxNumberReps);

for f = 1:length(uFreq)
    for a = 1:length(uAttn)
        ind = Freq == uFreq(f) & Attn == uAttn(a);
        
        idx = find(ind,maxNumberReps);
        
        R(a,f,:,:) = rTrial(:,idx); % -> Attn x Freq x Channel x Repetition
    end 
end

% compute mean over repetitions
mR = mean(R,4); % -> Attn x Freq x Channel
% mR = median(R,4); % -> Attn x Freq x Channel



uSL = baseSoundLevel - uAttn; % Attenuation -> Sound Level

nmstr = sprintf('%s - %s %ssec',data.cfg.experiment,metric,mat2str(aTimeWindow));
f = figure('Name',nmstr);
f.NumberTitle = 'off';
pos = layout.pos(1:32,:) + 1;
pos(:,2) = max(pos(:,2)) - pos(:,2) + 1;
mlp = max(pos);

t = tiledlayout(mlp(2),mlp(1));
t.TileSpacing = 'compact';
t.Padding = 'tight';


ax = [];
for ch = 1:size(R,3)
    imR = interp2(mR(:,:,ch),interpMultiplier,'makima');

    pidx = sub2ind(mlp,pos(ch,1),pos(ch,2));
    ax(ch) = nexttile(pidx);
    
    imagesc(uFreq,uSL,imR);
    
    title(data.label{ch});
    axis tight
end

sgtitle(nmstr,'interpreter','none','fontname','Consolas');
set(ax,'ydir','normal','xscale','log');
set(findobj(f,'-property','fontname'),'fontname','Consolas')




