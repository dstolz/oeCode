%% Tone Receptive Field Map

baseSoundLevel = 70; % dB
interpMultiplier = 3;
maxNumberReps = 10; % ensure equal number of stimulus presentations because baphy doesn't automatically stop

remMean = false;
% remMean = true;

window = [0 0.1]; %s
windowStepSize = 0.001; % s
windowDuration = 0.01; % s


metric = 'mean';
% metric = 'rms';
% metric = 'std';
% metric = 'max';
% metric = 'min';
% metric = 'maxabs';






% recover the frequency and attenuation data from trialinfo
Freq = fix(data.trialinfo);
Attn = round((data.trialinfo - Freq) * 100);

uFreq = unique(Freq);
uAttn = unique(Attn);

tvec = window(1):windowStepSize:window(2);


% Atten x Freq x Channel x Frame
mR = zeros(length(uAttn),length(uFreq),size(data.trial{1},1),length(tvec));

% first compute
for tm = 1:length(tvec)
    fprintf('Computing frame %d of %d ...',tm,length(tvec))

    thisTimeWindow = tvec(tm)+[0 windowDuration];
    
    % compute some metric over all trials
    switch metric
        case 'mean'
            rTrial = cellfun(@(a,b) mean(a(:,b>=thisTimeWindow(1)&b<thisTimeWindow(2)),2),data.trial,data.time,'uni',0);
        case 'rms'
            rTrial = cellfun(@(a,b) rms(a(:,b>=thisTimeWindow(1)&b<thisTimeWindow(2)),2),data.trial,data.time,'uni',0);
        case 'std'
            rTrial = cellfun(@(a,b) std(a(:,b>=thisTimeWindow(1)&b<thisTimeWindow(2)),0,2),data.trial,data.time,'uni',0);
        case 'max'
            rTrial = cellfun(@(a,b) max(a(:,b>=thisTimeWindow(1)&b<thisTimeWindow(2)),[],2),data.trial,data.time,'uni',0);
        case 'min'
            rTrial = cellfun(@(a,b) min(a(:,b>=thisTimeWindow(1)&b<thisTimeWindow(2)),[],2),data.trial,data.time,'uni',0);
        case 'maxabs'
            rTrial = cellfun(@(a,b) max(abs(a(:,b>=thisTimeWindow(1)&b<thisTimeWindow(2))),[],2),data.trial,data.time,'uni',0);
    end
    
    
    rTrial = cat(2,rTrial{:}); % -> channel x metric
    
    R = nan(length(uAttn),length(uFreq),size(rTrial,1),maxNumberReps);
    
    for i = 1:length(uFreq)
        for a = 1:length(uAttn)
            ind = Freq == uFreq(i) & Attn == uAttn(a);
            
            idx = find(ind,maxNumberReps);
            
            R(a,i,:,:) = rTrial(:,idx); % -> Attn x Freq x Channel x Repetition
        end
    end
    
    % compute mean over repetitions
    mR(:,:,:,tm) = mean(R,4); % -> Attn x Freq x Channel x Frame
    
        
    fprintf(' done\n')
end

if remMean
    mR = mR - mean(mR,4,'omitnan');
end

%% then plot


mav = squeeze(max(abs(mR),[],[1 2 4]));

nmstr = sprintf('%s - %s',data.cfg.experiment,metric);

f = figure('Name',nmstr,'Position',[1921 41 1920 963]);

clear FR
FR(length(tvec)) = struct('cdata',[],'colormap',[]);

uSL = baseSoundLevel - uAttn; % Attenuation -> Sound Level

for tm = 1:length(tvec)
    fprintf('Plotting frame %d of %d ...',tm,length(tvec))

    thisTimeWindow = tvec(tm)+[0 windowDuration];
    
    nmstr = sprintf('%s - %s %ssec',data.cfg.experiment,metric,mat2str(thisTimeWindow,3));
    
    f.Name = nmstr;
    f.NumberTitle = 'off';
    
    pos = layout.pos(1:end-2,:) + 1;
    pos(:,2) = max(pos(:,2)) - pos(:,2) + 1;
    pos = fix(pos);
    mlp = [7 4];
    
    t = tiledlayout(mlp(2),mlp(1));
    t.TileSpacing = 'compact';
    t.Padding = 'tight';
    

    
    ax = [];
    for ch = 1:size(mR,3)
        r = mR(:,:,ch,tm);
        r = r ./ mav(ch);
        
        imR = interp2(r,interpMultiplier,'makima');
        
        pidx = sub2ind(mlp,pos(ch,1),pos(ch,2));
        ax(ch) = nexttile(pidx);
        
        imagesc(uFreq,uSL,imR);
        
        title(data.label{ch});
        axis tight
    end
    
    set(ax,'ydir','normal','xscale','log');
    set(findobj(f,'-property','fontname'),'fontname','Consolas')
    
    
    %     cm = cell2mat(get(ax,'clim'));
    %     set(ax,'clim',[min(cm(:,1)) max(cm(:,2))]); % scale all plots together
    
    
    sgtitle(nmstr,'interpreter','none','fontname','Consolas');
    
    drawnow
    
    FR(tm) = getframe(f);
    
    fprintf(' done\n')
end

%% write video of animation

% note = '70-150Hz_BIP_abshilbert_highRes';
note = '2-150Hz_BIP';

vpath = 'G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG\';
vfn = sprintf('ToneRFanimation_%s_%s-%s',data.cfg.experiment,metric,note);

vffn = fullfile(vpath,vfn);

v = VideoWriter(vffn);

v.FrameRate = 10;
open(v);
for i= 1:length(FR)
    writeVideo(v,FR(i));
end
close(v);
