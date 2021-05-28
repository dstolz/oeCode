%% Stimulus response animation


interpMultiplier = 3;

remMean = false;
% remMean = true;

window = [0 0.25]; %s
windowStepSize = 0.001; % s
windowDuration = 0.001; % s


metric = 'mean';
% metric = 'rms';
% metric = 'std';
% metric = 'max';
% metric = 'min';
% metric = 'maxabs';






tvec = window(1):windowStepSize:window(2);


% Reps x Channel x Frame
mR = zeros(size(data.trial{1},1),length(tvec));

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
    
    
    R = cat(2,rTrial{:}); % -> channel x metric
    
    % compute mean over repetitions
    mR(:,tm) = mean(R,2); % -> Channel x Frame
    
        
    fprintf(' done\n')
end

if remMean
    mR = mR - mean(mR,2,'omitnan');
end

%% then plot


mav = max(abs(mR(:)));

mR = mR ./ mav;

nmstr = sprintf('%s - %s',data.cfg.experiment,metric);


clear FR
FR(length(tvec)) = struct('cdata',[],'colormap',[]);

fcfg = [];
fcfg.layout = layout;
fcfg.comment = 'no';
% tcfg.colorbar = 'EastOutside';
fcfg.shading = 'interp';
fcfg.interactive = 'no';

data_frame = data;
data_frame.time   = 0;
data_frame.dimord = 'chan_time';

for tm = 1:length(tvec)
    fprintf('Plotting frame %d of %d ...',tm,length(tvec))

    thisTimeWindow = tvec(tm)+[0 windowDuration];
    
    nmstr = sprintf('%s - %s %ssec',data.cfg.experiment,metric,mat2str(thisTimeWindow,3));
  
    
    data_frame.trial = mR(:,tm);
    
    
    ft_topoplotER(fcfg,data_frame);
    
    sgtitle(nmstr,'interpreter','none','fontname','Consolas');
    
    set(gca,'clim',[-1 1]*.5);
    
    set(gcf,'Name',nmstr,'Position',[2000 180 1063 673],'color','w');
        
    drawnow
    
    FR(tm) = getframe(gcf);
    
    close(gcf);
    fprintf(' done\n')
end

%% write video of animation

note = 'CSD_2-60Hz_1msbins_normalizedScale';


vpath = 'G:\My Drive\PROJECTS\VNS Targeted Platicity\ECoG\';
vfn = sprintf('ClickReponseAnimation_%s_%s-%s',data.cfg.experiment,metric,note);

vffn = fullfile(vpath,vfn);

v = VideoWriter(vffn);

v.FrameRate = 10;
open(v);
for i= 1:length(FR)
    writeVideo(v,FR(i));
end
close(v);
