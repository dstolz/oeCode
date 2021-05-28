function [trl,event] = ft_oe_trialfun_bhvr(cfg)
% [trl,event] = ft_oe_trialfun_bhvr(cfg)
% 
% cfg requires:
%    .path
%    .experiment
% 
% cfg optional:
%    .recording
%    .behaviormatfile
%    .channelmapfile
%
% Returns trl and event structures for use with FieldTrip
% cfg.trialdef.eventtype  ... name of an event
%    .trialdef.eventvalue ... value(s) of interest for the event
%    .trialdef.prestim    ... onset time (seconds) of window relative to
%                             event trigger (value > 0 means before trigger)
%    .trialdef.poststim   ... offset time (seconds) of window relative to event trigger
%    .trialdef.fsample    ... sampling rate to calculate event and trigger
%                             samples
%
%    .trialdef.skiptriggers ... trial triggers to skip, because Baphy is terrible
%


% Read OpenEphys event (TTL) data
% eventIDs seem to be TTL input channel # - 1
% All events are stored as on/off times synchronized with physiology data
%  event 0 : Lick IR TTL
%  event 1 : Shock TTL
%  event 3 : Trial TTL


if ~isfield(cfg.trialdef,'eventtype'),  cfg.trialdef.eventtype = [];     end
if ~isfield(cfg.trialdef,'eventvalue'), cfg.trialdef.eventvalue = [];    end
if ~isfield(cfg.trialdef,'eventid'),    cfg.trialdef.eventid = 3;        end
if ~isfield(cfg.trialdef,'prestim'),    cfg.trialdef.prestim  = 0;       end
if ~isfield(cfg.trialdef,'poststim'),   cfg.trialdef.poststim = 0.5;     end


if ~isfield(cfg,'behaviormatfile')
    d = dir(fullfile(cfg.path,cfg.experiment,'TrialData*.mat'));
    cfg.behaviormatfile = d(1).name;
end

if ~isfield(cfg,'recording')
    d = dir(fullfile(cfg.path,cfg.experiment,'Record Node*'));
    cfg.recording = d(1).name;
end

ffn = fullfile(cfg.path,cfg.experiment,cfg.recording,'all_channels.events');

[oeEventID,oeEventTimestamps,oeEventInfo] = load_open_ephys_data(ffn);

cfg.fsample = oeEventInfo.header.sampleRate;

eind = oeEventID == cfg.trialdef.eventid;
oeEventID(~eind) = [];
oeEventTimestamps(~eind) = [];


if rem(length(oeEventTimestamps),2)
    fprintf(2,'Odd number of triggers! Removing first trigger\n')
    oeEventTimestamps(1) = [];
    oeEventID(1) = [];
end

oeDuration        = oeEventTimestamps(2:2:end) - oeEventTimestamps(1:2:end);
oeEventTimestamps = oeEventTimestamps(1:2:end);
oeEventID         = oeEventID(1:2:end);
evOnS             = round(oeEventTimestamps.*cfg.fsample);



load(fullfile(cfg.path,cfg.experiment,cfg.behaviormatfile),'Data');
Data.trialData(:,2) = []; % get rid of anomalous empty data
%%%%%




% convert prestim and poststim parameters to samples
sprestim  = round(cfg.trialdef.prestim*cfg.fsample);
spoststim = round(cfg.trialdef.poststim*cfg.fsample);

trl(:,1) = evOnS-sprestim;
trl(:,2) = evOnS+spoststim-1;
trl(:,3) = -sprestim;


trialType     = [Data.trialData.TrialType]';
trialResponse = [Data.trialData.Response]';

event(size(Data.trialData)) = struct('type',[],'sample',[],'value',[],'timestamp',[],'trialData',[]);
for i = 1:size(trl,1)
    event(i).type   = 'Trial Type';
    event(i).sample = evOnS(i);
    event(i).value  = sprintf('%s\n%s',trialType(i),trialResponse(i));
    event(i).timestamp = oeEventTimestamps(i);
    event(i).duration = oeDuration(i);
    
    event(i).trialData = Data.trialData(i);
end




% event(size(oeEventID)) = struct('type',[],'sample',[],'value',[],'timestamp',[],'trialData',[]);
% ueid = unique(oeEventID);
% k = 1;
% for ie = 1:length(ueid)
%     idx = find(oeEventID == ueid(ie));
%     for i = 1:length(idx)
%         event(k).type   = eventTypes{ueid(ie)+1};
%         event(k).sample = evOnS(idx(i));
%         event(k).timestamp = oeEventTimestamps(idx(i));
%         event(k).duration = oeDuration(idx(i));
%         
%         if ueid(ie) == 3 % Trial code
%             event(k).value  = sprintf('%s\n%s',trialType(kk),trialResponse(kk));
%             event(k).trialData = Data.trialData(kk);
%             kk = kk + 1;
%         end
%         k = k + 1;
%     end
% end

