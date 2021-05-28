function [trl,event] = ft_oe_trialfun(cfg)
% [trl,event] = ft_oe_trialfun(cfg)
% 
% cfg requires:
%    .path
%    .experiment
% 
% cfg optional:
%    .recording
%    .baphymfile
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
if ~isfield(cfg.trialdef,'skiptriggers'),  cfg.trialdef.skiptriggers = [];     end

if ~isfield(cfg,'baphystimformat'), cfg.baphystimformat = '%*s%s%*[^\n]'; end

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

if ~isempty(cfg.trialdef.skiptriggers)
    fprintf(2,'Skipping triggers: %s\n',mat2str(cfg.trialdef.skiptriggers))
    oeEventID(cfg.trialdef.skiptriggers)         = [];
    oeEventTimestamps(cfg.trialdef.skiptriggers) = [];
end

if rem(length(oeEventTimestamps),2)
    fprintf(2,'Odd number of triggers! Removing last trigger\n')
    oeEventTimestamps(end) = [];
    oeEventID(end) = [];
end

oeDuration        = oeEventTimestamps(2:2:end) - oeEventTimestamps(1:2:end);
oeEventTimestamps = oeEventTimestamps(1:2:end);
oeEventID         = oeEventID(1:2:end);
evOnS             = round(oeEventTimestamps.*cfg.fsample);



% Read Baphy stimulus event values
% NOTE: VALUES FROM BAPHY ARE STORED IN AN M FILE
if ~isfield(cfg,'baphymfile') || isempty(cfg.baphymfile)
    d = dir(fullfile(cfg.path,cfg.experiment,'*.m'));
    cfg.baphymfile = fullfile(d.folder,d.name);
end
run(cfg.baphymfile); % creates 'exptevents'

clear globalparams

% INCORPORATE BAPHY EVENTS. BAPHY DOES NOT HAVE A WAY TO LABEL THE STIMLUS
% VALUES, SO JUST MAKE A AD HOC SOLUTION.
bind = startsWith({exptevents.Note},'Stim','IgnoreCase',true); %#ok<*NODEF>
exptevents(~bind) = [];

% baphy has highly unpredictable formatting
baphyVals = cellfun(@(a) textscan(a,cfg.baphystimformat,'delimiter',{',',' ','[',']'}, ... 
    'CollectOutput',true,'MultipleDelimsAsOne',1), ...
    {exptevents.Note},'uni',1);

baphyVals = baphyVals';

if isnumeric(baphyVals{1})
    baphyVals = cell2mat(baphyVals);
elseif iscellstr(baphyVals{1})
    baphyVals = cellfun(@char,baphyVals,'uni',0);
end





if length(baphyVals) ~= length(oeEventID)
    n = min(size(baphyVals,1),length(oeEventID));
    warning('Number of Baphy events = %d\nNumber of OpenEphys events = %d\nUsing first %d events\n', ...
        size(baphyVals,1),length(oeEventID),n);
    baphyVals(n+1:end,:) = [];
    evOnS(n+1:end) = [];
    oeEventTimestamps(n+1:end) = [];
    oeDuration(n+1:end) = [];
end


% convert prestim and poststim parameters to samples
sprestim  = round(cfg.trialdef.prestim*cfg.fsample);
spoststim = round(cfg.trialdef.poststim*cfg.fsample);

trl(:,1) = evOnS-sprestim;
trl(:,2) = evOnS+spoststim-1;
trl(:,3) = -sprestim;
if size(baphyVals,2) > 1
    trl(:,4) = baphyVals(:,1) + abs(baphyVals(:,2)/100); % hack to code for multiple stimulus parameters
elseif isnumeric(baphyVals)
    trl(:,4) = baphyVals;
elseif iscellstr(baphyVals)
    % ad hoc for baphy nonsense
    trl(:,4) = cell2mat(cellfun(@(a) textscan(a,'%*s%*d%d%*s','delimiter','_'),baphyVals));    
end


if iscellstr(baphyVals)
    baphyVals = string(baphyVals);
end


event(size(baphyVals)) = struct('type',[],'sample',[],'value',[],'timestamp',[]);
for i = 1:size(trl,1)
    event(i).type   = 'Stim';%baphyTypes{ind};
    event(i).sample = evOnS(i);
    event(i).value  = baphyVals(i,:);
    event(i).timestamp = oeEventTimestamps(i);
    event(i).duration = oeDuration(i);
end




