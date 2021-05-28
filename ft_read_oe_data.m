function [data,cfg] = ft_read_oe_data(cfg)
% [data,cfg] = ft_read_oe_data(cfg)
%
% cfg structure
%   .numchans   ... [1x1] expected number of channels to load. default = 32
%   .path       ... string path to data. default = pwd
%   .experiment ... string name of the experiment (typically a folder name
%                   in cfg.path). no default
%   .recording  ... string containing the specific recording. default is
%                   determined automatically.
%   .channelmapfile ... string filename specifying a channel map matlab
%                       file. The "Map" field in this structure is used to
%                       remap the channels if it is available. default = no
%                       file specified and no remapping takes place.
% 
% Output:
%   data    ... structure with trial data
%   cfg     ... cfg input structure is reflected back with additional
%               fields expected by the Fieldtrip Toolbox.
% 
% 
% Example:
%   cfg = [];
%   cfg.path = 'C:\Your\Path\to\ECoG_Data\';
%   cfg.experiment = 'Mangrove_2021-05-06_12-54-34_TORCS'; % OpenEphys data folder
%   data = ft_read_oe_data(cfg)
% 
% DJS 2021

if ~isfield(cfg,'numchans'), cfg.numchans = 32; end
if ~isfield(cfg,'path'), cfg.path = pwd; end
if ~isfield(cfg,'node'), cfg.node = '107'; end

if ~isfield(cfg,'recording')
    d = dir(fullfile(cfg.path,cfg.experiment,'Record Node*'));
    cfg.recording = d(1).name;
end

chFlag = 1; % sometimes "CH" is prepended to the file, and sometimes not ????
ffn = fullfile(cfg.path,cfg.experiment,cfg.recording,'*CH*.continuous');
d = dir(ffn);

if isempty(d)
    chFlag = 0;
    ffn = fullfile(cfg.path,cfg.experiment,cfg.recording,'*.continuous');
    d = dir(ffn);
end


for i = cfg.numchans:-1:1
    if chFlag
        efn = sprintf('%s_CH%d.continuous',cfg.node,i);
    else
        efn = sprintf('%s_%d.continuous',cfg.node,i);
    end
    ffn = fullfile(cfg.path,cfg.experiment,cfg.recording,efn);
    [data.trial{1}(i,:),timestamps,oeInfo] = load_open_ephys_data(ffn);    
    
    data.label{i,1} = sprintf('CH%02d',i);
    
    fprintf('\bdone\n')
end

data.fsample = oeInfo.header.sampleRate;


data.time{1} = linspace(timestamps(1),timestamps(end),size(data.trial{1},2)) - timestamps(1);
data.sampleinfo = 1;
data.sampleinfo(2) = data.sampleinfo + length(timestamps)-1;

firstIdx = round(timestamps(1).*data.fsample);

cfg.trl(:,[1 2]) = cfg.trl(:,[1 2]) - firstIdx + 1;


if isfield(cfg,'channelmapfile')
    fprintf('Remapping electrodes using "%s"\n',cfg.channelmapfile)
    load(cfg.channelmapfile)
    data.trial{1} = data.trial{1}(ChannelData.Map,:);
end



