%% NSMP: Neurophysiological substrates of movement preparation, TMS-EEG-EMG study
% ------------------------------------------------------------------------
% authors:  Dominika Sulcova, MSH, Hamburg
%           Mantosh Patnaik, UCLouvain, Brussels
% created:  January 2025   
% ------------------------------------------------------------------------
% project:  fill in project description
% 
% data:     1) raw EEG recordings (NeurOne, Bittium)
%           --> triggers:   1 - start of a trial
%                           2 - TMS stimulus onset
%                           3 - imperative cue (= bridge appears)
%                           4 - preparatory cue (= ball appears)
%           2) raw EMG recordings
%           --> triggers:   1 - TMS stimulus onset
% 
% script:
%   - part 1:   runs during the experimental session
%               encodes subject & session information
% 
%   - part 2:   processing of single-subject EEG data 
%               1) data import              
%                   - loads raw data from the subject folder --> dataset
%                   - crops close to trials
%                   - downsamples
%                   - DC + linear detrend on continuous data
%                   - saves for letswave
%               2) pre-process TEPs: letswave
%                   - interpolates channels if necessary
%                   - inputs missed trials
%                   - re-labels 'stimulation' events according to conditions
%                   - re-references to common average
%                   - epochs around TMS stimulus 
%                   - DC + linear detrend per epoch
%                   - removes TMS artifact via cubic interpolation 
%                   - concatenates datasets per condition
%                   - exports for EEGLAB
%               3) pre-process TEPs: EEGLAB
%                   - checks for bad trials
%                   - removes muscular artifact using SSP-SIR
%                   - corrects for baseline
%                   - quality check --> plots figure & saves
%                   - saves for letswave
%               4) pre-process pre-stimulus data
%                   - applies bandpass frequency filter
%                   - selects epochs from the pre-stimulus interval 
%                   - saves for letswave 
%               5) ICA
%                   - calculates ICA matrix for all datasets together
%                   - encodes discarded ICs
%                   - saves for letswave
%               6) extract TEP measures
%                   - calculates average TEP GFP and identifies peaks
%                   - plots average TEP with peak topographies
%                   - ...
%               7) extract beta measures 
%                   - ...
% 
%   - part 3:   processing of single-subject EMG data 
%               1) ...
% 
%   - part 4:   visualization of group-averaged data
%               1) ...
% 
% output:   output structures saved into the output MATLAB file 
%   - INFO   
%   - DATA
%   - MEASURES

%% ===================== PART 1: experimental session =====================
% directories
folder.output = uigetdir(pwd, 'Choose the OneDrive folder');        % output folder --> where the info structure is saved 
cd(folder.output)

% output
study = 'NSMP';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);

% load info structure
fprintf('loading the info structure...\n')
if exist(output_file) == 2
    output_vars = who('-file', output_file);
    if ismember('INFO', output_vars)
        load(output_file, 'INFO')
    else
        INFO = struct;
        save(output_file, 'INFO','-append')
    end
else
    INFO = struct;
    save(output_file, 'INFO')
end

% input subject & session info
prompt = {'subject', 'age:', 'sex:', 'handedness:', 'rMT:', 'session start:'};
dlgtitle = 'subject & session information';
dims = [1 35];
definput = {'S000', '00', 'female/male', 'right/left', '00', '00:00'};
info.session = inputdlg(prompt,dlgtitle,dims,definput);

% identify subject index
subject_idx = str2double(regexp(info.session{1}, '\d+', 'match', 'once'));

% update info structure
fprintf('encoding subject & session info...\n')
INFO(subject_idx).ID = info.session{1};
INFO(subject_idx).age = str2num(info.session{2});
INFO(subject_idx).sex = info.session{3};
INFO(subject_idx).handedness = info.session{4};
INFO(subject_idx).session.date = date;
INFO(subject_idx).session.start = info.session{6};
INFO(subject_idx).TMS.rMT = str2num(info.session{5});
INFO(subject_idx).TMS.intensity = ceil(INFO(subject_idx).TMS.rMT * 1.15); 

% input recording info
prompt = {'sampling rate (Hz):', 'reference:', 'ground:', 'triggers:', 'blocks:'};
dlgtitle = 'recording information';
dims = [1 50];
definput = {'20000', 'Fz', 'AFz', '1 - start, 2 - stimulation, 4 - imperative, 8 - preparatory', '1,2,3,4,5'};
info.recording = inputdlg(prompt,dlgtitle,dims,definput);

% update info structure
fprintf('encoding recording info...\n')
INFO(subject_idx).EEG.SR = str2num(info.recording{1});
INFO(subject_idx).EEG.ref = info.recording{2};
INFO(subject_idx).EEG.ground = info.recording{3};
for t = 1:4
    INFO(subject_idx).EEG.triggers(t).trigger = t;
    pattern = sprintf('%d - ([a-zA-Z]+)', t);
    INFO(subject_idx).EEG.triggers(t).label = regexp(info.recording{4}, pattern, 'tokens', 'once');
end
INFO(subject_idx).EEG.blocks = str2num(info.recording{5});

% save and continue
save(output_file, 'INFO','-append')
fprintf('done.\n\n')
clear output_vars info prompt dlgtitle dims definput pattern t

%% ====================== PART 2: EEG pre-processing ======================
% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % MATLAB toolboxes
folder.raw = uigetdir(pwd, 'Choose the input folder');              % raw data from one session --> hard-drive or local folder
folder.processed = uigetdir(pwd, 'Choose the data folder');         % processed data --> local folder 
folder.output = uigetdir(pwd, 'Choose the output folder');          % output folder --> OneDrive folder
cd(folder.processed)

% output
study = 'NSMP';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);

% load output structures
fprintf('loading the info structure...\n')
if exist(output_file) == 2
    output_vars = who('-file', output_file);
    if ismember('INFO', output_vars)
        load(output_file, 'INFO')
    else
        INFO = struct;
        save(output_file, 'INFO','-append')
    end
    if ismember('DATA', output_vars)
        load(output_file, 'DATA')
    else
        DATA = struct;
        save(output_file, 'DATA','-append')
    end
    if ismember('MEASURES', output_vars)
        load(output_file, 'MEASURES')
    else
        MEASURES = struct;
        save(output_file, 'MEASURES','-append')
    end
else
    INFO = struct; DATA = struct; MEASURES = struct;
    save(output_file, 'INFO', 'DATA', 'MEASURES')
end

% dataset
params.conditions = {'baseline' 'delay' 'no_TMS' 'catch'};

% current participant
prompt = {'current subject number:'};
dlgtitle = 'subject';
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);
subject_idx = str2num(input{1,1});

% visualization
figure_counter = 1;
clear output_vars prompt dlgtitle dims definput input

%% 1) import & pre-process continuous data
% ----- section input -----
params.suffix = {'crop' 'ds' 'dc'};
params.crop_margin = 5;
params.downsample = 20;
% ------------------------- 
fprintf('section 1: import & pre-process continuous data\n')

% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));

% identify import folder and subfolders
session_folders = dir(sprintf('%s\\%s', folder.raw, INFO(subject_idx).ID));
session_date = datetime(INFO(subject_idx).session.date, 'InputFormat', 'dd-MMM-yyyy', 'Format', 'yyyy-MM-dd');
session_date = char(session_date);
for a = 1:length(session_folders)
    if contains(session_folders(a).name, session_date)
        params.folder = sprintf('%s\\%s', session_folders(a).folder, session_folders(a).name);
    end
end
params.blocks = INFO(subject_idx).EEG.blocks;

% check if identified subfolders are avaiable 
data2import = dir(params.folder);
if ~isempty(data2import)
    % remove all datasets that are not labelled with numbers
    data_idx = logical([]);
    for b = 1:length(data2import)
        if isempty(str2num(data2import(b).name))
            data_idx(b) = true;
        else
            data2import(b).label = str2num(data2import(b).name);
            data_idx(b) = false;
        end
    end
    data2import(data_idx) = [];
    [~, file_idx] = sort([data2import.label]);
    data2import = data2import(file_idx);

    % check for correct labels
    file_idx = false(1, length(params.blocks));
    for b = 1:length(params.blocks)
        for c = 1:length(data2import)
            if params.blocks(b) == data2import(c).label
                file_idx(b) = true;
            end
        end
    end
    data2import = data2import(file_idx);

    % check number of datasets
    fprintf('%d datasets found in the directory.\n', length(data2import))
    if length(data2import) ~= length(params.blocks)
        error('ERROR: This does not match with expected number of datasets (%d)\n.Please verify manually.\n', length(params.blocks))
    end
else
    error('ERROR: no datasets found in the directory.\n')
end

% load datasets
fprintf('loading datasets:\n')
for d = 1:length(data2import)
    % provide update
    fprintf('--> block %d - ', d)

    % create the name
    params.name = sprintf('%s %s b%d', study, INFO(subject_idx).ID, d);
    
    % encode to the info structure
    if d == 1
        INFO(subject_idx).EEG.dataset(1).block = d;
    else
        INFO(subject_idx).EEG.dataset(end + 1).block = d;
    end
    INFO(subject_idx).EEG.dataset(end).subfolder = data2import(d).label;
    INFO(subject_idx).EEG.dataset(end).name = params.name;  

    % import the dataset
    [dataset.raw(d).header, dataset.raw(d).data, ~] = RLW_import_MEGA(data2import(d).folder, data2import(d).block);

    % rename in the header
    dataset.raw(d).header.name = params.name;
end  
fprintf('done.\n')

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% pre-process continuous data and save for letswave
fprintf('pre-processing:\n')
for d = 2:length(dataset.raw)
    % provide update
    fprintf('--> %s\n', INFO(subject_idx).EEG.dataset(d).name)

    % select data
    lwdata.header = dataset.raw(d).header;
    lwdata.data = dataset.raw(d).data; 

    % assign electrode coordinates
    fprintf('assigning electrode coordinates...\n')
    option = struct('filepath', sprintf('%s\\letswave 7\\res\\electrodes\\spherical_locations\\Standard-10-20-Cap81.locs', folder.toolbox), ...
        'suffix', '', 'is_save', 0);
    lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
    if d == 1
        INFO(subject_idx).EEG.processing(1).process = 'electrode coordinates assigned';
        INFO(subject_idx).EEG.processing(1).params.layout = 'standard 10-20-cap81';
        INFO(subject_idx).EEG.processing(1).suffix = [];
        INFO(subject_idx).EEG.processing(1).date = sprintf('%s', date);
    end

    % re-label and count events
    fprintf('checking events...\n') 
    event_idx = logical([]);
    for a = 1:length(lwdata.header.events)
        if isempty(str2num(lwdata.header.events(a).code))
            event_idx(a) = true;
        else
            event_idx(a) = false;
        end
    end
    lwdata.header.events(event_idx) = [];
    params.eventcodes = unique({lwdata.header.events.code});
    if length(params.eventcodes) == length(INFO(subject_idx).EEG.triggers)
        event_count = zeros(length(INFO(subject_idx).EEG.triggers), 1);
        for e = 1:length(lwdata.header.events)
            for a = 1:length(INFO(subject_idx).EEG.triggers)
                if strcmp(lwdata.header.events(e).code, num2str(INFO(subject_idx).EEG.triggers(a).trigger))
                    lwdata.header.events(e).code = INFO(subject_idx).EEG.triggers(a).label{1};
                    event_count(a) = event_count(a) + 1;
                end
            end
        end
    else
        error('ERROR: wrong number of triggers (%d) was found in the dataset!', length(params.eventcodes))
    end
    
    % update & encode
    fprintf('%d events in total were found in the dataset:\n%s - %d events\n%s - %d events\n%s - %d events\n%s - %d events\n', ...
        length(lwdata.header.events), ...
        INFO(subject_idx).EEG.triggers(1).label{1}, event_count(1), ...
        INFO(subject_idx).EEG.triggers(2).label{1}, event_count(2), ...
        INFO(subject_idx).EEG.triggers(3).label{1}, event_count(3), ...
        INFO(subject_idx).EEG.triggers(4).label{1}, event_count(4));
    INFO(subject_idx).EEG.dataset(d).trials = event_count(1);

    % crop data
    params.crop(1) = lwdata.header.events(1).latency - params.crop_margin;
    params.crop(2) = lwdata.header.events(end).latency + params.crop_margin;
    fprintf('cropping ...\n')
    option = struct('xcrop_chk', 1, 'xstart', params.crop(1), 'xend', params.crop(2), ...
        'suffix', params.suffix{1}, 'is_save', 0);
    lwdata = FLW_crop_epochs.get_lwdata(lwdata, option);
    if d == 1
        INFO(subject_idx).EEG.processing(2).process = 'continuous data cropped';
        INFO(subject_idx).EEG.processing(2).params.start = params.crop(1);
        INFO(subject_idx).EEG.processing(2).params.end = params.crop(2);
        INFO(subject_idx).EEG.processing(2).params.margin = params.crop_margin;
        INFO(subject_idx).EEG.processing(2).suffix = params.suffix{1};
        INFO(subject_idx).EEG.processing(2).date = sprintf('%s', date);
    end

    % downsample 
    fprintf('downsampling...\n')
    option = struct('x_dsratio', params.downsample, 'suffix', params.suffix{2}, 'is_save', 0);
    lwdata = FLW_downsample.get_lwdata(lwdata, option);
    if d == 1
        INFO(subject_idx).EEG.processing(3).process = sprintf('downsampled');
        INFO(subject_idx).EEG.processing(3).params.ratio = params.downsample;
        INFO(subject_idx).EEG.processing(3).params.fs_orig = 1/lwdata.header.xstep * params.downsample;
        INFO(subject_idx).EEG.processing(3).params.fs_final = 1/lwdata.header.xstep;
        INFO(subject_idx).EEG.processing(3).suffix = params.suffix{2};
        INFO(subject_idx).EEG.processing(3).date = sprintf('%s', date);
    end

    % remove DC + linear detrend continuous data
    fprintf('removing DC and applying linear detrend...')
    option = struct('linear_detrend', 1, 'suffix', params.suffix{3}, 'is_save', 1);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
    if d == 1
        INFO(subject_idx).EEG.processing(4).process = sprintf('DC + linear detrend on continuous data');
        INFO(subject_idx).EEG.processing(4).suffix = params.suffix{3};
        INFO(subject_idx).EEG.processing(4).date = sprintf('%s', date);
    end
    fprintf('done.\n')

    % update dataset
    dataset.raw(d).header = lwdata.header;
    dataset.raw(d).data = lwdata.data; 
end
fprintf('done.\n')

% save and continue
save(output_file, 'INFO','-append')
clear a b c d e session_folders session_date data2import data_idx file_idx  prompt dlgtitle dims definput event_idx event_count option lwdata 
fprintf('section 1 finished.\n\n')

%% 2) pre-process TEPs: letswave
% ----- section input -----
params.prefix = 'dc ds crop';
params.suffix = {'reref' 'dc' 'artifact' 'processed'};
params.interp_chans = 6;
params.epoch = [-1.5, 1.5];
params.artifact_interp = [-0.005 0.01];
params.ref = 'averef';
% -------------------------
fprintf('section 2: pre-process TEP data\n')

% load dataset if needed
if exist('dataset') ~= 1
    data2load = dir(sprintf('%s*%s*', params.prefix, INFO(subject_idx).ID));
    if length(data2load) == length(params.conditions) * 2
        dataset = reload_dataset(data2load, [], 'raw');
    else
        error('ERROR: Wrong number of available datasets to load! Check manually.')
    end
end

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% interpolate channels if needed
params.labels = {dataset.raw(1).header.chanlocs.labels};
prompt = {'indicate bad channels that need to be interpolated:'};
dlgtitle = 'channel interpolation';
dims = [1 120];
definput = {strjoin(params.labels,' ')};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if ~isempty(answer{1})
    % identify channels to interpolate
    chans2interpolate = split(answer{a}, ' ');

    % interpolate identified channels in all datsets
    for c = 1:length(chans2interpolate)
        if ~isempty(chans2interpolate{c})
            % provide update
            fprintf('interpolating channel %s\n', chans2interpolate{c})

            % indentify the channel to interpolate
            chan_n = find(strcmp(params.labels, chans2interpolate{c}));

            % calculate distances with other electrodes
            chan_dist = -ones(length(dataset.raw(1).header.chanlocs), 1);
            for b = setdiff(1:length(dataset.raw(1).header.chanlocs), chan_n)
                if dataset.raw(1).header.chanlocs(b).topo_enabled == 1
                    chan_dist(b) = sqrt((dataset.raw(1).header.chanlocs(b).X - dataset.raw(1).header.chanlocs(chan_n).X)^2 + ...
                        (dataset.raw(1).header.chanlocs(b).Y - dataset.raw(1).header.chanlocs(chan_n).Y)^2 + ...
                        (dataset.raw(1).header.chanlocs(b).Z - dataset.raw(1).header.chanlocs(chan_n).Z)^2);
                end
            end
            chan_dist((chan_dist==-1)) = max(chan_dist);
            [~,chan_idx] = sort(chan_dist);

            % identify neighbouring channels
            chan_idx = chan_idx(1:params.interp_chans);
            chans2use = params.labels;
            chans2use = chans2use(chan_idx);

            % cycle through all datasets
            for d = 1:length(dataset.raw)
                % select data
                lwdata.header = dataset.raw(d).header;
                lwdata.data = dataset.raw(d).data;
    
                % interpolate using the neighboring electrodes
                option = struct('channel_to_interpolate', chans2interpolate{c}, 'channels_for_interpolation_list', {chans2use}, ...
                    'suffix', '', 'is_save', 0);
                lwdata = FLW_interpolate_channel.get_lwdata(lwdata, option);
    
                % update dataset
                dataset.raw(d).header = lwdata.header;
                dataset.raw(d).data = lwdata.data;  
            end
            
            % encode
            if c == 1
                INFO(subject_idx).EEG.processing(5).process = sprintf('bad channels interpolated');
                INFO(subject_idx).EEG.processing(5).date = sprintf('%s', date);
            end
            INFO(subject_idx).EEG.processing(5).params.bad{c} = chans2interpolate{c};
            INFO(subject_idx).EEG.processing(5).params.chans_used{c} = strjoin(chans2use, ' ');  
        end
    end
else
    INFO(subject_idx).EEG.processing(5).process = sprintf('no channels interpolated');
    INFO(subject_idx).EEG.processing(5).date = sprintf('%s', date);
end

% input missed channels
answer;

% re-label 'stimulation' events according to conditions
fprintf('assigning condition labels to events: block ')

% epoch and sort into conditions
fprintf('pre-processing:\n')
for d = 1:length(dataset.raw)
    fprintf('block %d:\n', d)

    % select data
    lwdata.header = dataset.raw(d).header;
    lwdata.data = dataset.raw(d).data;

    % re-reference to common average
    fprintf('re-referencing...')
    option = struct('reference_list', {params.labels}, 'apply_list', {params.labels}, 'suffix', params.suffix{1}, 'is_save', 0);
    lwdata = FLW_rereference.get_lwdata(lwdata, option);
    if d == 1
        INFO(subject_idx).EEG.processing(6).process = sprintf('re-referenced to common average');
        INFO(subject_idx).EEG.processing(6).suffix = params.suffix{1};
        INFO(subject_idx).EEG.processing(6).date = sprintf('%s', date);
    end

    % epoch per condition
    fprintf('epoching ... ')
    lwdata_orig = lwdata;
    for c = 1:length(params.conditions) 
        % add letswave 7 to the top of search path
        addpath(genpath([folder.toolbox '\letswave 7']));

        % subset the data of this condition and epoch
        lwdata = lwdata_orig;
        option = struct('event_labels', params.conditions{c}, 'x_start', params.epoch(1), 'x_end', params.epoch(2), ...
            'x_duration', params.epoch(2)-params.epoch(1), 'suffix', sprints('%s b%d', params.conditions{c}, d), 'is_save', 0);
        lwdata = FLW_segmentation.get_lwdata(lwdata, option);
        if d == 1 && c == 1
            INFO(subject_idx).EEG.processing(7).process = sprintf('segmented to TEP epochs');
            INFO(subject_idx).EEG.processing(7).params.limits = params.epoch;
            INFO(subject_idx).EEG.processing(7).suffix = params.conditions;
            INFO(subject_idx).EEG.processing(7).date = sprintf('%s', date);
        end

        % remove DC + linear detrend epoch data
        fprintf('removing DC and applying linear detrend...')
        option = struct('linear_detrend', 1, 'suffix', params.suffix{2}, 'is_save', 0);
        lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
        if d == 1 && c == 1
            INFO(subject_idx).EEG.processing(8).process = sprintf('DC + linear detrend on epoched data');
            INFO(subject_idx).EEG.processing(8).suffix = params.suffix{2};
            INFO(subject_idx).EEG.processing(8).date = sprintf('%s', date);
        end

        % add letswave 6 to the top of search path
        addpath(genpath([folder.toolbox '\letswave 6']));

        % interpolate around TMS stimulus
        fprintf('interpolating TMS artifact...')
        [lwdata.header, lwdata.data, ~] = RLW_suppress_artifact_event(lwdata.header, lwdata.data, ...
            'xstart', params.artifact_interp(1), 'xend', params.artifact_interp(2), ...
            'event_code', params.conditions{c}, 'interp_method', 'pchip'); 
        lwdata.header.name = [params.suffix{3} ' ' lwdata.header.name];
        if d == 1 && c == 1
            INFO(subject_idx).EEG.processing(9).process = sprintf('TMS artifact interpolated');
            INFO(subject_idx).EEG.processing(9).params.limits = params.artifact_interp;
            INFO(subject_idx).EEG.processing(9).params.method = 'pchip';
            INFO(subject_idx).EEG.processing(9).suffix = params.suffix{3};
            INFO(subject_idx).EEG.processing(9).date = sprintf('%s', date);
        end

        % update in the dataset
        dataset.epoched((d - 1)*length(dataset.raw) + c).condition = params.conditions{c};
        dataset.epoched((d - 1)*length(dataset.raw) + c).header = lwdata.header;
        dataset.epoched((d - 1)*length(dataset.raw) + c).data = lwdata.data; 
        fprintf('done.\n')
    end
end

% concatenate datasets per condition & save for letswave
for c = 1:length(params.conditions)
    % select data from all matching datasets
    dataset.conditions(c).condition = params.conditions{c};
    dataset.conditions(c).header = 0;
    for d = 1:length(dataset.epoched)
        if strcmp(dataset.epoched(d).condition{1}, params.conditions{c})
            if dataset.conditions(c).header == 0
                dataset.conditions(c).header = dataset.epoched(d).header;
                dataset.conditions(c).data = dataset.epoched(d).data;
            else
                dataset.conditions(c).data = cat(1, dataset.conditions(c).data, dataset.epoched(d).data);
            end
        end
    end

    % adjust header
    dataset.conditions(c).header.name = dataset.conditions(c).header.name(1:end-3);
    dataset.conditions(c).header.datasize = size(dataset.conditions(c).data);

    % save
    header = dataset.conditions(c).header;
    save([dataset.conditions(c).header.name '.lw6'], 'header')
    data = dataset.conditions(c).data;
    save([dataset.conditions(c).header.name '.mat'], 'data')
end

% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% export in .set format for EEGLAB
fprintf('exporting for EEGLAB:\n')
for c = 1:length(params.conditions)
    fprintf('%s dataset ... ', params.conditions{c})
    
    % select data
    lwdata.header = dataset.conditions(c).header;
    lwdata.data = dataset.conditions(c).data;    
    
    % export 
    name = sprintf('%s %s %s', params.suffix{4}, INFO(subject_idx).ID, params.conditions{c});
    export_EEGLAB(lwdata, name, INFO(subject_idx).ID);
end
INFO(subject_idx).EEG.processing(10).process = sprintf('data exported for EEGLAB');
INFO(subject_idx).EEG.processing(10).params.format = '.set';
INFO(subject_idx).EEG.processing(10).suffix = params.suffix{4};
INFO(subject_idx).EEG.processing(10).date = sprintf('%s', date);
fprintf('done.\n')

% save and continue
save(output_file, 'INFO','-append')
clear params c d prompt dlgtitle dims definput answer chans2interpolate chan_n chan_dist chan_idx chans2use ...
    lwdata data header name
fprintf('section 2 finished.\n\n')

%% 3) pre-process TEPs: EEGLAB
% ----- section input -----
params.conditions = {'' ''};
params.suffix = {'processed' 'checked' 'sspsir_1' 'sspsir_2'};
params.baseline = [-0.25 -0.006]; 
params.time_range = [-0.006, 0.050];
% params.time_range = time_range/1000;              % uncomment this if the filter is not applied to appropriate time interval 
params.lwsave = true;                               % --> in that case a following change must be made in tesa_spssir function  
% -------------------------                         % in the calculation of the smoothening function:   if [timeRange(tidx,2) - timeRange(tidx,1)] < 1
fprintf('section 3: TEP pre-processing - EEGLAB\n')                                                     %   smoothLength = smoothLength / 1000;
                                                                                                        % end
% add eeglab to the top of search path and launch
addpath(fullfile(folder.toolbox, 'EEGLAB'));
eeglab 

% update figure counter
fig_all = findall(0, 'Type', 'figure');
figure_counter = length(fig_all) + 1;

% remove bad epochs
fprintf('removing bad epochs:\n')
for c = 1:length(params.conditions)
    fprintf('%s dataset ...', params.conditions{c})

    % load the dataset, re-refence
    name = sprintf('%s %s %s.set', params.suffix{1}, INFO(subject_idx).ID, params.conditions{c}); 
    EEG = pop_loadset('filename', name, 'filepath', folder.processed);
    EEG = pop_reref(EEG, []);
    EEG.filename = name;
    EEG.filepath = folder.processed;
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw  
    if a == 1 && b == 1
        INFO(subject_idx).processing(10).process = sprintf('re-referenced to common average');
        INFO(subject_idx).processing(10).date = sprintf('%s', date);
    end
    
    % visualize the average response
    figure(figure_counter); 
    pop_plottopo(pop_select(EEG, 'time', [-0.1 0.3]), [] , '', 0, 'ydir', 1, 'ylim', [-30 30]);
    sgtitle(sprintf('%s: %s', INFO(subject_idx).ID, params.conditions{c}))
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1])
    figure_counter = figure_counter + 1;
    
    % remove bad epochs
    pop_eegplot(EEG, 1, 1, 1);
    waitfor(gcf); 
    EEG = eeg_checkset(EEG);
    
    % encode bad epochs to info structure
    if a == 1 && b == 1
        INFO(subject_idx).processing(11).process = sprintf('bad trials discarded');
        INFO(subject_idx).processing(11).params.GUI = 'EEGLAB';
        INFO(subject_idx).processing(11).date = sprintf('%s', date);
    end
    for a = 1:length(ALLCOM)
        if contains(ALLCOM{a}, 'pop_rejepoch')
            match = regexp(ALLCOM{a}, '\[(.*?)\]', 'match');
            if ~isempty(match)
                discarded = str2num(match{1}(2:end-1)); 
            else
                discarded = []; 
            end
            break
        end
    end
    INFO(subject_idx).processing(11).params.discarded{c} = discarded;
      
    % save dataset
    pop_saveset(EEG, 'filename', EEG.filename, 'filepath', EEG.filepath);
end

% apply SSP-SIR individually to each dataset and save the filters
fprintf('applying SSP-SIR to the original dataset:\n')
for c = 1:length(params.conditions)
    fprintf('%s dataset ...', params.conditions{c})

    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
        
    % visual check - before SSP-SIR
    figure; 
    pop_timtopo(EEG, [-50  150], [5  15  30  45  60  100]);
    sgtitle(sprintf('%s - %s: original data', INFO(subject_idx).ID, params.conditions{c}))
    figure_name = sprintf('TEP no_filter %s %s', INFO(subject_idx).ID, params.conditions{c}); 
    savefig(sprintf('%s\\figures\\%s.fig', folder.output, figure_name))
    saveas(gcf, sprintf('%s\\figures\\%s.svg', folder.output, figure_name))

    % ask if the original ssp filter should be applied
    answer{counter} = questdlg('Do you want to apply original SSP-SIR', 'SSP-SIR', 'YES', 'NO', 'YES'); 

    % SSP-SIR - spherical model 
    [EEG, EEG_old] = pop_tesa_sspsir(EEG, 'artScale', 'manual', 'timeRange', time_range, 'PC', []);
    switch answer{counter}
        case 'YES'
            clear EEG_old
        case 'NO'
            EEG = EEG_old;
            clear EEG_old
    end  

    % baseline correct and save
    EEG = pop_rmbase(EEG, params.baseline, []);
    name = sprintf('%s %s %s.set', params.suffix{3}, INFO(subject_idx).ID, params.conditions{c}); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, c, 'setname', name, 'overwrite', 'on', 'gui', 'off');
    EEG.filename = name;
    eeglab redraw
    
    % extract information about rejected components
    sspsir(c) = EEG.SSPSIR;
end
fprint('done.\n')

% update info structure and save
INFO(subject_idx).EEG.processing(12).process = sprintf('muscular artifact filtered out');
INFO(subject_idx).EEG.processing(12).params.method = 'SSP-SIR';
INFO(subject_idx).EEG.processing(12).params.specifics = 'all filters applied to datasets of all conditions';
for c = 1:length(params.conditions)
    INFO(subject_idx).EEG.processing(12).params.PC(c).condition = params.conditions{c};
    INFO(subject_idx).EEG.processing(12).params.PC(c).PC = sspsir(c).PC;
end
INFO(subject_idx).EEG.processing(12).params.setings = rmfield(sspsir, 'PC');
INFO(subject_idx).EEG.processing(12).suffix = params.suffix([3,4]);
INFO(subject_idx).EEG.processing(12).date = sprintf('%s', date);
save(output_file, 'INFO','-append')

% apply SSP-SIR filters of other datasets
fprintf('applying each filter to all other datasets:\n')
for c = 1:length(params.conditions)
    fprintf('%s dataset ...', params.conditions{c})

    % select the dataset
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
            
    % SSP-SIR 
    other_datasets = 1:length(params.conditions);
    other_datasets = other_datasets(other_datasets ~= c);
    for a = other_datasets
        [EEG_trials] = SSP_SIR_trials(EEG, INFO(subject_idx).EEG.processing(12).params.setings(c).L_ave, ...
            INFO(subject_idx).EEG.processing(12).params.setings(a).topographies, ...
            INFO(subject_idx).EEG.processing(12).params.setings(a).filter, []);
        EEG.data =  EEG_trials.data;
    end
    
    % baseline correct and save
    EEG = pop_rmbase(EEG, params.baseline, []);
    name = sprintf('%s %s %s.set', params.suffix{4}, INFO(subject_idx).ID, params.conditions{c});
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, c, 'setname', name, 'overwrite', 'on', 'gui', 'off'); 
    eeglab redraw
    
    % visual check - after SSP-SIR
    figure; 
    pop_timtopo(EEG, [-50  150], [5  15  30  45  60  100]);
    sgtitle(sprintf('%s - %s: filtered data',  INFO(subject_idx).ID, params.conditions{c}))
    figure_name = sprintf('TEP SSPSIR %s %s', INFO(subject_idx).ID, params.conditions{c});
    savefig(sprintf('%s\\figures\\%s.fig', folder.output, figure_name))
    saveas(gcf, sprintf('%s\\figures\\%s.svg', folder.output, figure_name))
       
    % save dataset
    pop_saveset(EEG, 'filename', name, 'filepath', folder.processed);
end
fprint('done.\n')

% add letswave 7 to the top of search path
addpath(genpath([folder_toolbox '\letswave7-master']));

% export data and header
fprintf('exporting data back to letswave:\n')
for c = 1:length(params.conditions)
    fprintf('%s dataset ...', params.conditions{c})

    % select the dataset in EEGLAB
    EEG = ALLEEG(c);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, c);
    eeglab redraw
    
    % export in the letswave formate and save if required
    name = sprintf('%s %s %s', params.suffix{4}, INFO(subject_idx).ID, params.conditions{c});
    lwdata = export_lw(EEG, dataset.conditions(c).header, name);
    dataset.sspsir(c).condition = params.conditions{c};
    dataset.sspsir(c).header = lwdata.header;
    dataset.sspsir(c).data = lwdata.data;
    if params.lwsave
        CLW_save([], lwdata.header, lwdata.data);
    end
end
fprintf('done.\n')

% close EEGLAB and clear     
close eeglab
clear params c t e i name char_start char_end epochs_start epochs_end figure_name sspsir other_datasets ...
    data header lwdata ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG globalvars LASTCOM PLUGINLIST STUDY
fprintf('section 3 finished.\n\n')

%% 4) pre-process continuous data             
% - applies bandpass frequency filter
% - selects epochs from the pre-stimulus interval 
% - saves for letswave 
% ----- section input -----
% -------------------------
fprintf('section 4: pre-process continuous data\n')
fprintf('section 4 finished.\n\n')

%% 5) ICA
% - applies bandpass frequency filter
% - calculates ICA matrix for all datasets together
% - encodes discarded ICs
% - saves for letswave
% ----- section input -----
% -------------------------
fprintf('section 5: ICA\n')
fprintf('section 5 finished.\n\n')

%% 6) extract TEP measures
% - calculates average TEP GFP and identifies peaks
% - plots average TEP with peak topographies
% - ...
% ----- section input -----
% -------------------------
fprintf('section 6: extract TEP measures\n')
fprintf('section 6 finished.\n\n')

%% 7) extract beta measures
% ----- section input -----
% -------------------------
fprintf('section 7: extract beta measures\n')
fprintf('section 7 finished.\n\n')

%% ====================== PART 3: EMG pre-processing ======================

%% ======================== PART 4: visualization =========================

%% functions
function dataset = reload_dataset(data2load, conditions, fieldname)
% =========================================================================
% reloads pre-processed EEG data of a single subject for following 
% processing steps 
% input:    - data2load = list of datasets to load
%           - conditions = cell array with conditions
%           - fieldname
% =========================================================================  
% initiate output
dataset = struct;

% sort datasets
header_idx = logical([]);
data_idx = logical([]);
for d = 1:length(data2load)
    if contains(data2load(d).name, 'lw6') 
        header_idx(d) = true;
        data_idx(d) = false;
    elseif contains(data2load(d).name, 'mat') 
        header_idx(d) = false;
        data_idx(d) = true;
    end
end
headers = data2load(header_idx);
datas = data2load(data_idx);

% check conditions
if ~isempty(conditions)
    if length(conditions) == length(headers) 
    else
       conditions = []; 
       fprintf('WARNING: number of conditions does not match the number of found datasets.')
    end
end

% load 
if length(datas) == length(headers) 
    for d = 1:length(datas)
        % condition, if indicated
        if ~isempty(conditions)
            statement = sprintf('dataset.%s(d).condition = conditions{d};', fieldname);
            eval(statement) 
        end

        % header
        load(sprintf('%s\\%s', headers(d).folder, headers(d).name), '-mat')
        statement = sprintf('dataset.%s(d).header = header;', fieldname);
        eval(statement) 

        % data
        load(sprintf('%s\\%s', datas(d).folder, datas(d).name))
        statement = sprintf('dataset.%s(d).data = data;', fieldname);
        eval(statement) 
    end
else
    error('ERROR: Wrong number of available datasets to load! Check manually.')
end

end
function plot_ERP(input, varargin)
% =========================================================================
% plots an event-related potential
% input = structure with fields:    
%           data --> condition/electrode * sample
%           x --> vector with time samples
%           CI_upper --> condition/electrode * sample
%           CI_lower --> condition/electrode * sample
% varargins = name-value pairs: 
%           xlim --> 2-element vector (min, max)     
%           ylim --> 2-element vector (min, max) 
%           colours --> n*3 matrix of RGB values
%           shading --> 'on'(default)/'off'
%           alpha --> a float (default 0.2)           
%           plot_legend --> 'on'(default)/'off'
%           labels --> cell array with labels for the legend  
%           legend_loc --> legend location (default 'southeast')
%           eoi --> label of a channel to be highlighted
%           reverse --> 'on'/'off'(default) - flips y axis
% =========================================================================  
% set defaults
x_limits = [0,0];
y_limits = [0,0];
col = prism(size(input.data, 1));
shading = true;
alpha = 0.2;
plot_legend = true;
for c = 1:size(input.data, 1)
    labels{c} = sprintf('condition %d', c);
end
legend_loc = 'southeast';
highlight = false;
reverse = false;

% check for varargins
if ~isempty(varargin)
    % x limits
    a = find(strcmpi(varargin, 'xlim'));
    if ~isempty(a)
        x_limits = varargin{a + 1};
    end

    % y limits
    b = find(strcmpi(varargin, 'ylim'));
    if ~isempty(b)
        y_limits = varargin{b + 1};
    end

    % colours
    c = find(strcmpi(varargin, 'colours'));
    if ~isempty(c)
        col = varargin{c + 1};
    end

    % shading - default on
    d = find(strcmpi(varargin, 'shading'));
    if ~isempty(d) && strcmp(varargin{d + 1}, 'off')
        shading = false;
    end

    % alpha
    e = find(strcmpi(varargin, 'alpha'));
    if ~isempty(e)
        alpha = varargin{e + 1};
    end

    % legend - default on
    f = find(strcmpi(varargin, 'legend'));
    if ~isempty(f) && strcmp(varargin{f + 1}, 'off')
        plot_legend = false;
    end    

    % labels
    g = find(strcmpi(varargin, 'labels'));
    if ~isempty(g)
        labels = varargin{g + 1};
    end

    % legend location
    h = find(strcmpi(varargin, 'legend_loc'));
    if ~isempty(h) 
        legend_loc = varargin{h + 1};
    end  

    % highlighted channel - default off
    i = find(strcmpi(varargin, 'eoi'));
    if ~isempty(i)
        eoi = varargin{i + 1};
        highlight = true;
    end 

    % reverse y axis - default off
    r = find(strcmpi(varargin, 'reverse'));
    if ~isempty(r) && strcmp(varargin{r + 1}, 'on')
        reverse = true;
    end
end

% loop through datasets to plot
for t = 1:size(input.data, 1) 
    P(t) = plot(input.x, input.data(t, :), 'Color', col(t, :), 'LineWidth', 2);
    hold on
    if shading
        F(t) = fill([input.x fliplr(input.x)],[input.CI_upper(t, :) fliplr(input.CI_lower(t, :))], ...
            col(t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        hold on
    end
end

% check y limits
if y_limits(1) == 0 && y_limits(2) == 0
    y_limits = ylim;
end

% plot stimulus
line([0, 0], y_limits, 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '--')

% highlight channel if required
if highlight
    P(end + 1) = plot(input.x, input.data(eoi, :), 'Color', [0.9216    0.1490    0.1490], 'LineWidth', 3);
end

% plot legend if required
if plot_legend 
    legend(P, labels, 'Location', legend_loc, 'fontsize', 14)
    legend('boxoff');
end

% axes
box off;
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.TickDir = 'out'; 
ax.XColor = [0.5020    0.5020    0.5020]; 
ax.YColor = [0.5020    0.5020    0.5020]; 

% set x limits 
if x_limits(1) == 0 && x_limits(2) == 0
    xlim([input.x(1), input.x(end)]) 
else
    xlim(x_limits)
end

% referse y axis if required
if reverse
    set(gca, 'YDir', 'reverse');
end

% other parameters
xlabel('time (ms)')
ylabel('amplitude (\muV)')
set(gca, 'FontSize', 14)
ylim(y_limits)
set(gca, 'Layer', 'Top')
end
function export_EEGLAB(lwdata, filename, subj)
% =========================================================================
% exports data from letswave to EEGLAB format
% =========================================================================  
% dataset
EEG.setname = filename;
EEG.filename = [];
EEG.filepath = [];
EEG.subject = subj; 
EEG.session = 1;
    
% time properties
EEG.nbchan = lwdata.header.datasize(2);
EEG.trials = lwdata.header.datasize(1);
EEG.pnts = lwdata.header.datasize(6);
EEG.srate = 1/lwdata.header.xstep;
EEG.times = lwdata.header.xstart + (0:EEG.pnts-1)*lwdata.header.xstep;
EEG.xmin = EEG.times(1);
EEG.xmax = EEG.times(end);
EEG.data = permute(single(lwdata.data),[2,6,1,3,4,5]);
EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'SEEG_enabled');
EEG.chanlocs = rmfield(lwdata.header.chanlocs, 'topo_enabled');
    
% create events with appropriate latencies
EEG.event = lwdata.header.events;
if ~isempty(EEG.event)
    [EEG.event.type] = EEG.event.code;
    for e = 1:length(EEG.event)
        EEG.event(e).latency = (e-1)*EEG.pnts + 2001;
    end
    EEG.event = rmfield(EEG.event,'code');
end
    
% global tags
EEG.global_tags = lwdata.header.global_tags;
    
% create required empty fields
EEG.icawinv = [];
EEG.icaweights = [];
EEG.icasphere = [];
EEG.icaweights = [];
EEG.icaweights = [];
EEG.icaweights = [];
EEG.icaweights = [];
save([filename,'.set'], 'EEG');
end
function [lwdata] = export_lw(EEG, header, name)
% =========================================================================
% exports data from EEGLAB to letswave format
% ========================================================================= 
% transform data
data = [];
for t = 1:size(EEG.data, 3)
    for e = 1:size(EEG.data, 1)
        for i = 1:size(EEG.data, 2)
            data(t, e, 1, 1, 1, i) = EEG.data(e, i, t);
        end
    end
end
lwdata.data = data;    
    
% modify header
lwdata.header = header; 
lwdata.header.name = name;
lwdata.header.datasize = size(lwdata.data);
lwdata.header.chanlocs = lwdata.header.chanlocs(1:size(lwdata.data, 2));
lwdata.header.events = lwdata.header.events(1:size(lwdata.data, 1));
end
function [EEG_out] = SSP_SIR_trials(EEG_in, L, art_topographies, filt_ker, M)
% =========================================================================
% applies SSP-SIR using existing filters
% ========================================================================= 
%     L = param(c).L_ave;
%     EEG_in = EEG;
%     art_topographies = rmfield(param, {'PC', 'filter', 'L_ave'});
%     filt_ker = param(c).filter;
%     M = [];
%     EEG_out = EEG_in;
%     
%     % prepare individual filters
%     for t = 1:length(art_topographies)
%         statement = sprintf('P%d = eye(size(EEG_in.data,1)) - art_topographies(t).topographies*art_topographies(t).topographies'';', t);
%         eval(statement)
%     end
%     
%     % create final composite filter
%     statement = 'P = P1';
%     if length(art_topographies) > 1
%         for t = 2:length(art_topographies)
%             statement = [statement sprintf(' * P%d', t)];
%         end
%     end
%     statement = [statement ';'];
%     eval(statement);

EEG_out = EEG_in;
P = eye(size(EEG_in.data,1)) - art_topographies*art_topographies';

% apply filter
for i = 1:size(EEG_in.data,3)

    data = EEG_in.data(:,:,i);

    % suppressing the artifacts:
    data_clean = P*data;

    % performing SIR for the suppressed data:
    PL = P*L;

    if isempty (M)
        M = rank(data_clean) -  1 ;
    end

    tau_proj = PL*PL';
    [U,S,V] = svd(tau_proj);
    S_inv = zeros(size(S));
    S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
    tau_inv = V*S_inv*U';
    suppr_data_SIR = L*(PL)'*tau_inv*data_clean;

    % performing SIR for the original data:
    tau_proj = L*L';
    [U,S,V] = svd(tau_proj);
    S_inv = zeros(size(S));
    S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
    tau_inv = V*S_inv*U';
    orig_data_SIR = L*(L)'*tau_inv*data;

    if isempty(filt_ker)
        data_correct = suppr_data_SIR;
    else
        filt_ker_B = repmat(filt_ker,[size(suppr_data_SIR,1),1]);
        data_correct = filt_ker_B.*suppr_data_SIR + orig_data_SIR - filt_ker_B.*orig_data_SIR;
        filt_ker_B = [];
    end

    EEG_out.data(:,:,i) = data_correct;
end
end
