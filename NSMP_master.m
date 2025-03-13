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
%                           4 - preparatory cue (= ball appears)
%                           8 - imperative cue (= bridge appears)
%           2) raw EMG recordings
%           --> triggers:   1 - Start of Trial
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
triggers = split(info.recording{4}, ',');
for t = 1:length(triggers)
    INFO(subject_idx).EEG.triggers(t).trigger = str2double(regexp(triggers{t}, '\d+\.?\d*', 'match', 'once'));
    pattern = sprintf('- ([a-zA-Z]+)');
    INFO(subject_idx).EEG.triggers(t).label = regexp(triggers{t}, pattern, 'tokens', 'once');
end
INFO(subject_idx).EEG.blocks = str2num(info.recording{5});

% save and continue
save(output_file, 'INFO','-append')
fprintf('done.\n\n')
clear output_vars info prompt dlgtitle dims definput pattern t triggers

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
params.run_analysis = true;
% ------------------------- 
fprintf('section 1: import & pre-process continuous data\n')

% encode session and subject info if necessary
user_input = input('do you want to encode the session info now? (y/n): ', 's');
if lower(user_input) == 'y'
    % input subject & session info
    prompt = {'subject', 'age:', 'sex:', 'handedness:', 'rMT:', 'session date:', 'session start:'};
    dlgtitle = 'subject & session information';
    dims = [1 35];
    definput = {'S000', '00', 'female/male', 'right/left', '00', sprintf('%s', date), '00:00'};
    info.session = inputdlg(prompt,dlgtitle,dims,definput);

    % identify subject index
    subject_idx = str2double(regexp(info.session{1}, '\d+', 'match', 'once'));

    % update info structure
    fprintf('encoding subject & session info...\n')
    INFO(subject_idx).ID = info.session{1};
    INFO(subject_idx).age = str2num(info.session{2});
    INFO(subject_idx).sex = info.session{3};
    INFO(subject_idx).handedness = info.session{4};
    INFO(subject_idx).session.date = info.session{6};
    INFO(subject_idx).session.start = info.session{7};
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
    triggers = split(info.recording{4}, ',');
    for t = 1:length(triggers)
        INFO(subject_idx).EEG.triggers(t).trigger = str2double(regexp(triggers{t}, '\d+\.?\d*', 'match', 'once'));
        pattern = sprintf('- ([a-zA-Z]+)');
        INFO(subject_idx).EEG.triggers(t).label = regexp(triggers{t}, pattern, 'tokens', 'once');
    end
    INFO(subject_idx).EEG.blocks = str2num(info.recording{5});

    % save and continue
    save(output_file, 'INFO','-append')
    fprintf('done.\n\n')
end

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
fprintf('subject %d (%s): ', subject_idx, INFO(subject_idx).ID)
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
    file_idx = false(1, length(data2import));
    for b = 1:length(data2import)
        for c = 1:length(params.blocks)
            if data2import(b).label == params.blocks(c)
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
    [dataset.raw(d).header, dataset.raw(d).data, ~] = RLW_import_MEGA(data2import(d).folder, data2import(d).label);

    % rename in the header
    dataset.raw(d).header.name = params.name;
end  
fprintf('done.\n')

% check trigger and event nembers, ask for continuation
for a = 1:length(dataset.raw)
    % check total number of triggers
    fprintf('block %d:\n%d triggers found\n', a, length(dataset.raw(a).header.events));
    
    % check trigger labels
    triggers = unique([dataset.raw(a).header.events.code]');
    fprintf('triggers are labeled: ');
    for b = 1:length(triggers)
        fprintf('%s ', triggers(b))
    end
    fprintf('\n')

    % check total number of events
    if ismember('1', triggers)
        events = 0;
        for b = 1:length(dataset.raw(a).header.events)
            if str2double(dataset.raw(a).header.events(b).code) == 1
                events = events + 1;
            end
        end
    else
        error('ERROR: no trigger is labeled ''1''! Please check manually.')
    end
    fprintf('in total %d events found\n\n', events);
end
user_input = input('do you want to continue with the pre-processing? (y/n): ', 's');
if lower(user_input) == 'n'
    fprintf('script terminated by user.\n');
    return; 
end

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% pre-process continuous data and save for letswave
fprintf('pre-processing:\n')
for d = 1:length(dataset.raw)
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
    fprintf('removing DC and applying linear detrend...\n')
    option = struct('linear_detrend', 1, 'suffix', params.suffix{3}, 'is_save', 1);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
    if d == 1
        INFO(subject_idx).EEG.processing(4).process = sprintf('DC + linear detrend on continuous data');
        INFO(subject_idx).EEG.processing(4).suffix = params.suffix{3};
        INFO(subject_idx).EEG.processing(4).date = sprintf('%s', date);
    end
    fprintf('\n')

    % update dataset
    dataset.raw(d).header = lwdata.header;
    dataset.raw(d).data = lwdata.data; 
end
fprintf('done.\n')

% run event analysis if required
if params.run_analysis
    % initialize logfile
    fprintf('analysing event delays...\n')
    filename = sprintf('%s\\logfiles\\event_analysis_%s.txt', folder.output, INFO(subject_idx).ID); 
    initialize_logfile(INFO, subject_idx, filename);

    % loop through datasets
    for a = 1:length(dataset.raw)
        % launch output
        info = struct;        
        info.block = str2double(regexp(dataset.raw(a).header.name, 'b(.*)', 'tokens', 'once'));
    
        % mark trial starts
        idx_start = false(1, length(dataset.raw(a).header.events));
        for b = 1:length(dataset.raw(a).header.events)-2
            if strcmp(dataset.raw(a).header.events(b).code, 'start')
                idx_start(b) = true;
            end
        end
        trigger_start = find(idx_start);
        info.trials = sum(idx_start);
    
        % select baseline and delay trials
        idx_bl = false(1, length(dataset.raw(a).header.events));
        idx_delay = false(1, length(dataset.raw(a).header.events));
        for b = 1:length(trigger_start)
            if b ~= length(trigger_start)
                trigger_end = trigger_start(b+1) - 1;
            else
                trigger_end = length(dataset.raw(a).header.events);
            end
            if strcmp(dataset.raw(a).header.events(trigger_start(b) + 1).code, 'stimulation') ...
                    && strcmp(dataset.raw(a).header.events(trigger_start(b) + 2).code, 'preparatory') ...
                    && strcmp(dataset.raw(a).header.events(trigger_start(b) + 3).code, 'imperative')
                idx_bl(trigger_start(b):trigger_end) = true;
            elseif strcmp(dataset.raw(a).header.events(trigger_start(b) + 1).code, 'preparatory') ...
                    && strcmp(dataset.raw(a).header.events(trigger_start(b) + 2).code, 'stimulation') ...
                    && strcmp(dataset.raw(a).header.events(trigger_start(b) + 3).code, 'imperative')
                idx_delay(trigger_start(b):trigger_end) = true;
            end
        end
    
        % exract info about baseline trials
        events.baseline = dataset.raw(a).header.events(idx_bl);
        idx_bl_start = false(1, length(events.baseline));
        delay_bl.stimulus = [];
        delay_bl.preparatory = [];
        delay_bl.imperative = [];
        for c = 1:length(events.baseline)
            if strcmp(events.baseline(c).code, 'start')
                idx_bl_start(c) = true;
                delay_bl.stimulus(end + 1) = events.baseline(c + 1).latency - events.baseline(c).latency;
                delay_bl.preparatory(end + 1) = events.baseline(c + 2).latency - events.baseline(c + 1).latency;
                if strcmp(events.baseline(c + 3).code, 'imperative')
                    delay_bl.imperative(end + 1) = events.baseline(c + 3).latency - events.baseline(c + 2).latency;
                end
            end
        end
        info.baseline.trials = sum(idx_bl_start);
        info.baseline.delay.stimulus.mean = mean(delay_bl.stimulus);
        info.baseline.delay.stimulus.min = min(delay_bl.stimulus);
        info.baseline.delay.stimulus.max = max(delay_bl.stimulus);
        info.baseline.delay.preparatory.mean = mean(delay_bl.preparatory);
        info.baseline.delay.preparatory.min = min(delay_bl.preparatory);
        info.baseline.delay.preparatory.max = max(delay_bl.preparatory);
        info.baseline.delay.imperative.mean = mean(delay_bl.imperative);
        info.baseline.delay.imperative.min = min(delay_bl.imperative);
        info.baseline.delay.imperative.max = max(delay_bl.imperative);
    
        % exract info about delay trials
        events.delay = dataset.raw(a).header.events(idx_delay);
        idx_delay_start = false(1, length(events.delay));
        delay_delay.stimulus = [];
        delay_delay.preparatory = [];
        delay_delay.imperative = [];
        for c = 1:length(events.delay)
            if strcmp(events.delay(c).code, 'start')
                idx_delay_start(c) = true;
                delay_delay.preparatory(end + 1) = events.delay(c + 1).latency - events.delay(c).latency;
                delay_delay.stimulus(end + 1) = events.delay(c + 2).latency - events.delay(c + 1).latency;
                if strcmp(events.delay(c + 3).code, 'imperative')
                    delay_delay.imperative(end + 1) = events.delay(c + 3).latency - events.delay(c + 2).latency;
                end
            end
        end
        info.delay.trials = sum(idx_delay_start);
        info.delay.delay.preparatory.mean = mean(delay_delay.preparatory);
        info.delay.delay.preparatory.min = min(delay_delay.preparatory);
        info.delay.delay.preparatory.max = max(delay_delay.preparatory);
        info.delay.delay.stimulus.mean = mean(delay_delay.stimulus);
        info.delay.delay.stimulus.min = min(delay_delay.stimulus);
        info.delay.delay.stimulus.max = max(delay_delay.stimulus);
        info.delay.delay.imperative.mean = mean(delay_delay.imperative);
        info.delay.delay.imperative.min = min(delay_delay.imperative);
        info.delay.delay.imperative.max = max(delay_delay.imperative);
    
        % write to the text file
        logfile_entry(info, filename)
    end
    fprintf('done.\n')
end

% save and continue
save(output_file, 'INFO','-append')
clear a b c d e session_folders session_date data2import data_idx file_idx  prompt dlgtitle dims definput ...
    triggers events user_input event_idx event_count info option lwdata user_input input filename ...
    idx_bl idx_delay idx_start info trigger_start events idx_bl_start idx_delay_start delay_bl delay_delay trigger_end 
fprintf('section 1 finished.\n\n')

%% 2) pre-process TEPs 
% ----- section input -----
params.prefix = 'dc ds crop';
params.conditions = {'baseline' 'delay'};
params.suffix = {'ep' 'dc' 'art_interp' 'preprocessed'};
params.eventcode = 'stimulation';
params.epoch = [-1.5 1.5];
params.artifact_interp = [-0.005 0.01];
params.artifact_method = 'pchip';
% ------------------------- 
fprintf('section 2: TEP pre-processing\n')

% load dataset if needed
if exist('dataset') ~= 1
    fprintf('loading dataset...\n')
    data2load = dir(sprintf('%s*%s*', params.prefix, INFO(subject_idx).ID));
    dataset = reload_dataset(data2load, [], 'raw');
    fprintf('done.\n')
end

% pre-process 
for d = 1:length(dataset.raw)
    % provide update
    fprintf('***** recording block %d *****\n', str2double(regexp(dataset.raw(d).header.name, 'b(.*)', 'tokens', 'once')))

    % mark events of different conditions
    idx_bl = false(1, length(dataset.raw(d).header.events));
    idx_delay = false(1, length(dataset.raw(d).header.events));
    for e = 1:length(dataset.raw(d).header.events) 
        if strcmp(dataset.raw(d).header.events(e).code, 'start') ...
                && strcmp(dataset.raw(d).header.events(e + 1).code, 'stimulation') ...
                && strcmp(dataset.raw(d).header.events(e + 2).code, 'preparatory') ...
                && strcmp(dataset.raw(d).header.events(e + 3).code, 'imperative')
            idx_bl(e:e+3) = true;
        elseif strcmp(dataset.raw(d).header.events(e).code, 'start') ...
                && strcmp(dataset.raw(d).header.events(e + 1).code, 'preparatory') ...
                && strcmp(dataset.raw(d).header.events(e + 2).code, 'stimulation') ...
                && strcmp(dataset.raw(d).header.events(e + 3).code, 'imperative')
            idx_delay(e:e+3) = true;
        end
    end

    % loop through conditions
    for c = 1:length(params.conditions)
        % provide update
        fprintf('%s condition:\n', params.conditions{c})

        % add letswave 7 to the top of search path
        addpath(genpath([folder.toolbox '\letswave 7']));

        % subset data
        lwdata.header = dataset.raw(d).header;
        lwdata.data = dataset.raw(d).data; 
        if c == 1
            lwdata.header.events = lwdata.header.events(idx_bl);
        elseif c == 2
            lwdata.header.events = lwdata.header.events(idx_delay);
        end

        % segment
        fprintf('segmenting ...\n')
        option = struct('event_labels', params.eventcode, 'x_start', params.epoch(1), 'x_end', params.epoch(2), ...
            'x_duration', params.epoch(2)-params.epoch(1), 'suffix', params.suffix{1}, 'is_save', 0);
        lwdata = FLW_segmentation.get_lwdata(lwdata, option);
        if d == 1 && c == 1
            INFO(subject_idx).EEG.processing(5).process = 'segmented to epochs';
            INFO(subject_idx).EEG.processing(5).params.limits = params.epoch;
            INFO(subject_idx).EEG.processing(5).suffix = params.suffix{1};
            INFO(subject_idx).EEG.processing(5).date = sprintf('%s', date);
        end

        % add letswave 6 to the top of search path
        addpath(genpath([folder.toolbox '\letswave 6']));

        % interpolate TMS artifact
        fprintf('interpolating TMS artifact...\n')
        header = lwdata.header;
        data = lwdata.data;
        [header, data, ~] = RLW_suppress_artifact_event(header, data, ...
            'xstart', params.artifact_interp(1), 'xend', params.artifact_interp(2), ...
            'event_code', params.eventcode, 'interp_method', params.artifact_method); 
        header.name = [params.suffix{3} ' ' header.name];
        if d == 1 && c == 1
            INFO(subject_idx).EEG.processing(6).process = 'TMS artifact interpolated';
            INFO(subject_idx).EEG.processing(6).params.limits = params.artifact_interp;
            INFO(subject_idx).EEG.processing(6).params.method = params.artifact_method;
            INFO(subject_idx).EEG.processing(6).suffix = params.suffix{3};
            INFO(subject_idx).EEG.processing(6).date = sprintf('%s', date);
        end

        % remove DC + linear detrend
        fprintf('DC + linear detrend ...\n')
        [header, data, ~]=RLW_dc_removal(header, data, 'linear_detrend', 1);
        if d == 1 && c == 1
            INFO(subject_idx).EEG.processing(7).process = 'DC + linear detrend on epoched data';
            INFO(subject_idx).EEG.processing(7).suffix = params.suffix{2};
            INFO(subject_idx).EEG.processing(7).date = sprintf('%s', date);
        end

        % usave to temporary structure
        blocks_preprocessed((d-1)*length(params.conditions) + c).condition = params.conditions{c}; 
        blocks_preprocessed((d-1)*length(params.conditions) + c).header = header;
        blocks_preprocessed((d-1)*length(params.conditions) + c).data = data; 
    end 
end
fprintf('done.\n\n')

% concatenate by condition
fprintf('concatenating epochs by condition, saving...\n')
for c = 1:length(params.conditions)
    % identify datasets to concatenate
    idx = false(1, length(blocks_preprocessed));
    for b = 1:length(blocks_preprocessed)
        if strcmp(blocks_preprocessed(b).condition, params.conditions{c})
            idx(b) = true;
        end
    end

    % concatenate epochs
    for b = 1:length(blocks_preprocessed)
        if idx(b)
            if ~exist('subset')
                subset.header = blocks_preprocessed(b).header;
                subset.data = blocks_preprocessed(b).data;
            else
                subset.data = cat(1, subset.data, blocks_preprocessed(b).data);
                subset.header.events(end + 1 : end + length(blocks_preprocessed(b).header.events)) = blocks_preprocessed(b).header.events;
            end
        end
    end

    % adjust header
    subset.header.name = sprintf('%s %s %s %s', params.suffix{4}, study, INFO(subject_idx).ID, params.conditions{c}); 
    subset.header.datasize = size(subset.data);
    triggers_n = length(subset.header.events) / subset.header.datasize(1);
    event_counter = 1;
    for e = 1:triggers_n:length(subset.header.events)
        for a = 1:triggers_n
            subset.header.events(e - 1 + a).epoch = event_counter;
        end
        event_counter = event_counter + 1;
    end

    % save for letswave
    data = subset.data;
    header = subset.header;
    save(sprintf('%s.lw6', header.name), 'header')
    save(sprintf('%s.mat', header.name), 'data') 
    clear subset

    % encode
    if c == 1
        INFO(subject_idx).EEG.processing(8).process = sprintf('epochs concatenated by condition');
        INFO(subject_idx).EEG.processing(8).date = sprintf('%s', date);
    end
    INFO(subject_idx).EEG.processing(8).params(c).condition = params.conditions{c};
    INFO(subject_idx).EEG.processing(8).params(c).epochs = header.datasize(1);
    INFO(subject_idx).EEG.processing(8).suffix{c} = params.conditions{c};    
end
fprintf('done.\n')

% open letswave for visual check
fig_all = findall(0, 'Type', 'figure');
open = true;
for f = 1:length(fig_all)
    if contains(get(fig_all(f), 'Name'), 'Letswave', 'IgnoreCase', true)
        open = false;
        break;
    end
end
if open
    letswave
end
fprintf('\n')

% save and continue
save(output_file, 'INFO','-append')
clear a b c d e f data2load lwdata option idx idx_bl idx_delay header data blocks_preprocessed ...
    subset event_counter triggers_n fig_all open
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
function initialize_logfile(INFO, subject_idx, filename)
fileID = fopen(filename, 'w');
fprintf(fileID, sprintf('********** EVENT ANALYSIS **********\r\n')); 
fprintf(fileID, sprintf('subject: %s\r\n', INFO(subject_idx).ID));
fprintf(fileID, sprintf('date of experiment: %s\r\n', INFO(subject_idx).session.date));
fprintf(fileID, '\r\n');
fclose(fileID);
end
function logfile_entry(info, filename)
fileID = fopen(filename, 'a');
fprintf(fileID, sprintf('*** recording block %d ***\r\n', info.block));
fprintf(fileID, sprintf('trials in total: %d\r\n', info.trials));
fprintf(fileID, '\r\n');
fprintf(fileID, sprintf('baseline trials %d:\r\n', info.baseline.trials));
fprintf(fileID, sprintf('     - mean delay from start to stimulus: %ds (min: %ds; max: %ds)\r\n', ...
    info.baseline.delay.stimulus.mean, info.baseline.delay.stimulus.min, info.baseline.delay.stimulus.max));
fprintf(fileID, sprintf('     - mean delay from stimulus to preparatory cue: %ds (min: %ds; max: %ds)\r\n', ...
    info.baseline.delay.preparatory.mean, info.baseline.delay.preparatory.min, info.baseline.delay.preparatory.max));
fprintf(fileID, sprintf('     - mean delay from preparatory cue to imperative cue: %ds (min: %ds; max: %ds)\r\n', ...
    info.baseline.delay.imperative.mean, info.baseline.delay.imperative.min, info.baseline.delay.imperative.max));
fprintf(fileID, '\r\n');
fprintf(fileID, sprintf('delay trials %d:\r\n', info.delay.trials));
fprintf(fileID, sprintf('     - mean delay from start to preparatory cue: %ds (min: %ds; max: %ds)\r\n', ...
    info.delay.delay.preparatory.mean, info.delay.delay.preparatory.min, info.delay.delay.preparatory.max));
fprintf(fileID, sprintf('     - mean delay from preparatory cue to stimulus: %ds (min: %ds; max: %ds)\r\n', ...
    info.delay.delay.stimulus.mean, info.delay.delay.stimulus.min, info.delay.delay.stimulus.max));
fprintf(fileID, sprintf('     - mean delay from stimulus to imperative cue: %ds (min: %ds; max: %ds)\r\n', ...
    info.delay.delay.imperative.mean, info.delay.delay.imperative.min, info.delay.delay.imperative.max));
fprintf(fileID, '\r\n');
fclose(fileID);
end