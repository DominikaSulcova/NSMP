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
%               1) import & pre-process continuous data
%                       - assign electrode coordinates
%                       - check and re-label events
%                       - crop continuous data and save for future analysis
%                       - downsample by indicated ratio (default 20)
%                       - remove DC + linear detrend
%                       - optionally runs the analysis of trigger latencies
%               --> pre-processed data saved in a structure 'dataset'
%               2) pre-process TEPs 
%                       - segmented relative to TMS events (default [-1.5 1.5]s)
%                       - remove DC + linear detrend on single epochs
%                       - interpolate TMS artifact (cubic interpolation, default [-5 10]ms)
%                       - concatenate blocks by condition, save for letswave
%               ==> at this point, data should be visually inspected to
%               identify bad channels
%               3) export to EEGLAB
%                       - interpolate channels if needed (default 6 neighbor channels)
%                       - save as .set
%               4) SSP-SIR - removal of the muscular artifact 
%                       - re-reference to common average
%                       - apply SSP-SIR on merged data (spherical model) 
%                       - apply frequency filters - Hamming windowed sinc FIR filter
%                       (default bandpass [0.5 80]Hz, notch [48.5, 51.5]Hz)
%                       - save for letswave, plot output of filtering
%               ==> at this point, data should be visually inspected and bad trials
%               manually discarded (implemented in letswave 6)
%               5) ICA - removal of remaining artifacts 
%                       - load dataset with bad trials removed, encode 
%                       - compute ICA matrix, save for letswave
%                       - plot IC topographies and spectral content
%               ==> at this point, ICA should be conducted and artifactual ICs
%               manually discarded based on topography, timecourse, spectral content 
%               and MARA labelling (implemented in letswave 6)
%                       - encode and classify discarded ICs
%                       - extract and plot mean signal at electrode of interest
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
definput = {'20000', 'FCz', 'AFz', '1 - start, 2 - stimulation, 4 - preparatory, 8 - imperative', '1,2,3,4,5'};
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
    definput = {'20000', 'Fz', 'AFz', '1 - start, 2 - stimulation, 4 - preparatory, 8 - imperative', '1,2,3,4,5'};
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

    % save to the dataset
    dataset.preprocessed(c).condition = params.conditions{c};
    dataset.preprocessed(c).header = subset.header;
    dataset.preprocessed(c).data = subset.data;

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

%% 3) interpolate bad chanels, export to EEGLAB
% ----- section input -----
params.prefix = 'preprocessed';
params.conditions = {'baseline' 'delay'};
params.interp_chans = 6;
% -------------------------
fprintf('section 3: exporting to EEGLAB\n')

% load dataset if needed
if exist('dataset') ~= 1
    fprintf('loading dataset...\n')
    data2load = dir(sprintf('%s*%s*', params.prefix, INFO(subject_idx).ID));
    dataset = reload_dataset(data2load, [], 'preprocessed');
    fprintf('done.\n')
end

% interpolate channels if needed
params.labels = {dataset.preprocessed(1).header.chanlocs.labels};
prompt = {'channel(s) to interpolate:'};
definput = {strjoin(params.labels, ' ')};
dlgtitle = 'channel interpolation';
dims = [1 220];
answer = inputdlg(prompt,dlgtitle,dims,definput);
for a = 1:length(answer)
    if ~isempty(answer{a})
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
                chan_dist = -ones(length(dataset.preprocessed(1).header.chanlocs), 1);
                for b = setdiff(1:chan_n, chan_n)
                    if dataset.preprocessed(1).header.chanlocs(b).topo_enabled == 1
                        chan_dist(b) = sqrt((dataset.preprocessed(1).header.chanlocs(b).X - dataset.preprocessed(1).header.chanlocs(chan_n).X)^2 + ...
                            (dataset.preprocessed(1).header.chanlocs(b).Y - dataset.preprocessed(1).header.chanlocs(chan_n).Y)^2 + ...
                            (dataset.preprocessed(1).header.chanlocs(b).Z - dataset.preprocessed(1).header.chanlocs(chan_n).Z)^2);
                    end
                end
                chan_dist((chan_dist==-1)) = max(chan_dist);
                [~,chan_idx] = sort(chan_dist);
    
                % identify neighbouring channels
                chan_idx = chan_idx(1:params.interp_chans);
                chans2use = params.labels;
                chans2use = chans2use(chan_idx);
    
                % cycle through all datasets
                for d = 1:length(dataset.preprocessed)
                    % select data
                    lwdata.header = dataset.preprocessed(d).header;
                    lwdata.data = dataset.preprocessed(d).data;
        
                    % interpolate using the neighboring electrodes
                    option = struct('channel_to_interpolate', chans2interpolate{c}, 'channels_for_interpolation_list', {chans2use}, ...
                        'suffix', '', 'is_save', 0);
                    lwdata = FLW_interpolate_channel.get_lwdata(lwdata, option);
        
                    % update dataset
                    dataset.preprocessed(d).header = lwdata.header;
                    dataset.preprocessed(d).data = lwdata.data;  
                end
                
                % encode
                if c == 1
                    INFO(subject_idx).EEG.processing(9).process = sprintf('bad channels interpolated');
                    INFO(subject_idx).EEG.processing(9).date = sprintf('%s', date);
                end
                INFO(subject_idx).EEG.processing(9).params.bad{c} = chans2interpolate{c};
                INFO(subject_idx).EEG.processing(9).params.chans_used{c} = strjoin(chans2use, ' ');  
            end
        end
    else
        INFO(subject_idx).EEG.processing(9).process = sprintf('no channels interpolated');
        INFO(subject_idx).EEG.processing(9).date = sprintf('%s', date);
    end
end

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% save as .set
fprintf('saving for EEGLAB... ')
for d = 1:length(dataset.preprocessed)
    fprintf('. ')

    % select data
    lwdata.header = dataset.preprocessed(d).header;
    lwdata.data = dataset.preprocessed(d).data;

    % export 
    export_EEGLAB(lwdata, lwdata.header.name, INFO(subject_idx).ID);
end
fprintf('done.\n\n')

% save and continue
save(output_file, 'INFO','-append')
clear a b c d e prompt dlgtitle dims definput answer chans2interpolate chan_n chan_dist chan_idx chans2use lwdata option 
fprintf('section 3 finished.\n\n')

%% 4) SSP-SIR
% ----- section input -----
params.prefix = 'preprocessed';
params.suffix = {'preprocessed' 'sspsir' 'ffilt'};
params.conditions = {'baseline' 'delay'};
params.baseline = [-0.3 -0.006]; 
params.bandpass = [0.5, 80];
params.notch = [48.5, 51.5];
params.plot_output = true;
params.plot_toi = [-0.1 0.5];
params.plot_eoi = 'C3';
% -------------------------
fprintf('section 4: SSP-SIR\n')

% add eeglab to the top of search path and launch
addpath(fullfile(folder.toolbox, 'EEGLAB'));
eeglab

% load all datasets, re-reference 
fprintf('loading dataset: ')
for a = 1:length(dataset.preprocessed)
    fprintf('. ')
    % load dataset, re-reference  
    name = sprintf('%s %s %s %s.set', params.suffix{1}, study, INFO(subject_idx).ID, params.conditions{a}); 
    EEG = pop_loadset('filename', name, 'filepath', folder.processed);
    EEG.filename = name;
    EEG.filepath = folder.processed;
    EEG = pop_reref(EEG, []);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, a);
    eeglab redraw  
    if a == 1 
        INFO(subject_idx).EEG.processing(10).process = 're-referenced to common average';
        INFO(subject_idx).EEG.processing(10).date = sprintf('%s', date);
    end
end
fprintf('done.\n')

% merge the data
fprintf('merging ...\n')
merged_EEG = pop_mergeset(ALLEEG, 1:length(ALLEEG), 0);

% apply SSP-SIR - spherical model 
fprintf('applying SSP-SIR ...\n')
merged_EEG = pop_tesa_sspsir(merged_EEG, 'artScale', 'automatic', 'PC', []);
prompt = {'number of rejected PCs:'};
dlgtitle = 'SSP-SIR';
dims = [1 40];
definput = {''};
input = inputdlg(prompt,dlgtitle,dims,definput);

% encode 
INFO(subject_idx).EEG.processing(11).process = 'muscular artifact removed using SSP-SIR';
INFO(subject_idx).EEG.processing(11).params.method = 'SSP-SIR';
INFO(subject_idx).EEG.processing(11).params.leadfield = 'default spherical';
INFO(subject_idx).EEG.processing(11).params.PC_removed = str2num(input{1,1});
INFO(subject_idx).EEG.processing(11).suffix = params.suffix{2};
INFO(subject_idx).EEG.processing(11).date = sprintf('%s', date);

% split back to original datasets
idx_start = 1;
for b = 1:length(dataset.preprocessed)
    % identify number of epochs 
    n_epochs = ALLEEG(b).trials;

    % extract epochs and update ALLEEG
    EEG = pop_select(merged_EEG, 'trial', idx_start:(idx_start + n_epochs - 1));

    % save new dataset
    name = sprintf('%s %s %s %s', params.suffix{2}, study, INFO(subject_idx).ID, params.conditions{b}); 
    EEG.setname = name;
    EEG.filename = sprintf('%s.set', name);
    pop_saveset(EEG, 'filename', name, 'filepath', folder.processed);

    % update starting index 
    idx_start = idx_start + n_epochs;       
end
    
% apply frequency filters and save
fprintf('applying frequency filters: \n')
for b = 1:length(dataset.preprocessed)
    fprintf('. ')

    % load dataset
    name = sprintf('%s %s %s %s.set', params.suffix{2}, study, INFO(subject_idx).ID, params.conditions{b}); 
    EEG = pop_loadset('filename', name, 'filepath', folder.processed);
    [ALLEEG, EEG, b] = eeg_store(ALLEEG, EEG, b);
    eeglab redraw 

    % bandpass filter
    EEG = pop_eegfiltnew(EEG, 'locutoff', params.bandpass(1), 'hicutoff', params.bandpass(2));

    % notch filter
    EEG = pop_eegfiltnew(EEG, 'locutoff', params.notch(1), 'hicutoff', params.notch(2), 'revfilt', 1);

    % % check spectrum
    % pop_spectopo(EEG, 1, [], 'EEG', 'freqrange', [0 100]);

    % encode to info structure
    if b == 1
        INFO(subject_idx).EEG.processing(12).process = sprintf('frequency filters applied');
        INFO(subject_idx).EEG.processing(12).params.bandpass = params.bandpass;
        INFO(subject_idx).EEG.processing(12).params.notch = params.notch;
        INFO(subject_idx).EEG.processing(12).suffix = params.suffix{3};
        INFO(subject_idx).EEG.processing(12).date = sprintf('%s', date);
    end

    % update dataset and save for letswave
    name = sprintf('%s %s %s %s %s', params.suffix{3}, params.suffix{2},...
        study, INFO(subject_idx).ID, params.conditions{b}); 
    lwdata = export_lw(EEG, dataset.preprocessed(b).header, name);
    dataset.sspsir(b).header = lwdata.header;
    dataset.sspsir(b).data = lwdata.data;
    header = lwdata.header;
    data = lwdata.data;
    save([name '.lw6'], 'header')
    save([name '.mat'], 'data')
    fprintf('\n')
end

% update figure counter
fig_all = findall(0, 'Type', 'figure');
figure_counter = length(fig_all) + 1;

% plot the output figure
if params.plot_output
    % define common visual parameters
    params.labels = {dataset.sspsir(1).header.chanlocs.labels};
    visual.x = dataset.sspsir(1).header.xstart : dataset.sspsir(1).header.xstep : ...
        dataset.sspsir(1).header.xstep * dataset.sspsir(1).header.datasize(6) + dataset.sspsir(1).header.xstart - dataset.sspsir(1).header.xstep;
    visual.labels = {'original', 'filtered'};
    visual.xlim = params.plot_toi;
    visual.colors = [0.3333    0.4471    0.9020;
                    1.0000    0.0745    0.6510];
    visual.eoi  = find(strcmp(params.labels, params.plot_eoi));

    % launch the figure
    fig = figure(figure_counter);
    set(fig, 'units', 'normalized', 'outerposition', [0 0 1 1])
    hold on

    % plot data
    for a = 1:length(params.conditions)
        % define t value
        visual.t_value = tinv(0.975, size(dataset.sspsir(a).data, 1) - 1); 

        % select original data
        visual.data(1, :) = squeeze(mean(dataset.preprocessed(a).data(:, visual.eoi, 1, 1, 1, :), 1));  
        visual.sem(1, :) = squeeze(std(dataset.preprocessed(a).data(:, visual.eoi, 1, 1, 1, :), 0, 1)) / sqrt(size(dataset.preprocessed(a).data, 1)); 
        visual.CI_upper(1, :) = visual.data(1, :) + visual.t_value * visual.sem(1, :); 
        visual.CI_lower(1, :) = visual.data(1, :) - visual.t_value * visual.sem(1, :); 
            
        % select filtered data
        visual.data(2, :) = squeeze(mean(dataset.sspsir(a).data(:, visual.eoi, 1, 1, 1, :), 1));  
        visual.sem(2, :) = squeeze(std(dataset.sspsir(a).data(:, visual.eoi, 1, 1, 1, :), 0, 1)) / sqrt(size(dataset.sspsir(a).data, 1)); 
        visual.CI_upper(2, :) = visual.data(2, :) + visual.t_value * visual.sem(2, :); 
        visual.CI_lower(2, :) = visual.data(2, :) - visual.t_value * visual.sem(2, :); 
            
        % plot
        subplot(1, length(params.conditions), a)
        plot_ERP(visual, 'xlim', visual.xlim, 'colours', visual.colors, 'labels', visual.labels)
        title(sprintf('%s trials', params.conditions{a}))
    end

    % save and update figure counter
    saveas(fig, sprintf('%s\\figures\\SSPSIR_%s.png', folder.output, INFO(subject_idx).ID))
    figure_counter = figure_counter + 1;
end

% open letswave if not already open
addpath(genpath([folder.toolbox '\letswave 6']));
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
clear a b c f e fig_all name match prompt discarded answer data header lwdata data2load definput dims dlgtitle eoi fig idx_start input latency code ...
    n_epochs open tmpEEG tmpstr visual ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG merged_EEG globalvars LASTCOM PLUGINLIST STUDY GUI_handle GUI_OK
fprintf('section 4 finished.\nPlease check for bad trials now.\n\n')

%% 5) ICA
% ----- section input -----
params.prefix = 'ar ffilt sspsir';
params.conditions = {'baseline' 'delay'};
params.suffix = {'ar' 'ica' 'icfilt'};
params.ICA_comp = 40;
params.plot_toi = [-0.1 0.5];
params.plot_eoi = 'C3';
% -------------------------
fprintf('section 5: ICA\n')

% load dataset with bad trials removed
fprintf('updating dataset... ')
if exist('dataset') ~= 1
    % load previous dataset
    data2load = dir(sprintf('%s*%s*', params.prefix(4:end), INFO(subject_idx).ID));
    dataset = reload_dataset(data2load, [], 'sspsir');

    % append new dataset
    dataset_old = dataset;
    data2load = dir(sprintf('%s*%s*', params.prefix, INFO(subject_idx).ID));
    dataset = reload_dataset(data2load, [], 'checked');
    dataset.sspsir = dataset_old.sspsir;
else
    % append new dataset
    dataset_old = dataset;
    data2load = dir(sprintf('%s*%s*', params.prefix, INFO(subject_idx).ID));
    dataset = reload_dataset(data2load, [], 'checked');
    dataset.sspsir = dataset_old.sspsir;
end
fprintf('done.\n')

% encode bad trials
fprintf('encoding bad trials... ')
INFO(subject_idx).EEG.processing(13).process = sprintf('bad trials discarded');
INFO(subject_idx).EEG.processing(13).params.GUI = 'letswave';
INFO(subject_idx).EEG.processing(13).suffix = params.suffix{1};
INFO(subject_idx).EEG.processing(13).date = sprintf('%s', date);
for a = 1:length(params.conditions)
    % extract discarded expochs
    if ~isempty(dataset.checked(a).header.history(end).configuration)
        if ~isempty(dataset.checked(a).header.history(end).configuration.parameters.rejected_epochs)
            discarded = dataset.checked(a).header.history(end).configuration.parameters.rejected_epochs;
        end
    end

    % encode 
    INFO(subject_idx).EEG.processing(13).params.discarded{a} = discarded;
    INFO(subject_idx).EEG.processing(13).params.kept(a) = size(dataset.checked(a).data, 1);
end
fprintf('done.\n')

% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));
  
% compute ICA and save  
lwdataset = dataset.checked;
fprintf('computing ICA matrix:\n')
option = struct('ICA_mode', 2, 'algorithm', 1, 'num_ICs', params.ICA_comp, 'suffix', params.suffix{2}, 'is_save', 1);
lwdataset = FLW_compute_ICA_merged.get_lwdataset(lwdataset, option);
fprintf('done.\n')

% extract ICA parameters
matrix.mix = lwdataset(1).header.history(end).option.mix_matrix;
matrix.unmix = lwdataset(1).header.history(end).option.unmix_matrix;    
params.ICA_chanlocs = lwdataset(1).header.chanlocs;
for i = 1:size(matrix.mix, 2)
    params.ICA_labels{i} = ['IC',num2str(i)];
end
params.ICA_SR = 1/lwdataset(1).header.xstep;

% update dataset and adjust for letswave 6
dataset.ica = lwdataset;
for b = 1:length(params.conditions)
    dataset.ica(b).header.history(end).configuration.gui_info.function_name = 'LW_ICA_compute_merged';  
    dataset.ica(b).header.history(end).configuration.parameters = dataset.ica(b).header.history(end).option;  
    [dataset.ica(b).header.history(end).configuration.parameters.ICA_um] = dataset.ica(b).header.history(end).configuration.parameters.unmix_matrix; 
    [dataset.ica(b).header.history(end).configuration.parameters.ICA_mm] = dataset.ica(b).header.history(end).configuration.parameters.mix_matrix; 
    dataset.ica(b).header.history(end).configuration.parameters = rmfield(dataset.ica(b).header.history(end).configuration.parameters, {'unmix_matrix' 'mix_matrix'});
    header = dataset.ica(b).header;
    save(sprintf('%s.lw6', dataset.ica(b).header.name), 'header');
end

% unmix data
for b = 1:length(dataset.ica)
    for e = 1:size(dataset.ica(b).data, 1)
        dataset.unmixed(b).header = dataset.ica(b).header;
        dataset.unmixed(b).data(e, :, 1, 1, 1, :) = matrix.unmix * squeeze(dataset.ica(b).data(e, :, 1, 1, 1, :));        
    end
end

% update info structure
INFO(subject_idx).EEG.processing(14).process = 'ICA matrix computed';
INFO(subject_idx).EEG.processing(14).params.method = 'runica';
INFO(subject_idx).EEG.processing(14).params.components = params.ICA_comp;
INFO(subject_idx).EEG.processing(14).params.chanlocs = params.ICA_chanlocs;
INFO(subject_idx).EEG.processing(14).params.labels = params.ICA_labels;
INFO(subject_idx).EEG.processing(14).params.SR = params.ICA_SR;
INFO(subject_idx).EEG.processing(14).params.matrix = matrix;
INFO(subject_idx).EEG.processing(14).suffix = params.suffix{2};
INFO(subject_idx).EEG.processing(14).date = sprintf('%s', date);
 
% calculate PSD across all timepoints, components and trials 
fprintf('estimating spectral content of ICs...\n')
psd = [];
for b = 1:length(params.conditions)
    for c = 1:params.ICA_comp
        for e = 1:size(dataset.unmixed(b).data, 1)
            [psd(b, c, e, :), freq] = pwelch(squeeze(dataset.unmixed(b).data(e, c, 1, 1, 1, :)), ...
                [], [], [], INFO(subject_idx).EEG.processing(14).params.SR);  
        end
    end
end
INFO(subject_idx).EEG.processing(14).params.PSD = squeeze(mean(psd, [1, 3]));
fprintf('done.\n')

% plot IC spectral content (in 2 figures to keep the size readable)
for a = 1:2
    figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    components2plot = (a-1)*(params.ICA_comp/2) + [1:ceil(params.ICA_comp/2)];
    for f = 1:ceil(params.ICA_comp/2)
        % plot the topography
        matrix = INFO(subject_idx).EEG.processing(14).params.matrix.mix;
        subplot(ceil(length(components2plot)/3), 6, (f-1)*2 + 1);
        topoplot(double(matrix(:, components2plot(f))'), params.ICA_chanlocs, 'maplimits', [-4 4], 'shading', 'interp', 'whitebk', 'on', 'electrodes', 'off')
        set(gca,'color',[1 1 1]);
        title(params.ICA_labels{components2plot(f)})
    
        % plot the psd
        subplot(ceil(length(components2plot)/3), 6, (f-1)*2 + 2);
        plot(freq(1:max(find(freq <= 80))), log10(INFO(subject_idx).EEG.processing(14).params.PSD(components2plot(f), 1:max(find(freq <= 80)))));
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
    end
    saveas(gcf, sprintf('%s\\figures\\ICA_%s_%d.png', folder.output, INFO(subject_idx).ID, a));
end

% open letswave if not already open
fig_all = findall(0, 'Type', 'figure');
open = true;
for f = 1:length(fig_all)
    if contains(get(fig_all(f), 'Name'), 'Letswave', 'IgnoreCase', true)
        open = false;
        break;
    end
end
if open
    fprintf('opening letswave:\n')
    addpath(genpath([folder.toolbox '\letswave 6']));
    letswave
end

% identify outcome filenames and wait for datasets to appear
for a = 1:length(params.conditions)
    filenames{(a-1)*2 + 1} = sprintf('%s %s %s %s %s %s.lw6', params.suffix{3}, params.suffix{2}, params.prefix, study, INFO(subject_idx).ID, params.conditions{a});
    filenames{(a-1)*2 + 2} = sprintf('%s %s %s %s %s %s.mat', params.suffix{3}, params.suffix{2}, params.prefix, study, INFO(subject_idx).ID, params.conditions{a});
end
wait4files(filenames);
close(gcf)

% load filtered data
fprintf('loading filtered dataset...\n')
for a = 1:length(dataset.ica)
    data_name = sprintf('%s %s %s %s %s %s', params.suffix{3}, params.suffix{2}, params.prefix, study, INFO(subject_idx).ID, params.conditions{a});
    load(sprintf('%s.lw6', data_name), '-mat');
    dataset.ica_filtered(a).header = header; 
    load(sprintf('%s.mat', data_name));
    dataset.ica_filtered(a).data = data; 
end
fprintf('done.\n')

% ask for the input
prompt = {'blinks:', 'horizontal:', 'TMS:', 'muscles:', 'slow artifacts:'};
dlgtitle = 'ICA';  
dims = [1 60];
definput = {'', '', '', '', ''};
input = inputdlg(prompt,dlgtitle,dims,definput);

% encode ICA 
INFO(subject_idx).EEG.processing(15).process = 'artifactual ICs discarded';
INFO(subject_idx).EEG.processing(15).suffix = params.suffix{3};
INFO(subject_idx).EEG.processing(15).date = sprintf('%s', date);
INFO(subject_idx).EEG.processing(15).params.kept = params.ICA_comp - length([str2num(input{1}), str2num(input{2}), str2num(input{3}), str2num(input{4}), str2num(input{5})]);
INFO(subject_idx).EEG.processing(15).params.removed.blinks = str2num(input{1});
INFO(subject_idx).EEG.processing(15).params.removed.horizontal = str2num(input{2});
INFO(subject_idx).EEG.processing(15).params.removed.TMS = str2num(input{3});
INFO(subject_idx).EEG.processing(15).params.removed.muscles = str2num(input{4});
INFO(subject_idx).EEG.processing(15).params.removed.slow = str2num(input{5});
save(output_file, 'INFO','-append')

% launch the figure
fig_all = findall(0, 'Type', 'figure');
figure_counter = length(fig_all) + 1;
fig = figure(figure_counter);
screen_size = get(0, 'ScreenSize');
set(fig, 'Position', [screen_size(3)/4, screen_size(4)/4, screen_size(3) / 2.5, screen_size(4) / 2])

% define common visual parameters
params.labels = {dataset.ica(1).header.chanlocs.labels};
visual.x = dataset.ica(1).header.xstart : dataset.ica(1).header.xstep : ...
    dataset.ica(1).header.xstep * dataset.ica(1).header.datasize(6) + dataset.ica(1).header.xstart - dataset.ica(1).header.xstep;
visual.labels = params.conditions;
visual.xlim = params.plot_toi;
visual.colors = [0.9216    0.1490    0.1490;
                0    0.4471    0.7412];
visual.eoi  = find(strcmp(params.labels, params.plot_eoi));

% extract mean signal per condition
for a = 1:length(params.conditions)
    % subset data
    data2plot = squeeze(dataset.ica_filtered(a).data(:, visual.eoi, 1, 1, 1, :));

    % define t value
    visual.t_value = tinv(0.975, size(data2plot, 1) - 1); 

    % select original data
    visual.data(a, :) = mean(data2plot, 1)';  
    visual.sem(a, :) = std(data2plot, 0, 1)' / sqrt(size(data2plot, 1)); 
    visual.CI_upper(a, :) = visual.data(a, :) + visual.t_value * visual.sem(a, :); 
    visual.CI_lower(a, :) = visual.data(a, :) - visual.t_value * visual.sem(a, :);
end

% plot
plot_ERP(visual, 'xlim', visual.xlim, 'colours', visual.colors, 'labels', visual.labels)

% save and update counter
saveas(fig, sprintf('%s\\figures\\TEP_%s_%s.png', folder.output, params.plot_eoi, INFO(subject_idx).ID))
figure_counter = figure_counter + 1;

% ask for continuation
fprintf('section 5 finished.\n')
answer = questdlg('do you want to continue with next subject?', 'Continue?', 'YES', 'NO', 'YES'); 
if strcmp(answer, 'YES')
    subject_idx = subject_idx + 1;
    clear dataset
    fprintf('proceeding to subject %d.\n\n', subject_idx)
end
clear a b c d e f i data2load check discarded match prompt definput input input_old dims dlgtitle matrix fig_all open ...
    psd freq option lwdataset labels header dataset_old components2plot filenames data_name data2plot answer data fig screen_size visual

%% extract TEP measures
% - calculates average TEP GFP and identifies peaks
% - plots average TEP with peak topographies
% - ...
% ----- section input -----
% -------------------------
fprintf('section : extract TEP measures\n')
fprintf('section  finished.\n\n')

%% extract beta measures
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
        EEG.event(e).latency = (e-1)*EEG.pnts + EEG.xmin*(-1)*EEG.srate;
    end
    EEG.event = rmfield(EEG.event,'code');
end
    
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
function lwdata = export_lw(EEG, header, name)
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
%           legend --> 'on'(default)/'off'
%           labels --> cell array with labels for the legend  
%           legend_loc --> legend location (default 'southeast')
%           eoi --> label of a channel to be highlighted
%           reverse --> 'on'/'off'(default) - flips y axis
%           interpolated --> time window that was interpolated
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
interpolate = false;

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
        eoi_n = find(contains(input.chanlabels, eoi));
        highlight = true;
    end 

    % interpolated interval - default off
    j = find(strcmpi(varargin, 'interpolated'));
    if ~isempty(j)
        interpolate_toi = varargin{j + 1};
        interpolate = true;
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

% highlight channel if required
if highlight
    P(end + 1) = plot(input.x, input.data(eoi_n, :), 'Color', [0.9216    0.1490    0.1490], 'LineWidth', 4);
    hold on
end

% shade interpolated window if required
if interpolate
    interpolate_x = [interpolate_toi(1), interpolate_toi(2), interpolate_toi(2), interpolate_toi(1)];
    interpolate_y = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
    fill(interpolate_x, interpolate_y, [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% plot stimulus
line([0, 0], y_limits, 'Color', 'black', 'LineWidth', 2.5, 'LineStyle', '--')

% plot legend if required
if plot_legend 
    legend(P, labels, 'Location', legend_loc, 'fontsize', 14)
    legend('boxoff');
else
    legend('off')
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
xlabel('time (s)')
ylabel('amplitude (\muV)')
set(gca, 'FontSize', 14)
ylim(y_limits)
set(gca, 'Layer', 'Top')
end
function wait4files(filenames)
% =========================================================================
% waitForFiles pauses the script until all required files appear in the
% working working directory
% --> file names are specified in a cell array
% =========================================================================    
% loop to wait for files
while true
    allFilesExist = true;
    
    % check for each file
    for i = 1:length(filenames)
        if isempty(dir(filenames{i}))
            allFilesExist = false;
            break;
        end
    end
    
    % if all files exist, break the loop
    if allFilesExist
        break;
    end
    
    % pause for 2s
    pause(2);
end
end