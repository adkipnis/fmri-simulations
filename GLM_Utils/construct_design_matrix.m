function [design, Opts] = construct_design_matrix(Dirs, Opts, n)
    % - For sanity check: Create a design matrix that does not differentiate between pictures
    % design = struct('name', {'Picture'}, 'onset', {events.onset}, 'duration', {events.duration}, 'tmod', {0}, 'pmod', {''});

    % - Design matrix: Events (Check if design matrix file exists and create it if not)
        design_mat = fullfile(Dirs.subdir,'spm_design.mat');
        events = spm_load(Dirs.event_files{n});
        % events = spm_load(event_files{n}(1:end-1));

        % extract unique event names and remove missing rows
        stim_list = string(num2str(events.stim_id, '%10.6f'));
        unique_ids = unique(stim_list); 
        Opts.stim_id = unique_ids(unique_ids ~= "           NaN");

%                 % Alternatively: extract unique event names and rename missing entries to 'Baseline'
%                 stim_id = string(events.stim_id);
%                 stim_id(ismissing(stim_id(:,1))) = 'Baseline';
%                 events.stim_id = stim_id;
%                 stim_id = unique(string(events.stim_id));

        design = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
        for d = 1:length(Opts.stim_id)
            id = find(Opts.stim_id(d) == stim_list); % find rows that contain this stimulus
            design(d).name = char(Opts.stim_id(d));
            design(d).onset = events.onset(id); 
            design(d).duration = events.duration(id);
            design(d).tmod = 0; % time modulation
            design(d).pmod = {''}; % time modulation
            design(d).orth = 1;   
        end

        % Export design_mat and make filename for stim_id_file (we will use it later as a dictionary for our beta coefficients)
        save(design_mat,'design');
        
        % Collect design matrix metrics
        Opts.n_stim = length(Opts.stim_id);
        Opts.n_stim_betas = Opts.n_stim*(1 + sum(Opts.hrf_derivs));
        Opts.n_conf = length(Opts.confound_names);
        Opts.n_reg = Opts.n_stim_betas + Opts.n_conf + 1;
        
end