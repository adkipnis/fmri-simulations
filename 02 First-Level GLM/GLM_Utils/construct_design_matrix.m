function [Opts, Dirs] = construct_design_matrix(Dirs, Opts, n)

    % - Design matrix: Events (Check if design matrix file exists and create it if not)
        design_cond = fullfile(Dirs.results_dir,'spm_design_cond.mat');
        design_multi = fullfile(Dirs.results_dir,'spm_design_multi.mat');
        events = spm_load(Dirs.event_files{n});
        % events = spm_load(event_files{n}(1:end-1));

        % extract unique event names and remove missing rows
        stim_list = string(num2str(events.stim_id, '%10.6f'));
        unique_ids_tmp = unique(stim_list); 
        unique_ids = unique_ids_tmp(unique_ids_tmp ~= "           NaN");
        assert(length(unique_ids) == length(Opts.unique_conditions_to_include))
        
        % Save stim_ids and make copies for HRF derivatives
        Opts.stim_id = unique_ids(Opts.unique_conditions_to_include);
        Opts.stim_id_ext = Opts.stim_id; % for HRF_derivatives
        if sum(Opts.hrf_derivs) > 0
            for i = 1:sum(Opts.hrf_derivs)
                tmp = strcat(Opts.stim_id, '_hrf_', num2str(i));
                Opts.stim_id_ext = vertcat(Opts.stim_id_ext, tmp);
            end
            Opts.stim_id_ext  = sort(Opts.stim_id_ext);
        end
        
        % Construct design matrix for matlabbatch{1}.spm.stats.fmri_spec.sess(n).multi
        names = {};
        onsets = {};
        durations = {};
        
        for d = 1:length(Opts.stim_id)
            id = find(Opts.stim_id(d) == stim_list); % find rows that contain this stimulus
            names{1,d} = char(Opts.stim_id(d));
            onsets{1,d} = events.onset(id); 
            durations{1,d} = events.duration(id);
        end
        
        % Export design_mat and make filename for stim_id_file (we will use it later as a dictionary for our beta coefficients)
        save(design_multi,'names', 'onsets', 'durations');
        Dirs.design_multi = strcat(Dirs.results_dir, filesep, 'spm_design_multi.mat');
        
       
        
        % Collect design matrix metrics
        Opts.n_stim = length(Opts.stim_id);
        Opts.n_stim_betas = Opts.n_stim*(1 + sum(Opts.hrf_derivs));
        Opts.n_conf = length(Opts.confound_names);
        Opts.n_reg = Opts.n_stim_betas + Opts.n_conf + 1;
        
end