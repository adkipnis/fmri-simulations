function Dirs = save_confounds(Dirs, Opts, n)
% - Design matrix: Confounds
    Dirs.confounds_spm_file = fullfile(Dirs.subdir,'spm_confounds.txt');
    beta_id_file = fullfile(Dirs.subdir,'stim_ids.txt');
    
    if ~exist(Dirs.confounds_spm_file,'file') || Opts.rewrite
        confounds = spm_load(Dirs.confound_files{n});
        confounds_matrix = [];
        for c = 1:length(Opts.confound_names)
            confounds_matrix = horzcat(confounds_matrix, confounds.(Opts.confound_names{c}));
        end % add each confound vector as column
        confounds_matrix(isnan(confounds_matrix)) = 0; % replaces NaN values by zeros
        dlmwrite(Dirs.confounds_spm_file, confounds_matrix) % Save confound regressor matrix               
        T = table(vertcat(Opts.stim_id, Opts.confound_names', "intercept"), 'VariableNames', {'RegressorNames'}); % Concatenate the confound names with the stimulus IDs in a table
        writetable(T, beta_id_file); % Write table to .txt for later use as a beta coefficient dictionary
    end

end