function [fmri_spec_c, scans_per_run, Dirs] = concatenate_pooled_runs(fmri_spec, Dirs, Opts, i)
    % Initialize variables
    fmri_spec_c = fmri_spec;
    n_runs = length(fmri_spec_c.sess);
    confounds_matrix = [];
    scans = {};
    scans_per_run = [];
    
    % Remove structure fields
    fmri_spec_c = rmfield(fmri_spec_c, 'sess');
    fmri_spec_c = rmfield(fmri_spec_c, 'mask');
    
    % New folder for concatenated runs
    Dirs.results_dir_concat = char(fullfile(Dirs.outputdir, strcat('sub-', Dirs.sub_list{i}), strcat('ses-', Opts.session_type, '-concatenated')));
    if Opts.rewrite && exist(Dirs.results_dir_concat,'dir')
        rmdir(Dirs.results_dir_concat, 's');
    end
    spm_mkdir(Dirs.results_dir_concat);
    fmri_spec_c.dir = {Dirs.results_dir_concat};
    
    % Get design matrix for all runs
    % TODO
    
    % Get confound design matrix for all runs
    for run = 1:n_runs
        confound_txt_path = fmri_spec.sess(run).multi_reg;
        confounds_matrix_run = csvread(confound_txt_path{1});
        if length(confounds_matrix) == 0
            confounds_matrix = confounds_matrix_run;
        else
            confounds_matrix = vertcat(confounds_matrix, confounds_matrix_run);
        end
    end
    Dirs.confounds_spm_file = fullfile(Dirs.results_dir_concat,'spm_confounds.txt');
    dlmwrite(Dirs.confounds_spm_file, confounds_matrix) % Save confound regressor matrix   
    fmri_spec_c.sess(1).multi_reg = {char(Dirs.confounds_spm_file)};
    
    % Reduce brain mask list
    fmri_spec_c.mask(1) = fmri_spec.mask(1);
    
    % Set HPF
    fmri_spec_c.sess(1).hpf = fmri_spec.sess(1).hpf;
    
    % Stack paths to scans 
    for run = 1:n_runs
        scans_run = fmri_spec.sess(run).scans;
        scans_per_run = [scans_per_run, length(scans_run)];
        
        if length(scans) == 0
            scans = scans_run;
        else
            scans = vertcat(scans, scans_run);
        end
    end
    fmri_spec_c.sess(1).scans = scans;
    
end