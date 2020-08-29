%===============================================================
%     Post-fMRIPrep SPM12 first-level noise-regressors only GLM
%     (on concatenated runs)
%===============================================================

%% 1. Preparations
%----- Custom paths
clc
clear all
format long g
addpath('/home/alex/matlab/Toolboxes/spm12');
addpath('/home/alex/matlab/SPM Batchscripts/GLM_Utils/'); % path to utility functions for this script
Dirs.BIDSdir = '/home/alex/Datasets/ds001246/';

%------ Options
Opts = struct();
Opts.task = 'perception';
Opts.subtask = 'Test';
Opts.session_type = [Opts.task, Opts.subtask];
Opts = load_task_json(Opts, Dirs); % add task-related MRI specs to Opts
Opts.pool_inference = false; % Pool runs for model estimation and inference (don't change this, use *_pooled version instead!)
Opts.resmooth = false; % redo smoothing even if smoothed images exist
Opts.fwhm_s = 0; % Full-width at half maximum of Gaussian spatial high pass filter in mm - set to 0 if no spatial smoothing is wished
Opts.fwhm_t = Inf; % Full-width at half maximum of Gaussian temporal high pass filter in seconds - set to Inf if no hpf is wished
Opts.smooth_prefix = strcat('fwhm_',num2str(Opts.fwhm_s), '_'); % Add prefix to smoothed files
Opts.hrf_derivs = [0 0]; % Estimate coefficients for HRF derivatives (first and second)
Opts.thresh_desc = 'FWE'; % options: 'FWE', 'FDR', 'none' -> control FWER or FDR before including voxels to clusters
Opts.alpha_level = 0.05; % For statistical testing of clusters
Opts.test_noise_regressors = false; % perform additional F-Tests on groups of noise regressors
Opts.save_fitted_response = 'none'; % options: 'predicted', 'corrected', 'none'

%------ Conditions, confounds and directories
Opts.use_global_signals = 'none'; % (csf, white matter, global signal) options: 'none', 'basic', 'deriv', 'power', 'full'
Opts.use_motion_regressors = 'deriv'; % options: 'none', 'minimal, 'basic', 'deriv', 'power', 'full'
Opts.n_pcs = 6; % number of noise PCs to include as noise regressors
Opts.n_cos = 6; % number cosine bases for modeling drift
[Opts, Regs] = get_confound_names(Opts);
Dirs = parse_bids_base_res(Dirs.BIDSdir, Opts); % Parse BIDS directory
spm('Defaults','fMRI'); %Initialise SPM fmri
spm_jobman('initcfg');  %Initialise SPM batch mode



%% 2. Stage one: Regress out drift noise from concatenated runs
for i = 3 : Dirs.n_subs
    cd(Dirs.BIDSdir)
    Opts.rewrite = false; % overwrites previously saved outputs
    Dirs = parse_bids_sub(Dirs, Opts, i);
    r = 0;
    Opts.unique_conditions_to_include = logical(zeros(50,1)); % after sorting the unique conditions in your events file, this logical vector filters out unwanted conditions; set this to logical(ones(number_of_unique_conditions, 1)) if you want to include all

    for s = 1 : Dirs.n_ses  
        % 2.1 Create file lists
        Dirs = create_filelists_from_bids(Dirs, i, s);

        for n = 1 : Dirs.n_runs
            fprintf('Processing subject %d''s run #%d of session "%s".\n', i, n, Dirs.sesh_list{s})
            
            r = r+1;

            % 2.2 Add paths to Nifti files
            Dirs = add_nii_files(Dirs, Opts, i, s, n);

            % 2.3 Apply spatial smoothing (if speficifed in Opts)
            Dirs = run_smoothing(Dirs, Opts);
            if isfield(Dirs, 'run_scans_s')       
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cellstr(Dirs.run_scans_s);
            elseif isfield(Dirs, 'run_scans')  
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cellstr(Dirs.run_scans);
            end
            
            % Apply temporal smoothing
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = Opts.fwhm_t; % temporal high-pass filter (Inf disables hpf)

            % 2.4 Construct design matrix and add confounds
            [Opts, Dirs] = construct_design_matrix(Dirs, Opts, n);
            Dirs = save_confounds(Dirs, Opts, n);

            % 2.5 Model specification
            fprintf('Model specification...\n')
            
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % onset times in secs (seconds) or scans (TRs);
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = Opts.TR;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50; % Microtime resoultion (here: number of slices)
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25; % Microtime onset (here: middle slice)            
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.mask(r) = Dirs.mask_file; 
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {char(Dirs.confounds_spm_file)};  
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = Opts.hrf_derivs; 
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; 
        end
    end
    spm_specify_c = spm_specify;
    [spm_specify_c.matlabbatch{1}.spm.stats.fmri_spec, scans_per_run, Dirs] = concatenate_pooled_runs(spm_specify.matlabbatch{1}.spm.stats.fmri_spec, Dirs, Opts, i);
    spm_jobman('run', spm_specify_c.matlabbatch)
    cd(Dirs.results_dir_concat)
    spm_fmri_concatenate('SPM.mat', scans_per_run);
            
            
    %% 2.7 Model estimation for concatenated runs
    fprintf('Parameter estimation...\n') 
    spm_estimate = {};
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = {[Dirs.results_dir_concat filesep 'SPM.mat']};  
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1; % write_residuals                
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % ReML

    spm_jobman('run',spm_estimate.matlabbatch)
    
%% 3. Stage two: First Level GLM on residuals
    Opts.unique_conditions_to_include = logical(ones(50,1));
    Opts.rewrite = true; % overwrites previously saved outputs
    r = 0;
    for s = 1 : Dirs.n_ses  
        % 3.1 Create file lists
        Dirs = create_filelists_from_bids(Dirs, i, s);

        for n = 1 : Dirs.n_runs
            fprintf('Processing subject %d''s run #%d of session "%s".\n', i, n, Dirs.sesh_list{s})
            r = r+1;
            % 3.2 Add paths to Nifti files
            Dirs = add_nii_files(Dirs, Opts, i, s, n);
            cd(Dirs.results_dir)
            
            % 3.3 Apply smoothing (if speficifed in Opts)
            spm_specify_r = {};
            
            Dirs.input_residuals = residuals_as_GLM_input(spm_specify_c, scans_per_run, r);
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(Dirs.input_residuals);

            % Apply temporal smoothing
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = Opts.fwhm_t; % temporal high-pass filter (Inf disables hpf)

            % 3.4 Construct design matrix and add confounds
            [Opts, Dirs] = construct_design_matrix(Dirs, Opts, n);
            Dirs = save_confounds(Dirs, Opts, n);

            % 3.5 Model specification
            fprintf('Model specification...\n')
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.dir = {Dirs.results_dir};
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % onset times in secs (seconds) or scans (TRs);
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = Opts.TR;
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50; % Microtime resoultion (here: number of slices)
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25; % Microtime onset (here: middle slice)            
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.mask = Dirs.mask_file; 
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {char(Dirs.design_multi)};
%             spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {char(Dirs.confounds_spm_file)};  
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = Opts.hrf_derivs; 
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            spm_specify_r.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; 
            
            spm_jobman('run', spm_specify_r.matlabbatch)
            
            
            %% 3.6 Model estimation
            fprintf('Parameter estimation...\n') 
            spm_estimate_r = {};
            spm_estimate_r.matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = {[Dirs.results_dir filesep 'SPM.mat']};  
            spm_estimate_r.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1; % write_residuals                
            spm_estimate_r.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % ReML

            spm_jobman('run',spm_estimate_r.matlabbatch)
        end
    end
end
