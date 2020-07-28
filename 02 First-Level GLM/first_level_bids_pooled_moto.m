%===================================================================================
%     SPM12 First-level full GLM estimation post-fmriprep (w/ optional smoothing)
%===================================================================================

%% === 1. Preparations ===
% Custom parameters
clear all
format long g
addpath('/moto/home/ak4572/spm12');
addpath('/moto/home/ak4572/GLM_Utils/'); % path to utility functions for this script
Dirs.BIDSdir = '/moto/nklab/projects/ds001246/';

% Options
Opts = struct();
Opts.verbose = true;
Opts.session_type = 'perceptionTest';
Opts.space = 'T1w';
Opts.TR = 3; % Repetition time
Opts.rewrite = true; % overwrites previously saved outputs
Opts.pool_inference = true; % Pool runs for model estimation and inference (don't change this, use *_mono version instead!)
Opts.delete_estimates = false; % will delete estimates betas and residues if inference is pooled
Opts.t_smooth = false; % apply temporal smoothing
Opts.s_smooth = true; % apply spatial smoothing
Opts.resmooth = false; % redo smoothing even if smoothed images exist
Opts.fwhm_s = 3; % Full-width at half maximum of Gaussian spatial high pass filter in mm
Opts.fwhm_t = 128; % Full-width at half maximum of Gaussian temporal high pass filter in seconds
Opts.smooth_prefix = strcat('fwhm_',num2str(Opts.fwhm_s), '_'); % Add prefix to smoothed files
Opts.hrf_derivs = [0 0]; % Estimate coefficients for HRF derivatives (first and second)
Opts.alpha_level = 0.05; % For statistical testing

% Confounds and directories
Opts.test_noise_regressors = true; % perform additional F-Tests on groups of noise regressors
Opts.use_global_signals = 'none'; % (csf, white matter, global signal) options: 'none', 'basic', 'deriv', 'power', 'full'
Opts.use_motion_regressors = 'full'; % options: 'none', 'minimal, 'basic', 'deriv', 'power', 'full'
Opts.n_pcs = 6; % number of noise PCs to include as noise regressors
Opts.n_cos = 6; % number cosine bases for modeling drift
[Opts.confound_names, Regs] = get_confound_names(Opts);
Dirs = parse_bids_base(Dirs.BIDSdir, Opts); % Parse BIDS directory

% Initialize SPM
spm('Defaults','fMRI'); %Initialise SPM fmri
spm_jobman('initcfg');  %Initialise SPM batch mode

%% === 2. Create SPM Batch structure
for i = 1 : Dirs.n_subs
    if Opts.verbose, fprintf('Processing all of subject %d''s runs session-type "%s".\n', i, Opts.session_type), end
    Dirs = parse_bids_sub(Dirs, Opts, i);
    r = 0; % successive run-counter (across sessions)
    spm_specify = {};
    
    for s = 1 : Dirs.n_ses  
        % 2.1 Create file lists
        Dirs = create_filelists_from_bids(Dirs, Opts, i, s);

        for n = 1 : Dirs.n_runs
            
            r = r+1;
            
            % 2.2 Add paths to Nifti files
            Dirs = add_nii_files(Dirs, Opts, i, s, n);

            % 2.3 Apply smoothing (if speficifed in Opts)
            [spm_specify, Dirs] = run_smoothing(spm_specify, Dirs, Opts, r);

            % 2.4 Construct design matrix and add confounds
            [design, Opts, Dirs] = construct_design_matrix(Dirs, Opts, n);
            Dirs = save_confounds(Dirs, Opts, n);

            % 2.5 Model specification
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.dir = {Dirs.subject_results_dir};
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % onset times in secs (seconds) or scans (TRs);
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = Opts.TR;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50; % Microtime resoultion (here: number of slices)
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25; % Microtime onset (here: middle slice)            
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.mask = Dirs.mask_file; 
%             spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond = design;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {char(Dirs.design_multi)};
%             spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = {''};
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {char(Dirs.confounds_spm_file)};  
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = Opts.hrf_derivs; 
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; 
        end
    end
    
    spm_jobman('run', spm_specify.matlabbatch)
       
        
    %% 2.6 Model estimation
    if Opts.verbose, fprintf('Parameter estimation...\n'), end 
    spm_estimate = {};
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = {[Dirs.subject_results_dir filesep 'SPM.mat']};  
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1; % write_residuals                
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % ReML

    spm_jobman('run',spm_estimate.matlabbatch)

    %% 2.7 Signal extra SS F-Test pipeline (Contrasts, Results and Exports)
    if Opts.verbose, fprintf('F-Testing...\n'), end 

    % Set contrast name
    contrast_name = 'Effects-of-interest';

    % Specify F-constrasts and get default SPM results
    F_test_pipeline(Dirs.subject_results_dir, contrast_name, Opts, r);

    % Save SPM results
    save(strcat('xSPM_', contrast_name,'.mat'),'xSPM')
    save(strcat('TabDat_', contrast_name,'.mat'),'TabDat')

    % Write nice table for each contrast analysis
    better_results_table(TabDat, Dirs.subject_results_dir, contrast_name, Opts);

    % Save cluster map and significance map as nifti file
    make_inferential_maps(xSPM, Dirs.subject_results_dir, 'FWE');
    make_inferential_maps(xSPM, Dirs.subject_results_dir, 'FDR');

    %% 2.8 Optional: Noise regressors extra SS F-Test pipelines

    if Opts.test_noise_regressors         
        contrast_names = fieldnames(Regs);
        regressor_lists = struct2cell(Regs);
        print = false; % we only want to print out the last ps table (otherwise we get recursive bloated versions thereof)

        % Loop over all sets of regressors
        for reg = 1:length(regressor_lists)
            % Set contrast name
            contrast_name = contrast_names{reg,1};
            regressor_list = regressor_lists{reg,1};

            % Print only the last results table (containing all previous results)
            if reg == length(regressor_lists)
                print = true;
            end

            % Specify F-constrasts and get default SPM results
            F_test_noise_pipeline(Dirs.subject_results_dir, contrast_name, regressor_list, Opts, print, r);

            % Save SPM results
            save(strcat('xSPM_', contrast_name,'_.mat'),'xSPM')
            save(strcat('TabDat_', contrast_name,'_.mat'),'TabDat')

            % Write nice table for each contrast analysis
            better_results_table(TabDat, Dirs.subject_results_dir, contrast_name, Opts);

            % Save cluster map and significance map as nifti file
            make_inferential_maps(xSPM, Dirs.subject_results_dir, 'FWE');
        end
    end
    
    %% Clean up
    if Opts.delete_estimates
        cd(Dirs.subject_results_dir)
        delete Beta*.*
        delete Res*.*
    end
end
