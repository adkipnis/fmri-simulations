%===================================================================================
%     SPM12 First-level full GLM estimation post-fmriprep (w/ optional smoothing)
%===================================================================================

%% === 1. Preparations ===
% Custom parameters
clear all
format long g
addpath('/home/alex/matlab/Toolboxes/spm12');
addpath('/home/alex/matlab/SPM Batchscripts/GLM_Utils/'); % path to utility functions for this script
Dirs.BIDSdir = '/home/alex/Datasets/ds001246/';

% Options
Opts = struct();
Opts.verbose = true;
Opts.session_type = 'perceptionTest';
Opts.space = 'T1w';
Opts.TR = 3; % Repetition time
Opts.rewrite = true; % overwrites previously saved outputs
Opts.pool_inference = false; % Pool runs for model estimation and inference (don't change this, use *_multi version instead!)
Opts.t_smooth = false; % apply temporal smoothing
Opts.s_smooth = true; % apply spatial smoothing
Opts.resmooth = false; % redo smoothing even if smoothed images exist
Opts.fwhm_s = 3; % Full-width at half maximum of Gaussian spatial high pass filter in mm
Opts.fwhm_t = 128; % Full-width at half maximum of Gaussian temporal high pass filter in seconds
Opts.smooth_prefix = strcat('fwhm_',num2str(Opts.fwhm_s), '_'); % Add prefix to smoothed files
Opts.hrf_derivs = [0 0]; % Estimate coefficients for HRF derivatives (first and second)

% Confounds and directories
Opts.use_global_signals = 'none'; % (csf, white matter, global signal) options: 'none', 'basic', 'deriv', 'power', 'full'
Opts.use_motion_regressors = 'full'; % options: 'none', 'minimal, 'basic', 'deriv', 'power', 'full'
Opts.n_pcs = 6; % number of noise PCs to include as noise regressors
Opts.n_cos = 6; % number cosine bases for modeling drift
Opts.confound_names = get_confound_names(Opts);
Dirs = parse_bids_base(Dirs.BIDSdir, Opts); % Parse BIDS directory

% Initialize SPM
spm('Defaults','fMRI'); %Initialise SPM fmri
spm_jobman('initcfg');  %Initialise SPM batch mode

%% === 2. Create SPM Batch structure
for i = 1 %: Dirs.n_subs
    Dirs = parse_bids_sub(Dirs, Opts, i);
    
    for s = 1 : Dirs.n_ses  
        % 2.1 Create file lists
        Dirs = create_filelists_from_bids(Dirs, Opts, i, s);

        for n = 1 : Dirs.n_runs
            if Opts.verbose, fprintf('Processing subject %d''s run #%d of session "%s".\n', i, n, Dirs.sesh_list{s}), end

            % 2.2 Add paths to Nifti files
            Dirs = add_nii_files(Dirs, Opts, i, s, n);

            % 2.3 Apply smoothing (if speficifed in Opts)
            spm_specify = {};
            [spm_specify, Dirs] = run_smoothing(spm_specify, Dirs, Opts, 1);

            % 2.4 Construct design matrix and add confounds
            [design, Opts, Dirs] = construct_design_matrix(Dirs, Opts, n);
            Dirs = save_confounds(Dirs, Opts, n);

            % 2.5 Model specification
            if Opts.verbose, fprintf('Model specification...\n'), end
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.dir = {Dirs.results_dir};
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % onset times in secs (seconds) or scans (TRs);
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = Opts.TR;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50; % Microtime resoultion (here: number of slices)
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25; % Microtime onset (here: middle slice)            
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.mask = Dirs.mask_file; 
%             spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = design;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {char(Dirs.design_multi)};
%             spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = {''};
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {char(Dirs.confounds_spm_file)};  
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = Opts.hrf_derivs; 
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; 
            
            spm_jobman('run', spm_specify.matlabbatch)


            %% 2.7 Model estimation
            if Opts.verbose, fprintf('Parameter estimation...\n'), end 
            spm_estimate = {};
            spm_estimate.matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = {[Dirs.results_dir filesep 'SPM.mat']};  
            spm_estimate.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1; % write_residuals                
            spm_estimate.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % ReML

            spm_jobman('run',spm_estimate.matlabbatch)

            %% 2.8 F-Contrasts (Just for sanity check)
            if Opts.verbose, fprintf('F-Test...\n'), end 
            spm_contrasts = {};
            spm_contrasts.matlabbatch{1}.spm.stats.con.spmmat = {[Dirs.results_dir filesep 'SPM.mat']}; 
            spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'Effects-of-interest';
            spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights =  effects_of_interest(Opts); 
            spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
            spm_contrasts.matlabbatch{1}.spm.stats.con.delete = 0;

            spm_jobman('run',spm_contrasts.matlabbatch);

            %% 2.9 Default Results
            if Opts.verbose, fprintf('Results...\n'), end
            spm_results = {};
            spm_results.matlabbatch{1}.spm.stats.results.spmmat = {[Dirs.results_dir filesep 'SPM.mat']}; 
            spm_results.matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
            spm_results.matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
            spm_results.matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
            spm_results.matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
            spm_results.matlabbatch{1}.spm.stats.results.conspec.extent = 0;
            spm_results.matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            spm_results.matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
%             spm_results.matlabbatch{1}.spm.stats.results.conspec.mask.image.name = Dirs.mask_file;
%             spm_results.matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
            spm_results.matlabbatch{1}.spm.stats.results.conspec.titlestr = {'Effects of interest'};
            spm_results.matlabbatch{1}.spm.stats.results.units = 1;
            spm_results.matlabbatch{1}.spm.stats.results.print = true;
            
            spm_jobman('run',spm_results.matlabbatch);
            
            %% 2.10 Custom Results
            
            % Save SPM results output
            save('xSPM.mat','xSPM')
            save('TabDat.mat','TabDat')
            
            % Write nice table for each contrast analysis
            better_results_table(Dirs.results_dir, spm_results, TabDat);
            
            % Save cluster map and significance map as nifti file
            make_inferential_maps(Dirs.results_dir, 'FWE');
            make_inferential_maps(Dirs.results_dir, 'FDR');
            
        end
    end
end
