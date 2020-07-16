%===================================================================================
%     SPM12 First-level full GLM estimation post-fmriprep (w/ optional smoothing)
%===================================================================================

%% === 1. Preparations ===
% 1.1 Custom parameters
clear all
addpath('/moto/home/ak4572/spm12');
BIDSdir = '/moto/nklab/projects/ds001246/';
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
Opts.verbose = true;
Opts.session_type = 'perceptionTest';
Opts.space = 'T1w';
Opts.TR = 3; % Repetition time
Opts.rewrite = true; % overwrites previously saved outputs
Opts.t_smooth = true; % apply temporal smoothing
Opts.s_smooth = true; % apply spatial smoothing
Opts.resmooth = true; % redo smoothing even if smoothed images exist
Opts.fwhm_s = 3; % Full-width at half maximum of Gaussian spatial high pass filter in mm
Opts.fwhm_t = 128; % Full-width at half maximum of Gaussian temporal high pass filter in seconds
Opts.n_pcs = 6; % number of noise PCs to include as noise regressors
Opts.n_cos = 4; % number cosine bases for modeling drift
Opts.use_motion_regressors = true;
Opts.hrf_derivs = [1 0]; % Estimate coefficients for HRF derivatives (first and second)

% Confounds and directories
Opts.confound_names = get_confound_names({'global_signal', 'dvars', 'framewise_displacement'}, Opts);
Dirs = parse_bids(Dirs.BIDSdir, Opts.session_type, Opts.s_smooth); % Parse BIDS directory
Opts.n_ses = sum(Dirs.sesh_filt);

% Initialize SPM
spm('Defaults','fMRI'); %Initialise SPM fmri
spm_jobman('initcfg');  %Initialise SPM batch mode

%% === 2. Create SPM Batch structure
for i = 1 %: length(Dirs.sub_list)
    try
        for s = 1 %: Opts.n_ses
            
            % 2.1 Create file lists
            Dirs = create_filelists_from_bids(Dirs, Opts, i, s);
        
            for n = 1 %: length(Dirs.nii_files)
                if Opts.verbose, fprintf('Processing subject %d''s run #%d of session "%s".\n', i, n, Dirs.sesh_list{s}), end
                
                % 2.2 Add paths to Nifti files
                Dirs = add_nii_files(Dirs, Opts, i, s, n);

                % 2.3 Apply smoothing (if speficifed in Opts)
                spm_specify = {};
                spm_specify = run_smoothing(Dirs, Opts);
                        
                % 2.4 Construct design matrix and add confounds
                [design, Opts] = construct_design_matrix(Dirs, Opts, n);
                Dirs = save_confounds(Dirs, Opts, n);
                
                % 2.6 Model specification
                if Opts.verbose, fprintf('Model specification...\n'), end
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.dir = {Dirs.subdir};
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % onset times in secs (seconds) or scans (TRs);
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = Opts.TR;
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50; % Microtime resoultion (here: number of slices)
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25; % Microtime onset (here: middle slice)            
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.mask = {Dirs.mask_file{1}}; 
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = design;
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
%                 spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = {''};
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
                spm_estimate.matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = {[Dirs.subdir filesep 'SPM.mat']};  
                spm_estimate.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1; % write_residuals                
                spm_estimate.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % ReML
                
                spm_jobman('run',spm_estimate.matlabbatch)
                
                
                %% 2.8 F-Contrasts (Just for sanity check)
                if Opts.verbose, fprintf('F-Test...\n'), end 
                spm_contrasts = {};
                spm_contrasts.matlabbatch{1}.spm.stats.con.spmmat = {[Dirs.subdir filesep 'SPM.mat']}; 
                spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'Effects of interest';
                spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights =  effects_of_interest(Opts); %eye(Opts.n_reg)
                spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
                spm_contrasts.matlabbatch{1}.spm.stats.con.delete = 0;
                
                spm_jobman('run',spm_contrasts.matlabbatch);

                %% 2.9 Results
                if Opts.verbose, fprintf('Results...\n'), end
                spm_results = {};
                spm_results.matlabbatch{1}.spm.stats.results.spmmat = {[Dirs.subdir filesep 'SPM.mat']}; 
                spm_results.matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
                spm_results.matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
                spm_results.matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
                spm_results.matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
                spm_results.matlabbatch{1}.spm.stats.results.conspec.extent = 0;
                spm_results.matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
                spm_results.matlabbatch{1}.spm.stats.results.conspec.mask = {Dirs.mask_file{1}}; 
                spm_results.matlabbatch{1}.spm.stats.results.units = 1;
                spm_results.matlabbatch{1}.spm.stats.results.export{1}.ps = true;
                
                spm_jobman('run',spm_results.matlabbatch);
                %% Finally: Run matlabbatch jobs
                spm_jobman('run',matlabbatch);

            end
        
            
        end
    catch
        continue
    end
end