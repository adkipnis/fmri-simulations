%===============================================================
%     Post-fMRIPrep SPM12 first-level noise-regressors only GLM 
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
Opts.rewrite = true; % overwrites previously saved outputs
Opts.resmooth = false; % redo smoothing even if smoothed images exist
Opts.fwhm_s = 0; % Full-width at half maximum of Gaussian spatial high pass filter in mm - set to 0 if no spatial smoothing is wished
Opts.fwhm_t = Inf; % Full-width at half maximum of Gaussian temporal high pass filter in seconds - set to Inf if no hpf is wished
Opts.smooth_prefix = strcat('fwhm_',num2str(Opts.fwhm_s), '_'); % Add prefix to smoothed files
Opts.hrf_derivs = [0 0]; % Estimate coefficients for HRF derivatives (first and second)
Opts.thresh_desc = 'FWE'; % options: 'FWE', 'FDR', 'none' -> control FWER or FDR before including voxels to clusters
Opts.alpha_level = 0.05; % For statistical testing of clusters
Opts.save_fitted_response = 'none'; % options: 'predicted', 'corrected', 'none'

%------ Conditions, confounds and directories
Opts.unique_conditions_to_include = logical(zeros(50,1)); % after sorting the unique conditions in your events file, this logical vector filters out unwanted conditions; set this to logical(ones(number_of_unique_conditions, 1)) if you want to include all
Opts.test_noise_regressors = false; % perform additional F-Tests on groups of noise regressors
Opts.use_global_signals = 'none'; % (csf, white matter, global signal) options: 'none', 'basic', 'deriv', 'power', 'full'
Opts.use_motion_regressors = 'deriv'; % options: 'none', 'minimal, 'basic', 'deriv', 'power', 'full'
Opts.n_pcs = 6; % number of noise PCs to include as noise regressors
Opts.n_cos = 6; % number cosine bases for modeling drift
[Opts, Regs] = get_confound_names(Opts);
Dirs = parse_bids_base_res(Dirs.BIDSdir, Opts); % Parse BIDS directory


%% 2. Create SPM Batch structure
spm('Defaults','fMRI'); %Initialise SPM fmri
spm_jobman('initcfg');  %Initialise SPM batch mode

for i = 3 : Dirs.n_subs
    Dirs = parse_bids_sub(Dirs, Opts, i);
    
    for s = 1 : Dirs.n_ses  
        % 2.1 Create file lists
        Dirs = create_filelists_from_bids(Dirs, i, s);

        for n = 1 : Dirs.n_runs
            fprintf('Processing subject %d''s run #%d of session "%s".\n', i, n, Dirs.sesh_list{s})
spm_specify = {};

            % 2.2 Add paths to Nifti files
            Dirs = add_nii_files(Dirs, Opts, i, s, n);

            % 2.3 Apply spatial smoothing (if speficifed in Opts)
            Dirs = run_smoothing(Dirs, Opts, n);
            if isfield(Dirs, 'run_scans_s')       
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(Dirs.run_scans_s);
            elseif isfield(Dirs, 'run_scans')  
                spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(Dirs.run_scans);
            end
            
            % Apply temporal smoothing
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = Opts.fwhm_t; % temporal high-pass filter (Inf disables hpf)


            % 2.4 Construct design matrix and add confounds
            [Opts, Dirs] = construct_design_matrix(Dirs, Opts, n);
            Dirs = save_confounds(Dirs, Opts, n);

            % 2.5 Model specification
            fprintf('Model specification...\n')
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.dir = {Dirs.results_dir};
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % onset times in secs (seconds) or scans (TRs);
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = Opts.TR;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50; % Microtime resoultion (here: number of slices)
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25; % Microtime onset (here: middle slice)            
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.mask = Dirs.mask_file; 
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {char(Dirs.confounds_spm_file)};  
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = Opts.hrf_derivs; 
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            spm_specify.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; 
            
            spm_jobman('run', spm_specify.matlabbatch)
            
            
            %% 2.7 Model estimation
            fprintf('Parameter estimation...\n') 
            spm_estimate = {};
            spm_estimate.matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = {[Dirs.results_dir filesep 'SPM.mat']};  
            spm_estimate.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1; % write_residuals                
            spm_estimate.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % ReML

            spm_jobman('run',spm_estimate.matlabbatch)
            
            

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
                    F_test_noise_pipeline(Dirs.results_dir, contrast_name, regressor_list, Opts, print);

                    % Save SPM results
                    save(strcat('xSPM_', contrast_name,'_.mat'),'xSPM')
                    save(strcat('TabDat_', contrast_name,'_.mat'),'TabDat')

                    % Write nice table for each contrast analysis
                    better_results_table(TabDat, Dirs.results_dir, contrast_name, Opts);

                    % Save cluster map and significance map as nifti file
                    make_cluster_maps(xSPM, Dirs.results_dir, Opts);
                end
            end
        end
    end
end
