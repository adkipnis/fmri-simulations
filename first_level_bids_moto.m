%===================================================================================
%     SPM12 First-level full GLM estimation post-fmriprep (w/ optional smoothing)
%===================================================================================

%% === 1. Preparations ===
% 1.1 Custom parameters
clear all
addpath('/moto/home/ak4572/spm12');
BIDSdir = '/moto/nklab/projects/ds001246/';
session_type = 'perceptionTest';
space = 'T1w';
TR = 3; % Repetition time

rewrite = true; % overwrites previously saved outputs
s_smooth = true; % apply spatial smoothing
t_smooth = false; % apply temporal smoothing
fwhm_s = 3; % Full-width at half maximum of Gaussian spatial high pass filter in mm
fwhm_t = 128; % Full-width at half maximum of Gaussian temporal high pass filter in seconds
n_pcs = 6; % number of noise PCs to include as noise regressors
n_cos = 4; % number cosine bases for modeling drift
confound_names = {'dvars', 'framewise_displacement'}; 



% 1.2 Directories
% Add first 6 noise PCs and model first 4 cosine bases for drift
if n_pcs > 0
    for k = 1:n_pcs 
        tmp1(k) = strcat('a_comp_cor_0', string(k-1)); end %python indexing starts with 0
else
    tmp1 = [];
end

if n_cos > 0
    for k = 1:n_cos
        tmp2(k) = strcat('cosine0', string(k-1)); end
else 
    tmp2 = [];
end

confound_names_tmp = [confound_names, tmp1, tmp2];
confound_names = confound_names_tmp;
clear confound_names_tmp tmp1 tmp2

% - Parse BIDS directory
BIDS = spm_BIDS(BIDSdir);
sub_list = spm_BIDS(BIDS,'subjects'); 

% - Path to sessions
sesh_list_full = spm_BIDS(BIDS,'sessions');
sesh_filt = contains(sesh_list_full,session_type); % filter out unneeded sessions
sesh_list = sesh_list_full(sesh_filt);
n_ses = sum(sesh_filt);

% - Output directories
if s_smooth
    outputdir = strcat(BIDSdir, 'derivatives',filesep,'SPM_s',filesep);
else
    outputdir = strcat(BIDSdir, 'derivatives',filesep,'SPM',filesep);
end
spm_mkdir(outputdir,strcat('sub-', char(sub_list)));

% - Initialize SPM
spm('Defaults','fMRI'); %Initialise SPM fmri
spm_jobman('initcfg');  %Initialise SPM batch mode

%% === 2. Create SPM Batch structure
for i = 1 : length(sub_list)
    try
        for s = 1 : n_ses
        
            %% 2.1 create file lists
            BIDSdir_raw = string(strcat(BIDSdir, 'sub-', sprintf('%02s', string(i)), filesep,'ses-', sesh_list(s), filesep, 'func'));
            BIDSdir_prep = string(strcat(BIDSdir, 'derivatives',filesep,'fmriprep',filesep, 'sub-', sprintf('%02s', string(i)), filesep,'ses-', sesh_list(s), filesep, 'func'));
            raw_files = dir(BIDSdir_raw);
            deriv_files = dir(BIDSdir_prep);
            file_names_r = extractfield(raw_files, 'name');
            file_locs_r = extractfield(raw_files, 'folder');
            file_names_p = extractfield(deriv_files, 'name');
            file_locs_p = extractfield(deriv_files, 'folder');

            % 2.1.1 Get the all preprocessed functional images and corresponding masks for subject i, the relevant session
            nii_filt = logical(contains(file_names_p, 'preproc').*contains(file_names_p, space).*contains(file_names_p, 'nii.gz')); %only list all files containing these substrings
            nii_files = strcat(file_locs_p(nii_filt), filesep, file_names_p(nii_filt));
            nii_filt_mask = logical(contains(file_names_p, 'mask').*contains(file_names_p, space).*contains(file_names_p, 'nii.gz')); %only list all files containing these substrings
            mask_files = strcat(file_locs_p(nii_filt_mask), filesep, file_names_p(nii_filt_mask));

            % 2.1.2 Corresponding confound regressors
            confound_filt = logical(contains(file_names_p, 'confounds_regressors.tsv'));
            confound_files = strcat(file_locs_p(confound_filt), filesep, file_names_p(confound_filt));

            % 2.1.3 Corresponding event files
            event_filt = logical(contains(file_names_r, 'events.tsv'));
            event_files = strcat(file_locs_r(event_filt), filesep, file_names_r(event_filt));


            % 2.1.4 Get corresponding metadata, including TR:
            %meta_locs_all = spm_BIDS(BIDS,'metadata','sub',sub_list(i),'run',run_list(n),'task','perception','type','bold');
            %meta_locs = meta_locs_all(sesh_filt(2:end));



            %% 2.2 Loop over runs
         
            for n = 1 : length(nii_files)

                fprintf('Processing subject %d''s run #%d of session "%s".\n', i, n, sesh_list{s}) 
                matlabbatch = {}; 
                
                % 2.2.1 Directories
                % - Get path to nii_files
                nii_file_gz = nii_files(n);
                nii_file = erase(nii_file_gz,".gz");
                [d, f, e] = fileparts(nii_file{1});
                nii_file_s = {strcat(d, filesep, 's_', f, e)};
                
                mask_file_gz = mask_files(n);
%                 % different mask file
%                 mask_file_gz = {'/home/alex/Documents/10. Semester - CU/Master Thesis/Python Code/ds001246/derivatives/FSL/training_glm_and_test_glm/single_trial_extraction_by_volumes/sub-01/sub-01_brain_mask.nii.gz'}
                mask_file = erase(mask_file_gz,".gz");

                % - This is silly but necessary because reading nii.gz files is not yet implemented in SPM 12
                if ~exist(nii_file{1},'file'), gunzip(nii_file_gz) 
                end
                if ~exist(mask_file{1},'file'), gunzip(mask_file_gz) 
                end

                % - (Re-)set output dirs where you save SPM.mat
                subdir=char(fullfile(outputdir,strcat('sub-', sub_list{i}), sesh_list{s}, strcat('run-', sprintf('%02s', string(n)))));
                if rewrite && exist(subdir,'dir')
                    rmdir(subdir, 's');
                end
                spm_mkdir(subdir);
                run_scans = spm_select('Expand',nii_file); % create list with path to nifti file for every sample/scan


                % 2.2.2 Spatial Smoothing
                if s_smooth             
                    if ~exist(nii_file_s{1},'file') || rewrite
                        fprintf('Spatial smoothing...\n') 
                        smooth.matlabbatch{1}.spm.spatial.smooth.data = cellstr(run_scans);
                        smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm_s fwhm_s fwhm_s]; % spatial filter width
                        smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0; % data type: same
                        smooth.matlabbatch{1}.spm.spatial.smooth.im = 0; % implicit masking: off
                        smooth.matlabbatch{1}.spm.spatial.smooth.prefix = 's_';
                        spm_jobman('run',smooth.matlabbatch) % mellow down
                        run_scans_s = spm_select('Expand',nii_file_s); % create list with path to smoothed nifti file for every sample/scan
                        matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(run_scans_s);
                    end
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(run_scans);
                end
                
                fprintf('Model specification...\n')
                 % 2.2.3 Basic parameters and pointers
                matlabbatch{1}.spm.stats.fmri_spec.dir = {subdir};
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % onset times in secs (seconds) or scans (TRs);
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
%                 matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16; % Microtime onset (relevant if Slice-Timing Correction has been performed)
%                 matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; % Microtime resolution (see above)
                                
                matlabbatch{1}.spm.stats.fmri_spec.mask = {mask_file{1}};

                % 2.2.4 Design matrix
                % - For sanity check: Create a design matrix that does not differentiate between pictures
                % design = struct('name', {'Picture'}, 'onset', {events.onset}, 'duration', {events.duration}, 'tmod', {0}, 'pmod', {''});

                % - Design matrix: Events (Check if design matrix file exists and create it if not)
                design_mat = strcat(subdir,filesep,'sub-', sub_list{i},'-',sesh_list{s},'-run-', string(n),'_spm_design.mat');

                if exist(design_mat,'file') && ~rewrite
                    load(design_mat);
                else
                    events = spm_load(event_files{n});

                    % extract unique event names and remove missing rows
                    stim_id = unique(string(events.stim_id));
                    stim_id = rmmissing(stim_id); 

    %                 % Alternatively: extract unique event names and rename missing entries to 'Baseline'
    %                 stim_id = string(events.stim_id);
    %                 stim_id(ismissing(stim_id(:,1))) = 'Baseline';
    %                 events.stim_id = stim_id;
    %                 stim_id = unique(string(events.stim_id));

                    design = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                    for d = 1:length(stim_id)
                        id = find(stim_id(d)==string(events.stim_id)); % find rows that contain this stimulus
                        design(d).name = char(stim_id(d));
                        design(d).onset = events.onset(id); 
                        design(d).duration = events.duration(id);
                        design(d).tmod = 0; % time modulation
                        design(d).pmod = {''}; % time modulation
                        design(d).orth = 1;
                    end

    %                 for d = 1:length(events.stim_id)-2 % we want to skip the first and last row of the events file
    %                     design(d).name = char(string(events.stim_id(d+1)));
    %                     design(d).onset = events.onset(d+1);
    %                     design(d).duration = events.duration(d+1);
    %                     design(d).tmod = 0; % time modulation
    %                     design(d).pmod = {''}; % time modulation
    %                     design(d).orth = 1;
    %                 end

                    % Export design_mat and make filename for stim_id_file (we will use it later as a dictionary for our beta coefficients)
                    save(design_mat,'design');
                    beta_id_file=strcat(subdir,filesep,'sub-', sub_list{i},'-',sesh_list{s},'-run-', string(n),'_stim_ids.txt');


                end

                matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = design;
                matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = {''};

                % - Design matrix: Confounds
                confounds_spm_file=strcat(subdir,filesep,'sub-', sub_list{i},'-',sesh_list{s},'-run-', string(n),'_spm_confounds.txt');

                if ~exist(confounds_spm_file,'file') || rewrite
                    confounds=spm_load(confound_files{n});
                    confounds_matrix = [];
                    for c = 1:length(confound_names)
                        confounds_matrix = horzcat(confounds_matrix, confounds.(confound_names{c}));
                    end % add each confound vector as column
                    confounds_matrix(isnan(confounds_matrix)) = 0; % replaces NaN values by zeros
                    dlmwrite(confounds_spm_file,confounds_matrix) % Save confound regressor matrix               
                    T = table(vertcat(stim_id, confound_names', "intercept"), 'VariableNames', {'RegressorNames'}); % Concatenate the confound names with the stimulus IDs in a table
                    writetable(T, beta_id_file); % Write table to .txt for later use as a beta coefficient dictionary
                end

                matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {char(confounds_spm_file)};
                
                if t_smooth
                    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = fwhm_t; % temporal high-pass filter 
                end

                % 2.2.5 Default Model Specifications
                matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % keep derivatives out due to multi-condition design
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

                %spm_jobman('run',matlabbatch)


                %% 2.2.6 Model estimation (Default)
                fprintf('Parameter estimation...\n') 
                matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = {[subdir filesep 'SPM.mat']};
                % write_residuals
                matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
                % method
                matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

                % 2.2.7 Contrasts (Just for sanity check)
    %             Nregr = 14; % first: picture viewing, 12 confound regressors, and model intercept
    %             matlabbatch{3}.spm.stats.con.spmmat = {[subdir filesep 'SPM.mat']};
    %             matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Picture viewing';
    %             matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = zeros(1, Nregr);
    %             matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights(1) = 1;
    %             matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    %             matlabbatch{3}.spm.stats.con.delete = 0;


                % Finally: Run matlabbatch jobs
                spm_jobman('run',matlabbatch);

            end
        
            
        end
    catch
        continue
    end
end