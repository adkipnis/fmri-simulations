function Dirs = create_filelists_from_bids(Dirs, Opts, i, s)
            Dirs.func_dir_raw = string(fullfile(Dirs.sub_dir, Dirs.sesh_list(s), 'func'));            
            Dirs.func_dir_prep = string(fullfile(Dirs.BIDSdir, 'derivatives','fmriprep', strcat('sub-', sprintf('%02s', string(i))), Dirs.sesh_list(s), 'func'));
            Dirs.raw_files = dir(Dirs.func_dir_raw);
            Dirs.deriv_files = dir(Dirs.func_dir_prep);
            Dirs.file_names_r = extractfield(Dirs.raw_files, 'name');
            Dirs.file_locs_r = extractfield(Dirs.raw_files, 'folder');
            Dirs.file_names_p = extractfield(Dirs.deriv_files, 'name');
            Dirs.file_locs_p = extractfield(Dirs.deriv_files, 'folder');

            % 2.1.1 Get the all preprocessed functional images and corresponding masks for subject i, the relevant session
            nii_filt = logical(contains(Dirs.file_names_p, 'preproc').*contains(Dirs.file_names_p, Opts.space).*contains(Dirs.file_names_p, 'nii.gz')); %only list all files containing these substrings
            Dirs.nii_files = strcat(Dirs.file_locs_p(nii_filt), filesep, Dirs.file_names_p(nii_filt));
            Dirs.n_runs = length(Dirs.nii_files);
            nii_filt_mask = logical(contains(Dirs.file_names_p, 'mask').*contains(Dirs.file_names_p, Opts.space).*contains(Dirs.file_names_p, 'nii.gz')); %only list all files containing these substrings
            Dirs.mask_files = strcat(Dirs.file_locs_p(nii_filt_mask), filesep, Dirs.file_names_p(nii_filt_mask));
            
            % 2.1.2 Corresponding confound regressors
            confound_filt = logical(contains(Dirs.file_names_p, 'confounds_regressors.tsv'));
            Dirs.confound_files = strcat(Dirs.file_locs_p(confound_filt), filesep, Dirs.file_names_p(confound_filt));

            % 2.1.3 Corresponding event files
            event_filt = logical(contains(Dirs.file_names_r, 'events.tsv'));
            Dirs.event_files = strcat(Dirs.file_locs_r(event_filt), filesep, Dirs.file_names_r(event_filt));
            
            % 2.1.4 Get corresponding metadata, including TR:
            %meta_locs_all = spm_BIDS(BIDS,'metadata','sub',sub_list(i),'run',run_list(n),'task','perception','type','bold');
            %meta_locs = meta_locs_all(sesh_filt(2:end));
            
%             % - (Re-)set output dirs session-wise output dirs
%             Dirs.session_results_dir = char(fullfile(Dirs.outputdir, strcat('sub-', Dirs.sub_list{i}), Dirs.sesh_list{s}, 'session-results'));
%             if Opts.rewrite && exist(Dirs.session_results_dir,'dir')
%                 rmdir(Dirs.session_results_dir, 's');
%             end
%             spm_mkdir(Dirs.session_results_dir);

end