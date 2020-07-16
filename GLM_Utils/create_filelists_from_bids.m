function Dirs = create_filelists_from_bids(Dirs, Opts, i, s)
            Dirs.BIDSdir_raw = string(strcat(Dirs.BIDSdir, 'sub-', sprintf('%02s', string(i)), filesep,'ses-', Dirs.sesh_list(s), filesep, 'func'));            
            Dirs.BIDSdir_prep = string(strcat(Dirs.BIDSdir, 'derivatives',filesep,'fmriprep',filesep, 'sub-', sprintf('%02s', string(i)), filesep,'ses-', Dirs.sesh_list(s), filesep, 'func'));
            Dirs.raw_files = dir(Dirs.BIDSdir_raw);
            Dirs.deriv_files = dir(Dirs.BIDSdir_prep);
            Dirs.file_names_r = extractfield(Dirs.raw_files, 'name');
            Dirs.file_locs_r = extractfield(Dirs.raw_files, 'folder');
            Dirs.file_names_p = extractfield(Dirs.deriv_files, 'name');
            Dirs.file_locs_p = extractfield(Dirs.deriv_files, 'folder');

            % 2.1.1 Get the all preprocessed functional images and corresponding masks for subject i, the relevant session
            nii_filt = logical(contains(Dirs.file_names_p, 'preproc').*contains(Dirs.file_names_p, Opts.space).*contains(Dirs.file_names_p, 'nii.gz')); %only list all files containing these substrings
            Dirs.nii_files = strcat(Dirs.file_locs_p(nii_filt), filesep, Dirs.file_names_p(nii_filt));
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

end