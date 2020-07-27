function Dirs = parse_bids_sub(Dirs, Opts, i)

    Dirs.sub_dir = fullfile(Dirs.BIDSdir, strcat('sub-', sprintf('%02s', string(i))));
    sub_dir_list = dir(Dirs.sub_dir);
    Dirs.sesh_list_full = extractfield(sub_dir_list, 'name');
    Dirs.sesh_filt = contains(Dirs.sesh_list_full, Opts.session_type); % filter out unneeded sessions
    Dirs.sesh_list = Dirs.sesh_list_full(Dirs.sesh_filt);
    Dirs.n_ses = sum(Dirs.sesh_filt);

    % - (Re-)set subject-wise output dirs
    if Opts.pool_inference
        Dirs.subject_results_dir = char(fullfile(Dirs.outputdir, strcat('sub-', Dirs.sub_list{i}), strcat(Opts.session_type, '-results')));
        if Opts.rewrite && exist(Dirs.subject_results_dir,'dir')
            rmdir(Dirs.subject_results_dir, 's');
        end
        spm_mkdir(Dirs.subject_results_dir);
    end
end

