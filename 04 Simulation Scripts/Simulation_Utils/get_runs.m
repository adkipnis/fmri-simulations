function Dirs = get_runs(Dirs, s)
    ses_dir = string(fullfile(Dirs.sub_dir, Dirs.sesh_list(s)));            
    ses_dir_contents = dir(ses_dir);
    ses_dir_fnames = extractfield(ses_dir_contents, 'name');
    run_filt = logical(contains(ses_dir_fnames, 'run'));
    Dirs.n_runs = sum(run_filt);
end