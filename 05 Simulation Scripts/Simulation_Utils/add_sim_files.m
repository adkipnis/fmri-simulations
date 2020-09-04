function Dirs = add_sim_files(Dirs, i, s, n)
    Dirs.output_dir = char(fullfile(Dirs.outputdir, strcat('sub-', Dirs.sub_list{i}), Dirs.sesh_list{s}, strcat('run-', sprintf('%02s', string(n)))));
    Dirs.output_dir_pendant = char(fullfile(Dirs.GLM_results, strcat('sub-', Dirs.sub_list{i}), Dirs.sesh_list{s}, strcat('run-', sprintf('%02s', string(n)))));
    Dirs.input_files = dir(Dirs.output_dir);
    input_file_names = extractfield(Dirs.input_files, 'name');
    input_file_locs = extractfield(Dirs.input_files, 'folder');
    res_filt = logical(contains(input_file_names, 'Res_perm'));
    Dirs.res_perm_files = strcat(input_file_locs(res_filt), filesep, input_file_names(res_filt));
    sig_filt = logical(contains(input_file_names, 'Signal'));
    Dirs.signal_file = strcat(input_file_locs(sig_filt), filesep, input_file_names(sig_filt));
    Dirs.n_permutations = sum(res_filt);
end