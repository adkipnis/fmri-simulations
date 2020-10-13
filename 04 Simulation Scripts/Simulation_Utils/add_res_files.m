function Dirs = add_res_files(Dirs, Opts, i, s, n)
    Dirs.output_dir = char(fullfile(Dirs.outputdir, strcat('sub-', Dirs.sub_list{i}), Dirs.sesh_list{s}, strcat('run-', sprintf('%02s', string(n)))));
    if Opts.rewrite && exist(Dirs.output_dir,'dir')
        rmdir(Dirs.output_dir, 's');
    end
    spm_mkdir(Dirs.output_dir);
    
    Dirs.input_dir = char(fullfile(Dirs.GLM_results, strcat('sub-', Dirs.sub_list{i}), Dirs.sesh_list{s}, strcat('run-', sprintf('%02s', string(n)))));
    Dirs.input_files = dir(Dirs.input_dir);
    input_file_names = extractfield(Dirs.input_files, 'name');
    input_file_locs = extractfield(Dirs.input_files, 'folder');
    res_filt = logical(contains(input_file_names, 'Res_'));
    Dirs.res_files = strcat(input_file_locs(res_filt), filesep, input_file_names(res_filt));
    Dirs.n_res = sum(res_filt);
    if Dirs.n_res ~= 178
        warning('Wrong amount of Residual files found');
    end
  
end