function mask = get_mask(Dirs)
    input_file_names = extractfield(Dirs.input_files, 'name');
    input_file_locs = extractfield(Dirs.input_files, 'folder');
    mask_filt = logical(contains(input_file_names, 'mask'));
    Dirs.mask_file = strcat(input_file_locs(mask_filt), filesep, input_file_names(mask_filt));
    mask = niftiread(Dirs.mask_file{1});
end