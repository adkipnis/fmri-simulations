function Dirs = parse_bids_sub(Dirs, Opts, i)

Dirs.sub_dir = fullfile(Dirs.BIDSdir, strcat('sub-', sprintf('%02s', string(i))));
sub_dir_list = dir(Dirs.sub_dir);
Dirs.sesh_list_full = extractfield(sub_dir_list, 'name');
Dirs.sesh_filt = contains(Dirs.sesh_list_full, Opts.session_type); % filter out unneeded sessions
Dirs.sesh_list = Dirs.sesh_list_full(Dirs.sesh_filt);
Dirs.n_ses = sum(Dirs.sesh_filt);
end

