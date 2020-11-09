function Dirs = parse_bids_base_name(Dirs, name)
    Dirs.BIDS = spm_BIDS(Dirs.BIDSdir);
    Dirs.sub_list = unique({Dirs.BIDS.subjects.name});
    for i = 1:length(Dirs.sub_list)
        Dirs.sub_list{i} = erase(Dirs.sub_list{i}, 'sub-');
    end
    Dirs.n_subs = length(Dirs.sub_list);
    suffix = '';
    Dirs.outputdir = fullfile(Dirs.BIDSdir, 'derivatives',strcat(name, suffix));
    spm_mkdir(Dirs.outputdir,strcat('sub-', char(Dirs.sub_list)));
end