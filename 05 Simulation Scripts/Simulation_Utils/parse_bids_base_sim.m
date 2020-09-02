function Dirs = parse_bids_base_sim(BIDSdir)
    Dirs.BIDSdir = BIDSdir;
    Dirs.BIDS = spm_BIDS(Dirs.BIDSdir);
    Dirs.sub_list = spm_BIDS(Dirs.BIDS,'subjects'); 
    Dirs.n_subs = length(Dirs.sub_list);
    suffix = '';
    Dirs.outputdir = fullfile(Dirs.BIDSdir, 'derivatives',strcat('Sim', suffix));
    spm_mkdir(Dirs.outputdir,strcat('sub-', char(Dirs.sub_list)));
end