function Dirs = parse_bids_base(BIDSdir, Opts)
    Dirs.BIDSdir = BIDSdir;
    Dirs.BIDS = spm_BIDS(Dirs.BIDSdir);
    Dirs.sub_list = spm_BIDS(Dirs.BIDS,'subjects'); 
    Dirs.n_subs = length(Dirs.sub_list);
    if Opts.s_smooth
        Dirs.outputdir = fullfile(Dirs.BIDSdir, 'derivatives',strcat('SPM_', string(Opts.fwhm_s)));
    else
        Dirs.outputdir = fullfile(Dirs.BIDSdir, 'derivatives','SPM');
    end
    spm_mkdir(Dirs.outputdir,strcat('sub-', char(Dirs.sub_list)));
    
end