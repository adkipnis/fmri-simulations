function Directories = parse_bids(BIDSdir, session_type, s_smooth)
    Directories.BIDSdir = BIDSdir
    Directories.BIDS = spm_BIDS(Directories.BIDSdir);
    Directories.sub_list = spm_BIDS(Directories.BIDS,'subjects'); 
    Directories.sesh_list_full = spm_BIDS(Directories.BIDS,'sessions'); % - Path to sessions
    Directories.sesh_filt = contains(Directories.sesh_list_full, session_type); % filter out unneeded sessions
    Directories.sesh_list = Directories.sesh_list_full(Directories.sesh_filt);
    if s_smooth
        Directories.outputdir = fullfile(Directories.BIDSdir, 'derivatives','SPM_s');
    else
        Directories.outputdir = fullfile(Directories.BIDSdir, 'derivatives','SPM');
    end
    spm_mkdir(Directories.outputdir,strcat('sub-', char(Directories.sub_list)));
end