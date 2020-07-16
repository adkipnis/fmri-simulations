function Dirs = add_nii_files(Dirs, Opts, i, s, n)
    % 2.2.1 Directories
    % - Get path to nii_files
    nii_file_gz = Dirs.nii_files(n);
    nii_file = erase(nii_file_gz,".gz");
    [d, f, e] = fileparts(nii_file{1});
    Dirs.nii_file_s = {strcat(d, filesep, 's_', f, e)}; % s_ signifies smoothed

    mask_file_gz = Dirs.mask_files(n);
    Dirs.mask_file = erase(mask_file_gz,".gz");

    % - This is silly but necessary because reading nii.gz files is not yet implemented in SPM 12
    if ~exist(nii_file{1},'file')
        gunzip(nii_file_gz) 
    end
    
    if ~exist(Dirs.mask_file{1},'file')
        gunzip(mask_file_gz) 
    end
    
    % - (Re-)set output dirs where you save SPM.mat
    Dirs.subdir = char(fullfile(Dirs.outputdir,strcat('sub-', Dirs.sub_list{i}), Dirs.sesh_list{s}, strcat('run-', sprintf('%02s', string(n)))));
    if Opts.rewrite && exist(Dirs.subdir,'dir')
        rmdir(Dirs.subdir, 's');
    end
    spm_mkdir(Dirs.subdir);
    Dirs.run_scans = spm_select('Expand',nii_file); % create list with path to nifti file for every sample/scan
    Dirs.run_scans_s = spm_select('Expand', Dirs.nii_file_s); % create list with path to smoothed nifti file for every sample/scan
  
end