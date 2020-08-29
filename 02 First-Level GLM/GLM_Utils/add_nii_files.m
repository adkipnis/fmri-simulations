function Dirs = add_nii_files(Dirs, Opts, i, s, n)
    % 2.2.1 Directories
    % - Get path to nii_files
    nii_file_gz = Dirs.nii_files(n);
    nii_file = erase(nii_file_gz,".gz");
    [d, f, e] = fileparts(nii_file{1});
    Dirs.nii_file_s = cellstr(strcat(d, filesep, Opts.smooth_prefix, f, e)); 

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
    Dirs.results_dir = char(fullfile(Dirs.outputdir, strcat('sub-', Dirs.sub_list{i}), Dirs.sesh_list{s}, strcat('run-', sprintf('%02s', string(n)))));
    if Opts.rewrite && exist(Dirs.results_dir,'dir')
        rmdir(Dirs.results_dir, 's');
    end
    spm_mkdir(Dirs.results_dir);
    Dirs.run_scans = spm_select('Expand', nii_file); % create list with path to nifti file for every sample/scan
    
  
end