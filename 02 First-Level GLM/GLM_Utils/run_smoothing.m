function [spm_specify, Dirs] = run_smoothing(spm_specify, Dirs, Opts, r)
    
    % Spatial
    if Opts.fwhm_s > 0         
        if ~exist(Dirs.nii_file_s{1},'file') || Opts.resmooth
            if Opts.verbose, fprintf('Spatial smoothing...\n'), end 
            spm_smooth.matlabbatch{1}.spm.spatial.smooth.data = cellstr(Dirs.run_scans);
            spm_smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [Opts.fwhm_s, Opts.fwhm_s, Opts.fwhm_s]; % spatial filter width
            spm_smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0; % data type: same
            spm_smooth.matlabbatch{1}.spm.spatial.smooth.im = 0; % implicit masking: off
            spm_smooth.matlabbatch{1}.spm.spatial.smooth.prefix =  Opts.smooth_prefix;
            spm_jobman('run',spm_smooth.matlabbatch) % mellow down
        end
        Dirs.run_scans_s = spm_select('Expand', Dirs.nii_file_s); % create list with path to smoothed nifti file for every sample/scan
        spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cellstr(Dirs.run_scans_s);
    else
        spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cellstr(Dirs.run_scans);
    end

    % Temporal
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = Opts.fwhm_t; % temporal high-pass filter (Inf disables hpf)

end