function level_1_GLM(nii_file, i, s, n, p, snr, Dirs, Opts)
    % Set some Dirs
    Dirs.results_dir = char(fullfile(Dirs.outputdir, ['sub-', ...
        Dirs.sub_list{i}], Dirs.sesh_list{s}, ['run-', sprintf('%02s', ...
        string(n))], ['GLM_Data_perm_', Opts.sim_type, '_', sprintf('%04s', ...
        string(p)), '_snr_', char(string(snr))]));
    if Opts.rewrite && exist(Dirs.results_dir,'dir')
        rmdir(Dirs.results_dir, 's');
    end
    spm_mkdir(Dirs.results_dir);
    Dirs.run_scans = spm_select('Expand', nii_file); % create list with path to nifti file for every sample/scan
    Dirs.design_multi = fullfile(Dirs.input_dir, 'spm_design_multi.mat');

    % Model specification
    spm_specify = struct();
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(Dirs.run_scans);
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = Inf;
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.dir = {Dirs.results_dir};
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % onset times in secs (seconds) or scans (TRs);
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = Opts.TR;
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50; % Microtime resoultion (here: number of slices)
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25; % Microtime onset (here: middle slice)            
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(Dirs.input_dir, 'mask.nii')}; 
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {char(Dirs.design_multi)};
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; 
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    spm_specify.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; 
    spm_jobman('run', spm_specify.matlabbatch)

    % Model estimation
    spm_estimate = struct();
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.spmmat(1) = {[Dirs.results_dir filesep 'SPM.mat']};  
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1; % write_residuals                
    spm_estimate.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % ReML
    spm_jobman('run',spm_estimate.matlabbatch)
end