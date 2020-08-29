function run_scans = residuals_as_GLM_input(spm_specify_c, scans_per_run, r)
    % initialization
    run_scans = {};
    res_dir = spm_specify_c.matlabbatch{1}.spm.stats.fmri_spec.dir;
    cumulative_scans_per_run = cumsum(scans_per_run);

    % Interval for residuals
    ending_scan = cumulative_scans_per_run(r);
    if r>1
        starting_scan = ending_scan - scans_per_run(r-1) + 1;
    else
        starting_scan = 1;
    end

    % Create scan_names
    for scan_num = starting_scan:ending_scan
       run_scan = char(fullfile(res_dir, strcat('Res_', sprintf('%04s', string(scan_num)), '.nii')));
       run_scans = vertcat(run_scans, run_scan);
    end
end
    