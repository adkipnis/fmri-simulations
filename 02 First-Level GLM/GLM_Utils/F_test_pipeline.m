function F_test_pipeline(results_dir, contrast_name, Opts, r)
    
    if Opts.pool_inference == false
        r = 1;
    elseif ~exist('r','var') && Opts.pool_inference == true
        error("Specify number of runs in F-Test!")
    end
    
    if Opts.test_noise_regressors
        spm_get_defaults('cmdline',true);
    else
        spm_get_defaults('cmdline',false);
    end
    %% F-Contrasts
    spm_contrasts = {};
    spm_contrasts.matlabbatch{1}.spm.stats.con.spmmat = {[results_dir filesep 'SPM.mat']}; 
    spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = contrast_name;
    
    if strcmp(contrast_name, 'Effects-of-interest') 
        spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights =  effects_of_interest(Opts, r); 
    elseif strcmp(contrast_name, 'Full-model')  
        spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights =  eye(Opts.n_reg); 
    else
        error("Specify a viable F-contrast scheme ('Effects-of-interest' or 'Full-model') or change code.")
    end
    spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    spm_contrasts.matlabbatch{1}.spm.stats.con.delete = 0;
  
    spm_jobman('run',spm_contrasts.matlabbatch);

    %% Default Results
    if Opts.verbose, fprintf('Results...\n'), end
    spm_results = {};
    spm_results.matlabbatch{1}.spm.stats.results.spmmat = {[results_dir filesep 'SPM.mat']}; 
    spm_results.matlabbatch{1}.spm.stats.results.conspec.titlestr = contrast_name;
    spm_results.matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
    spm_results.matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
    spm_results.matlabbatch{1}.spm.stats.results.conspec.thresh = Opts.alpha_level;
    spm_results.matlabbatch{1}.spm.stats.results.conspec.extent = 0;
    spm_results.matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
    spm_results.matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
%             spm_results.matlabbatch{1}.spm.stats.results.conspec.mask.image.name = Dirs.mask_file;
%             spm_results.matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0;
    spm_results.matlabbatch{1}.spm.stats.results.conspec.titlestr = {contrast_name};
    spm_results.matlabbatch{1}.spm.stats.results.units = 1;
    spm_results.matlabbatch{1}.spm.stats.results.export{1}.ps = true;

    spm_jobman('run',spm_results.matlabbatch);

end