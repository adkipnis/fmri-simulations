function F_test_noise_pipeline(results_dir, contrast_name, regressor_list, Opts, print, r)
    
    if Opts.pool_inference == false
        r = 1;
    elseif ~exist('r','var') && Opts.pool_inference == true
        error("Specify number of runs in F-Test!")
    end
    
    if ~print
        spm_get_defaults('cmdline',true);
    else
        spm_get_defaults('cmdline',false);
    end
    
    if strcmp(Opts.thresh_desc, 'FDR')
        spm_get_defaults('stats.topoFDR',false);
     else
        spm_get_defaults('stats.topoFDR',true);
    end

    %% 2.7 F-Contrasts (Just for sanity check)
    spm_contrasts = {};
    spm_contrasts.matlabbatch{1}.spm.stats.con.spmmat = {[results_dir filesep 'SPM.mat']}; 
    spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = contrast_name;
    spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights =  effects_of_noise(Opts, regressor_list, r); 
    spm_contrasts.matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    spm_contrasts.matlabbatch{1}.spm.stats.con.delete = 0;
  
    spm_jobman('run',spm_contrasts.matlabbatch);

    %% Default Results
    fprintf('Results...\n')
    spm_results = {};
    spm_results.matlabbatch{1}.spm.stats.results.spmmat = {[results_dir filesep 'SPM.mat']}; 
    spm_results.matlabbatch{1}.spm.stats.results.conspec.titlestr = contrast_name;
    spm_results.matlabbatch{1}.spm.stats.results.conspec.contrasts = Inf;
    spm_results.matlabbatch{1}.spm.stats.results.conspec.threshdesc = Opts.thresh_desc;
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