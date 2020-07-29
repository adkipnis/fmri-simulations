function make_cluster_maps(xSPM, results_dir, Opts)
    % Partially derived from https://github.com/solleo/myspm/edit/master/myspm_result.m
    
    con_title = xSPM.title;
    correction_type = Opts.thresh_desc;
    %% Get cluster data

    XYZ     = xSPM.XYZ;
    Z       = round(spm_clusters(XYZ));
    num     = max(Z);
    [clusSize, ni] = sort(histc(Z,1:num), 2, 'descend');
    n       = size(ni);
    n(ni)   = 1:num;
    Z       = n(Z); % renumbering by cluster size

    
    %% Map significant clusters and include cluster size info
    % Assign cluster size to voxel values instead of cluster ID
    Z_tmp = Z;
    for i = 1:length(clusSize)
        Z_tmp(Z==i) = clusSize(i);
    end
    
%     fname_sigclus=['sized_significant_clusters_', con_title, '_', correction_type];
%     fname_sigclus = strrep(fname_sigclus,' ','');
%     spm_write_filtered(Z_tmp, XYZ, xSPM.DIM, xSPM.M,...
%     sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
%     [fname_sigclus,'.nii']);
%     fprintf('Created cluster map in "%s".\n', [results_dir, filesep, [fname_sigclus,'.nii']])
    
    %% Map significant clusters in a binary form

    sigs = repmat(1, 1, length(Z));   
    
    fname_sig=['significant_clusters_', con_title, '_', correction_type];
    fname_sig = strrep(fname_sig,' ','');
    spm_write_filtered(sigs, XYZ, xSPM.DIM, xSPM.M,...
    sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
    [fname_sig,'.nii']);
    fprintf('Created significance map in "%s".\n', [results_dir, filesep, [fname_sig,'.nii']])
    
end