function make_inferential_maps(Dirs)
    % Partially derived from https://github.com/solleo/myspm/edit/master/myspm_result.m

    load([Dirs.results_dir,'/SPM.mat']);
    load([Dirs.results_dir,'/xSPM.mat']);
    con_title = SPM.xCon(1).name;
    
    % Get cluster data
    XYZ     = xSPM.XYZ;
    Z       = round(spm_clusters(XYZ));
    num     = max(Z);
    [clusSize, ni] = sort(histc(Z,1:num), 2, 'descend');
    n       = size(ni);
    n(ni)   = 1:num;
    Z       = n(Z); % renumbering by cluster size
    
    
    % Assign cluster size to voxel values instead of cluster ID
    Z_tmp = Z;
    for i = 1:length(clusSize)
        Z_tmp(Z==i) = clusSize(i);
    end
    
    
    % Sanity check
    if (min(xSPM.Z)<xSPM.u) || (min(clusSize)<xSPM.k)
    error('Non-significant cluster included!')
    end
    
    
    % Save cluster map
    fname_sigclus=['sigclus_', con_title];
    fname_sigclus = strrep(fname_sigclus,' ','');
    spm_write_filtered(Z_tmp, XYZ, xSPM.DIM, xSPM.M,...
    sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
    [fname_sigclus,'.nii']);
    fprintf('Created cluster map in "%s".\n', [Dirs.results_dir, filesep, [fname_sigclus,'.nii']])
    
    
    % Save significance map
    sigs = repmat(1, 1, length(Z));   
    fname_sig=['sig_', con_title];
    fname_sig = strrep(fname_sig,' ','');
    spm_write_filtered(sigs, XYZ, xSPM.DIM, xSPM.M,...
    sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
    [fname_sig,'.nii']);
    fprintf('Created significance map in "%s".\n', [Dirs.results_dir, filesep, [fname_sig,'.nii']])
    
end