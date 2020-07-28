function generate_predicted_image(SPM, xSPM, Opts, Ic)

    if ~exist('Ic','var')    
        Ic = 1;
    end
    
    image_name = SPM.xCon(Ic).name;
    
    if ~strcmp(image_name, xSPM.title)
        error("Contrast name mismatch, check provided SPM and xSPM.")
    end
    
    %%%

    % Specify what spm_graph is supposed to do:
    % get fitted responses, use predicted time series (alternative: adjusted, i.e. with noise predictors regressed out)
    % Based on 'scans' (i.e., on voxel coordinates) and using the Contrast Index Ic.
    Opts.time_series_specs = struct('def', 'Fitted responses', ...
        'spec', struct('Ic', Ic, 'predicted', 1, 'x', struct('scan', true))); 

    
    
    Y = [];
    y = [];
    XYZ = xSPM.XYZ;
    for v=1:length(XYZ)
        [Y(:,v),y(:,v),~,~,~] = spm_graph(SPM, XYZ(:,v), Opts.time_series_specs); 
    end
    % Y  - fitted   data for the selected voxel
    % y  - adjusted data for the selected voxel

    % Save predicted responses
    for scan = 1:SPM.nscan         
        spm_write_filtered(Y(1,:), XYZ, xSPM.DIM, xSPM.M,...
        sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
        [image_name,'_Pred_', sprintf('%04s', string(scan)), '.nii']);
    end
end