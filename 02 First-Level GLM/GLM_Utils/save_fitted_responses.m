function save_fitted_responses(SPM, xSPM, Opts, Ic)

    if ~exist('Ic','var')    
        Ic = 1;
    end
    
    image_name = SPM.xCon(Ic).name;
    
    if ~strcmp(image_name, xSPM.title)
        error("Contrast name mismatch, check provided SPM and xSPM.")
    end
    
    %%%

    Y = [];
    XYZ = xSPM.XYZ;

    beta  = spm_data_read(SPM.Vbeta,'xyz',XYZ);

    if strcmp(Opts.save_fitted_response,'predicted')

        % fitted (predicted) data (Y = X1*beta)
        %--------------------------------------------------------------
        Y = SPM.xX.xKXs.X*SPM.xCon(Ic).c*pinv(SPM.xCon(Ic).c)*beta;
   
    elseif strcmp(Opts.save_fitted_response,'corrected')

        % fitted (corrected)  data (Y = X1o*beta)
        %--------------------------------------------------------------
        Y = spm_FcUtil('Yc',SPM.xCon(Ic),SPM.xX.xKXs,beta);

    end

    % Save predicted responses
    for scan = 1:SPM.nscan         
        spm_write_filtered(Y(scan,:), XYZ, xSPM.DIM, xSPM.M,...
        sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
        [image_name,'_Pred_', sprintf('%04s', string(scan)), '.nii']);
    end
end