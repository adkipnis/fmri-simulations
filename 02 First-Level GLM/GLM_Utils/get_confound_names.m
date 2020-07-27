function confound_names = get_confound_names(Opts)
    
    confound_names = {};
    
    % Set global signals confound names
    if ~strcmp(Opts.use_global_signals, 'none')
        global_signals = ["csf", "white_matter", "global_signal"];
        deriv_global_signals = strcat(global_signals, '_derivative1');
        power_global_signals = strcat(global_signals, '_power2');
        power_deriv_global_signals = strcat(deriv_global_signals, '_power2');
        
        confound_names = [confound_names, global_signals];
        
        % Add standard confounds
        if strcmp(Opts.use_global_signals, 'deriv') % add their first derivatives only
            confound_names = [confound_names, deriv_global_signals];
        elseif strcmp(Opts.use_global_signals, 'power') % add their powers only
            confound_names = [confound_names, power_global_signals];
        elseif strcmp(Opts.use_global_signals, 'full') % add their first derivatives, powers for both the basic and the derivative parameters
            confound_names = [confound_names, deriv_global_signals, power_global_signals, power_deriv_global_signals];
        end    
        confound_names = cellstr(confound_names);
    end
        
        
    % Set motion confound names
    if ~strcmp(Opts.use_motion_regressors, 'none')
        bulk_motion_params = ["framewise_displacement", "dvars"];
        basic_motion_params = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"];
        deriv_motion_params = strcat(basic_motion_params, '_derivative1');
        power_motion_params = strcat(basic_motion_params, '_power2');
        power_deriv_motion_params = strcat(deriv_motion_params, '_power2');
        
        % Basic motion parameters
        if strcmp(Opts.use_motion_regressors, 'minimal')
            confound_names = [confound_names, bulk_motion_params]; % add six basic parameters        
        else
            confound_names = [confound_names, bulk_motion_params, basic_motion_params]; % add six basic parameters   
        end
        
        % Additional motion parameters
        if strcmp(Opts.use_motion_regressors, 'deriv') % add their first derivatives only
            confound_names = [confound_names, deriv_motion_params];
        elseif strcmp(Opts.use_motion_regressors, 'power') % add their powers only
            confound_names = [confound_names, power_motion_params];
        elseif strcmp(Opts.use_motion_regressors, 'full') % add their first derivatives, powers for both the basic and the derivative parameters
            confound_names = [confound_names, deriv_motion_params, power_motion_params, power_deriv_motion_params];
        end    
        confound_names = cellstr(confound_names);
    end
    
    % Set anatomical noise PCs

    
    if Opts.n_pcs > 0
        for k = 1:Opts.n_pcs 
            anatomical_PCs(k) = string(strcat('a_comp_cor_', sprintf('%02s', string(k-1)))); %python indexing starts with 0
        end 
        confound_names = [confound_names, anatomical_PCs];
    end
    
%     if Opts.n_pcs > Opts.n_cos
%        fprintf("fMRIPrep does high-pass filtering before running CompCor. Therefore, when using CompCor regressors, the corresponding cosine_XX regressors should also be included in the design matrix.")
%        fprintf("Too few Cosine basis functions: setting n_cos = n_pcs.")
%        Opts.n_cos = Opts.n_cos;
%     end
    
    % Set discrete cosine-basis functions
    if Opts.n_cos > 0
        if Opts.t_smooth == true
            fprintf('Warning: You should not apply temporal filtering AND model DCTFs together! Skipping discrete cosine-basis functions...')
        else
            for k = 1:Opts.n_cos
                cosine_basis_functions(k) = string(strcat('cosine', sprintf('%02s', string(k-1))));
            end
            confound_names = [confound_names, cosine_basis_functions];
        end
    end
end