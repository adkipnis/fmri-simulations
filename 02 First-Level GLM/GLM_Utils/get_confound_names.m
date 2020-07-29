function [Opts, Regs] = get_confound_names(Opts)
    Regs = struct();
    Opts.confound_names = {};
    
    
    % Set global signals confound names
    if ~strcmp(Opts.use_global_signals, 'none')
        global_signals = ["csf", "white_matter", "global_signal"];
        deriv_global_signals = strcat(global_signals, '_derivative1');
        power_global_signals = strcat(global_signals, '_power2');
        power_deriv_global_signals = strcat(deriv_global_signals, '_power2');
        
        Regs.global_signals = global_signals;
        
        % Add standard confounds
        if strcmp(Opts.use_global_signals, 'deriv') % add their first derivatives only
            Regs.deriv_global_signals = deriv_global_signals;
        elseif strcmp(Opts.use_global_signals, 'power') % add their powers only
            Regs.power_global_signals = power_global_signals;
        elseif strcmp(Opts.use_global_signals, 'full') % add their first derivatives, powers for both the basic and the derivative parameters
            Regs.deriv_global_signals = deriv_global_signals;
            Regs.power_global_signals = power_global_signals;
            Regs.power_deriv_global_signals = power_deriv_global_signals;
        end    
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
            Regs.bulk_motion_params = bulk_motion_params;
        else
            Regs.bulk_motion_params = bulk_motion_params;
            Regs.basic_motion_params = basic_motion_params;
        end
        
        % Additional motion parameters
        if strcmp(Opts.use_motion_regressors, 'deriv') % add their first derivatives only
            Regs.deriv_motion_params = deriv_motion_params;
        elseif strcmp(Opts.use_motion_regressors, 'power') % add their powers only
            Regs.power_motion_params = power_motion_params;
        elseif strcmp(Opts.use_motion_regressors, 'full') % add their first derivatives, powers for both the basic and the derivative parameters
            Regs.deriv_motion_params = deriv_motion_params;
            Regs.power_motion_params = power_motion_params;
            Regs.power_deriv_motion_params = power_deriv_motion_params;
        end    
    end
    
    % Set anatomical noise PCs

    
    if Opts.n_pcs > 0
        for k = 1:Opts.n_pcs 
            anatomical_PCs(k) = string(strcat('a_comp_cor_', sprintf('%02s', string(k-1)))); %python indexing starts with 0
        end 
        Regs.anatomical_PCs = anatomical_PCs;
    end
    
%     if Opts.n_pcs > Opts.n_cos
%        fprintf("fMRIPrep does high-pass filtering before running CompCor. Therefore, when using CompCor regressors, the corresponding cosine_XX regressors should also be included in the design matrix.")
%        fprintf("Too few Cosine basis functions: setting n_cos = n_pcs.")
%        Opts.n_cos = Opts.n_cos;
%     end
    
    % Set discrete cosine-basis functions
    if Opts.n_cos > 0
        if Opts.fwhm_t < Inf
            warning('You should not apply temporal filtering AND model DCTFs together! Skipping discrete cosine-basis functions...')
            Opts.n_cos = 0;
        else
            for k = 1:Opts.n_cos
                cosine_basis_functions(k) = string(strcat('cosine', sprintf('%02s', string(k-1))));
            end
        Regs.cosine_basis_functions = cosine_basis_functions;
        end
    end
    
    % Gather all regressor names
    C = struct2cell(Regs);
    Opts.confound_names = horzcat(C{:});
end