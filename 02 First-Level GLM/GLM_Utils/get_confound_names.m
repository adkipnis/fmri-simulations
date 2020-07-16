function confound_names = get_confound_names(confound_names_init, Opt)
    
    confound_names = confound_names_init;
    
    % Motion parameters
    if Opt.use_motion_regressors
        motion_confound_names = {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'};
        confound_names = [confound_names, motion_confound_names];
    end
    
    % number of aCompCors
    if Opt.n_pcs > 0
        for k = 1:Opt.n_pcs 
            tmp1(k) = strcat('a_comp_cor_0', string(k-1)); %python indexing starts with 0
        end 
    else
        tmp1 = [];
    end
    
    % number of Cosine Basis functions
    if Opt.n_cos > 0
        for k = 1:Opt.n_cos
            tmp2(k) = strcat('cosine0', string(k-1));
        end
    else 
        tmp2 = [];
    end
    
    confound_names_tmp = [confound_names, tmp1, tmp2];
    confound_names = confound_names_tmp;
end