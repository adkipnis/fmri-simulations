function basis = effects_of_interest(Opts, n)
    
    if Opts.pool_inference == false
        n = 1;
    elseif ~exist('n','var') && Opts.pool_inference == true
        error("Specify number of runs in F-Test!")
    end
    
    
    %%%
    
    basis_init = zeros(Opts.n_reg, Opts.n_reg-1);
    step = 1 + sum(Opts.hrf_derivs);
    
    basis_tmp = basis_init;
    for i = 1:step:Opts.n_stim_betas
        basis_tmp(i,i) = 1;
    end
    
    
    %%%
    
    if Opts.pool_inference == true
        basis_multi_tmp = repmat(basis_tmp, 1, n);
        basis = [basis_multi_tmp, zeros(Opts.n_reg, n)];
    else
        basis = [basis_tmp, zeros(Opts.n_reg, n)];
    end
    
end