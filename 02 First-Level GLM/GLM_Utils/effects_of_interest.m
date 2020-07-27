function basis = effects_of_interest(Opts)
    
    basis = zeros(Opts.n_reg, Opts.n_reg);
    step = 1 + sum(Opts.hrf_derivs);

    for i = 1:step:Opts.n_stim_betas
        basis(i,i) = 1;
    end
        
end