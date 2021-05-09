function [yp, e, A] = ar_n_bare_metal(y, n)
% AR(n) function with LS estimator for arrays of timeseries
% Largely derived from MATLAB's ar() function

% Args:
%    y (double): n_timepoints x n_channels
%    n (int): natural number (how many previous time points are used to
%    predict next time point)
% Returns:
%    yp (double): size(y), predicted timeseries
%    e (double): size(y), prediction errors
%    A (double): 1 x n x n_channels, parameter estimates

    [samples, channels] = size(y);
    R1 = zeros(n+1, n+1, channels);
    jj = (n+1:samples);
    
    for chan = 1:channels
        yy = y(:, chan); % time series at channel chan
        phi = zeros(length(jj), n); 
        
        % gather AR model predictors (channel values at k1 = 1, ..., n
        % time points before the regressed value)
        for k1 = 1:n
            phi(:,k1) = -yy(jj-k1); % design matrix
        end
        
        % Get upper triangular matrix R in A = Q*R,
        % where A is the design matrix concatenated with the criterium
        R1_k = triu(qr([R1(:,:,chan); [phi,yy(jj)]]));
        [nRr,nRc] = size(R1_k);
        R1(:,:,chan) = R1_k(1:min(nRr,nRc),:);
    end
    
    % get covariance estimates for the first n predictors
    R1 = -R1;
    covR = R1(1:n,1:n,:);
    
    
    if n > 1
        % get precision matrix for each channel
        P = zeros(size(covR));
        for chan = 1 : channels
            P(:,:,chan) = pinv(covR(:,:,chan));
        end
        
        % LS estimation of AR(n) parameters per channel
        A = zeros(1, n, size(P,3));
        for chan = 1 : channels
          A(:,:,chan) = (P(:,:,chan) * R1(1:n, n+1, chan)).';
        end
    else
        P = covR.^-1;
        A = P.*R1(1:n,n+1,:);
    end
    
    e = zeros(size(y));
    yp = zeros(size(y));
    for chan = 1:channels
       % get prediction errors using estimated AR(n) parameters for each
       % channel
       tt = filter([1 A(:,:,chan)], 1, y(:,chan)); 
       tt(1:n) = zeros(n,1);
       e(:,chan) = tt; 
       yp(:,chan) = y(:,chan) - tt;
    end
end