function [yp, e, A] = ar_n_bare_metal(y, n)
% AR(N) function with LS estimator for arrays of timeseries
% Largely derived from MATLAB's ar() function

% Args:
%    y (double): n_timepoints x n_channels
%    n (int): natural number
% Returns:
%    yp (double): size(y), predicted timeseries
%    e (double): size(y), prediction errors
%    A (double): 1 x n x n_channels, parameter estimates

    [Ncap, Ne] = size(y);
    R1 = zeros(n+1,n+1,Ne);
    jj = (n+1:Ncap);
    
    for kexp = 1:Ne
        yy = y(:,kexp);
        phi = zeros(length(jj),n);

        for k1 = 1:n
            phi(:,k1) = -yy(jj-k1);
        end

        R1_k = triu(qr([R1(:,:,kexp); [phi,yy(jj)]]));
        [nRr,nRc] = size(R1_k);
        R1(:,:,kexp) = R1_k(1:min(nRr,nRc),:);
    end
    
    R1 = -R1;
    covR = R1(1:n,1:n,:);
    
    if n > 1
        P = zeros(size(covR));
        for i = 1 : Ne
            P(:,:,i) = pinv(covR(:,:,i));
        end

        A = zeros(1, n, size(P,3));
        for i = 1 : Ne
          A(:,:,i) = (P(:,:,i) * R1(1:n, n+1, i)).';
        end
    else
        P = covR.^-1;
        A = P.*R1(1:n,n+1,:);
    end
    
    e = zeros(size(y)); yp = zeros(size(y));
    for kexp = 1:Ne
       tt = filter([1 A(:,:,kexp)], 1, y(:,kexp));
       tt(1:n) = zeros(n,1);
       e(:,kexp) = tt; 
       yp(:,kexp) = y(:,kexp) - tt;
    end

end