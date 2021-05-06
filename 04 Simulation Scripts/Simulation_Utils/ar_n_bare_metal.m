function [yp, e, A] = ar_n_bare_metal(y, n, average)
% AR(N) function with LS estimator for arrays of timeseries
% Largely derived from MATLAB's ar() function

% Args:
%    y (double): n_timepoints x n_channels
%    n (int): natural number
% Returns:
%    yp (double): size(y), predicted timeseries
%    e (double): size(y), prediction errors
%    A (double): 1 x n x n_channels, parameter estimates

    if ~exist('average','var') || isempty(average)
        average = false;
    end

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
        if average
            covR = mean(covR,3);
            R1 = mean(R1,3);
            P = pinv(covR);
            A = repmat((P * R1(1:n,n+1)).', 1,1,size(y,2));
        else
            P = zeros(size(covR));
            A = zeros(1, n, size(P,3));
            for i = 1 : Ne
                P(:,:,i) = pinv(covR(:,:,i));
                A(:,:,i) = (P(:,:,i) * R1(1:n, n+1, i)).';
            end
        end
    else
        if average
            P = mean(covR,3).^-1;
            A = repmat((P * R1(1:n,n+1)), 1,1,size(y,2));
        else
            P = covR.^-1;
            A = P.*R1(1:n,n+1,:);
        end
    end
    
    e = zeros(size(y)); yp = zeros(size(y));
    for kexp = 1:Ne
       tt = filter([1 A(:,:,kexp)], 1, y(:,kexp));
       e(:,kexp) = tt;
       tt(1:n) = zeros(n,1);
       yp(:,kexp) = y(:,kexp) - tt;
    end

end