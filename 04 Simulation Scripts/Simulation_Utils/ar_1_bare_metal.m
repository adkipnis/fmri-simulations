function [y_pred, phi] = ar_1_bare_metal(y)
    T = size(y,1);
    y0 = y(1,:);
    y_temp = y(2:end-1, :);
    q = movprod(y_temp,2);
    q_m = mean(q);
    p = y_temp.^2;
    p_m = mean(p);
    phi = q_m./p_m;
    y_pred = zeros(T,size(y, 2));
    y_pred(1,:) = y0;
    for i = 2:T
        y_pred(i,:) = y(i-1,:).*phi;
    end
end