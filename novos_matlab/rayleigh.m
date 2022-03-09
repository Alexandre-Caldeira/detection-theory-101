function [Y,Rm] = rayleigh(y,tj,M)
%rayleigh modified rayleigh test
%   TODO explain
    y = reshape(y(1:tj*M), tj, M);
    Y = fft(y);

    for i = 1:size(Y,2)
        Y(i) = sort(Y(i));
    end

    ang = unwrap(angle(Y));

    s = zeros(1,tj); 
    c = zeros(1,tj);

    for ii = 1:tj
        for jj = 1:M
            s(ii) = s(ii) +jj*sin(ang(ii,jj)); 
            c(ii) = c(ii) +jj*cos(ang(ii,jj));
        end
%         s(ii) = s(ii)/M;
%         c(ii) = c(ii)/M;
    end

    

    Rm = sqrt(c.^2 + s.^2)./(M^(3/2));
%     (M*sqrt(M));

end

