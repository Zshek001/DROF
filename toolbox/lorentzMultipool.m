function y = lorentzMultipool(p, x)
% Multi-pool Lorentzian function
% Huabin Zhang   Email: huabinz@connect.hku.hk
% Input
%   p: exchange pool parameter, p(1) for Z0, p(2-4) for water, p(5-7) for pool_1, ...
%   x: frequency offset, in ppm
% Output
%   y: Z-spectra

pn = (length(p)-1)/3;
y = p(1);
dw_0 = p(4); % B0 shift / water peak shift
for i = 1:pn
    A = p(3*i-1);
    G = p(3*i);
    dw = p(3*i+1);
    if i == 1
        y = y - A*G^2/4./(G^2/4+(x-dw_0).^2);
    else
        y = y - A*G^2/4./(G^2/4+(x-dw_0-dw).^2);
    end
end

end