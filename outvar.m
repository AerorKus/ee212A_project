function nvar = outvar(num,den)
[r,p,K] = residue(num,den);
R = size(r,1);
R2 = size(K,1);
    if R2 > 1
        disp('Cannot continue...')
        return;
    end

    if R2 == 1
        nvar = K^2;
    else
        nvar = 0;
    end

for k = 1:R
    for m = 1:R
        integral = r(k)*conj(r(m))/(1-p(k)*conj(p(m)));
        nvar = nvar+integral;
    end
end