function [ dHp ] = getdHp( x )

load cmbparams.mat

top = -(omega_m + omega_b).*exp(-x) - 2 .* (omega_r).*exp(-2.*x) + 2.*omega_v.*exp(2.*x);
bottom = sqrt((omega_m + omega_b).*exp(-x) + (omega_r).*exp(-2.*x) + omega_v.*exp(2.*x));

dHp = (H0/2) .* (top./bottom);

end

