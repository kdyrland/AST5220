function [ H ] = getH( x )

load cmbparams.mat

a = exp(x);
H = H0 .* sqrt((omega_b+omega_m).*a.^(-3) + (omega_r).*a.^(-4) + omega_v);

end

