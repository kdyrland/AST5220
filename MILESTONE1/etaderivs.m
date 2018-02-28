function [ dnda ] = etaderivs( x )

load cmbparams.mat

a = exp(x);
dnda = c ./ (a.^2.*getH(x));


end

