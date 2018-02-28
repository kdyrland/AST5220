function [ Hp ] = getHp( x )

load cmbparams.mat

a = exp(x);
Hp = a .* getH(x);

end

