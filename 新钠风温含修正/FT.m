function [ZZ] = FT(XX,YY,PP,JJ)

% RV RT J
% FV Na Nu Modify

x = XX;
y =YY;
j =JJ;
PxT = PP;

p00 =PxT(1,j);
p10 =PxT(2,j);
p01 =PxT(3,j);
p20 =PxT(4,j);
p11 =PxT(5,j);
p02 =PxT(6,j);
p30 =PxT(7,j);
p21 =PxT(8,j);
p12 =PxT(9,j);
p03 =PxT(10,j);
p40 =PxT(11,j);
p31 =PxT(12,j);
p22 =PxT(13,j);
p13 =PxT(14,j);
p04 =PxT(15,j);
p50 =PxT(16,j);
p41 =PxT(17,j);
p32 =PxT(18,j);
p23 =PxT(19,j);
p14 =PxT(20,j);
p05 =PxT(21,j);

FT = p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.^2 + p30.*x.^3 ...
    + p21.*x.^2.*y + p12.*x.*y.^2 + p03.*y.^3 + p40.*x.^4 + p31.*x.^3.*y ...
    + p22.*x.^2.*y.^2 + p13.*x.*y.^3 + p04.*y.^4 + p50.*x.^5 + p41.*x.^4.*y ...
    + p32.*x.^3.*y.^2 + p23.*x.^2.*y.^3 + p14.*x.*y.^4 + p05.*y.^5;

ZZ = FT;
end
