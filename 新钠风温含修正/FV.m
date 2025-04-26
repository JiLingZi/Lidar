function [ZZ] = FV(XX,YY,PP,JJ)

% RV RT J
% FV Na Nu Modify

x = XX;
y =YY;
j =JJ;
PxV = PP;

p00 =PxV(1,j);
p10 =PxV(2,j);
p01 =PxV(3,j);
p20 =PxV(4,j);
p11 =PxV(5,j);
p02 =PxV(6,j);
p30 =PxV(7,j);
p21 =PxV(8,j);
p12 =PxV(9,j);
p03 =PxV(10,j);
p40 =PxV(11,j);
p31 =PxV(12,j);
p22 =PxV(13,j);
p13 =PxV(14,j);
p04 =PxV(15,j);
p50 =PxV(16,j);
p41 =PxV(17,j);
p32 =PxV(18,j);
p23 =PxV(19,j);
p14 =PxV(20,j);
p05 =PxV(21,j);

FV = p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.^2 + p30.*x.^3 ...
    + p21.*x.^2.*y + p12.*x.*y.^2 + p03.*y.^3 + p40.*x.^4 + p31.*x.^3.*y ...
    + p22.*x.^2.*y.^2 + p13.*x.*y.^3 + p04.*y.^4 + p50.*x.^5 + p41.*x.^4.*y ...
    + p32.*x.^3.*y.^2 + p23.*x.^2.*y.^3 + p14.*x.*y.^4 + p05.*y.^5;

ZZ = FV;
end
