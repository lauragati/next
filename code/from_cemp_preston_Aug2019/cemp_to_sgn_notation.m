function [fn, An, fl, Al, Gl, h, C, Ql, R, Gn, Qn, Qln] = cemp_to_sgn_notation(fk, fpibar, Apibar, fxi, Axi, Sxi, SIG, hpibar, H, R)
% 22 June 2019
k_l = size(fxi,1);

fn = [fk; fpibar]; 
An = [zeros(1,k_l);Apibar];
fl = fxi;
Al = Axi;
Gl = Sxi;
h = hpibar;
C = H';
Ql = SIG;
% R = R; R remains the same

% Fortunately these guys are zero for CEMP
% Gn = zeros(k_n); Qn = zeros(k_n); Qln = zeros(k_l,k_n);
Gn = 0; Qn = 0; Qln = 0;

