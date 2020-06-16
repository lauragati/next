% computes the matrices I named "s" in front of fb, st, fa and st in
% equations A9 and A10 (e.g. materials25, appendix)
% THE MAIN CODE SMAT
function [s1, s2, s3, s4, s5] = smat(param,hx,knowTR,mpshock)
% this_code = mfilename;
% max_no_inputs = nargin(this_code);
% if nargin < max_no_inputs-1 
%     knowTR = 1; % knowing the TR is the default
%     mpshock = 1;% and having a monpol shock is the default.
% end

kapp = param.kapp;
bet  = param.bet;
alph = param.alph;
sig  = param.sig;
nx = length(hx);

s1 = [sig,1-bet,-sig*bet];
s2 = sig*[1,0,0]*(eye(nx)-bet*hx)^(-1);
s3 = [(1-alph)*bet,kapp*alph*bet,0];
s4 = [0,0,1]*(eye(nx)-alph*bet*hx)^(-1);
s5 = [0,0,0]; % if you don't wanna include a monpol shock.
if mpshock==1
    s5 = [0,1,0]; % if you wanna include a monpol shock.
end


psi_pi = param.psi_pi;
psi_x  = param.psi_x;
s1_TR = [sig-sig*bet*psi_pi, 1-bet-sig*bet*psi_x,0];

% Impose condition (*): input s1_old for s1, which really is the
% assumption of agents knowing the Taylor rule:
if knowTR==1
s1 = s1_TR;
end

% % When you call smat the first time, it should display what info assumption
% % about the Taylor rule you're using
% [ST] = dbstack;
% caller = ST(2).name;
% if strcmp(caller,'materials26')==1
%     star = sum(s1==s1_TR);
%     if star > 0
%         disp('Agents imposed to know Taylor rule')
%     else
%         disp('Agents DO NOT know Taylor rule')
%         
%     end
% end

