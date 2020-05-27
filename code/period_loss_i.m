function loss = period_loss_i(param,hx,i,fa,fb,ssim,L)
lamx = param.lamx;
sig  = param.sig;
kapp = param.kapp;
[s1, s2, s3, s4] = smat(param,hx);

fas = s3*fa+s4*ssim;
fbs = s1*fb +s2*ssim;

% Notes 27 May 2020, double-checked with Mathematica, materials31.nb
loss = (kapp^2+lamx)*sig^2.*i.^2 -(2*sig*(kapp^2 +lamx)*fbs +2*sig*kapp*fas).*i ...
    +2*kapp*fas.*fbs + fas.^2 + (kapp^2+lamx)*fbs.^2 -L; 
