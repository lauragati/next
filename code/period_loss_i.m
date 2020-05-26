function loss_scalar = period_loss_i(param,hx,i,fa,fb,ssim,L)
lamx = param.lamx;
sig  = param.sig;
kapp = param.kapp;
[s1, s2, s3, s4] = smat(param,hx);


s3fa = s3*fa;
s1fb = s1*fb;
s4s  = s4*ssim;
s2s  = s2*ssim;

loss =(-1).*L+kapp.^2.*s1fb.^2+lamx.*s1fb.^2+2.*kapp.^2.*s1fb.*s2s+2.* ...
  lamx.*s1fb.*s2s+kapp.^2.*s2s.^2+lamx.*s2s.^2+2.*kapp.*s1fb.*s3fa+ ...
  2.*kapp.*s2s.*s3fa+s3fa.^2+2.*kapp.*s1fb.*s4s+2.*kapp.*s2s.*s4s+ ...
  2.*s3fa.*s4s+s4s.^2+i.*((-2).*kapp.^2.*s1fb.*sig+(-2).*lamx.* ...
  s1fb.*sig+(-2).*kapp.^2.*s2s.*sig+(-2).*lamx.*s2s.*sig+(-2).* ...
  kapp.*s3fa.*sig+(-2).*kapp.*s4s.*sig)+i.^2.*(kapp.^2.*sig.^2+ ...
  lamx.*sig.^2);

loss_scalar = max(max(abs(loss)));
