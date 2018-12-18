
%-Wiktor Olszowy
%-adjusted from spm_get_bf.m

dt      = 0.1
fMRI_T  = 16
[bf, p] = spm_hrf(dt, [], fMRI_T);

%-Canonical hemodynamic response function
dp      = 1;
p(6)    = p(6) + dp;
D       = (bf(:,1) - spm_hrf(dt,p,fMRI_T))/dp;
bf      = [bf D(:)];
p(6)    = p(6) - dp;

dp      = 0.01;
p(3)    = p(3) + dp;
D       = (bf(:,1) - spm_hrf(dt,p,fMRI_T))/dp;
bf      = [bf D(:)];
