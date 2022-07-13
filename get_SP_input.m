function [SP_dirind,SP_ef] = get_SP_input(cdist,efmon,refvecs)


% Calculate efield according to current distribution
efx = cdist'*efmon;
efx(efx==0) = nan;
efx = reshape(efx,[],3);

% Get components
comp = efx*refvecs;
[mxtmp,SP_x] = max(comp,[],2);
SP_x(isnan(mxtmp)) = nan;
nvx = size(SP_x,1);

% Get direction indices
SP_dirind = sub2ind([size(refvecs,2),nvx],SP_x',1:nvx);

% Calculate normvectors
SP_ef = (vecnorm(efx,2,2))';