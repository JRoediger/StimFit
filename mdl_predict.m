function [meanpdf,ncomp] = mdl_predict(SP_dirind,SP_ef,MDL_efvals,MDL_predsd,MDL_validind,MDL_range)

% Get overlap between valid models and sample directions as indices for
% both, Sample and Model matrices
[sp_bool,mdl_ind] = ismember(SP_dirind,MDL_validind);                                           % Check which components were calculated and get locations
mdl_ind(mdl_ind==0) = [];

% Get number of valid voxels
ncomp = sum(sp_bool);

% Get indices for sample patient in MDL
[~,sp_loc] = min(abs(squeeze(MDL_efvals(:,mdl_ind))-SP_ef(sp_bool)));                      % Get indices for sample efield on grid                                  % calculate differences in xvals (vector magnitudes) to find index with minimum difference
mu = MDL_range(sp_loc);                                                                     % Convert to mean predicted values
                                                               
sd = MDL_predsd(sub2ind(size(MDL_predsd),sp_loc,mdl_ind));                                                                          % get standard deviation grid
pdf_tmp = nan(length(MDL_range),length(sd));

% Calculate probability density function over MDL_range for each voxel i

parfor i = 1:length(sd)
    pdf_tmp(:,i) = pdf(gmdistribution(mu(i),sd(i)),MDL_range);
end

meanpdf = mean(pdf_tmp,2);     % Average PDFs