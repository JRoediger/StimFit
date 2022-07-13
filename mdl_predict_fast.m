function [meanpdf,ncomp] = mdl_predict_fast(SP_dirind,SP_ef,MDL_efvals,MDL_predsd,MDL_validind,MDL_pdf,MDL_pdfrange)

% Get overlap between valid models and sample directions as indices for
% both, Sample and Model matrices

[sp_bool,mdl_ind] = ismember(SP_dirind,MDL_validind);                                           % Check which components were calculated and get locations
mdl_ind(mdl_ind==0) = [];

% Get number of valid voxels
ncomp = sum(sp_bool);

% Get indices for mu in sample data (requires 50% of total function time)
[~,muind] = min(abs(squeeze(MDL_efvals(:,mdl_ind))-SP_ef(sp_bool)));                      % Get indices for sample efield on grid                                  % calculate differences in xvals (vector magnitudes) to find index with minimum difference

% Get standard deviations from model
sd = MDL_predsd(sub2ind(size(MDL_predsd),muind,mdl_ind));                                                                          % get standard deviation grid

% Get indices of std from all distributions (requires 25% of total function time)
[~,sdind] = min(abs(repmat(sd,size(MDL_pdfrange))-MDL_pdfrange));

indPDF = sub2ind([size(MDL_predsd,1),size(MDL_pdfrange,1)],muind,sdind);

PDF_all = MDL_pdf(indPDF,:); %Requires 20% of total function time

meanpdf = (mean(PDF_all,1))';     % Average PDFs


