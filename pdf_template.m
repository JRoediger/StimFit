function MDL_pdf = pdf_template(range_mu,range_sd)

range_mu = range_mu(:);
range_sd = range_sd(:);

MDL_pdf = nan(length(range_mu),length(range_sd),length(range_mu));


for i_mu = 1:length(range_mu)
    for i_sd = 1:length(range_sd)
        MDL_pdf(i_mu,i_sd,:) = pdf(gmdistribution(range_mu(i_mu),range_sd(i_sd)),range_mu);
    end
end

