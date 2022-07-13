function stimfit_main(options)

% StimFit Main function launched by StimFit UI

%% Check if selected folder contains patients with electrode reconstructions

list = dir(options.folder);
elefound = zeros(length(list),1);
singlefound = 0;

for i = 3:length(list)
    elefound(i) = exist([list(i).folder,filesep,list(i).name,filesep,'ea_reconstruction.mat'],'file')==2;
end
patindex = find(elefound);
npat = sum(elefound);

if npat > 1
    ea_dispt([num2str(npat), ' electrode reconstructions found.'])
elseif npat == 1
    ea_dispt('1 electrode reconstruction found.')
else
    % single folder was selected
    singlefound = exist([options.folder,filesep,'ea_reconstruction.mat'],'file')==2;
    if singlefound
        npat = 1;
        ea_dispt('1 electrode reconstruction found.')
    else
        error('No electrode reconstructions found... Please reconstruct electrodes first!')
    end
end

%% Prepare Model
ea_dispt('Load and prepare model...')

load(options.modelpth)
[~,options.mdlname] = fileparts(options.modelpth);

if options.fastpredict
    range_mu = model.raMDL_range;
    rr = log([min(model.raMDL_predsd,[],'all'),max(model.raMDL_predsd,[],'all')]);
    model.raPDFrange = exp(rr(1):diff(rr)/(options.nbins-1):rr(2));
    model.raMDL_pdf = pdf_template(range_mu,model.raPDFrange);
    rr = log([min(model.trMDL_predsd,[],'all'),max(model.trMDL_predsd,[],'all')]);
    model.trPDFrange = exp(rr(1):diff(rr)/(options.nbins-1):rr(2));
    model.trMDL_pdf = pdf_template(range_mu,model.trPDFrange);
    range_mu = model.sMDL_range;
    rr = log([min(model.sMDL_predsd,[],'all'),max(model.sMDL_predsd,[],'all')]);
    model.sPDFrange = exp(rr(1):diff(rr)/(options.nbins-1):rr(2));
    model.sMDL_pdf = pdf_template(range_mu,model.sPDFrange);
else
    if options.calcoptdistr
        warning(['You are about to start the optimizer without the "fastpredict" option... This might result in long execution times and will be stopped after ',num2str(options.optimizer.MaxTime/60),' minutes.'])
    end
    
    model.raMDL_pdf = nan;
    model.raPDFrange = nan;
    model.trMDL_pdf = nan;
    model.trPDFrange = nan;
    model.sMDL_pdf = nan;
    model.sPDFrange = nan;
end


%% Prepare parallel computing

if options.optimizer.nCores>1 && strcmp(options.optimizer.Solver,'mult')
    pp = gcp('nocreate');
    if isempty(pp)
        poolsize = 0;
    else
        poolsize = pp.NumWorkers;
    end
    if poolsize>1 && options.optimizer.nCores <=1
        warning('Found an interactive session. Cancelling existing session to recruit the specified number of workers...')
        delete(pp)
    end
    
    if options.optimizer.nCores~=poolsize
        try
            parpool(options.optimizer.nCores)
        catch
            initnewpar(options.optimizer.nCores)
        end
    end
end


%% Process patients

for i = 1:npat
    if singlefound
        patientfolder = options.folder;
        options.figname = 'Live Tracking';
    else
        patientfolder = [list(patindex(i)).folder,filesep,list(patindex(i)).name];
        ea_dispt(['Analyzing patient ', list(patindex(i)).name])
        options.figname = ['Live Tracking - ',list(patindex(i)).name];
    end
    
    % Generate E-field templates if necessary
    if ~exist([patientfolder,filesep,'stimulations',filesep,model.settings.spacedef,filesep,model.settings.emodel],'dir')==7
        ea_dispt('No suited E-field templates found... Generating E-fields first... this may take a while.')
        try
            generate_EF(patientfolder,options)
        catch
            warning('Generating E-field templates failed!')
            continue
        end
    end
    
    try
        [optimresults(i).eff_comb,optimresults(i).eff_ra,optimresults(i).eff_tr,optimresults(i).eff_se,optimresults(i).ncomp,optimresults(i).perc,optimresults(i).current,...
            monrevresults(i).eff_comb, monrevresults(i).eff_ra,monrevresults(i).eff_tr,monrevresults(i).eff_se,monrevresults(i).ncomp,monrevresults(i).coords,options] = stimpredict(patientfolder,model,options,options.tremorweight);
    catch
        warning('Prediction failed!')
        continue
    end
    
    if options.export
        outpth = [patientfolder,filesep,'predictions'];
        stimpredict_export(optimresults(i),monrevresults(i),options,outpth)
    end
    %     if options.vizz
    %         stimpredict_vizz(monrevresults(i),options,model,outpth);
    %     end
    
end


ea_dispt('Done')

end


%% Functions

function [eff_comb,eff_ra,eff_tr,eff_se,eff_ncomp,finalperc,totalcurrent,monrev_comb,monrev_ra,monrev_tr,monrev_se,monrev_ncomp,contcoords,options] = stimpredict(patientfolder,model,options,tremorweight)

% Loads efield templates and executes predictions of monopolar review
% and/or optimal setting


%% Initialize variables
refvecs = model.refvecs;

raMDL_efvals = model.raMDL_efvals;
raMDL_predsd = model.raMDL_predsd;
raMDL_validind = model.raMDL_validind;
raMDL_range = model.raMDL_range';
raMDL_pdf = reshape(model.raMDL_pdf,[],size(model.raMDL_pdf,3));
raMDL_pdfrange = model.raPDFrange';

trMDL_efvals = model.trMDL_efvals;
trMDL_predsd = model.trMDL_predsd;
trMDL_validind = model.trMDL_validind;
trMDL_range = model.trMDL_range';
trMDL_pdf = reshape(model.trMDL_pdf,[],size(model.trMDL_pdf,3));
trMDL_pdfrange = model.trPDFrange';

sMDL_efvals = model.sMDL_efvals;
sMDL_predsd = model.sMDL_predsd;
sMDL_validind = model.sMDL_validind;
sMDL_range = model.sMDL_range';
sMDL_pdf = reshape(model.sMDL_pdf,[],size(model.sMDL_pdf,3));
sMDL_pdfrange = model.sPDFrange';

fastpredict = options.fastpredict;

eff_comb = nan;
eff_ra = nan;
eff_tr = nan;
eff_se = nan;
eff_ncomp = nan;
totalcurrent = nan;
monrev_comb = nan;
monrev_ra = nan;
monrev_tr = nan;
monrev_se = nan;
monrev_ncomp = nan;
contcoords = nan;

%% Create predictions folder

savepth = [patientfolder,filesep,'predictions',filesep,options.mdlname];
if ~exist([patientfolder,filesep,'predictions'],'dir')
    mkdir([patientfolder,filesep,'predictions'])
end
if ~exist(savepth,'dir')
    mkdir(savepth)
end


%% Load patient data

sidename = {'right','left'};

% Check Version
vatname = ['c00_a',sprintf('%02d',model.settings.maxamp*10),'_vat_evec',model.settings.dims{1},'_',sidename{1},'.nii'];
try
    stimpth = [patientfolder,filesep,'stimulations',filesep,'MNI_ICBM_2009b_NLIN_ASYM',filesep,model.settings.spacedef,filesep,model.settings.emodel];
    spm_vol([stimpth,filesep, vatname]);
catch
    stimpth = [patientfolder,filesep,'stimulations',filesep,model.settings.spacedef,filesep,model.settings.emodel];
    spm_vol([stimpth,filesep, vatname]);
end

% Determine number of contacts
load([patientfolder,filesep,'ea_reconstruction.mat'],'reco')
nconts = length(reco.mni.coords_mm{1,1});

ea_dispt('Loading efields...')
notfound = [];
if exist(stimpth,'dir')
    for side = 1:2
        for c = 1:nconts
            %             try
            for dim = 1:length(model.settings.dims)
                
                % Load efields (remember to switch PL and PM contacts
                % on left side)
                if strcmp(sidename{side},'left') && contains(reco.props(side).elmodel,'Directed')
                    if  c==3 || c==6
                        vatname = ['c',sprintf('%02d',c),'_a',sprintf('%02d',model.settings.maxamp*10),'_vat_evec',model.settings.dims{dim},'_',sidename{side},'.nii'];
                    elseif c==4 || c==7
                        vatname = ['c',sprintf('%02d',c-2),'_a',sprintf('%02d',model.settings.maxamp*10),'_vat_evec',model.settings.dims{dim},'_',sidename{side},'.nii'];
                    else
                        vatname = ['c',sprintf('%02d',c-1),'_a',sprintf('%02d',model.settings.maxamp*10),'_vat_evec',model.settings.dims{dim},'_',sidename{side},'.nii'];
                    end
                else
                    vatname = ['c',sprintf('%02d',c-1),'_a',sprintf('%02d',model.settings.maxamp*10),'_vat_evec',model.settings.dims{dim},'_',sidename{side},'.nii'];
                end
                
                V = spm_vol([stimpth,filesep, vatname]);
                Vimg = spm_read_vols(V);
                Vimg = Vimg(1:model.settings.sf:end,1:model.settings.sf:end,1:model.settings.sf:end);
                
                if dim == 1 && side==2
                    Vimg = Vimg*-1;
                end
                ef(side,c,:,dim) = Vimg(:);
            end
            %             catch
            %                 notfound = [notfound; c,model.settings.maxamp,side];
            %             end
        end
        
    end
else
    notfound = 1; % generate stims
end

if notfound == 1
    ea_dispt('No efield simulation folder found... continuing generating efields. This may take a while...')
elseif size(notfound,1) > 1
    ea_dispt([num2str(size(notfound,1)), '/', num2str(nconts*2), ' simulations not found... continuing generating efields. This may take a while...'])
end

ef = ef/5;
efsample = ef(:,:,model.vxsample,:);
ind = isnan(squeeze(sum(sum(sum(efsample,1),2),4)));
efsample = efsample(:,:,~ind,:);

nsides = size(efsample,1);
nconts = size(efsample,2);
options.nconts = nconts;
options.nsides = nsides;
options.ele = reco.props(1).elmodel;

if contains(options.ele,'Directed')
    options.dir = 1;
else
    options.dir = 0;
end

if options.vizz
    [fax,txt,fig_cont,fig_main] = openvizzfig(options.calcmonrev,options.calcoptdistr,options,sidename);
else
    fax = [];
    txt = [];
    fig_cont = [];
    fig_main = [];
end

% Get monopolar review data
monrevfound = 0;
if exist([savepth,filesep,'monrev.mat'],'file')==2
    load([savepth,filesep,'monrev'])
    if exist('amplitudes','var')
        if isequal(options.ampl,amplitudes)
            monrevfound = 1;
            disp('Monopolar reviews found in patient folder. They will serve as starting points for the optimizer.')
            if options.vizz
                txt.String = 'Monopolar reviews found in patient folder. They will serve as starting points for the optimizer.';
                
                for side = 1:2
                    Cdata = squeeze(monrev_comb(:,:,side));
                    imagesc(fax{side},Cdata)
                    caxis(fax{side},[0,max(Cdata,[],'all')])
                    set(fax{side},'YDir','normal')
                    title(fax{side},['Monopolar Settings (', sidename{side},')'])
                    xlabel(fax{side},'amplitude [mA]')
                    ylabel(fax{side},'contact')
                    
                    [~,maxc] = max(mean(Cdata,2));
                    contfig = fig_main;
                    tmp2 = squeeze(contfig(:,:,2));
                    tmp2(fig_cont{maxc}) = 0;
                    tmp3 = squeeze(contfig(:,:,3));
                    tmp3(fig_cont{maxc}) = 0;
                    contfig = cat(3,contfig(:,:,1),tmp2,tmp3);
                    axes(fax{side+2});
                    imshow(contfig)
                    drawnow
                end
                
            end
            
        else
            warning('Monopolar review settings differ from the once found in path... Continue to generate and overwrite new monopolar reviews')
            keyboard
        end
    end
end


%% Predict Monopolar Review

% Sollte noch parallelisiert werden
if options.calcmonrev && ~monrevfound

    ea_dispt('')
    pb = CmdLineProgressBar('Predicting monopolar review... ');
    
    amplitudes = options.ampl;
    
    nampl = length(amplitudes);
    
    % Preallocate
    monrev_comb = nan(nconts,nampl,nsides);
    monrev_ra = nan(nconts,nampl,nsides);
    monrev_tr = nan(nconts,nampl,nsides);
    monrev_se = nan(nconts,nampl,nsides);
    monrev_ncomp = nan(3,nconts,nampl,nsides);
    contcoords = nan(nconts,3,2);
    totit = nsides*nconts*nampl;
    it = 1;
    
    for side = 1:nsides
        efs = squeeze(efsample(side,:,:,:));
        efs = reshape(efs,nconts,[]);
        
        if side ==2 && options.dir
            efs(3:4,:) = flipud(efs(3:4,:));                     % Switch indices for PL and PM contacts on left side
            efs(6:7,:) = flipud(efs(6:7,:));
        end
        
        for ai = 1:nampl
            for ci = 1:nconts
                
                cdist = zeros(nconts,1);
                cdist(ci) = options.ampl(ai);
                [monrev_comb(ci,ai,side),monrev_ra(ci,ai,side),monrev_tr(ci,ai,side),monrev_se(ci,ai,side),monrev_ncomp(:,ci,ai,side)]...
                    = checkscore(cdist,efs,refvecs,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_range,raMDL_pdf,raMDL_pdfrange,...
                    trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_range,trMDL_pdf,trMDL_pdfrange,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_range,sMDL_pdf,sMDL_pdfrange,tremorweight,fastpredict);
                pb.print(it,totit)
                
                if options.vizz
                    Cdata = squeeze(monrev_comb(:,:,side));
                    imagesc(fax{side},Cdata)
                    caxis(fax{side},[0,max(Cdata,[],'all')])
                    set(fax{side},'YDir','normal')
                    title(fax{side},['Monopolar Settings (', sidename{side},')'])
                    xlabel(fax{side},'amplitude [mA]')
                    ylabel(fax{side},'contact')
                    
                    contfig = fig_main;
                    tmp2 = squeeze(contfig(:,:,2));
                    tmp2(fig_cont{ci}) = 0;
                    tmp3 = squeeze(contfig(:,:,3));
                    tmp3(fig_cont{ci}) = 0;
                    contfig = cat(3,contfig(:,:,1),tmp2,tmp3);
                    axes(fax{side+2});
                    imshow(contfig)
                    
                    txt.String = ['Predicted monopolar review: ',num2str(it),'/',num2str(totit)];
                    pause(0.01)
                end
                
                it = it+1;
            end
        end
        
        if options.vizz
            [~,maxc] = max(mean(Cdata,2));
            contfig = fig_main;
            tmp2 = squeeze(contfig(:,:,2));
            tmp2(fig_cont{maxc}) = 0;
            tmp3 = squeeze(contfig(:,:,3));
            tmp3(fig_cont{maxc}) = 0;
            contfig = cat(3,contfig(:,:,1),tmp2,tmp3);
            axes(fax{side+2});
            imshow(contfig)
            pause(0.01)
        end
        
        if side == 2
            contcoords(:,:,side) = ea_flip_lr_nonlinear(reco.mni.coords_mm{side});
        else
            contcoords(:,:,side) = reco.mni.coords_mm{side};
        end
    end
    toc
    save([savepth,filesep,'monrev'],'monrev_comb','monrev_ra','monrev_tr','monrev_se','monrev_ncomp','amplitudes')
end

%% Run Optimization function

finalperc = nan(nconts,options.optimizer.nIter,nsides);

if options.calcoptdistr

    trackglob = options.vizz | options.optimizer.Track;
    
    % Define optimization parameters
    lb = zeros(1,nconts);                                               % Minimum current per contact
    ub = ones(1,nconts)*options.optimizer.MaxCurrCont;                % Maximum current per contact
    A = ones(1,nconts);
    b = options.optimizer.MaxCurrTot;                                 % Maximum total current
    %     Aeq = [];
    %     beq = [];
    
    optimizeroptions = optimoptions('fmincon','Display', options.optimizer.Display.local,'CheckGradients',options.optimizer.CheckGradients,...
        'OptimalityTolerance',options.optimizer.OptimalityTolerance,'DiffMinChange',options.optimizer.DiffMinChange,'MaxIterations',options.optimizer.MaxIterations,...
        'PlotFcn',options.optimizer.PlotFcn,'StepTolerance',options.optimizer.StepTolerance,'ConstraintTolerance',options.optimizer.ConstraintTolerance);
    
    % Preallocate
    finalcurrent = nan(nconts,options.optimizer.nIter,2);
    eff_comb = nan(options.optimizer.nIter,2);
    eff_ra = nan(options.optimizer.nIter,2);
    eff_tr = nan(options.optimizer.nIter,2);
    eff_se = nan(options.optimizer.nIter,2);
    eff_ncomp = nan(3,options.optimizer.nIter,2);
    %     allcurs = nan(nconts,options.optimizer.Rerun,nconts,options.optimizer.nIter,2);
    %     allfvals = nan(nconts,options.optimizer.Rerun,options.optimizer.nIter,2);
    %     allorigins = nan(nconts,nconts,options.optimizer.nIter,2);
    
    it = 1;
    
    checkval = cell(4,nsides);
    
    for side = 1:nsides
        efs = squeeze(efsample(side,:,:,:));
        efs = reshape(efs,nconts,[]);
        
        if side ==2 && options.dir                      % Switch indices for PL and PM contacts on left side
            efs(3:4,:) = flipud(efs(3:4,:));
            efs(6:7,:) = flipud(efs(6:7,:));
        end
        
        if trackglob
            global checkit_mot checkit_se checkit_distmot checkit_distse counter optcounter
            checkit_mot = [];
            checkit_se = [];
            checkit_distmot = [];
            checkit_distse = [];
            counter = 1;
            optcounter = 0;
        end
        
        for i = 1:options.optimizer.nIter
            ea_dispt(['Running gradient descent: ',num2str(it),'/',num2str(nsides*options.optimizer.nIter)])
            
            maxsideeff = options.maxse*(1-options.optimizer.PercDiff/100)^(i-1);
            
            % Create starting points
            if exist('monrev_se','var')
                % Get monopolar settings as starting points for multi
                % and single start
                
                combtmp = monrev_comb;
                combtmp(monrev_se>options.maxse) = nan;
                [~,start_ind] = max(combtmp,[],2);
                
                startmulti = diag(ones(nconts,1)).*options.ampl(start_ind(:,:,side));
            else
                % Create artificial starting points
                if nconts==8        % assume directional leads
                    startmulti = ([diag(ones(8,1))*3;[0,1,1,1,0,0,0,0];[0,0,0,0,1,1,1,0]])*2.5/3;     % set 10 starting points with 2.5mA at each monopolar setting and for circular stimulations
                else
                    startmulti = diag(ones(nconts,1))*2.5;                                             % set n starting points with 2.5mA at each monopolar setting
                end
            end
            
            cdist = startmulti(end,:)';
            
            % Define ojective function
            f = @(cdist)generatescore(cdist,efs,refvecs,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_range,raMDL_pdf,raMDL_pdfrange,...
                trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_range,trMDL_pdf,trMDL_pdfrange,tremorweight,fastpredict,trackglob);
            
            % Define nonlinear constraint function for side-effects
            nonlcon = @(cdist)sideeffconstr(cdist,efs,refvecs,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_range,sMDL_pdf,sMDL_pdfrange,maxsideeff,fastpredict,trackglob,options.vizz,fax,fig_cont,fig_main,side);
            
            % Create optimization problem
            problem = createOptimProblem('fmincon','x0',cdist,'objective',f,'lb',lb,'ub',ub,'Aineq',A,'bineq',b,'nonlcon',nonlcon,'options',optimizeroptions);
            
            if strcmp('single',options.optimizer.Solver)
                % Run single gradiest descent
                
                current = fmincon(problem);
                
            elseif strcmp('global',options.optimizer.Solver)
                % Run global optimizer (does not support parallel computing
                % and might therefore lead to long execution times)
                
                % HAVENT PLAYED AROUND MUCH WITH THE SETTINGS YET SO MAYBE ROOM FOR IMPROVEMENT!
                
                gs = GlobalSearch('Display',options.optimizer.Display.global,'MaxTime',options.optimizer.MaxTime);
                
                [current,fval,exitflag,output,solutions] = run(gs,problem);
                
            elseif strcmp('mult',options.optimizer.Solver)
                
                % Run multi-start optimizer (default)
                startmulti = startmulti+10^-6;                          % add a little bit of current to each contact because 0 needs to be excluded
                tpoints = CustomStartPointSet(startmulti);
                
                ms =  MultiStart('UseParallel', options.optimizer.nCores>1,'MaxTime',options.optimizer.MaxTime,'Display',options.optimizer.Display.mult);

                [~,fval,~,~,solutions] = run(ms,problem,tpoints);

                % get best rounded solution
                allcomb = nan(size(solutions,2),1);
                allside = nan(size(solutions,2),1);
                for si = 1:size(solutions,2)
                    [perctmp,curtmp] = roundnxclude(solutions(si).X,options.optimizer.MinPercCont);
                    disttmp = perctmp*curtmp/100;
                    [allcomb(si),~,~,allside(si)] = checkscore(disttmp,efs,refvecs,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_range,raMDL_pdf,raMDL_pdfrange,...
                        trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_range,trMDL_pdf,trMDL_pdfrange,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_range,sMDL_pdf,sMDL_pdfrange,tremorweight,fastpredict);
                end
                alldiff = [(fval+allcomb)*1000,maxsideeff-allside];
                [rankscore,rankind] =  sort(sum(alldiff,2),'descend');
                
                if rankscore(1)<-5
                    warning('Rounding currents to adequate percentages changed the predicted results')
                end
                
                current = solutions(rankind(1)).X;
            end
            
            [finalperc(:,i,side),totalcurrent(i,side)] = roundnxclude(current,options.optimizer.MinPercCont);
            finalcurrent(:,i,side) = finalperc(:,i,side)*totalcurrent(i,side)/100;
            [eff_comb(i,side),eff_ra(i,side),eff_tr(i,side),eff_se(i,side),eff_ncomp(:,i,side)] = checkscore(finalcurrent(:,i,side),efs,refvecs,raMDL_efvals,...
                raMDL_predsd,raMDL_validind,raMDL_range,raMDL_pdf,raMDL_pdfrange,trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_range,trMDL_pdf,trMDL_pdfrange,...
                sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_range,sMDL_pdf,sMDL_pdfrange,tremorweight,fastpredict);
            it = it+1;
            
        end
        if options.optimizer.Track
            checkval{1,side} = checkit_mot;
            checkval{2,side} = checkit_se;
            checkval{3,side} = checkit_distmot;
            checkval{4,side} = checkit_distse;
            checkval{5,side} = startmulti;
            save('checkval_test2_2','checkval')
        end
        toc
    end
end

% finalcurrent = squeeze(finalcurrent);
% save('optres_testdata_ra','optimresults')

end

%% Functions

function eff_comb = generatescore(cdist,efmon,refvecs,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_range,raMDL_pdf,raMDL_pdfrange,...
    trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_range,trMDL_pdf,trMDL_pdfrange,tremorweight,fastpredict,track)

[SP_dirind,SP_ef] = get_SP_input(cdist,efmon,refvecs);

if fastpredict
    [meanpdf_ra,~] = mdl_predict_fast(SP_dirind,SP_ef,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_pdf,raMDL_pdfrange);
    [meanpdf_tr,~] = mdl_predict_fast(SP_dirind,SP_ef,trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_pdf,trMDL_pdfrange);
else
    [meanpdf_ra,~] = mdl_predict(SP_dirind,SP_ef,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_range);
    [meanpdf_tr,~] = mdl_predict(SP_dirind,SP_ef,trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_range);
end

% Get results from PDF and combine predictions
[~,maxind] = max(meanpdf_ra,[],1);
eff_ra = raMDL_range(maxind);
[~,maxind] = max(meanpdf_tr,[],1);
eff_tr = trMDL_range(maxind);

eff_comb = eff_ra*(1-tremorweight) + eff_tr*tremorweight;

if track
    global checkit_mot counter checkit_distmot
    checkit_mot(counter) = eff_comb;
    checkit_distmot(counter,:) = cdist;
    counter = counter+1;
end

eff_comb = -eff_comb;

end

function [eff_comb,eff_ra,eff_tr,eff_se,ncomp] = checkscore(cdist,efmon,refvecs,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_range,raMDL_pdf,raMDL_pdfrange,...
    trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_range,trMDL_pdf,trMDL_pdfrange,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_range,sMDL_pdf,sMDL_pdfrange,tremorweight,fastpredict)

[SP_dirind,SP_ef] = get_SP_input(cdist,efmon,refvecs);

if fastpredict
    [meanpdf_ra,ncomp(1)] = mdl_predict_fast(SP_dirind,SP_ef,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_pdf,raMDL_pdfrange);
    [meanpdf_tr,ncomp(2)] = mdl_predict_fast(SP_dirind,SP_ef,trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_pdf,trMDL_pdfrange);
    [meanpdf_se,ncomp(3)] = mdl_predict_fast(SP_dirind,SP_ef,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_pdf,sMDL_pdfrange);
else
    [meanpdf_ra,ncomp(1)] = mdl_predict(SP_dirind,SP_ef,raMDL_efvals,raMDL_predsd,raMDL_validind,raMDL_range);
    [meanpdf_tr,ncomp(2)] = mdl_predict(SP_dirind,SP_ef,trMDL_efvals,trMDL_predsd,trMDL_validind,trMDL_range);
    [meanpdf_se,ncomp(3)] = mdl_predict(SP_dirind,SP_ef,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_range);
end

% Normalize probabilities
meanpdf_ra = meanpdf_ra/100;
meanpdf_ra = meanpdf_ra/sum(meanpdf_ra);
meanpdf_tr = meanpdf_tr/100;
meanpdf_tr = meanpdf_tr/sum(meanpdf_tr);
meanpdf_se = meanpdf_se/10;
meanpdf_se = meanpdf_se/sum(meanpdf_se);

[~,maxind] = max(meanpdf_ra,[],1);
eff_ra = raMDL_range(maxind);
[~,maxind] = max(meanpdf_tr,[],1);
eff_tr = trMDL_range(maxind);

eff_comb = eff_ra*(1-tremorweight) + eff_tr*tremorweight;
eff_se = logit2prob(sMDL_range'*meanpdf_se);

end

function [c,ceq] = sideeffconstr(cdist,efmon,refvecs,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_range,sMDL_pdf,sMDL_pdfrange,maxsideeff,fastpredict,track,vizz,fax,fig_cont,fig_main,side)


[SP_dirind,SP_ef] = get_SP_input(cdist,efmon,refvecs);

if fastpredict
    [meanpdf_se,~] = mdl_predict_fast(SP_dirind,SP_ef,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_pdf,sMDL_pdfrange);
else
    [meanpdf_se,~] = mdl_predict(SP_dirind,SP_ef,sMDL_efvals,sMDL_predsd,sMDL_validind,sMDL_range);
end

% Normalize probabilities
meanpdf_se = meanpdf_se/10;
meanpdf_se = meanpdf_se/sum(meanpdf_se);

sideff = logit2prob(sMDL_range'*meanpdf_se);

if track
    global checkit_se checkit_distse counter
    checkit_se(counter-1) = sideff;
    checkit_distse(counter-1,:) = cdist;
    
    % Start visualization after parameterspace exploration
    if vizz && counter>10
        global optcounter checkit_mot
        
        % Check if current setting is better then optimal setting
        if ~optcounter
            optmot = nan;
            optdist = zeros(1,size(checkit_distse,2));
        else
            optmot = checkit_mot(optcounter)*100;
            optdist = checkit_distse(optcounter,:);
        end
        
        if checkit_se(counter-1)>maxsideeff*1.05
            currentmot = nan;
        else
            currentmot = checkit_mot(counter-1)*100;
        end
        
        if optmot < currentmot
            optcounter = counter-1;
            optdist = checkit_distse(optcounter,:);
            optmot = checkit_mot(optcounter)*100;
        elseif isnan(optmot) && ~isnan(currentmot)
            optcounter = counter-1;
            optdist = checkit_distse(optcounter,:);
            optmot = checkit_mot(optcounter)*100;
            addpoints(fax{side+6},counter-10,optmot);
            ylim(fax{side+4},optmot+[-0.1 0.1])
        end
        
        % Update figure (plot)
        addpoints(fax{side+6},counter-9,optmot);
        xticks(fax{side+4},'auto')
        yticks(fax{side+4},'auto')
        xlim(fax{side+4},[1 counter-9])
        limtmp= ylim(fax{side+4});
        limtmp(2) = optmot+0.1;
        ylim(fax{side+4},limtmp)
        
        % Update figure (contacts)
        totcurr = sum(optdist);
        contfig = fig_main;
        tmp2 = squeeze(contfig(:,:,2));
        tmp3 = squeeze(contfig(:,:,3));
        
        for ci = 1:length(optdist)
            ciperc = optdist(ci)/totcurr;
            tmp2(fig_cont{ci}) = round(255-255*ciperc);
            tmp3(fig_cont{ci}) = round(255-255*ciperc);
        end
        
        contfig = cat(3,contfig(:,:,1),tmp2,tmp3);
        axes(fax{side+8});
        imshow(contfig)
        
        drawnow
        
    end
end
c = sideff-maxsideeff;
ceq = [];


end

function stimpredict_export(optimresults,monrevresults,options,outpth)

% Export results
ea_dispt('Exporting results...')

if ~exist(outpth,'dir')
    mkdir(outpth)
end
if ~exist([outpth,filesep,options.mdlname],'dir')
    mkdir([outpth,filesep,options.mdlname])
end

if options.calcoptdistr
    
    cnames1 = {'Side','Improvement_Combined','Improvement_Rig/Akin','Improvement_Tremor','Side-effect_probability','Amplitude'};
    cnames2 = {'Contact','Percent'};
    
    for side = 1:options.nsides
        save([outpth,filesep,options.mdlname,filesep,[options.sidename{side},'_optimresults']],'optimresults')
        overviewtable(side,:) = table(side,optimresults.eff_comb(:,side)*100,optimresults.eff_ra(:,side)*100,optimresults.eff_tr(:,side)*100,round(optimresults.eff_se(:,side)*100,2),optimresults.current(:,side),'VariableNames',cnames1);
        writetable(overviewtable,[outpth,filesep,options.mdlname,filesep,'optimresults.xlsx'],'Sheet','Overview')
    end
    for side = 1:options.nsides
        for i = 1:options.optimizer.nIter
            itertable = table(flipud((1:options.nconts)'),flipud(optimresults.perc(:,i,side)),'VariableNames',cnames2);
            writetable(itertable,[outpth,filesep,options.mdlname,filesep,'optimresults.xlsx'],'Sheet',[options.sidename{side},'_Iteration_', num2str(i)],'Range','A3')
            writecell({'Amplitude',optimresults.current(i,side)},[outpth,filesep,options.mdlname,filesep,'optimresults.xlsx'],'Sheet',[options.sidename{side},'_Iteration_', num2str(i)],'Range','A1')
        end
    end
end

if options.calcmonrev
    
    for side = 1:options.nsides
        
        comb = monrevresults.eff_comb(:,:,side);
        ra = monrevresults.eff_ra(:,:,side);
        tr = monrevresults.eff_tr(:,:,side);
        sideff = monrevresults.eff_se(:,:,side);
        conts = ones(size(comb)).*(1:options.nconts)';
        amps = ones(size(comb)).*options.ampl;
        
        cnames = {'Contact','Amplitude','Improvement_Combined','Improvement_Rig/Akin','Improvement_Tremor','Side-effect_probability'};
        outputtable = table(conts(:),amps(:),comb(:),ra(:),tr(:),sideff(:),'VariableNames',cnames);
        
        pth = [outpth,filesep,options.mdlname,filesep,options.sidename{side},'Hemisphere_monrev.xlsx'];
        %         warning_fileexists(pth)
        writetable(outputtable,pth)
        
        pth = [outpth,filesep,options.mdlname,filesep,options.sidename{side},'_monrevresults.mat'];
        %         warning_fileexists(pth)
        save(pth,'outputtable')
        
    end
    
end

end
%
% function stimpredict_vizz(monrevresults,options,model,outpth)
%
% if options.calcoptdistr
%     % Nice output figure
% end
%
% if options.calcmonrev
%     % Figure
%     h = figure;
%     set(h, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%
%     for side = 1:options.nsides
%         subplot(2,2,side*2-1)
%         b1 = bar3(monrevresults.eff_comb(:,:,side));
%         for k = 1:length(b1)
%             zdata = b1(k).ZData;
%             b1(k).CData = zdata;
%             b1(k).FaceColor = 'interp';
%         end
%         view(2)
%         title(['Motor improvement - ',options.sidename{side}])
%         xlabel('Amplitude (mA)')
%         ylabel('Contacts')
%         xticklabels(num2cell((0.5:0.5:5)'))
%         colorbar
%
%         subplot(2,2,side*2)
%         b1 = bar3(monrevresults.sideff(:,:,side));
%         for k = 1:length(b1)
%             zdata = b1(k).ZData;
%             b1(k).CData = zdata;
%             b1(k).FaceColor = 'interp';
%         end
%         view(2)
%         title(['Side-effects - ',options.sidename{side}])
%         xlabel('Amplitude (mA)')
%         ylabel('Contacts')
%         xticklabels(num2cell((0.5:0.5:5)'))
%         colorbar
%     end
%
%     if options.export
%         saveas(h,[outpth,filesep,model.options.modelpth,filesep,'Monrevpredict.png'])
%     end
%
% end
%
% end
%
%
