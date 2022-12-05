function generate_EF(ptpth)

% Function to generate E-field templates in patient folder based on lead-DBS SimBio/Fieldtrip pipeline

efpth = [fileparts(mfilename('fullpath')),filesep,'EF_model',filesep];
% Load atlas
load([efpth,'atlas',filesep,'atlas_index.mat'])

% Load reco
load([ptpth,filesep,'ea_reconstruction.mat'])
elstruct = reco.mni;
elstruct.stretchfactor = 0.75;

% Load template of S struct
load([efpth,'Stemplate.mat'])
amp = 5;                                % Default amplitude to generate Efield templates

% Load electrode model
if strcmp(reco.props(1).elmodel,'Medtronic B33005')
    load([efpth,'elmodel',filesep,'medtronic_b33005_elspec.mat'])
elseif strcmp(reco.props(1).elmodel,'Boston Scientific Vercise Directed')
    load([efpth,'elmodel',filesep,'boston_vercise_directed_elspec.mat'])
else
    warning('Electrode model not found')
end


for c = 1:8
    stimset = [amp,0,0,0,0,0,0,0,0];
    stimset(c+1) = -100;
    S = jr_modifySstruct(stimset,stimset,Stemplate,0);
    S.label = ['c',num2str(c-1,'%02d'),'_a',num2str(amp*10,'%02d')];
    
    for side = 1:2  
        execute_SimBioFieldtrip(c,atlases,efpth,elspec,elstruct,ptpth,reco,S,Stemplate,stimset,structures,side)
    end
    
end


% Helpers
function S = jr_modifySstruct(rightsets,leftsets,S,volcur)

rightsets(isnan(rightsets)) = 0;
leftsets(isnan(leftsets)) = 0;
sidelabel = {'R','L'};

% Stim OFF
for side = 1:2
    for s = 1:4
        f1 = [sidelabel{side},'s',num2str(s)];
        for k = 0:7
            f2 = ['k',num2str(k+8*(side-1))];
            S.(f1).(f2).perc = 0;
            S.(f1).(f2).pol = 0;
            S.(f1).(f2).imp = 1;
        end
        S.(f1).amp = 0;
        S.(f1).va = volcur;             % mA
        S.(f1).case.perc = 100;
        S.(f1).case.pol = 2;
        S.(f1).case.imp = 1;
    end
    
end
S.activecontacts{1,1} = [0,0,0,0,0,0,0,0];
S.activecontacts{1,2} = [0,0,0,0,0,0,0,0];
S.amplitude{1,1} = [0,0,0,0];
S.amplitude{1,2} = [0,0,0,0];


% Stim ON

S.Rs1.amp = rightsets(1);
S.Ls1.amp = leftsets(1);
S.active = [1,1];
S.amplitude{1,1} = [rightsets(1),0,0,0];
S.amplitude{1,2} = [leftsets(1),0,0,0];
S.activecontacts{1,1} = abs(rightsets(2:end))>0;
S.activecontacts{1,2} = abs(leftsets(2:end))>0;

rsum = sum(rightsets(2:end));
lsum = sum(leftsets(2:end));

S.Rs1.case.perc = abs(-rsum);
S.Ls1.case.perc = abs(-lsum);
if rsum<0
    S.Rs1.case.pol = 2;
elseif rsum>0
    S.Rs1.case.pol = 1;
else
    S.Rs1.case.pol = 0;
end

if lsum<0
    S.Ls1.case.pol = 2;
elseif lsum>0
    S.Ls1.case.pol = 1;
else
    S.Ls1.case.pol = 0;
end


for cont = 1:8
    
    kr = ['k' num2str(cont-1)];
    S.Rs1.(kr).perc = abs(rightsets(cont+1));
    
    if rightsets(cont+1) > 0
        S.Rs1.(kr).pol = 2;
    elseif rightsets(cont+1) < 0
        S.Rs1.(kr).pol = 1;
    else
        S.Rs1.(kr).pol = 0;
    end
    
    kl = ['k' num2str(cont+7)];
    S.Ls1.(kl).perc = abs(leftsets(cont+1));
    
    if leftsets(cont+1) > 0
        S.Ls1.(kl).pol = 2;
    elseif leftsets(cont+1) < 0
        S.Ls1.(kl).pol = 1;
    else
        S.Ls1.(kl).pol = 0;
    end
    
end

