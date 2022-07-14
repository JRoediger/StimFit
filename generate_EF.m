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
load([efpth,'elmodel',filesep,'Boston_vercise_directed_elspec.mat'])


for c = 1:8
    stimset = [amp,0,0,0,0,0,0,0,0];
    stimset(c+1) = -100;
    S = jr_modifySstruct(stimset,stimset,Stemplate,0);
    S.label = ['c',num2str(c-1,'%02d'),'_a',num2str(amp*10,'%02d')];
    
    for side = 1:2
        
        if c == 1
            ea_dispt('Headmodel needs to be re-calculated. This may take a while...');
            
            cnt=1;
            mesh.tet=[];
            mesh.pnt=[];
            mesh.tissue=[];
            mesh.tissuelabel={'gray','white','contacts','insulation'};
            
            [elfv,ntissuetype,Y,electrode]=ea_buildelfv(elspec,elstruct,side);
            [mesh.tet,mesh.pnt,activeidx,wmboundary,centroids,mesh.tissue,conts]=ea_mesh_electrode([],elfv,ntissuetype,electrode,[],S,side,electrode.numel,Y,elspec);
            
            % replace wmboundary:
            tess = mesh.tet(:,1:4);
            tess = sort(tess,2);
            
            % all faces
            faces=[tess(:,[1 2 3]);tess(:,[1 2 4]); ...
                tess(:,[1 3 4]);tess(:,[2 3 4])];
            
            faces = sortrows(faces);                                        % find replicate faces
            k = find(all(diff(faces)==0,2));
            faces([k;k+1],:) = [];                                          % delete the internal (shared) faces
            
            wmboundary = unique(faces(:))';
            
            meshregions=mesh.tet(:,5);
            mesh.tet=mesh.tet(:,1:4);
            
            ea_dispt('Creating volume conductor...');
            mesh.pnt=mesh.pnt/1000; % in meter
            mesh.unit='m';
            
            try
                vol=ea_ft_headmodel_simbio(mesh,'conductivity',[0.2 0.2 10^8 10^-16]);
            catch % reorder elements so not to be degenerated.
                tmesh=mesh;
                tmesh.tet=tmesh.tet(:,[1 2 4 3]);
                vol=ea_ft_headmodel_simbio(tmesh,'conductivity',[0.2 0.2 10^8 10^-16]);
            end
            mesh.pnt=mesh.pnt*1000; % convert to mm
            mesh.unit='mm';
            
            if ~exist([ptpth,filesep,'current_headmodel',filesep,'template'],'dir')
                mkdir([ptpth,filesep,'current_headmodel',filesep,'template']);
            end
            save([ptpth,filesep,'current_headmodel',filesep,'template',filesep,'headmodel',num2str(side),'.mat'],'vol','mesh','centroids','wmboundary','elfv','meshregions','conts','-v7.3');
            
        else
            ea_dispt('Loading headmodel...');
            load([ptpth,filesep,'current_headmodel',filesep,'template',filesep,'headmodel',num2str(side),'.mat']);
            activeidx = jr_activeidx(S,side,conts,elspec);                  % Load contactidx instead of searching for activecontacts
        end
        
        switch side
            case 1
                sidec='R';
                cnts={'k0','k1','k2','k3','k4','k5','k6','k7'};
            case 2
                sidec='L';
                cnts={'k8','k9','k10','k11','k12','k13','k14','k15'};
        end
        
        if ~isfield(S, 'sources')
            S.sources=1:4;
        end
        
        for source=S.sources
            stimsource=S.([sidec,'s',num2str(source)]);
            
            for cnt=1:length(cnts)
                
                U(cnt)=(stimsource.(cnts{cnt}).perc/100)*stimsource.amp;
                if stimsource.(cnts{cnt}).pol==1
                    U(cnt)=U(cnt)*-1;
                end
            end
            
            Acnt=find(U); % active contact
            if ~isempty(Acnt)
                
                volts=U(U~=0);
                
                % calculate voltage distribution based on dipole
                ea_dispt('Calculating voltage distribution...');
                
                unipolar=1;
                ix=[];
                voltix=[];
                
                cnt=1;
                for ac=Acnt
                    ix=[ix;activeidx(source).con(ac).ix];
                    voltix=[voltix;repmat(U(ac),length(activeidx(source).con(ac).ix),1),...
                        repmat(cnt,length(activeidx(source).con(ac).ix),1)];
                    cnt=cnt+1;
                end
                
                voltix(:,1)=voltix(:,1)/1000; % from mA to A
                
                potential = ea_apply_dbs(vol,ix,voltix,unipolar,0,wmboundary); % output in V. 4 indexes insulating material.
                % save('results','mesh','vol','ix','voltix','unipolar','constvol','wmboundary','potential3v','potential3ma','gradient3v','gradient3ma');
                
                voltix=voltix(:,1); % get rid of index column
                
                ea_dispt('Calculating E-Field...');
                gradient{source} = ea_calc_gradient(vol,potential); % output in V/m.
                
                % stuff by Till to get high EF values for active electrodes
                % this can be adjusted by assigning all tetrahedar belonging to the
                % active electrode a new value:
                % gradient{source}(elec_tet_ix,:) = new_value;
                elec_tet_ix = sub2ind(size(mesh.pnt),vertcat(ix,ix,ix),vertcat(ones(length(ix),1),ones(length(ix),1).*2,ones(length(ix),1).*3));
                elec_tet_ix = find(sum(ismember(mesh.tet,elec_tet_ix),2)==4);
                
                % gradient{source}(elec_tet_ix,:) = repmat(max(gradient{source}),[length(elec_tet_ix),1]); %assign maximum efield value
                tmp = sort(abs(gradient{source}),'descend');
                gradient{source}(elec_tet_ix,:) = repmat(mean(tmp(1:ceil(length(tmp(:,1))*0.001),:)),[length(elec_tet_ix),1]); % choose mean of highest 0.1% as new efield value
                clear tmp
                
            else % empty source..
                gradient{source}=zeros(size(vol.tet,1),3);
            end
            
        end
        
        gradient=gradient{1}+gradient{2}+gradient{3}+gradient{4}; % combined gradient from all sources.
        
        vol.pos=vol.pos*1000; % convert back to mm.
        vat.pos=mean(cat(3,vol.pos(vol.tet(:,1),:),vol.pos(vol.tet(:,2),:),vol.pos(vol.tet(:,3),:),vol.pos(vol.tet(:,4),:)),3); % midpoints of each pyramid
        vat.ET=sqrt(sum(gradient'.^2,1)); % vol.cond(vol.tissue).*ngrad; would be stromstaerke.
        
        % Extract voltage gradients in all dimensions
        vat.Vx=gradient(:,1)';
        vat.Vy=gradient(:,2)';
        vat.Vz=gradient(:,3)';
        
        % remove electrode
        vat = jr_remove_electrode(vat,elstruct,mesh,side,elspec);
        
        % Flip left to right (remember to multiply the voltage gradient in
        % x-direction by -1 for the flipped side... Might make sense to already
        % implement this here since we are already flipping the coordinates
        if side == 2
            vat.pos = ea_flip_lr_nonlinear(vat.pos);
        end
        
        ea_dispt('Calculating interpolant on scattered FEM mesh data...');
        F=scatteredInterpolant(vat.pos(:,1),vat.pos(:,2),vat.pos(:,3),vat.ET','linear','none');
        Fx=scatteredInterpolant(vat.pos(:,1),vat.pos(:,2),vat.pos(:,3),vat.Vx','linear','none');
        Fy=scatteredInterpolant(vat.pos(:,1),vat.pos(:,2),vat.pos(:,3),vat.Vy','linear','none');
        Fz=scatteredInterpolant(vat.pos(:,1),vat.pos(:,2),vat.pos(:,3),vat.Vz','linear','none');
        
        %% Define window and resolution of output nii
        ea_dispt('Converting to equispaced image data...');
        
        res=100;
        gv=cell(3,1);
        
        % Predifined area of 40*40*40mm with resolution of 0.4mm (you should be
        % able to modifiy this without any problems if necessary)
        vxs = 0.4;
        org = [-7.6, -32.4, -26];
        nvox = 100;
        for dim=1:3
            gv{dim}=linspace(org(dim),org(dim)+nvox*vxs-vxs,res);
        end
        
        ea_dispt('Creating nifti header for export...');
        Vvat.mat = [vxs,0,0,org(1);0,vxs,0,org(2);0,0,vxs,org(3);0,0,0,1];
        Vvat.dim=[nvox,nvox,nvox];
        Vvat.dt=[16,0];
        Vvat.n=[1 1];
        Vvat.descrip='lead dbs - vat';
        
        if ~exist([ptpth,filesep,'stimulations',filesep,'template'],'file')
            mkdir([ptpth,filesep,'stimulations',filesep,'template']);
        end
        
        ea_dispt('Filling data with values from interpolant...');
        eeg = F(gv);
        eeg(isnan(eeg))=0;
        
        egx = Fx(gv);
        egx(isnan(eeg))=0;
        
        egy = Fy(gv);
        egy(isnan(eeg))=0;
        
        egz = Fz(gv);
        egz(isnan(eeg))=0;
        
        
        %% Binarized (VTA)
        % I did not use them but might be useful to export VTAs here for sanity checks etc
        ea_dispt('Calculating output file data...');
        thresh = 200;       % Threshold for VTA
        eg=eeg;
        eg=eg>thresh;       % this is the VTA
        
        [xg,yg,zg] = meshgrid(gv{1},gv{2},gv{3});
        
        % I think the maximum radius of the VTA is calculated here... Never
        % actually looked at it
        XYZmax=[max(yg(eg>0)),max(xg(eg>0)),max(zg(eg>0))]; % x and y need to be permuted here (should be correct but wouldnt matter anyways since only serves to calc radius)
        try
            radius=ea_pdist([XYZmax;dpvx]);
        end
        
        %eg=smooth3(eg,'gaussian',[25 25 25]);
        ea_dispt('Calculating volume...');
        
        vatvolume=sum(eg(:))*vxs^3; % returns volume of vat in mm^3
        S.volume(side)=vatvolume;
        
        %% Export stuff
        ea_dispt('Writing files...');
        
        stimname = 'elmov_mA_Homogeneous_20_20';                            % Those are the efield settings StimFit was trained with
        
        % determine stimulation name:
        Vvate=Vvat; Vvatvx=Vvat; Vvatvy=Vvat; Vvatvz=Vvat;
        addname = [S.label,'_'];
        
        if ~exist([ptpth,filesep,'stimulations',filesep,'template',filesep,stimname],'file')
            mkdir([ptpth,filesep,'stimulations',filesep,'template',filesep,stimname]);
        end
        
        switch side
            case 1
                sidename = 'right';
            case 2
                sidename = 'left';
        end
        
        Vvat.fname=[ptpth,filesep,'stimulations',filesep,'template',filesep,stimname,filesep,addname,'vat_',sidename,'.nii'];
        Vvate.fname=[ptpth,filesep,'stimulations',filesep,'template',filesep,stimname,filesep,addname,'vat_efield_',sidename,'.nii'];
        Vvatvx.fname=[ptpth,filesep,'stimulations',filesep,'template',filesep,stimname,filesep,addname,'vat_evecx_',sidename,'.nii'];
        Vvatvy.fname=[ptpth,filesep,'stimulations',filesep,'template',filesep,stimname,filesep,addname,'vat_evecy_',sidename,'.nii'];
        Vvatvz.fname=[ptpth,filesep,'stimulations',filesep,'template',filesep,stimname,filesep,addname,'vat_evecz_',sidename,'.nii'];
        
        Vvate.img=eeg; %permute(eeg,[2,1,3]);
        Vvate.dt=[16,0];
        ea_write_nii(Vvate);
        
        Vvat.img=eg; %permute(eg,[1,2,3]);
        ea_write_nii(Vvat);
        
        Vvatvx.img=egx; %permute(eeg,[2,1,3]);
        Vvatvx.dt=[16,0];
        ea_write_nii(Vvatvx);
        
        Vvatvy.img=egy; %permute(eeg,[2,1,3]);
        Vvatvy.dt=[16,0];
        ea_write_nii(Vvatvy);
        
        Vvatvz.img=egz; %permute(eeg,[2,1,3]);
        Vvatvz.dt=[16,0];
        ea_write_nii(Vvatvz);
        
        ea_dispt('');
        
        
    end
    
end


%% Additional functions

function activeidx = jr_activeidx(S,side,conts,elspec)

% Load contact locations and assign active contact
if side == 1
    sidec = 'R';
else
    sidec = 'L';
end
for con = 1:8
    for source=1:4
        if S.([sidec,'s',num2str(source)]).amp % then this active contact could be from this source since source is active
            if S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).perc % current captured contact is from this source
                activeidx(source).con(con).ix=conts.contactidx{con};
                activeidx(source).con(con).pol=S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).pol;
                activeidx(source).con(con).perc=S.([sidec,'s',num2str(source)]).(['k',num2str(con+8*(side-1)-1)]).perc;
            end
        end
    end
end

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

function potential = ea_apply_dbs(vol,elec,val,unipolar,constvol,boundarynodes)
if constvol
    if unipolar
        dirinodes = [boundarynodes, elec'];
    else
        dirinodes = elec;
    end
    
    rhs = zeros(length(vol.pos),1);
    dirival = zeros(size(vol.pos,1),1);
    dirival(elec) = val(:,1);
else
    
    if unipolar
        dirinodes = boundarynodes;
    else
        dirinodes = 1;
    end
    dirival = zeros(size(vol.pos,1),1);
    
    rhs = zeros(size(vol.pos,1),1);
    uvals=unique(val(:,2));
    if unipolar && length(uvals)==1
        elec_center_id = ea_find_elec_center(elec,vol.pos);
        rhs(elec_center_id) = val(1,1);
    else
        
        for v=1:length(uvals)
            elec_center_id = ea_find_elec_center(elec(val(:,2)==uvals(v)),vol.pos);
            thesevals=val(val(:,2)==uvals(v),1);
            rhs(elec_center_id) = thesevals(1);
        end
        
        %warning('Bipolar constant current stimulation currently not implemented!');
    end
end

[stiff, rhs] = ea_dbs(vol.stiff,rhs,dirinodes,dirival);

potential = ea_sb_solve(stiff,rhs);

function gradient = ea_calc_gradient(vol,potential)
normal = cross(vol.pos(vol.tet(:,4),:)-vol.pos(vol.tet(:,3),:),vol.pos(vol.tet(:,3),:)-vol.pos(vol.tet(:,2),:));
gradient = repmat(potential(vol.tet(:,1))./sum(normal.*(vol.pos(vol.tet(:,1),:)-(vol.pos(vol.tet(:,2),:)+vol.pos(vol.tet(:,3),:)+vol.pos(vol.tet(:,4),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,1),:)-vol.pos(vol.tet(:,4),:),vol.pos(vol.tet(:,4),:)-vol.pos(vol.tet(:,3),:));
gradient = gradient + repmat(potential(vol.tet(:,2))./sum(normal.*(vol.pos(vol.tet(:,2),:)-(vol.pos(vol.tet(:,3),:)+vol.pos(vol.tet(:,4),:)+vol.pos(vol.tet(:,1),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,2),:)-vol.pos(vol.tet(:,1),:),vol.pos(vol.tet(:,1),:)-vol.pos(vol.tet(:,4),:));
gradient = gradient + repmat(potential(vol.tet(:,3))./sum(normal.*(vol.pos(vol.tet(:,3),:)-(vol.pos(vol.tet(:,4),:)+vol.pos(vol.tet(:,1),:)+vol.pos(vol.tet(:,2),:))/3),2),1,3).*normal;
normal = cross(vol.pos(vol.tet(:,3),:)-vol.pos(vol.tet(:,2),:),vol.pos(vol.tet(:,2),:)-vol.pos(vol.tet(:,1),:));
gradient = gradient + repmat(potential(vol.tet(:,4))./sum(normal.*(vol.pos(vol.tet(:,4),:)-(vol.pos(vol.tet(:,1),:)+vol.pos(vol.tet(:,2),:)+vol.pos(vol.tet(:,3),:))/3),2),1,3).*normal;

function center_id = ea_find_elec_center(elec, pos)

center = mean(pos(elec,:));

dist_center = sqrt(sum((pos(elec,:)-repmat(center,length(elec),1)).^2,2));
[dist, elec_id] = min(dist_center);
center_id = elec(elec_id);

function [stiff,rhs] = ea_dbs(stiff,rhs,dirinodes,dirival)

diagonal = diag(stiff);
stiff = stiff + stiff';
rhs = rhs - stiff*dirival;
stiff(dirinodes,:) = 0.0;
stiff(:,dirinodes) = 0.0;
diagonal = -diagonal;
%diagonal(:) = 0;
diagonal(dirinodes) = 1.0;
stiff = stiff + spdiags(diagonal(:),0,length(diagonal),length(diagonal));
rhs(dirinodes) = dirival(dirinodes);

function x = ea_sb_solve(sysmat,vecb)

% SB_SOLVE
%
% $Id: sb_solve.m 8776 2013-11-14 09:04:48Z roboos $
try
    L = ichol(sysmat);
catch
    alpha = max(sum(abs(sysmat),2)./diag(sysmat))-2;
    L = ichol(sysmat, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
end

%scalen
[~,x]=evalc('pcg(sysmat,vecb,10e-10,5000,L,L'',vecb)');

%% SimBio/Fieldtrip Functions

function vol = ea_ft_headmodel_simbio(geom, varargin)

% FT_HEADMODEL_SIMBIO creates a volume conduction model of the head
% using the finite element method (FEM) for EEG. This function takes
% as input a volumetric mesh (hexahedral or tetrahedral) and
% returns as output a volume conduction model which can be used to
% compute leadfields.
%
% This implements
%       ...
%
% Use as
%   vol = ft_headmodel_simbio(geom,'conductivity', conductivities, ...)
%
% The geom is given as a volumetric mesh, using ft_datatype_parcellation
%   geom.pos = vertex positions
%   geom.tet/geom.hex = list of volume elements
%   geom.tissue = tissue assignment for elements
%   geom.tissuelabel = labels correspondig to tissues
%
% Required input arguments should be specified in key-value pairs and have
% to include
%   conductivity   = vector containing tissue conductivities using ordered
%                    corresponding to geom.tissuelabel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To run this on Windows the following packages are necessary:
%
% Microsoft Visual C++ 2008 Redistributable
%
% Intel Visual Fortran Redistributables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

% $Id: ft_headmodel_simbio.m 8445 2013-09-03 10:01:42Z johvor $


% get the optional arguments
conductivity    = ea_ft_getopt(varargin, 'conductivity');

% start with an empty volume conductor
geom = ea_ft_datatype_parcellation(geom);
vol = [];
if isfield(geom,'pos')
    vol.pos = geom.pos;
else
    error('Vertex field is required!')
end

if isfield(geom,'tet')
    vol.tet = geom.tet;
elseif isfield(geom,'hex')
    vol.hex = geom.hex;
else
    error('Connectivity information is required!')
end

if isfield(geom,'tissue')
    vol.tissue = geom.tissue;
else
    error('No element indices declared!')
end

if isempty(conductivity)
    error('No conductivity information!')
end

if length(conductivity) >= length(unique(vol.tissue))
    vol.cond = conductivity;
else
    keyboard
    error('Wrong conductivity information!')
end

if ~isfield(geom,'tissuelabel')
    numlabels = size(unique(geom.tissue),1);
    vol.tissuelabel = {};
    ulabel = unique(labels);
    for i = 1:numlabels
        vol.tissuelabel{i} = num2str(ulabel(i));
    end
else
    vol.tissuelabel = geom.tissuelabel;
end

vol.stiff = ea_sb_calc_stiff(vol);
vol.type = 'simbio';

function val = ea_ft_getopt(opt, key, default, emptymeaningful)

% ea_ft_getopt gets the value of a specified option from a configuration structure
% or from a cell-array with key-value pairs.
%
% Use as
%   val = ea_ft_getopt(s, key, default, emptymeaningful)
% where the input values are
%   s               = structure or cell-array
%   key             = string
%   default         = any valid MATLAB data type
%   emptymeaningful = boolean value (optional, default = 0)
%
% If the key is present as field in the structure, or as key-value
% pair in the cell-array, the corresponding value will be returned.
%
% If the key is not present, ea_ft_getopt will return an empty array.
%
% If the key is present but has an empty value, then the emptymeaningful
% flag specifies whether the empty value or the default value should
% be returned. If emptymeaningful==true, then an empty array will be
% returned. If emptymeaningful==false, then the specified default will
% be returned.
%
% See also FT_SETOPT, FT_CHECKOPT

% Copyright (C) 2011-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ea_ft_getopt.m 7123 2012-12-06 21:21:38Z roboos $

if nargin<3
    default = [];
end

if nargin < 4
    emptymeaningful = 0;
end

if isa(opt, 'struct') || isa(opt, 'config')
    % get the key-value from the structure
    fn = fieldnames(opt);
    if ~any(strcmp(key, fn))
        val = default;
    else
        val = opt.(key);
    end
    
elseif isa(opt, 'cell')
    % get the key-value from the cell-array
    if mod(length(opt),2)
        error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
    end
    
    % the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
    keys = opt(1:2:end);
    vals = opt(2:2:end);
    
    % the following may be faster than cellfun(@ischar, keys)
    valid = false(size(keys));
    for i=1:numel(keys)
        valid(i) = ischar(keys{i});
    end
    
    if ~all(valid)
        error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
    end
    
    hit = find(strcmpi(key, keys));
    if isempty(hit)
        % the requested key was not found
        val = default;
    elseif length(hit)==1
        % the requested key was found
        val = vals{hit};
    else
        error('multiple input arguments with the same name');
    end
    
elseif isempty(opt)
    % no options are specified, return default
    val = default;
end % isstruct or iscell or isempty

if isempty(val) && ~isempty(default) && ~emptymeaningful
    % use the default value instead of the empty input that was specified:
    % this applies for example if you do functionname('key', []), where
    % the empty is meant to indicate that the user does not know or care
    % what the value is
    val = default;
end

function parcellation = ea_ft_datatype_parcellation(parcellation, varargin)

% FT_DATATYPE_PARCELLATION describes the FieldTrip MATLAB structure for parcellated
% cortex-based data and atlases. A parcellation can either be indexed or probabilistic
% (see below).
%
% A parcellation describes the tissue types for each of the surface elements.
% Parcellations are often, but not always labeled. A parcellatoin can be used to
% estimate the activity from MEG data in a known region of interest. A surface-based
% atlas is basically a very detailled parcellation with an anatomical label for each
% vertex.
%
% An example of a surface based Brodmann parcellation looks like this
%
%              pos: [8192x3]         positions of the vertices forming the cortical sheet
%              tri: [16382x3]        triangles of the cortical sheet
%         coordsys: 'ctf'            the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'             the units in which the coordinate system is expressed
%         brodmann: [8192x1 uint8]   values from 1 to N, the value 0 means unknown
%    brodmannlabel: {Nx1 cell}
%
% An alternative representation of this parcellation is
%
%              pos: [8192x3]           positions of the vertices forming the cortical sheet
%              tri: [16382x3]          triangles of the cortical sheet
%         coordsys: 'ctf'              the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'               the units in which the coordinate system is expressed
%  Brodmann_Area_1: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  Brodmann_Area_2: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  Brodmann_Area_3: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  ...
%
% The examples above demonstrate that a parcellation can be indexed, i.e. consisting of
% subsequent integer numbers (1, 2, ...) or probabilistic, consisting of real numbers
% ranging from 0 to 1 that represent probabilities between 0% and 100%. An extreme case
% is one where the probability is either 0 or 1, in which case the probability can be
% represented as a binary or logical array.
%
% The only difference to the source data structure is that the parcellation structure
% contains the additional fields xxx and xxxlabel. See FT_DATATYPE_SOURCE for further
% details.
%
% Required fields:
%   - pos
%
% Optional fields:
%   - tri, coordsys, unit
%
% Deprecated fields:
%   - none
%
% Obsoleted fields:
%   - none
%
% Revision history:
% (2012/latest) The initial version was defined in accordance with the representation of
% a voxel-based segmentation.
%
% See also ea_ft_datatype, FT_DATATYPE_SOURCE, ea_ft_datatype_segmentation

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_datatype_parcellation.m 10213 2015-02-11 19:38:33Z roboos $

% get the optional input arguments, which should be specified as key-value pairs
version           = ea_ft_getopt(varargin, 'version', 'latest');
parcellationstyle = ea_ft_getopt(varargin, 'parcellationstyle');  % can be indexed or probabilistic

if strcmp(version, 'latest')
    parcelversion = '2012';
    sourceversion = 'latest';
    clear version
else
    parcelversion = version;
    sourceversion = version;
    clear version
end

if isempty(parcellation)
    return;
end

switch parcelversion
    case '2012'
        
        if isfield(parcellation, 'pnt')
            parcellation.pos = parcellation.pnt;
            parcellation = rmfield(parcellation, 'pnt');
        end
        
        % convert the inside/outside fields, they should be logical rather than an index
        if isfield(parcellation, 'inside')
            parcellation = ea_fixinside(parcellation, 'logical');
        end
        
        dim = size(parcellation.pos,1);
        
        % make a list of fields that represent a parcellation
        fn = fieldnames(parcellation);
        fn = setdiff(fn, 'inside'); % exclude the inside field from any conversions
        sel = false(size(fn));
        for i=1:numel(fn)
            sel(i) = isnumeric(parcellation.(fn{i})) && numel(parcellation.(fn{i}))==dim;
        end
        % only consider numeric fields of the correct size
        fn = fn(sel);
        
        % determine whether the style of the input fields is probabilistic or indexed
        [indexed, probabilistic] = ea_determine_segmentationstyle(parcellation, fn, dim);
        
        % ignore the fields that do not contain a parcellation
        sel = indexed | probabilistic;
        fn            = fn(sel);
        indexed       = indexed(sel);
        probabilistic = probabilistic(sel);
        
        if ~any(probabilistic) && ~any(indexed)
            % rather than being described with a tissue label for each vertex
            % it can also be described with a tissue label for each surface or volme element
            for i = 1:length(fn)
                fname = fn{i};
                switch fname
                    case 'tri'
                        dim = size(parcellation.tri,1);
                    case 'hex'
                        dim = size(parcellation.hex,1);
                    case 'tet'
                        dim = size(parcellation.tet,1);
                end
            end
            [indexed, probabilistic] = ea_determine_segmentationstyle(parcellation, fn, dim);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that the parcellation is internally consistent
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if any(probabilistic)
            parcellation = ea_fixsegmentation(parcellation, fn(probabilistic), 'probabilistic');
        end
        
        if any(indexed)
            parcellation = ea_fixsegmentation(parcellation, fn(indexed), 'indexed');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert the parcellation to the desired style
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(parcellationstyle, 'indexed') && any(probabilistic)
            parcellation  = convert_segmentationstyle(parcellation, fn(probabilistic), [dim 1], 'indexed');
        elseif strcmp(parcellationstyle, 'probabilistic') && any(indexed)
            parcellation  = convert_segmentationstyle(parcellation, fn(indexed), [dim 1], 'probabilistic');
        end % converting converting to desired style
        
    otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('unsupported version "%s" for parcellation datatype', parcelversion);
end

% the parcellation is a speciat type of volume structure, so ensure that it also fulfills the requirements for that
parcellation = ea_ft_datatype_source(parcellation, 'version', sourceversion);

function [indexed, probabilistic] = ea_determine_segmentationstyle(segmentation, fn, dim)

% DETERMINE_SEGMENTATIONSTYLE is a helper function that determines the type of segmentation
% contained in each of the fields. It is used by ea_ft_datatype_segmentation and
% ft_datatype_parcellation.
%
% See also FIXSEGMENTATION, CONVERT_SEGMENTATIONSTYLE

indexed       = false(size(fn));
probabilistic = false(size(fn));

% determine for each of the fields whether it is probabilistic, indexed or something else
for i=1:numel(fn)
    if numel(segmentation.(fn{i}))~=prod(dim)
        % this does not look like a segmentation
        continue
    elseif strcmp(fn{i}, 'anatomy')
        % this should not be interpreted as segmentation, also not when it is a uint8 or uint16 representation
        continue
    else
        if isfield(segmentation, [fn{i} 'label'])
            % the xxxlabel field exists, which only makes sense for an indexed representation
            probabilistic(i) = false;
            indexed(i)       = true;
        else
            % this looks like a segmentation
            tmp = segmentation.(fn{i});
            tmp = tmp(:);       % convert to vector
            sel = isnan(tmp);   % find NaN values
            if any(sel)
                % note that the the finding and removing of NaNs have been separated to speed up the code
                tmp = tmp(~sel);  % remove NaN values
            end
            clear sel
            probabilistic(i) =  islogical(tmp) || all(tmp>=-0.001 & tmp<=1.001); % allow some roundoff error
            indexed(i)       = ~islogical(tmp) && all(abs(tmp - round(tmp))<1000*eps);
            
            if probabilistic(i) && indexed(i)
                % the xxxlabel does not exist, so treat it as a probabilistic representation
                probabilistic(i) = true;
                indexed(i)       = false;
            end
        end
    end
end % for each of the fields

function source = ea_ft_datatype_source(source, varargin)

% FT_DATATYPE_SOURCE describes the FieldTrip MATLAB structure for data that is
% represented at the source level. This is typically obtained with a beamformer of
% minimum-norm source reconstruction using FT_SOURCEANALYSIS.
%
% An example of a source structure obtained after performing DICS (a frequency
% domain beamformer scanning method) is shown here
%
%           pos: [6732x3 double]       positions at which the source activity could have been estimated
%        inside: [6732x1 logical]      boolean vector that indicates at which positions the source activity was estimated
%           dim: [xdim ydim zdim]      if the positions can be described as a 3D regular grid, this contains the
%                                       dimensionality of the 3D volume
%     cumtapcnt: [120x1 double]        information about the number of tapers per original trial
%          time: 0.100                 the latency at which the activity is estimated (in seconds)
%          freq: 30                    the frequency at which the activity is estimated (in Hz)
%           pow: [6732x120 double]     the estimated power at each source position
%     powdimord: 'pos_rpt'             defines how the numeric data has to be interpreted,
%                                       in this case 6732 dipole positions x 120 repetitions (i.e. trials)
%           cfg: [1x1 struct]          the configuration used by the function that generated this data structure
%
% Required fields:
%   - pos
%
% Optional fields:
%   - time, freq, pow, coh, eta, mom, ori, cumtapcnt, dim, transform, inside, cfg, dimord, other fields with a dimord
%
% Deprecated fields:
%   - method, outside
%
% Obsoleted fields:
%   - xgrid, ygrid, zgrid, transform, latency, frequency
%
% Historical fields:
%   - avg, cfg, cumtapcnt, df, dim, freq, frequency, inside, method,
%   outside, pos, time, trial, vol, see bug2513
%
% Revision history:
%
% (2014) The subfields in the avg and trial fields are now present in the
% main structure, e.g. source.avg.pow is now source.pow. Furthermore, the
% inside is always represented as logical vector.
%
% (2011) The source representation should always be irregular, i.e. not
% a 3-D volume, contain a "pos" field and not contain a "transform".
%
% (2010) The source structure should contain a general "dimord" or specific
% dimords for each of the fields. The source reconstruction in the avg and
% trial substructures has been moved to the toplevel.
%
% (2007) The xgrid/ygrid/zgrid fields have been removed, because they are
% redundant.
%
% (2003) The initial version was defined
%
% See also ea_ft_datatype, FT_DATATYPE_COMP, FT_DATATYPE_DIP, FT_DATATYPE_FREQ,
% FT_DATATYPE_MVAR, FT_DATATYPE_RAW, FT_DATATYPE_SOURCE, FT_DATATYPE_SPIKE,
% FT_DATATYPE_TIMELOCK, FT_DATATYPE_VOLUME

% Copyright (C) 2013-2014, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_datatype_source.m 10265 2015-03-04 12:20:41Z jansch $

% FIXME: I am not sure whether the removal of the xgrid/ygrid/zgrid fields
% was really in 2007

% get the optional input arguments, which should be specified as key-value pairs
version = ea_ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest') || strcmp(version, 'upcoming')
    version = '2014';
end

if isempty(source)
    return;
end

% old data structures may use latency/frequency instead of time/freq. It is
% unclear when these were introduced and removed again, but they were never
% used by any fieldtrip function itself
if isfield(source, 'frequency')
    source.freq = source.frequency;
    source      = rmfield(source, 'frequency');
end
if isfield(source, 'latency')
    source.time = source.latency;
    source      = rmfield(source, 'latency');
end

switch version
    case '2014'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that it has individual source positions
        source = ea_fixpos(source);
        
        % ensure that it is always logical
        source = ea_fixinside(source, 'logical');
        
        % remove obsolete fields
        if isfield(source, 'method')
            source = rmfield(source, 'method');
        end
        if isfield(source, 'transform')
            source = rmfield(source, 'transform');
        end
        if isfield(source, 'xgrid')
            source = rmfield(source, 'xgrid');
        end
        if isfield(source, 'ygrid')
            source = rmfield(source, 'ygrid');
        end
        if isfield(source, 'zgrid')
            source = rmfield(source, 'zgrid');
        end
        
        if isfield(source, 'avg') && isstruct(source.avg)
            % move the average fields to the main structure
            fn = fieldnames(source.avg);
            for i=1:length(fn)
                dat = source.avg.(fn{i});
                if isequal(size(dat), [1 size(source.pos,1)])
                    source.(fn{i}) = dat';
                else
                    source.(fn{i}) = dat;
                end
                clear dat
            end % j
            source = rmfield(source, 'avg');
        end
        
        if isfield(source, 'inside')
            % the inside is by definition logically indexed
            probe = find(source.inside, 1, 'first');
        else
            % just take the first source position
            probe = 1;
        end
        
        if isfield(source, 'trial') && isstruct(source.trial)
            npos = size(source.pos,1);
            
            % concatenate the fields for each trial and move them to the main structure
            fn = fieldnames(source.trial);
            
            for i=1:length(fn)
                % some fields are descriptive and hence identical over trials
                if strcmp(fn{i}, 'csdlabel')
                    source.csdlabel = dat;
                    continue
                end
                
                % start with the first trial
                dat    = source.trial(1).(fn{i});
                datsiz = ea_getdimsiz(source, fn{i});
                nrpt   = datsiz(1);
                datsiz = datsiz(2:end);
                
                
                if iscell(dat)
                    datsiz(1) = nrpt; % swap the size of pos with the size of rpt
                    val  = cell(npos,1);
                    indx = find(source.inside);
                    for k=1:length(indx)
                        val{indx(k)}          = nan(datsiz);
                        val{indx(k)}(1,:,:,:) = dat{indx(k)};
                    end
                    % concatenate all data as {pos}_rpt_etc
                    for j=2:nrpt
                        dat = source.trial(j).(fn{i});
                        for k=1:length(indx)
                            val{indx(k)}(j,:,:,:) = dat{indx(k)};
                        end
                        
                    end % for all trials
                    source.(fn{i}) = val;
                    
                else
                    % concatenate all data as pos_rpt_etc
                    val = nan([datsiz(1) nrpt datsiz(2:end)]);
                    val(:,1,:,:,:) = dat(:,:,:,:);
                    for j=2:length(source.trial)
                        dat = source.trial(j).(fn{i});
                        val(:,j,:,:,:) = dat(:,:,:,:);
                    end % for all trials
                    source.(fn{i}) = val;
                    
                    %         else
                    %           siz = size(dat);
                    %           if prod(siz)==npos
                    %             siz = [npos nrpt];
                    %           elseif siz(1)==npos
                    %             siz = [npos nrpt siz(2:end)];
                    %           end
                    %           val = nan(siz);
                    %           % concatenate all data as pos_rpt_etc
                    %           val(:,1,:,:,:) = dat(:);
                    %           for j=2:length(source.trial)
                    %             dat = source.trial(j).(fn{i});
                    %             val(:,j,:,:,:) = dat(:);
                    %           end % for all trials
                    %           source.(fn{i}) = val;
                    
                end
            end % for each field
            
            source = rmfield(source, 'trial');
            
        end % if trial
        
        % ensure that it has a dimord (or multiple for the different fields)
        source = ea_fixdimord(source);
        
        
    case '2011'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that it has individual source positions
        source = ea_fixpos(source);
        
        % remove obsolete fields
        if isfield(source, 'xgrid')
            source = rmfield(source, 'xgrid');
        end
        if isfield(source, 'ygrid')
            source = rmfield(source, 'ygrid');
        end
        if isfield(source, 'zgrid')
            source = rmfield(source, 'zgrid');
        end
        if isfield(source, 'transform')
            source = rmfield(source, 'transform');
        end
        
        % ensure that it has a dimord (or multiple for the different fields)
        source = ea_fixdimord(source);
        
    case '2010'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that it has individual source positions
        source = ea_fixpos(source);
        
        % remove obsolete fields
        if isfield(source, 'xgrid')
            source = rmfield(source, 'xgrid');
        end
        if isfield(source, 'ygrid')
            source = rmfield(source, 'ygrid');
        end
        if isfield(source, 'zgrid')
            source = rmfield(source, 'zgrid');
        end
        
        % ensure that it has a dimord (or multiple for the different fields)
        source = ea_fixdimord(source);
        
    case '2007'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ensure that it has individual source positions
        source = ea_fixpos(source);
        
        % remove obsolete fields
        if isfield(source, 'dimord')
            source = rmfield(source, 'dimord');
        end
        if isfield(source, 'xgrid')
            source = rmfield(source, 'xgrid');
        end
        if isfield(source, 'ygrid')
            source = rmfield(source, 'ygrid');
        end
        if isfield(source, 'zgrid')
            source = rmfield(source, 'zgrid');
        end
        
    case '2003'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(source, 'dimord')
            source = rmfield(source, 'dimord');
        end
        
        if ~isfield(source, 'xgrid') || ~isfield(source, 'ygrid') || ~isfield(source, 'zgrid')
            if isfield(source, 'dim')
                minx = min(source.pos(:,1));
                maxx = max(source.pos(:,1));
                miny = min(source.pos(:,2));
                maxy = max(source.pos(:,2));
                minz = min(source.pos(:,3));
                maxz = max(source.pos(:,3));
                source.xgrid = linspace(minx, maxx, source.dim(1));
                source.ygrid = linspace(miny, maxy, source.dim(2));
                source.zgrid = linspace(minz, maxz, source.dim(3));
            end
        end
        
    otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('unsupported version "%s" for source datatype', version);
end

function source = ea_fixpos(source)
if ~isfield(source, 'pos')
    if isfield(source, 'xgrid') && isfield(source, 'ygrid') && isfield(source, 'zgrid')
        source.pos = ea_grid2pos(source.xgrid, source.ygrid, source.zgrid);
    elseif isfield(source, 'dim') && isfield(source, 'transform')
        source.pos = ea_dim2pos(source.dim, source.transform);
    else
        error('cannot reconstruct individual source positions');
    end
end

function pos = ea_dim2pos(dim, transform)
[X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos = [X(:) Y(:) Z(:)];
pos = ea_ft_warp_apply(transform, pos, 'homogenous');

function [warped] = ea_ft_warp_apply(M, input, method, tol)

% ea_ft_warp_apply performs a 3D linear or nonlinear transformation on the input
% coordinates, similar to those in AIR 3.08. You can find technical
% documentation on warping in general at http://bishopw.loni.ucla.edu/AIR3
%
% Use as
%   [warped] = ea_ft_warp_apply(M, input, method, tol)
% where
%   M        vector or matrix with warping parameters
%   input    Nx3 matrix with coordinates
%   warped   Nx3 matrix with coordinates
%   method   string describing the warping method
%   tol      (optional) value determining the numerical precision of the
%             output, to deal with numerical round off imprecisions due to
%             the warping
%
% The methods 'nonlin0', 'nonlin2' ... 'nonlin5' specify a
% polynomial transformation. The size of the transformation matrix
% depends on the order of the warp
%   zeroth order :  1 parameter  per coordinate (translation)
%   first  order :  4 parameters per coordinate (total 12, affine)
%   second order : 10 parameters per coordinate
%   third  order : 20 parameters per coordinate
%   fourth order : 35 parameters per coordinate
%   fifth  order : 56 parameters per coordinate (total 168)
% The size of M should be 3xP, where P is the number of parameters
% per coordinate. Alternatively, you can specify the method to be
% 'nonlinear', where the order will be determined from the size of
% the matrix M.
%
% If the method 'homogeneous' is selected, the input matrix M should be
% a 4x4 homogenous transformation matrix.
%
% If the method 'sn2individual' or 'individual2sn' is selected, the input
% M should be a structure based on nonlinear (warping) normalisation parameters
% created by SPM8 for alignment between an individual structural MRI and the
% template MNI brain.  These options call private functions of the same name.
% M will have subfields like this:
%     Affine: [4x4 double]
%         Tr: [4-D double]
%         VF: [1x1 struct]
%         VG: [1x1 struct]
%      flags: [1x1 struct]
%
% If any other method is selected, it is assumed that it specifies
% the name of an auxiliary function that will, when given the input
% parameter vector M, return an 4x4 homogenous transformation
% matrix. Supplied functions in the warping toolbox are translate,
% rotate, scale, rigidbody, globalrescale, traditional, affine,
% perspective.
%
% See also FT_WARP_OPTIM, FT_WARP_ERROR

% Copyright (C) 2000-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ea_ft_warp_apply.m 10132 2015-01-27 16:08:29Z johzum $

if nargin<4
    tol = [];
end

if nargin<3 && all(size(M)==4)
    % no specific transformation mode has been selected
    % it looks like a homogenous transformation matrix
    method = 'homogeneous';
elseif nargin<3
    % the default method is 'nonlinear'
    method = 'nonlinear';
end

if size(input,2)==2
    % convert the input points from 2D to 3D representation
    input(:,3) = 0;
    input3d = false;
else
    input3d = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear warping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(strcmp(method, {'nonlinear', 'nonlin0', 'nonlin1', 'nonlin2', 'nonlin3', 'nonlin4', 'nonlin5'}))
    x = input(:,1);
    y = input(:,2);
    z = input(:,3);
    s = size(M);
    
    if s(1)~=3
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin0') && s(2)~=1
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin1') && s(2)~=4
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin2') && s(2)~=10
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin3') && s(2)~=20
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin4') && s(2)~=35
        error('invalid size of nonlinear transformation matrix');
    elseif strcmp(method, 'nonlin5') && s(2)~=56
        error('invalid size of nonlinear transformation matrix');
    end
    
    if s(2)==1
        % this is a translation, which in a strict sense is not the 0th order nonlinear transformation
        xx = M(1,1) + x;
        yy = M(2,1) + y;
        zz = M(3,1) + z;
    elseif s(2)==4
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z;
    elseif s(2)==10
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z;
    elseif s(2)==20
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z;
    elseif s(2)==35
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z + M(1,21)*x.*x.*x.*x + M(1,22)*x.*x.*x.*y + M(1,23)*x.*x.*x.*z + M(1,24)*x.*x.*y.*y + M(1,25)*x.*x.*y.*z + M(1,26)*x.*x.*z.*z + M(1,27)*x.*y.*y.*y + M(1,28)*x.*y.*y.*z + M(1,29)*x.*y.*z.*z + M(1,30)*x.*z.*z.*z + M(1,31)*y.*y.*y.*y + M(1,32)*y.*y.*y.*z + M(1,33)*y.*y.*z.*z + M(1,34)*y.*z.*z.*z + M(1,35)*z.*z.*z.*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z + M(2,21)*x.*x.*x.*x + M(2,22)*x.*x.*x.*y + M(2,23)*x.*x.*x.*z + M(2,24)*x.*x.*y.*y + M(2,25)*x.*x.*y.*z + M(2,26)*x.*x.*z.*z + M(2,27)*x.*y.*y.*y + M(2,28)*x.*y.*y.*z + M(2,29)*x.*y.*z.*z + M(2,30)*x.*z.*z.*z + M(2,31)*y.*y.*y.*y + M(2,32)*y.*y.*y.*z + M(2,33)*y.*y.*z.*z + M(2,34)*y.*z.*z.*z + M(2,35)*z.*z.*z.*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z + M(3,21)*x.*x.*x.*x + M(3,22)*x.*x.*x.*y + M(3,23)*x.*x.*x.*z + M(3,24)*x.*x.*y.*y + M(3,25)*x.*x.*y.*z + M(3,26)*x.*x.*z.*z + M(3,27)*x.*y.*y.*y + M(3,28)*x.*y.*y.*z + M(3,29)*x.*y.*z.*z + M(3,30)*x.*z.*z.*z + M(3,31)*y.*y.*y.*y + M(3,32)*y.*y.*y.*z + M(3,33)*y.*y.*z.*z + M(3,34)*y.*z.*z.*z + M(3,35)*z.*z.*z.*z;
    elseif s(2)==56
        xx = M(1,1) + M(1,2)*x + M(1,3)*y + M(1,4)*z + M(1,5)*x.*x + M(1,6)*x.*y + M(1,7)*x.*z + M(1,8)*y.*y + M(1,9)*y.*z + M(1,10)*z.*z + M(1,11)*x.*x.*x + M(1,12)*x.*x.*y + M(1,13)*x.*x.*z + M(1,14)*x.*y.*y + M(1,15)*x.*y.*z + M(1,16)*x.*z.*z + M(1,17)*y.*y.*y + M(1,18)*y.*y.*z + M(1,19)*y.*z.*z + M(1,20)*z.*z.*z + M(1,21)*x.*x.*x.*x + M(1,22)*x.*x.*x.*y + M(1,23)*x.*x.*x.*z + M(1,24)*x.*x.*y.*y + M(1,25)*x.*x.*y.*z + M(1,26)*x.*x.*z.*z + M(1,27)*x.*y.*y.*y + M(1,28)*x.*y.*y.*z + M(1,29)*x.*y.*z.*z + M(1,30)*x.*z.*z.*z + M(1,31)*y.*y.*y.*y + M(1,32)*y.*y.*y.*z + M(1,33)*y.*y.*z.*z + M(1,34)*y.*z.*z.*z + M(1,35)*z.*z.*z.*z + M(1,36)*x.*x.*x.*x.*x + M(1,37)*x.*x.*x.*x.*y + M(1,38)*x.*x.*x.*x.*z + M(1,39)*x.*x.*x.*y.*y + M(1,40)*x.*x.*x.*y.*z + M(1,41)*x.*x.*x.*z.*z + M(1,42)*x.*x.*y.*y.*y + M(1,43)*x.*x.*y.*y.*z + M(1,44)*x.*x.*y.*z.*z + M(1,45)*x.*x.*z.*z.*z + M(1,46)*x.*y.*y.*y.*y + M(1,47)*x.*y.*y.*y.*z + M(1,48)*x.*y.*y.*z.*z + M(1,49)*x.*y.*z.*z.*z + M(1,50)*x.*z.*z.*z.*z + M(1,51)*y.*y.*y.*y.*y + M(1,52)*y.*y.*y.*y.*z + M(1,53)*y.*y.*y.*z.*z + M(1,54)*y.*y.*z.*z.*z + M(1,55)*y.*z.*z.*z.*z + M(1,56)*z.*z.*z.*z.*z;
        yy = M(2,1) + M(2,2)*x + M(2,3)*y + M(2,4)*z + M(2,5)*x.*x + M(2,6)*x.*y + M(2,7)*x.*z + M(2,8)*y.*y + M(2,9)*y.*z + M(2,10)*z.*z + M(2,11)*x.*x.*x + M(2,12)*x.*x.*y + M(2,13)*x.*x.*z + M(2,14)*x.*y.*y + M(2,15)*x.*y.*z + M(2,16)*x.*z.*z + M(2,17)*y.*y.*y + M(2,18)*y.*y.*z + M(2,19)*y.*z.*z + M(2,20)*z.*z.*z + M(2,21)*x.*x.*x.*x + M(2,22)*x.*x.*x.*y + M(2,23)*x.*x.*x.*z + M(2,24)*x.*x.*y.*y + M(2,25)*x.*x.*y.*z + M(2,26)*x.*x.*z.*z + M(2,27)*x.*y.*y.*y + M(2,28)*x.*y.*y.*z + M(2,29)*x.*y.*z.*z + M(2,30)*x.*z.*z.*z + M(2,31)*y.*y.*y.*y + M(2,32)*y.*y.*y.*z + M(2,33)*y.*y.*z.*z + M(2,34)*y.*z.*z.*z + M(2,35)*z.*z.*z.*z + M(2,36)*x.*x.*x.*x.*x + M(2,37)*x.*x.*x.*x.*y + M(2,38)*x.*x.*x.*x.*z + M(2,39)*x.*x.*x.*y.*y + M(2,40)*x.*x.*x.*y.*z + M(2,41)*x.*x.*x.*z.*z + M(2,42)*x.*x.*y.*y.*y + M(2,43)*x.*x.*y.*y.*z + M(2,44)*x.*x.*y.*z.*z + M(2,45)*x.*x.*z.*z.*z + M(2,46)*x.*y.*y.*y.*y + M(2,47)*x.*y.*y.*y.*z + M(2,48)*x.*y.*y.*z.*z + M(2,49)*x.*y.*z.*z.*z + M(2,50)*x.*z.*z.*z.*z + M(2,51)*y.*y.*y.*y.*y + M(2,52)*y.*y.*y.*y.*z + M(2,53)*y.*y.*y.*z.*z + M(2,54)*y.*y.*z.*z.*z + M(2,55)*y.*z.*z.*z.*z + M(2,56)*z.*z.*z.*z.*z;
        zz = M(3,1) + M(3,2)*x + M(3,3)*y + M(3,4)*z + M(3,5)*x.*x + M(3,6)*x.*y + M(3,7)*x.*z + M(3,8)*y.*y + M(3,9)*y.*z + M(3,10)*z.*z + M(3,11)*x.*x.*x + M(3,12)*x.*x.*y + M(3,13)*x.*x.*z + M(3,14)*x.*y.*y + M(3,15)*x.*y.*z + M(3,16)*x.*z.*z + M(3,17)*y.*y.*y + M(3,18)*y.*y.*z + M(3,19)*y.*z.*z + M(3,20)*z.*z.*z + M(3,21)*x.*x.*x.*x + M(3,22)*x.*x.*x.*y + M(3,23)*x.*x.*x.*z + M(3,24)*x.*x.*y.*y + M(3,25)*x.*x.*y.*z + M(3,26)*x.*x.*z.*z + M(3,27)*x.*y.*y.*y + M(3,28)*x.*y.*y.*z + M(3,29)*x.*y.*z.*z + M(3,30)*x.*z.*z.*z + M(3,31)*y.*y.*y.*y + M(3,32)*y.*y.*y.*z + M(3,33)*y.*y.*z.*z + M(3,34)*y.*z.*z.*z + M(3,35)*z.*z.*z.*z + M(3,36)*x.*x.*x.*x.*x + M(3,37)*x.*x.*x.*x.*y + M(3,38)*x.*x.*x.*x.*z + M(3,39)*x.*x.*x.*y.*y + M(3,40)*x.*x.*x.*y.*z + M(3,41)*x.*x.*x.*z.*z + M(3,42)*x.*x.*y.*y.*y + M(3,43)*x.*x.*y.*y.*z + M(3,44)*x.*x.*y.*z.*z + M(3,45)*x.*x.*z.*z.*z + M(3,46)*x.*y.*y.*y.*y + M(3,47)*x.*y.*y.*y.*z + M(3,48)*x.*y.*y.*z.*z + M(3,49)*x.*y.*z.*z.*z + M(3,50)*x.*z.*z.*z.*z + M(3,51)*y.*y.*y.*y.*y + M(3,52)*y.*y.*y.*y.*z + M(3,53)*y.*y.*y.*z.*z + M(3,54)*y.*y.*z.*z.*z + M(3,55)*y.*z.*z.*z.*z + M(3,56)*z.*z.*z.*z.*z;
    else
        error('invalid size of nonlinear transformation matrix');
    end
    
    warped = [xx yy zz];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % linear warping using homogenous coordinate transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(method, 'homogenous') || strcmp(method, 'homogeneous')
    if all(size(M)==3)
        % convert the 3x3 homogenous transformation matrix (corresponding with 2D)
        % into a 4x4 homogenous transformation matrix (corresponding with 3D)
        M = [
            M(1,1) M(1,2)  0  M(1,3)
            M(2,1) M(2,2)  0  M(2,3)
            0      0       0  0
            M(3,1) M(3,2)  0  M(3,3)
            ];
    end
    
    %warped = M * [input'; ones(1, size(input, 1))];
    %warped = warped(1:3,:)';
    
    % below achieves the same as lines 154-155
    warped = [input ones(size(input, 1),1)]*M(1:3,:)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % using external function that returns a homogeneous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif exist(method, 'file') && ~isa(M, 'struct')
    % get the homogenous transformation matrix
    H = feval(method, M);
    warped = ea_ft_warp_apply(H, input, 'homogeneous');
    
elseif strcmp(method, 'sn2individual') && isa(M, 'struct')
    % use SPM structure with parameters for an inverse warp
    % from normalized space to individual, can be non-linear
    warped = sn2individual(M, input);
    
elseif strcmp(method, 'individual2sn') && isa(M, 'struct')
    % use SPM structure with parameters for a warp from
    % individual space to normalized space, can be non-linear
    %error('individual2sn is not yet implemented');
    warped = individual2sn(M, input);
else
    error('unrecognized transformation method');
end

if ~input3d
    % convert from 3D back to 2D representation
    warped = warped(:,1:2);
end

if ~isempty(tol)
    if tol>0
        warped = fix(warped./tol)*tol;
    end
end

function [source] = ea_fixinside(source, opt)

% FIXINSIDE ensures that the region of interest (which is indicated by the
% field "inside") is consistently defined for source structures and volume
% structures. Furthermore, it solves backward compatibility problems.
%
% Use as
%   [source] = fixinside(source, 'logical');
% or
%   [source] = fixinside(source, 'index');

% Copyright (C) 2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: fixinside.m 9663 2014-06-22 07:06:19Z roboos $


if nargin<2
    opt = 'logical';
end

if ~isfield(source, 'inside')
    if isfield(source, 'pos')
        % assume that all positions are inside the region of interest
        source.inside  = [1:size(source.pos,1)]';
        source.outside = [];
    elseif isfield(source, 'dim')
        source.inside  = [1:prod(source.dim)]';
        source.outside = [];
    end
end

if ~isfield(source, 'inside')
    % nothing to do
    return;
end

% determine the format
if isa(source.inside, 'logical')
    logicalfmt = 1;
elseif all(source.inside(:)==0 | source.inside(:)==1)
    source.inside = logical(source.inside);
    logicalfmt = 1;
else
    logicalfmt = 0;
end

if ~logicalfmt && strcmp(opt, 'logical')
    % convert to a logical array
    if ~isfield(source, 'outside')
        source.outside = [];
    end
    inside(source.inside)  = (1==1);  % true
    inside(source.outside) = (1==0);  % false
    source.inside = inside(:);
    if isfield(source, 'outside')
        source = rmfield(source, 'outside');
    end
elseif logicalfmt && strcmp(opt, 'index')
    % convert to a vectors with indices
    tmp = source.inside;
    source.inside  = find( tmp(:));
    source.outside = find(~tmp(:));
else
    % nothing to do
end

function [data] = ea_fixdimord(data)

% FIXDIMORD ensures consistency between the dimord string and the axes
% that describe the data dimensions. The main purpose of this function
% is to ensure backward compatibility of all functions with data that has
% been processed by older FieldTrip versions
%
% Use as
%   [data] = fixdimord(data)
% This will modify the data.dimord field to ensure consistency.
% The name of the axis is the same as the name of the dimord, i.e. if
% dimord='freq_time', then data.freq and data.time should be present.
%
% The default dimensions in the data are described by
%  'time'
%  'freq'
%  'chan'
%  'chancmb'
%  'refchan'
%  'subj'
%  'rpt'
%  'rpttap'
%  'pos'
%  'ori'
%  'rgb'
%  'comp'
%  'voxel'

% Copyright (C) 2009-2014, Robert Oostenveld, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: fixdimord.m 9972 2014-11-19 08:09:34Z roboos $

% if nargin<2, keepsourcedimord = 0; end
%
% if any(ea_ft_datatype(data, {'source', 'volume'})) && isfield(data, 'dimord') && ~keepsourcedimord
%   % the old source data representation does not have a dimord, whereas the new source data representation does have a dimord
%   warning(sprintf('removing dimord "%s" from source representation data', data.dimord));
%   data = rmfield(data, 'dimord');
%   return
% else
%   % it is ok
%   return
% end

if ~isfield(data, 'dimord')
    if ea_ft_datatype(data, 'raw')
        % it is raw data, which does not have a dimord -> this is ok
        return
    elseif ea_ft_datatype(data, 'comp')
        % it is component data, which resembles raw data -> this is ok
        return
    elseif ea_ft_datatype(data, 'volume')
        % it is volume data, which does not have a dimord -> this is ok
        return
    else
        fn = fieldnames(data);
        sel = true(size(fn));
        for i=1:length(fn)
            sel(i) = contains(fn{i}, 'dimord');
        end
        df = fn(sel);
        
        if isempty(df)
            if ea_ft_datatype(data, 'source') || ea_ft_datatype(data, 'parcellation')
                % it is old-style source data -> this is ok
                % ft_checkdata will convert it to new-style
                return
            else
                error('the data does not contain a dimord, but it also does not resemble raw or component data');
            end
        end
        
        % use this function recursively on the XXXdimord fields
        for i=1:length(df)
            data.dimord = data.(df{i});
            data = fixdimord(data);
            data.(df{i}) = data.dimord;
            data = rmfield(data, 'dimord');
        end
        % after the recursive call it should be ok
        return
    end
end

if strcmp(data.dimord, 'voxel')
    % this means that it is position
    data.dimord = 'pos';
end

dimtok = tokenize(data.dimord, '_');
if strncmp('{pos_pos}', data.dimord, 9)
    % keep these together for bivariate source structures
    dimtok = {'{pos_pos}', dimtok{3:end}};
end

for i=1:length(dimtok)
    switch dimtok{i}
        case {'tim' 'time' 'toi' 'latency'}
            dimtok{i} = 'time';
            
        case {'frq' 'freq' 'foi' 'frequency'}
            dimtok{i} = 'freq';
            
        case {'sgn' 'label' 'chan'}
            dimtok{i} = 'chan';
            
        case {'rpt' 'trial'}
            dimtok{i} = 'rpt';
            
        case {'subj' 'subject'}
            dimtok{i} = 'subj';
            
        case {'comp'}
            % don't change, it is ok
            
        case {'sgncmb' 'labelcmb' 'chancmb'}
            dimtok{i} = 'chancmb';
            
        case {'rpttap'}
            % this is a 2-D field, coding trials and tapers along the same dimension
            % don't change, it is ok
            
        case {'refchan'}
            % don't change, it is ok
            
        case {'ori'}
            % don't change, it is ok
            
        case {'rgb'}
            % don't change, it is ok
            
        case {'voxel' 'vox' 'repl' 'wcond'}
            % these are used in some fieldtrip functions, but are not considered standard
            warning_once('unexpected dimord "%s"', data.dimord);
            
        case {'pos'}
            % this is for source data on a 3-d grid, a cortical sheet, or unstructured positions
            
        case {'{pos}' '{pos}_rpt' '{pos}_rpttap'}
            % this is for source data on a 3-d grid, a cortical sheet, or unstructured positions
            % the data itself is represented in a cell-array, e.g. source.mom or source.leadfield
            
        case {'{pos_pos}'}
            % this is for bivariate source data on a 3-d grid, a cortical sheet, or unstructured positions
            
        otherwise
            error(sprintf('unexpected dimord "%s"', data.dimord));
            
    end % switch dimtok
end % for length dimtok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data, 'tim'),         data.time      = data.tim         ; data = rmfield(data, 'tim')        ; end
if isfield(data, 'toi'),         data.time      = data.toi         ; data = rmfield(data, 'toi')        ; end
if isfield(data, 'latency'),     data.time      = data.latency     ; data = rmfield(data, 'latency')    ; end
if isfield(data, 'frq'),         data.freq      = data.frq         ; data = rmfield(data, 'frq')        ; end
if isfield(data, 'foi'),         data.freq      = data.foi         ; data = rmfield(data, 'foi')        ; end
if isfield(data, 'frequency'),   data.freq      = data.frequency   ; data = rmfield(data, 'frequency')  ; end
if isfield(data, 'sgn'),         data.label     = data.sgn         ; data = rmfield(data, 'sgn')        ; end
if isfield(data, 'chan'),        data.label     = data.chan        ; data = rmfield(data, 'chan')       ; end
% if isfield(data, 'trial'),         data.rpt     = data.trial         ; data = rmfield(data, 'trial')        ; end  % DO NOT CONVERT -> this is an exception
if isfield(data, 'subject'),     data.subj      = data.subject     ; data = rmfield(data, 'subject')    ; end
if isfield(data, 'sgncmb'),      data.labelcmb  = data.sgncmb      ; data = rmfield(data, 'sgncmb')     ; end
if isfield(data, 'chancmb'),     data.labelcmb  = data.chancmb     ; data = rmfield(data, 'chancmb')    ; end

% ensure that it is a column
if isfield(data, 'label')
    data.label = data.label(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if isfield(data, 'trial')
%   mat = data.trial;
% elseif isfield(data, 'individual')
%   mat = data.individual;
% elseif isfield(data, 'avg')
%   mat = data.avg;
% elseif isfield(data, 'crsspctrm')
%   mat = data.crsspctrm;
% elseif isfield(data, 'powspctrm')
%   mat = data.powspctrm;
% elseif isfield(data, 'fourierspctrm')
%   mat = data.fourierspctrm;
% end
%
% add the descriptive axis for each dimension
% for i=1:length(dimtok)
%   if isfield(data, dimtok{i})
%     % the dimension is already described with its own axis
%     % data = setfield(data, dimtok{i}, getfield(data, dimtok{i}));
%   else
%     % add an axis to the output data
%     data = setfield(data, dimtok{i}, 1:size(mat,i));
%   end
% end

% undo the tokenization
data.dimord = dimtok{1};
for i=2:length(dimtok)
    data.dimord = [data.dimord '_' dimtok{i}];
end

function [type, dimord] = ea_ft_datatype(data, desired)

% ea_ft_datatype determines the type of data represented in a FieldTrip data
% structure and returns a string with raw, freq, timelock source, comp,
% spike, source, volume, dip.
%
% Use as
%   [type, dimord] = ea_ft_datatype(data)
%   [status]       = ea_ft_datatype(data, desired)
%
% See also FT_DATATYPE_COMP FT_DATATYPE_FREQ FT_DATATYPE_MVAR
% ea_ft_datatype_segmentation FT_DATATYPE_PARCELLATION FT_DATATYPE_SOURCE
% FT_DATATYPE_TIMELOCK FT_DATATYPE_DIP FT_DATATYPE_HEADMODEL
% FT_DATATYPE_RAW FT_DATATYPE_SENS FT_DATATYPE_SPIKE FT_DATATYPE_VOLUME

% Copyright (C) 2008-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ea_ft_datatype.m 10064 2014-12-22 14:30:50Z roboos $

if nargin<2
    desired = [];
end

% determine the type of input data, this can be raw, freq, timelock, comp, spike, source, volume, dip, segmentation, parcellation
israw          =  isfield(data, 'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && ~isfield(data,'trialtime');
isfreq         = (isfield(data, 'label') || isfield(data, 'labelcmb')) && isfield(data, 'freq') && ~isfield(data,'trialtime') && ~isfield(data,'origtrial'); %&& (isfield(data, 'powspctrm') || isfield(data, 'crsspctrm') || isfield(data, 'cohspctrm') || isfield(data, 'fourierspctrm') || isfield(data, 'powcovspctrm'));
istimelock     =  isfield(data, 'label') && isfield(data, 'time') && ~isfield(data, 'freq') && ~isfield(data,'timestamp') && ~isfield(data,'trialtime') && ~(isfield(data, 'trial') && iscell(data.trial)); %&& ((isfield(data, 'avg') && isnumeric(data.avg)) || (isfield(data, 'trial') && isnumeric(data.trial) || (isfield(data, 'cov') && isnumeric(data.cov))));
iscomp         =  isfield(data, 'label') && isfield(data, 'topo') || isfield(data, 'topolabel');
isvolume       =  isfield(data, 'transform') && isfield(data, 'dim') && ~isfield(data, 'pos');
issource       =  isfield(data, 'pos');
isdip          =  isfield(data, 'dip');
ismvar         =  isfield(data, 'dimord') && contains(data.dimord, 'lag');
isfreqmvar     =  isfield(data, 'freq') && isfield(data, 'transfer');
ischan         = ea_check_chan(data);
issegmentation = ea_check_segmentation(data);
isparcellation = ea_check_parcellation(data);

if ~isfreq
    % this applies to a freq structure from 2003 up to early 2006
    isfreq = all(isfield(data, {'foi', 'label', 'dimord'})) && contains(data.dimord, 'frq');
end

% check if it is a spike structure
spk_hastimestamp  = isfield(data,'label') && isfield(data, 'timestamp') && isa(data.timestamp, 'cell');
spk_hastrials     = isfield(data,'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && isfield(data, 'trialtime') && isa(data.trialtime, 'numeric');
spk_hasorig       = isfield(data,'origtrial') && isfield(data,'origtime'); % for compatibility
isspike           = isfield(data, 'label') && (spk_hastimestamp || spk_hastrials || spk_hasorig);

% check if it is a sensor array
isgrad = isfield(data, 'label') && isfield(data, 'coilpos') && isfield(data, 'coilori');
iselec = isfield(data, 'label') && isfield(data, 'elecpos');

if isspike
    type = 'spike';
elseif israw && iscomp
    type = 'raw+comp';
elseif istimelock && iscomp
    type = 'timelock+comp';
elseif isfreq && iscomp
    type = 'freq+comp';
elseif israw
    type = 'raw';
elseif iscomp
    type = 'comp';
elseif isfreqmvar
    % freqmvar should conditionally go before freq, otherwise the returned ea_ft_datatype will be freq in the case of frequency mvar data
    type = 'freqmvar';
elseif isfreq
    type = 'freq';
elseif ismvar
    type = 'mvar';
elseif isdip
    % dip should conditionally go before timelock, otherwise the ea_ft_datatype will be timelock
    type = 'dip';
elseif istimelock
    type = 'timelock';
elseif issegmentation
    % a segmentation data structure is a volume data structure, but in general not vice versa
    % segmentation should conditionally go before volume, otherwise the returned ea_ft_datatype will be volume
    type = 'segmentation';
elseif isvolume
    type = 'volume';
elseif isparcellation
    % a parcellation data structure is a source data structure, but in general not vice versa
    % parcellation should conditionally go before source, otherwise the returned ea_ft_datatype will be source
    type = 'parcellation';
elseif issource
    type = 'source';
elseif ischan
    % this results from avgovertime/avgoverfreq after timelockstatistics or freqstatistics
    type = 'chan';
elseif iselec
    type = 'elec';
elseif isgrad
    type = 'grad';
else
    type = 'unknown';
end

if nargin>1
    % return a boolean value
    switch desired
        case 'raw'
            type = any(strcmp(type, {'raw', 'raw+comp'}));
        case 'timelock'
            type = any(strcmp(type, {'timelock', 'timelock+comp'}));
        case 'freq'
            type = any(strcmp(type, {'freq', 'freq+comp'}));
        case 'comp'
            type = any(strcmp(type, {'comp', 'raw+comp', 'timelock+comp', 'freq+comp'}));
        case 'volume'
            type = any(strcmp(type, {'volume', 'segmentation'}));
        case 'source'
            type = any(strcmp(type, {'source', 'parcellation'}));
        case 'sens'
            type = any(strcmp(type, {'elec', 'grad'}));
        otherwise
            type = strcmp(type, desired);
    end % switch
end

if nargout>1
    % FIXME this should be replaced with getdimord in the calling code
    % also return the dimord of the input data
    if isfield(data, 'dimord')
        dimord = data.dimord;
    else
        dimord = 'unknown';
    end
end

function [res] = ea_check_chan(data)

if ~isstruct(data) || any(isfield(data, {'time', 'freq', 'pos', 'dim', 'transform'}))
    res = false;
elseif isfield(data, 'dimord') && any(strcmp(data.dimord, {'chan', 'chan_chan'}))
    res = true;
else
    res = false;
    fn = fieldnames(data);
    for i=1:numel(fn)
        if isfield(data, [fn{i} 'dimord']) && any(strcmp(data.([fn{i} 'dimord']), {'chan', 'chan_chan'}))
            res = true;
            break;
        end
    end
end

function [res] = ea_check_segmentation(volume)
res = false;

if ~isfield(volume, 'dim')
    return
end

if isfield(volume, 'pos')
    return
end

if any(isfield(volume, {'seg', 'csf', 'white', 'gray', 'skull', 'scalp', 'brain'}))
    res = true;
    return
end

fn = fieldnames(volume);
isboolean = [];
cnt = 0;
for i=1:length(fn)
    if isfield(volume, [fn{i} 'label'])
        res = true;
        return
    else
        if (islogical(volume.(fn{i})) || isnumeric(volume.(fn{i}))) && isequal(size(volume.(fn{i})),volume.dim)
            cnt = cnt+1;
            if islogical(volume.(fn{i}))
                isboolean(cnt) = true;
            else
                isboolean(cnt) = false;
            end
        end
    end
end
if ~isempty(isboolean)
    res = all(isboolean);
end

function [res] = ea_check_parcellation(source)
res = false;

if ~isfield(source, 'pos')
    return
end

fn = fieldnames(source);
fb = false(size(fn));
npos = size(source.pos,1);
for i=1:numel(fn)
    % for each of the fields check whether it might be a logical array with the size of the number of sources
    tmp = source.(fn{i});
    fb(i) = numel(tmp)==npos && islogical(tmp);
end
if sum(fb)>1
    % the presence of multiple logical arrays suggests it is a parcellation
    res = true;
end

if res == false      % check if source has more D elements
    check = 0;
    for i = 1: length(fn)
        fname = fn{i};
        switch fname
            case 'tri'
                npos = size(source.tri,1);
                check = 1;
            case 'hex'
                npos = size(source.hex,1);
                check = 1;
            case 'tet'
                npos = size(source.tet,1);
                check = 1;
        end
    end
    if check == 1   % check if elements are labelled
        for i=1:numel(fn)
            tmp = source.(fn{i});
            fb(i) = numel(tmp)==npos && islogical(tmp);
        end
        if sum(fb)>1
            res = true;
        end
    end
end

fn = fieldnames(source);
for i=1:length(fn)
    if isfield(source, [fn{i} 'label']) && isnumeric(source.(fn{i}))
        res = true;
        return
    end
end

function [stiff, diinsy, cols, sysmat] = ea_sb_calc_stiff(vol)

% SB_CALC_STIFF
%
% $Id: sb_calc_stiff.m 8776 2013-11-14 09:04:48Z roboos $

if(~(size(vol.pos,2)==3))
    if(size(vol.pos,1)==3)
        node = vol.pos';
        warning('Dimension of vol.pos should be #nodes x 3!')
    else
        error('vol.pos has wrong dimension!')
    end
else
    node = vol.pos;
end
npnt = size(node,1);
npnt = int32(npnt);

if isfield(vol,'tet')
    if size(vol.tet,1) == 4
        mele = size(vol.tet,1);
        elem = vol.tet;
    elseif size(vol.tet,2) == 4
        mele = size(vol.tet,2);
        elem = vol.tet';
    else
        error('vol.tet has wrong dimensions!')
    end
    elem = [elem; zeros(4,size(elem,2))];
elseif isfield(vol,'hex')
    if size(vol.hex,1) == 8
        mele = size(vol.hex,1);
        elem = vol.hex;
    elseif size(vol.hex,2) == 8
        mele = size(vol.hex,2);
        elem = vol.hex';
    else
        error('vol.hex has wrong dimensions!')
    end
else
    error('Could not find connectivity information!')
end

if min(min(elem(1:mele,:))) == 0
    elem = elem + 1;
    warning('Numbering of nodes in vol.tet/vol.hex must start at 1 (Fortran numbering)!')
elseif min(min(elem(1:mele,:))) < 0
    error('No negative indices for conectivity information allowed!')
end

if isfield(vol,'cond') && isfield(vol,'tissue') && isfield(vol,'tissuelabel')
    if length(vol.tissuelabel) == length(vol.cond)
        if length(vol.tissue) == size(elem,2)
            cond = zeros(size(elem,2),1);
            numlabels = length(vol.tissuelabel);
            for i=1:numlabels
                cond(vol.tissue == i) = vol.cond(i);
            end
        else
            error('Dimensions of vol.tet or vol.hex and vol.tissue do not fit!');
        end
    else
        error('Dimensions of vol.cond and entries of vol.tissuelabel do not fit!');
    end
end

mele = int32(mele);
elem = int32(elem);

% check whether the nodes have right orientation

if isfield(vol,'tet')
    if ~ea_sb_test_ori(node,elem(1:4,:)')
        error('Elements have wrong orientation, consider exchanging node 3 and 4');
        return;
    end
elseif isfield(vol,'hex')
    if ~ea_sb_test_ori(node,elem')
        error('Elements have wrong orientation or are degenerated');
        return
    end
end

try
    [diinsy,cols,sysmat] = ea_calc_stiff_matrix_val_wrapper(node,elem,cond,mele);
    ea_delete([pwd, filesep, 'fort.6']);
catch err
    if ispc && strcmp(err.identifier,'MATLAB:invalidMEXFile')
        error('Error executing mex-file. Microsoft Visual C++ 2008 Redistributables and Intel Visual Fortran Redistributables are required.')
    else
        rethrow(err)
    end
end
npnt = double(npnt);
diinsy = double(diinsy);
cols = double(cols);
rows = ea_sb_sparse_to_mat(diinsy);
stiff = sparse(rows,cols,sysmat,npnt,npnt,length(sysmat));

function err = ea_sb_test_ori(pnt,elem)
err = 1;
if(size(elem,2) == 4)
    det = sum(cross(pnt(elem(:,2),:)-pnt(elem(:,1),:),pnt(elem(:,4),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,3),:)-pnt(elem(:,1),:)),2);
    if length(find(det <= 0)) > 0
        err = 0;
    end
elseif(size(elem,2) == 8)
    det1 = sum(cross(pnt(elem(:,6),:)-pnt(elem(:,1),:),pnt(elem(:,8),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,5),:)-pnt(elem(:,1),:)),2);
    det2 = sum(cross(pnt(elem(:,3),:)-pnt(elem(:,1),:),pnt(elem(:,6),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,2),:)-pnt(elem(:,1),:)),2);
    det3 = sum(cross(pnt(elem(:,8),:)-pnt(elem(:,1),:),pnt(elem(:,3),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,4),:)-pnt(elem(:,1),:)),2);
    det4 = sum(cross(pnt(elem(:,8),:)-pnt(elem(:,3),:),pnt(elem(:,6),:)-pnt(elem(:,3),:),2).*(pnt(elem(:,7),:)-pnt(elem(:,3),:)),2);
    if (length(find(det1 <= 0)) > 0 || length(find(det2 <= 0)) > 0 || length(find(det3 <= 0)) > 0 || length(find(det4 <= 0)) > 0)
        err = 0;
    end
else
    error('Invalid number of nodes per element!');
end

function rows = ea_sb_sparse_to_mat(diinsy)

% SB_SPARSE_TO_MAT
%
% $Id: sb_sparse_to_mat.m 8776 2013-11-14 09:04:48Z roboos $

rows = zeros(max(diinsy),1);
rows(diinsy) = 1;
rows = [1;rows];
rows = cumsum(rows);
rows = rows(1:end-1);
