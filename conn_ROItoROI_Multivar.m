% CONN_BATCH_HUMANCONNECTOMEPROJECT batch processing script for the Human Connectome Project resting-state dataset (HCP; http://www.humanconnectome.org/; tested on Q6 497-subjects release)

% This script assumes that the dataset is already available in your system (set the CONNECTOMEpath variable in this script to the appropriate folder)

% The script will create a new conn_HCP project and:
%  Load the preprocessed "clean" ICA+FIX functional series (four sessions)
%   Apply structural segmentation, functional ART scrubbing, and functional smoothing
%   Apply default aCompCor denoising and band-pass filtering
%   Compute seed-to-voxel and ROI-to-ROI bivariate correlation connectivity measures for all CONN default ROIs (atlas + dmn)

% Default settings (see code for details):
%    Run each individual subject in a separate parallel stream (edit RUNPARALLEL and NJOBS variables in this script to modify these settings)
%    Assumes user has write permission into connectome data folders (edit COPYFILES variable in this script when users have only read-permissions)
%    Run all subjects in dataset (edit NSUBJECTS in this script to process instead only a subset of subjects)
%
%   Margo: How I use this code, I let it copy 1 subject from the L Disk to
%   Scratch, let it run and copy everything back. This is way faster
% note: before running this script it is recommended to test it first with just a few subjects to make sure everything is working
% as expected. To do so, set NSUBJECTS=4; in the DEFAULT SETTINGS section
% of this script

%% DEFAULT SETTINGS: EDIT THE LINES BELOW (minimally set CONNECTOMEpath to the actual location of your connectome data)

addpath('/opt/amc/matlab/toolbox/conn-18b')
addpath('/opt/amc/matlab/toolbox/spm12')
TARGETpath='/scratch/mslomp/Conn';        % target folder for conn project (default current folder) % ! ! try to have this also as current folder!!
CONNECTOMEpath={path}; %this is the location where your data is stored
SAVEPATH={path}; %on the L disk (storage) %Change this to your folder
RUNPARALLEL=false; %used to be false                                                           % run in parallel using computer cluster
NSUBJECTS=[]; %used to be []                                                 % number of subjects to include in your project (leave empty for all subjects)
NJOBS=[50];                                                                   % number of parallel jobs to submit (leave empty for one job per subject)
COPYFILES=true;                                                            % true/false: set to true if you do not have write-permissions into connectome data folders
                                                                            %   This will create a local copy (in same folder as your conn project) of the structural/functional data
                                                                            %   where any post-processed files will also be store
OVERWRITE=false;                                                            % overwrites files if they exist in target folder (unzipped files and/or files in local-copy folder)
                                                                            %   Set to false if you have already unzipped / copied-to-local-folder your data and would like to skip this step

%% FINDS STRUCTURAL/FUNCTIONAL/REALIGNMENT FILES
tic

subs=dir(regexprep([CONNECTOMEpath '/ProcessedData/'],'%s.*$','*')); %%  < Other wise, run:
subs=subs([subs.isdir]>0);
subs={subs.name};
subs=subs(cellfun(@(s)all(s>='0'&s<='9'),subs)); %now 'subs' is a cell array with all subject numbers present in the data folder

NSUBJECTS=1; %run one subject at the time
if isempty(NJOBS), NJOBS=NSUBJECTS; end
NJOBS=min(NSUBJECTS,NJOBS);
if isempty(dir(fullfile(TARGETpath,'/Results_ROItoROI_72pps')))
mkdir(fullfile(TARGETpath,'/Results_ROItoROI_72pps'));
end

%n=1
for n=1:numel(subs)%n%=1
tic
clear FUNCTIONAL_FILE* REALIGNMENT_FILE* STRUCTURAL_FILE* HABENULA_SEED* CONTROL_CM* CONTROL_DM* FMC_SEED* NAC_SEED* AMY_SEED* HPC_SEED* INS_SEED* HYP_SEED* VTA_SEED*;
name = sprintf('connHCP_%s',subs{n});
if isempty(dir(fullfile(TARGETpath,name)))
mkdir(fullfile(TARGETpath,name));
end
fprintf('Locating subject %s files\n',subs{n});
    t1=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear','T1w_restore.2.nii.gz');                                        % STRUCTURAL VOLUME
    f1=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear/Results/','rfMRI_REST1_LR','rfMRI_REST1_LR_hp2000_clean.nii.gz');   % FUNCTIONAL VOLUME (1/4)
    f2=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear/Results/','rfMRI_REST2_LR','rfMRI_REST2_LR_hp2000_clean.nii.gz');   % FUNCTIONAL VOLUME (2/4)
    f3=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear/Results/','rfMRI_REST1_RL','rfMRI_REST1_RL_hp2000_clean.nii.gz');   % FUNCTIONAL VOLUME (3/4)
    f4=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear/Results/','rfMRI_REST2_RL','rfMRI_REST2_RL_hp2000_clean.nii.gz');   % FUNCTIONAL VOLUME (4/4)
    r1=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear/Results/','rfMRI_REST1_LR','Movement_Regressors_dt.txt');              % REALIGNMENT FILE (1/4)
    r2=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear/Results/','rfMRI_REST2_LR','Movement_Regressors_dt.txt');              % REALIGNMENT FILE (2/4)
    r3=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear/Results/','rfMRI_REST1_RL','Movement_Regressors_dt.txt');              % REALIGNMENT FILE (3/4)
    r4=fullfile(CONNECTOMEpath,subs{n},'MNINonLinear/Results/','rfMRI_REST2_RL','Movement_Regressors_dt.txt');              % REALIGNMENT FILE (4/4)
    s1=fullfile(CONNECTOMEpath,subs{n},'VNET_Segmentation/SeedHb_MNI_Bi.nii.gz'); %load in the seed 
   

%Load the ROIs
    s2=fullfile(TARGETpath, '/ROIs/Cluster_pgACC_Bi.nii.gz'); 
    s3=fullfile(TARGETpath, '/ROIs/50Threshold_Insula_Bilateral.nii.gz');
    s4=fullfile(TARGETpath, 'ROIs/50Threshold_Bi_Accumbens.nii.gz');
    s5=fullfile(TARGETpath, '/ROIs/50Threshold_Bi_Amygdala.nii.gz');
    s6=fullfile(TARGETpath, '/ROIs/50Threshold_Bi_Hippocampus.nii.gz');
    s7=fullfile(CONNECTOMEpath,subs{n},'HypothalamusSegmentation/Hypo_Seg_MNI_Vol0.nii.gz'); %load in  the hpothalamus seed
    s8=fullfile(TARGETpath, '/ROIs/0.075Threshold_Bi_VTA.nii.gz'); % VTA

    cT1=fullfile(TARGETpath,'ControlSeeds/72subs_HbA1c/Ctr_CM_Bi.nii'); % load in the control seed (bilateral, CM) 
    cT2=fullfile(TARGETpath,'ControlSeeds/72subs_HbA1c/Ctr_DM_Bi.nii'); % load in the control seed (bilateral, DM)
    
    if isempty(dir(t1)), error('file %s not found',t1); end
    if isempty(dir(f1)), error('file %s not found',f1); end
    if isempty(dir(f2)), error('file %s not found',f2); end
    if isempty(dir(r1)), error('file %s not found',r1); end
    if isempty(dir(r2)), error('file %s not found',r2); end
    if isempty(dir(s1)), error('file %s not found',s1); end 
    if isempty(dir(s2)), error('file %s not found',s2); end 
    if isempty(dir(s3)), error('file %s not found',s3); end 
    if isempty(dir(s4)), error('file %s not found',s4); end 
    if isempty(dir(s6)), error('file %s not found',s6); end 
    if isempty(dir(s7)), error('file %s not found',s7); end 
    if isempty(dir(s8)), error('file %s not found',s8); end %
    if isempty(dir(cT1)), error('file %s not found',cT1); end %control thalamus seed 
    if isempty(dir(cT2)), error('file %s not found',cT2); end %control thalamus seed 

    if COPYFILES

        if ~isempty(dir(fullfile([SAVEPATH 'LocalCopyDataFiles'],subs{n}))) %if conn alreadu hase files, move to next step
            if isempty(dir(fullfile(TARGETpath,'LocalCopyDataFiles',subs{n})))
            disp('Folder already excists, moving if from L disk to Scratch...')
            movefile(fullfile([SAVEPATH 'LocalCopyDataFiles'],subs{n}),fullfile(TARGETpath,'LocalCopyDataFiles') ); %move the Local Copy Files from this subject back to the L disk
            else
            disp('Files already present in Scratch folder, moving on...')
            end
        else        
        fprintf('Copying files to local folder\n');
        [ok,nill]=mkdir(TARGETpath,'LocalCopyDataFiles');
        [ok,nill]=mkdir(fullfile(TARGETpath,'LocalCopyDataFiles'),subs{n});
        end

        t1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'structural.nii.gz');  if OVERWRITE||isempty(dir(t1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',t1,t1b)); end; t1=t1b;
        f1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional1.nii.gz'); if OVERWRITE||isempty(dir(f1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f1,f1b)); end; f1=f1b;
        f2b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional2.nii.gz'); if OVERWRITE||isempty(dir(f2b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f2,f2b)); end; f2=f2b;
        f3b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional3.nii.gz'); if OVERWRITE||isempty(dir(f3b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f3,f3b)); end; f3=f3b;
        f4b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional4.nii.gz'); if OVERWRITE||isempty(dir(f4b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f4,f4b)); end; f4=f4b;
        r1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Movement1.txt'); if OVERWRITE||isempty(dir(r1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r1,r1b)); end; r1=r1b;
        r2b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Movement2.txt'); if OVERWRITE||isempty(dir(r2b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r2,r2b)); end; r2=r2b;
        r3b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Movement3.txt'); if OVERWRITE||isempty(dir(r3b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r3,r3b)); end; r3=r3b;
        r4b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Movement4.txt'); if OVERWRITE||isempty(dir(r4b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r4,r4b)); end; r4=r4b;
        s1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'SeedHb.nii.gz'); if OVERWRITE||isempty(dir(s1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',s1,s1b)); end; s1=s1b;
        s2b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Cluster_pgACC_Bi.nii.gz'); if OVERWRITE||isempty(dir(s2b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',s2,s2b)); end; s2=s2b;       
        s3b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Insula.nii.gz'); if OVERWRITE||isempty(dir(s3b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',s3,s3b)); end; s3=s3b; 
        s4b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'NAc.nii.gz'); if OVERWRITE||isempty(dir(s4b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',s4,s4b)); end; s4=s4b; 
        s5b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n}, 'AMY.nii.gz'); if OVERWRITE||isempty(dir(s5b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',s5,s5b)); end; s5=s5b; 
        s6b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n}, 'HPC.nii.gz'); if OVERWRITE||isempty(dir(s6b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',s6,s6b)); end; s6=s6b; 
        s7b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n}, 'HYP.nii.gz'); if OVERWRITE||isempty(dir(s7b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',s7,s7b)); end; s7=s7b; 
        s8b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n}, 'VTA.nii.gz'); if OVERWRITE||isempty(dir(s8b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',s8,s8b)); end; s6=s6b;   
        cT1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Ctr_Seed_CM.nii'); if OVERWRITE||isempty(dir(cT1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',cT1,cT1b)); end; cT1=cT1b;
        cT2b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Ctr_Seed_DM.nii'); if OVERWRITE||isempty(dir(cT2b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',cT2,cT2b)); end; cT2=cT2b;
        end
toc
    fprintf('Unzipping files\n');
    if OVERWRITE||isempty(dir(regexprep(t1,'\.gz$',''))), gunzip(t1); end; t1=regexprep(t1,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(f1,'\.gz$',''))), gunzip(f1); end; f1=regexprep(f1,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(f2,'\.gz$',''))), gunzip(f2); end; f2=regexprep(f2,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(f3,'\.gz$',''))), gunzip(f3); end; f3=regexprep(f3,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(f4,'\.gz$',''))), gunzip(f4); end; f4=regexprep(f4,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(s1,'\.gz$',''))), gunzip(s1); end; s1=regexprep(s1,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(s2,'\.gz$',''))), gunzip(s2); end; s2=regexprep(s2,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(s3,'\.gz$',''))), gunzip(s3); end; s3=regexprep(s3,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(s4,'\.gz$',''))), gunzip(s4); end; s4=regexprep(s4,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(s5,'\.gz$',''))), gunzip(s5); end; s5=regexprep(s5,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(s6,'\.gz$',''))), gunzip(s6); end; s6=regexprep(s6,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(s7,'\.gz$',''))), gunzip(s7); end; s7=regexprep(s7,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(s8,'\.gz$',''))), gunzip(s8); end; s8=regexprep(s8,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(cT1,'\.gz$',''))), gunzip(cT1); end; cT1=regexprep(cT1,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(cT2,'\.gz$',''))), gunzip(cT2); end; cT2=regexprep(cT2,'\.gz$','');
    r1b=regexprep(r1,'\.txt$','\.deg.txt'); if OVERWRITE||isempty(dir(r1b)), [ok,nill]=system(sprintf('mv ''%s'' ''%s''',r1,r1b)); end; r1=r1b; % note: angles in degrees
    r2b=regexprep(r2,'\.txt$','\.deg.txt'); if OVERWRITE||isempty(dir(r2b)), [ok,nill]=system(sprintf('mv ''%s'' ''%s''',r2,r2b)); end; r2=r2b; 
    r3b=regexprep(r3,'\.txt$','\.deg.txt'); if OVERWRITE||isempty(dir(r3b)), [ok,nill]=system(sprintf('mv ''%s'' ''%s''',r3,r3b)); end; r3=r3b; 
    r4b=regexprep(r4,'\.txt$','\.deg.txt'); if OVERWRITE||isempty(dir(r4b)), [ok,nill]=system(sprintf('mv ''%s'' ''%s''',r4,r4b)); end; r4=r4b; 

    %row is always 1 cus running 1 subject at the time
    STRUCTURAL_FILE{1,1}=t1;
    FUNCTIONAL_FILE(1,1:4)={f1,f2,f3,f4};
    REALIGNMENT_FILE(1,1:4)={r1,r2,r3,r4};
    HABENULA_SEED{1,1}=s1;
    FMC_SEED{1,1}=s2;
    INS_SEED{1,1}=s3;
    NAC_SEED{1,1}=s4;
    AMY_SEED{1,1}=s5;
    HPC_SEED{1,1}=s6; 
    HYP_SEED{1,1}=s7;
    VTA_SEED{1,1}=s8;
    CONTROL_CM(1,1)={cT1};
    CONTROL_DM(1,1)={cT2};
    
nsessions=4;
fprintf('%d subjects, %d sessions\n',NSUBJECTS,nsessions);
toc
%% CREATES CONN BATCH STRUCTURE

clear batch;
batch.filename=fullfile(TARGETpath,[name '.mat']); % 
if RUNPARALLEL

batch.parallel.N=NJOBS;                             % number of parallel processing batch jobs
batch.parallel.profile ='Slurm';
end                                                     % note: use default parallel profile (defined in GUI Tools.GridSettings)

% CONN Setup                                           
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS;
batch.Setup.RT=0.72;                                    % TR (seconds)
batch.Setup.conditions.names={'rest'};                  % single condition (aggregate across all sessions)
for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions,      batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);     % Point to functional volumes for each subject/session
for nsub=1:NSUBJECTS,for nses=1:nsessions,                  batch.Setup.functionals{nsub}{nses}=FUNCTIONAL_FILE(nsub,nses);end; end 
batch.Setup.structurals=STRUCTURAL_FILE;                % Point to anatomical volumes for each subject

batch.Setup.voxelresolution=1;                          % default 2mm isotropic voxels analysis space

batch.Setup.covariates.names={'realignment' };
batch.Setup.covariates.files{1}=repmat({{}},[NSUBJECTS,1]);      
for nsub=1:NSUBJECTS,for nses=1:nsessions,                  batch.Setup.covariates.files{1}{nsub}{nses}=REALIGNMENT_FILE(nsub,nses);end; end 
batch.Setup.rois.names={'SeedHb', 'FrontalMedialCortex', 'Insula', 'NAc', 'AMY', 'HPC','HYP', 'VTA', 'Ctr_DM', 'Ctr_CM';}; %

batch.Setup.rois.files=repmat({{}},[NSUBJECTS,4]);     % Point to habenula Seeds for each subject
for nsub=1:NSUBJECTS,   batch.Setup.rois.files(nsub,:)=[HABENULA_SEED(nsub,1),FMC_SEED(nsub,1),INS_SEED(nsub,1),NAC_SEED(nsub,1), AMY_SEED(nsub,1), HPC_SEED(nsub,1),HYP_SEED(nsub,1), VTA_SEED(nsub,1),CONTROL_DM(nsub,1),CONTROL_CM(nsub,1)];  end ;
                   
batch.Setup.rois.mask(1:10) = 1; %mask rois with GM

batch.Setup.analyses=[1];                             % ROI-to-ROI 
batch.Setup.overwrite='Yes';                            % 1: ROI-to-ROI; 2: Seed-to-voxel; 3: Voxel-to-voxel; 4: Dynamic FC)                            
batch.Setup.done=1;

batch.Setup.preprocessing.steps={'structural_segment','functional_art'};  % Run additional preprocessing steps: segmentation,ART,smoothing

% CONN Denoising                                    
batch.Denoising.filter=[0.01, 0.10];                    % frequency filter (band-pass values, in Hz)
batch.Denoising.confounds.names={'realignment','scrubbing', 'White Matter', 'CSF' ,'Effect of rest'}; %use motion regressors, white matter and csf ase confound
batch.Denoising.done=1; %make 1 to perform              % use default denoising step
batch.Denoising.overwrite='Yes';

% CONN Analysis                                         % Default options (uses all ROIs in conn/rois/ as connectivity sources); see conn_batch for additional options 
batch.Analysis.done=1;
batch.Analysis.overwrite='Yes';
batch.Analysis.type=1; 
batch.Analysis.sources={'SeedHb', 'Ctr_DM', 'Ctr_CM'}; 
          %use Hb and thalamus control rois as seeds

          %Changed to 3 ipv 4: bivariate regression ipv multivariate
batch.Analysis.measure=4; %    1 = 'correlation (bivariate), 2 = 'correlation (semipartial)', 3 =regression (bivariate), 4 = 'regression (multivariate)';                           %correlation (bivariate) = 1, multivariate regression (4)
batch.Analysis.outputfiles=[0,1,1,0,1,0]; %

%%Quality Assurance plots (QA)(optional)

%%batch.QA.plots = {'QA_NORM functional','QA_NORM rois', 'QA_REG functional','QA_DENOISE histogram'};

%% RUNS CONN BATCH STRUCTURE
toc
conn_batch(batch);
toc
fprintf('Done running Conn batch for subject %s... \n ',subs{n});



%% Move results to subject specific folder

if isdir([SAVEPATH '/Results_ROItoROI_72pps/connHCP_' subs{n}]); rmdir([SAVEPATH '/Results_ROItoROI_72pps/connHCP_' subs{n} '/'],'s');end; %is there is already a folder excisting, delete

movefile(fullfile(TARGETpath,name), [SAVEPATH '/Results_ROItoROI_72pps']);

movefile(fullfile(TARGETpath,[name '.mat']), [SAVEPATH '/Results_ROItoROI_72pps']); %move .mat file


%% Move the LocalCopyDataFiles subject folder back to the L disk

fprintf('Moving LocalCopy files back to the L disk \n')

if ~isempty(fullfile([SAVEPATH 'LocalCopyDataFiles'], subs{n}))

delete(fullfile([SAVEPATH 'LocalCopyDataFiles'], subs{n}, '*'));    

end

movefile(fullfile(TARGETpath,'LocalCopyDataFiles',subs{n}), [SAVEPATH 'LocalCopyDataFiles']); %move the Local Copy Files from this subject back to the L disk



fprintf('Done, moving on to next subject... \n')



end %loop through all subjects






%% CONN Display

%conn

%conn('load',fullfile(pwd,'conn_HCP.mat'));

%conn gui_results


