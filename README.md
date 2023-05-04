# fMRI_HCP_Habenula
Scripts used in "Investigating Habenula Functional Connectivity and Reward-Related Activity in Obesity using Human Connectome Project Data" 

Authors: M.Slomp 1,2,3,4, I.G.S. de Lange 1,4, J.D.Mul 2,5, A.Schrantee 2,6*, S.E.la Fleur 1,2,3,4*
1. Amsterdam UMC location University of Amsterdam, Endocrinology Laboratory, Dept Laboratory Medicine, Meibergdreef 9, Amsterdam, The Netherlands
2. Amsterdam Neuroscience, Amsterdam, The Netherlands
3. Amsterdam Gastroenterology Endocrinology & Metabolism, Amsterdam, The Netherlands
4. Metabolism and Reward Group, Royal Netherlands Academy of Arts and Sciences, Netherlands Institute of Neuroscience, Amsterdam, Netherlands. 
5. Brain Plasticity group, Swammerdam Institute for Life Sciences, Faculty of Science, University of Amsterdam, Amsterdam, The Netherlands.
6. Department of Radiology and Nuclear Medicine, Amsterdam UMC, University of Amsterdam, Meibergdreef 9, Amsterdam, Netherlands. 
*contributed equally

Contact info: m.slomp@amsterdamumc.nl

## Data used: HCP S1200 Release.
You can downdload this dataset on: https://db.humanconnectome.org/app/template/Login.vm

## Habenula segmentation: 
Used V-Net automated habenula segmentation using only T1w scans from Kim, Joo-won, and Junqian Xu. "Automated Human Habenula Segmentation from T1-weighted Magnetic Resonance Images using V-Net." bioRxiv (2022): 2022-01, see https://github.com/joowon-kim/hb_seg_vnet

Then, resulting habenula seeds were thresholded as in Ely, Benjamin A., et al. "Resting‐state functional connectivity of the human habenula in healthy individuals: Associations with subclinical depression." Human brain mapping 37.7 (2016): 2369-2384.

   Used script: Segmentationscript_Vnet.sh

## Running CONN
When all seeds and ROIs are achieved, you can run CONN for the fMRI analysis. I adapted this script: https://github.com/alfnie/conn/blob/master/conn_batch_humanconnectomeproject.m

  Used script: conn_ROItoROI_Multivar.m
  Slightly adapted per dataset
  
## Resting state fMRI: 2nd level analysis in R
I combine all output from CONN in 1 excel as input to R for 2nd level analysis.
   Used script: RScript_ROItoROI_300pps-BMI-RestingState_clean.R
   
  *for the exploratory dataset (n=72), the same is used but different input, and BMI * HbA1c instead of BMI + HbA1c.

## Task fMRI: 2nd level analysis in R
  Used script: R_MixedModel_Script_tfMRI_final.R






