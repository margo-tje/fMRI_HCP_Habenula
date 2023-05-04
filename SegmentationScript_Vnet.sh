#!/bin/bash

FSLDIR=/opt/amc/fsl-6.0.0/
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH


#Make sure that in this folder, there is only folders with the subjects number
cd {foldername}

ls -d */ > subjects.txt#
subjects=subjects.txt
subjects_nos=`cat $subjects`
 

#VNet segmentation is already runned; only need thresholding and downsampling. The segmentation itself is binary; so transform it to 2mm; then threshold to find the best fit volume wise
######################
#the seed is made, you need to adjust it's size (thresholding)
folder={foldername}
folder2='VNET_Segmentation'
file='/hb_seg'
ext='.nii.gz'
acpc='/MNINonLinear/xfms/acpc_dc2standard.nii.gz'
name='/SeedHb'
thres='_thres.nii.gz'
output='_MNI.nii.gz'
fsl='/opt/amc/fsl-6.0.0/data/standard/MNI152_T1_2mm.nii.gz'
bin='_bin'
save='/Old'
output1='/SeedHb_thresholded'
output2='/SeedHb_MNI'
Voltext1='SeedHb_MNI'
Voltext2='_Volume.txt'
twomm='_2mm.nii.gz'

declare -a side=("_left" "_right") #runs through L and R sides

for sub in $subjects_nos ;  do
mkdir $folder$sub$folder2$save
for sid in ${side[@]};do
repeat=1 
threshold=0.5 #every subject starts with threshold = 0.5
echo $sub$sid

#transform to MNI space and 2mm resolution 
applywarp -i $folder$sub$folder2$file$sid$ext -o $folder$sub$folder2$name$sid$twomm -r /opt/amc/fsl-6.0.0/data/standard/MNI152_T1_2mm.nii.gz -w $folder$sub$acpc
#fslmaths $folder$sub$folder2$name$sid$twomm -bin $folder$sub$folder2$name$sid$bin$twomm #ADDED > can weg?

while [ $repeat -eq 1 ]; do #find the best threshold for 2mm image

#threshold (start with 50% and then 5% steps up or down)
fslmaths $folder$sub$folder2$name$sid$twomm -thr $threshold -bin $folder$sub$folder2$name$sid$threshold$thres

#volume and mean signal 
anatomic_vol=`fslstats $folder$sub$folder2$file$sid$ext -V` 
mni_vol=`fslstats $folder$sub$folder2$name$sid$threshold$thres -V`

anatomic_volume=$(echo $anatomic_vol | cut -d' ' -f 2)     #only take the second volume (mm3), the 							first one is voxels
mni_volume=$(echo $mni_vol | cut -d' ' -f 2) 

#check if MNI volume is within range of target volume (+- 4 mm3)
anatomic_min=$(echo $anatomic_volume -4 | bc -l)  #-4 for min
anatomic_max=$(echo $anatomic_volume +4 | bc -l)  #+4 for max

echo "MNI mask of $sub $sid should be between $anatomic_min and $anatomic_max" 

if (( $(echo "$mni_volume > $anatomic_max" | bc -l) ));then
 echo "MNI Mask is $mni_volume , so it is too big"
 if [ $threshold == 1.0000 ]; then #if threshold = 1, then take that
  echo "Threshold is $threshold"
  echo $mni_volume > $folder$sub$folder2/$Voltext1$sid$Voltext2 # save MNI volume at .text file
  mv $folder$sub$folder2$name$sid$threshold$thres $folder$sub$folder2$output2$sid$ext
  #mv $folder$sub$folder2$name$sid$threshold$output $folder$sub$folder2$output2$sid$ext #rename final 
  if (($sid = '_L'));then
   name1=$(echo $folder$sub$folder2$output2$sid$ext)
  else
   name2=$(echo $folder$sub$folder2$output2$sid$ext)
  fi
 repeat=0
 else #is threshold is other than 1
  mv $folder$sub$folder2$name$sid$threshold$thres 	$folder$sub$folder2$save$name$sid$threshold$thres #move file to 'Old' folder
  threshold=$(echo $threshold+0.000025 | bc -l)
  echo "New threshold is $threshold"
  repeat=1
 fi
elif (( $(echo "$mni_volume < $anatomic_min" | bc -l) ));then
 mv $folder$sub$folder2$name$sid$threshold$thres 	$folder$sub$folder2$save$name$sid$threshold$thres #move file to 'Old' folder
 echo "MNI Mask is $mni_volume , so it is too small"
 threshold=$(echo $threshold-0.000025 | bc -l)
 echo "New threshold is $threshold"
 repeat=1
else
 echo "MNI Mask is $mni_volume , so it is good"
 echo $mni_volume > $folder$sub$folder2$Voltext1$sid$Voltext2 # save MNI volume at .text file
 mv $folder$sub$folder2$name$sid$threshold$thres $folder$sub$folder2$output2$sid$ext
 #mv $folder$sub$folder2$name$sid$threshold$output $folder$sub$folder2$output2$sid$ext #rename the final files 
 if (($sid = '_L'));then
 name1=$(echo $folder$sub$folder2$output2$sid$ext)
 else
 name2=$(echo $folder$sub$folder2$output2$sid$ext)
 fi
 repeat=0
fi
 
done

done

done

#Add left & right together
for sub in $subjects_nos ;  do
cd $folder$sub$folder2
fslmaths SeedHb_MNI_left -add SeedHb_MNI_right SeedHb_MNI_Bi #make bilateral again
done




