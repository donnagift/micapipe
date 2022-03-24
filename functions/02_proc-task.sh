#!/bin/bash
# Task fMRI processing adapted and modified from the resting state processing module (rs module details below)
# Written by Donna Gift Cabalo March 2022
#Resting-state module information
  # Written by Casey Paquola and Reinder Vos De Wael (Oct 2018).
  # and a tiny bit from Sara (Feb 2019)...
  # and a whole lot from Sara (August 2019)
  # and incorporation to mica-pipe by Raul (August-September 2020)
  # and addition of a bunch of fancy flags by Jessica (October-November 2020)

umask 003
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
changeTopupConfig=$8
changeIcaFixTraining=$9
thisMainScan=${10}
thisPhase=${11}
smooth=${12}
taskScanStr=${13}
fmri_pe=${14}
fmri_rpe=${15}
performNSR=${16}
performGSR=${17}
noFIX=${18}
sesAnat=${19}
regAffine=${20}
fmri_acq=${21}
dropTR=${22}
PROC=${23}
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

if [[ "$sesAnat" != FALSE  ]]; then
  sesAnat=${sesAnat/ses-/}
  BIDSanat="${subject}_ses-${sesAnat}"
  dir_anat="${out}/${subject}/ses-${sesAnat}/anat"
  dir_volum="${dir_anat}/volumetric"
  dir_conte69="${dir_anat}/surfaces/conte69"
  T1nativepro="${dir_anat}/${BIDSanat}_space-nativepro_t1w.nii.gz"
  T1nativepro_brain="${dir_anat}/${BIDSanat}_space-nativepro_t1w_brain.nii.gz"
  T1nativepro_mask="${dir_anat}/${BIDSanat}_space-nativepro_t1w_brain_mask.nii.gz"
  dir_freesurfer="${dir_surf}/${subject}_ses-${sesAnat}"
  T1freesurfr="${dir_freesurfer}/mri/T1.mgz"
else
  BIDSanat="${idBIDS}"
  dir_anat="${proc_struct}"
fi
T1_seg_subcortex="${dir_volum}/${BIDSanat}_space-nativepro_t1w_atlas-subcortical.nii.gz"
T1_seg_cerebellum="${dir_volum}/${BIDSanat}_space-nativepro_t1w_atlas-cerebellum.nii.gz"

### CHECK INPUTS: task fmri , phase encoding, structural proc, topup and ICA-FIX files
Info "Inputs:"
Note "Topup Config     :" "$changeTopupConfig"
Note "ICA fix training :" "$changeIcaFixTraining"
Note "fMRI task scan   :" "$taskScanStr"
Note "Phase scan       :" "$fmri_pe"
Note "Reverse Phase    :" "$fmri_rpe"
Note "Smoothing        :" "$smooth"
Note "Perform NSR      :" "$performNSR"
Note "Perform GSR      :" "$performGSR"
Note "No FIX           :" "$noFIX"
Note "Longitudinal ses :" "$sesAnat"
Note "fmri acq         :" "$fmri_acq"

#------------------------------------------------------------------------------#
if [[ "$taskScanStr" == FALSE ]]; then
  Error "No task fMRI string was provided"; exit;
else
    Info "Using user provided main scan: ${subject_bids}/func/${idBIDS}_${taskScanStr}"
    mainScan=$(ls "${subject_bids}/func/${idBIDS}_${taskScanStr}".nii* 2>/dev/null)
    taskScanJson=$(ls "${subject_bids}/func/${idBIDS}_${taskScanStr}".json 2>/dev/null)

fi
# If no json is found search at the top BIDS directory
if [[ -z ${taskScanJson} ]]; then taskScanJson="${BIDS}/task-rest_bold.json"; fi

#------------------------------------------------------------------------------#
# Phase encoding
N_mainPhase=${#bids_mainPhase[@]}
N_revPhase=${#bids_reversePhase[@]}
if [ "$N_mainPhase" -gt 1 ] || [ "$N_revPhase" -gt 1 ]; then
    if [[ "$thisPhase" == "DEFAULT" ]]; then
        Error "Found multiple phase reversal runs in BIDS rawdata directory! Please specify which run should be processed using flag -phaseReversalRun"; exit;
    elif [ "$thisPhase" -gt "$N_mainPhase" ] || [ "$thisPhase" -gt "$N_revPhase" ]; then
        Warning "Specified run number ($thisPhase) is greater than number of phase reversal scans scans found ($N_mainPhase and $N_revPhase). Using first filename in list as default";
        mainPhaseScan=${bids_mainPhase[$thisPhase-1]}
        reversePhaseScan=${bids_reversePhase[$thisPhase-1]}
    else
        Info "Found $N_mainPhase and $N_revPhase phase reversal scans, processing specified scan # $thisPhase"
        mainPhaseScan=${bids_mainPhase[$thisPhase-1]}
        reversePhaseScan=${bids_reversePhase[$thisPhase-1]}
    fi
else
    mainPhaseScan=${bids_mainPhase[0]}
    reversePhaseScan=${bids_reversePhase[0]}
    if [[ "$thisPhase" == "DEFAULT" ]]; then
        Info "No run number specified for phase reversals and did not find more than one phase reversal scan - all good!"
    else
        if [ "$thisPhase" -gt "$N_mainPhase" ] || [ "$thisPhase" -gt "$N_revPhase" ]; then
            Warning "Specified run number ($thisPhase) is greater than number of phase reversal scans scans found ($N_mainPhase and $N_revPhase). Using first filename in list as default"; fi
    fi
fi

# Manually defined Phase scan and reverse phase scan
if [[ "$fmri_pe" != DEFAULT ]] && [[ -f "$fmri_pe" ]]; then mainPhaseScan="$fmri_pe"; fi
if [[ "$fmri_rpe" != DEFAULT ]] && [[ -f "$fmri_rpe" ]]; then reversePhaseScan="$fmri_rpe"; fi

# Check inputs
if [ ! -f "$mainScan" ]; then Error "Couldn't find $id main task fMRI scan : \n\t ls ${mainScan} Warning: If file exist in the directory, check if input for taskScanStr is spelled correctly"; exit; fi #Last check to make sure file exists
if [ ! -f "$taskScanJson" ]; then Error "Couldn't find $id main task fMRI scan json file: \n\t ls ${taskScanJson}"; exit; fi #Last check to make sure file exists
if [ -z "$mainPhaseScan" ]; then  Warning "Subject $id doesn't have acq-APse_bold: TOPUP will be skipped"; fi
if [ -z "$reversePhaseScan" ]; then Warning "Subject $id doesn't have acq-PAse_bold: TOPUP will be skipped"; fi

# Check requirements: Structural nativepro scan and freesurfer, and post_structural
if [ ! -f "$T1nativepro" ]; then Error "Subject $id doesn't have T1_nativepro: run -proc_structural"; exit; fi
if [ ! -f "$T1freesurfr" ]; then Error "Subject $id doesn't have a T1 in freesurfer space: <SUBJECTS_DIR>/${idBIDS}/mri/T1.mgz"; exit; fi
if [ ! -f "$T1_seg_cerebellum" ]; then Error "Subject $id doesn't have cerebellar segmentation:\n\t\t ls ${T1_seg_cerebellum} \n\t\tRUN -post_structural"; exit; fi
if [ ! -f "$T1_seg_subcortex" ]; then Error "Subject $id doesn't have subcortical segmentation:\n\t\t ls ${T1_seg_subcortex} \n\t\t -post_structural"; exit; fi

# Check topup input
if [[ ${changeTopupConfig} == "DEFAULT" ]]; then
    Info "Will use default config file for TOPUP: ${topupConfigFile}"
else
    topupConfigFile=${changeTopupConfig}
    Info "Will use specified config file for TOPUP: ${topupConfigFile}"
fi

# Check FIX: run or no?
if [[ "$noFIX" -eq 1 ]]; then
    Info "ICA-FIX will be skipped! Consider performing nuisance signal regression with <-regress_WM_CSF> or <-GSR>"

    # Check ICA-FIX Training input
    if [[ ! ${changeIcaFixTraining} == "DEFAULT" ]]; then
        Error "If ICA-FIX is skipped, <-icafixTraining> must remain empty"; exit; fi
else
    Info "ICA-FIX pipeline will be run!"

    # Check ICA-FIX Training input
    if [[ ${changeIcaFixTraining} == "DEFAULT" ]]; then
        Info "Will use default training file for ICA-FIX: ${icafixTraining}"
    else
        icafixTraining=${changeIcaFixTraining}
        Info "Will use specified training file for ICA-FIX: ${icafixTraining}"
    fi
fi

# Check smoothing
if [[ $smooth == 1 ]]; then
    Info "Smoothing of native surface timeseries will be performed using workbench command"
else
    Info "Smoothing of native surface timeseries will be performed using FreeSurfer tools (default)"
fi

# Check nuisance signal regression
if [[ $performNSR == 1 ]]; then
    Info "White matter and CSF signals will be regressed from processed timeseries"
elif [[ $performGSR == 1 ]]; then
    Info "Global, white matter and CSF signals will be regressed from processed timeseries"
else
    Info "Global, white matter and CSF signal regression will not be performed (default)"
fi

# gettin dat from taskScanJson exit if Not found
readoutTime=$(grep TotalReadoutTime "${taskScanJson}" | awk -F ' ' '{print $2}' | awk -F ',' '{print $1}')
RepetitionTime=$(grep RepetitionTime "${taskScanJson}" | awk -F ' ' '{print $2}' | awk -F ',' '{print $1}')
if [[ -z "$readoutTime" ]]; then Warning "readoutTime is missing in $taskScanJson, if TOPUP was selected it will likely FAIL"; fi
if [[ -z "$RepetitionTime" ]]; then Error "RepetitionTime is missing in $taskScanJson $RepetitionTime"; exit; fi

#------------------------------------------------------------------------------#
Title "Task fMRI processing\n\t\tmicapipe $Version, $PROC "
micapipe_software
bids_print.variables-taskfmri
Info "Saving temporal dir: $nocleanup"
Info "ANTs will use $threads threads"
Info "wb_command will use $OMP_NUM_THREADS threads"
# taskfmri directories
if [[ ${fmri_acq} == "TRUE" ]]; then
  fmri_tag=$(echo $mainScan | awk -F ${idBIDS}_ '{print $2}' | cut -d'.' -f1); fmri_tag=${fmri_tag/_bold/}
  tagMRI="${proc_taskfmri}"
  proc_taskfmri="$subject_dir/func/${fmri_tag}"
  Info "Outputs will be stored in:"
  Note "fMRI path:" "${proc_taskfmri}"

else
  tagMRI="taskfmri"
  proc_taskfmri="$subject_dir/func/${fmri_tag}"
fi
Note "tagMRI:      " "${tagMRI}"
#	Timer
aloita=$(date +%s)
Nsteps=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_proc-taskfmri_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Define directories
export SUBJECTS_DIR="$dir_surf"
taskfmri_volum="${proc_taskfmri}/$taskScanStr/volumetric"   # volumetricOutputDirectory
taskfmri_surf="${proc_taskfmri}/$taskScanStr/surfaces"      # surfaceOutputDirectory
taskfmri_ICA="$proc_taskfmri/$taskScanStr/ICA_MELODIC"      # ICAOutputDirectory

# Make directories - exit if processing directory already exists (to prevent deletion of existing files at the end of this script).
for x in "$taskfmri_surf" "$taskfmri_volum"; do
    [[ ! -d "${x}" ]] && mkdir -p "${x}"
done

#------------------------------------------------------------------------------#
# Begining of the REAL processing
# Scans to process
toProcess=($mainScan $reversePhaseScan $mainPhaseScan)
tags=(mainScan reversePhaseScan mainPhaseScan)
singleecho="${taskfmri_volum}/${idBIDS}"_space-taskfmri_desc-singleecho.nii.gz

# scan for registration to T1
# taskfmri4reg="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_mainPhaseAlignedTopup_mean.nii.gz"

# Processing single.
if [[ ! -f "${singleecho}" ]]; then
    # Loop over all scans for everything before motion correction across scans.
    for i in {0..2}; do
        # Get basic parameters
        rawNifti=${toProcess[$i]}
        tag=${tags[$i]}
        Info "Processing TAG-$tag scan, readout time: $readoutTime ms"
        # IF FILE NOT FOUND DON'T RUN
        if [[ ! -z "${rawNifti}" ]] && [[ -f "${rawNifti}" ]]; then
              Note "RAWNIFTI:" "$rawNifti"

              # Drop first five TRs and reorient to standard
              if [ "$tag" == "mainScan" ] && [ "$dropTR" == "TRUE" ]; then
                  Do_cmd nifti_tool -cbl -prefix "${tmp}/${tag}_trDrop.nii.gz" -infiles "$rawNifti"'[5..$]'
                  Do_cmd 3dresample -orient LPI -prefix "${tmp}/${tag}_reorient.nii.gz" -inset "${tmp}/${tag}_trDrop.nii.gz"
                  Do_cmd fslreorient2std "${tmp}/${tag}_reorient.nii.gz" "${tmp}/${tag}_reorient.nii.gz"
              else
                  Do_cmd 3dresample -orient LPI -prefix "${tmp}/${tag}_reorient.nii.gz" -inset "$rawNifti"
                  Do_cmd fslreorient2std "${tmp}/${tag}_reorient.nii.gz" "${tmp}/${tag}_reorient.nii.gz"
              fi

              # Remove slices to make an even number of slices in all directions (requisite for topup).
              #dimensions=`fslhd ${tmp}/${tag}_reorient.nii.gz | grep -E "^dim[1-3]" | awk '{print $NF}'`
              #newDimensions=`echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}'`
              #Do_cmd fslroi ${tmp}/${tag}_reorient.nii.gz ${tmp}/${tag}_sliceCut.nii.gz `echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /'` # I removed -1 -1

              # Skipping fslroi step. Rename files for simplicity
              mv "${tmp}/${tag}_reorient.nii.gz" "${tmp}/${tag}_sliceCut.nii.gz"

              # Motion correction within scans
              Do_cmd fslmaths "${tmp}/${tag}_sliceCut.nii.gz" -Tmean "${tmp}/${tag}_sliceCutMean.nii.gz"
              Do_cmd 3dvolreg -Fourier -twopass -base "${tmp}/${tag}_sliceCutMean.nii.gz" \
                              -zpad 4 -prefix "${tmp}/${tag}_mc.nii.gz" \
                              -1Dfile "${taskfmri_volum}/${idBIDS}_space-taskfmri_${tag}.1D" \
                              "${tmp}/${tag}_sliceCut.nii.gz"
              Do_cmd fslmaths "${tmp}/${tag}_mc.nii.gz" -Tmean "${tmp}/${tag}_mcMean.nii.gz"
        fi
    done

    # Calculate motion outliers with FSL
    if [[ ! -f "${taskfmri_volum}/${idBIDS}_space-taskfmri_singleecho.1D" ]]; then
        Do_cmd fsl_motion_outliers -i "${tmp}/mainScan_sliceCut.nii.gz" \
                                   -o "${taskfmri_volum}/${idBIDS}_space-taskfmri_spikeRegressors_FD.1D" \
                                   -s "${taskfmri_volum}/${idBIDS}_space-taskfmri_metric_FD.1D" --fd
        Do_cmd mv "${taskfmri_volum}/${idBIDS}_space-taskfmri_mainScan.1D ${taskfmri_volum}/${idBIDS}_space-taskfmri_singleecho.1D"; ((Nsteps++))
    else
        Info "Subject ${id} has a singleecho.1D with motion outliers"; ((Nsteps++))
    fi

    # If ONLY Reverse phase scan is provided mainScan will be the mainPhaseScan
    if [[ -f "$reversePhaseScan" ]] && [[ ! -f "$mainPhaseScan" ]]; then
          main_pe="NO-main_pe"
          Warning "reversePhaseScan was found but NO mainPhaseScan, using mainScan as mainPhaseScan"
          mainPhaseScan="${tmp}/mainPhaseScan_mc.nii.gz"
          Do_cmd cp "${tmp}/mainScan_mc.nii.gz" "$mainPhaseScan"
          Do_cmd cp "${tmp}/mainScan_mcMean.nii.gz" "${tmp}/singleecho_mainPhaseAlignedMean.nii.gz"
    fi

    # Only do distortion correction if field maps were provided, if not then rename the scan to distortionCorrected (just to make the next lines of code easy).
    if [ -z "${mainPhaseScan}" ] || [ -z "${reversePhaseScan}" ]; then
        Warning "No AP or PA acquisition was found, TOPUP will be skip!!!!!!!"
        export statusTopUp="NO"
        Do_cmd mv -v "${tmp}/mainScan_mc.nii.gz" "${singleecho}"; ((Nsteps++))
    else
        if [[ ! -f "${taskfmri_volum}/TOPUP.txt" ]] && [[ ! -f "${singleecho}" ]]; then
            mainPhaseScanMean=$(find "$tmp"    -maxdepth 1 -name "*mainPhaseScan*_mcMean.nii.gz")
            mainPhaseScan=$(find "$tmp"        -maxdepth 1 -name "*mainPhaseScan*_mc.nii.gz")
            reversePhaseScanMean=$(find "$tmp" -maxdepth 1 -name "*reversePhaseScan*_mcMean.nii.gz")
            reversePhaseScan=$(find "$tmp"     -maxdepth 1 -name "*reversePhaseScan*_mc.nii.gz")
            mainScan=$(find "$tmp"             -maxdepth 1 -name "*mainScan*_mc.nii.gz")

            Do_cmd flirt -in "$reversePhaseScanMean" -ref "$tmp"/mainScan_mcMean.nii.gz -omat "$tmp"/singleecho_tmpXfmSecondary.omat
            Do_cmd flirt -in "$reversePhaseScan" -ref "$tmp"/mainScan_mcMean.nii.gz -applyxfm -init "$tmp"/singleecho_tmpXfmSecondary.omat -out "$tmp"/singleecho_secondaryPhaseAligned.nii.gz
            Do_cmd fslmaths "$tmp"/singleecho_secondaryPhaseAligned.nii.gz -Tmean "$tmp"/singleecho_secondaryPhaseAlignedMean.nii.gz

            if [[ "$main_pe" != "NO-main_pe" ]]; then
                Do_cmd flirt -in "$mainPhaseScanMean" -ref "$tmp"/mainScan_mcMean.nii.gz -omat "$tmp"/singleecho_tmpXfmMain.omat
                Do_cmd flirt -in "$mainPhaseScan" -ref "$tmp"/mainScan_mcMean.nii.gz -applyxfm -init "$tmp"/singleecho_tmpXfmMain.omat -out "$tmp"/singleecho_mainPhaseAligned.nii.gz
                Do_cmd fslmaths "$tmp"/singleecho_mainPhaseAligned.nii.gz -Tmean "$tmp"/singleecho_mainPhaseAlignedMean.nii.gz
            fi

            # Distortion correction
            echo -e "0 1 0 ${readoutTime} \n0 -1 0 ${readoutTime}" > "$tmp"/singleecho_topupDataIn.txt
            Info "topup datain:\n$(cat "${tmp}"/singleecho_topupDataIn.txt)"
            Do_cmd fslmerge -t "${tmp}/singleecho_mergeForTopUp.nii.gz" "${tmp}/singleecho_mainPhaseAlignedMean.nii.gz" "${tmp}/singleecho_secondaryPhaseAlignedMean.nii.gz"
            Do_cmd topup --imain="${tmp}/singleecho_mergeForTopUp.nii.gz" --datain="${tmp}/singleecho_topupDataIn.txt" --config="${topupConfigFile}" --out="$tmp/singleecho_topup"
            Do_cmd applytopup --imain="${mainScan}" --inindex=1 --datain="${tmp}/singleecho_topupDataIn.txt" --topup="${tmp}/singleecho_topup" --method=jac --out="${singleecho}"

            # # taskfmri for registration
            # Do_cmd applytopup --imain="${tmp}/singleecho_mainPhaseAligned.nii.gz" --inindex=1 --datain="${tmp}/singleecho_topupDataIn.txt" --topup="${tmp}/singleecho_topup" --method=jac --out="${tmp}/singleecho_mainPhaseAlignedTopup.nii.gz"
            # Do_cmd fslmaths "${tmp}/singleecho_mainPhaseAlignedTopup.nii.gz" -Tmean "$taskfmri4reg"

            # Check if it worked
            if [[ ! -f "${singleecho}" ]]; then Error "Something went wrong with TOPUP check ${tmp} and log:\n\t\t${dir_logs}/proc_taskfmri.txt"; exit; fi
            export statusTopUp="YES"; ((Nsteps++))
        else
            Info "Subject ${id} has singleecho in fmrispace with TOPUP"; export statusTopUp="YES"; ((Nsteps++))
        fi
    fi
else
      Info "Subject ${id} has a singleecho_fmrispace processed"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
Info "!!!!!  goin str8 to ICA-FIX yo  !!!!!"

fmri_mean="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_mean.nii.gz"
fmri_HP="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_HP.nii.gz"
fmri_brain="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_brain.nii.gz"
fmri_mask="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_brain_mask.nii.gz"

if [[ ! -f "$fmri_mask" ]] || [[ ! -f "$fmri_brain" ]]; then
    Info "Generating a taskfmri binary mask"
    # Calculates the mean taskfmri volume
    Do_cmd fslmaths "$singleecho" -Tmean "$fmri_mean"

    # Creates a mask from the motion corrected time series
    Do_cmd bet "$fmri_mean" "${fmri_brain}" -m -n

    # masked mean taskfmri time series
    Do_cmd fslmaths "$fmri_mean" -mul "$fmri_mask" "$fmri_brain"
    if [[ -f "${fmri_mask}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a binary mask of the taskfmri"; ((Nsteps++))
fi

# High-pass filter - Remove all frequencies EXCEPT those in the range
if [[ ! -f "$fmri_HP" ]]; then
    Info "High pass filter"
    Do_cmd 3dTproject -input "${singleecho}" -prefix "$fmri_HP" -passband 0.01 666
        if [[ -f "${fmri_HP}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has High-pass filter"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# run MELODIC for ICA-FIX
melodic_IC="${taskfmri_ICA}/filtered_func_data.ica/melodic_IC.nii.gz"
fmri_filtered="${taskfmri_ICA}/filtered_func_data.nii.gz"

# melodic will run ONLY no FIX option is selected
if [[ "$noFIX" -eq 0 ]] && [[ ! -f "${melodic_IC}" ]]; then
    [[ ! -d "${taskfmri_ICA}" ]] && Do_cmd mkdir -p "${taskfmri_ICA}"
    Info "Running melodic"
    Do_cmd cp "$fmri_HP" "$fmri_filtered"
    Do_cmd melodic --in="${fmri_filtered}" \
          --tr="${RepetitionTime}" \
          --nobet \
          --mask="${fmri_mask}" \
          --bgthreshold=3 \
          --mmthresh=0.5 \
          --report \
          --Oall \
          --outdir="${taskfmri_ICA}/filtered_func_data.ica" \
          --Omean="${taskfmri_ICA}/mean_func.nii.gz"
    if [[ -f "${melodic_IC}" ]]; then export statusMel="YES"; else export statusMel="FAILED"; fi
else
    Info "Subject ${id} has MELODIC outputs"; export statusMel="YES"
fi
if [[ "$noFIX" -eq 1 ]]; then export statusMel="NO"; fi
#------------------------------------------------------------------------------#
fmri_in_T1nativepro="${proc_struct}/${idBIDS}_space-nativepro_desc-${tagMRI}_bold.nii.gz"
T1nativepro_in_fmri="${taskfmri_volum}/${idBIDS}_space-${tagMRI}_t1w.nii.gz"
str_taskfmri_affine="${dir_warp}/${idBIDS}_taskfmri_from-${tagMRI}_to-nativepro_mode-image_desc-affine_"
mat_taskfmri_affine="${str_taskfmri_affine}0GenericAffine.mat"
t1bold="${proc_struct}/${idBIDS}_space-nativepro_desc-t1wbold.nii.gz"

str_taskfmri_SyN="${dir_warp}/${idBIDS}_taskfmri_from-nativepro_${tagMRI}_to-${tagMRI}_mode-image_desc-SyN_"
SyN_taskfmri_affine="${str_taskfmri_SyN}0GenericAffine.mat"
SyN_taskfmri_warp="${str_taskfmri_SyN}1Warp.nii.gz"
SyN_taskfmri_Invwarp="${str_taskfmri_SyN}1InverseWarp.nii.gz"

# Registration to native pro
Nreg=$(ls "$mat_taskfmri_affine" "$fmri_in_T1nativepro" "$T1nativepro_in_fmri" 2>/dev/null | wc -l )
if [[ "$Nreg" -lt 3 ]]; then
    # if [[ -f "$taskfmri4reg" ]]; then
    #     Do_cmd fslmaths "$taskfmri4reg" -mul "$fmri_mask" "$fmri_brain"
    # fi
    if [[! -f "${t1bold}" ]]; then
        Info "Creating a synthetic BOLD image for registration"
        # Inverse T1w
        Do_cmd ImageMath 3 "${tmp}/${id}_t1w_nativepro_NEG.nii.gz" Neg "$T1nativepro"
        # Dilate the T1-mask
        #Do_cmd ImageMath 3 "${tmp}/${id}_t1w_mask_dil-2.nii.gz" MD "$T1nativepro_mask" 2
        # Masked the inverted T1w
        Do_cmd ImageMath 3 "${tmp}/${id}_t1w_nativepro_NEG_brain.nii.gz" m "${tmp}/${id}_t1w_nativepro_NEG.nii.gz" "$T1nativepro_mask"
        # Match histograms values acording to taskfmri
        Do_cmd ImageMath 3 "${tmp}/${id}_t1w_nativepro_NEG-rescaled.nii.gz" HistogramMatch "${tmp}/${id}_t1w_nativepro_NEG_brain.nii.gz" "$fmri_brain"
        # Smoothing
        Do_cmd ImageMath 3 "$t1bold" G "${tmp}/${id}_t1w_nativepro_NEG-rescaled.nii.gz" 0.35
    else
        Info "Subject ${id} has a synthetic BOLD image for registration"
    fi

    Info "Registering fmri space to nativepro"

    # Affine from taskfmri to t1-nativepro
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$T1nativepro_brain" -m "$fmri_brain" -o "$str_taskfmri_affine" -t a -n "$threads" -p d
    Do_cmd antsApplyTransforms -d 3 -i "$t1bold" -r "$fmri_brain" -t ["$mat_taskfmri_affine",1] -o "${tmp}/T1bold_in_fmri.nii.gz" -v -u int

    if [[ ${regAffine}  == "FALSE" ]]; then
        # SyN from T1_nativepro to t1-nativepro
        Do_cmd antsRegistrationSyN.sh -d 3 -m "${tmp}/T1bold_in_fmri.nii.gz" -f "$fmri_brain" -o "$str_taskfmri_SyN" -t s -n "$threads" -p d #-i "$mat_taskfmri_affine"
        export reg="Affine+SyN"
        transformsInv="-t ${SyN_taskfmri_warp} -t ${SyN_taskfmri_affine} -t [${mat_taskfmri_affine},1]" # T1nativepro to taskfmri
        transform="-t ${mat_taskfmri_affine} -t [${SyN_taskfmri_affine},1] -t ${SyN_taskfmri_Invwarp}"  # taskfmri to T1nativepro
        xfmat="-t ${SyN_taskfmri_affine} -t [${mat_taskfmri_affine},1]" # T1nativepro to taskfmri only lineal
    elif [[ ${regAffine}  == "TRUE" ]]; then
        export reg="Affine"
        transformsInv="-t [${mat_taskfmri_affine},1]"
        transform="-t ${SyN_taskfmri_affine}"
        xfmat="-t [${mat_taskfmri_affine},1]"
    fi

    # fmri to t1-nativepro
    Do_cmd antsApplyTransforms -d 3 -i "$fmri_brain" -r "$t1bold" "${transform}" -o "$fmri_in_T1nativepro" -v -u int
    # t1-nativepro to fmri
    Do_cmd antsApplyTransforms -d 3 -i "$T1nativepro_brain" -r "$fmri_brain" "${transformsInv}" -o "${T1nativepro_in_fmri}" -v -u int

    if [[ -d "${taskfmri_ICA}/filtered_func_data.ica" ]]; then Do_cmd cp "${T1nativepro_in_fmri}" "${taskfmri_ICA}/filtered_func_data.ica/t1w2fmri_brain.nii.gz"; fi
    if [[ -f "${SyN_taskfmri_Invwarp}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a taskfmri volume and transformation matrix in T1nativepro space"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Register taskfmri to Freesurfer space with Freesurfer
fmri2fs_dat="${dir_warp}/${idBIDS}_from-${tagMRI}_to-fsnative_bbr.dat"
if [[ ! -f "${fmri2fs_dat}" ]] ; then
  Info "Registering fmri to FreeSurfer space"
    Do_cmd bbregister --s "$BIDSanat" --mov "$fmri_mean" --reg "${fmri2fs_dat}" --o "${dir_warp}/${idBIDS}_from-${tagMRI}_to-fsnative_bbr_outbbreg_FIX.nii.gz" --init-fsl --bold
    if [[ -f "${fmri2fs_dat}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a dat transformation matrix from fmri to Freesurfer space"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# run ICA-FIX IF melodic ran succesfully
fix_output="${taskfmri_ICA}/filtered_func_data_clean.nii.gz"
fmri_processed="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_clean.nii.gz"

# Run if fmri_clean does not exist
if [[ "$noFIX" -eq 0 ]]; then
    if [[ ! -f "${fmri_processed}" ]] ; then
          if  [[ -f "${melodic_IC}" ]] && [[ -f $(which fix) ]]; then
              if [[ ! -f "${fix_output}" ]] ; then
                    Info "Getting ICA-FIX requirements"
                    Do_cmd mkdir -p "${taskfmri_ICA}"/{reg,mc}
                    # FIX requirements - https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide
                    # $fmri_filtered                                                                                 preprocessed 4D data
                    # $melodic_IC                                                                                    melodic (command-line program) full output directory
                    Do_cmd cp "${taskfmri_volum}/${idBIDS}_space-taskfmri_singleecho.1D" "${taskfmri_ICA}/mc/prefiltered_func_data_mcf.par"   # motion parameters created by mcflirt
                    Do_cmd cp "$fmri_mask" "${taskfmri_ICA}/mask.nii.gz"                                                 # valid mask relating to the 4D data
                    Do_cmd cp "${taskfmri_ICA}/filtered_func_data.ica/mean.nii.gz" "${taskfmri_ICA}/mean_func.nii.gz"      # temporal mean of 4D data
                    middleSlice=$(mrinfo "$fmri_filtered" -size | awk -F ' ' '{printf "%.0f\n", $4/2}')
                    Do_cmd fslroi "$fmri_filtered" "${taskfmri_ICA}/reg/example_func.nii.gz" "$middleSlice" 1          # example middle image from 4D data
                    Do_cmd cp "$T1nativepro_brain" "${taskfmri_ICA}/reg/highres.nii.gz"                                  # brain-extracted structural

                    # REQUIRED by FIX - reg/highres2example_func.mat                                               # FLIRT transform from structural to functional space
                    if [[ ! -f "${taskfmri_ICA}/reg/highres2example_func.mat" ]]; then
                        # Get transformation matrix T1native to taskfmri space (ICA-FIX requirement)
                        Do_cmd antsApplyTransforms -v 1 -o Linear["$tmp/highres2example_func.mat",0] "${xfmat}"
                        # Transform matrix: ANTs (itk binary) to text
                        Do_cmd ConvertTransformFile 3 "$tmp/highres2example_func.mat" "$tmp/highres2example_func.txt"

                        # Fixing the transformations incompatibility between ANTS and FSL
                        tmp_ants2fsl_mat="$tmp/itk2fsl_highres2example_func.mat"
                        # Transform matrix: ITK text to matrix (FSL format)
                        Do_cmd lta_convert --initk "$tmp/highres2example_func.txt" --outfsl "$tmp_ants2fsl_mat" --src "$T1nativepro" --trg "$fmri_brain"
                        # apply transformation with FSL
                        Do_cmd flirt -in "$T1nativepro" -out "$tmp/t1w2fmri_brain_ants2fsl.nii.gz" -ref "$fmri_brain" -applyxfm -init "$tmp_ants2fsl_mat"
                        # correct transformation matrix
                        Do_cmd flirt -in "$tmp/t1w2fmri_brain_ants2fsl.nii.gz" -ref "$T1nativepro_in_fmri" -omat "$tmp/ants2fsl_fixed.omat" -cost mutualinfo -searchcost mutualinfo -dof 6
                        # concatenate the matrices to fix the transformation matrix
                        Do_cmd convert_xfm -concat "$tmp"/ants2fsl_fixed.omat -omat "${taskfmri_ICA}/reg/highres2example_func.mat" "$tmp_ants2fsl_mat"
                    else Info "Subject ${id} has reg/highres2example_func.mat for ICA-FIX"; fi

                    Info "Running ICA-FIX"
                    Do_cmd fix "$taskfmri_ICA" "$icafixTraining" 20 -m -h 100

                    # Replace file if melodic ran correctly - Change single-echo files for clean ones
                    if [[ -f "$fix_output" ]]; then
                        yes | Do_cmd cp -rf "$fix_output" "$fmri_processed"
                        export statusFIX="YES"
                    else
                        Error "FIX failed, but MELODIC ran log file:\n\t $(ls "${dir_logs}"/proc_taskfmri_*.txt)"; exit
                    fi
              else
                    Info "Subject ${id} has filtered_func_data_clean from ICA-FIX already"
                    cp -rf "$fix_output" "$fmri_processed"; export statusFIX="YES"
              fi
          else
              Warning "!!!!  Melodic Failed and/or FIX was not found, check the software installation !!!!
                             If you've installed FIX try to install required R packages and re-run:
                             'kernlab','ROCR','class','party','e1071','randomForest'"
              Do_cmd cp -rf "$fmri_HP" "$fmri_processed"
              export statusFIX="NO"
          fi
    else
        Info "Subject ${id} has singleecho_fmrispace_clean"; export statusFIX="YES"
    fi
    json_taskfmri "${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_clean.json"
else
    # Skip FIX processing but rename variables anyways for simplicity
    Info "Further processing will be performed on distorsion corrected images."
    cp -rf "${fmri_HP}" "$fmri_processed"
    if [[ "$noFIX" -eq 1 ]]; then export statusFIX="NO"; fi
    json_taskfmri "${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_clean.json"
fi


#------------------------------------------------------------------------------#
global_signal="${taskfmri_volum}/${idBIDS}_space-taskfmri_global.txt"
if [[ ! -f "${global_signal}" ]] ; then
    Info "Calculating tissue-specific and global signals changes"
    tissues=(CSF GM WM)
    for idx in {0..2}; do
        tissue=${tissues[$idx]}
        tissuemap="${dir_anat}/${BIDSanat}_space-nativepro_t1w_brain_pve_${idx}.nii.gz"
        tissue_series="${taskfmri_volum}/${idBIDS}_space-taskfmri_pve_${tissue}.txt"
        if [[ ! -f "${tissue_series}" ]] ; then
            Do_cmd antsApplyTransforms -d 3 -i "$tissuemap" -r "$fmri_mean" "${transformsInv}" -o "${tmp}/${idBIDS}_space-taskfmri_${tissue}.nii.gz" -v -u int
            Do_cmd fslmaths "${tmp}/${idBIDS}_space-taskfmri_${tissue}.nii.gz" -thr 0.9 "${tmp}/${idBIDS}_space-taskfmri_${tissue}.nii.gz"
            Do_cmd fslmeants -i "$fmri_processed" -o "$tissue_series" -m "${tmp}/${idBIDS}_space-taskfmri_${tissue}.nii.gz" -w
        else
             Info "Subject ${idBIDS} has $tissue time-series"
        fi
    done
    # Global signal from brain mask
    Do_cmd fslmeants -i "$fmri_processed" -o "$global_signal" -m "${fmri_mask}" -w
    if [[ -f "${global_signal}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has Global time-series"; ((Nsteps++))
fi

# Motion confound
spikeRegressors="${taskfmri_volum}/${idBIDS}_space-taskfmri_spikeRegressors_REFRMS.1D"
if [[ ! -f "$spikeRegressors" ]] ; then
    Do_cmd fsl_motion_outliers -i "$fmri_processed" -o "$spikeRegressors" -s "${taskfmri_volum}/${idBIDS}_space-taskfmri_metric_REFRMS.1D" --refmse --nomoco
    if [[ -f "$spikeRegressors" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a spike Regressors from fsl_motion_outliers"; ((Nsteps++))
fi

# Register to surface
# If three surfaces are found skipp this step
Nsurf=$(ls "${taskfmri_surf}/${idBIDS}"_taskfmri_space-fsnative_?h.mgh \
            "${taskfmri_surf}/${idBIDS}"_taskfmri_space-fsnative_?h_10mm.mgh \
            "${taskfmri_surf}/${idBIDS}"_taskfmri_space-fsaverage5_?h_10mm.mgh \
            "${taskfmri_surf}/${idBIDS}"_taskfmri_space-conte69-32k_?h_10mm.mgh 2>/dev/null | wc -l)

if [ "$Nsurf" -lt 8 ]; then
for hemisphere in lh rh; do
    HEMI=$(echo "${hemisphere/h/}" | tr [:lower:] [:upper:])
    Info "Mapping volumetric timeseries to native surface ${hemisphere}"
    vol2surfTS="${taskfmri_surf}/${idBIDS}"_taskfmri_space-fsnative_${hemisphere}.mgh
    if [[ ! -f "$vol2surfTS" ]] ; then

          # Map the non high-passed volumetric timeseries to the surface so we can compute tSNR
          Do_cmd mri_vol2surf \
              --mov "$singleecho" \
              --reg "$fmri2fs_dat" \
              --projfrac-avg 0.2 0.8 0.1 \
              --trgsubject "$BIDSanat" \
              --interp trilinear \
              --hemi "${hemisphere}" \
              --out "${taskfmri_surf}/${idBIDS}"_taskfmri_space-fsnative_"${hemisphere}"_NoHP.mgh

          # Map processed timeseries to surface
          Do_cmd mri_vol2surf \
              --mov "$fmri_processed "\
              --reg "$fmri2fs_dat" \
              --projfrac-avg 0.2 0.8 0.1 \
              --trgsubject "$BIDSanat" \
              --interp trilinear \
              --hemi "${hemisphere}" \
              --out "$vol2surfTS"

          if [[ -f "$vol2surfTS" ]] ; then ((Nsteps++)); fi
    else
        Info "Subject ${id} volumetric timeseries have been mapped to ${HEMI} cortical surface"; ((Nsteps++))
    fi

    # Convert native timeseries to gifti
    Do_cmd mri_convert "${taskfmri_surf}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}.mgh" "${tmp}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}.func.gii"

    # Apply smoothing on native surface
    out_surf_native="${taskfmri_surf}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}_10mm.mgh"
    if [[ ! -f "$out_surf_native" ]] ; then
          if [[ "$smooth" == 1 ]] ; then
            Do_cmd wb_command -metric-smoothing \
                "${dir_freesurfer}/surf/${hemisphere}.midthickness.surf.gii"  \
                "${tmp}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}.func.gii" \
                10 \
                "${tmp}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}_10mm.func.gii"
            Do_cmd mri_convert "${tmp}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}_10mm.func.gii" "$out_surf_native"
          else
            Do_cmd mri_surf2surf \
                --hemi "${hemisphere}" \
                --srcsubject "$BIDSanat" \
                --sval "${taskfmri_surf}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}.mgh" \
                --trgsubject "$BIDSanat" \
                --tval "$out_surf_native" \
                --fwhm-trg 10
          fi
    if [[ -f "$out_surf_native" ]] ; then ((Nsteps++)); fi
    else
        Info "Subject ${id} has native timeseries smoothed on ${HEMI} surface"; ((Nsteps++))
    fi

    # Register to fsa5 and smooth
    out_surf_fsa5="${taskfmri_surf}/${idBIDS}_taskfmri_space-fsaverage5_${hemisphere}.mgh"
    if [[ ! -f "$out_surf_fsa5" ]] ; then
         Do_cmd mri_surf2surf \
            --hemi "${hemisphere}" \
            --srcsubject "$BIDSanat" \
            --sval "${taskfmri_surf}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}.mgh" \
            --trgsubject fsaverage5 \
            --tval "$out_surf_fsa5"
         if [[ -f "$out_surf_fsa5" ]] ; then ((Nsteps++)); fi
    else
         Info "Subject ${id} has timeseries mapped to ${HEMI} fsa5"; ((Nsteps++))
    fi

    out_surf_fsa5_sm="${taskfmri_surf}/${idBIDS}_taskfmri_space-fsaverage5_${hemisphere}_10mm.mgh"
    if [[ ! -f "$out_surf_fsa5_sm" ]] ; then
         Do_cmd mri_surf2surf \
            --hemi "${hemisphere}" \
            --srcsubject "$BIDSanat" \
            --sval "${taskfmri_surf}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}.mgh" \
            --trgsubject fsaverage5 \
            --tval "$out_surf_fsa5_sm" \
            --fwhm-trg 10
         if [[ -f "$out_surf_fsa5_sm" ]] ; then ((Nsteps++)); fi
    else
         Info "Subject ${id} has smoothed timeseries mapped to ${HEMI} fsa5"; ((Nsteps++))
    fi

    # Register to conte69 and smooth
    out_surf="${taskfmri_surf}/${idBIDS}_taskfmri_space-conte69-32k_${hemisphere}_10mm.mgh"
    if [[ ! -f "$out_surf" ]] ; then
          # Register to conte69
          Do_cmd wb_command -metric-resample \
              "${tmp}/${idBIDS}_taskfmri_space-fsnative_${hemisphere}.func.gii" \
              "${dir_conte69}/${BIDSanat}_${hemisphere}_sphereReg.surf.gii" \
              "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMI}.sphere.32k_fs_LR.surf.gii" \
              ADAP_BARY_AREA \
              "${tmp}/${idBIDS}_taskfmri_space-conte69-32k_${hemisphere}.func.gii" \
              -area-surfs \
              "${dir_freesurfer}/surf/${hemisphere}.midthickness.surf.gii" \
              "${dir_conte69}/${BIDSanat}_space-conte69-32k_desc-${hemisphere}_midthickness.surf.gii"
          # Apply smooth on conte69
          Do_cmd wb_command -metric-smoothing \
              "${util_surface}/fsaverage.${HEMI}.midthickness_orig.32k_fs_LR.surf.gii" \
              "${tmp}/${idBIDS}_taskfmri_space-conte69-32k_${hemisphere}.func.gii" \
              10 \
              "${tmp}/${idBIDS}_taskfmri_space-conte69-32k_${hemisphere}_10mm.func.gii"

          Do_cmd mri_convert "${tmp}/${idBIDS}_taskfmri_space-conte69-32k_${hemisphere}.func.gii" "${taskfmri_surf}/${idBIDS}_taskfmri_space-conte69-32k_${hemisphere}.mgh"
          Do_cmd mri_convert "${tmp}/${idBIDS}_taskfmri_space-conte69-32k_${hemisphere}_10mm.func.gii" "${out_surf}"
          if [[ -f "$out_surf" ]] ; then ((Nsteps++)); fi
    else
          Info "Subject ${id} has a singleecho fmri2fs ${hemisphere} on conte69-32k 10mm surface"; ((Nsteps++))
    fi
done
else
  Info "Subject ${id} has a singleecho fmri2fs on all surfaces: native, native_fwhm-10mm, fsa5, fsa5_fwhm-10mm, and conte69fwhm_10mm"; Nsteps=$((Nsteps+10))
fi

#------------------------------------------------------------------------------#
#                           S U B C O R T E X
# Subcortical segmentation (nativepro) to taskfmri space
taskfmri_subcortex="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_subcortical.nii.gz"
timese_subcortex="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_timeseries_subcortical.txt"

if [[ ! -f "$timese_subcortex" ]] ; then
      Info "Getting subcortical timeseries"
      Do_cmd antsApplyTransforms -d 3 -i "$T1_seg_subcortex" -r "$fmri_mean" -n GenericLabel  "${transformsInv}" -o "$taskfmri_subcortex" -v -u int
      # Extract subcortical timeseries
      # Output: ascii textis this correct? file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
      Do_cmd mri_segstats --i "$fmri_processed" --seg "$taskfmri_subcortex" --exclude 0 --exclude 16 --avgwf "$timese_subcortex"
      if [[ -f "$timese_subcortex" ]] ; then ((Nsteps++)); fi
else
      Info "Subject ${id} has taskfmri subcortical time-series"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
#                           C E R E B E L L U M
taskfmri_cerebellum="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_cerebellum.nii.gz"
timese_cerebellum="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_timeseries_cerebellum.txt"
stats_cerebellum="${taskfmri_volum}/${idBIDS}_space-taskfmri_desc-singleecho_cerebellum_roi_stats.txt"

if [[ ! -f "$timese_cerebellum" ]] ; then
      Info "Getting cerebellar timeseries"
      Do_cmd antsApplyTransforms -d 3 -i "$T1_seg_cerebellum" -r "$fmri_mean" -n GenericLabel "${transformsInv}" -o "$taskfmri_cerebellum" -v -u int
      # Extract cerebellar timeseries (mean, one ts per segemented structure, exluding nuclei because many are too small for our resolution)
      # Output: ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
      Do_cmd mri_segstats --i "$fmri_processed" --seg "$taskfmri_cerebellum" --exclude 0 --avgwf "$timese_cerebellum"
      3dROIstats -mask "$taskfmri_cerebellum" -nzsum "$taskfmri_cerebellum" > "$stats_cerebellum"
      if [[ -f "$timese_cerebellum" ]] ; then ((Nsteps++)); fi
else
      Info "Subject ${id} has taskfmri cerebellar time-series"; ((Nsteps++))
fi

# -----------------------------------------------------------------------------------------------
# QC: taskfmri processing Input files
if [[ ${fmri_acq} == "FALSE" ]]; then QC_proc-taskfmri; fi

#------------------------------------------------------------------------------#
# run post-taskfmri
cleanTS="${taskfmri_surf}/${idBIDS}_taskfmri_space-conte69-32k_desc-timeseries_clean.txt"
if [[ ! -f "$cleanTS" ]] ; then
    Info "Running Task fMRI post processing"
    labelDirectory="${dir_freesurfer}/label/"
    Do_cmd python "$MICAPIPE"/functions/03_FC.py "$idBIDS" "${proc_taskfmri}/$taskScanStr" "$labelDirectory" "$util_parcelations" "$dir_volum" "$performNSR" "$performGSR" "taskfmri"
    if [[ -f "$cleanTS" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has post-processed conte69 time-series"; ((Nsteps++))
fi


#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 21 ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
Title "Task fMRI processing and post processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m:
\tSteps completed : $(printf "%02d" "$Nsteps")/21
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/proc_taskfmri_*.txt)"
if [[ ${fmri_acq} == "FALSE" ]]; then
    grep -v "${id}, ${SES/ses-/}, proc_taskfmri" "${out}/micapipe_processed_sub.csv" > "${tmp}/tmpfile" && mv "${tmp}/tmpfile" "${out}/micapipe_processed_sub.csv"
    #echo "${id}, ${SES/ses-/}, proc_taskfmri, ${status}, $(printf "%02d" "$Nsteps")/21, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipe_processed_sub.csv"
fi
cleanup "$tmp" "$nocleanup" "$here"
