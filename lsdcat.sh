#! /bin/bash


export PYTHONPATH=${PYTHONPATH}:/Users/mxhf/work/MPE/hetdex/src/lsdcat/lib
export PATH=${PATH}:/Users/mxhf/work/MPE/hetdex/src/lsdcat
export PATH=${PATH}:/Users/mxhf/work/MPE/hetdex/src/lsdcat/tools

INCUBE=$1
EXPMAP=$2

# PSF Parameters
field_id=01
p0=1.
p1=0.

# PATH + filenames - note how the filenames are generated here from field_id
input_cube=$INCUBE
exp_map=$EXPMAP

input_cube_base=`basename ${input_cube}`

echo Field ID: $field_id
echo p0: $p0
echo p1: $p1


floor_output=${output_dir}c${input_cube_base}
floor_com="python ./floor.py ${input_cube} -200 200 ${floor_output}"
echo ${floor_com}
${floor_com}


# median filtering
floor_output_base=`basename ${floor_output} .fits`
med_filt_output=${output_dir}median_filtered_${floor_output_base}
med_filt_com="median-filter-cube.py ${floor_output_base} --signalHDU=0 --varHDU=0 --num_cpu=4 --width=151 \
     --output=${med_filt_output}"
echo ${med_filt_com}
${med_filt_com}

# applying effective noise
#apply_eff_noise_com="apply_eff_noise.py ${med_filt_output} ${effnoise_file} --NHDU=1 --blowup --rsexp \
#    --output=${med_filt_output}_effnoised.fits"
#echo ${apply_eff_noise_com}
#${apply_eff_noise_com}

#spatial cross-correlation
med_filt_output_base=`basename ${med_filt_output} .fits`
spat_cced_out=${output_dir}spat_cced_${med_filt_output_base}
lsd_cc_spat_com="lsd_cc_spatial.py --input=${med_filt_output} --SHDU=0 \
    --threads=4 --gaussian --lambda0=4500 -p0=${p0} -p1=${p1} --output=${spat_cced_out}"
echo ${lsd_cc_spat_com}
${lsd_cc_spat_com}

# spectral cross-correlation
spat_cced_out_base=`basename ${spat_cced_out}`
spec_cced_out=${output_dir}spec_cced_${spat_cced_out_base}
lsd_cc_spec_com="lsd_cc_spectral.py --input=${spat_cced_out} --threads=4 --FWHM=330 --SHDU=0 --NHDU=1 \
    --output=${spec_cced_out}"
echo ${lsd_cc_spec_com}
${lsd_cc_spec_com}

# creation of S/N cube
spec_cced_out_base=`basename ${spec_cced_out}`
sn_out=${output_dir}sn_${spec_cced_out_base}
sn_com="s2n-cube.py --input=${spec_cced_out_base}  \
    --output=${sn_out} --clobber"
echo ${sn_com}
${sn_com}



sn_out_base=`basename ${sn_out}`
#cat_com="lsd_cat.py -t 1.5 -i ${sn_out_base} --clobber --tabvalues ID,X_PEAK_SN,Y_PEAK_SN,Z_PEAK_SN,RA_PEAK_SN,DEC_PEAK_SN,LAMBDA_PEAK_SN,DETSN_MAX -e ${exp_map}"
cat_com="lsd_cat.py -t 1.5 -i ${sn_out_base} --clobber --tabvalues ID,X_PEAK_SN,Y_PEAK_SN,Z_PEAK_SN,RA_PEAK_SN,DEC_PEAK_SN,LAMBDA_PEAK_SN,DETSN_MAX "
echo ${cat_com}
${cat_com}

