#!/bin/bash

SITE=http://www.bruzual.org/bc03/Original_version_2003
IMF=(chabrier salpeter)
IMF_SHORT=(chab salp)
IMF_REP=(ch sa)
RESOLUTIONS=(hr lr)
METAL=(m22 m32 m42 m52 m62 m72)
METAL_REP=(z0001 z0004 z004 z008 z02 z05)

for RES in ${RESOLUTIONS[@]}; do
    mkdir -p ssp.${RES}
done

mkdir -p tmp_download_bc03
cd tmp_download_bc03
for ((i=0; i<${#IMF[@]}; ++i)); do
    wget ${SITE}/bc03.models.padova_1994_${IMF[i]}_imf.tar.gz -O archive.tar.gz
    tar -xvzf archive.tar.gz

    for RES in ${RESOLUTIONS[@]}; do
        for ((m=0; m<${#METAL[@]}; ++m)); do
            OFILE=bc03/models/Padova1994/${IMF[i]}/bc2003_${RES}_${METAL[m]}_${IMF_SHORT[i]}_ssp.ised_ASCII
            gzip -df ${OFILE}.gz
            mv ${OFILE} ../ssp.${RES}/bc03_${RES}_${IMF_REP[i]}_${METAL_REP[m]}.ised_ASCII
        done
    done

    rm -r bc03
done

cd ..
rm -r tmp_download_bc03
