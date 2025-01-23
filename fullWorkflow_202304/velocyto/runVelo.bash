mkdir -p velocyto
module load R
module load python
rdsFile=$1

Rscript runPerSample.R $1

filename="${rdsFile##*/}"
filename="${filename%.*}"

python runVelo.py $filename.h5ad
