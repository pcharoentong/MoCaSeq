nextflow run roland-rad-lab/MoCaSeq -r human-pipeline-nextflow-2  -profile charliecloud,slurm -work-dir /gpfs/scratch/pn29ya/${USER}/${USER}/test/mocaseq/work --output_base /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/output --genome_build.human GRCh38.p12 --custom_config_version mocaseq-lrz --custom_config_base https://raw.githubusercontent.com/roland-rad-lab/MoCaSeq/human-pipeline-nextflow-2/conf --input https://raw.githubusercontent.com/roland-rad-lab/MoCaSeq/human-pipeline-nextflow-2/input/mocaseq/PCSI_0683.tsv