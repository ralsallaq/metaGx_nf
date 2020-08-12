# metaGx_nf
A number of pipelines that are managed by nextflow for conducting various analyses on high-throughput NGS 
# pipelines
There are four ![pipelines](metaGx_pipelines.tiff) implemented in this repository and all can be run using the same main code "metaGx_develope.nf" and choosing a mode for the analyses:
1. Read-based CDS annotation: by aligning reads to a protein database (--mode read-annot)
2. Assembly-based CDS annotation (--mode assembly)
3. Antimicrobial resistance (AMR) annotation (--mode resistome)
4. Specific-gene discovery (--mode read-gene)

# Example: run assembly-based annotation on paired reads
prepare a sample information CSV file (e.g. sampleInfo.csv) with the following headings: 
R1=absolute path to forward reads
R2=absolute path to reverse reads
sname=sample unique identifier that will be used through the pipeline
## case 1 perform quality trimming first: 
./nextflow run metaGx_develope.nf -w /path/to/scratch/directory --sampleInfo sampleInfo.csv --mode assembly -resume
## case 2 the paths in sampleInfo.csv points to quality-trimmed files 
./nextflow run metaGx_develope.nf -w /path/to/scratch/directory --sampleInfo sampleInfo.csv --mode assembly --skipQtrim true -resume
