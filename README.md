# metaGx_nf
A number of pipelines that are managed by nextflow for conducting various analyses on high-throughput NGS 
# pipelines
There are four ![pipelines](metaGx_pipelines.tiff) implemented in this repository and all can be run using the same main code "metaGx_develope.nf" and choosing a mode for the analyses:
1. Read-based CDS annotation: by aligning reads to a protein database (--mode read-annot)
2. Assembly-based CDS annotation (--mode assembly)
3. Antimicrobial resistance (AMR) annotation (--mode resistome)
4. Specific-gene discovery (--mode read-gene)

# Example: run assembly-based annotation on paired reads
