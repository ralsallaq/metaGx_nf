#!/usr/bin/env nextflow

/* 
 * Run using the following command: 
 * 
 *   nextflow run metaGx.nf 
 * 
 * 
 * The following parameters can be provided as command line options  
 * replacing the prefix `params.` with `--` e.g.:   
 * 
 *   nextflow run metaGx.nf --sampleInfo ./sampleInfo.csv
 *   Ramzi Alsallaq
 * 
 */
/* set parameters to default values */
params.sampleInfo = "./sampleInfo.csv"
params.scratchDir = "/scratch_space/ralsalla/metaGx_nf/"
params.hgenomeSource = "$PWD/humanG/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
params.hgenomeIndexSource = "$PWD/humanG/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz"
//can be "EBV" or "MT" or "both"
params.humanGenomeFilter = "both" 
params.humanGDir = "$PWD/humanG/"
params.adaptersF = "$PWD/adapters.fasta"

params.butyRateCSV = "/research/home/margogrp/ralsalla/butyrate_EC_pathways.csv"
params.bileAcidCSV = "/research/home/margogrp/ralsalla/bileSalt_EC_pathways.csv"

// resistome params
params.megaResDBAnnotCSV = "/research/home/margogrp/ralsalla/databaseFiles/megares_v2.00/megares_drugs_annotations_v2.00.csv"
params.megaResDBFNA = "/research/home/margogrp/ralsalla/databaseFiles/megares_v2.00/megares_drugs_database_v2.00.fasta"
params.resistomeGCC = 80 //resistome gene coverage threshold
params.rarefactionMin=5
params.rarefactionMax=100
params.rarefactionSkip=5
params.numIterationsForSampling=1
params.megaResIndDir = "$PWD/megaRes/"

// kraken params
params.krakenDBDir = "/scratch_space/ralsalla/databases/kraken2_db/"

params.outD = "$PWD/analysis/"
params.mode = null // can be "assembly" for assembly-based annotation, "resistome" for resistome analysis, "read-gene" for read-based specific-gene annotation, or "read-annot" for read-based annotation
params.skipQtrim = false // if true then params.sampleInfo would be considered as the file encompassing information on user provided quality trimmed files 
params.skipTaxClassify = false // if true then kraken step will be skipped 
// read-gene params
params.geneEC = null 
params.geneName = null //must always have a value in read-gene mode
params.geneRetrieve = true // if false then provide seed faa file
params.geneSeedFAA = null //must have a value if geneRetrieve=false
params.geneExtensionDBF = "/research/home/margogrp/ralsalla/metaG/databases/uniref90.fasta.gz"
params.blosum80 = "/home/ralsalla/mmseqs/matrices/blosum80.out"
params.blosum62 = "/home/ralsalla/mmseqs/matrices/blosum62.out"

params.useDiamond = true

// default mode to use
def default_mode = "assembly"

// variable to use throughout the workflow
def runMode

// check if CLI arg was passed; if so used that instead
if(params.mode == null){
    runMode = default_mode
    } else if(params.mode == "resistome" | params.mode == "read-gene" | params.mode == "read-annot" | params.mode == "assembly") {
        runMode = params.mode
} else {
    log.error("No valid mode specified. Valid modes are: 'assembly', 'resistome', 'read-gene', or 'read-annot' ")
    exit 1
}


//make sure everything is fine if the read-gene mode is selected
if(runMode == "read-gene"){
    if(params.geneName == null) {
        log.error("No valid gene name specified. In the read-gene run mode a gene name must be specified")
        exit 1
        }
    if(! params.geneRetrieve) {
        if(params.geneSeedFAA == null) {
            log.error("No valid gene seed FAA file specified. In the read-gene run mode when 'geneRetrieve' is false a gene seed AA fasta file must be specified")
            exit 1
            } else {gene_seedFAA_ch=Channel.fromPath(params.geneSeedFAA)
                    geneName_ch = Channel.from(params.geneName)
                    geneName_ch.merge(gene_seedFAA_ch).set{gene_input_ch}
                    }
    } else {//gene_seedFAA_ch=Channel.from(params.geneSeedFAA)
            //geneName_ch = Channel.from(params.geneName)
            //geneName_ch.merge(gene_seedFAA_ch).set{gene_input_ch}
            gene_input_ch = Channel.empty() //so process setGeneSeed will have an empty input channel
            }
}  

//gene_input_ch.subscribe {println "value: $it"}


// set val("${params.geneName}"), file("${params.geneName}.seed.fasta") into refu_ch 

/* containers */
container_bbtools = "bryce911/bbtools:latest"
container_fastqc = "biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1"
container_cutadapt = "kfdrc/cutadapt:latest"
container_bwa = "biocontainers/bwa:v0.7.17-3-deb_cv1"
container_samtools = "biocontainers/samtools:v1.9-4-deb_cv1"
container_kraken2 = "quay.io/biocontainers/kraken2:2.0.8_beta--pl526hc9558a2_1"
container_refu = "ralsallaq/refu:v0.0"
container_diamond = "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0"
container_mmseqs2 = "soedinglab/mmseqs2:latest"
container_famli = "quay.io/fhcrc-microbiome/famli:v1.5"
container_spades = "biocontainers/spades:v3.13.0dfsg2-2-deb_cv1"
// semi-official prokka (found to be with the best up to date packages)
container_prokka = "bionivid/prokka:latest"
container_eggnog_mapper = "golob/eggnog-mapper:2xx__bcw.0.3.1A"
container_bedtools = "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
container_seqtk = "biocontainers/seqtk:v1.3-1-deb_cv1"
container_picard = "broadinstitute/picard:latest"

log.info """\
METAGX-nf 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
runMode: $runMode
sampleInfo : $params.sampleInfo
outD : $params.outD
scratchDir : $params.scratchDir
adaptersF : $params.adaptersF
butyRateCSV : $params.butyRateCSV
bileAcidCSV : $params.bileAcidCSV
megaResDBFNA : $params.megaResDBFNA
megaResDBAnnotCSV : $params.megaResDBAnnotCSV
resistomeGCC : $params.resistomeGCC
rarefactionMin : $params.rarefactionMin
rarefactionMax : $params.rarefactionMax
rarefactionSkip : $params.rarefactionSkip
numIterationsForSampling : $params.numIterationsForSampling
megaResIndDir : $params.megaResIndDir
resistomeGCC : $params.resistomeGCC
krakenDBDir : $params.krakenDBDir
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
"""
/* process input parameter */
/* create channel from CSV file and create two channels from that*/
Channel
    .fromPath(params.sampleInfo)
    .splitCsv(header: true, sep: ",")
    .into { rawInput1; rawInput2; rawInput3}

rawInput1.map{ sample -> [ sample.sname, file(sample.R1), file(sample.R2)] }.set{noqtsplit_ch1}
rawInput2.map{ sample -> [ sample.sname, file(sample.R1), file(sample.R2)] }.set{noqtsplit_ch2}
rawInput3.map{ sample -> [ sample.sname, file(sample.R1), file(sample.R2)] }.set{noqtsplit_ch3}
// prep adapters fasta file
adaptersF_name = file(params.adaptersF).name
adaptersF_path = file(params.adaptersF).parent
hgenomeSource = file(params.hgenomeSource)
hgenomeIndexSource = file(params.hgenomeIndexSource)

geneExDBFAA=file(params.geneExtensionDBF)
blosum80 = file(params.blosum80)
blosum62 = file(params.blosum62)
kraken2DB = file(params.krakenDBDir)


// check that megaRes DB and annotations are the same
if((params.megaResDBFNA.split("drugs").size()>1 & params.megaResDBAnnotCSV.split("drugs").size()>1) | (params.megaResDBFNA.split("modified").size()>1 & params.megaResDBAnnotCSV.split("modified").size()>1)) {
    megaResDBFNA = file(params.megaResDBFNA)
    megaResDBAnnotCSV = file(params.megaResDBAnnotCSV)
    } else {
    log.error("The assigned megaRes database and annotations differ")
    exit 1
}

/*
// if you uncomment this you will see that the two channels have the same order of files
// One has to make sure that when combining the two results they are combined on the same order, how?
noqtsplit_ch1.subscribe { println "value: $it" }
noqtsplit_ch2.subscribe { println "value: $it" }
*/

/* echo "${queryF.getSimpleName()}" > ${queryF.getSimpleName()}.sname*/


// no need to wget the ftp files just use them directly as nextflow recognize them as files

if(params.skipQtrim == false){

    process filterHumGenome { 
        label 'io_mem'
        publishDir "${params.humanGDir}/fltred", mode: 'copy'
    
        input:
        file hgenome from hgenomeSource
    
        output:
        file 'hgenome.filtered.fna' into FilthumG_ch
    
        """
        module load python/2.7.15-rhel7
        filterHumanG.py \
        -i ${hgenome} \
        -o hgenome.filtered.fna \
        -fltr ${params.humanGenomeFilter}
        """
    }
    
    process bwaIndexHG {
        container "${container_bwa}"
        label 'io_mem'
        publishDir "${params.humanGDir}/bwaInd/", mode: 'copy'
    
        input:
        file hgenomeFltr from FilthumG_ch
    
        output:
        file('hgenome.bwa_index.tgz') into humGBWAIndx 
        file('human_genome_bwa.*') into FiltIndxhumG_ch
    
        """
        bwa index -a bwtsw -p "human_genome_bwa" ${hgenomeFltr}
        wait
        #create an archive for the bwa index for the filtered genome 
        tar czvf hgenome.bwa_index.tgz human_genome_bwa.* 
        """
    }
    
    process verifyAndLinkReads {
        container "${container_bbtools}"
        label 'multithread'
    
        input:
        set sname, file(R1), file(R2) from noqtsplit_ch1
        
        output:
        set sname, file("${sname}.reads.fq.gz") into vpair_ch
        
        """
        reformat.sh in1=${R1} in2=${R2} vpair out=${sname}.reads.fq.gz  
        """
    }
    
    
    // vpair_ch.subscribe { println "value: $it" }
    
    
    vpair_ch.into { vpair_ch1; vpair_ch2; vpair_ch3}
    
    
    process chckFQ {
        container "${container_fastqc}"
        label 'multithread'
        publishDir "${params.outD}/qualityChecks/", mode: 'copy', pattern: '*.html' 
    
        input:
        set sname, file(intrlF) from vpair_ch1
        //output:
        //file "${intrlF.getSimpleName()}_fastqc.html"
    
        """
        fastqc -t 4 ${intrlF} 
        """
    }
    
    
    process removeAdapters {
        container "${container_cutadapt}"
        label 'multithread'
        publishDir "${params.outD}/noAdpt/", mode: 'copy'
        
        input:
        set sname, file(intrlF) from vpair_ch2
        file adaptersF_path
    
        output:
        set sname, file("${sname}.noAdpt.fq.gz") into noAdpt_ch
    
        """
        cutadapt --interleaved -j ${task.cpus}  -a file:$adaptersF_path/$adaptersF_name -A file:$adaptersF_path/$adaptersF_name -o ${sname}.noAdpt.fq.gz ${intrlF}
        """
    }
    
    noAdpt_ch.into { noAdpt_ch1; noAdpt_ch2 }
    
    process filterArtifacts {
        container "${container_bbtools}"
        label 'multithread'
        publishDir "${params.outD}/noArtiFacts/", mode: 'copy'
        
        input:
        set sname, file(noAdptF) from noAdpt_ch1
    
        output:
        set sname, file("${sname}.noArtiF.fq.gz"), file("${sname}.noArtiF.stats.txt") into noArtiF_ch
        set sname, file("${sname}.noArtiF.fq.gz") into noArtiF_ch1, noArtiF_ch2
    
        """
        bbduk.sh in=${noAdptF} out=${sname}.noArtiF.fq.gz ref=/bbmap/resources/phix174_ill.ref.fa.gz,/bbmap/resources/phix_adapters.fa.gz,/bbmap/resources/sequencing_artifacts.fa.gz,artifacts,phix  k=31 hdist=1 ordered cardinality stats=${sname}.noArtiF.stats.txt 
        """
    }
    
    // noArtiF_ch.into { noArtiF_ch1; noArtiF_ch2 }
    
    
    // buffer it so it will not be processed as separate files
    FiltIndxhumG_ch.flatMap { it -> it[0].parent/it[0].getSimpleName() }.buffer( size:1).set{idxbase_ch}
    // get a file and then its path
    idxbase_parent = idxbase_ch.first().parent
    
    
    process decontByMap1 {
        container "${container_bwa}"
        label 'multithread'
        // publishDir "${params.outD}/decontaminated/", mode: 'copy'
        
        input:
        set sname, file(noArtiF) from noArtiF_ch1
        path idxbase_parent
    
        output:
        set sname, file("${sname}.alignment.sam") into alnToHuman_ch
    
       // script:
       // def idxbase = hgenomBWAIndexFs.toList()[0].getSimpleName()
        
        """
        #align noArtiFile to human genome
        bwa mem -t ${task.cpus} -o ${sname}.alignment.sam -p ${idxbase_parent}/human_genome_bwa ${noArtiF} 
        """
    }
    
    
    process decontByMap2 {
        container "${container_samtools}"
        label 'multithread'
        publishDir "${params.outD}/decontaminated/", mode: 'copy'
        
        input:
        set sname, file(samalnToH) from alnToHuman_ch
    
        output:
        set sname, file("${sname}.decontaminated.fq.gz") into noHuman_ch
    
        """
        #convert to bam
        samtools view -bh -o ${sname}.alignment.bam ${samalnToH}
        #extract unaligned pairs (i.e. non-human sequences)
        samtools fastq  --threads ${task.cpus}  -f 12 ${sname}.alignment.bam  > ${sname}.decontaminated.fq
        #create a compressed file
        gzip -c ${sname}.decontaminated.fq > ${sname}.decontaminated.fq.gz
        """
    }
    
    noHuman_ch.into { noHuman_ch1; noHuman_ch2 }
    
    process qualityTrim {
        container "${container_cutadapt}"
        label 'multithread'
        publishDir "${params.outD}/qtrimmed/", mode: 'copy'
        
        input:
        set sname, file(decontF) from noHuman_ch1
    
        output:
        set sname, file("${sname}.qtrimmed.fq.gz") into qtrimmed_ch1, qtrimmed_ch2, qtrimmed_ch3, qtrimmed_ch4, qtrimmed_ch5, qtrimmed_ch6, qtrimmed_ch7, qtrimmed_ch8
    
        """
        cutadapt --interleaved -j ${task.cpus} -q 10 -m 50 -A XXX -o ${sname}.qtrimmed.fq.gz ${decontF}
        """
    }
    
    
    vpair_ch3.join(noAdpt_ch2).set{comb_ch1}
    comb_ch1.join(noArtiF_ch2).set{comb_ch2}
    comb_ch2.join(noHuman_ch2).set{comb_ch3}
    comb_ch3.join(qtrimmed_ch1).set{comb_ch4}
    // comb_ch4.subscribe  {println "value: $it"} 
    
    
    process trackReads {
        label 'io_mem'
        //publishDir "${params.outD}/", mode: 'copy'
        
        input:
        set sname, file(inpF), file(noAdptF), file(noArtiF), file(decontF), file(qtrimF) from comb_ch4
    
        output:
        file result 
    
        """
        inputN=\$(echo \$(zcat $inpF|wc -l)/4|bc)
        noAdptN=\$(echo \$(zcat $noAdptF|wc -l)/4|bc)
        filtN=\$(echo \$(zcat $noArtiF|wc -l)/4|bc)
        decontN=\$(echo \$(zcat $decontF|wc -l)/4|bc)
        qtrimN=\$(echo \$(zcat $qtrimF|wc -l)/4|bc)
    
        echo "#sample  #input  #noAdpt  #filtered  #decontaminated  #qtrimmed"|column -t > result 
        echo "${sname} \$inputN \$noAdptN \$filtN \$decontN \$qtrimN" |column -t >> result
    
        """
    }
    
    //collect result files into one file
    result
        .collectFile(name: "${params.outD}/trackReads.out", keepHeader:true, skip:1)
        .println { file -> "sample\tinput\tnoAdpt\tfilt\tdecont\tqtrim:\n ${file.text}" }
    
    /* Now that we have qtrimmed reads */
    
    process splitFiles {
        container "${container_bbtools}"
        label 'multithread'
    
        input:
        set sname, file(qtrimmedF) from qtrimmed_ch2
    
        output:
        set sname, file("${sname}.qt.R1.fq.gz"), file("${sname}.qt.R2.fq.gz") into qtsplit_ch

        // only execute process if params.skipTaxClassify is set to false
        when:
        !params.skipTaxClassify
    
        """
        reformat.sh in=${qtrimmedF}  out1=${sname}.qt.R1.fq.gz out2=${sname}.qt.R2.fq.gz
        """
    }

    process taxaClassifyQt {
        container "${container_kraken2}"
        label 'multithread'
        publishDir "${params.outD}/kraken2/", mode: 'copy'
        
        input:
        //mix the two channels : qtrimmed and noqtrimmed to account for either choice by the user
        set sname, file(R1), file(R2) from qtsplit_ch 
        path kraken2DB
    
        output:
        set sname, file("${sname}.krakenReport"), file("${sname}.krakenClass"), file("${sname}.krakenLabels") into kraken_ch
    
        // only execute process if params.skipTaxClassify is set to false
        when:
        !params.skipTaxClassify
    
        """
        kraken2 --use-names  --threads ${task.cpus} --db ${kraken2DB}  --gzip-compressed --paired  ${R1} ${R2} --report ${sname}.krakenReport > ${sname}.krakenClass
        #produce krakenLabels
        grep '^C' ${sname}.krakenClass | sed 's/^C\t//g' >  ${sname}.krakenLabels
        """
    }

} else {//skip qtrim need to consider sampleInfo file as a quality trimmed splitted files in noqtsplit_ch2/noqtsplit_ch3 channels

    process verifyAndLinkNQtrimmedReads {
        container "${container_bbtools}"
        label 'multithread'
    
        input:
        set sname, file(R1), file(R2) from noqtsplit_ch2
        
        output:
        set sname, file("${sname}.qtrimmed.fq.gz") into qtrimmed_ch1, qtrimmed_ch2, qtrimmed_ch3, qtrimmed_ch4, qtrimmed_ch5, qtrimmed_ch6, qtrimmed_ch7, qtrimmed_ch8
        
        """
        reformat.sh in1=${R1} in2=${R2} vpair out=${sname}.qtrimmed.fq.gz  
        """
    }


    process taxaClassifyNQt {
        container "${container_kraken2}"
        label 'multithread'
        publishDir "${params.outD}/kraken2/", mode: 'copy'
        
        input:
        //mix the two channels : qtrimmed and noqtrimmed to account for either choice by the user
        set sname, file(R1), file(R2) from noqtsplit_ch3 
        path kraken2DB
    
        output:
        set sname, file("${sname}.krakenReport"), file("${sname}.krakenClass"), file("${sname}.krakenLabels") into kraken_ch
    
        // only execute process if params.skipTaxClassify is set to false
        when:
        !params.skipTaxClassify
    
        """
        kraken2 --use-names  --threads ${task.cpus} --db ${kraken2DB}  --gzip-compressed --paired  ${R1} ${R2} --report ${sname}.krakenReport > ${sname}.krakenClass
        #produce krakenLabels
        grep '^C' ${sname}.krakenClass | sed 's/^C\t//g' >  ${sname}.krakenLabels
        """
    }
}



if(runMode == "resistome") 
{ // starting the block for resistome analyses

    process bwaIndexMegaRes {
        container "${container_bwa}"
        label 'io_mem'
        publishDir "${params.megaResIndDir}/", mode: 'copy'
    
        input:
        file megaResDBFNA 
    
        output:
        file('megaRes_bwa.*') into indxMegaRes_ch
    
        """
        bwa index -a bwtsw -p "megaRes_bwa" ${megaResDBFNA}
        """
    }

     // buffer it so it will not be processed as separate files
     indxMegaRes_ch.flatMap { it -> it[0].parent/it[0].getSimpleName() }.buffer( size:1).set{megaResIdxbase_ch}
     // get a file and then its path
     megaResIdxbase_parent = megaResIdxbase_ch.first().parent
    
    process alnReadsToMegaRes {
        container "${container_bwa}"
        label 'multithread'

        input:
        set sname, file(qtrimmedF) from qtrimmed_ch3
        path megaResIdxbase_parent

        output:
        set sname, file("${sname}.megaResalgn.sam") into alnToMegaRes_ch1, alnToMegaRes_ch2, alnToMegaRes_ch3
        
        """
        bwa mem -t ${task.cpus} -o ${sname}.megaResalgn.sam -p ${megaResIdxbase_parent}/megaRes_bwa ${qtrimmedF} 
        """

    }


    process resistomeAnalysis {
        // container "${container_cutadapt}"
        label 'multithread'
        publishDir "${params.outD}/AMR/", mode: 'copy'
        
        input:
        set sname, file(samalnToMegaRes) from alnToMegaRes_ch1
    
        output:
        set sname, file("${sname}.gene.tsv"), file("${sname}.group.tsv"), file("${sname}.class.tsv"), file("${sname}.mechanism.tsv") into amr_ch
    
        """
        resistome \
             -ref_fp ${params.megaResDBFNA} \
             -annot_fp ${params.megaResDBAnnotCSV} \
             -sam_fp ${samalnToMegaRes} \
             -gene_fp ${sname}.gene.tsv \
             -group_fp ${sname}.group.tsv \
             -class_fp ${sname}.class.tsv \
             -mech_fp ${sname}.mechanism.tsv \
             -t ${params.resistomeGCC}
        """
    }

    process rareFactionAnalysis {
        // container "${container_cutadapt}"
        label 'multithread'
        publishDir "${params.outD}/AMR/", mode: 'copy'
        
        input:
        set sname, file(samalnToMegaRes) from alnToMegaRes_ch2
    
        output:
        set sname, file("${sname}.gene.rareFD.tsv"), file("${sname}.group.rareFD.tsv"), file("${sname}.class.rareFD.tsv"), file("${sname}.mech.rareFD.tsv") into amr_rareF_ch
    
        """
        rarefaction \
             -ref_fp ${params.megaResDBFNA} \
             -sam_fp ${samalnToMegaRes} \
             -annot_fp ${params.megaResDBAnnotCSV} \
             -gene_fp ${sname}.gene.rareFD.tsv \
             -group_fp ${sname}.group.rareFD.tsv \
             -class_fp ${sname}.class.rareFD.tsv \
             -mech_fp ${sname}.mech.rareFD.tsv \
             -min ${params.rarefactionMin} \
             -max ${params.rarefactionMax} \
             -skip ${params.rarefactionSkip} \
             -samples ${params.numIterationsForSampling} \
             -t ${params.resistomeGCC} 
        """
    }

    process snpFinder {
        // container "${container_cutadapt}"
        label 'multithread'
        publishDir "${params.outD}/AMR/", mode: 'copy'
        
        input:
        set sname, file(samalnToMegaRes) from alnToMegaRes_ch3
    
        output:
        set sname, file("${sname}.snp.tsv") into amr_snp_ch
    
        """
         snpfinder \
             -amr_fp ${params.megaResDBFNA} \
             -sampe ${samalnToMegaRes} \
             -out_fp ${sname}.snp.tsv

        """
    }


}//end of resistome analyses 

/*===================================*/


else if(runMode == "read-gene") 
{ // starting the block for read-based specific gene annotation


    process setGeneSeed {
        label 'io_mem'
        publishDir "${params.outD}/gene_discover/", mode: 'copy'
        
        input:
        set geneName, seedF from gene_input_ch 
    
        output:
        set geneName, file("${geneName}.seed.fasta") into notrefu_ch 

        // only execute process if geneRetrieve is set to false; might be redundant as gene_input_ch is empty when params.geneRetrieve=true
        //when:
        //!params.geneRetrieve
    
        """
        cp ${seedF} ${geneName}.seed.fasta
        """
    }


    process getGeneSeeds {
        container "${container_refu}"
        label 'io_mem'
        publishDir "${params.outD}/gene_discover/", mode: 'copy'
    
        output:
        set val("${params.geneName}"), file("${params.geneName}.seed.fasta") into refu_ch 

        // only execute process if geneRetrieve is set to true; might be redundant as process setGeneSeed will return empty channel when params.geneRetrieve=true
        //when:
        //params.geneRetrieve
    
        """
        getEnzFromUniprot.py -ec ${params.geneEC} -g ${params.geneName} --reviewed -o ${params.geneName}.seed.fasta 
        """
    }


    process alnToGeneExDB {
        container "${container_mmseqs2}"
        label 'multithread'
        publishDir "${params.outD}/gene_discover/${geneName}/", mode: 'copy'
        
        input:
        set geneName, file(seedF) from refu_ch.mix(notrefu_ch)
        file geneExDBFAA
        file blosum62 
        
        output:
        set geneName, file("result_${geneName}.m8") into extDBAln_ch
        
        """
        mkdir tmp
        #search the geneExDBFAA (database) for hits for the seedF at probability of substitution of 62
        mmseqs easy-search ${seedF}  ${geneExDBFAA} result_${geneName}.m8 tmp --sub-mat ${blosum62} --threads ${task.cpus} \
        --start-sens 1 --sens-steps 3 -s 7 -e 0.00001 --remove-tmp-files true \
        --format-output "query,target,qlen,tlen,qstart,qend,tstart,tend,gapopen,mismatch,pident,evalue,bits"  

        """
    }


    process getExtendedFAA {
        container "${container_seqtk}"
        label 'io_mem'
        publishDir "${params.outD}/gene_discover/${geneName}/", mode: 'copy'
        
        input:
        set geneName, file(blastLikeAlnF) from extDBAln_ch
        set geneName, file(seedF) from refu_ch
        file geneExDBFAA
        
        output:
        set geneName, file("${geneName}.extended.faa") into geneExtFAA_ch
        
        """
        # get target ids
        cat "${blastLikeAlnF}" | awk '{print \$2}' > target_ids
        # remove duplicates
        awk '!seen[\$1]++'  target_ids >  target_ids_nodup
        # use seqtk to subseq the database
        seqtk subseq ${geneExDBFAA} target_ids_nodup > ${geneName}.extended.faa 
		
        """
    }


    process getRepFAA {
        container "${container_mmseqs2}"
        label 'multithread'
        publishDir "${params.outD}/gene_discover/${geneName}/", mode: 'copy'
        
        input:
        set geneName, file(extFAA) from geneExtFAA_ch
        
        output:
        set geneName, file("${geneName}_clstr_rep_seq.fasta"), file("${geneName}_clstr_all_seqs.fasta"), file("${geneName}_clstr_cluster.tsv") into clustGene_ch
        
        """
        mkdir tmp
        #search the geneExDBFAA (database) for hits for the seedF at probability of substitution of 62
        #cluster to 80% of sequence identity, and asking for the representative sequence that cover at least 90% of members 
        #and the members cover at least 90% of it, this seems comparable if not better than CD-HIT
        
        mmseqs easy-cluster ${extFAA} ${geneName}_clstr tmp --min-seq-id 0.8 -c 0.90 --cov-mode 0 --remove-tmp-files true

        """
    }


    clustGene_ch.combine(qtrimmed_ch4).into{findHits_ch1; findHits_ch2}

    findHits_ch1.subscribe{println it}

    process findHits {
        container "${container_mmseqs2}"
        label 'multithread'
        publishDir "${params.outD}/gene_discover/${geneName}/", mode: 'copy'
        
        input:
        set geneName, file(clstrRepFAA), file(allFAA), file(clstrTab), sname, file(qtrimmedF) from findHits_ch2
        file blosum80 
        
        output:
        set geneName, sname, file("${geneName}_hits_in_${sname}.m8") into geneHits_ch
        
        """
        mkdir tmp
        #mmseqs easy-search ${qtrimmedF} ${clstrRepFAA} ${geneName}_hits_in_${sname}.m8 tmp --sub-mat ${blosum80} --threads ${task.cpus} --remove-tmp-files true \
        #   -s 5.5  --format-output "query,target,alnlen,qlen,tlen,qstart,qend,tstart,tend,pident,nident,mismatch,qcov,tcov,evalue,bits"  
        echo "the format of the output is set to be able to run famli fltr next"
        mmseqs easy-search ${qtrimmedF} ${clstrRepFAA} ${geneName}_hits_in_${sname}.m8 tmp --sub-mat ${blosum80} --threads ${task.cpus} --remove-tmp-files true \
        -s 5.5  --format-output "query,target,alnlen,qlen,tlen,qstart,qend,tstart,tend,pident,nident,mismatch,qcov,tcov,evalue,bits"  
        #--start-sens 1 --sens-steps 3 -s 7  --format-output "query,target,alnlen,qlen,tlen,qstart,qend,tstart,tend,pident,nident,mismatch,qcov,tcov,evalue,bits"  

        """
    }


    process filterHits {
        container "${container_famli}"
        label 'multithread'
        publishDir "${params.outD}/gene_discover/${geneName}/", mode: 'copy'
        errorStrategy 'ignore'
        
        input:
        set geneName, sname, file(hitsInsample) from geneHits_ch
        
        output:
        set geneName, sname, file("${sname}.json.gz") into filtrdHits_ch
        
        """
        if [[ `cat ${hitsInsample} | wc -l` -gt 0 ]]; then
        famli filter --input ${hitsInsample} --output ${sname}.json.gz  --threads ${task.cpus} \
        --output-aln ${sname}.aln --logfile ${sname}.log --qseqid-ix 0 --sseqid-ix 1 \
        --sstart-ix 7 --send-ix 8 --bitscore-ix 13 --slen-ix 4 --sd-mean-cutoff 1.0 --strim-5 5 --strim-3 5
        else
        echo "no hits to be filtered! reporting empty file"
        touch ${sname}.json.gz
        fi
        """
    }


    process filterHitsCSV {
        label 'io_limited'
        publishDir "${params.outD}/gene_discover/${geneName}/", mode: 'copy'
        errorStrategy 'ignore'
        
        input:
        set geneName, sname, file(hitsInjason) from filtrdHits_ch
        
        output:
        file("${geneName}_hits_in_${sname}.csv")
        
        """
        module load python/2.7.15-rhel7
        if [[ `cat ${hitsInjason} | wc -l` -gt 0 ]]; then
        filterHitsToCSV.py -j ${hitsInjason} -o ${geneName}_hits_in_${sname}.csv 
        else
        echo "no hits to report; empty file found! reporting empty file"
        touch ${geneName}_hits_in_${sname}.csv
        fi
        """
    }

}//end of read-gene analyses

/*===================================*/

else if(runMode == "assembly")
{ // starting the block for assembly-based analyses

    process assembleToContigs {
        container "${container_spades}"
        label 'multithread'
        publishDir "${params.outD}/assembly-based/scaffolds/", mode: 'copy'
        
        input:
        set sname, file(qtrimmedF) from qtrimmed_ch5
        
        output:
        set sname, file("${sname}_scaffolds.fasta") into spades_ch
        
        """
        spades.py --meta --only-assembler -t ${task.cpus} --phred-offset 33 -k 25,55,95  --12 ${qtrimmedF}  -o spades_out.${sname} 
        wait
        cp spades_out.${sname}/scaffolds.fasta ${sname}_scaffolds.fasta

        """
    }

    process predictCDS {
        //container "${container_prokka}"
        label 'multithread'
        publishDir "${params.outD}/assembly-based/CDS/", mode: 'copy'
        errorStrategy 'ignore'
        
        input:
        set sname, file(scaffoldFNA) from spades_ch
        
        output:
        //set sname, file("prokka_Annot.${sname}/*") into prokka_ch1, prokka_ch2, prokka_ch3
        set sname, file("prokka_Annot.${sname}/${sname}.prokka.faa") into prokka_ch1, prokka_ch2, prokka_ch3
        
        """
        set +eu
        activate=/hpcf/apps/conda2/install/4.1.1/bin/activate
        echo "activate prokka env to get prokka"
        source \$activate prokkaLatest
        prokka --outdir prokka_Annot.$sname --force --centre metaG --norna --notrna --metagenome --prefix ${sname}.prokka --cpus ${task.cpus} ${scaffoldFNA} 

        """
    }


    qtrimmed_ch6.join(prokka_ch1).set{qtrim_prokka_ch}


    process getEggNogData {
        container "${container_eggnog_mapper}"
        label 'io_mem'
        //errorStrategy 'retry'
        //maxRetries 4

        output:
        path("eggNogData") into eggnogData_ch

        """
        mkdir eggNogData
        download_eggnog_data.py -y --data_dir eggNogData    
        """
    }
 
    //convert to file so we do not need to emit within the process
    eggNogData_path = eggnogData_ch.first()
 

    process annotateCDS {
        container "${container_eggnog_mapper}"
        label 'multithread'
        publishDir "${params.outD}/assembly-based/CDS_annot/", mode: 'copy'
        errorStrategy 'ignore'
        //errorStrategy 'retry'
        //maxRetries 4
        
        input:
        set sname, file(cdsFAA) from prokka_ch2
        path eggNogData_path
        
        output:
        set sname, file("${sname}.emapper.annotations"), file("${sname}.emapper.annotations.reformat") into eggnog_ch1, eggnog_ch2
        
        """
        emapper.py --resume -i ${cdsFAA} -m diamond --data_dir ${eggNogData_path} --cpu ${task.cpus} -o ${sname}
        wait
        grep -v '^#' ${sname}.emapper.annotations |awk 'BEGIN{FS="\t"}{print \$0}'>${sname}.emapper.annotations.reformat
        """
    }
    
    

    process mapReadsToCDS {
        container "${container_diamond}" //<--bwa only does DNA-DNA alignment
        label 'multithread'
        //publishDir "${params.outD}/assembly-based/CDS/readsOntoCDS/", mode: 'copy'
        
        input:
        set sname, file(qtrimmedF), file(cdsFAA) from qtrim_prokka_ch
        
        output:
        set sname, file("${sname}.mapOnTo.CDS.sam") into readsOntoCDS_sam_ch1, readsOntoCDS_sam_ch2
        
        """
        #MMSEQS2 support for SAM format is still beta and breaks
        #mkdir tmp
        #mmseqs createdb ${cdsFAA}  targetDB --dbtype 0 --shuffle 1 --createdb-mode 0 --id-offset 0 --compressed 0 -v 3
        #mmseqs search -a --alignment-mode 4  --threads ${task.cpus} -s 2 -c 0.9 --cov-mode 2 
        #echo "the format of the output is set to SAM format the AS tag contains the raw score, NM is the miss match count."
        #mmseqs easy-search ${qtrimmedF} ${cdsFAA} ${sname}.mapOnTo.CDS.sam tmp --threads ${task.cpus} --alignment-mode 4 -s 2 \
        #-e 0.00001 -c 0.9 --cov-mode 1 --remove-tmp-files true --format-mode 1 
        diamond makedb --in ${cdsFAA} -d prokka
        diamond blastx --threads ${task.cpus} --query ${qtrimmedF} --db prokka.dmnd  --outfmt 101  -o ${sname}.mapOnTo.CDS.sam
        """
    }


    readsOntoCDS_sam_ch1.join(prokka_ch3).set{dmndSam_prokka_ch}

    process samTobam {
        container "${container_samtools}"
        label 'multithread'
        publishDir "${params.outD}/assembly-based/readsOntoCDS/", mode: 'copy'
        
        input:
        set sname, file(samAln), file(cdsFAA) from dmndSam_prokka_ch
        
        output:
        set sname, file("${sname}.mapOnTo.CDS.bam"), file("${sname}.mapOnTo.CDS.sorted.bam") into readsOntoCDS_bam_ch
        
        """
        samtools faidx ${cdsFAA}
        wait
        faiF=`ls *.fai`
        samtools view -bh -t \$faiF -o ${sname}.mapOnTo.CDS.bam  ${samAln}
        samtools sort -@ ${task.cpus} ${sname}.mapOnTo.CDS.bam -o ${sname}.mapOnTo.CDS.sorted.bam
        """
    }


//    not necessary
//    process dropDuplicates {
//        container "${container_picard}"
//        label 'multithread'
//        publishDir "${params.outD}/assembly-based/readsOntoCDS/", mode: 'copy'
//        
//        input:
//        set sname, file(unsortedBam), file(sortedBam) from readsOntoCDS_bam_ch
//        
//        output:
//        set sname, file("${sname}.mapOnTo.CDS.markedup.bam"), file("${sname}.mapOnTo.CDS.markedup.metrics") into readsOntoCDS_nodup_ch
//        
//        """
//        java -jar /usr/picard/picard.jar MarkDuplicates INPUT=${sortedBam}  OUTPUT=${sname}.mapOnTo.CDS.markedup.bam METRICS_FILE=${sname}.mapOnTo.CDS.markedup.metrics AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=60000 REMOVE_DUPLICATES=TRUE
//
//        """
//    }


    process quantifyCDSCov {
        container "${container_bedtools}"
        label 'multithread'
        publishDir "${params.outD}/assembly-based/readsOntoCDS/", mode: 'copy'
        
        input:
        set sname, file(unsortedBam), file(sortedBam) from readsOntoCDS_bam_ch
        
        output:
        set sname, file("${sname}.cds_coverage.hist") into cdsCovHisto_ch
        
        """
        bedtools genomecov -ibam ${sortedBam}  > ${sname}.cds_coverage.hist
        """
    }


    eggnog_ch1.join(cdsCovHisto_ch).set{annot_cov_ch}

    process quantifyCDS {
        //conda "envs/python.yml" //define conda env for this step to get necessary modules loaded
        label 'io_limited'
        publishDir "${params.outD}/assembly-based/readsOntoCDS/", mode: 'copy'
        
        input:
        set sname, file(eggNogAnnot), file(eggNogAnnotReformat), file(covHisto) from annot_cov_ch
        
        
        output:
        set sname, file("${sname}.cds_annot.Wcov.csv") into cdsAnnotCov_ch
        
        """
        module load python/2.7.15-rhel7
        addCovToAnnot.py ${eggNogAnnotReformat} ${covHisto} ${sname}.cds_annot.Wcov.csv
        """
    }




}//end of assembly-based analyses

/*===================================*/

else if(runMode == "read-annot")
{ // starting the block for read-based annotation


    process buildDiamondDB {
        container "${container_diamond}"
        label 'io_mem'
        //publishDir "${params.outD}/uniref/", mode: 'copy'
        
        input:
        file geneExDBFAA
        
        output:
        file('targetDB.dmnd') into targetDmnd_ch

        // only execute process if useDiamond is set to true
        when:
        params.useDiamond
        
        """
        diamond makedb --in ${geneExDBFAA}  -d targetDB.dmnd
        """
    }

    process buildMMseqDB {
        container "${container_mmseqs2}"
        label 'io_mem'
        //publishDir "${params.outD}/uniref/", mode: 'copy'
        
        input:
        file geneExDBFAA

        // only execute process if useDiamond is set to false
        when:
        ! params.useDiamond
        
        output:
        file('targetDBMMseq*') into targetMMseq_ch
        
        """
        mmseqs createdb ${geneExDBFAA}  targetDBMMseq --dbtype 0 --shuffle 1 --createdb-mode 0 --id-offset 0 --compressed 0 -v 3
        """
    }

    // get the path to the db files so it will not be processed as separate files
    if(!params.useDiamond) { 
    targetDmnd_ch.mix(targetMMseq_ch).flatMap { it -> it[0].parent/it[0].getSimpleName() }.buffer( size:1).into {targetDB_ch1; targetDB_ch2}
    targetDB_parent = targetDB_ch1.first().parent
                            } else { targetDB_parent = targetDmnd_ch.first()} 

    process alnReadsToDBWdmnd {
       // container "${container_mmseqs2}"
        container "${container_diamond}"
        label 'multithread'
        publishDir "${params.outD}/read_based_annotation/", mode: 'copy'
        
        input:
        set sname, file(qtrimmedF) from qtrimmed_ch7
        path targetDB_parent 
        
        output:
        set sname, file("all_hits_in_${sname}.d6") into sampleHitsDmnD_ch

        // only execute process if useDiamond is set to true
        when:
        params.useDiamond
        
        """
        #diamond way
        diamond blastx --query ${qtrimmedF} --out all_hits_in_${sname}.d6 --threads ${task.cpus} \
        --db ${targetDB_parent} --outfmt 6 \
        qseqid sseqid length qlen slen qstart qend sstart send pident mismatch gapopen  evalue bitscore \
        --min-score 20 --query-cover 95 --id 80 --top 10 --block-size 5 --query-gencode 11 --unal 0

        """
    }

    process alnReadsToDBWmmseq {
        container "${container_mmseqs2}"
        label 'multithread'
        publishDir "${params.outD}/read_based_annotation/", mode: 'copy'
        maxForks = 4
        
        input:
        set sname, file(qtrimmedF) from qtrimmed_ch8
        path targetDB_parent 
        
        output:
        set sname, file("all_hits_in_${sname}.m8") into sampleHitsMMseq_ch

        // only execute process if useDiamond is set to false
        when:
        ! params.useDiamond
        
        """
        #mmseqs2 way
        mkdir tmp
        mmseqs createdb ${qtrimmedF} queryDB --dbtype 0 --shuffle 1 --createdb-mode 0 --id-offset 0 --compressed 0 -v 3 
        wait

        mmseqs search queryDB ${targetDB_parent}/targetDBMMseq alnDB tmp \
                -a  --alignment-mode 4 --threads ${task.cpus} \
                -s 2 --cov-mode 2 -c 0.95 -e 0.00001 \
                --remove-tmp-files true  
        wait

        echo "the format of the output is set to be able to run famli fltr next"
        mmseqs convertalis queryDB ${targetDB_parent}/targetDBMMseq alnDB all_hits_in_${sname}.m8 \
                --format-output "query,target,alnlen,qlen,tlen,qstart,qend,tstart,tend,pident,mismatch,gapopen,evalue,bits,qcov,tcov"

        """
    }

    process filterHitAnnot {
        container "${container_famli}"
        label 'multithread'
        publishDir "${params.outD}/read_based_annotation/filtered", mode: 'copy'
        errorStrategy 'ignore'
        maxForks = 4
        
        input:
        set sname, file(hitsInsample) from sampleHitsDmnD_ch.mix(sampleHitsMMseq_ch)
        
        output:
        set sname, file("${sname}.json.gz") into filtrdHits_ch
        
        """
        famli filter --input ${hitsInsample} --output ${sname}.json.gz  --threads ${task.cpus} \
        --output-aln ${sname}.aln --logfile ${sname}.log --qseqid-ix 0 --sseqid-ix 1 \
        --sstart-ix 7 --send-ix 8 --bitscore-ix 13 --slen-ix 4 --sd-mean-cutoff 1.0 --strim-5 5 --strim-3 5

        """
    }

    process filterHitAnnotCSV {
        label 'io_limited'
        publishDir "${params.outD}/read_based_annotation/filtered", mode: 'copy'
        errorStrategy 'ignore'
        maxForks = 4
        
        input:
        set sname, file(hitsInjason) from filtrdHits_ch
        
        output:
        set sname, file("hits_in_${sname}.csv") into filtrdHitsCSV_ch1, filtrdHitsCSV_ch2
        
        """
        module load python/2.7.15-rhel7
        filterHitsToCSV.py -j ${hitsInjason} -o hits_in_${sname}.csv
        """
    }


    process addAnnotFromID {
        label 'io_mem'
        publishDir "${params.outD}/read_based_annotation/filtered", mode: 'copy'
        errorStrategy 'ignore'
        
        input:
        //set sname, file("hitsInCSV"), file("hitsInFaa") from filtrdHitsComb_ch
        set sname, file("hitsInCSV") from filtrdHitsCSV_ch2
        
        output:
        set sname, file("hits_in_${sname}_annot.csv") into filtrdHits_annot_ch 
        
        """
        module load python/2.7.15-rhel7
        addAnnotToFilteredHits_fromUniProtAccess.py -i ${hitsInCSV} -o  hits_in_${sname}_annot.csv
		
        """
    }

}//end of read-annot

/*===================================*/

else
{ // raise error
    log.error("No valid mode specified. Valid modes are: 'assembly', 'resistome', 'read-gene', or 'read-annot' ")
    exit 1
}
