// WINK (What's In my Nanopore reads, with Kraken2) pipeline
// real-time monitoring of fastq_pass folder, execute kraken+bracken (optionally kaiju), shiny app on top

if( !nextflow.version.matches('>=19.08') ) {
    println "This workflow requires Nextflow version 19.08 or greater and you are running version $nextflow.version"
    exit 1
}

/* 
 * pipeline input parameters 
 */
 params.fastq_pass = "fastq_pass"
 params.results = "${workflow.launchDir}/results-wink"
 params.kraken_db = "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz"
 params.weakmem = false
 params.taxlevel = "S" //level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
 params.help = false

def helpMessage() {
log.info """
        ===========================================
         W I N K - What's In my Nanopore reads, with Kraken2

         real-time monitoring of a Nanopore run + taxonomic and abundance information

         Parameters:
        -------------------------------------------
         --fastq_pass       : the folder where basecalled reads are saved during a run, must contain barcodes
         --results          : where the results will go
         --kraken_db        : 
         --taxlevel         : 
         """
         .stripIndent()
}

if (params.help){
    helpMessage()
    exit(0)
}

log.info """
        ===========================================
         W I N K - What's In my Nanopore reads, with Kraken2

         real-time monitoring of a Nanopore run + taxonomic and abundance information

         Used parameters:
        -------------------------------------------
         --fastq_pass       : ${params.fastq_pass}
         --results          : ${params.results}
         --kraken_db        : ${params.kraken_db}
         --taxlevel         : ${params.taxlevel}

         Runtime data:
        -------------------------------------------
         Running with profile:   ${workflow.profile}
         Container:              ${workflow.container}
         Running as user:        ${workflow.userName}
         Launch dir:             ${workflow.launchDir}
         Base dir:               ${baseDir}
         """
         .stripIndent()


// used to touch the files at start so that all available fastq files are processed 
 Channel
    .fromPath( "${params.fastq_pass}/**.fastq", type: 'file')
    //.ifEmpty { error "Can not find any barcode* folders in ${params.fastq_pass}" }
    //.view()
    .set { fastq_ch }


Channel
    .watchPath("${params.fastq_pass}/**.fastq", 'create,modify')
    //.view()
    .set { watch_ch }

// touch when program starts so that all fastq are processed
process touch {
    input:
        file x from fastq_ch.collect()
    script:
    """
    touch \$(realpath $x)
    """
}

process watch {
    publishDir "${params.results}/latest-fastq", mode: 'copy', overwrite: true, pattern: '*.fastq'
    publishDir "${params.results}/latest-stats", mode: 'copy', overwrite: true, pattern: '*.txt'
    
    tag "new reads: ${x}"
    echo true

    input:
        file x from watch_ch.flatten()
    output:
        file '*.fastq' into fastq_latest_ch
        file '*stats.txt'

    script:
    """
    dir=\$(dirname \$(realpath $x))
    barcodename=\$(basename \$dir)

    # to get first and last time stamps, use only these files, sorted by mtime
    # when read in R, they are POSIXct already
    get-times.sh min \$(ls -d \$dir/* | head -n 1) > firsttime.txt
    get-times.sh max \$(ls -d \$dir/* | tail -n 1) > lasttime.txt

    cat \$dir/*.fastq > \$barcodename.fastq

    # make stats file, adding also  min and max time
    seqkit stats -a \$barcodename.fastq | sed '/^file/d' | tr -d ',' | paste -d " " - firsttime.txt lasttime.txt > \$barcodename-stats.txt

    """
}

