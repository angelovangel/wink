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
 params.kraken_gz = "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz"
 params.kraken_store = "$HOME/db/kraken"
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
         --kraken_gz        : path to kraken2 database, as tgz file (ftp or absolute path)
         --kraken_store     : path to permanently store the kraken2 database, will be used in subsequent runs
         --taxlevel         : taxonomic level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
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
         --kraken_gz        : ${params.kraken_gz}
         --kraken_store     : ${params.kraken_store}
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

// dir to store and watch merged fastq
 latestfastqdir = file("${workflow.launchDir}/fastq-latest")
 if( !latestfastqdir.exists() ) {
 latestfastqdir.mkdir()
 }

// used to touch the files at start so that all available fastq files are processed 
 Channel
    .fromPath( "${params.fastq_pass}/**.fastq", type: 'file')
    //.ifEmpty { error "Can not find any barcode* folders in ${params.fastq_pass}" }
    //.view()
    .set { fastq_ch }


Channel
    .watchPath("${params.fastq_pass}/**.fastq", 'create,modify')
    //.collate( 10 ) // how often will the process execute
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
    //nothing to publish, this ends here
    tag "new reads detected: ${x}"

    // flatten is used as it emits each file as a single item (does not wait), try collate here?
    input:
        file x from watch_ch.flatten()

    script:
    """
    dir=\$(dirname \$(realpath $x))
    barcodename=\$(basename \$dir)

    # this is executed for each new bunch of reads, leads to blowing up the storage! 
    # so cat directly to latestfastqdir
    cat \$dir/*.fastq > $latestfastqdir/\$barcodename.fastq

    """
}

// use collate to determine how often kraken will be run - 
// 
Channel
    .watchPath("${latestfastqdir}/*.fastq", 'create,modify')
    .map { file -> tuple(file.simpleName, file) }
    .into { merged_fastq_ch1; merged_fastq_ch2 }

process seqkit {
    publishDir "${params.results}/latest-stats", mode: 'copy', overwrite: true, pattern: '*stats.txt'
    tag "working on: ${x}"
    //echo true

    input:
        tuple filename, file(x) from merged_fastq_ch1
    output:
        file '*stats.txt'

    script:
    """
    # get first and last time stamps from fastq headers
    # when read in R, they are POSIXct already

    head $x | get-times.sh min > firsttime.txt
    tail $x | get-times.sh max > lasttime.txt

    seqkit stats -a ${x} | sed '/^file/d' | tr -d ',' | paste -d " " - firsttime.txt lasttime.txt > ${filename}-stats.txt

    """ 
}

if(params.kraken_gz){
    Channel
        .of( "${params.kraken_gz}" )
        .set { kraken_gz_ch }
} else {
        kraken_gz_ch = Channel.empty()
}

process krakenDB {
    storeDir "${params.kraken_store}"

input:
    path kraken_file from kraken_gz_ch

output:
    path "**", type: 'dir' into kraken_db_ch

script:
"""
tar -xf $kraken_file
"""
}

kraken_channel = merged_fastq_ch2.combine(kraken_db_ch) //list with 3 elements

process kraken2 {
    container 'aangeloo/kraken2:latest'
    tag "working on: ${sample_id}"
    //echo true
    publishDir "${params.results}/samples", mode: 'copy', overwrite: true, pattern: '*.{report,tsv}'
    
    input:
        //path db from kraken_db_ch
        tuple sample_id, file(x), file(y) from kraken_channel
    
    output:
        file("*report") into kraken2mqc_ch // both kraken2 and the bracken-corrected reports are published and later used in pavian?
        tuple sample_id, file("*bracken.tsv") into bracken2dt_ch
        file("*bracken.tsv") into bracken2summary_ch
    
    script:
    //def single = x instanceof Path
    //def kraken_input = single ? "\"${ x }\"" : "--paired \"${ x[0] }\"  \"${ x[1] }\""
    def memory = params.weakmem ? "--memory-mapping" : ""  // use --memory-mapping to avoid loading db in ram on weak systems
    def rlength = 250
    
        """
        kraken2 \
            -db $y \
            $memory \
            --report ${sample_id}_kraken2.report \
            ${x} \
            > kraken2.output

        bracken \
            -d $y \
            -r $rlength \
            -i ${sample_id}_kraken2.report \
            -l ${params.taxlevel} \
            -o ${sample_id}_bracken.tsv
        """

}