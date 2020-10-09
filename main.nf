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
    .set { fastq_ch }

process touch {
    input:
        file x from fastq_ch.collect()
    script:
    """
    touch \$(realpath $x)
    """
}

/*
This is the key to the whole pipeline - 
*/
// 
Channel
    .watchPath("${params.fastq_pass}/**.fastq", 'create,modify')
    .map { file ->
        //def key = file.name.toString().tokenize('_').get(2)
        def key = file.getParent()
        return key
     }
    .distinct()
    //.view()
    .set { watch_fastq_pass }

process watch {
    //nothing to publish, this ends here
    publishDir "${params.results}/latest-stats", mode: 'copy', overwrite: true, pattern: '*stats.txt'
    tag "new reads detected: ${x}"

    // flatten is used as it emits each file as a single item (does not wait), try collate here?
    input:
        path x from watch_fastq_pass
    output:
        file '*stats.txt'

    script:
    """
    # this is executed for each new bunch of reads, leads to blowing up the storage! 
    # so cat directly to latestfastqdir

    cat ${x}/*.fastq > $latestfastqdir/${x}.fastq

    # do this in the same process, otherwise comes to file blocking from cat and clashes
    head $latestfastqdir/${x}.fastq | get-times.sh min > firsttime.txt
    tail $latestfastqdir/${x}.fastq | get-times.sh max > lasttime.txt

    seqkit stats -a $latestfastqdir/${x}.fastq | \
    sed '/^file/d' | tr -d ',' | \
    paste -d " " - firsttime.txt lasttime.txt > ${x}-stats.txt


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

Channel
    .watchPath("${params.fastq_pass}/**.fastq", 'create,modify')
    .map { file ->
        def key = file.getParent().getFileName()//get parent folder name, i.e. barcode01
        def filename = file.simpleName
        return tuple(key, filename, file)
     }
    .set { watch_fastq_pass_2 }

Channel
    .watchPath("${latestfastqdir}/*.fastq", 'create,modify')
    .map { file -> tuple(file.simpleName, file) }
    .set { merged_fastq_ch }

// this is now barcode, filename, file, krakendb
kraken_ch_2 = watch_fastq_pass_2.combine(kraken_db_ch)
//kraken_ch_2.view()

process kraken {
    container 'aangeloo/kraken2:latest'
    tag "working on: ${filename} from ${key}"
    publishDir "${params.results}/scratch-kraken", mode: 'copy', pattern: '*.tsv'
    
    input:
        tuple key, filename, file(fastq), file(db) from kraken_ch_2
    
    output:
        //file("*report") // both kraken2 and the bracken-corrected reports are published and later used in pavian?
        file("*bracken.tsv")
    
    script:
    def memory = params.weakmem ? "--memory-mapping" : ""  // use --memory-mapping to avoid loading db in ram on weak systems
    def rlength = 250
    
        """
        kraken2 \
            -db $db \
            $memory \
            --report ${filename}_k2.report \
            ${fastq} \
            > kraken2.output

        bracken \
            -d $db \
            -r $rlength \
            -i ${filename}_k2.report \
            -l ${params.taxlevel} \
            -o ${filename}_bracken.tsv
        """

}
