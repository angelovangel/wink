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
 params.kraken_db = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20200919.tar.gz"
 params.kraken_store = "$HOME/db/kraken"
 params.weakmem = false
 params.skip_kraken = false
 params.taxlevel = "S" //level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
 params.help = false

/*
assign store dir dynamically
*/
kraken_dbname = file("${params.kraken_db}").getSimpleName()
// even this method exists! getSimpleName()
curr_kraken_store = "${params.kraken_store}/${kraken_dbname}"
println("Will use ${curr_kraken_store} as kraken_store dir")

/*
 dir to store and watch merged fastq and scratch for intermediate bracken results
*/

 latestfastqdir = file("${params.results}/latest-fastq")
 if( !latestfastqdir.exists() ) {
 latestfastqdir.mkdirs() // attention - mkdir() vs mkdirs(), the latter creates parents if they do not exist
 }

 scratchdir = file("${params.results}/scratch")
 if( !scratchdir.exists() ) {
 scratchdir.mkdirs()
 }

def helpMessage() {
log.info """
        ===========================================
         W I N K - What's In my Nanopore reads, with Kraken2

         real-time monitoring of a Nanopore run + taxonomic and abundance information

         Parameters:
        -------------------------------------------
         --fastq_pass       : the folder where basecalled reads are saved during a run, must contain barcodes
         --results          : where the results will go
         --kraken_db        : path to kraken2 database, as tgz file (ftp or absolute path), no quotes
         --kraken_store     : path to permanently store the kraken2 database, will be used in subsequent runs
         --taxlevel         : taxonomic level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
         --weakmem          : do not load database in memory, use on weak machines
         --skip_kraken      : skip kraken2 classification, just run statistics
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
         --kraken_store     : ${params.kraken_store}
         --taxlevel         : ${params.taxlevel}

         Runtime data:
        -------------------------------------------
         Running with profile:   ${workflow.profile}
         Container:              ${workflow.container}
         Running as user:        ${workflow.userName}
         Launch dir:             ${workflow.launchDir}
         Base dir:               ${baseDir}
         kraken db store dir     ${curr_kraken_store}
         """
         .stripIndent()




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

//
//This is the key to the whole pipeline - 

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

    # make sure the files are cat in the right order
    find -L ${x} -type f -name '*.fastq' | sort -V | xargs cat > $latestfastqdir/${x}.fastq


    # do this in the same process, otherwise comes to file blocking from cat and clashes
    head $latestfastqdir/${x}.fastq | get-times.sh min > firsttime.txt
    tail $latestfastqdir/${x}.fastq | get-times.sh max > lasttime.txt

    seqkit stats -a $latestfastqdir/${x}.fastq | \
    sed '/^file/d' | tr -d ',' | \
    paste -d " " - firsttime.txt lasttime.txt > ${x}-stats.txt


    """
}

// some work to enable kraken db to be one of:
// 1. http:// --> does not exist as a file on the file system
// 2. downloaded tgz file --> exists and isFile() returns true (method for file obj)
// 3. a ready to use directory --> exists and isDirectory() returns true (method for file obj)
// 1 and 2 are covered by the db_ch_file

// set a var to track if an existing dir was provided. In case it is, emit empty kraken_gz_ch (skip DBPrep)
// and emit the directory to kraken_db_ch
doKrakenDBPrep = true
File kraken_dir_exists = new File("${params.kraken_db}") 

if( kraken_dir_exists.isDirectory() ) {
    println "will use existing kraken db directory"
    doKrakenDBPrep = false
    Channel
        .fromPath("${params.kraken_db}", type: 'dir')
        .set { db_ch_dir }
}

if(!params.skip_kraken){
    Channel
        .fromPath( "${params.kraken_db}" )
        .set { db_ch_file }
} else {
        db_ch_file = Channel.empty()
}

process krakenDBPrep {
    storeDir "${curr_kraken_store}"

    when:
        !params.skip_kraken
        doKrakenDBPrep
    input:
        path kraken_file from db_ch_file
    output:
        path "**", type: 'dir' into kraken_db_ch

    script:
    dbname = kraken_file.baseName
    """
    mkdir -p $dbname && tar -xzf $kraken_file -C $dbname
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

// this is now barcode, filename, file, krakendb
// kraken_db_path solves the issue where the tar archive may be with or without a leading directory
kraken_db_path = doKrakenDBPrep ? kraken_db_ch.flatten().last() : db_ch_dir
kraken_ch = watch_fastq_pass_2.combine(kraken_db_path)
//kraken_ch_2.view()

//exit 1

process kraken {
    //container 'aangeloo/kraken2:latest'
    tag "working on: ${filename}; barcode: ${key}"
    publishDir "${scratchdir}/${key}", mode: 'copy', pattern: '*.tsv'
    
    when:
        !params.skip_kraken
    input:
        tuple key, filename, file(fastq), file(db) from kraken_ch
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

        # get # of unclassified reads from _k2.report and add to _bracken.tsv, they are under 'kraken_assigned_reads', name: unclassified
        
        unclassified=\$(head -n 1 ${filename}_k2.report | cut -f 2)
        echo -e unclassified '\t\t\t' \$unclassified '\t' 0 '\t' \$unclassified '\t' >> ${filename}_bracken.tsv
        """

}


Channel
    .watchPath("${scratchdir}/*", 'create,modify') // simple, just watch for directory change, mtime updates when new files are written
    .map { n -> 
            def barcode = n.getFileName()
            def barcodedir = n
            return tuple(barcode, barcodedir)
     }
     //.unique() // can not use unique() here, because it will execute once per barcode and then stop
     // distinct() works, but dangerous (if only one sample!), what about interleave with a fake?
    .set { combine_ch } 
    //combine_ch.view()

process combine {
    tag "working on: ${barcode}"
    publishDir "${params.results}/latest-bracken", mode: 'copy', pattern: '*tsv'
    
    when:
        !params.skip_kraken
    input:
        tuple barcode, path(barcodedir) from combine_ch
    output:
        file '*.tsv'
    
    script:
    
    """
    bracken-combine.R ${barcodedir} ${barcode}-bracken.tsv
    """

}
