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
 params.kraken_db = "" // this is a db FOLDER prepared manually, gives more flexibiliy for custom db and simplifies things
 params.weakmem = false
 params.skip_kraken = false
 params.taxlevel = "S" //level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
 params.help = false


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
         --kraken_db        : path to kraken2 database (the kraken db folder)
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
    //
    publishDir "${params.results}/latest-stats", mode: 'copy', overwrite: true, pattern: '*stats.txt'
    tag "new reads detected: ${x}"
    containerOptions "--volume '${latestfastqdir}:${latestfastqdir}'"

    // flatten is used as it emits each file as a single item (does not wait), try collate here?
    input:
        path x from watch_fastq_pass
    output:
        file '*stats.txt'

    script:
    """
    # this is executed for each new bunch of reads, leads to blowing up the storage! 
    # so cat directly to latestfastqdir
    # BUT latestfastqdir has to be mounted in docker explicitly!

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


// kraken db is a ready to use directory -

if( !params.skip_kraken) {
Channel
    .fromPath("${params.kraken_db}", type: 'dir')
    .ifEmpty { error "Can not find ${params.kraken_db}" }
    .set { db_ch_dir }
} else {
    db_ch_dir = Channel.empty()
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
kraken_ch = watch_fastq_pass_2.combine(db_ch_dir)

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
