


ref_fasta = file("s3://uwlm-chemistry-data/reference/swiss.HUMAN.20141123.amyloid.v2.fasta", checkIfExists: true)
publish_path = 's3://uwlm-chemistry-data/crux_output/'

raw_file_ch = Channel.fromPath( "${params.input_path}/*.raw" )

process msconvert {
    cpus 2
    memory "2GB"
    container 'docker.labmed.uw.edu/docker-images/msconvert:045dd857'
    label 'msconvert'

    input:
      file(input) from raw_file_ch
    output:
      file("*.ms2") into ms2_file_ch

    script:
    """
    WINEDEBUG=-all
    wine msconvert \
      --ms2 --inten64 --zlib \
      --filter "peakPicking true 1-2" \
      --filter "msLevel 1-2" \
      ${input}
    """

    
}


process crux {
    cpus 8
    memory "8GB"
    container 'docker.labmed.uw.edu/docker-images/crux:91681dfd'
    label 'crux'
    publishDir publish_path, overwrite: true

    input:
        path ref_fasta
        path("inputs/*") from ms2_file_ch
    output:
        path("inputs/")
    
    script:
    """
    cruxpipeline.py --search-engine comet ${ref_fasta} inputs/
    """
}
