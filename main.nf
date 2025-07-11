#!/usr/bin/env nextflow


params.in = "/mnt/c/Users/Vanessa Fleischer/NextFlow/NextFlow_ABI/exam/hepatitis/"
params.se_glob = "*.{fa,fasta}"      // For all FastQ, we'll filter later for true single-end
params.accession = "M21012"
params.out = 'results/'

//------------- processes ---------------//

process download_reference {

    conda 'bioconda::entrez-direct=24.0'
    
    input: 
        val accession
    
    output:
        path ref_file

    script:
        ref_file = "${accession}.fa"

        """
        esearch -db nucleotide -query "$accession"\
        | efetch -format fasta \
            > "$ref_file"
        """
}

process mergeFasta {

    input:
    path fasta_files
    path ref_ch

    output:
    path 'merged.fasta'

    script:
    """
    cat ${ref_ch} ${fasta_files.join(' ')} > merged.fasta
    """
}

process mafft_aligner {

    conda 'bioconda::mafft=7.525'

    input:
        path merged_fasta


    output:
        path "aligned.fasta"

    script:
        """
        mafft "$merged_fasta" > aligned.fasta
        """

}

process trimal {
    publishDir 'results/cleaned_alignment', mode: 'copy'

    conda 'bioconda::trimal=1.5.0'

    input:
        path aligned_fasta

    output:
        path "${aligned_fasta.simpleName}_trimmed.fasta"
        path "${aligned_fasta.simpleName}_report.html"

    script:
        """
        trimal -in $aligned_fasta \
               -out ${aligned_fasta.simpleName}_trimmed.fasta \
               -htmlout ${aligned_fasta.simpleName}_report.html \
               -automated1
        """
}


workflow {
    def ref_ch = download_reference(params.accession)

    def fasta_files = channel
    .fromPath(
        "$params.in/$params.se_glob")
    .collect()

    mergeFasta (fasta_files, ref_ch)
    .set { merged_fasta }

    mafft_aligner(merged_fasta)
        .set { aligned_fasta }

    trimal(aligned_fasta)
}
