/*
 * Copyright (c) 2013-2019, Centre for Genomic Regulation (CRG).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 *
 */


/*
 * Proof of concept of a RNAseq pipeline implemented with Nextflow
 *
 * Authors:
 * - Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * - Emilio Palumbo <emiliopalumbo@gmail.com>
 * - Evan Floden <evanfloden@gmail.com>
 */


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"

transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)

log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 transcriptome: ${params.transcriptome}
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch; read_pairs2_ch }

Channel
    .fromPath(params.transcriptome)
    .ifEmpty { exit 1, "Transcript fasta file is unreachable: ${params.transcriptome}"  }
    .set { tx_fasta_ch  }

process index {
    tag "$transcriptome.simpleName"

    input:
    file transcriptome from tx_fasta_ch

    output:
    file 'index' into index_ch

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


process quant {
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy'

    input:
    file index from index_ch.collect()
    set pair_id, file(reads) from read_pairs_ch

    output:
    file(pair_id) into quant_ch

    script:
    """
    salmon quant --validateMappings --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

process fastqc {
    tag "FASTQC on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    set sample_id, file(reads) from read_pairs2_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}


process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file('*') from quant_ch.mix(fastqc_ch).collect()

    output:
    file('multiqc_report.html')

    script:
    """
    multiqc .
    """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
