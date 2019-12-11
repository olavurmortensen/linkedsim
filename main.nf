#!/usr/bin/env nextflow

params.n_samples = null
params.region = null
params.mutation_rate = null
params.reference = null
params.vcf = null
params.barcodes = null
params.processors = null
params.CF = null
params.CR = null
params.N_FP = null
params.Mu_F = null
params.SR = null
params.Fast_mode = null
params.Seq_error = null
params.Error_rate = null
params.Path_Seq_qual = null
params.Path_Barcode_qual = null
params.Mu_IS = null
params.Std_IS = null
params.Path_barcodepool = null
params.outdir = null

assert params.n_samples != null, 'Input parameter "n_samples" cannot be unasigned.'
assert params.region != null, 'Input parameter "region" cannot be unasigned.'
assert params.mutation_rate != null, 'Input parameter "mutation_rate" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.vcf != null, 'Input parameter "vcf" cannot be unasigned.'
assert params.processors != null, 'Input parameter "processors" cannot be unasigned.'
assert params.CF != null, 'Input parameter "CF" cannot be unasigned.'
assert params.CR != null, 'Input parameter "CR" cannot be unasigned.'
assert params.N_FP != null, 'Input parameter "N_FP" cannot be unasigned.'
assert params.Mu_F != null, 'Input parameter "Mu_F" cannot be unasigned.'
assert params.SR != null, 'Input parameter "SR" cannot be unasigned.'
assert params.Fast_mode != null, 'Input parameter "Fast_mode" cannot be unasigned.'
assert params.Seq_error != null, 'Input parameter "Seq_error" cannot be unasigned.'
assert params.Error_rate != null, 'Input parameter "Error_rate" cannot be unasigned.'
assert params.Path_Seq_qual != null, 'Input parameter "Path_Seq_qual" cannot be unasigned.'
assert params.Path_Barcode_qual != null, 'Input parameter "Path_Barcode_qual" cannot be unasigned.'
assert params.Mu_IS != null, 'Input parameter "Mu_IS" cannot be unasigned.'
assert params.Std_IS != null, 'Input parameter "Std_IS" cannot be unasigned.'
assert params.Path_barcodepool != null, 'Input parameter "Path_barcodepool" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

reference = file(params.reference, checkIfExists: true)
vcf = file(params.vcf, checkIfExists: true)
Path_Seq_qual = file(params.Path_Seq_qual, checkIfExists: true)
Path_Barcode_qual = file(params.Path_Barcode_qual, checkIfExists: true)
Path_barcodepool = file(params.Path_barcodepool, checkIfExists: true)
outdir = file(params.outdir)


process make_small_fasta {
    output:
    file "*.fa" into ref_fa_mutate_ch, ref_fa_hap_ch

    script:
    """
    samtools faidx $reference $params.region > "reference_region.fa"
    """
}

process mutate_fasta {
    input:
    file fasta from ref_fa_mutate_ch
    
    output:
    file "*.fa" into mutated_fa_hap_ch
    
    script:
    """
    mutate_fasta.py $fasta $params.mutation_rate > mutated.fa
    """
}

sample_id_ch = Channel.from(1..params.n_samples)

make_hap_data = sample_id_ch.combine(ref_fa_hap_ch).combine(mutated_fa_hap_ch)

process make_hap {
    input:
    set sample, file(ref_fa), file(mutated_fa) from make_hap_data
    
    output:
    set sample, file("hap1.fa"), file("hap2.fa") into fasta_hap_config_ch, fasta_hap_run_ch
    
    script:
    """
    make_hap.py $ref_fa $mutated_fa > hap1.fa
    make_hap.py $ref_fa $mutated_fa > hap2.fa
    """
}

process make_config {
    input:
    set sample, file(hap1), file(hap2) from fasta_hap_config_ch

    output:
    set sample, file("config"), file(hap1), file(hap2) into data_run_ch
    
    script:
    config_file = "config/sample${sample}.config"
    """
    mkdir config
    echo Path_Fastahap1=$hap1 >> $config_file
    echo Path_Fastahap2=$hap2 >> $config_file
    echo processors=$params.processors >> $config_file
    echo CF=$params.CF >> $config_file
    echo CR=$params.CR >> $config_file
    echo N_FP=$params.N_FP >> $config_file
    echo Mu_F=$params.Mu_F >> $config_file
    echo SR=$params.SR >> $config_file
    echo Fast_mode=$params.Fast_mode >> $config_file
    echo Seq_error=$params.Seq_error >> $config_file
    echo Error_rate=$params.Error_rate >> $config_file
    echo Path_Seq_qual=$Path_Seq_qual >> $config_file
    echo Path_Barcode_qual=$Path_Barcode_qual >> $config_file
    echo Mu_IS=$params.Mu_IS >> $config_file
    echo Std_IS=$params.Std_IS >> $config_file
    echo Path_barcodepool=$Path_barcodepool >> $config_file
    echo Hap=2 >> $config_file
    """
}

process run_sim {
    publishDir "$outdir/sample$sample", mode: 'copy'

    input:
    set sample, file("config"), file(hap1), file(hap2) from data_run_ch
    
    output:
    file "*fastq.gz"
    
    script:
    """
    LRTK-SIM.py $config
    mv config/*/*fastq.gz .
    """
}












