rule detect_differential_isoforms:
    input:
        group_1='data/splicing_order_stringent_version/{sample}_group_2.txt',
        group_2='data/splicing_order_stringent_version/{sample}_group_2.txt',
        abundance='data/splicing_order_stringent_version/replicates_separate_esp_format/{sample}_GM19209_interm_counts_per_allele.hac.all_introns.min10reads.filterND.stringent.3_isoforms.esp.txt',
    output:
        DTU='results/detect_differential_isoforms_interm_counts_per_allele/{sample}/differential_transcripts.tsv',
    log:
        'logs/detect_differential_isoforms_espresso_filtered/{sample}/differential_transcripts.log',
    params:
        out_dir='results/detect_differential_isoforms_espresso/{sample}',
        adj_pvalue=config['adj_pvalue'],
        delta_proportion=config['delta_proportion'],
    conda:
        '../envs/rmats_long.yaml'
    threads: 1
    shell:
        'python workflow/resources/rMATS-long/scripts/detect_differential_isoforms.py'
        ' --abundance {input.abundance}'
        ' --out-dir {params.out_dir}'
        ' --group-1 {input.group_1}'
        ' --group-2 {input.group_2}'
        ' --num-threads {threads}'
        ' --adj-pvalue {params.adj_pvalue}'
        ' --delta-proportion {params.delta_proportion}'
        ' &> {log}'


rule classify_isoform_differences:
    input:
        gtf='workflow/resources/references/ensembl.Homo_sapiens.GRCh38.86.gtf',
        updated_gtf='data/espresso/LCL_espresso_updated.gtf',
    output:
        isoform_differences='results/classify_isoform_differences/isoform_differences_{main_tx_id}_to_{second_tx_id}.tsv'
    log:
        'logs/classify_isoform_differences/{main_tx_id}_to_{second_tx_id}.log'
    conda:
        '../envs/rmats_long.yaml'
    shell:
        'python workflow/resources/rMATS-long/scripts/classify_isoform_differences.py'
        ' --main-transcript-id {wildcards.main_tx_id}'
        ' --second-transcript-id {wildcards.second_tx_id}'
        ' --gencode-gt {input.gtf}'
        ' --updated-gtf {input.updated_gtf}'
        ' --out-tsv {output.isoform_differences}'
