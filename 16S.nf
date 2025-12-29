#!/bin/env nextflow
nextflow.enable.dsl=1

/*
 * Releases
 */

// user: heweihuang
// version: v1.0.1  whhe(2025-12-23) 16S pipeline
/*
 * print usage
 */
params.h = false
params.help = false
if(params.help || params.h){
log.info ''
log.info '16S.nf'
log.info '=================================================================================================================================='
log.info 'methlylation sample for rawdata fastq.gz'
log.info 'Usage:'
log.info '  nextflow run BS_panel.nf --i /rawdata input path [*fq.gz]/ --t cpu_number --o /outut path/ '
log.info ''
log.info 'Options:'
log.info '      --help                      Show this message and exit.'
log.info '      --i                 <str>   input (fastq.gz/fq.gz) path                              [./]'
log.info '      --id                <str>   the sample name                                          [out]'
log.info '      --o                 <str>   the output dir                                           [./]'
log.info '      --groups_type       <str>   groups type'
log.info '  QC_options:'
log.info '      --qc                <str>   Whether the rawdata is qc filtered                       [true]'
log.info '      --qc_cpu            <int>   cpu number for each QC job                               [4]'
log.info '      --adapter           <str>   adapter sequences file                                   [/fuer2/01.software/2.fastp/adapters/TruSeq3-PE-2.fa]'
log.info '      --clean_rate        <int>   threshold of clean data rate                             [80]'
log.info '  Map_options:'
log.info '      --m                 <str>   Whether the reads is mapping the reference database      [true]'
log.info '      --map_cpu           <int>   cpu number for each rm rRNA job                          [8]'
log.info '      --mapping_rate      <int>   threshold of rRNA mapping rate                           [80]'
log.info '      --target_bed        <str>   Capture range of the panel                               [/fuer2/01.Pipeline/03.pipeline_test/01.ST_DNA/db/bed/regioncount.bed]'
log.info '      --expan_bed         <str>   panel bed expan range                                    [/fuer2/01.Pipeline/03.pipeline_test/01.ST_DNA/db/bed/regioncount.expan.bed]'
exit 1
}

outdirAbs              = ""

if( params.d[0] == "." ){
    outdirAbs = "${launchDir}/${params.d}"
}else{
    outdirAbs = "${params.d}"
}

Channel
    .fromFilePairs("$params.i/*{_1,_2,_R1,_R2}*{fq,fastq}*",size:2)
    .ifEmpty { error "Cannot find fq/fq.gz/fastq/fastq.gz files: $params.i" }
    .set{raw_reads}

process merge{
    tag "${libName}"
    cpus "${params.qc_cpus}"
   
    input:
    set val(libNameRaw), file(r_reads) from raw_reads
    output:
    //set val(libName), file("*.join.fastq") into join_fq
    // 输出具体的文件名，而不是通配符
    file("${libName}.join.fastq") into join_fq
    script:
       libName = libNameRaw.split(/[_]/)[0]

   """
   ${params.join_paired} -m fastq-join -j 8 -p 10 -f ${r_reads[0]} -r ${r_reads[1]} -o joined/${libName}
   cp  joined/${libName}/fastqjoin.join.fastq ${libName}.join.fastq
   """
   }

process QC{
    tag "merge_QC"
    cpus "${params.qc_cpus}"

    input:
    file(join_files) from join_fq.collect()
    output:
    file("filter.fa") into clean_fq
    script:
    def list = join_files.name.collect { it.replace('.join.fastq', '') }.join(',')
    def list_fq = join_files.name.collect{ it }.join(',')  
   """
   (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.meger_split_fq}"' -i '"${list_fq}"' --sample_id '"${list}"' -o ./ -m '"${params.map_txt}"' -q 19 --barcode_type '"'not-barcoded'"' --phred_offset=33')
   ${params.vsearch} --derep_fulllength seqs.fna --sizeout --minuniquesize 8 --output qc_tmp.1.fa
   ${params.vsearch} --uchime_ref qc_tmp.1.fa -db ${params.fa}  --nonchimeras qc_tmp.2.fa
   ${params.vsearch} --fastx_filter qc_tmp.2.fa -fastaout filter.fa  --fastq_minlen 200 --fastq_maxlen 1000 --threads ${params.qc_cpus}
   """
   }

 
process Map_OTU{
    tag "Map_OTU"
    cpus "${params.map_cpus}"
    input:
    file(filter_fa) from clean_fq
    output:
    file("tax_table.txt") into otu_out
    file("result/otu_table_mc2_w_tax_no_pynast_failures.biom") into biom_out1, biom_out2, biom_out3
    file("result/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta") into map_fa

    """
    rm -fr result
    export PATH="/data_center_01/home/heweihuang/miniconda3/envs/qiime/bin:\$PATH"
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.pick_otu}"' -o result -i filter.fa -r '"${params.fa}"'  -p '"${params.params_txt}"'')
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.biom}"' convert -i result/otu_table_mc2_w_tax_no_pynast_failures.biom -o tax_table.txt --to-tsv --header-key taxonomy --output-metadata-id '"'Consensus Lineage'"'')
    """
}


process diversity{
    tag "alpha-beta"
    cpus "${params.diver_cpus}"
    errorStrategy='ignore'
    input:
    file(biom) from biom_out1
    output:
    file("result/bdiv_even*/unweighted_unifrac_dm.txt") into unifrac_dm

    """
    export PATH="/data_center_01/home/heweihuang/miniconda3/envs/qiime/bin:\$PATH"    
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.biom}"' summarize-table -i otu_table_mc2_w_tax_no_pynast_failures.biom > sample.depth')
    min_depth=`awk 'NR>15' sample.depth  | awk '{print \$2}' |sort |awk 'NR==1{print int(\$0*0.85)}'`
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.diversity}"' -i otu_table_mc2_w_tax_no_pynast_failures.biom -o result -m '"${params.map_txt}"' -t '"${params.ref_tre}"' -e '"\${min_depth}"' -p '"${params.diver_par}"'') 
    """
 }


process tree{
    tag "tree"
    cpus "${params.tree_cpus}"
    errorStrategy='ignore'
    input:
    file(map_file) from map_fa
    file(biom_file) from biom_out2
    output:
    file("otu_table_0.001.fraction.txt") into tree_txt
    file("rep_set_0.001_fraction.tre") into tree_tre

    """
    export PATH="/data_center_01/home/heweihuang/miniconda3/envs/qiime/bin:\$PATH"
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.filter_otus}"' -i otu_table_mc2_w_tax_no_pynast_failures.biom --min_count_fraction 0.001 -o otu_table_0.001.fraction.biom')
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.filter_fa}"' -f rep_set_aligned_pfiltered.fasta -b otu_table_0.001.fraction.biom -o rep_set_0.001_fraction.fa')
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.make_phylogeny}"' -i rep_set_0.001_fraction.fa -o rep_set_0.001_fraction.tre')
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.biom}"' convert -i otu_table_0.001.fraction.biom -o otu_table_0.001.fraction.txt --to-tsv --header-key taxonomy')
    """
 }

process diff{
    tag "diff"
    cpus "${params.diff_cpus}"
    errorStrategy='ignore'
    input:
    file(biom2) from biom_out3
    output:
    file("diff_otus.txt") into diff_out

    """
    export PATH="/data_center_01/home/heweihuang/miniconda3/envs/qiime/bin:\$PATH"
    (bash -c 'for var in \$(env | grep -iE "SLURM|PMI|MPI|OMPI" | cut -d= -f1); do unset \$var; done; exec '"${params.diff_abun}"' -i otu_table_mc2_w_tax_no_pynast_failures.biom -o diff_otus.txt -m '"${params.map_txt}"' -c SampleType -x normal -y tumor -a DESeq2_nbinom')
    """
 }

process plot{
    tag "plot"
    cpus "${params.plot_cpus}"
    errorStrategy='ignore'
    input:
    file(f1) from otu_out
    file(f2) from tree_txt
    file(f3) from tree_tre
    file(f4) from unifrac_dm
    output:
    file("*.png") into plot_out

    """
    list=`awk 'NR>1{print \$1}' ${params.map_txt}|xargs|sed 's/ /,/g'`
    Rscript ${params.venn} tax_table.txt venn.pdf \${list}
    Rscript ${params.heatmap} tax_table.txt heatmap.pdf 100 
    python ${params.otu_phy} -i tax_table.txt -o ./
    Rscript ${params.ggtree_heatmap} otu_table_0.001.fraction.txt rep_set_0.001_fraction.tre circos_output.pdf heatmap_output.pdf
    awk -F '\\t' 'BEGIN{OFS="\\t"}{if(NR==1){if(\$1==""){\$1="Sample"}print \$0}else{print \$0}}' unweighted_unifrac_dm.txt  > unweighted_unifrac_dm.plot.txt
    python ${params.upgma} -i unweighted_unifrac_dm.plot.txt -o my_cluster.png
    """
 }

workflow.onComplete {
    if( workflow.success ) {
        File f = new File("./pipe.Done")
        f.write("Pipeline completed at: $workflow.complete"+"\n") 
        f.append("Execution status:      ${ workflow.success ? 'OK' : 'failed' }" + "\n")
        f.append("User:                  $workflow.userName" + "\n")
        f.append("Launch time:           ${workflow.start.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Ending time:           ${workflow.complete.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Duration:              $workflow.duration" + "\n")
        f.append("Total CPU-Hours:       ${workflow.stats.computeTimeFmt ?: '-'}" + "\n")
        f.append("Tasks stats:           Succeeded ${workflow.stats.succeedCountFmt}; Cached ${workflow.stats.cachedCountFmt}; Ignored ${workflow.stats.ignoredCountFmt}; Failed ${workflow.stats.failedCountFmt}"  + "\n")
    }else{
        File f = new File("./pipe.Failed")
        f.write("Pipeline failed at: $workflow.complete"+"\n")
        f.append("Execution status:   ${ workflow.success ? 'OK' : 'failed' }" + "\n")
        f.append("User:               $workflow.userName" + "\n")
        f.append("Launch time:        ${workflow.start.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Ending time:        ${workflow.complete.format('yyyy-MMM-dd HH:mm:ss')}" + "\n")
        f.append("Duration:           $workflow.duration" + "\n")
        f.append("Total CPU-Hours:    ${workflow.stats.computeTimeFmt ?: '-'}" + "\n")
        f.append("Tasks stats:        Succeeded ${workflow.stats.succeedCountFmt}; Cached ${workflow.stats.cachedCountFmt}; Ignored ${workflow.stats.ignoredCountFmt}; Failed ${workflow.stats.failedCountFmt}" + "\n")
        f.append("ERROR message:\n" + "  ${workflow.errorMessage}" + "\n")
        f.append("ERROR report:\n" + "  ${workflow.errorReport}" + "\n")
    }
}
