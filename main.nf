#!/usr/bin/env nextflow
/*
 * Authors:
 * - Daniel Cook <danielecook@gmail.com>
 * - Ye Wang <yewangfaith@gmail.com>
 */

/*
    Globals
*/

// Define contigs here!
CONTIG_LIST = ["I", "II", "III", "IV", "V", "X", "MtDNA"]
contig_list = Channel.from(CONTIG_LIST)

/* 
    ======
    Params
    ======
*/

date = new Date().format( 'yyyyMMdd' )
params.out = "concordance-${date}"
params.debug = false
params.tmpdir = "tmp/"
params.email = ""
params.alignments = "bam"

reference_handle = file("${params.reference}/*")

fq_concordance_script = file("fq_concordance.R")

// Debug
if (params.debug == true) {
    println """

        *** Using debug mode ***

    """
    params.bamdir = "${params.out}/bam"
    params.fq_file = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.fq_file_prefix = "${workflow.projectDir}/test_data"

} else {
    // The SM sheet that is used is located in the root of the git repo
    params.bamdir = "(required)"
    params.fq_file = "SM_sample_sheet.tsv"
    params.fq_file_prefix = null;
}

File fq_file = new File(params.fq_file);

/* 
    =======================
    Filtering configuration
    =======================
*/

min_depth=3
qual=30
mq=40
dv_dp=0.5

/* 
    ==
    UX
    ==
*/

param_summary = '''

┌─┐┌─┐┌┐┌┌─┐┌─┐┬─┐┌┬┐┌─┐┌┐┌┌─┐┌─┐  ┌┐┌┌─┐
│  │ │││││  │ │├┬┘ ││├─┤││││  ├┤───│││├┤ 
└─┘└─┘┘└┘└─┘└─┘┴└──┴┘┴ ┴┘└┘└─┘└─┘  ┘└┘└  
                                                         
''' + """

    parameters              description                    Set/Default
    ==========              ===========                    =======

    --debug                 Set to 'true' to test          ${params.debug}
    --cores                 Regular job cores              ${params.cores}
    --out                   Directory to output results    ${params.out}
    --fq_file               fastq file (see help)          ${params.fq_file}
    --fq_file_prefix        fastq file (see help)          ${params.fq_file_prefix}
    --reference             Reference Genome               ${params.reference}
    --bamdir                Location for bams              ${params.bamdir}
    --tmpdir                A temporary directory          ${params.tmpdir}
    --email                 Email to be sent results       ${params.email}

    HELP: http://andersenlab.org/dry-guide/pipeline-concordance/

"""


/*
    Fetch fastq files and additional information.
*/
if (params.fq_file_prefix) {
println "Using fq prefix"
println params.fq_file_prefix
//fq_file_prefix = fq_file.getParentFile().getAbsolutePath();
fqs = Channel.from(fq_file.collect { it.tokenize( '\t' ) })
             .map { SM, ID, LB, fq1, fq2, seq_folder -> ["${SM}", ID, LB, file("${params.fq_file_prefix}/${fq1}"), file("${params.fq_file_prefix}/${fq2}"), seq_folder] }
             .view()

} else {
fqs = Channel.from(fq_file.collect { it.tokenize( '\t' ) })
         .map { SM, ID, LB, fq1, fq2, seq_folder -> [SM, ID, LB, file("${fq1}"), file("${fq2}"), seq_folder] }
}

Channel.fromPath("${params.reference}/*")
       .into { reference1; reference2; reference3 }

/* 
    =========
    Alignment
    =========
*/

process perform_alignment {

    cpus params.cores

    tag { ID }

    input:
        set SM, ID, LB, file(fq1), file(fq2), seq_folder from fqs
        file("${params.genome}/*") from reference1.collect()
    output:
        set val(ID), file("${ID}.bam"), file("${ID}.bam.bai") into fq_bam_set
        set val(SM), file("${ID}.bam"), file("${ID}.bam.bai") into SM_aligned_bams

    
    """
        bwa mem -t ${task.cpus} -R '@RG\\tID:${ID}\\tLB:${LB}\\tSM:${SM}' ${params.genome}/${params.genome}.fa.gz $fq1 $fq2 | \\
        sambamba view --nthreads=${task.cpus} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${task.cpus} --show-progress --tmpdir=${params.tmpdir} --out=${ID}.bam /dev/stdin
        sambamba index --nthreads=${task.cpus} ${ID}.bam

        if [[ ! \$(samtools view ${ID}.bam | head -n 10) ]]; then
            exit 1;
        fi
    """
}

fq_bam_set.into { fq_cov_bam; fq_stats_bam; fq_idx_stats_bam }

process coverage_fq {

    tag { ID }

    label 'bam_coverage'

    input:
        set val(ID), file("${ID}.bam"), file("${ID}.bam.bai") from fq_cov_bam
    output:
        file("${ID}.coverage.tsv") into fq_coverage


    """
        bam coverage ${ID}.bam > ${ID}.coverage.tsv
    """
}

process coverage_fq_merge {

    publishDir "${params.out}/fq", mode: 'copy', overwrite: true

    input:
        file fq_set from fq_coverage.toSortedList()

    output:
        file("fq_coverage.full.tsv")
        file("fq_coverage.tsv")

    """
        echo -e 'fq\\tcontig\\tstart\\tend\\tproperty\\tvalue' > fq_coverage.full.tsv
        cat ${fq_set} >> fq_coverage.full.tsv

        cat <(echo -e 'fq\\tcoverage') <( cat fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > fq_coverage.tsv
    """
}

process fq_idx_stats {
    
    tag { ID }

    input:
        set val(ID), file("${ID}.bam"), file("${ID}.bam.bai") from fq_idx_stats_bam
    output:
        file fq_idxstats into fq_idxstats_set

    """
        samtools idxstats ${ID}.bam | awk '{ print "${ID}\\t" \$0 }' > fq_idxstats
    """
}

process fq_combine_idx_stats {

    publishDir "${params.out}/fq", mode: 'copy', overwrite: true

    input:
        file("?.stat.txt") from fq_idxstats_set.toSortedList()

    output:
        file("fq_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > fq_bam_idxstats.tsv
        cat *.stat.txt >> fq_bam_idxstats.tsv
    """

}

process fq_bam_stats {

    tag { ID }

    input:
        set val(ID), file("${ID}.bam"), file("${ID}.bam.bai") from fq_stats_bam

    output:
        file 'bam_stat' into fq_bam_stat_files

    """
        cat <(samtools stats ${ID}.bam | grep ^SN | cut -f 2- | awk '{ print "${ID}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_fq_bam_stats {

    publishDir "${params.out}/fq", mode: 'copy', overwrite: true

    input:
        file("*.stat.txt") from fq_bam_stat_files.toSortedList()

    output:
        file("fq_bam_stats.tsv")

    """
        echo -e "ID\\tvariable\\tvalue\\tcomment" > fq_bam_stats.tsv
        cat *.stat.txt >> fq_bam_stats.tsv
    """
}

process merge_bam {

    cpus params.cores

    tag { SM }

    input:
        set SM, file(bam), file(index) from SM_aligned_bams.groupTuple()

    output:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") into merged_bam_set
        file("${SM}.duplicates.txt") into duplicates_file

    """
    ls -al 1>&2
    count=`ls $bam | wc -l`
    if [ "\${count}" -eq "1" ]; then
        ln -s *.bam ${SM}.merged.bam
        ln -s *.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${task.cpus} --show-progress ${SM}.merged.bam $bam
        sambamba index --nthreads=${task.cpus} ${SM}.merged.bam
    fi
    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${task.cpus} ${SM}.bam
    """

}

merged_bam_set.into {
    subsampled_bams_for_coverage;
    bams_idxstats;
    bams_stats;
    merged_bams_for_coverage;
    merged_bams_individual;
    merged_bams_union
}

process subsampled_coverage_SM {

    tag { SM }

    label 'bam_coverage'

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from subsampled_bams_for_coverage

    output:
        file("${SM}.coverage.tsv") into subsampled_SM_coverage

    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}

process sampled_coverage_SM_merge {

    publishDir "${params.out}/strain", mode: 'copy', overwrite: true

    input:
        file(sm_set) from subsampled_SM_coverage.collect()

    output:
        file("SM_coverage_subsampled.tsv")

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6 | sort) > SM_coverage_subsampled.tsv
    """
}

process SM_idx_stats {
    
    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bams_idxstats
    output:
        file('bam_idxstats.txt') into bam_idxstats_set

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > bam_idxstats.txt
    """
}

process SM_combine_idx_stats {

    publishDir "${params.out}/strain", mode: 'copy', overwrite: true

    input:
        file("*.stat.txt") from bam_idxstats_set.toSortedList()

    output:
        file("SM_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
        cat *.stat.txt | sort >> SM_bam_idxstats.tsv
    """

}

process SM_bam_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bams_stats

    output:
        file('bam_stat.txt') into SM_bam_stat_files

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat.txt
    """
}

process combine_SM_bam_stats {

    publishDir "${params.out}/strain", mode: 'copy', overwrite: true

    input:
        file("?.stat.txt") from SM_bam_stat_files.toSortedList()

    output:
        file("SM_bam_stats.tsv")

    """
        echo -e "ID\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
        cat *.stat.txt | sort >> SM_bam_stats.tsv
    """
}

process format_duplicates {

    publishDir "${params.out}/duplicates", mode: 'copy', overwrite: true

    input:
        file(duplicates_set) from duplicates_file.toSortedList()

    output:
        file("bam_duplicates.tsv")

    """
        echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > bam_duplicates.tsv
        for i in ${duplicates_set.join(" ")}; do
            f=\$(basename \${i})
            cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> bam_duplicates.tsv
        done;
    """
}

process coverage_SM {

    publishDir "${params.out}/coverage_SM", mode: 'copy', overwrite: true

    tag { SM }

    label 'bam_coverage'

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_for_coverage

    output:
        file("${SM}.coverage.tsv") into SM_coverage
        file("${SM}.1mb.coverage.tsv") into SM_1mb_coverage
        file("${SM}.100kb.coverage.tsv") into SM_100kb_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
        bam coverage --window=1000000 ${SM}.bam > ${SM}.1mb.coverage.tsv
        bam coverage --window=100000 ${SM}.bam > ${SM}.100kb.coverage.tsv
    """
}

process coverage_SM_merge {

    publishDir "${params.out}/strain", mode: 'copy', overwrite: true

    input:
        file(sm_set) from SM_coverage.toSortedList()

    output:
        file("SM_coverage.full.tsv")
        file("SM_coverage.tsv") into SM_coverage_merged

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6 | sort) > SM_coverage.tsv
    """
}

SM_coverage_merged.into {
    for_concordance;
    for_strain_list
}

process coverage_bins_merge {

    publishDir "${params.out}/strain", mode: 'copy', overwrite: true

    input:
        file(mb) from SM_1mb_coverage.toSortedList()
        val kb_100 from SM_100kb_coverage.toSortedList()

    output:
        file("SM_coverage.mb.tsv.gz")

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.mb.tsv
        cat ${mb.join(" ")} >> SM_coverage.mb.tsv
        gzip SM_coverage.mb.tsv
    """
}

process call_variants_individual {

    cpus params.cores

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_individual
        file("${params.genome}/*") from reference2.collect()

    output:
        file("${SM}.individual.sites.tsv") into individual_sites

    """
    # Perform individual-level calling
    contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
    echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${task.cpus} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${params.genome}/${params.genome}.fa.gz ${SM}.bam | bcftools call --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${SM}.{}.individual.vcf.gz"
    order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".individual.vcf.gz" }'`
    
    # Output variant sites
    bcftools concat \${order} -O v| bcftools view -O z > ${SM}.individual.vcf.gz
    bcftools index ${SM}.individual.vcf.gz
    rm \${order}

    bcftools view -M 2 -m 2 -O v ${SM}.individual.vcf.gz | \\
    bcftools filter --include 'FORMAT/DP > 3' | \\
    egrep '(^#|1/1)' | \\
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' > ${SM}.individual.sites.tsv
    """
}

process merge_variant_list {

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        file(sites) from individual_sites.toSortedList()

    output:
        file("sitelist.tsv.gz") into gz_sitelist
        file("sitelist.tsv.gz") into sitelist_stat
        file("sitelist.tsv.gz.tbi") into gz_sitelist_index


    """
        echo ${sites}
        cat ${sites.join(" ")} | sort -k1,1 -k2,2n | uniq > sitelist.tsv
        bgzip sitelist.tsv -c > sitelist.tsv.gz && tabix -s1 -b2 -e2 sitelist.tsv.gz
    """
}

union_vcf_set = merged_bams_union.combine(gz_sitelist).combine(gz_sitelist_index)

process call_variants_union {

    cpus params.cores

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from union_vcf_set
        file("${params.genome}/*") from reference3.collect()

    output:
        file("${SM}.union.vcf.gz") into union_vcf_to_list

    """
        contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
        echo \${contigs} | \\
        tr ' ' '\\n' | \\
        xargs --verbose -I {} -P ${task.cpus} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,INFO/AD,SP --fasta-ref ${params.genome}/${params.genome}.fa.gz ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`

        # Output variant sites
        bcftools concat \${order} -O v | \\
        # vk geno het-polarization - | \\
        bcftools filter -O u --threads ${task.cpus} --set-GTs . --include "QUAL >= ${qual} || FORMAT/GT == '0/0'" |  \\
        bcftools filter -O u --threads ${task.cpus} --set-GTs . --include "FORMAT/DP > ${min_depth}" | \\
        bcftools filter -O u --threads ${task.cpus} --set-GTs . --include "INFO/MQ > ${mq}" | \\
        bcftools filter -O u --threads ${task.cpus} --set-GTs . --include "(FORMAT/AD[*:1])/(FORMAT/DP) >= ${dv_dp} || FORMAT/GT == '0/0'" | \\
        bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """
}

union_vcf_to_list.into { union_vcf_to_list1; union_vcf_to_list2 }

process generate_union_vcf_list {

    cpus 1 

    publishDir "${params.out}/variation", mode: 'copy'

    input:
       file(vcf_set) from union_vcf_to_list1.toSortedList()

    output:
       file("union_vcfs.txt") into union_vcfs
       file("*.tbi") into union_vcfs_indexes

    """
        echo ${vcf_set.join(" ")} | tr ' ' '\\n' > union_vcfs.txt
        for file in ${vcf_set.join(" ")}; do
            tabix -fp vcf \$file
        done
    """
}

union_vcfs_in = union_vcfs.combine(contig_list)

process merge_union_vcf_chromosome {

    cpus params.cores

    tag { chrom }

    input:
        set file(union_vcfs:"union_vcfs.txt"), val(chrom) from union_vcfs_in
        file(union_vcfs_2) from union_vcf_to_list2.collect()
        file(indexes) from union_vcfs_indexes.collect()

    output:
        file("${chrom}.merged.raw.vcf.gz") into raw_vcf

    """
        ls -al 1>&2
        bcftools merge --regions ${chrom} -O z -m all ${union_vcfs_2.join(" ")} > ${chrom}.merged.raw.vcf.gz
        bcftools index ${chrom}.merged.raw.vcf.gz
    """
}

// Generate a list of ordered files.

process concatenate_union_vcf {

    cpus params.cores

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        file(merge_vcf) from raw_vcf.toSortedList()

    output:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") into raw_vcf_concatenated

    """
        bcftools concat --threads ${task.cpus} -O z ${merge_vcf.join(" ")}  > merged.raw.vcf.gz
        bcftools index --threads ${task.cpus} merged.raw.vcf.gz
    """
}

process filter_union_vcf {

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") from raw_vcf_concatenated

    output:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") into filtered_vcf
    
    """
        bcftools view merged.raw.vcf.gz | \\
        bcftools filter --mode=+x --include 'F_MISSING  <= 0.1' - | \\
        bcftools view -m1 - |\\
        bcftools view -O z - > concordance.vcf.gz
        bcftools index concordance.vcf.gz
    """
}

filtered_vcf.into {
    filtered_vcf_gtcheck;
    filtered_vcf_stat;
    filtered_vcf_pairwise;
    het_check_vcf;
    strain_pairwise_vcf;
    npr1_allele
}

process calculate_gtcheck {

    publishDir "${params.out}/concordance", mode: 'copy'

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from filtered_vcf_gtcheck

    output:
        file("gtcheck.tsv") into gtcheck
        file("gtcheck.tsv") into pairwise_compare_gtcheck

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > gtcheck.tsv
        bcftools gtcheck -H -G 1 concordance.vcf.gz | egrep '^CN' | cut -f 2-6 >> gtcheck.tsv
    """
}

process stat_tsv {

    publishDir "${params.out}/vcf", mode: 'copy'

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from filtered_vcf_stat

    output:
        file("concordance.stats") into filtered_stats

    """
        bcftools stats --verbose concordance.vcf.gz > concordance.stats
    """
}

process process_concordance_results {

    publishDir "${params.out}/concordance", mode: "copy"

    input:
        file 'gtcheck.tsv' from gtcheck
        file 'filtered.stats.txt' from filtered_stats
        file 'SM_coverage.tsv' from for_concordance

    output:
        file("concordance.pdf")
        file("concordance.png")
        file("xconcordance.pdf")
        file("xconcordance.png")
        file("isotype_groups.tsv") into isotype_groups_ch
        file("isotype_count.txt")
        file("WI_metadata.tsv")

    """
    # Run concordance analysis
    process_concordance.R ${params.debug}
    """
}

isotype_groups_ch.into { isotype_groups; for_combined_final}

process generate_isotype_groups {

    input:
        file("isotype_groups.tsv") from isotype_groups

    output:
        file("pairwise_groups.txt") into pairwise_groups

    """
    cat isotype_groups.tsv | awk '{ curr_strain = \$2; curr_group = \$1; if (group_prev == curr_group) { print prev_strain "," curr_strain "\t" \$1 "\t" \$3 } ; prev_strain = \$2; group_prev = \$1; }' > pairwise_groups.txt
    """

}

pairwise_groups_input = pairwise_groups.splitText( by:1 )

// Look for diverged regions among isotypes.
process pairwise_variant_compare {

    publishDir "${params.out}/concordance/pairwise/within_group", mode: 'copy', overwrite: true

    tag { pair }

    input:
        val(pair_group) from pairwise_groups_input
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from filtered_vcf_pairwise

    output:
        file("${group}.${isotype}.${pair.replace(",","_")}.png")
        file("${group}.${isotype}.${pair.replace(",","_")}.tsv")

    script:
        pair_group = pair_group.trim().split("\t")
        pair = pair_group[0]
        group = pair_group[1]
        isotype = pair_group[2]

    """
        bcftools query -f '%CHROM\t%POS[\t%GT]\n' -s ${pair} concordance.vcf.gz > out.tsv
        Rscript --vanilla `which plot_pairwise.R` 
        mv out.png ${group}.${isotype}.${pair.replace(",","_")}.png
        mv out.tsv ${group}.${isotype}.${pair.replace(",","_")}.tsv
    """
}

process heterozygosity_check {

    cpus params.cores

    publishDir "${params.out}/concordance", mode: "copy"

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from het_check_vcf

    output:
        file("heterozygosity.tsv")

    """
        bcftools query -l concordance.vcf.gz | xargs --verbose -I {} -P ${task.cpus} sh -c "bcftools query -f '[%SAMPLE\t%GT\n]' --samples={} concordance.vcf.gz | grep '0/1' | uniq -c >> heterozygosity.tsv"
    """

}

// The belows are new processes for futher checking

process strain_pairwise_list {

    //executor 'local'

    //publishDir "${params.out}/concordance/pairwise/between_strains", mode: "copy"

    input:
        file("SM_coverage.tsv") from for_strain_list

    output:
        file("strain_pairwise_list.tsv") into strain_pairwise

    """
        # generate strain level pairwise comparison list
        cat SM_coverage.tsv | cut -f1 | sed '1d' > raw_strain.tsv

        for i in `cat raw_strain.tsv` ; do
            for j in `cat raw_strain.tsv` ; do
                if [ "\$i" '<' "\$j" ] ; then
                    echo \${i}-\${j}
                fi
            done
        done | tr "-" "\t" > strain_pairwise_list.tsv
    """
}

strain_pairwise.splitText( by:1 )
               .set { new_strain_pairwise }


process query_between_group_pairwise_gt {

    publishDir "${params.out}/variation", mode: 'copy', overwrite: true

    cpus params.cores

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from strain_pairwise_vcf

    output:
        file("out_gt.tsv") into gt_pairwise

    """
        # query the genotype for pairwise comparison
        bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n' concordance.vcf.gz -H | sed 's/\\:GT//g' | sed 's/\\[[0-9]*\\]//g' | sed 's/\\#//g' | csvtk rename -t -f 1 -n CHROM > out_gt.tsv
    """
}

process between_group_pairwise {

    publishDir "${params.out}/concordance/pairwise/between_group", mode: 'copy', overwrite: true, pattern: '*.png'

    tag "${sp1}_${sp2}"

    input:
        val(pair_group) from new_strain_pairwise
        file("out_gt.tsv") from gt_pairwise

    output:
        file("${sp1}-${sp2}.tsv") into between_group_pairwise_out
        file("${sp1}-${sp2}-distribution.tsv") into cutoff_distribution
        file("${sp1}-${sp2}.disconcordance.png") optional true
        file("${sp1}-${sp2}.hist.png") optional true
    
    script:
        pair_group = pair_group.trim().split("\t")
        sp1 = pair_group[0]
        sp2 = pair_group[1]

    """            
        csvtk cut -t -f CHROM,POS,${sp1},${sp2} out_gt.tsv > ${sp1}-${sp2}.queried.tsv
        process_strain_pairwise.R ${sp1} ${sp2} ${sp1}-${sp2}.queried.tsv
        mv condition_results.tsv ${sp1}-${sp2}.tsv
        mv for_distribution.tsv ${sp1}-${sp2}-distribution.tsv
        rm ${sp1}-${sp2}.queried.tsv
    """
}

process npr1_allele_check {

    cpus params.cores

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        set file("concordance.vcf.gz"), file("concordance.vcf.gz.csi") from npr1_allele

    output:
        file("npr1_allele_strain.tsv") into npr1_out

    """
        echo -e 'problematic_strain\\tgt' > npr1_allele_strain.tsv
        bcftools view --threads ${params.cores} -t X:4768788 concordance.vcf.gz | bcftools query -f '[%SAMPLE\\t%GT\\n]' | awk '\$2 != "1/1"' >> npr1_allele_strain.tsv
    """
}

process merge_betweengroup_pairwise_output {

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        file(bg_pairwise) from between_group_pairwise_out.toSortedList()

    output:
        file("merge_betweengroup_pairwise_output.tsv") into combine_pairwise_results_ch


    """
        echo ${bg_pairwise}
        echo -e 'pairwise\\tconcordant_bin_gt_70\\tmax_discordant_bin_count_lt_3\\tmean_discordant_bin_count_lt_2.5\\tno_bin_lt_0.9\\tsuspected_introgress' > merge_betweengroup_pairwise_output.tsv
        cat ${bg_pairwise.join(" ")} | cut -f 1,2,3,4,5,6 >> merge_betweengroup_pairwise_output.tsv
    """
}

process cutoff_distribution {

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        file(cutoff_val) from cutoff_distribution.toSortedList()

    output:
        file("cutoff_distribution.tsv")


    """
        echo ${cutoff_val}
        ls -al 1>&2
        echo -e 'pairwise\\tconcordant_bin_gt_70\\tmax_discordant_bin_count_lt_3\\tmean_discordant_bin_count_lt_2.5\\tmin_bin' > cutoff_distribution.tsv
        cat ${cutoff_val.join(" ")} | cut -f 1,2,3,4,5 >> cutoff_distribution.tsv
    """
}

process combine_pairwise_results {

    publishDir "${params.out}/concordance", mode: 'copy', overwrite: true

    input:
        file("isotype_groups.tsv") from for_combined_final
        file("merge_betweengroup_pairwise_output.tsv") from combine_pairwise_results_ch
        file("npr1_allele_strain.tsv") from npr1_out

    output:
        file("new_isotype_groups.tsv")

    """
        merge_groups_info.R isotype_groups.tsv merge_betweengroup_pairwise_output.tsv npr1_allele_strain.tsv
    """
}

workflow.onComplete {

    user="whoami".execute().text

    summary = """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    User: ${user}
    """

    println summary

    // mail summary
    ['mail', '-s', 'concordance-nf', params.email].execute() << summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << summary
        outlog << "\n--------pyenv-------\n"
        outlog << "pyenv versions".execute().text
        outlog << "--------ENV--------"
        outlog << "ENV".execute().text
        outlog << "--------brew--------"
        outlog << "brew list".execute().text
        outlog << "--------R--------"
        outlog << "Rscript -e 'devtools::session_info()'".execute().text
    }
}

