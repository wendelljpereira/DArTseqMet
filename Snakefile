configfile: "config.yaml"

rule all:
    input:
        expand("methylated_sites/{generated_bed_file_of_DNA_methylations}", generated_bed_file_of_DNA_methylations=config["generated_bed_file_of_DNA_methylations"])

rule barcode_removing:
    input:
        expand("{barcodes_files}", barcodes_files=config["barcodes_files"]),
        expand("{adapters_file}", adapters_file=config["adapters_file"]),
        expand("{fastq}.FASTQ", fastq=config["fastq"])
    threads: config["available_cores"]
    output:
        expand("{fastq}.without_barcodes_and_adapters.fq", fastq=config["fastq"])
    shell:
        """
        # defines the input files

        barcode_file1={input[0]}
        barcode_file2={input[1]}
        adapters_file={input[2]}

        # execute using the two sets of barcodes

        while read FILE barcode_sequence barcode_size
        do
        trimmomatic SE -threads {threads} $FILE.FASTQ $FILE.without_barcodes_and_adapters.fq ILLUMINACLIP:$adapters_file:2:30:6 HEADCROP:$barcode_size MINLEN:20
        done < "$barcode_file1"

        while read FILE barcode_sequence barcode_size
        do
        trimmomatic SE -threads {threads} $FILE.FASTQ $FILE.without_barcodes_and_adapters.fq ILLUMINACLIP:$adapters_file:2:30:6 HEADCROP:$barcode_size MINLEN:20
        done < "$barcode_file2"

        """

rule fastqc:
    input:
        expand("{fastq}.without_barcodes_and_adapters.fq", fastq=config["fastq"])
    output:
        expand("fastqc_after_trimming/{fastq}.without_barcodes_and_adapters_fastqc.html", fastq=config["fastq"]),
        expand("fastqc_after_trimming/{fastq}.without_barcodes_and_adapters_fastqc.zip", fastq=config["fastq"])
    shell:
        """
        fastqc *.fq
        mv *fastqc.zip fastqc_after_trimming
        mv *fastqc.html fastqc_after_trimming
        """

rule combining_samples:
    input:
        expand("{fastq}.without_barcodes_and_adapters.fq", fastq=config["fastq"])
    output:
        expand("{combining_samples}", combining_samples=config["combining_samples"])
    shell:
        "cat {input} > {output}"

rule mapping:
    input:
        expand("{reference_genome}", reference_genome=config["reference_genome"]),
        expand("{combining_samples}", combining_samples=config["combining_samples"])
    conda: "bowtie2.yaml" 
    threads: config["available_cores"]
    output:
        expand("mapping/{mapping}.fq", mapping=config["mapping"]["fastq"]),
        expand("mapping/{mapping}.sam", mapping=config["mapping"]["sam"])
    shell:
        """

            base=$(basename {input[0]} ".fa")

            bowtie2-build -f --threads {threads} {input[0]} $base

            bowtie2 -p {threads} -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --score-min L,0,0 -x $base --un {output[0]} -U {input[1]} -S {output[4]}

            bowtie2 -p {threads} -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-2,0 -x $base --un {output[1]} -U {output[0]} -S {output[5]}

            bowtie2 -p {threads} -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-4,0 -x $base --un {output[2]} -U {output[1]} -S {output[6]}

            bowtie2 -p {threads} -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-6,0 -x $base --un {output[3]} -U {output[2]} -S {output[7]}

        """

rule extracts_single_mapping_reads:
    input:
        expand("mapping/{mapping}.sam", mapping=config["mapping"]["sam"])
    threads: config["available_cores"]
    output:
        expand("mapping/{mapping}.single_map.sam", mapping=config["mapping"]["sam"])
    shell:
        """

            for sample in `ls mapping/*.sam`
            do
            base=$(basename $sample \".sam\")

            samtools view -@ {threads} -h mapping/${{base}}.sam | grep -v \"XS:i:\" > mapping/${{base}}.single_map.sam
            done

        """

rule sam_to_bam:
    input:
        expand("mapping/{mapping}.single_map.sam", mapping=config["mapping"]["sam"]),
    threads: config["available_cores"]
    output:
        expand("mapping/{mapping}.single_map.bam", mapping=config["mapping"]["sam"])
    shell:
        """

            for sample in `ls mapping/*single_map.sam`
            do

            base=$(basename $sample \".sam\")

            samtools view -@ {threads} -b mapping/${{base}}.sam > mapping/${{base}}.bam
            samtools sort -@ {threads} mapping/${{base}}.bam > mapping/tmp.bam

            mv mapping/tmp.bam mapping/${{base}}.bam
            done

        """

rule merge_bam:
    input:
        expand("mapping/{mapping}.single_map.bam", mapping=config["mapping"]["sam"])
    threads: config["available_cores"]
    output:
        expand("mapping/{combine_bam}", combine_bam=config["combine_bam"])
    shell:
        """

            samtools merge -@ {threads} -f {output[0]} {input[0]} {input[1]} {input[2]} {input[3]}
            samtools sort -@ {threads} {output[0]} > mapping/temp.bam

            mv mapping/temp.bam {output[0]}

        """

rule bam_to_bed:
    input:
        expand("mapping/{combine_bam}", combine_bam=config["combine_bam"])
    output:
        expand("bed_files/{bed_name}", bed_name=config["bam_to_bed"]["bed_name"])
    shell:
        """
            #conversion from bam to bed

            #Deletes intermediate files

            #rm mapping/*.sam
            #rm mapping/*map.bam

            bamToBed -i {input[0]} > {output[0]}
            bedtools sort -i {output[0]} > tmp.bed
            mv tmp.bed {output[0]}

        """

rule find_the_find_methylation_sites_position:
    input:
        expand("{reference_genome}", reference_genome=config["reference_genome"]),
        expand("{restriction_site_sequences}", restriction_site_sequences=config["restriction_site_sequences"])
    output:
       expand("{position_of_RE_sites}", position_of_RE_sites=config["position_of_RE_sites"])
    shell:
        "Rscript finding_restriction_site_2.0.R --out1 {output[0]} {input[0]} {input[1]}"

rule defines_the_position_of_sampled_sites:
    input:
        expand("bed_files/{bed_name}", bed_name=config["bam_to_bed"]["bed_name"]),
        expand("{position_of_RE_sites}", position_of_RE_sites=config["position_of_RE_sites"])
    output:
        expand("sample_sites_generation/{defines_the_position_of_sampled_sites}", defines_the_position_of_sampled_sites=config["defines_the_position_of_sampled_sites"])
    shell:
        """

            cat {input[0]} {input[1]} > bed_files/enz_sites_sample_reads.bed

            bedtools sort -i bed_files/enz_sites_sample_reads.bed > bed_files/tmp.bed

            mv bed_files/tmp.bed sample_sites_generation/enz_sites_sample_reads.bed

            bedtools cluster -d -1 -s -i sample_sites_generation/enz_sites_sample_reads.bed > {output[0]}

            Rscript split_clusters.R --out1 {output[1]} {output[0]}

        """

rule finds_unique_positions:
    input:
        expand("sample_sites_generation/{defines_the_position_of_sampled_sites}", defines_the_position_of_sampled_sites=config["defines_the_position_of_sampled_sites"])[1]
    output:
        expand("multiple_bed/{finds_unique_positions}", finds_unique_positions=config["finds_unique_positions"])[0:2],
        expand("bed_files/{finds_unique_positions}", finds_unique_positions=config["finds_unique_positions"])[2]
    shell:
        """

            cd multiple_bed
            rm -rf ./clusters_without_pstI/
            rm -rf ./good_clusters/

            awk '{{print >> $7; close($7)}}' ./../{input[0]} 

            Rscript ./../find_unique_positions_from_splited_clusters.R

            cd good_clusters/

            ls |  while read filename;  do cat $filename;  done >./../good_clusters.bed

            cd ..

            cd clusters_without_pstI

            ls |  while read filename;  do cat $filename;  done >./../clusters_without_pstI.bed

            cd ./../../

            cat {output[0]} {output[1]} > {output[2]}

            bedtools sort -i {output[2]} > temp.bed
            mv temp.bed {output[2]}

            bedtools cluster -d -1 -s -i bed_files/msdartseq_positions.bed > bed_files/tmp.bed
            mv bed_files/tmp.bed bed_files/msdartseq_positions.bed

            bedtools sort -i {output[2]} > temp.bed
            mv temp.bed {output[2]}

            awk 'BEGIN{{OFS="\t"}}; {{name="MS-DArT_site_"NR; print $1,$2,$3,name,$5,$6,$8}}' {output[2]} > temp.bed
            mv temp.bed {output[2]}
        """

rule trimmomatic_individual_libraries:
    input:
        expand("{barcodes_files}", barcodes_files=config["barcodes_files"]),
        expand("{adapters_file}", adapters_file=config["adapters_file"]),
        expand("{fastq}.FASTQ", fastq=config["fastq"])
    threads: config["available_cores"]
    output:
        expand("{fastq}.without_barc_and_adapt_plus_qc.fq", fastq=config["fastq"])
    shell:
        """

            # defines the input files

            barcode_file1={input[0]}
            barcode_file2={input[1]}
            adapters_file={input[2]}

            # execute using the two sets of barcodes

            while read FILE barcode_sequence barcode_size
            do
            trimmomatic SE -threads {threads} $FILE.FASTQ $FILE.without_barc_and_adapt_plus_qc.fq ILLUMINACLIP:$adapters_file:2:30:6 SLIDINGWINDOW:5:25 HEADCROP:$barcode_size MINLEN:20
            done < "$barcode_file1"

            while read FILE barcode_sequence barcode_size
            do
            trimmomatic SE -threads {threads} $FILE.FASTQ $FILE.without_barc_and_adapt_plus_qc.fq ILLUMINACLIP:$adapters_file:2:30:6 SLIDINGWINDOW:5:25 HEADCROP:$barcode_size MINLEN:20
            done < "$barcode_file2"
        """

rule fastqc_individual_libraries:
    input:
        expand("{fastq}.without_barc_and_adapt_plus_qc.fq", fastq=config["fastq"])
    output:
        expand("fastqc_after_trimming/{fastq}.without_barc_and_adapt_plus_qc_fastqc.html", fastq=config["fastq"]),
        expand("fastqc_after_trimming/{fastq}.without_barc_and_adapt_plus_qc_fastqc.zip", fastq=config["fastq"])
    shell:
        "fastqc *_qc.fq;"
        "mv *fastqc.zip fastqc_after_trimming;"
        "mv *fastqc.html fastqc_after_trimming"

rule mapping_individual_libraries:
    input:
        expand("{reference_genome}", reference_genome=config["reference_genome"]),
        expand("{fastq}.without_barc_and_adapt_plus_qc.fq", fastq=config["fastq"])
    conda: "bowtie2.yaml" 
    threads: config["available_cores"]
    output:
        expand("mapping/{fastq}_{map}.fq", fastq=config["fastq"], map=config["mapping"]["fastq"]),
        expand("mapping/{fastq}{sam_samples}.sam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    shell:
        """

            base=$(basename {input[0]} ".fa")

            bowtie2-build -f --threads {threads} {input[0]} $base

        # 0 mismatch
        for sample in `ls ./*_qc.fq`
        do

        base=$(basename $sample ".without_barc_and_adapt_plus_qc.fq")
        echo $base
        bowtie2 -p {threads} -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --score-min L,0,0 -x Egrandis_297_v2.0.softmasked --un mapping/${{base}}_no_mapped_mm0.fq -U $sample -S mapping/${{base}}_mm0.sam

        done

        #1 mismatch
        for sample in `ls mapping/*_no_mapped_mm0.fq`
        do
        base=$(basename $sample "_no_mapped_mm0.fq")

        echo $base
        bowtie2 -p {threads} -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-2,0 -x Egrandis_297_v2.0.softmasked --un mapping/${{base}}_no_mapped_mm1.fq -U mapping/${{base}}_no_mapped_mm0.fq -S mapping/${{base}}_mm1.sam
        done

        #2 mismatches
        for sample in `ls mapping/*_no_mapped_mm1.fq`
        do
        base=$(basename $sample "_no_mapped_mm1.fq")
        echo $base
        bowtie2 -p {threads} -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-4,0 -x Egrandis_297_v2.0.softmasked --un mapping/${{base}}_no_mapped_mm2.fq -U mapping/${{base}}_no_mapped_mm1.fq -S mapping/${{base}}_mm2.sam
        done

        #3 mismatches
        for sample in `ls mapping/*_no_mapped_mm2.fq`
        do
        base=$(basename $sample "_no_mapped_mm2.fq")

        echo $base
        bowtie2 -p {threads} -q -k 2 -N 0 -L 15 -i S,1,0.50 -R 3 --no-1mm-upfront --end-to-end --no-unal --mp 2 --np 10 --rdg 10 --rfg 10 --score-min L,-6,0 -x Egrandis_297_v2.0.softmasked --un mapping/${{base}}_no_mapped_mm3.fq -U mapping/${{base}}_no_mapped_mm2.fq -S mapping/${{base}}_mm3.sam
        done

        """

rule extracts_single_mapping_reads_individual_libraries:
    input:
        expand("mapping/{fastq}{sam_samples}.sam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    threads: config["available_cores"]
    output:
        expand("mapping/{fastq}{sam_samples}.single_map.sam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    shell:
        """
            for sample in `ls mapping/*.sam`
            do
            base=$(basename $sample \".sam\")

            samtools view -@ {threads} -h mapping/${{base}}.sam | grep -v \"XS:i:\" > mapping/${{base}}.single_map.sam
            done

        """

rule sam_to_bam_individual_libraries:
    input:
        expand("mapping/{fastq}{sam_samples}.single_map.sam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    threads: config["available_cores"]
    output:
        expand("mapping/{fastq}{sam_samples}.single_map.bam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    shell:
        """

            for sample in `ls mapping/*single_map.sam`
            do
            base=$(basename $sample \".sam\")
            samtools view -@ {threads} -b mapping/${{base}}.sam > mapping/${{base}}.bam
            samtools sort -@ {threads} mapping/${{base}}.bam > mapping/tmp.bam
            mv mapping/tmp.bam mapping/${{base}}.bam
            done

        """

rule merge_bam_individual_libraries:
    input:
        expand("mapping/{fastq}{sam_samples}.single_map.bam", fastq=config["fastq"], sam_samples=config["mapping"]["sam_samples"])
    threads: config["available_cores"]
    output:
        expand("mapping/{fastq}_combined.bam", fastq=config["fastq"])
    shell:
        """

            for sample in `ls mapping/*_mm0.single_map.bam`
            do
            base=$(basename $sample "_mm0.single_map.bam")

            samtools merge -@ {threads} -f mapping/${{base}}_combined.bam mapping/${{base}}_mm0.single_map.bam mapping/${{base}}_mm1.single_map.bam mapping/${{base}}_mm2.single_map.bam mapping/${{base}}_mm3.single_map.bam

            samtools sort -@ {threads} mapping/${{base}}_combined.bam > tmp.bam
            mv tmp.bam mapping/${{base}}_combined.bam

            done

        """

rule prepare_featureCounts_input:
    input:
        expand("bed_files/{finds_unique_positions}", finds_unique_positions=config["finds_unique_positions"])[2]
    output:
        expand("bed_files/{finds_unique_positions}", finds_unique_positions=config["featurecounts"])[0]
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}; {{print $4,$1,$2,$3,$6}}' {input[0]} > {output[0]}
        """

rule featureCounts:
    input:
        expand("bed_files/{finds_unique_positions}", finds_unique_positions=config["featurecounts"])[0],
        expand("mapping/{fastq}_combined.bam", fastq=config["fastq"])
    threads: config["available_cores"]
    output:
        expand("counts/{finds_unique_positions}.tst", finds_unique_positions=config["featurecounts"])[1]
    shell:
        """
        declare -a list_of_bams
        list_of_bams=({input})
        list_of_bams="${{list_of_bams[@]:1}}"

        featureCounts -f -F SAF -R CORE -s 1 -T {threads} -O -a {input[0]} -o {output[0]} $list_of_bams

        sed '1d' {output[0]} > tmp_counts
        mv tmp_counts {output[0]}
        """

rule counts_correction:
    input:
        expand("bed_files/{finds_unique_positions}", finds_unique_positions=config["finds_unique_positions"])[2],
        expand("counts/{finds_unique_positions}.tst", finds_unique_positions=config["featurecounts"])[1]
    output:
        expand("bed_files/{counts_correction}", counts_correction=config["counts_correction"])[0],
        expand("counts/{counts_correction}", counts_correction=config["counts_correction"])[1]
    shell:
        "Rscript counts_correction.R --out1 {output[0]} --out2 {output[1]} {input[0]} {input[1]};"

rule selects_true_counts:
    input:
        expand("counts/{counts_correction}", counts_correction=config["counts_correction"])[1]
    params:
        samples_from_ms_insenstive=expand("{samples_from_ms_insenstive}", samples_from_ms_insenstive=config["samples_from_ms_insenstive"]),
    output:
        expand("true_sites/{output_bigger_than_0}.txt", output_bigger_than_0=config["marks_with_msp_outputs"])
    shell:
        """
        Rscript selects_true_counts.R --ms_insensitive_samples {params.samples_from_ms_insenstive} --out1 {output[0]} {input[0]}
        """

rule sampled_site_position:
    input:
        expand("bed_files/{counts_correction}", counts_correction=config["counts_correction"])[0],
        "enzymes_positions_of_RE_sites_sorted.bed",
        expand("counts/{counts_correction}", counts_correction=config["counts_correction"])[1]
    output:
        expand("position_of_the_sampled_sites/{methylatio_site_output}.bed", methylatio_site_output=config["sampled_site_position"])[0],
        expand("position_of_the_sampled_sites/{methylatio_site_output}.bed", methylatio_site_output=config["sampled_site_position"])[1],
        expand("{methylatio_site_output}.csv", methylatio_site_output=config["sampled_site_position"])[2]
    shell:
        "Rscript marks_closest_restriction_site_search.R --out1 {output[0]} --out2 {output[1]} --out3 {output[2]} {input[0]} {input[1]} {input[2]}"

rule sampled_site_position_part_2:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}.bed", methylatio_site_output=config["sampled_site_position"])[1]
    output:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.tst", methylatio_site_output=config["sampled_site_position"])[1],
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["sampled_site_position"])[1]
    shell:
        """

        bedtools merge -d -1 -s -c 6,2,3,4 -o distinct,collapse,collapse,collapse -delim '|' -i {input[0]} > {output[0]}

        awk 'BEGIN{{OFS=\"\t\"}}; {{print $1,$2,$3,$7,0,$4}}' {output[0]} > tempfile & mv tempfile {output[1]}

        """

rule executes_edgeR:
    input:
        expand("true_sites/{output_bigger_than_0}.txt", output_bigger_than_0=config["marks_with_msp_outputs"]),
        expand("{methylatio_site_output}.csv", methylatio_site_output=config["sampled_site_position"])[2],
    params:
        samples_from_ms_insenstive=expand("{samples_from_ms_insenstive}", samples_from_ms_insenstive=config["samples_from_ms_insenstive"]),
        Grouping_var=expand("{Grouping_var}", Grouping_var=config["Grouping_var"]),
        prefix=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["prefix"]),
        min_msp=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["min_msp"]),
        fdr=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["fdr"]),
        log_fold_change=expand("{edgeR_params}", edgeR_params=config["edgeR_with_DArTCounts_params"]["log_fold_change"]),
    output:
        expand("{prefix_edgeR}_DE_marks.txt", prefix_edgeR=config["edgeR_with_DArTCounts_params"]["prefix"], groups=config["edgeR_with_DArTCounts_params"]),
    script:
        "new_edgeR_with_DArTCounts.R"

rule executes_DEseq2:
    input:
        expand("true_sites/{output_bigger_than_0}.txt", output_bigger_than_0=config["marks_with_msp_outputs"]),
        expand("{methylatio_site_output}.csv", methylatio_site_output=config["sampled_site_position"])[2]
    params:
        samples_from_ms_insenstive=expand("{samples_from_ms_insenstive}", samples_from_ms_insenstive=config["samples_from_ms_insenstive"]),
        Grouping_var=expand("{Grouping_var}", Grouping_var=config["Grouping_var"]),
        prefix=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["prefix"]),
        min_msp=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["min_msp"]),
        fdr=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["fdr"]),
        log_fold_change=expand("{deseq_params}", deseq_params=config["DEseq2_with_DArTCounts_params"]["log_fold_change"]),
    output:
        expand("{prefix_deseq}_DE_marks.txt", prefix_deseq=config["DEseq2_with_DArTCounts_params"]["prefix"])
    script:
        "new_deseq2_with_DArTCounts.R"

rule compares_edgeR_and_DEseq2:
    input:
        expand("{prefix_edgeR}_DE_marks.txt", prefix_edgeR=config["edgeR_with_DArTCounts_params"]["prefix"]),
        expand("{prefix_deseq}_DE_marks.txt", prefix_deseq=config["DEseq2_with_DArTCounts_params"]["prefix"]),
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["sampled_site_position"])[1]
    output:
        expand("compares_edgeR_and_DEseq2/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[0],
        expand("compares_edgeR_and_DEseq2/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1],
        expand("images/compares_edgeR_and_DEseq2/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[2]
    script:
        "edger_vs_deseq2.R"

rule generated_bed_file_of_DNA_methylations:
    input:
        expand("position_of_the_sampled_sites/{methylatio_site_output}_merged.bed", methylatio_site_output=config["sampled_site_position"])[1],
        expand("compares_edgeR_and_DEseq2/{edger_vs_deseq2}", edger_vs_deseq2=config["edger_vs_deseq2_output"])[1]
    output:
        expand("methylated_sites/{generated_bed_file_of_DNA_methylations}", generated_bed_file_of_DNA_methylations=config["generated_bed_file_of_DNA_methylations"])
    script:
        "generated_bed_file_of_DNA_methylations.R"