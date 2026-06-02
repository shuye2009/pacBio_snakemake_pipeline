# =============================================================================
# MethBat rules for methylation analysis
# =============================================================================

import os

SCRIPTS_DIR = os.path.join(workflow.basedir, "scripts")

# =============================================================================
# Convert_bed_to_tsv: Convert BED file to MethBat-compatible TSV format
# 
# MethBat requires a TSV file with headers (chrom, start, end, cpg_label)
# rather than a standard BED file. This rule adds the required header.
# =============================================================================
rule Convert_region_bed_to_tsv:
    input:
        bed=config["methbat"]["regions"],
    output:
        tsv=METHBAT_DIR + "/regions.tsv",
    run:
        import os
        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)
        with open(input.bed, "r") as infile, open(output.tsv, "w") as outfile:
            outfile.write("chrom\tstart\tend\tcpg_label\n")
            for line in infile:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                chrom = fields[0]
                start = fields[1]
                end = fields[2]
                label = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
                outfile.write(f"{chrom}\t{start}\t{end}\t{label}\n")


# =============================================================================
# Methbat_profile: Generate methylation profile for a single sample
# 
# Creates a per-region methylation summary for predefined genomic regions.
# Used for cohort-based methylation analysis.
# =============================================================================
rule Methbat_region_profile:
    input:
        bed=config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
        regions=METHBAT_DIR + "/regions.tsv",
    output:
        profile=METHBAT_DIR + "/profiles_region/{sample}.region.profile.tsv",
    wildcard_constraints:
        sample="|".join([s.replace(".", r"\.") for s in config["samples"]["case"] + config["samples"]["control"]])
    params:
        input_prefix=config["directory"]["output"] + "/pb_cpg_tools/{sample}",
    log:
        config["directory"]["output"] + "/logs/methbat/profile_{sample}.log",
    shell:
        """
        set +u
        source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
        conda activate methbat
        set -u
        mkdir -p $(dirname {output.profile})
        mkdir -p $(dirname {log})
        methbat profile \
            --input-prefix {params.input_prefix} \
            --input-regions {input.regions} \
            --output-region-profile {output.profile} \
            2>&1 | tee {log}
        conda deactivate
        """


# =============================================================================
# Create_profile_collection: Create collection file for cohort analysis
# 
# Generates a TSV file listing all sample profiles with their case/control
# labels. This collection is used by methbat build for cohort analysis.
# =============================================================================
rule Create_region_profile_collection:
    input:
        profiles=expand(
            METHBAT_DIR + "/profiles_region/{sample}.region.profile.tsv",
            sample=ALL_SAMPLES
        ),
    output:
        collection=METHBAT_DIR + "/region_profile_collection.tsv",
    params:
        profile_dir=METHBAT_DIR + "/profiles_region",
        case_samples=config["samples"]["case"],
        control_samples=config["samples"]["control"],
    run:
        import os
        os.makedirs(os.path.dirname(output.collection), exist_ok=True)
        with open(output.collection, "w") as f:
            f.write("identifier\tfilename\tlabels\n")
            for sample in params.case_samples:
                profile = f"{params.profile_dir}/{sample}.region.profile.tsv"
                f.write(f"{sample}\t{profile}\tcase\n")
            for sample in params.control_samples:
                profile = f"{params.profile_dir}/{sample}.region.profile.tsv"
                f.write(f"{sample}\t{profile}\tcontrol\n")


# =============================================================================
# Methbat_build: Build cohort methylation profile from individual profiles
# 
# Aggregates individual sample profiles into a single cohort profile that
# can be used for case vs control comparison.
# =============================================================================
rule Methbat_region_build:
    input:
        collection=METHBAT_DIR + "/region_profile_collection.tsv",
    output:
        cohort_profile=METHBAT_DIR + "/region_cohort.profile.tsv",
    log:
        config["directory"]["output"] + "/logs/methbat/build.log",
    shell:
        """
        set +u
        source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
        conda activate methbat
        set -u
        mkdir -p $(dirname {output.cohort_profile})
        mkdir -p $(dirname {log})
        methbat build \
            --input-collection {input.collection} \
            --output-profile {output.cohort_profile} \
            2>&1 | tee {log}
        conda deactivate
        """


# =============================================================================
# Methbat_compare: Compare methylation between case and control cohorts
# 
# Performs statistical comparison of methylation levels at each region
# between case and control groups. Outputs comparison statistics.
# =============================================================================
rule Methbat_region_compare:
    input:
        cohort_profile=METHBAT_DIR + "/region_cohort.profile.tsv",
    output:
        comparison=METHBAT_DIR + "/region_cohort_comparison.tsv",
    log:
        config["directory"]["output"] + "/logs/methbat/compare.log",
    shell:
        """
        set +u
        source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
        conda activate methbat
        set -u
        mkdir -p $(dirname {output.comparison})
        mkdir -p $(dirname {log})
        methbat compare \
            --input-profile {input.cohort_profile} \
            --output-comparison {output.comparison} \
            --baseline-category control \
            --compare-category case \
            2>&1 | tee {log}
        conda deactivate
        """


# =============================================================================
# Create_collection: Create collection file for signature analysis
# 
# Generates a TSV file listing all sample BED file prefixes with their
# case/control labels. Used by methbat signature for DMR discovery.
# =============================================================================
rule Create_collection:
    input:
        beds=expand(
            config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
            sample=ALL_SAMPLES
        ),
    output:
        collection=METHBAT_BASE + "/collection.tsv",
    params:
        output_dir=config["directory"]["output"] + "/pb_cpg_tools",
        case_samples=config["samples"]["case"],
        control_samples=config["samples"]["control"],
    run:
        import os
        os.makedirs(os.path.dirname(output.collection), exist_ok=True)
        with open(output.collection, "w") as f:
            f.write("identifier\tfilename\tlabels\n")
            for sample in params.case_samples:
                prefix = f"{params.output_dir}/{sample}"
                f.write(f"{sample}\t{prefix}\tcase\n")
            for sample in params.control_samples:
                prefix = f"{params.output_dir}/{sample}"
                f.write(f"{sample}\t{prefix}\tcontrol\n")


# =============================================================================
# Methbat_signature: Identify differentially methylated regions (DMRs)
# 
# Discovers genomic regions with significantly different methylation patterns
# between case and control samples. Does not require predefined regions.
# Requires at least 3 samples total.
# =============================================================================
rule Methbat_signature:
    input:
        collection=METHBAT_BASE + "/collection.tsv",
    output:
        regions=METHBAT_BASE + "/signature.signature_regions.bed",
        stats=METHBAT_BASE + "/signature.signature_stats.tsv",
    params:
        output_prefix=METHBAT_BASE + "/signature",
        threads=config["methbat"]["threads"],
        min_delta=config["methbat"]["min_delta"],
        min_zscore=config["methbat"]["min_zscore"],
        min_sample_frac=config["methbat"]["min_sample_frac"],
    log:
        config["directory"]["output"] + "/logs/methbat/signature.log",
    shell:
        """
        set +u
        source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
        conda activate methbat
        set -u
        mkdir -p $(dirname {params.output_prefix})
        mkdir -p $(dirname {log})
        methbat signature \
            --threads {params.threads} \
            --baseline-category control \
            --compare-category case \
            --input-collection {input.collection} \
            --output-prefix {params.output_prefix} \
            --min-delta {params.min_delta} \
            --min-zscore {params.min_zscore} \
            --min-sample-frac {params.min_sample_frac} \
            2>&1 | tee {log}
        conda deactivate
        """


# =============================================================================
# HAPLOTYPE-SPECIFIC METHYLATION ANALYSIS
# 
# The following rules perform allele-specific (haplotype-resolved) methylation
# analysis using hap1 and hap2 BED files from pb-CpG-tools. These rules are
# only executed when phasing is enabled.
# =============================================================================

if config.get("phasing", {}).get("enabled", False):

    # =========================================================================
    # Symlink_haplotype_beds: Make haplotype files compatible with MethBat
    # 
    # MethBat expects input prefixes to resolve to *.combined.bed.gz
    # pb-CpG-tools outputs *.hap1.bed.gz without the "combined" part.
    # =========================================================================
    rule Symlink_haplotype_beds:
        input:
            bed=config["directory"]["output"] + "/pb_cpg_tools/{sample}.{haplotype}.bed.gz",
            tbi=config["directory"]["output"] + "/pb_cpg_tools/{sample}.{haplotype}.bed.gz.tbi",
        output:
            bed_symlink=config["directory"]["output"] + "/pb_cpg_tools/{sample}.{haplotype}.combined.bed.gz",
            tbi_symlink=config["directory"]["output"] + "/pb_cpg_tools/{sample}.{haplotype}.combined.bed.gz.tbi",
        wildcard_constraints:
            haplotype="hap1|hap2"
        shell:
            """
            cd $(dirname {output.bed_symlink})
            ln -sf $(basename {input.bed}) $(basename {output.bed_symlink})
            ln -sf $(basename {input.tbi}) $(basename {output.tbi_symlink})
            """


    # =========================================================================
    # Methbat_region_profile_haplotype: Generate methylation profile per haplotype
    # =========================================================================
    rule Methbat_region_profile_haplotype:
        input:
            bed_symlink=config["directory"]["output"] + "/pb_cpg_tools/{sample}.{haplotype}.combined.bed.gz",
            regions=METHBAT_DIR + "/regions.tsv",
        output:
            profile=METHBAT_DIR + "/profiles_region_haplotype/{sample}.{haplotype}.region.profile.tsv",
        params:
            input_prefix=config["directory"]["output"] + "/pb_cpg_tools/{sample}",
        wildcard_constraints:
            haplotype="hap1|hap2"
        log:
            config["directory"]["output"] + "/logs/methbat/profile_{sample}_{haplotype}.log",
        shell:
            """
            set +u
            source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
            conda activate methbat
            set -u
            mkdir -p $(dirname {output.profile})
            mkdir -p $(dirname {log})
            methbat profile \
                --input-prefix {params.input_prefix}.{wildcards.haplotype} \
                --input-regions {input.regions} \
                --output-region-profile {output.profile} \
                2>&1 | tee {log}
            conda deactivate
            """


    # =========================================================================
    # Create_haplotype_collection: Create collection for haplotype-specific analysis
    # 
    # Creates separate collections for hap1 and hap2 BED files.
    # =========================================================================
    rule Create_haplotype_collection:
        input:
            bed_symlinks=expand(
                config["directory"]["output"] + "/pb_cpg_tools/{sample}.{{haplotype}}.combined.bed.gz",
                sample=ALL_SAMPLES
            ),
        output:
            collection=METHBAT_BASE + "/collection_{haplotype}.tsv",
        params:
            output_dir=config["directory"]["output"] + "/pb_cpg_tools",
            case_samples=config["samples"]["case"],
            control_samples=config["samples"]["control"],
        wildcard_constraints:
            haplotype="hap1|hap2"
        run:
            import os
            os.makedirs(os.path.dirname(output.collection), exist_ok=True)
            hap = wildcards.haplotype
            with open(output.collection, "w") as f:
                f.write("identifier\tfilename\tlabels\n")
                for sample in params.case_samples:
                    prefix = f"{params.output_dir}/{sample}.{hap}"
                    f.write(f"{sample}_{hap}\t{prefix}\tcase\n")
                for sample in params.control_samples:
                    prefix = f"{params.output_dir}/{sample}.{hap}"
                    f.write(f"{sample}_{hap}\t{prefix}\tcontrol\n")


    # =========================================================================
    # Methbat_signature_haplotype: DMR discovery per haplotype
    # 
    # Identifies differentially methylated regions separately for each haplotype.
    # =========================================================================
    rule Methbat_signature_haplotype:
        input:
            collection=METHBAT_BASE + "/collection_{haplotype}.tsv",
        output:
            regions=METHBAT_BASE + "/signature_{haplotype}.signature_regions.bed",
            stats=METHBAT_BASE + "/signature_{haplotype}.signature_stats.tsv",
        params:
            output_prefix=METHBAT_BASE + "/signature_{haplotype}",
            threads=config["methbat"]["threads"],
            min_delta=config["methbat"]["min_delta"],
            min_zscore=config["methbat"]["min_zscore"],
            min_sample_frac=config["methbat"]["min_sample_frac"],
        wildcard_constraints:
            haplotype="hap1|hap2"
        log:
            config["directory"]["output"] + "/logs/methbat/signature_{haplotype}.log",
        shell:
            """
            set +u
            source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
            conda activate methbat
            set -u
            mkdir -p $(dirname {params.output_prefix})
            mkdir -p $(dirname {log})
            methbat signature \
                --threads {params.threads} \
                --baseline-category control \
                --compare-category case \
                --input-collection {input.collection} \
                --output-prefix {params.output_prefix} \
                --min-delta {params.min_delta} \
                --min-zscore {params.min_zscore} \
                --min-sample-frac {params.min_sample_frac} \
                2>&1 | tee {log}
            conda deactivate
            """


    # =========================================================================
    # Create_region_profile_collection_haplotype: Collection for haplotype cohort analysis
    # =========================================================================
    rule Create_region_profile_collection_haplotype:
        input:
            profiles=expand(
                METHBAT_DIR + "/profiles_region_haplotype/{sample}.{{haplotype}}.region.profile.tsv",
                sample=ALL_SAMPLES
            ),
        output:
            collection=METHBAT_DIR + "/region_profile_collection_{haplotype}.tsv",
        params:
            profile_dir=METHBAT_DIR + "/profiles_region_haplotype",
            case_samples=config["samples"]["case"],
            control_samples=config["samples"]["control"],
        wildcard_constraints:
            haplotype="hap1|hap2"
        run:
            import os
            os.makedirs(os.path.dirname(output.collection), exist_ok=True)
            hap = wildcards.haplotype
            with open(output.collection, "w") as f:
                f.write("identifier\tfilename\tlabels\n")
                for sample in params.case_samples:
                    profile = f"{params.profile_dir}/{sample}.{hap}.region.profile.tsv"
                    f.write(f"{sample}_{hap}\t{profile}\tcase\n")
                for sample in params.control_samples:
                    profile = f"{params.profile_dir}/{sample}.{hap}.region.profile.tsv"
                    f.write(f"{sample}_{hap}\t{profile}\tcontrol\n")


    # =========================================================================
    # Methbat_region_build_haplotype: Build cohort profile per haplotype
    # =========================================================================
    rule Methbat_region_build_haplotype:
        input:
            collection=METHBAT_DIR + "/region_profile_collection_{haplotype}.tsv",
        output:
            cohort_profile=METHBAT_DIR + "/region_cohort_{haplotype}.profile.tsv",
        wildcard_constraints:
            haplotype="hap1|hap2"
        log:
            config["directory"]["output"] + "/logs/methbat/build_{haplotype}.log",
        shell:
            """
            set +u
            source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
            conda activate methbat
            set -u
            mkdir -p $(dirname {output.cohort_profile})
            mkdir -p $(dirname {log})
            methbat build \
                --input-collection {input.collection} \
                --output-profile {output.cohort_profile} \
                2>&1 | tee {log}
            conda deactivate
            """


    # =========================================================================
    # Methbat_region_compare_haplotype: Compare methylation per haplotype
    # =========================================================================
    rule Methbat_region_compare_haplotype:
        input:
            cohort_profile=METHBAT_DIR + "/region_cohort_{haplotype}.profile.tsv",
        output:
            comparison=METHBAT_DIR + "/region_cohort_comparison_{haplotype}.tsv",
        wildcard_constraints:
            haplotype="hap1|hap2"
        log:
            config["directory"]["output"] + "/logs/methbat/compare_{haplotype}.log",
        shell:
            """
            set +u
            source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
            conda activate methbat
            set -u
            mkdir -p $(dirname {output.comparison})
            mkdir -p $(dirname {log})
            methbat compare \
                --input-profile {input.cohort_profile} \
                --output-comparison {output.comparison} \
                --baseline-category control \
                --compare-category case \
                2>&1 | tee {log}
            conda deactivate
            """


    # =========================================================================
    # Allele-Specific Methylation (ASM) Analysis Rules
    #
    # These rules analyze differences in methylation between haplotypes
    # and compare ASM patterns between irradiated (case) and control groups.
    # =========================================================================

    # =========================================================================
    # Create_ASM_collection: Create collection for joint haplotype segmentation
    #
    # Creates separate collection files for case and control groups.
    # Each entry is a sample with its prefix (without haplotype suffix).
    # Methbat joint-segment will automatically look for .hap1 and .hap2 files.
    # =========================================================================
    rule Create_ASM_collection:
        input:
            hap1_beds=lambda wildcards: expand(
                config["directory"]["output"] + "/pb_cpg_tools/{sample}.hap1.combined.bed.gz",
                sample=config["samples"]["case"] if wildcards.group == "case" else config["samples"]["control"]
            ),
            hap2_beds=lambda wildcards: expand(
                config["directory"]["output"] + "/pb_cpg_tools/{sample}.hap2.combined.bed.gz",
                sample=config["samples"]["case"] if wildcards.group == "case" else config["samples"]["control"]
            ),
        output:
            collection=METHBAT_BASE + "/asm_{group}_collection.tsv",
        params:
            output_dir=config["directory"]["output"] + "/pb_cpg_tools",
            samples=lambda wildcards: config["samples"]["case"] if wildcards.group == "case" else config["samples"]["control"],
        wildcard_constraints:
            group="case|control"
        run:
            import os
            os.makedirs(os.path.dirname(output.collection), exist_ok=True)
            with open(output.collection, "w") as f:
                f.write("identifier\tfilename\tlabels\n")
                for sample in params.samples:
                    prefix = f"{params.output_dir}/{sample}"
                    f.write(f"{sample}\t{prefix}\t{wildcards.group}\n")


    # =========================================================================
    # Methbat_joint_segment_ASM: Joint segmentation for ASM discovery
    #
    # Performs joint segmentation separately for case and control groups
    # to identify regions with allele-specific methylation patterns.
    # Outputs: asm_{group}.meth_regions.bed, asm_{group}.asm.bedgraph, 
    #          asm_{group}.combined_methyl.bedgraph
    # =========================================================================
    rule Methbat_joint_segment_ASM:
        input:
            collection=METHBAT_BASE + "/asm_{group}_collection.tsv",
        output:
            regions=METHBAT_BASE + "/asm_{group}.meth_regions.bed",
            asm_bedgraph=METHBAT_BASE + "/asm_{group}.asm.bedgraph",
            combined_bedgraph=METHBAT_BASE + "/asm_{group}.combined_methyl.bedgraph",
        params:
            output_prefix=METHBAT_BASE + "/asm_{group}",
            threads=config["methbat"]["threads"],
        wildcard_constraints:
            group="case|control"
        log:
            config["directory"]["output"] + "/logs/methbat/joint_segment_asm_{group}.log",
        shell:
            """
            set +u
            source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
            conda activate methbat
            set -u
            mkdir -p $(dirname {params.output_prefix})
            mkdir -p $(dirname {log})
            methbat joint-segment \
                --threads {params.threads} \
                --input-collection {input.collection} \
                --output-prefix {params.output_prefix} \
                2>&1 | tee {log}
            conda deactivate
            """


    # =========================================================================
    # Merge_ASM_regions: Merge case and control ASM regions
    #
    # Creates a union of regions from both case and control joint-segmentation.
    # All samples will be profiled against this merged region set to ensure
    # consistent entries for methbat build.
    # =========================================================================
    rule Merge_ASM_regions:
        input:
            case_bed=METHBAT_BASE + "/asm_case.meth_regions.bed",
            control_bed=METHBAT_BASE + "/asm_control.meth_regions.bed",
        output:
            merged_tsv=METHBAT_BASE + "/asm_merged_regions.tsv",
        run:
            import os
            os.makedirs(os.path.dirname(output.merged_tsv), exist_ok=True)
            
            # Collect all unique regions
            regions = set()
            for bed_file in [input.case_bed, input.control_bed]:
                with open(bed_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        fields = line.split("\t")
                        chrom = fields[0]
                        start = int(fields[1])
                        end = int(fields[2])
                        label = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
                        # only include AlleleSpecificMethylation regions
                        if(label == "AlleleSpecificMethylation"):
                            regions.add((chrom, start, end, label))
            
            # Sort and write
            sorted_regions = sorted(regions, key=lambda x: (x[0], x[1], x[2]))
            with open(output.merged_tsv, "w") as f:
                f.write("chrom\tstart\tend\tcpg_label\n")
                for chrom, start, end, label in sorted_regions:
                    f.write(f"{chrom}\t{start}\t{end}\t{label}\n")
            
            print(f"Merged {len(sorted_regions)} unique regions from case and control")


    # =========================================================================
    # Methbat_profile_ASM: Profile methylation at ASM segments
    #
    # Generates methylation profiles at the merged ASM regions.
    # All samples use the same merged region set for consistent profiling.
    # Uses the original combined BED files which contain both haplotypes.
    # =========================================================================
    rule Methbat_profile_ASM:
        input:
            bed=config["directory"]["output"] + "/pb_cpg_tools/{sample}.combined.bed.gz",
            regions=METHBAT_BASE + "/asm_merged_regions.tsv",
        output:
            profile=METHBAT_BASE + "/profiles_asm/{sample}.profile.tsv",
        params:
            input_prefix=config["directory"]["output"] + "/pb_cpg_tools/{sample}",
        wildcard_constraints:
            sample="|".join([s.replace(".", r"\.") for s in config["samples"]["case"] + config["samples"]["control"]])
        log:
            config["directory"]["output"] + "/logs/methbat/profile_asm_{sample}.log",
        shell:
            """
            set +u
            source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
            conda activate methbat
            set -u
            mkdir -p $(dirname {output.profile})
            mkdir -p $(dirname {log})
            
            methbat profile \
                --input-prefix {params.input_prefix} \
                --input-regions {input.regions} \
                --output-region-profile {output.profile} \
                2>&1 | tee {log}
            conda deactivate
            """


    # =========================================================================
    # Create_ASM_profile_collection: Collection for ASM cohort analysis
    #
    # Creates a collection with case/control labels for group comparison.
    # Each profile contains hap1/hap2 data internally for ASM analysis.
    # =========================================================================
    rule Create_ASM_profile_collection:
        input:
            profiles=expand(
                METHBAT_BASE + "/profiles_asm/{sample}.profile.tsv",
                sample=ALL_SAMPLES
            ),
        output:
            collection=METHBAT_BASE + "/asm_profile_collection.tsv",
        params:
            profile_dir=METHBAT_BASE + "/profiles_asm",
            case_samples=config["samples"]["case"],
            control_samples=config["samples"]["control"],
        run:
            import os
            os.makedirs(os.path.dirname(output.collection), exist_ok=True)
            with open(output.collection, "w") as f:
                f.write("identifier\tfilename\tlabels\n")
                for sample in params.case_samples:
                    profile = f"{params.profile_dir}/{sample}.profile.tsv"
                    f.write(f"{sample}\t{profile}\tcase\n")
                for sample in params.control_samples:
                    profile = f"{params.profile_dir}/{sample}.profile.tsv"
                    f.write(f"{sample}\t{profile}\tcontrol\n")


    # =========================================================================
    # Methbat_build_ASM: Build cohort profile for ASM analysis
    # =========================================================================
    rule Methbat_build_ASM:
        input:
            collection=METHBAT_BASE + "/asm_profile_collection.tsv",
        output:
            cohort_profile=METHBAT_BASE + "/asm_cohort.profile.tsv",
        log:
            config["directory"]["output"] + "/logs/methbat/build_asm.log",
        shell:
            """
            set +u
            source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
            conda activate methbat
            set -u
            mkdir -p $(dirname {output.cohort_profile})
            mkdir -p $(dirname {log})
            methbat build \
                --input-collection {input.collection} \
                --output-profile {output.cohort_profile} \
                2>&1 | tee {log}
            conda deactivate
            """


    # =========================================================================
    # ASM_haplotype_comparison: Extract significant ASM regions from cohort
    #
    # ASM (hap1 vs hap2 difference) is computed during profiling and stored
    # in avg_abs_meth_deltas column. This rule extracts regions with
    # significant allele-specific methylation based on delta threshold.
    # =========================================================================
    rule ASM_haplotype_comparison:
        input:
            cohort_profile=METHBAT_BASE + "/asm_cohort.profile.tsv",
        output:
            comparison=METHBAT_BASE + "/asm_haplotype_comparison.tsv",
        params:
            min_delta=config["methbat"]["min_delta"],
        log:
            config["directory"]["output"] + "/logs/methbat/asm_haplotype_comparison.log",
        run:
            import pandas as pd
            import os
            
            os.makedirs(os.path.dirname(output.comparison), exist_ok=True)
            
            # Load cohort profile, skip comment lines
            df = pd.read_csv(input.cohort_profile, sep='\t', comment='#')
            print(f"Loaded {len(df)} entries from cohort profile")
            
            # Get ALL rows (aggregate across samples)
            all_df = df[df['data_category'] == 'ALL'].copy()
            print(f"Found {len(all_df)} aggregate (ALL) entries")
            
            # Filter by avg_abs_meth_deltas (hap1 vs hap2 difference)
            if 'avg_abs_meth_deltas' in all_df.columns:
                sig_asm = all_df[all_df['avg_abs_meth_deltas'] >= params.min_delta].copy()
                print(f"Found {len(sig_asm)} regions with avg_abs_meth_deltas >= {params.min_delta}")
                
                sig_asm.to_csv(output.comparison, sep='\t', index=False)
            else:
                print("Warning: avg_abs_meth_deltas column not found")
                all_df.head(0).to_csv(output.comparison, sep='\t', index=False)
            
            with open(log[0], 'w') as f:
                f.write(f"Loaded {len(df)} entries from cohort profile\n")
                f.write(f"Found {len(all_df)} aggregate (ALL) entries\n")
                if 'avg_abs_meth_deltas' in all_df.columns:
                    f.write(f"Found {len(sig_asm)} regions with avg_abs_meth_deltas >= {params.min_delta}\n")


    # =========================================================================
    # Methbat_compare_ASM_groups: Compare case vs control ASM regions
    #
    # Compares ASM regions between irradiated (case) and control groups
    # by computing overlaps between the joint-segmented regions.
    # Creates a confusion matrix showing overlaps by summary_label.
    # =========================================================================
    rule Methbat_compare_ASM_groups:
        input:
            case_bed=METHBAT_BASE + "/asm_case.meth_regions.bed",
            control_bed=METHBAT_BASE + "/asm_control.meth_regions.bed",
        output:
            comparison=METHBAT_BASE + "/asm_group_comparison.tsv",
            confusion_matrix=METHBAT_BASE + "/asm_group_confusion_matrix.tsv",
            case_specific=METHBAT_BASE + "/asm_case_specific.bed",
            control_specific=METHBAT_BASE + "/asm_control_specific.bed",
            shared=METHBAT_BASE + "/asm_shared.bed",
        params:
            script=os.path.join(SCRIPTS_DIR, "compare_asm_groups.py"),
            min_jaccard=config["methbat"].get("min_jaccard", 0.5),
        log:
            config["directory"]["output"] + "/logs/methbat/compare_asm_groups.log",
        shell:
            """
            module load bedtools/2.27.1
            mkdir -p $(dirname {output.comparison})
            mkdir -p $(dirname {log})
            python {params.script} \
                --case-bed {input.case_bed} \
                --control-bed {input.control_bed} \
                --output-comparison {output.comparison} \
                --output-confusion-matrix {output.confusion_matrix} \
                --output-case-specific {output.case_specific} \
                --output-control-specific {output.control_specific} \
                --output-shared {output.shared} \
                --min-jaccard {params.min_jaccard} \
                2>&1 | tee {log}
            """
