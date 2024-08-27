configfile: "config.yaml"

container: "envs/diadumenernaseq.sif"

rule report:
    input:
        correlation_heatmap=expand("data/processed/DESeq2/DESeq2_{level}/correlation_heatmap_{level}_combined.png", level=config["deseq2level"]),
        blastn_barplot=expand("data/processed/blastn/{transcriptome}_{db}_barplot.png", transcriptome=config["transcriptome"], db=config["blastndb"]), # For blast to work properly, I downloaded the nt_euk database into the data/processed/blastn/db folder, and also downloaded the taxdb into the same folder. For some reason, species search required that I move taxdb.bti and taxdb.btd to the working directory from which the blast command is run i.e. the data/processed/blastn folder
        volcano=expand("data/processed/DESeq2/DESeq2_{level}/{level}_volcano_plot.png", level=config["deseq2level"]),
        diff_heatmap=expand("data/processed/DESeq2/DESeq2_{level}/DESeq2_{level}_heatmap_TMM_combined.png", level=config["deseq2level"]),
        fastqc=expand("data/processed/fastqc/{sample}.merged.SE_trimmed_fastqc.html", sample=config["samples"]),
        dexseq=expand("data/processed/dexseq/{transcriptome}.dexseq.results.dat", transcriptome=config["transcriptome"]),


# Read in the assembled transcriptome
rule trinity_assembly:
    output:
        expand("data/raw/{transcriptome}.fasta", transcriptome=config["transcriptome"])

# Reduce transcript redundancy
rule cdhit:
    input:
        fasta=rules.trinity_assembly.output,
    output:
        expand("data/raw/{transcriptome}.cdhit.fasta", transcriptome=config["transcriptome"])
    threads: 30
    params:
        seq_id=config["cdhit_seq_id_threshold"],
        memory=config["cdhit_memory"],
        word_size=config["cd_hit_word_size"]
    log:
        logfile=expand("logs/cdhit/{transcriptome}.log", transcriptome=config["transcriptome"])
    shell:
        "cd-hit-est -i {input.fasta} -o {output} -c {params.seq_id} -n {params.word_size} -T {threads} -M {params.memory} > {log.logfile} 2>&1"

rule extract_orfs:
    input:
        fasta=rules.cdhit.output
    output:
        pep="data/processed/transdecoder/longest_orfs.pep"
    params:
        outputdir="data/processed/transdecoder"
    log:
        logfile="logs/extract_orfs/longest_orfs.log"
    shell:
        """
        mkdir -p {params.outputdir} &&

        cd {params.outputdir} &&

        TransDecoder.LongOrfs -t ../../../{input.fasta} > ../../../{log.logfile} 2>&1 &&

        mv $(basename ../../../{input.fasta}).transdecoder_dir/* . &&

        rm -r $(basename ../../../{input.fasta}).transdecoder_dir*
        """

rule get_ref_genome:
    output:
        accession=expand("data/references/{genome_accession}.fasta", genome_accession=config["genome_accession"]),
        genome=expand("data/references/{genome}.fasta", genome=config["genome"])
    params:
        genome=config["genome"],
        genome_accession=config["genome_accession"],
        outdir="data/references"
    shell:
        """
        mkdir -p {params.outdir} && 

        cd {params.outdir} && 

        datasets download genome accession {params.genome_accession} --include genome --no-progressbar && 
        
        unzip ncbi_dataset.zip && 
        
        mv ncbi_dataset/data/* . && 
        
        mv {params.genome_accession}/* . && 

        mv {params.genome_accession}*.fna {params.genome_accession}.fasta && 

        ln -s {params.genome_accession}.fasta {params.genome}.fasta && 
        
        rm -r ncbi_dataset* {params.genome_accession}/

        """


rule make_blastn_db_genome:
    input:
        fasta=expand("data/references/{genome}.fasta", genome=config["genome"]),
    output:
        expand("data/references/{genome}.fasta.njs", genome=config["genome"])
    params:
        db=expand("data/references/{genome}", genome=config["blastndb"])
    threads: 30
    log:
        logfile=expand("logs/{genome}_makedb.log", genome=config["genome"])
    shell:
        """

        mkdir -p data/references && 

        cd data/references && 

        makeblastdb -in ../../{input.fasta} -dbtype nucl -parse_seqids -taxid 1789172 -logfile ../../{log.logfile}

        """

# For rule below to work, user must provide the input files in the data/references folder

rule make_blastp_uniprot_db:
    input:
        fasta=expand("data/references/{blastpdb}/{blastpdb}_{blastpdb_taxid}.fasta", blastpdb=config["blastpdb"], blastpdb_taxid=config["blastpdb_taxid"]),
        taxonmap="data/references/prot.accession2taxid.FULL.gz",
        taxonnodes="data/references/nodes.dmp",
        taxonnames="data/references/names.dmp"
    output:
        expand("data/references/{blastpdb}.dmnd", blastpdb=config["blastpdb"])
    params:
        db=expand("data/references/{blastpdb}", blastpdb=config["blastpdb"])
    threads: 30
    log:
        logfile=expand("logs/{blastpdb}_{blastpdb_taxid}_makedb.log", blastpdb=config["blastpdb"], blastpdb_taxid=config["blastpdb_taxid"])
    shell:
        "diamond makedb --in {input.fasta} -d {params.db} --threads {threads} --taxonmap {input.taxonmap} --taxonnodes {input.taxonnodes} --taxonnames {input.taxonnames} > {log.logfile} 2>&1"


rule blastp_uniprot:
    input:
        pep=rules.extract_orfs.output.pep,
        dbfile=rules.make_blastp_uniprot_db.output
    output:
        blastout="data/processed/blastp/extract_orfs_longest_orfs.blastp.outfmt6"
    params:
        db=expand("data/references/{blastpdb}", blastpdb=config["blastpdb"]),
        evalue=config["evalue"]
    threads: 30
    log:
        logfile="logs/blastp_uniprot/extract_orfs2.log"
    shell:
        """

        mkdir -p tmp &&
        
        diamond blastp -d {params.db} -q {input.pep} --max-target-seqs 1 -e {params.evalue} --threads {threads} -t tmp -o {output.blastout} --outfmt 6 > {log.logfile} 2>&1 &&

        rm -r tmp
        
        """


rule salmon:
    input:
        fasta=rules.cdhit.output,
        samples="data/raw/salmon_samples.txt"
    output:
        quant=expand("data/processed/salmon/{transcriptome}.cdhit.fasta.gene_trans_map", transcriptome=config["transcriptome"]),
        multiqc=touch("data/processed/salmon.done")
    params:
        outputdir="data/processed/salmon",
        tmp_file=expand("data/raw/{transcriptome}.cdhit.fasta.gene_trans_map", transcriptome=config["transcriptome"])
    threads: 30
    log:
        logfile=expand("logs/salmon/{transcriptome}.log", transcriptome=config["transcriptome"])
    shell:
        """
        mkdir -p {params.outputdir} && 

        cd {params.outputdir} &&

        $TRINITY_HOME/util/align_and_estimate_abundance.pl --seqType fq --samples_file ../../../{input.samples} --transcripts ../../../{input.fasta} --est_method salmon --trinity_mode --prep_reference --thread_count {threads} > ../../../{log.logfile} 2>&1 &&

        mv ../../../{params.tmp_file} .

        """

rule counts_matrix:
    input:
        samples=rules.salmon.input.samples,
        gene_trans_map=rules.salmon.output
    output:
        # quants="data/processed/salmon/quant_files.list",
        counts="data/processed/salmon/trinity.{level}.counts.matrix",
        TPM="data/processed/salmon/trinity.{level}.TPM.not_cross_norm",
        TMM="data/processed/salmon/trinity.{level}.TMM.EXPR.matrix"
    params:
        salmondir="data/processed/salmon",
        prefix="trinity"
    log:
        logfile="logs/counts_matrix/{level}.log"
    shell:
        """
        cd {params.salmondir} &&

        find *_rep_* -name "quant.sf" > quant_files.list && 

        $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method salmon --out_prefix {params.prefix} --name_sample_by_basedir --quant_files quant_files.list --gene_trans_map ../../../{input.gene_trans_map} > ../../../{log.logfile} 2>&1

        """


rule DESeq2:
    input:
        counts="data/processed/salmon/trinity.{level}.counts.matrix",
        samples=rules.counts_matrix.input.samples
    output:
        deseq2="data/processed/DESeq2/trinity.{level}.counts.matrix.16_vs_rt.DESeq2.DE_results"
    params:
        deseq2dir="data/processed/DESeq2"
    log:
        logfile="logs/deseq2/{level}.log"
    shell:
        "$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.counts} --samples {input.samples} --method DESeq2 --output {params.deseq2dir} > {log.logfile} 2>&1"


rule extract_diff_expr_transcripts:
    input:
        matrix="data/processed/salmon/trinity.{level}.TMM.EXPR.matrix",
        samples="data/raw/salmon_samples.txt",
        transcripts=rules.DESeq2.output.deseq2
    output:
        # expand("data/processed/DESeq2_{level}/diffExpr.P{FDR}_C{foldchange}.matrix", FDR=config["FDR"], foldchange=config["foldchange"], allow_missing=True),
        transcripts="data/processed/DESeq2/DESeq2_{level}/differentially_expressed_{level}s.txt",
        normalized_heatmap="data/processed/DESeq2/DESeq2_{level}/DESeq2_{level}_heatmap_TMM_combined.png",
        vertical_heatmap="data/processed/DESeq2/DESeq2_{level}/DE_TPM+LogFold_{level}_heatmap_vertical.png",
        normalized_heatmap_pdf="data/processed/DESeq2/DESeq2_{level}/DESeq2_{level}_heatmap_TMM_combined.pdf"
    params:
        fdr=config["FDR"],
        foldchange=config["foldchange"],
        out_dir="data/processed/DESeq2/DESeq2_{level}",
        # lev="{level}"
    log:
        diffExprlog="logs/extract_diff_expr_{level}.log",
        # clusterslog=expand("log/define_clusters_2_{level}.log", level=config["deseq2level"], allow_missing=True),
        verticalHeatLog="logs/verticalHeatmap_{level}.log",
        normalizedHeatLog="logs/normalizedHeatmap_{level}.log"
    shell:
        """
        mkdir -p {params.out_dir} &&     

        cp {input.transcripts} {params.out_dir} &&  

        cd {params.out_dir} && 

        $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../../../../{input.matrix} --samples ../../../../{input.samples} -P {params.fdr} -C {params.foldchange} > ../../../../{log.diffExprlog} 2>&1 &&

        ../../../../scripts/verticalHeatmapDEPipeline.sh diffExpr.P{params.fdr}_C{params.foldchange}.matrix > ../../../../{log.verticalHeatLog} 2>&1 && 

        ../../../../scripts/makeNormalizedHeatmapDEPipeline.sh {params.fdr} {params.foldchange} > ../../../../{log.normalizedHeatLog} 2>&1 && 

        cd ../../../../ 

        """

rule add_gene_symbols:
    input:
        diffexpr=rules.extract_diff_expr_transcripts.output.transcripts,
        blastout=rules.blastp_uniprot.output
    output:
        withsymbols="data/processed/DESeq2/DESeq2_{level}/differentially_expressed_{level}s_with_symbols.txt"
    params:
        folder="data/processed/DESeq2/DESeq2_{level}",
        level="{level}"
    log:
        logfile="logs/add_gene_symbols/{level}.log"
    shell:
        """
        cd {params.folder} && 
        
        ../../../../scripts/add_gene_symbols.sh {params.level} {input.blastout} > ../../../../{log.logfile} 2>&1 && 

        cd ../../../../
        """


rule volcano_plot:
    input:
        dge_results=rules.DESeq2.output,
        gene_symbols=rules.add_gene_symbols.output,
        trinity_genes="data/processed/DESeq2/DESeq2_{level}/differentially_expressed_{level}s.txt"
    output:
        volcano_plot="data/processed/DESeq2/DESeq2_{level}/{level}_volcano_plot.png",
        enrichr="data/processed/DESeq2/DESeq2_{level}/{level}_enrichr_kegg_2021_human.png",
        enrichr_table="data/processed/DESeq2/DESeq2_{level}/{level}_enrichr_kegg_2021_human.txt",
        volcano_pdf="data/processed/DESeq2/DESeq2_{level}/{level}_volcano_plot.pdf",
        enrichr_pdf="data/processed/DESeq2/DESeq2_{level}/{level}_enrichr_kegg_2021_human.pdf"
    log:
        logfile="logs/volcano_plot_{level}.log"
    script:
        "scripts/volcanoplot.R"


rule correlate:
    input:
        counts="data/processed/salmon/trinity.{level}.counts.matrix",
        samples="data/raw/salmon_samples.txt"
    output:
        heatmap="data/processed/DESeq2/DESeq2_{level}/correlation_heatmap_{level}_combined.png"
    log:
        logfile="logs/correlate_{level}.log"
    shell:
        """
        cd data/processed/DESeq2/DESeq2_{wildcards.level} && 

        bash ../../../../scripts/correlateReplicates.sh {input.samples} > ../../../../{log.logfile} 2>&1 

        """

rule blastn_all_db_setup:
    container: None
    output:
        directory("data/processed/blastn/db")
    params:
        sys_db_path=config["blastndb_system_path"],
        local_db_path="data/processed/blastn",
    shell:
        """
        mkdir -p {params.local_db_path} &&

        mkdir -p {params.local_db_path}/db &&

        cd {params.local_db_path}/db &&

        cp -u {params.sys_db_path}/* .

        """


rule blastn_all:
    input:
        query=rules.cdhit.output,
        db_dir=rules.blastn_all_db_setup.output
    output:
        expand("data/processed/blastn/{transcriptome}.outfmt7", transcriptome=config["transcriptome"]),
    params:
        qcov=config["qcov_hsp_perc"],
        db=config["blastndb"], # This is the actual blast db
        outdir="data/processed/blastn"
    threads: 30
    log:
        logfile=expand("logs/blastn/{transcriptome}.log", transcriptome=config["transcriptome"])
    shell:
        """

        cd {params.outdir} &&

        # If taxdb.bti and taxdb.btd files are not present in the current directory, mv them from the db directory

        if [ ! -f taxdb.bti ] || [ ! -f taxdb.btd ]; then

            mv db/taxdb.bti . &&

            mv db/taxdb.btd .

        fi &&

        blastdbcmd -info -db db/{params.db} > ../../../{log.logfile} 2>&1 &&

        blastn -query ../../../{input.query} -db db/{params.db} -out ../../../{output} -qcov_hsp_perc {params.qcov} -outfmt "7 saccver evalue pident qcovs staxids sscinames scomnames sblastnames" -num_threads {threads} >> ../../../{log.logfile} 2>&1

        """

rule blastn_all_data_extraction:
    input:
        blastn=rules.blastn_all.output,
        cdhit=rules.cdhit.output
    output:
        extract="data/processed/blastn/{transcriptome}_{db}_extracted.tsv",
        barplot_png="data/processed/blastn/{transcriptome}_{db}_barplot.png",
        barplot_pdf="data/processed/blastn/{transcriptome}_{db}_barplot.pdf"
    params:
        barplot_file_name="data/processed/blastn/{transcriptome}_{db}_barplot"
    log:
        logfile="logs/blastn/{transcriptome}_{db}_extracted.log"
    shell:
        """
        scripts/extract_blastn_all.sh {input.blastn} {input.cdhit} {output.extract} {params.barplot_file_name} > {log.logfile} 2>&1

        """

rule make_blastx_nr_db:
    input:
        fasta=expand("data/references/{blastxdb}.fasta", blastxdb=config["blastxdb"]),
        taxonmap="data/references/prot.accession2taxid.FULL.gz",
        taxonnodes="data/references/nodes.dmp",
        taxonnames="data/references/names.dmp"
    output:
        expand("data/references/{blastxdb}.dmnd", blastxdb=config["blastxdb"])
    params:
        db=expand("data/references/{blastxdb}", blastxdb=config["blastxdb"])
    threads: 30
    log:
        logfile=expand("logs/{blastxdb}_makedb.log", blastxdb=config["blastxdb"])
    shell:
        "diamond makedb --in {input.fasta} -d {params.db} --threads {threads} --taxonmap {input.taxonmap} --taxonnodes {input.taxonnodes} --taxonnames {input.taxonnames} > {log.logfile} 2>&1"

rule supertranscripts:
    input:
        fasta=rules.cdhit.output,
    output:
        super_fasta="data/processed/dexseq/{transcriptome}.fasta",
        super_gtf="data/processed/dexseq/{transcriptome}.gtf"
    params:
        out_dir="data/processed/dexseq/",
        prefix="{transcriptome}"
    threads: 30
    log:
        log="logs/supertranscript_{transcriptome}.log"
    shell:
        """
        mkdir -p {params.out_dir} &&  
        
        cd {params.out_dir} &&

        $TRINITY_HOME/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
        --trinity_fasta ../../../{input.fasta} \
        --out_prefix {params.prefix} \
        > ../../../{log.log} 2>&1
        """



# Run FASTQC on the raw data
rule fastqc:
    input:
        "data/raw/{sample}.merged.SE_trimmed.fq.gz"
    output:
        html="data/processed/fastqc/{sample}.merged.SE_trimmed_fastqc.html",
        zip="data/processed/fastqc/{sample}.merged.SE_trimmed_fastqc.zip", # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
        # multiqc = touch("data/processed/{sample}fastqc.done")
    params:
        extra = "--quiet",
        outdir = "data/processed/fastqc"
    log:
        "logs/fastqc/{sample}.log"
    threads: 30
    resources:
        mem_mb = 16000
    shell:
        "fastqc --threads {threads} --memory {resources.mem_mb} --outdir {params.outdir} {input} > {log} 2>&1"


rule dexseq:
    input:
        fasta=rules.supertranscripts.output.super_fasta,
        gtf=rules.supertranscripts.output.super_gtf,
        samples="data/raw/salmon_samples.txt"
    output:
        dex_out="data/processed/dexseq/{transcriptome}.dexseq.results.dat",
    params:
        out_dir="data/processed/dexseq/",
        prefix="{transcriptome}",
        image="envs/trinity.simg"
    threads: 30
    log:
        "logs/dexseq_{transcriptome}.log" 
    shell:
        """
        cd {params.out_dir} &&

        ../../../scripts/dexseq_wrapper.pl \
        --genes_fasta ../../../{input.fasta} \
        --genes_gtf ../../../{input.gtf} \
        --samples_file ../../../{input.samples} \
        --out_prefix {params.prefix} \
        --aligner STAR \
        --CPU {threads} \
        > ../../../{log} 2>&1

        """

