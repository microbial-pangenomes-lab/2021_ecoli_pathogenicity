from snakemake.utils import validate
import pandas as pd
import os

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

##### functions #####

def _read_samples(infile):
  m = pd.read_csv(infile, sep='\t', index_col=0)
  return set(m.index)

##### rules #####

rule prepare_clinical:
  input: config["samples"]
  output:
    config["clinical"],
    config["clinical_reduced"],
    config["clinical_full"]
  shell:
    "python workflow/scripts/prepare_clinical_data.py {input} {output}"

rule impute_clinical:
  input: config["clinical_full"]
  output: config["clinical_imputed"]
  conda: "../envs/mice.yaml"
  shell:
    "Rscript workflow/scripts/impute.R {input} {output}"

rule poppunk:
  input: config["poppunk_input"]
  output: config["poppunk"]
  params:
    indir=config["poppunk_db"],
    outdir=config["poppunk_dir"]
  threads: 32
  conda: "../envs/poppunk.yaml"
  log: "out/logs/poppunk.log"
  shell:
     """
     poppunk --assign-query --ref-db {params.indir} \
             --q-files {input} --output {params.outdir} \
             --threads {threads} && \
     tail -n+2 {params.outdir}/*_clusters.csv | \
     sed 's/data\/fastas\///g' | sed 's/.fasta//g' | \
     awk -F',' '{{print $1"\\t"$2}}' > {output}
     """

rule prepare_lineages:
  input: config["typing"]
  output: config["lineage"],
  shell:
    "tail -n+2 {input} | awk -F',' '{{print $1\"\\t\"$2}}' > {output}"

rule mash_sketch:
  input: config["mash_input"]
  output: config["sketches"]
  params: config["sketches_base"]
  threads: 5
  conda: "../envs/mash.yaml"
  log: "out/logs/mash_sketch.log"
  shell:
    "mash sketch -p {threads} -s 10000 -o {params} -l {input}"

rule distance:
  input: config["sketches"]
  output: config["distances"]
  threads: 5
  conda: "../envs/mash.yaml"
  log: "out/logs/distance.log"
  shell:
    "mash dist -p {threads} {input} {input} | square_mash > {output}"

rule similarity:
  input: config["tree"]
  output: config["similarities"]
  conda: "../envs/pyseer.yaml"
  log: "out/logs/similarity.log"
  shell:
    "python3 workflow/scripts/phylogeny_distance.py --calc-C {input} > {output}"

rule pangenome:
  input: config["panaroo_input"]
  output:
    config["pangenome"],
    config["pangenome_roary"],
    config["structural"]
  params: config["panaroo_dir"]
  threads: 40
  conda: "../envs/panaroo.yaml"
  log: "out/logs/pangenome.log"
  shell:
    "panaroo -t {threads} -i {input} -o {params} --clean-mode strict"

rule pangenome_ref:
  input: config["panaroo_ref_input"]
  output:
    config["pangenome_ref"],
    config["pangenome_ref_roary"]
  params: config["panaroo_ref_dir"]
  threads: 40
  conda: "../envs/panaroo.yaml"
  log: "out/logs/pangenome_ref.log"
  shell:
    "panaroo -t {threads} -i {input} -o {params} --clean-mode strict"

rule unitigs:
  input: config["unitigs_input"]
  output: 
    config["unitigs"],
    config["unitigs_rtab"]
  params: config["unitigs_dir"]
  threads: 40
  conda: "../envs/unitig-counter.yaml"
  log: "out/logs/unitigs.log"
  shell:
    "rm -rf {params} && unitig-counter -strains {input} -output {params} -nb-cores {threads}"

rule univariate:
  input: config["clinical_imputed"]
  output: config["univariate"]
  params: config["univariate_dir"]
  log: "out/logs/univariate.log"
  conda: "../envs/pyseer.yaml"
  shell:
    "python3 workflow/scripts/univariate_analysis.py {input} {params} > {output}"

rule multivariate:
  input: config["univariate"]
  output: os.path.join(config["univariate_dir"], 'multivariate.done')
  params: config["univariate_dir"]
  log: "out/logs/multivariate.log"
  conda: "../envs/r-mass.yaml"
  shell:
    """
    Rscript workflow/scripts/multivariate_analysis.R {params}/deces.tsv {params}/deces.multi.tsv
    Rscript workflow/scripts/multivariate_analysis.R {params}/choc.tsv {params}/choc.multi.tsv
    Rscript workflow/scripts/multivariate_analysis.R {params}/passage_en_rea.tsv {params}/passage_en_rea.multi.tsv
    Rscript workflow/scripts/multivariate_analysis.R {params}/deces_septicoli.tsv {params}/deces_septicoli.multi.tsv
    Rscript workflow/scripts/multivariate_analysis.R {params}/choc_septicoli.tsv {params}/choc_septicoli.multi.tsv
    Rscript workflow/scripts/multivariate_analysis.R {params}/passage_en_rea_septicoli.tsv {params}/passage_en_rea_septicoli.multi.tsv
    Rscript workflow/scripts/multivariate_analysis.R {params}/deces_colibafi.tsv {params}/deces_colibafi.multi.tsv
    Rscript workflow/scripts/multivariate_analysis.R {params}/choc_colibafi.tsv {params}/choc_colibafi.multi.tsv
    Rscript workflow/scripts/multivariate_analysis.R {params}/passage_en_rea_colibafi.tsv {params}/passage_en_rea_colibafi.multi.tsv
    touch {output}
    """

rule prepare_pyseer:
  input:
    variants=config["unitigs_input"],
    phenotypes=config["clinical_imputed"],
    similarity=config["similarities"],
    distances=config["distances"],
    lineage=config["lineage"]
  output:
    phenotypes=os.path.join(config["association_inputs"], "phenotypes.tsv"),
    similarity=os.path.join(config["association_inputs"], "similarity.tsv"),
    distances=os.path.join(config["association_inputs"], "distances.tsv"),
    lineage=os.path.join(config["association_inputs"], "lineages.tsv")
  params: config["association_inputs"]
  log: "out/logs/prepare_pyseer.log"
  shell:
    "python3 workflow/scripts/prepare_pyseer.py {input} {params}"

rule sift:
  input: config["snps_reference_faa"]
  output: "out/snps/sift/sift.done"
  params:
    db=config["uniref50"],
    outdir=config["sift"]
  threads: 36
  conda: "../envs/sift4g.yaml"
  log: "out/logs/sift.log"
  shell:
    """
    sift4g -q {input} -d {params.db} --out {params.outdir} \
           -t {threads} --sub-results \
           --median-threshold 3.25
    Rscript workflow/scripts/install_siftr.R || true
    for i in $(ls {params.outdir} | grep aligned.fasta); \
    do \
      echo $i; \
      Rscript workflow/scripts/sift_scorer.R {params.outdir}/$i {params.outdir}/$(basename $i .aligned.fasta).tsv || true; \
    done
    touch {output}
    """

rule prepare_rare_variants:
  input:
    vcf=expand("out/snps/{sample}/snps.vcf.gz",
               sample=_read_samples(config["unitigs_input"])),
    sift="out/snps/sift/sift.done"
  output: config["rare_snps"]
  params: config["sift"]
  threads: 36
  conda: "../envs/bcftools.yaml"
  log: "out/logs/prepare_rare_variants.log"
  shell:
    """
    bcftools merge {input.vcf} -0 -O z --threads {threads} > out/snps/merged.vcf.gz
    bcftools norm -m - -O z --threads {threads} out/snps/merged.vcf.gz > out/snps/norm.vcf.gz
    bcftools view -Q 0.01 out/snps/norm.vcf.gz -O z --threads {threads} > out/snps/filtered.vcf.gz
    python workflow/scripts/vcf2deleterious.py out/snps/filtered.vcf.gz {params} | bgzip > {output}
    bcftools index {output}
    rm out/snps/merged.vcf.gz out/snps/norm.vcf.gz out/snps/filtered.vcf.gz
    """

rule get_snps:
  input:
    ref=config["snps_reference"],
    ctgs=os.path.join(config["fixed_fastas"], "{sample}.fasta")
  output:
    "out/snps/{sample}/snps.vcf.gz"
  params:
    "out/snps/{sample}"
  conda: "../envs/nucmer.yaml"
  log: "out/logs/snippy_{sample}.log"
  shell:
    "snippy --force --outdir {params} --ref {input.ref} --ctgs {input.ctgs} --cpus 1 --ram 8"

rule prepare_regions:
  input: config["snps_reference_gff"],
  output: config["regions"]
  log: "out/logs/prepare_regions.log"
  shell:
    """
    grep CDS {input} | \
    awk -F '\t' '{{print $1":"$4"-"$5"\t"$9}}' | \
    sed 's/.1:/:/g' | awk -F ';' '{{print $1}}' | \
    sed 's/ID=//g' | awk '{{print $2"\t"$1}}' > {output} 
    """

rule heritability:
  input:
    expand('out/associations/{target}/heritability.tsv',
           target=config["targets"]),
    expand('out/associations/{target}/heritability_lineages.tsv',
           target=config["targets"]),
    expand('out/associations/{target}/heritability_septicoli.tsv',
           target=config["targets"]),
    expand('out/associations/{target}/heritability_septicoli_lineages.tsv',
           target=config["targets"]),
    expand('out/associations/{target}/heritability_colibafi.tsv',
           target=config["targets"]),
    expand('out/associations/{target}/heritability_colibafi_lineages.tsv',
           target=config["targets"]),
    expand('out/associations/{target}/heritability.ci.tsv',
           target=config["targets"]),

rule unitigs2covariance:
  input:
    config["unitigs_rtab"]
  output:
    os.path.join(config["association_inputs"], "unitigs_covariance.tsv")
  log: "out/logs/unitigs2covariance.log"
  shell:
    "python workflow/scripts/unitigs2covar.py {input} > {output}"

rule lineages2covariance:
  input:
    os.path.join(config["association_inputs"], "lineages.tsv")
  output:
    os.path.join(config["association_inputs"], "lineages_covariance.tsv")
  log: "out/logs/lineages2covariance.log"
  shell:
    "python workflow/scripts/lineage2covar.py {input} > {output}"

rule run_heritability:
  input:
    phenotypes=os.path.join(config["association_inputs"], "phenotypes.tsv"),
    similarity=os.path.join(config["association_inputs"], "similarity.tsv"),
    lineages=os.path.join(config["association_inputs"], "lineages_covariance.tsv"),
  output:
    h_ci="out/associations/{target}/heritability.ci.tsv",
    h="out/associations/{target}/heritability.tsv",
    h_lineages="out/associations/{target}/heritability_lineages.tsv",
    h_septicoli="out/associations/{target}/heritability_septicoli.tsv",
    h_septicoli_l="out/associations/{target}/heritability_septicoli_lineages.tsv",
    h_colibafi="out/associations/{target}/heritability_colibafi.tsv",
    h_colibafi_l="out/associations/{target}/heritability_colibafi_lineages.tsv",
  params:
    covariates=lambda wildcards: config["covariates"][wildcards.target]    
  conda: "../envs/limix.yaml"
  log: "out/logs/heritability_{target}.log"
  shell:
    """
    python workflow/scripts/estimate_heritability.py {input.phenotypes} {input.similarity} \
          -p {wildcards.target} \
          --use-covariates {params.covariates} \
    > {output.h}
    python workflow/scripts/prepare_fiesta.py {input.phenotypes} {input.similarity} \
	  -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --prefix /tmp/{wildcards.target}
    grep normal {output.h} | awk '{{print $3}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h} | awk '{{print $3}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimating" | grep -v '#' > {output.h_ci}
    grep normal {output.h} | awk '{{print $4}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h} | awk '{{print $4}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt -x /tmp/{wildcards.target}_covariates.txt \
          -v /tmp/{wildcards.target}_vectors.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
 
    python workflow/scripts/estimate_heritability.py {input.phenotypes} {input.lineages} \
          -p {wildcards.target} \
          --use-covariates {params.covariates} \
    > {output.h_lineages} 
    python workflow/scripts/prepare_fiesta.py {input.phenotypes} {input.lineages} \
	  -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --prefix /tmp/{wildcards.target}
    grep normal {output.h_lineages} | awk '{{print $3}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_lineages} | awk '{{print $3}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    grep normal {output.h_lineages} | awk '{{print $4}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_lineages} | awk '{{print $4}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt -x /tmp/{wildcards.target}_covariates.txt \
          -v /tmp/{wildcards.target}_vectors.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    
    python workflow/scripts/estimate_heritability.py {input.phenotypes} {input.similarity} \
          -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --study septicoli \
    > {output.h_septicoli} 
    python workflow/scripts/prepare_fiesta.py {input.phenotypes} {input.similarity} \
	  -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --prefix /tmp/{wildcards.target}
    grep normal {output.h_septicoli} | awk '{{print $3}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_septicoli} | awk '{{print $4}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    grep normal {output.h_septicoli} | awk '{{print $4}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_septicoli} | awk '{{print $4}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt -x /tmp/{wildcards.target}_covariates.txt \
          -v /tmp/{wildcards.target}_vectors.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    
    python workflow/scripts/estimate_heritability.py {input.phenotypes} {input.lineages} \
          -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --study septicoli \
    > {output.h_septicoli_l} 
    python workflow/scripts/prepare_fiesta.py {input.phenotypes} {input.lineages} \
	  -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --prefix /tmp/{wildcards.target}
    grep normal {output.h_septicoli_l} | awk '{{print $3}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_septicoli_l} | awk '{{print $3}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    grep normal {output.h_septicoli_l} | awk '{{print $4}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_septicoli_l} | awk '{{print $4}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt -x /tmp/{wildcards.target}_covariates.txt \
          -v /tmp/{wildcards.target}_vectors.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    
    python workflow/scripts/estimate_heritability.py {input.phenotypes} {input.similarity} \
          -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --study colibafi \
    > {output.h_colibafi}
    python workflow/scripts/prepare_fiesta.py {input.phenotypes} {input.similarity} \
	  -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --prefix /tmp/{wildcards.target}
    grep normal {output.h_colibafi} | awk '{{print $3}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_colibafi} | awk '{{print $3}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    grep normal {output.h_colibafi} | awk '{{print $4}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_colibafi} | awk '{{print $4}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt -x /tmp/{wildcards.target}_covariates.txt \
          -v /tmp/{wildcards.target}_vectors.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    
    python workflow/scripts/estimate_heritability.py {input.phenotypes} {input.lineages} \
          -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --study colibafi \
    > {output.h_colibafi_l}
    python workflow/scripts/prepare_fiesta.py {input.phenotypes} {input.lineages} \
	  -p {wildcards.target} \
          --use-covariates {params.covariates} \
          --prefix /tmp/{wildcards.target}
    grep normal {output.h_colibafi_l} | awk '{{print $3}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_colibafi_l} | awk '{{print $3}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    grep normal {output.h_colibafi_l} | awk '{{print $4}}' > /tmp/{wildcards.target}.estimates.txt
    grep normal {output.h_colibafi_l} | awk '{{print $4}}' >> /tmp/{wildcards.target}.estimates.txt
    python albi/albi.py -k /tmp/{wildcards.target}_values.txt -x /tmp/{wildcards.target}_covariates.txt \
          -v /tmp/{wildcards.target}_vectors.txt \
          -f /tmp/{wildcards.target}.estimates.txt \
    | grep -v "Estimat" | grep -v '#' >> {output.h_ci}
    """

rule pyseer:
  input:
    expand('out/associations/{target}/unitigs_filtered.tsv',
           target=config["targets"])

rule pyseer_rare:
  input:
    expand('out/associations/{target}/rare_filtered.tsv',
           target=config["targets"])

rule pyseer_nc:
  input:
    expand('out/associations/{target}/nc_unitigs_filtered.tsv',
           target=config["targets"])

rule run_pyseer:
  input:
    unitigs=config["unitigs"],
    gpa=config["pangenome"],
    struct=config["structural"],
    phenotypes=os.path.join(config["association_inputs"], "phenotypes.tsv"),
    similarity=os.path.join(config["association_inputs"], "similarity.tsv"),
    distances=os.path.join(config["association_inputs"], "distances.tsv"),
    lineages=os.path.join(config["association_inputs"], "lineages.tsv")
  output:
    unitigs="out/associations/{target}/unitigs.tsv",
    unitigs_f="out/associations/{target}/unitigs_filtered.tsv",
    gpa="out/associations/{target}/gpa.tsv",
    gpa_f="out/associations/{target}/gpa_filtered.tsv",
    struct="out/associations/{target}/struct.tsv",
    struct_f="out/associations/{target}/struct_filtered.tsv"
  params:
    covariates=lambda wildcards: config["covariates"][wildcards.target]    
  threads: 2
  conda: "../envs/pyseer.yaml"
  log: "out/logs/pyseer_{target}.log"
  shell:
    """
    pyseer --phenotypes {input.phenotypes} \
           --phenotype-column {wildcards.target} \
           --kmers {input.unitigs} \
           --similarity {input.similarity} \
           --lmm --uncompressed \
           --covariates {input.phenotypes} \
           --use-covariates {params.covariates} \
           --output-patterns out/associations/{wildcards.target}/unitigs_patterns.txt \
           --cpu {threads} \
           --lineage --lineage-clusters {input.lineages} \
           --lineage-file out/associations/{wildcards.target}/unitigs_lineage.txt \
           --distances {input.distances} \
           > {output.unitigs} && \
    cat <(head -1 {output.unitigs}) <(LC_ALL=C awk -v pval=$(python workflow/scripts/count_patterns.py --threshold out/associations/{wildcards.target}/unitigs_patterns.txt) '$4<pval {{print $0}}' {output.unitigs}) > {output.unitigs_f}
    pyseer --phenotypes {input.phenotypes} \
           --phenotype-column {wildcards.target} \
           --pres {input.gpa} \
           --similarity {input.similarity} \
           --lmm --uncompressed \
           --covariates {input.phenotypes} \
           --use-covariates {params.covariates} \
           --output-patterns out/associations/{wildcards.target}/gpa_patterns.txt \
           --cpu {threads} \
           --lineage --lineage-clusters {input.lineages} \
           --lineage-file out/associations/{wildcards.target}/gpa_lineage.txt \
           --distances {input.distances} \
           > {output.gpa} && \
    cat <(head -1 {output.gpa}) <(LC_ALL=C awk -v pval=$(python workflow/scripts/count_patterns.py --threshold out/associations/{wildcards.target}/gpa_patterns.txt) '$4<pval {{print $0}}' {output.gpa}) > {output.gpa_f}
    pyseer --phenotypes {input.phenotypes} \
           --phenotype-column {wildcards.target} \
           --pres {input.struct} \
           --similarity {input.similarity} \
           --lmm --uncompressed \
           --covariates {input.phenotypes} \
           --use-covariates {params.covariates} \
           --output-patterns out/associations/{wildcards.target}/struct_patterns.txt \
           --cpu {threads} \
           --lineage --lineage-clusters {input.lineages} \
           --lineage-file out/associations/{wildcards.target}/struct_lineage.txt \
           --distances {input.distances} \
           > {output.struct} && \
    cat <(head -1 {output.struct}) <(LC_ALL=C awk -v pval=$(python workflow/scripts/count_patterns.py --threshold out/associations/{wildcards.target}/struct_patterns.txt) '$4<pval {{print $0}}' {output.struct}) > {output.struct_f}
    """

rule run_pyseer_rare:
  input:
    snps=config["rare_snps"],
    regions=config["regions"],
    phenotypes=os.path.join(config["association_inputs"], "phenotypes.tsv"),
    similarity=os.path.join(config["association_inputs"], "similarity.tsv"),
    lineages=os.path.join(config["association_inputs"], "lineages.tsv")
  output:
    rare="out/associations/{target}/rare.tsv",
    rare_f="out/associations/{target}/rare_filtered.tsv",
  params:
    covariates=lambda wildcards: config["covariates"][wildcards.target]    
  conda: "../envs/pyseer.yaml"
  log: "out/logs/pyseer_rare_{target}.log"
  shell:
    """
    pyseer --phenotypes {input.phenotypes} \
           --phenotype-column {wildcards.target} \
           --vcf {input.snps} \
           --burden {input.regions} \
           --similarity {input.similarity} \
           --lmm \
           --covariates {input.phenotypes} \
           --use-covariates {params.covariates} \
           --output-patterns out/associations/{wildcards.target}/rare_patterns.txt \
           --cpu 1 \
           > {output.rare} && \
    cat <(head -1 {output.rare}) <(LC_ALL=C awk -v pval=$(python workflow/scripts/count_patterns.py --threshold out/associations/{wildcards.target}/rare_patterns.txt) '$4<pval {{print $0}}' {output.rare}) > {output.rare_f}
    """

rule run_pyseer_nc:
  input:
    unitigs=config["unitigs"],
    gpa=config["pangenome"],
    struct=config["structural"],
    phenotypes=os.path.join(config["association_inputs"], "phenotypes.tsv"),
    similarity=os.path.join(config["association_inputs"], "similarity.tsv"),
    distances=os.path.join(config["association_inputs"], "distances.tsv"),
    lineages=os.path.join(config["association_inputs"], "lineages.tsv")
  output:
    unitigs="out/associations/{target}/nc_unitigs.tsv",
    unitigs_f="out/associations/{target}/nc_unitigs_filtered.tsv",
    gpa="out/associations/{target}/nc_gpa.tsv",
    gpa_f="out/associations/{target}/nc_gpa_filtered.tsv",
    struct="out/associations/{target}/nc_struct.tsv",
    struct_f="out/associations/{target}/nc_struct_filtered.tsv"
  threads: 2
  conda: "../envs/pyseer.yaml"
  log: "out/logs/pyseer_nc_{target}.log"
  shell:
    """
    pyseer --phenotypes {input.phenotypes} \
           --phenotype-column {wildcards.target} \
           --kmers {input.unitigs} \
           --similarity {input.similarity} \
           --lmm --uncompressed \
           --output-patterns out/associations/{wildcards.target}/nc_unitigs_patterns.txt \
           --cpu {threads} \
           > {output.unitigs} && \
    cat <(head -1 {output.unitigs}) <(LC_ALL=C awk -v pval=$(python workflow/scripts/count_patterns.py --threshold out/associations/{wildcards.target}/nc_unitigs_patterns.txt) '$4<pval {{print $0}}' {output.unitigs}) > {output.unitigs_f}
    pyseer --phenotypes {input.phenotypes} \
           --phenotype-column {wildcards.target} \
           --pres {input.gpa} \
           --similarity {input.similarity} \
           --lmm --uncompressed \
           --lineage --lineage-clusters {input.lineages} \
           --lineage-file out/associations/{wildcards.target}/nc_lineage.txt \
           --distances {input.distances} \
           --output-patterns out/associations/{wildcards.target}/nc_gpa_patterns.txt \
           --cpu {threads} \
           > {output.gpa} && \
    cat <(head -1 {output.gpa}) <(LC_ALL=C awk -v pval=$(python workflow/scripts/count_patterns.py --threshold out/associations/{wildcards.target}/nc_gpa_patterns.txt) '$4<pval {{print $0}}' {output.gpa}) > {output.gpa_f}
    pyseer --phenotypes {input.phenotypes} \
           --phenotype-column {wildcards.target} \
           --pres {input.struct} \
           --similarity {input.similarity} \
           --lmm --uncompressed \
           --output-patterns out/associations/{wildcards.target}/nc_struct_patterns.txt \
           --cpu {threads} \
           > {output.struct} && \
    cat <(head -1 {output.struct}) <(LC_ALL=C awk -v pval=$(python workflow/scripts/count_patterns.py --threshold out/associations/{wildcards.target}/nc_struct_patterns.txt) '$4<pval {{print $0}}' {output.struct}) > {output.struct_f}
    """

rule map_back:
  input:
    expand("out/associations/{target}/mapped/{sample}.txt",
           target=config["targets"],
           sample=_read_samples(config["unitigs_input"]))
  output:
    "out/associations/mapped.done"
  shell:
    "touch {output}"

rule map_back_nc:
  input:
    expand("out/associations/{target}/nc_mapped/{sample}.txt",
           target=config["targets"],
           sample=_read_samples(config["unitigs_input"]))
  output:
    "out/associations/mapped_nc.done"
  shell:
    "touch {output}"

rule run_map_back:
  input:
    passing="out/associations/{target}/unitigs_filtered.tsv",
    fasta=os.path.join(config["fixed_fastas"], "{sample}.fasta"),
    gff=os.path.join(config["gffs"], "{sample}.gff"),
    pangenome=config["pangenome_roary"]
  output:
    "out/associations/{target}/mapped/{sample}.txt"
  params:
    "/tmp/map_back_{target}_{sample}"
  conda: "../envs/pyseer.yaml"
  log: "out/logs/map_back_{target}_{sample}.log"
  shell:
    """
    mkdir -p {params} && \
    python workflow/scripts/map_back.py {input.passing} {input.fasta} --tmp-prefix {params} --gff {input.gff} --print-details --roary {input.pangenome} > {output}
    """

rule run_map_back_nc:
  input:
    passing="out/associations/{target}/nc_unitigs_filtered.tsv",
    fasta=os.path.join(config["fixed_fastas"], "{sample}.fasta"),
    gff=os.path.join(config["gffs"], "{sample}.gff"),
    pangenome=config["pangenome_roary"]
  output:
    "out/associations/{target}/nc_mapped/{sample}.txt"
  params:
    "/tmp/map_back_nc_{target}_{sample}"
  conda: "../envs/pyseer.yaml"
  log: "out/logs/map_back_nc_{target}_{sample}.log"
  shell:
    """
    mkdir -p {params} && \
    python workflow/scripts/map_back.py {input.passing} {input.fasta} --tmp-prefix {params} --gff {input.gff} --print-details --roary {input.pangenome} > {output}
    """

rule map_summary:
  input:
    mapped=expand("out/associations/{target}/mapped.tsv",
                  target=config["targets"]),
    summary=expand("out/associations/{target}/summary.tsv",
                   target=config["targets"])

rule map_summary_nc:
  input:
    mapped=expand("out/associations/{target}/nc_mapped.tsv",
                  target=config["targets"]),
    summary=expand("out/associations/{target}/nc_summary.tsv",
                   target=config["targets"])

rule run_map_summary:
  input:
    phenotypes=os.path.join(config["association_inputs"], "phenotypes.tsv"),
    filtered="out/associations/{target}/unitigs_filtered.tsv",
    pangenome=config["pangenome"],
    mapped="out/associations/mapped.done"
  output:
    mapped="out/associations/{target}/mapped.tsv",
    summary="out/associations/{target}/summary.tsv"
  log: "out/logs/map_summary_{target}.log"
  shell:
    """
    echo -e "strain\\tunitig\\tcontig\\tstart\\tend\\tstrand\\tupstream\\tgene\\tdownstream" > {output.mapped}
    cat out/associations/{wildcards.target}/mapped/*.txt >> {output.mapped}
    python workflow/scripts/mapped_summary.py {output.mapped} \
           {input.phenotypes} {wildcards.target} {input.filtered} \
           --pangenome {input.pangenome} \
           --length 30 --minimum-hits 9 --maximum-genes 10 \
           --unique \
           > {output.summary}
    """

rule run_map_summary_nc:
  input:
    phenotypes=os.path.join(config["association_inputs"], "phenotypes.tsv"),
    filtered="out/associations/{target}/nc_unitigs_filtered.tsv",
    pangenome=config["pangenome"],
    mapped="out/associations/mapped_nc.done"
  output:
    mapped="out/associations/{target}/nc_mapped.tsv",
    summary="out/associations/{target}/nc_summary.tsv"
  log: "out/logs/map_summary_nc_{target}.log"
  shell:
    """
    echo -e "strain\\tunitig\\tcontig\\tstart\\tend\\tstrand\\tupstream\\tgene\\tdownstream" > {output.mapped}
    cat out/associations/{wildcards.target}/mapped/*.txt >> {output.mapped}
    python workflow/scripts/mapped_summary.py {output.mapped} \
           {input.phenotypes} {wildcards.target} {input.filtered} \
           --pangenome {input.pangenome} \
           --length 30 --minimum-hits 9 --maximum-genes 10 \
           --unique \
           > {output.summary}
    """

rule annotate_summary:
  input:
    expand("out/associations/{target}/annotated_summary.tsv",
           target=config["targets"])

rule run_annotate_summary:
  input:
    summary="out/associations/{target}/summary.tsv",
    pangenome=config["pangenome_roary"],
    genes=config["pangenome_genes"]
  output:
    "out/associations/{target}/annotated_summary.tsv"
  params:
    emapper_data=config["emapper"],
    emapper_base="out/associations/{target}/summary",
    sample="out/associations/{target}/sample.faa",
    annotations="out/associations/{target}/summary.emapper.annotations"
  threads: 8 
  conda: "../envs/eggnog-mapper.yaml"
  log: "out/logs/annotate_{target}.log"
  shell:
    """
    python workflow/scripts/sample_pangenome.py {input.pangenome} {input.genes} \
           --groups {input.summary} > {params.sample} && \
    emapper.py -i {params.sample} -o {params.emapper_base} \
               --cpu {threads} --target_orthologs one2one --go_evidence all \
               --tax_scope Bacteria --pfam_realign realign --override \
               --data_dir {params.emapper_data} || touch {output} && \
    python workflow/scripts/enhance_summary.py {input.summary} {params.annotations} \
    > {output}
    """

rule simulations:
  input:
    expand("out/sim/sim_{number}.done",
           number=config["simulations"])

rule generate_simulation:
  input: "config/bacgwasim_{number}.yaml"
  output:
    done="out/sim/sim_{number}.done",
    msa="out/sim/sim_{number}/simulations/genSim/genSim.fasta",
    tree="out/sim/sim_{number}/simulations/genSim/phylogeny.nwk",
    vcf="out/sim/sim_{number}/simulations/genSim/sims.vcf",
    p0="out/sim/sim_{number}/simulations/phenSim/0/phenSim.phen",
    p1="out/sim/sim_{number}/simulations/phenSim/1/phenSim.phen",
    p2="out/sim/sim_{number}/simulations/phenSim/2/phenSim.phen",
    p3="out/sim/sim_{number}/simulations/phenSim/3/phenSim.phen",
    p4="out/sim/sim_{number}/simulations/phenSim/4/phenSim.phen",
    p5="out/sim/sim_{number}/simulations/phenSim/5/phenSim.phen",
    p6="out/sim/sim_{number}/simulations/phenSim/6/phenSim.phen",
    p7="out/sim/sim_{number}/simulations/phenSim/7/phenSim.phen",
    p8="out/sim/sim_{number}/simulations/phenSim/8/phenSim.phen",
    p9="out/sim/sim_{number}/simulations/phenSim/9/phenSim.phen"
  threads: 24
  log: "out/logs/simulation_{number}.log"
  conda: "../envs/bacgwasim.yaml"
  shell:
    """
    cd BacGWASim && python BacGWASim_runner.py --config ../{input} && touch ../{output.done}
    """

rule similarity_simulation:
  input: "out/sim/sim_{number}/simulations/genSim/phylogeny.nwk"
  output: "out/sim/sim_{number}/simulations/genSim/similarity.tsv.gz"
  threads: 1
  log: "out/logs/similarity_simulation_{number}.log"
  conda: "../envs/pyseer.yaml"
  shell:
    """
    python workflow/scripts/phylogeny_distance.py {input} --calc-C | gzip > {output}
    """

rule distance_simulation:
  input: "out/sim/sim_{number}/simulations/genSim/genSim.fasta"
  output: "out/sim/sim_{number}/simulations/genSim/distance.tsv.gz"
  threads: 12
  log: "out/logs/distance_simulation_{number}.log"
  conda: "../envs/mash.yaml"
  shell:
    """
    mash dist -p {threads} -i {input} {input} | square_mash | gzip > {output}
    """

rule phenotypes_simulation:
  input:
    expand("out/sim/sim_{number}/simulations/phenSim/{round}/phenotype.tsv",
           number=config["simulations"],
           round=range(10))

rule run_phenotypes_simulation:
  input: "out/sim/sim_{number}/simulations/phenSim/{round}/phenSim.phen"
  output:
    all="out/sim/sim_{number}/simulations/phenSim/{round}/phenotype.tsv",
    thousand="out/sim/sim_{number}/simulations/phenSim/{round}/phenotype.1000.tsv",
    fivethousand="out/sim/sim_{number}/simulations/phenSim/{round}/phenotype.5000.tsv"
  threads: 1
  log: "out/logs/phenotypes_simulation_{number}_{round}.log"
  shell:
    """
    python workflow/scripts/sim2phen.py {input} > {output.all}
    python workflow/scripts/sim2phen.py {input} --sample 0.1 --seed {wildcards.round} > {output.thousand}
    python workflow/scripts/sim2phen.py {input} --sample 0.5 --seed {wildcards.round} > {output.fivethousand}
    """

rule pyseer_simulation:
  input:
    expand("out/sim/sim_{number}/associations/{round}/vcf.tsv",
           number=config["simulations"],
           round=range(10))

rule run_pyseer_simulation:
  input:
    all="out/sim/sim_{number}/simulations/phenSim/{round}/phenotype.tsv",
    thousand="out/sim/sim_{number}/simulations/phenSim/{round}/phenotype.1000.tsv",
    fivethousand="out/sim/sim_{number}/simulations/phenSim/{round}/phenotype.5000.tsv",
    distance="out/sim/sim_{number}/simulations/genSim/distance.tsv.gz",
    vcf="out/sim/sim_{number}/simulations/genSim/sims.vcf"
  output:
    vcf="out/sim/sim_{number}/associations/{round}/vcf.tsv",
    vcf1000="out/sim/sim_{number}/associations/{round}/vcf.1000.tsv",
    vcf5000="out/sim/sim_{number}/associations/{round}/vcf.5000.tsv",
    patterns="out/sim/sim_{number}/associations/{round}/patterns.txt", 
    patterns1000="out/sim/sim_{number}/associations/{round}/patterns.1000.txt",
    patterns5000="out/sim/sim_{number}/associations/{round}/patterns.5000.txt" 
  threads: 4
  log: "out/logs/pyseer_simulation_{number}_{round}.log"
  conda: "../envs/pyseer.yaml"
  shell:
    """
    pyseer --vcf {input.vcf} \
           --phenotypes {input.all} \
           --distance {input.distance} \
           --output-patterns {output.patterns} \
           --cpu {threads} \
           --block_size 1000 \
           --max-dimensions 4 \
           --min-af 0.001 \
           --max-af 0.999 \
           > {output.vcf} &&
    pyseer --vcf {input.vcf} \
           --phenotypes {input.thousand} \
           --distance {input.distance} \
           --output-patterns {output.patterns1000} \
           --cpu {threads} \
           --block_size 1000 \
           --max-dimensions 4 \
           --min-af 0.01 \
           --max-af 0.99 \
           > {output.vcf1000} &&
    pyseer --vcf {input.vcf} \
           --phenotypes {input.fivethousand} \
           --distance {input.distance} \
           --output-patterns {output.patterns5000} \
           --cpu {threads} \
           --block_size 1000 \
           --max-dimensions 4 \
           --min-af 0.002 \
           --max-af 0.998 \
           > {output.vcf5000}
    """
