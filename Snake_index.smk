configfile: "config.yaml"
include: "utils.py"

prefix = config["prefix"]

rule all:
    input:
        expand("/mnt/project2/tmp70T/gth/ICGC/coverage_10X/{prefix}.bam.bai", prefix=config["prefix"]),

################################## for BED FILE or CHRx ###################################################
#pour extraire un chromosome ou une sequence particuliere BED (chr1 par exemple)
# On peut desormais l'indexer 
rule indexation:
    input:
        bam = "/mnt/project2/tmp70T/gth/ICGC/coverage_10X/{prefix}.bam"
    output:
        bam_index = "/mnt/project2/tmp70T/gth/ICGC/coverage_10X/{prefix}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """