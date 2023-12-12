configfile: "config.yaml"
include: "utils.py"

prefix = config["prefix"]
bedfile = config["bedfile"]
chr_value = config["chr"]
is_ARN = config["isARN"]
max_mem = config["max_mem"]
patient = config["patient"]
igh_data = config["igh_data"]   

#Pour lancer le fichier en ligne de commande :
#snakemake -s snakefile.smk --cores 8 --config prefix=4cd74e8e5c0f2ec4fd155ccad999f585 isARN=yes bedfile=chemin
#snakemake -s snakefile.smk --cores 8 --config prefix=4cd74e8e5c0f2ec4fd155ccad999f585 chr=chr1
#snakemake -s snakefile.smk --cores 8 --config prefix=4cd74e8e5c0f2ec4fd155ccad999f585 chr=chr1 isARN=yes
#snakemake -s Snakefile.smk --cores 8 --config prefix=PCAWG.b22c5a18-f6e0-11e4-8ef4-b5a7fb958afb.STAR.v1 isARN=yes bedfile=bed/hg38ig.bed
#snakemake -s Snakefile.smk --cores 8 --config prefix=e7cf3d7837e1be92a4be335629152069 chr=chr1 patient= DO52700 
#prefix : (nom du fichier sans extension .bam)
#si aucun fichier BED ne rient ecrire
class MonException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


rule all:
    input:
        expand("bam/{prefix}.bam",prefix=config["prefix"]),
        expand("bam/{prefix}.txt",prefix=config["prefix"]),
        expand("fastq/{prefix}.txt", prefix=config["prefix"]),
        expand("fastq/{prefix}_R1.fastq",prefix=config["prefix"]),
        expand("fastq/{prefix}_R2.fastq",prefix=config["prefix"]),
        expand("fastq/{prefix}_single.fastq",prefix=config["prefix"]), 
        expand("fastq/{prefix}.tx2t",prefix=config["prefix"]),
        expand("reportfalco/{prefix}",prefix=config["prefix"]),
        expand("reportfalco/{prefix}_R1",prefix=config["prefix"]),
        expand("reportfalco/{prefix}_R2",prefix=config["prefix"]),
        expand("reportfalco/{prefix}_R1/summary.txt",prefix=config["prefix"]),
        expand("trimfastq/{prefix}_R1_trim.fastq",prefix=config["prefix"]),
	    expand("trimfastq/{prefix}_R2_trim.fastq",prefix=config["prefix"]),
        expand("trimfastq/{prefix}.txt",prefix=config["prefix"]),
        expand("sam/{prefix}_PEM.sam",prefix=config["prefix"]),
        expand("sambam/{prefix}_PEM_sam.bam",prefix=config["prefix"]),
        expand("sambam/{prefix}.txt",prefix=config["prefix"]),
        expand("bamfiltered/{prefix}_PEM_sam.bam",prefix=config["prefix"]),
        expand("bamfiltered/{prefix}.txt",prefix=config["prefix"]),
        expand("sambamsort/{prefix}_PEM_sam_sort.bam",prefix=config["prefix"]),
        expand("sambamsort/{prefix}.txt",prefix=config["prefix"]),
        expand("sambamsort/{prefix}.txtp",prefix=config["prefix"]),
        expand("sambamsortNoDup/{patient}/{prefix}_PEM_sam_sort.bam", prefix=prefix,patient=patient),
        expand("sambamsortNoDup/{patient}/{prefix}.txt", prefix=prefix,patient=patient),
        expand("sambamsortNoDup/{patient}/{prefix}.txtp", prefix=prefix,patient=patient),
        expand("sambamchorbed/{patient}/{prefix}_PEM_sam_sort_extracted.bam",prefix=prefix,patient=patient),
        expand("sambamchorbed/{patient}/{prefix}.txt",prefix=prefix,patient=patient),
        expand("fastq/{patient}/{prefix}_PEM_sam_sort_extracted_single.fastq", prefix=prefix,patient=patient),
        expand("fastq/{patient}/{prefix}_PEM_sam_sort_extracted.txt", prefix=prefix,patient=patient),
        expand("trinity/{patient}/trinity_{prefix}/{prefix}_PEM_sam_sort_extracted_fasqR.txt", prefix=prefix,patient=patient),
        expand("blastn/{patient}/{prefix}_PEM_sam_sort_extracted_filtered.fasta", prefix=prefix,patient=patient),
        expand("blastn/{patient}/{prefix}_PEM_sam_sort_extracted_filtered.txt", prefix=prefix,patient=patient),
        expand("igblast_report/{patient}/{prefix}.IgBLAST_out.txt", prefix=prefix,patient=patient),
        expand("igblast_report/{patient}/{prefix}.txt", prefix=prefix,patient=patient),



rule sorted:
    input:
        bamfile = "/mnt/project2/tmp70T/gth/ICGC/coverage_10X/{prefix}.bam"
    output:
        single = "bam/{prefix}.bam",
        succes_file = "bam/{prefix}.txt"
    threads:
        8
    shell:
        """
        samtools sort -n {input.bamfile} -o {output.single}
        echo {wildcards.prefix} ok first sorted > {output.succes_file} 
        """

#conversion bam en fastq
#il est important de preciser que les paires d'ADN doivent etre separe
#loption -@ permet de preciser le nombre de cpu a utiliser
rule decompression:
    input:
        bamfile = "bam/{prefix}.bam",
        succes_file = "bam/{prefix}.txt"
    output:
        single = "fastq/{prefix}_single.fastq",
        succes_file = "fastq/{prefix}.txt"
    threads:
        8
    shell:
        """
        samtools bam2fq -@ {threads} {input.bamfile} > {output.single} 
        echo {wildcards.prefix} ok > {output.succes_file} 
        """

rule decuplate:
    input:
        fqsingle = "fastq/{prefix}_single.fastq",
        succes_file = "fastq/{prefix}.txt"
    output:
        r1 = "fastq/{prefix}_R1.fastq",
        r2 = "fastq/{prefix}_R2.fastq",
        succes_file = "fastq/{prefix}.tx2t"
    run:
        shell("cat {input.fqsingle} | grep '^@.*/1$' -A 3 --no-group-separator > {output.r1}")
        shell("cat {input.fqsingle} | grep '^@.*/2$' -A 3 --no-group-separator > {output.r2}")
        #shell("cat {input.fqsingle} | tee >(cut -f 1-4 | tr '\\t' '\\n' > {output.r1}) | cut -f 5-8 | tr '\\t' '\\n' > {output.r2}")
        shell("echo {wildcards.prefix} decuplate > {output.succes_file}")

#on verifie la qualite des fichiers generes avec falco
rule testquatility:
    input:
        fqr1 = "fastq/{prefix}_R1.fastq",
        fqr2 = "fastq/{prefix}_R2.fastq",
        fqsingle = "fastq/{prefix}_single.fastq",
        succes_file = "fastq/{prefix}.txt"
    output:
        falrep1 = directory("reportfalco/{prefix}_R1"),
        falrep2 = directory("reportfalco/{prefix}_R2"),
        falrep = directory("reportfalco/{prefix}"),
        falreport1 =  "reportfalco/{prefix}_R1/summary.txt"
    shell:
        """
        falco {input.fqr1} -o {output.falrep1}
        falco {input.fqr2} -o {output.falrep2}
        falco {input.fqsingle} -o {output.falrep}
        """


#si les fichiers ne sont pas tous de bonnes qualites
#il faut faire le trimming (avec fastp)
#il faut faire attention aux ARN, si la qualité est bonne, pas de trimming
rule trimming:
    input:
        fqr1_in = "fastq/{prefix}_R1.fastq",
        fqr2_in = "fastq/{prefix}_R2.fastq",
        falrep1 = "reportfalco/{prefix}_R1",
        falreport1 =  "reportfalco/{prefix}_R1/summary.txt"
    output:
        fqr1_trim = "trimfastq/{prefix}_R1_trim.fastq",
        fqr2_trim = "trimfastq/{prefix}_R2_trim.fastq",
        trimdir = directory("trimfastq/{prefix}"),
        succes_file = "trimfastq/{prefix}.txt"
    run:
        quality = "Per base sequence quality"
        adaptor = "Adapter Content"
        #print({input.falrep})
        try:
            shell("mkdir {output.trimdir}")
            if config['isARN']=='no':
                shell("fastp -i {input.fqr1_in} -I {input.fqr2_in} -o {output.fqr1_trim} -O {output.fqr2_trim}") 
                shell("echo {wildcards.prefix} ok > {output.succes_file}")
            elif config['isARN']=='yes' and report_falco_val(quality,input.falreport1)== 'FAIL' or report_falco_val(adaptor,input.falreport1)=='FAIL':
                shell("trim_galore --paired {input.fqr1_in} {input.fqr2_in} --output_dir {output.trimdir}")
                shell("echo {wildcards.prefix} ok > {output.succes_file}")
            elif config['isARN']=='yes':
                shell("cp {input.fqr1_in} {output.fqr1_trim}")
                shell("cp {input.fqr2_in} {output.fqr2_trim}")
                shell("echo {wildcards.prefix} no trimming > {output.succes_file}")
            elif config['isARN']!='yes' and config['isARN']!='no':
                raise MonException("isARN expected value (yes or no) but have " + str(config['isARN']))
        except MonException as e:
            print("\033[91mError when we try to execute trimmimg :", e.message, "\033[0m")

#Lorsque tout est correct
#convertir les fastq en sam (alignement)
# il faut a ce niveau toujours preciser le genome de reference (pour notre cas hg38)
# loption -t permet de preciser le nombre de cpu a utiliser
# il est important de noter que lagniment se fait en paired end mode
# pour generer lindex STAR :
# STAR --runMode genomeGenerate --genomeDir /mnt/project3/hgayap/pipeline/refgenomestar --genomeFastaFiles /mnt/project3/hgayap/pipeline/refgenome/hg38.fa --sjdbGTFfile /mnt/project3/hgayap/pipeline/refgenome/hg38.refGene.gtf
rule align:
    input:
        ref = "/mnt/project3/hgayap/pipeline/refgenome/hg38.fa",
        index = "/mnt/project3/hgayap/pipeline/refgenomestar/",
        fqr1_trim_in = "trimfastq/{prefix}_R1_trim.fastq",
        fqr2_trim_in = "trimfastq/{prefix}_R2_trim.fastq",
        succes_file_in = "trimfastq/{prefix}.txt"
    output:
        paired_al ="sam/{prefix}_PEM.sam"
    threads: 8
    run:
        if config['isARN']=='yes':
            shell("STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.fqr1_trim_in} {input.fqr2_trim_in} --outFileNamePrefix aligned_reads_")
            shell("mv aligned_reads_Aligned.out.sam {output.paired_al}")
        else:
            shell("bwa mem -t {threads} {input.ref} {input.fqr1_trim_in} {input.fqr2_trim_in} > {output.paired_al}")


# on doit ensuite faire la conversion des .sam en .bam 
# lobjectif est de pouvoir indexer les fichiers aligner pour les visualer correctement
rule sam_to_bam:
    input:
        paired_al_in="sam/{prefix}_PEM.sam"
    output:
        sambam="sambam/{prefix}_PEM_sam.bam",
        succes_file="sambam/{prefix}.txt"
    threads: 8
    shell:
        """
        samtools view -@ {threads} -bS -o {output.sambam} {input.paired_al_in} 
        echo {wildcards.prefix} ok > {output.succes_file} 
        """

rule filter:
    input:
        sambam = "sambam/{prefix}_PEM_sam.bam",
        succes_file="sambam/{prefix}.txt"
    
    output:
        bamfilt = "bamfiltered/{prefix}_PEM_sam.bam",
        succes_file="bamfiltered/{prefix}.txt"
    shell:
        """
        samtools view -b -F 0x4 -f 0x2 -q 20 {input.sambam} > {output.bamfilt}
        echo {wildcards.prefix} filtered > {output.succes_file} 
        """


# Il faut trier car lordre dans lequel les lectures sont alignees peuvent etre desordonnee
rule sort:
    input:
        bamfilt = "bamfiltered/{prefix}_PEM_sam.bam",
        succes_file="bamfiltered/{prefix}.txt"
    output:
        bamsort="sambamsort/{prefix}_PEM_sam_sort.bam",
        succes_file="sambamsort/{prefix}.txt"
    threads: 8
    shell:
        """
        samtools sort -@ {threads} {input.bamfilt} -o {output.bamsort}
        echo {wildcards.prefix} ok > {output.succes_file} 
        """


# On peut desormais l'indexer 
rule indexation:
    input:
        bamsort_in = "sambamsort/{prefix}_PEM_sam_sort.bam",
        succes_file = "sambamsort/{prefix}.txt"
    output:
        succes_file = "sambamsort/{prefix}.txtp"
    shell:
        """
        samtools index {input.bamsort_in}
        echo {wildcards.prefix} Indexed success > {output.succes_file}
        """

rule dupMark:
    input:
        bamsort_in = "sambamsort/{prefix}_PEM_sam_sort.bam",
        succes_file = "sambamsort/{prefix}.txtp"
    
    output:
        bamsort_out = "sambamsortNoDup/{patient}/{prefix}_PEM_sam_sort.bam",
        succes_file = "sambamsortNoDup/{patient}/{prefix}.txt"

    shell:
        """
        samtools rmdup {input.bamsort_in} {output.bamsort_out}
        echo {wildcards.prefix} duplicate removed success > {output.succes_file}
        """
# On peut desormais l'indexer 
rule indexation2:
    input:
        bamsort_in = "sambamsortNoDup/{patient}/{prefix}_PEM_sam_sort.bam",
        succes_file = "sambamsortNoDup/{patient}/{prefix}.txt"
    output:
        succes_file = "sambamsortNoDup/{patient}/{prefix}.txtp"
    shell:
        """
        samtools index {input.bamsort_in}
        echo {wildcards.prefix} Indexed success > {output.succes_file}
        """


#la visualisation avec loutils igv peut etre possible desormais

################################## for BED FILE or CHRx ###################################################
#pour extraire un chromosome ou une sequence particuliere BED (chr1 par exemple)
rule ch_extraction:
    input:
        bamsort_in = "sambamsortNoDup/{patient}/{prefix}_PEM_sam_sort.bam",
        succes_file = "sambamsortNoDup/{patient}/{prefix}.txtp"
    output:
        sambamch = "sambamchorbed/{patient}/{prefix}_PEM_sam_sort_extracted.bam",
        succes_file = "sambamchorbed/{patient}/{prefix}.txt"
    threads: 8
    run:
        if config["bedfile"] != 'EMPTY':
            shell("samtools view -@ {threads} -b -L {bedfile} {input.bamsort_in} > {output.sambamch}")
            shell("echo {wildcards.prefix} Extract successfully > {output.succes_file}")
        elif is_chrom(chr_value):
            shell("samtools view -@ {threads} -b -h {input.bamsort_in} {chr_value} > {output.sambamch}")
            shell("echo {wildcards.prefix} Extract successfully > {output.succes_file}")


rule decompression2:
    input:
        bamfile = "sambamchorbed/{patient}/{prefix}_PEM_sam_sort_extracted.bam",
        succes_file = "sambamchorbed/{patient}/{prefix}.txt"
    output:
        single = "fastq/{patient}/{prefix}_PEM_sam_sort_extracted_single.fastq",
        succes_file = "fastq/{patient}/{prefix}_PEM_sam_sort_extracted.txt"
    threads:
        8
    shell:
        """
        samtools fastq -@ {threads} {input.bamfile} > {output.single}
        echo {wildcards.prefix} ok > {output.succes_file} 
        """

#conversion bam en fastq
rule assembling:
    input:
        single = "fastq/{patient}/{prefix}_PEM_sam_sort_extracted_single.fastq",
        succes_file = "fastq/{patient}/{prefix}_PEM_sam_sort_extracted.txt"
    output:
        trin_dir = directory("trinity/{patient}/trinity_{prefix}"),
        trin_file = "trinity/{patient}/trinity_{prefix}.Trinity.fasta",
        succes_file = "trinity/{patient}/trinity_{prefix}/{prefix}_PEM_sam_sort_extracted_fasqR.txt"
    threads:
        8
    run:
        if is_ARN == 'yes':
            shell("Trinity --seqType fq --max_memory {max_mem} --single {input.single} --CPU {threads} --output {output.trin_dir} --trimmomatic --full_cleanup --no_bowtie")
            shell("echo {wildcards.prefix} ok > {output.succes_file}")
        else:
            shell("echo {wildcards.prefix} ok no ARN sorry > {output.trin_file}") #To change when ADN 
            shell("echo {wildcards.prefix} ok > {output.succes_file}")

#Filtrer lARN
rule filtering:
    input:
        trin_dir = "trinity/{patient}/trinity_{prefix}",
        trin_file = "trinity/{patient}/trinity_{prefix}.Trinity.fasta",
        succes_file = "trinity/{patient}/trinity_{prefix}/{prefix}_PEM_sam_sort_extracted_fasqR.txt"
    output:
        blastfile = "blastn/{patient}/{prefix}_PEM_sam_sort_extracted_filtered.fasta",
        succes_file = "blastn/{patient}/{prefix}_PEM_sam_sort_extracted_filtered.txt"
    threads: 8
    run:
        if is_ARN == 'yes':
            shell("blastn -db {igh_data}/IGHV -query {input.trin_file} -out {output.blastfile} -num_threads {threads} -outfmt 6")
            shell("echo {wildcards.prefix} ok > {output.succes_file}")
        else:
            shell("echo {wildcards.prefix} ok no ARN Sorry > {output.blastfile}") #To change
            shell("echo {wildcards.prefix} ok > {output.succes_file}")


#cette règle prend un fichier fasta en entrée, exécute l'outil "igblastn" 
#pour analyser la séquence et génère un rapport IgBLAST
rule seqAnalys:
    input:
        trin_dir = "trinity/{patient}/trinity_{prefix}",
        trin_file = "trinity/{patient}/trinity_{prefix}.Trinity.fasta",
        succes_file = "trinity/{patient}/trinity_{prefix}/{prefix}_PEM_sam_sort_extracted_fasqR.txt"
    output:
        igblast_file = "igblast_report/{patient}/{prefix}.IgBLAST_out.txt",
        succes_file = "igblast_report/{patient}/{prefix}.txt"
    threads:
        8
    run:
        ARN_is = config["isARN"]
        if  ARN_is == 'yes':
            shell("igblastn -germline_db_V {igh_data}/IGHV -num_alignments_V 3 -germline_db_J {igh_data}/IGHJ -num_alignments_J 3 -germline_db_D {igh_data}/IGHD -num_alignments_D 3 -organism human -query {input.trin_file} -show_translation -auxiliary_data {igh_data}/human_gl.aux > {output.igblast_file}")
            shell("echo {wildcards.prefix} ok > {output.succes_file}")
        else:
            shell("echo {wildcards.prefix} ok no ARN Sorry > {output.igblast_file}") #To change when ADN side is ok
            shell("echo {wildcards.prefix} ok > {output.succes_file}")


################################## END for BED FILE ###################################################




