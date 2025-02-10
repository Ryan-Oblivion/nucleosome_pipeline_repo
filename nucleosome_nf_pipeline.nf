

//modules

/*
FASTQC='fastqc/0.11.7'
FASTP='fastp/0.22.0'
BWA='bwa/0.7.17'
SAMTOOLS='samtools/1.16'
BEDTOOLS='bedtools/2.30.0'
DEEPTOOLS='deeptools/3.5.1'
MACS2='macs2/2.1.1'
R_LAN='r/4.3.2'
*/

// making a parameter for the user to input a reference genome

//params.mouse_ref_mm10='/gpfs/home/rj931/projects/chip_seq_practice/nextflow_workstation/mouse_genome_grcm38/mm10.fa.gz'


// process for downloading the mouse reference genome and gtf


// enable dsl 2
nextflow.enable.dsl=2
/*
process mouse_ref_genome_gtf {
    cache true
    publishDir './ref_genome_dir', mode:'copy', pattern:'*.fa.gz'
    publishDir './gtf_genome_dir', mode:'copy', pattern:'Mus_musculus*.gtf.gz'

    input:

    output:
    path '*.fa.gz', emit: ref_genome
    path 'Mus_musculus*.gtf.gz', emit: gtf_file

    """
    #!/bin/env bash

    # This is the mouse reference genome
    curl -O https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz

    # this is GRCm38.p6 version of the mouse genome
    #curl -o Mus_musculus.GRCm38.p6.fasta https://www.ebi.ac.uk/ena/browser/api/embl/GCA_000001635.8?download=true&gzip=true && gzip -f Mus_musculus.GRCm38.p6.fasta

    # this is the mouse gtf
    curl -O https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz


    """

}


process fastp_se { 
    cache true
    //executor 'slurm'
    memory '30 GB'
    
    publishDir './fastp_filter_output', mode:'copy', pattern:'*'

    input:
    tuple val(control), path(control_reads)
    tuple val(treatment), path(treatment_reads)
    val(filenames_control)
    val(filenames_treatment)
    
    output:

    tuple path("${filenames_control}_filt.fastq"), path("${filenames_treatment}_filt.fastq"), emit: filtered_fastq_c_t

    // output the .html to a channel

    path "*.html", emit: fastp_html
    
    script:

    """
    #!/bin/env bash

    module load $FASTP
    module load $FASTQC


    # this is single end use
    fastp \
    -i ${control_reads} \
    -o $filenames_control'_filt.fastq' \
    -h $filenames_control'.html' \
    --qualified_quality_phred=15 \
    --dedup \
    --dup_calc_accuracy=4 

    fastp \
    -i ${treatment_reads} \
    -o $filenames_treatment'_filt.fastq' \
    -h $filenames_treatment'.html' \
    --qualified_quality_phred=15 \
    --dedup \
    --dup_calc_accuracy=4 

    """
}






// process for making the ref index

process bwa_index {
    cache true
    //memory '100 GB'

    publishDir './ref_indices', mode:'copy', pattern:'*'

    input:

    path ref

    

    output:
    path '*', emit: ref_index_files

    """
    #!/bin/env bash

    module load $BWA

    bwa index \
    -p $ref \
    -a bwtsw \
    $ref
    """
}


//process for alignment of filtered reads to reference genome.
// alignment of the filtered reads to the reference genome

process bwa_alignment {

    cache true

    publishDir './Chip_sam_storage', mode:'copy', pattern:'*.sam'

    input:

    tuple path(control_reads), path(treatment_reads)
    path ref
    path index_files
    path gtf
    tuple val(control_name), val(treatment_name)

    output:

    tuple path("${sam_control_file}"), path("${sam_treatment_file}"), emit: sam_files_tuple_ch


    script:

    sai_control_output = "${control_name}.sai"
    sam_control_file = "${control_name}.sam"

    sai_treatment_output = "${treatment_name}.sai"
    sam_treatment_file = "${treatment_name}.sam"

    """
    #!/bin/env bash

    module load $BWA

    # aligning the control files

    bwa aln \
    "${ref}" \
    "${control_reads}" \
    > "${sai_control_output}"


    bwa samse \
    "${ref}" \
    "${sai_control_output}" \
    "${control_reads}" \
    > "${sam_control_file}"


    # now aligning the treatment files

    bwa aln \
    "${ref}" \
    "${treatment_reads}" \
    > "${sai_treatment_output}"


    bwa samse \
    "${ref}" \
    "${sai_treatment_output}" \
    "${treatment_reads}" \
    > "${sam_treatment_file}"

    """
}


// I should use samtools to order the coordinates and then output that as a bam file

process samtools_convert {

    publishDir './Sorted_bam_files', mode: 'copy', pattern: '*.bam'
    publishDir './Sorted_bam_files', mode: 'copy', pattern: '*.bai'
    publishDir './Sorted_bai_files', mode: 'copy', pattern: '*.bai'
    publishDir './Bigwig_files', mode: 'copy', pattern: '*.bigwig'


    input:

    tuple path(sam_control), path(sam_treatment) 
    tuple val(control_basename), val(treatment_basename)
    
    output:

    tuple path("${out_control_file}"), path("${out_treatment_file}"), emit: sorted_bam_tuple_ch
    path '*.bai', emit: sorted_bai_ch
    tuple path("${out_control_bigwig}"), path("${out_treatment_bigwig}"), emit: bigwig_c_t_tuple
    

    script:

    out_control_file = "${control_basename}_sorted.bam"
    out_treatment_file = "${treatment_basename}_sorted.bam"

    out_control_bigwig = "${control_basename}_normalized.bigwig"
    out_treatment_bigwig = "${treatment_basename}_normalized.bigwig"

    """
    #!/bin/env bash

    module unload $MACS2
    module load $SAMTOOLS
    module load $DEEPTOOLS

    samtools sort \
    -O bam \
    -o "${out_control_file}" \
    "${sam_control}" 

    samtools sort \
    -O bam \
    -o "${out_treatment_file}" \
    "${sam_treatment}" 

    # I need to index the bam files for use in bamCoverage

    samtools index \
    "${out_control_file}"

    samtools index \
    "${out_treatment_file}" 


    # If I need fast random access, I can use samtools index and use the coordinate sorted bam file as input
    # the output will be a .bai file.

    # maybe make a normalized bigwig file using rpm

    # using deeptools I can convert bed to bigwig and normalize by rpm or other methods

    # making the control bigwig file first
    # using the table that shows the effective genome size for a mouse genome grcm39
    bamCoverage \
    -b "${out_control_file}" \
    -o "${out_control_bigwig}" \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2654621783

    # now doing the same for treatment bam to bigwig 
    bamCoverage \
    -b "${out_treatment_file}" \
    -o "${out_treatment_bigwig}" \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2654621783

    


    ######### NOT USING THIS SECTION #############################
    #scaling_factor_control=\$(samtools idxstats "{out_control_file}" | '{sum+=\$3}END{print sum/1000000}')
    
    #bedtools genomecov \
    #-ibam "{out_control_file}" \
    #-scale scaling_factor_control \
    #-bg \
    #> "{out_control_bedgraph}"
    

    # now converting the bedgraph to bigwig
    
    
    #scaling_factor_treatment=\$(samtools idxstats "{out_treatment_file}" | '{sum+=\$3}END{print sum/1000000}')


    #bedtools genomecov \
    #-ibam "{out_treatment_file}" \
    #-scale scaling_factor_treatment \
    #-bg \
    #> "{out_treatment_bedgraph}"


    """
}


// now to use macs2 to call the peaks

process peak_calling_mouse {

    publishDir './macs_outputs', mode: 'copy', pattern:'*'

    input:
    // i am doing it this way because i didnt make separate processes for each treatment and control groups
    tuple path(control_bam), path(treatment_bam)
    tuple val(control_name), val(treatment_name)

    // inputting the bigwig files to make the peak count table with deeptools
    //tuple path(control_bigwig), path(treatment_bigwig)
    //tuple val(control_bw_basename), val(treatment_bw_basename)
    
    output:

    path "*", emit: peak_files

    script:
    // macs2 does not allow for out file name
    experiment_name = "${control_name}_vs_${treatment_name}"

    //peak_count_table_out = "${control_bw_basename}_vs_${treatment_bw_basename}"

    """
    #!/bin/env bash
    # module unload \$DEEPTOOLS
    module unload $MACS2
    module unload python
    module load $MACS2
    

     # the -g flag is showing that the genome size is comming from the mouse
    # the -t flag takes the treatment group of experiments
    # the -c flag takes the control group of experiments
    # the -f flag specifies that I am using bam files if BAM is set
    # -B flag is whether or not to save extended fragment pileup
    # --trackline tells macs to include trackline with gedgraph files
    # --SPMR saves signal per million reads for fragment pileup profiles
    # -q flag is the minimum FDR cutoff for peak detection
    # --cutoff-analysis will give a summary of number of total peaks that can be called by different thresholds. so we can make a decision
    # -n flag is the name of the experiment that will be used as prefix for output file names
    macs2 callpeak \
    -t "${treatment_bam}" \
    -c "${control_bam}" \
    -f BAM \
    -g 'mm' \
    -B \
    --trackline \
    --SPMR \
    -q 0.05 \
    --cutoff-analysis \
    --outdir '.' \
    -n "${experiment_name}"


    ########### DOING THIS IN THE MAKE PEAK COUNT TABLE PROCESS ################
    # I need to unload macs2 and then load deeptools for it to work
    #module unload \$MACS2
    #module unload \$DEEPTOOLS
    #module load \$DEEPTOOLS

    # now I want to create a peak count table using deeptools command multiBigwigSummary BED-file

    #multiBigwigSummary BED-file \
    #-b "\${control_bigwig}" "\${treatment_bigwig}" \
    #-o "\${peak_count_table_out}" \
    #--BED "*.bed"
    
    """

}
*/

process peak_calling_jan {

    // many ways to use conda in the processes of nextflow. I will use the path to an existing conda envrionment
    // but it may be better to use the path to a conda yml config file later on for portability
    conda '/lustre/fs4/home/rjohnson/miniconda3/envs/macs2_rj'
    // this is the conda yml env file below
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/macs2_rj_env.yml'

    publishDir './macs_outputs', mode: 'copy', pattern:'*'

    input:
    

    // these are the bam files that I put in this process in the workflow section
    path(bam_sub100_files_ch) 
    val(bam_sub100_names_ch) //, path(bam_20u), path(bam_100u), path(bam_500u)

    // these below are the names i gave them to make it easier to work with the files
    // with input channels you can call them different names here but they retain the order so adding '_name' just makes it easier for me to know i am using the name and not the file
    
    //tuple val(bam_4u_name) //, val(bam_20u_name), val(bam_100u_name), val(bam_500u_name)
    
    // inputting the bigwig files to make the peak count table with deeptools
    //tuple path(control_bigwig), path(treatment_bigwig)
    //tuple val(control_bw_basename), val(treatment_bw_basename)
    
    output:

    path "*", emit: peak_files
    path "*summits.bed", emit: summits_bed_file // this is used to make the peak count table

    script:
    // macs2 does not allow for out file name
    //experiment_name = "${control_name}_vs_${treatment_name}"

    //peak_count_table_out = "${control_bw_basename}_vs_${treatment_bw_basename}"

    experiment_name = "${bam_sub100_names_ch}"

    """
    #!/bin/env bash
    # module unload \$DEEPTOOLS
    #module unload \$MACS2
    #module unload python
    #module load \$MACS2
    

     # the -g flag is showing that the genome size is comming from the mm mouse or hs from human
     # the --nolambda flag so that MACS2 does not attempt to automatically build the local λ background from your treatment alone.
    # the -t flag takes the treatment group of experiments
    # the -c flag takes the control group of experiments
    # the -f flag specifies that I am using bam files if BAM is set
    # -B flag is whether or not to save extended fragment pileup
    # --trackline tells macs to include trackline with gedgraph files
    # --SPMR saves signal per million reads for fragment pileup profiles
    # -q flag is the minimum FDR cutoff for peak detection
    # --cutoff-analysis will give a summary of number of total peaks that can be called by different thresholds. so we can make a decision
    # -n flag is the name of the experiment that will be used as prefix for output file names
    macs2 callpeak \
    -t "${bam_sub100_files_ch}" \
    -f BAM \
    -g 'hs' \
    --nolambda \
    -B \
    --trackline \
    --SPMR \
    -q 0.05 \
    --outdir '.' \
    -n "${experiment_name}"


    ########### DOING THIS IN THE MAKE PEAK COUNT TABLE PROCESS ################
    # I need to unload macs2 and then load deeptools for it to work
    #module unload \$MACS2
    #module unload \$DEEPTOOLS
    #module load \$DEEPTOOLS

    # now I want to create a peak count table using deeptools command multiBigwigSummary BED-file

    #multiBigwigSummary BED-file \
    #-b "\${control_bigwig}" "\${treatment_bigwig}" \
    #-o "\${peak_count_table_out}" \
    #--BED "*.bed"
    
    """

}






process sem_process {

    publishDir './CAD/GM12878/CHOMP/SEM_outputs', mode: 'copy', pattern: '*'
    conda '/ru-auth/local/home/risc_soft2/miniconda3/envs/SEM_env'

    // using the bamtools or samtools conda environment to look at the bam file
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools_rj_env.yml'
    cache true


    input:
    tuple path(genome_dict), path(genome_fai), path(genome_fa), path(genome_xml) 
    //path(genome_files)
    path(bam_files)
    val(bam_names)

    output:

    path "${out_file_name}*", emit: sem_out_files_ch

    script:

    out_file_name = "${bam_names}"

    """
    #!/bin/env bash

    echo 'this is the dict: ${genome_dict}; this is the fai: ${genome_fai}'
    # this shows that the correct files are in the correct variable name. just remove the comment hashtag to check your self. in the .command.log file

    # i am using the tool sem from yenlab https://github.com/YenLab/SEM?tab=readme-ov-file#readme
    
    #### parameters used ####
    # --Xmx50G;  this is requesting memory and without this sem didn't run. so this helped even though i requested 50gb in the nextflow config file. why?
    # --out; takes the output file prefix name
    # --geninfo; takes the genome info file (fa or fai not sure yet) the example has the fai file
    # --expt; takes the input file name (the bam file)
    # --format; you need to specify the type of file; its either you use this or the expt
    #########################

    # checking the format of the input files
    #head -n 5 '${genome_fai}'
    #tail -n 5 '${genome_fai}'

    # this shows that the bam files and the genome file have the same chr naming convention
    #samtools view '${bam_files}' | head -n 10
    

    sem \
    -Xmx50G \
    --out '${out_file_name}' \
	--geninfo '${genome_fai}' \
	--expt '${bam_files}' \
    --numClusters 4 \
    --threads 8
    

    """
}





process make_peak_count_table {

    publishDir './macs_outputs', mode: 'copy', pattern:'*'

    input:

    // inputting the bigwig files to make the peak count table with deeptools
    tuple path(control_bigwig), path(treatment_bigwig)
    tuple val(control_bw_basename), val(treatment_bw_basename)

    tuple file(lambda_bdg), file(analysis_txt), file(model_r), file(peak_narrowpeak), file(peaks_xls), file(summits_bed), file(pileup_bdg)

    output:

    path "${peak_count_table_out}*", emit: peak_count_table_ch

    script:

    peak_count_table_out = "${control_bw_basename}_vs_${treatment_bw_basename}"
    bed_file = "${summits_bed}"

    """
    # I need to unload macs2 and then load deeptools for it to work
    module unload $MACS2
    module load $DEEPTOOLS

    # now I want to create a peak count table using deeptools command multiBigwigSummary BED-file

    multiBigwigSummary BED-file \
    -b "${control_bigwig}" "${treatment_bigwig}" \
    -o "${peak_count_table_out}_peak_counts.npz" \
    --outRawCounts "${peak_count_table_out}_peak_counts.tsv" \
    --BED "${bed_file}"

    """

}



// I need to create the plots with ngsplot
/*
process ngs_plot {

    publishDir './Visualizations', mode: 'copy', pattern: '*.pdf'

    input:

    path bam
    path bai

    output:
    
    path "*.pdf"

    script:


    """
    #!/bin/env bash

    module load $R_LAN
    # setting up the environment and tool to use ngsplot with the mouse genome

    # first downloading the mouse genome from their google drive
    # i got the file_id and put it equal to the link for google drive provided elsewhere

    #curl -o ngsplotdb_mm10_75_3.00.tar.gz \
    # "https://drive.google.com/uc?export=download&id=0B5hDZ2BucCI6NXNzNjZveXdadU0"
    


    #git clone  https://github.com/shenlab-sinai/ngsplot.git

    #ngsplotdb.py install -y ngsplotdb_mm10_75_3.00.tar.gz


    # This is for the treatment file, which I put second in the channel in this nf script
    # creating an enrichment plot of the transcription start site (tss) using ngsplot
     ngs.plot.r \
     -G mm10 \
     -R tss \
     -C "${bam[1]}" \
     -O klf4.tss \
     -T klf4 \
     -L 3000 \
     -FL 50 \
     -RR 70


    # This is for the control file which I put first in the nf script
    # creating the same graphs for the input gene

    ngs.plot.r \
     -G mm10 \
     -R tss \
     -C "${bam[0]}" \
     -O input.tss \
     -T Input \
     -L 3000 \
     -FL 50 \
     -CO blue2 \
     -RR 70

    """

}
*/






// sorting the bam files to see if ngs plot takes them

process sort_with_samtools {
    // for some reason having the memory at 11GB here ensures that the large bam files get sorted into one file without multiple temp files instead
    memory '200 GB'
    publishDir './sorted_indexed_bams', mode:'copy', pattern:'*'

    conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/samtools_rj_env.yml'
    cache true

    input:

    path (bams)
    val (bam_names)

    output:

    path "${bam_names}.sorted.bam", emit: sorted_bam_ch
    path "*.bai", emit: indexed_bam_ch

    script:

    sorted_bams_out = "${bam_names}.sorted.bam"
    //indexed_bams_out = "${bam_names}.indexed.bam.bai"
   
    """
    #!/bin/env bash

    samtools sort \
    -o "${sorted_bams_out}" \
    -O bam \
    "${bams}"

    samtools index \
    -b \
    "${bams}" 



    """
}




// making ngs_plot process for jan

process ngs_plot_jan {

    // since i cant fix the zip file error thing yet i will just tell the process to ignore any errors and continue
    //errorStrategy 'ignore'


    // lets add more memory to see if it completes faster
    //memory '150 GB'
    // this was the old yml file with ngsplot but I found and made a new one
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/ngsplot_rj_env.yml'
    
    // better ngsplot yml file with hg38 genome but no other genome
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/ngsplot_w_hg38_rj_env.yml'
    // for some reason the yml file is not good so i have to give the actual conda environment i made
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/ngsplot_db_hg38'
    // these other two conda environments do not use ngsplot properly
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/base-r_rj_env.yml'
    //conda '/lustre/fs4/home/rjohnson/conda_env_files_rj_test/ngsplot_rj_env.yml'
    publishDir './Visualizations', mode: 'copy', pattern: '*.pdf'

    input:

    path(bam_files)
    val(bam_names)
    
    // dont need any of this anymor but will keep it if i find a new problem with the new yml ngsplot file
    //path(ngsplot_tar_file)
    //path(ngsplot_r_update)
    //path(bai_files)
    //path bai

    output:
    
    path "*.pdf"
    
    script:


    """
    #!/usr/bin/env bash

    # for ngsplot you need to download the most recent ngsplot zip file from their google drive; link on their github
    # ill see if i can just put ngsplot zip file in the bin
    ########################
    # I got it to work by downloading the zip file from ngs plot google dirve
    # https://drive.google.com/drive/folders/0B1PVLadG_dCKN1liNFY0MVM1Ulk?resourcekey=0-8RNVZduOvcSGsafKYa0AkA
    # the version i used is 2..61

    
    # a package r-shortread will be needed with this new setup
    # i can probably just comment all of this out since it is saved in the yml file
    # using this version of r 4.2.3
    # a package r-shortread will be needed with this new setup
    #conda install conda-forge::r-base=4.2.3 -y 
    #conda config --add channels bioconda
    #conda config --add channels conda-forge
    #conda config --set channel_priority strict
    #conda install bioconductor-shortread -y

    # i want to use hg38
    # conda install bioconda::r-ngsplotdb-hg38

    # now need BSgenome
    #conda install bioconductor-bsgenome -y 

    # this ensures i can use the zip feature without any problems
    #conda install conda-forge::zip -y
    
    # lets try the next steps with a conda install of domc and catools so it is saved in the exported yml file
    #conda install conda-forge::r-domc -y
    #conda install conda-forge::r-catools -y
    
    ###### dont need this below since the above works ######## installing it using conda will save it in the yml file ######
    # is a bioconductior package so update biocmanager (domc)
    # Rscript -e 'BiocManager::install("doMC")' 
    # Rscript -e 'BiocManager::install("caTools")'

    # these are the last dependencies that needed to be set up for ngsplot to work
    ##################################


    ############# Dont need this anymore since i have a new yml file ################
    # first extract the ngsplot tar file
    # or maybe i have to put the tar file in a channel
    #tar -xvf "\${ngsplot_tar_file}"
    # i needed to then remove the old ngs.plot.r file and add the new one that removed the zip commands that kept causing an error
    #rm ngsplot/bin/ngs.plot.r
    #cp "\${ngsplot_r_update}" ngsplot/bin

    # i want to see where ngs plot is located and set the environment variable
    #pwd
    #ls -l
    ###################################################

    # ensuring ngsplot will work anywhere in this environment
    #export NGSPLOT=\$(pwd)/ngsplot
    
    # when using this second export command it works and can find ngs.plot.r in the terminal
    #export PATH=\$(pwd)/ngsplot/bin:\$PATH

    

    #I dont need this yet
    #ngsplotdb.py install -y ngsplotdb_mm10_75_3.00.tar.gz

    #ngsplotdb.py install -y  ngsplotdb_hg19_75_3.00.tar.gz   


    # Debugging
    #echo "PATH: \$PATH"
    #which zip
    # This is for the treatment file, which I put second in the channel in this nf script
    # creating an enrichment plot of the transcription start site (tss) using ngsplot
    
    


    ngs.plot.r \
     -G hg38 \
     -R tss \
     -C "${bam_files}" \
     -O "${bam_names}.tss" \
     -T "${bam_names}" \
     -L 6000 \
     -FL 50 \
     -RR 70 \
     -SC global \
     -P 12 \
     -MQ 20 \
     -SE 0 \
     -LOW 0 \
     -GO total



    """

}

/*
workflow original_workflow{

    // i want to put this in the workflow section to see if it can still work
    params.mouse_ref_mm10='/gpfs/home/rj931/projects/chip_seq_practice/nextflow_workstation/mouse_genome_grcm38/mm10.fa.gz'

    mm10_ref=Channel.fromPath(params.mouse_ref_mm10, checkIfExists: true)

    mouse_ref_genome_gtf()
    //getting the ref and gtf channels
    mouse_ref_genome = mouse_ref_genome_gtf.out.ref_genome

    mouse_gtf_genome = mouse_ref_genome_gtf.out.gtf_file

    mouse_ref_genome.view()
    mouse_gtf_genome.view()


    // this can be over ridden by setting the --se_reads parameter for single end reads
    // this is NOT control but igg. this experiment does not have the control just the treatment
    params.se_reads_control = '../chip_fastqs/control/chip_*.fastq'
    params.se_reads_treatment = '../chip_fastqs/treatment/chip_*.fastq'
    Channel.fromPath(params.se_reads_control, checkIfExists: true)
            .map{file -> tuple('control', file)}
            .set{SE_reads_control }
    SE_reads_control.view()
    
    Channel.fromPath(params.se_reads_treatment, checkIfExists: true)
            .map{file -> tuple('control', file)}
            .set{SE_reads_treatment}
    
    SE_reads_control
        .map { file -> file[1].baseName}
        .set {filenames_control}

    filenames_control.view()

    SE_reads_treatment
        .map { file -> file[1].baseName}
        .set {filenames_treatment}
    
    
    // I want to try getting the srr files using nextflow
    //params.sra_txt = '../SRR_Acc_List.txt'
    //SRR_txt_list = Channel.fromPath(params.sra_txt, checkIfExists: true).splitText().collect() //.view()
    //SRR_txt_list.view()
    //ids = SRR_txt_list
    //println ids
    // now that I have my SRR list, I can put it into the channel that works with sra

    //sra_acc = Channel.fromSRA(SRR_txt_list, apiKey:'a5b95701789e19a261c98a7fcdb920f46f08', cache:true)

    //sra_acc.view()
    //Channel.fromSRA()

    //control_fastq =  channel.fromPath('../chip_fastqs/goat_igg.fastq')

    fastp_se(SE_reads_control, SE_reads_treatment, filenames_control, filenames_treatment )

    filt_files_tuple = fastp_se.out.filtered_fastq_c_t
    filt_files_tuple.view()
    
    
    
    // get the base names of the filtered files
    filt_files_tuple
        .map {file -> file.baseName}
        .set {filt_basenames}
    filt_basenames.view()

    
    html_fastp_file = fastp_se.out.fastp_html.collect() 

    //filt_files_tuple.view()
    html_fastp_file.view()

    //bwa_index(mouse_ref_genome)
    
    // other version of mouse ref from above
    bwa_index(mm10_ref)

    // now taking those reference index files and giving them their own channel

    index_files = bwa_index.out.ref_index_files

    index_files.view()
        //.tap {amb_index, ann_index, bwt_index, pac_index, sa_index}

   


    //bwa_alignment(filt_files_tuple, mouse_ref_genome, index_files, mouse_gtf_genome, filt_basenames)

    // using the mm10 genome
    bwa_alignment(filt_files_tuple, mm10_ref, index_files, mouse_gtf_genome, filt_basenames)

    sam_output_files_ch = bwa_alignment.out.sam_files_tuple_ch

    sam_output_files_ch.view()

    

    // give the basename of the sam files.

    sam_output_files_ch
        .map{file -> file.baseName}
        .set{sam_basenames}

    samtools_convert(sam_output_files_ch, sam_basenames)

    sorted_bams = samtools_convert.out.sorted_bam_tuple_ch
    sorted_bams.view()

    big_wig_tuple_ch = samtools_convert.out.bigwig_c_t_tuple

    bai_files = samtools_convert.out.sorted_bai_ch
    // making the basename again from the sorted bams channel

    sorted_bams
        .map{file -> file.baseName}
        .set{bams_basename}

    big_wig_tuple_ch
        .map {file -> file.baseName}
        .set {big_wig_basename}

    peak_calling_mouse(sorted_bams, bams_basename )

    called_peaks_ch = peak_calling_mouse.out.peak_files
    called_peaks_ch.view()


    make_peak_count_table(big_wig_tuple_ch, big_wig_basename, called_peaks_ch )

    peak_table_ch = make_peak_count_table.out.peak_count_table_ch
    peak_table_ch.view()


    ngs_plot(sorted_bams, bai_files)
}*/
/*
params.bam_4u='/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_4U_TR1_TR2_merged_sub100.bam'
params.bam_20u='/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_20U_TR1_TR2_merged_sub100.bam'
params.bam_100u='/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_100U_TR1_TR2_merged_sub100.bam'
params.bam_500u='/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_500U_TR1_TR2_merged_sub100.bam'
*/

workflow {

    // first get the bam files to use as input

    params.ref_genome_files = file('/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/*')
    
    //making a  channel for the reference genome files
    Channel.value(params.ref_genome_files)
            .set{ref_genome_files} // genome_dict, genome_fai, genome_fa, genome_xml i looked at .view to see what order the files entered then i put them in a channel correctly
    
    // the files should be in this path /rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/
    // only run macs2 on the four bam files and look at a comparison of all the peaks 

    // you can call each of these parameters in the command line to change the file from these defualts
    //params.bam_4u = '/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_4U_TR1_TR2_merged_sub100.bam'
    //params.bam_20u = '/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_20U_TR1_TR2_merged_sub100.bam'
    //params.bam_100u = '/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_100U_TR1_TR2_merged_sub100.bam'
    //params.bam_500u = '/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_500U_TR1_TR2_merged_sub100.bam'
    
    //params.bam_test = file('/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_500U_TR1_TR2_merged.bam')
    
    // I want to get all four of Jan's files in one parameter
    params.bams_sub100 = file('/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_*_sub100.bam')
    
    //CHOMP_GM12878_500U_TR1_TR2_merged.bam
    // using the peak_calling_jan process

    // making channels that has a name and the file attached
    

    // i have to make the params into a channel

    channel.value(params.bams_sub100)
            .flatten()
            .set{bam_sub100_files_ch}
    //bam_sub100_files_ch.view()
    
    // getting the basename of each file that matched the glob pattern
    
    channel.value(params.bams_sub100)
            .map{file -> file.baseName}
            .flatten()
            .set{bam_sub100_names_ch}
    //bam_sub100_names_ch.view()
    
    //test_bam = Channel.value(params.bam_test)
    //test_name = test_bam.map{file -> file.baseName}
    // bam_4u might have zero reads mapped to it
    peak_calling_jan(bam_sub100_files_ch, bam_sub100_names_ch)

    //peak_calling_jan.out.peak_files.view()
    //peak_calling_jan.out.summits_bed_file.view()

    // checking to see if having the bai file would help
    // it didnt work but having this set up is still good
    /*
    
    params.bai_files = file('/rugpfs/fs0/risc_lab/scratch/jhickling/CAD/CAD-C/GM_CADC/GM_CAD_Titration/CHOMP_GM12878_*_sub100.bam.bai')

    bai_files_ch = Channel.value(params.bai_files).flatten()

    bai_files_ch.view()

    Channel.value(params.bai_files)
            .map{file -> file.baseName}
            .flatten()
            .set{bai_names_ch}
    */  

    // lets use the summits.bed files instead from the peak_calling_jan process

    //peak_calling_jan.out.summits_bed_file.view()

    
    summits_bed_ch = peak_calling_jan.out.summits_bed_file
    summits_bed_ch.view()
    summits_bed_ch.map{file -> file.baseName}
                .flatten()
                .set{summits_name_ch}

    // let's try to do ngs plot 
    

    //lets sort the bam files then see if that works in ngs plot

    sort_with_samtools(bam_sub100_files_ch, bam_sub100_names_ch)

    sort_with_samtools.out.sorted_bam_ch.view()

    bam_sorted_files = sort_with_samtools.out.sorted_bam_ch

    bam_sorted_files.flatten()
                    .map{file -> file.baseName}
                    .set{sorted_bam_names_ch}

    

    // i need to put the ngsplot file in a channel

    ngsplot_tar=file('bin/ngsplot-2.61.tar')
    ngsplot_tar_file_ch = Channel.value(ngsplot_tar)

    // i need to get the updated ngs.plot.r script I made from the bin in the project directory
    ngs_plot_r_script = file('bin/ngs.plot.r')
    ngs_plot_r_update = Channel.value(ngs_plot_r_script)

    // keeping this for the future but not using it since i have a better ngsplot yml file
    //ngs_plot_jan(bam_sorted_files, sorted_bam_names_ch, ngsplot_tar_file_ch, ngs_plot_r_update)
    ngs_plot_jan(bam_sorted_files, sorted_bam_names_ch)
    
    // I need to come back to this and figure out how to get ngsplot to work
    
    //ngs_plot_jan(summits_bed_ch, summits_name_ch)
    
    //make sem_process

    ref_genome_files.view()
    

    sem_process(ref_genome_files, bam_sub100_files_ch, bam_sub100_names_ch)
    
}



