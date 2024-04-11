#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.run_type}" == "r2d2" ]] || [[ "${params.run_type}" == "raven" ]] || [[ "${params.run_type}" == "studio" ]]; 
      then
        cd ${params.image_folder}
        if [[ ! -f diffbind-3.12.0.sif ]] ;
          then
            singularity pull diffbind-3.12.0.sif docker://index.docker.io/mpgagebioinformatics/diffbind:3.12.0
        fi
        if [[ ! -f macs-2.2.9.1.sif ]] ;
          then
            singularity pull macs-2.2.9.1.sif docker://index.docker.io/mpgagebioinformatics/macs:2.2.9.1
        fi
    fi


    if [[ "${params.run_type}" == "local" ]] ; 
      then
        docker pull mpgagebioinformatics/diffbind:3.12.0
        docker pull mpgagebioinformatics/macs:2.2.9.1
    fi

    """

}

process make_samplesheet {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val bowtie_out
    val macs2_out
    val samplesheet

  when:
      (  ! file("${params.project_folder}/diffbind_sample_sheet.csv").exists() )

  script:
    """
#!/usr/local/bin/python
import pandas as pd
import numpy as np
import openpyxl
import os

bowtie_out='${bowtie_out}'
macs2_out='${macs2_out}'
tmp = '${params.project_folder}/tmp/'
sampleFile= '${params.samplestable}'
sampleCSV= '${samplesheet}'

if not os.path.exists(tmp):
    os.makedirs(tmp)

book = openpyxl.load_workbook(sampleFile)
samples_in = pd.read_excel(sampleFile, sheet_name= 'samples', engine="openpyxl")
for col in samples_in.columns:
    if samples_in[col].dtype == 'object':
        samples_in[col] = samples_in[col].str.replace('\\n', '').str.replace('"', '')
if 'input' in book.sheetnames:
    input_match = pd.read_excel(sampleFile, sheet_name = 'input', engine = 'openpyxl')
    if len(input_match) >= 1:
        for col in input_match.columns:
            if input_match[col].dtype == 'object':
                input_match[col] = input_match[col].str.replace('\\n', '').str.replace('"', '')
        if not pd.isnull(input_match['Input Sample']).all():
            control_input = True
        else:
            control_input = False
    else:
        input_match = pd.DataFrame(columns = ['ChIP Sample', 'Input Sample'])
        control_input = False
else:
    input_match = pd.DataFrame(columns = ['ChIP Sample', 'Input Sample'])
    control_input = False
samples_in
input_match.head()
# subset samples to only chip samples, i.e. remove samples from samples_in which are in "input_match['Input Sample']"
samples_out = samples_in[~samples_in['Sample'].isin(input_match['Input Sample'])]
samples_out = samples_out[['Sample', 'Group', 'Replicate']]
samples_out.rename(columns = {'Sample': 'SampleID', 'Group': 'Treatment'}, inplace = True)
samples_out['bamReads'] = bowtie_out + samples_out['Treatment'] + '.Rep_'+ samples_out['Replicate'].astype(str) + '.md.bam'
if control_input:
    samples_in['SampleID'] = samples_in['Group'] + '.Rep_'+ samples_in['Replicate'].astype(str) 
    samples_dict = samples_in.set_index('Sample').to_dict()['SampleID']
    input_match.rename(columns = {'ChIP Sample': 'SampleID', 'Input Sample': 'ControlID'}, inplace = True)
    # merge input match table
    samples_out = pd.merge(samples_out, input_match, how = 'outer', on = 'SampleID')
    samples_out['bamControl'] = np.nan
    for index, sample in samples_out.iterrows():
        if not pd.isnull(sample['ControlID']):
            samples_out.loc[index, 'bamControl'] = bowtie_out + samples_dict[sample['ControlID']] + '.md.bam'
    with open(tmp + '/input.txt', 'w') as the_file:
        the_file.write('yes\\n')
else:
    with open(tmp+ '/input.txt', 'w') as the_file:
        the_file.write('no\\n')
samples_out
samples_out['Peaks'] = macs2_out + samples_out['Treatment'] + '.Rep_'+ samples_out['Replicate'].astype(str) + '_peaks.xls'
samples_out['PeakCaller'] = 'macs'
samples_out
samples_out.to_csv(sampleCSV, sep=',', index = False)


    """
}


process diffbind_R {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.project_folder}/diffbind3_output/consensus_peaks.tsv").exists() ) 
  input:
    val samplesheet
  
  script:
    """
#!/usr/bin/Rscript

print(Sys.getenv("R_LIBS_USER")) # update
print("loading libraries")
library(DiffBind)
library(openxlsx)

TMP = '/workdir/tmp/'

dir.create("/workdir/diffbind3_output/", recursive = TRUE, showWarnings = FALSE)
setwd("/workdir/diffbind3_output/")

samples = read.csv("${samplesheet}")
# check how many peaks are in each samples
for(s in 1:nrow(samples)){
  filein=gsub('_peaks.xls', '_summits.bed', samples[s,'Peaks'])
  con <- file(gsub('${params.project_folder}', '/workdir/', filein))
  print(con)
  samples[s, 'nPeaks'] = length(readLines(con))
  close(con)
}
# remove samples without peaks
samples = subset(samples, nPeaks > 0)
# load peakset
hist<-dba(sampleSheet=samples)
# start QC pdf
pdf("general_QC.pdf")
plot(hist)
# count reads
hist<-dba.count(hist, bParallel=FALSE)
plot(hist)
# check fraction of reads in peaks and save results
info <- dba.show(hist)
print(info)
if ('FRiP' %in% names(info)){
  libsizes <- cbind(LibReads=info[, 'Reads'], FRiP=info[, 'FRiP'], PeakReads=round(info[, 'Reads'] * info[,'FRiP']))
  rownames(libsizes) <- info[, 'ID']
  write.xlsx(libsizes, "libsizes.xlsx", row.names = TRUE)
}
dba.plotPCA(hist, attributes=DBA_TREATMENT, label=DBA_ID)
# normalize data
hist<-dba.normalize(hist)
plot(hist)
# close QC pdf
dev.off()
# add contrast if there are at least 2 groups with at least 2 replicates
if (sum(table(samples[, 'Treatment']) > 1) > 1) {
  hist<-dba.contrast(hist, categories=DBA_TREATMENT, minMembers = 2)
  # differential binding analysis
  hist<-dba.analyze(hist, method=DBA_DESEQ2)
  # get list of all comparisons
  comparisons <- dba.show(hist, bContrasts=T)
  for(comp in 1:nrow(comparisons)){
    comp_name = paste0(comparisons[comp, 'Group'], '.vs.', comparisons[comp, 'Group2'])
    pdf(paste0(comp_name, '.pdf'))
    # get significant peaks
    hist_ <- dba.report(hist, contrast = comp)
    
    # various QC plots for differential binding analysis
    if (!is.null(hist_) && length(hist_@ranges) > 1){
      dba.plotPCA(hist, contrast=comp, method=DBA_DESEQ2, attributes=DBA_TREATMENT, label=DBA_ID)
      plot(hist,contrast=comp)
      dba.plotVolcano(hist, contrast = comp)
      pvals <- dba.plotBox(hist, contrast = comp)
      # heatmap values
      corvals <- dba.plotHeatmap(hist, correlations = FALSE, scale = 'row', contrast = comp)
      write.table(corvals, paste0('normalized_exression_of_significantly_changed_peaks.', comp_name, '.tsv'), sep = '\t', quote = F, row.names = F)
      write.xlsx(corvals, paste0('normalized_exression_of_significantly_changed_peaks.', comp_name, '.xlsx'), row.names = F)
      # save significant peaks
      write.table(hist_, paste0("significant.",comp_name , ".tsv"), sep="\t", row.names=FALSE, quote = FALSE) 
      write.xlsx(hist_, paste0("significant.",comp_name , ".xlsx"), row.names = FALSE)
    } else {
      print(paste("No significant values in", comp_name))
    }
    dba.plotMA(hist, method=DBA_DESEQ2, contrast = comp)
    dba.plotMA(hist, bXY=TRUE, contrast = comp)
    # save all peaks
    hist_all = as.data.frame(dba.report(hist, th = 1, contrast = comp))
    write.table(hist_all, paste0("all_peaks.",comp_name,".tsv"), sep="\t", row.names=FALSE, quote = FALSE)
    write.xlsx(hist_all, paste0("all_peaks.",comp_name,".xlsx"))
    write.table(hist_all[,1:3], paste0("all_peaks.",comp_name,".bed"), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
    dev.off()
  
  }
} else {
  print("
        
        Not enough replicates for differential analysis
        
        ")
}
consensus = as.data.frame(hist[['binding']])
consensus[, 'CHR'] = factor(consensus[,'CHR'], levels = 1:length(hist[['chrmap']]))
consensus[,'chrom_names'] = rep(hist[['chrmap']], as.numeric(table(consensus[, 'CHR'])))
write.xlsx(consensus, 'consensus_peaks.xlsx', row.names = FALSE)
write.table(consensus, 'consensus_peaks.tsv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
save.image(paste0(TMP,'/diffbind.Rdata'))
sessionInfo()

    """
}


process upload_paths {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
  """
    rm -rf upload.txt

    cd ${params.project_folder}/diffbind3_output/
    for f in \$(ls *.{pdf,xlsx}) ; do echo "diffbind \$(readlink -f \${f})" >>  upload.txt_ ; done  
    uniq upload.txt_ upload.txt 
    rm upload.txt_
  """
}


workflow images {
  main:
    get_images()
}

workflow samplesheet{
    // if ( 'bowtie_out' in params.keySet() ) {
    //     bowtie_out="${params.bowtie_out}"
    // }else {
    //     bowtie_out= '${params.project_folder}/bowtie2_output/'
    // }

    // if ( 'macs2_out' in params.keySet() ) {
    //     macs2_out="${params.macs2_out}"
    // }else {
    //     macs2_out= '${params.project_folder}/macs2_output/'
    // }
    def bowtie_out = params.containsKey('bowtie_out') ? params.bowtie_out : "${params.project_folder}/bowtie2_output/"
    def macs2_out = params.containsKey('macs2_out') ? params.macs2_out : "${params.project_folder}/macs2_output/"
    def samplesheet = params.containsKey('samplesheet') ? params.samplesheet : "${params.project_folder}/diffbind_sample_sheet.csv"

    make_samplesheet(bowtie_out,macs2_out,samplesheet)
}

workflow {
    def samplesheet = params.containsKey('samplesheet') ? params.samplesheet : "${params.project_folder}/diffbind_sample_sheet.csv"

    diffbind_R(samplesheet)
}

workflow upload {
  main:
    upload_paths()
}