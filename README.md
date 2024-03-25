# nf-diffbind


```
PROFILE=raven
nextflow run nf-diffbind -params-file ./nf-diffbind/params.slurm.json -entry images -profile ${PROFILE}  && \
nextflow run nf-diffbind -params-file ./nf-diffbind/params.slurm.json -entry samplesheet -profile ${PROFILE}  && \
nextflow run nf-diffbind -params-file ./nf-diffbind/params.slurm.json -profile ${PROFILE}
```