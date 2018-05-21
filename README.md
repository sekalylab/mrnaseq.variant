## Variant calling based on RNA-Seq
RNA-Seq variant calling pipeline according to [GATK Best practices](https://www.broadinstitute.org/gatk/guide/article?id=3891).

#### Run mRNA-Seq variant pipeline
![picture](img/mrnaseq-variant.pipeline.png)
```bash
bash mRNA.variant_master.sh -d {raw_directory}

arguments:
d=[d]irectory with raw data (directory; required)  
g=directory with the reference [g]enome  
    accepted values: GRCh98, Mmul_8  
e=[e]mail address
```


