# Piano di analisi e struttura validazione TSO500

## Coverage

### Coverage DNA

1. **SNV-Indel: StitchedRealigned**
I bam sono in `/data/novaseq_results/research/validation/bam/dna/snv`

```
ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210715_A01423_0008_AH35CWDRXY/Logs_Intermediates/StitchedRealigned/*/*.bam* ./

ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210729_A01423_0009_AH33WGDRXY/Logs_Intermediates/StitchedRealigned/*/*.bam* ./

ln -s /data/novaseq_results/211022_A01423_0010_AHGYFYDRXY_sequencerOutput/Logs_Intermediates/StitchedRealigned/*/*.bam* ./
```

Aggiunti bam della run 211111_A01423_0011_AHH2Y2DRXY in `/data/novaseq_results/research/validation/bam/dna/snv/211111_A01423_0011_AHH2Y2DRXY`

```
ln -s /data/novaseq_results/211111_A01423_0011_AHH2Y2DRXY/Logs_Intermediates/StitchedRealigned/*/*.bam* ./
```


I risultati sono in `/data/novaseq_results/research/validation/bam/results/coverage/dna/snv`



2. **CNV: DnaRealignment**
I bam sono in `/data/novaseq_results/research/validation/bam/dna/cnv`

```
ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210715_A01423_0008_AH35CWDRXY/Logs_Intermediates/DnaRealignment/*/*.bam* ./

ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210729_A01423_0009_AH33WGDRXY/Logs_Intermediates/DnaRealignment/*/*.bam* ./

ln -s /data/novaseq_results/211022_A01423_0010_AHGYFYDRXY_sequencerOutput/Logs_Intermediates/DnaRealignment/*/*.bam* ./
```
Aggiunti bam della run 211111_A01423_0011_AHH2Y2DRXY in `/data/novaseq_results/research/validation/bam/dna/cnv/211111_A01423_0011_AHH2Y2DRXY`


```
ln -s /data/novaseq_results/211111_A01423_0011_AHH2Y2DRXY/Logs_Intermediates/DnaRealignment/*/*.bam* ./
```


Scritp utilizzato:

`/data/novaseq_results/research/validation/bam/results/coverage/dna/cnv/cnv_coverage.sh`


I risultati sono in `/data/novaseq_results/research/validation/bam/results/coverage/dna/cnv`

### Coverage RNA

1. **RNA MarkDuplicates**

I bam sono in `/data/novaseq_results/research/validation/bam/rna`

```
ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210715_A01423_0008_AH35CWDRXY/Logs_Intermediates/RnaMarkDuplicates/*/*.bam* ./

ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210729_A01423_0009_AH33WGDRXY/Logs_Intermediates/RnaMarkDuplicates/*/*.bam* ./

ln -s /data/novaseq_results/211022_A01423_0010_AHGYFYDRXY_sequencerOutput/Logs_Intermediates/RnaMarkDuplicates/*/*.bam* ./
```

Aggiunti bam della run 211111_A01423_0011_AHH2Y2DRXY in `/data/novaseq_results/research/validation/bam/rna/211111_A01423_0011_AHH2Y2DRXY`

```
ln -s /data/novaseq_results/211111_A01423_0011_AHH2Y2DRXY/Logs_Intermediates/RnaMarkDuplicates/*/*.bam* ./
```

I risultati sono in `/data/novaseq_results/research/validation/bam/results/coverage/rna`

Scritp utilizzato:

`/data/novaseq_results/research/validation/bam/results/coverage/rna/rna_coverage.sh`




## Vcf annotation

### Vcf DNA

1. **SNV-Indel**

I vcf sono in `/data/novaseq_results/research/validation/vcf/dna/snv`

```
ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210715_A01423_0008_AH35CWDRXY/Results/*/*/*DNA_Merged*.vcf ./

ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210729_A01423_0009_AH33WGDRXY/Results/*/*/*DNA_Merged*.vcf ./

ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210729_A01423_0009_AH33WGDRXY/Results/*/*DNA_Merged*.vcf ./

ln -s /data/novaseq_results/211022_A01423_0010_AHGYFYDRXY_sequencerOutput/Results/*/*/*DNA_Merged*.vcf ./

ln -s /data/novaseq_results/211022_A01423_0010_AHGYFYDRXY_sequencerOutput/Results/*/*DNA_Merged*.vcf ./
```


Aggiunta la run 211111_A01423_0011_AHH2Y2DRXY in `/data/novaseq_results/research/validation/vcf/dna/snv/211111_A01423_0011_AHH2Y2DRXY`

```
ln -s /data/novaseq_results/211111_A01423_0011_AHH2Y2DRXY/Results/*_DNA*/*MergedSmallVariants.genome.vcf ./
```

I risultati sono in `/data/novaseq_results/research/validation/vcf/results/dna/snv`

Scritp utilizzato:

`/data/novaseq_results/research/validation/vcf/results/dna/snv/snv_annotation.sh`


2. **CNV**

I vcf sono in `/data/novaseq_results/research/validation/vcf/dna/cnv`

```
ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210715_A01423_0008_AH35CWDRXY/Results/*/*_DNA/*_DNA_CopyNumberVariants.vcf ./

ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210729_A01423_0009_AH33WGDRXY/Results/*/*_DNA/*_DNA_CopyNumberVariants.vcf ./

ln -s /data/novaseq_results/211022_A01423_0010_AHGYFYDRXY_sequencerOutput/Results/*/*_DNA/*_DNA_CopyNumberVariants.vcf ./
```

Aggiunta la run 211111_A01423_0011_AHH2Y2DRXY in `/data/novaseq_results/research/validation/vcf/dna/cnv/211111_A01423_0011_AHH2Y2DRXY`

```
ln -s /data/novaseq_results/211111_A01423_0011_AHH2Y2DRXY/Results/*_DNA*/*CopyNumberVariants.vcf ./
```

I risultati sono in `/data/novaseq_results/research/validation/vcf/results/dna/cnv`

Scritp utilizzato:

`/data/novaseq_results/research/validation/vcf/results/dna/cnv/cnv_annotation.sh`



### Vcf RNA

I vcf sono in `/data/novaseq_results/research/validation/vcf/rna`

```
ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210715_A01423_0008_AH35CWDRXY/Results/*/*_RNA/*.vcf ./

ln -s /data/novaseq/Diagnostic/NovaSeq/Results/210729_A01423_0009_AH33WGDRXY/Results/*/*_RNA/*.vcf ./

ln -s /data/novaseq_results/211022_A01423_0010_AHGYFYDRXY_sequencerOutput/Results/*/*_RNA/*.vcf ./
```

Aggiunta la run 211111_A01423_0011_AHH2Y2DRXY in `/data/novaseq_results/research/validation/vcf/rna/211111_A01423_0011_AHH2Y2DRXY`

```
ln -s /data/novaseq_results/211111_A01423_0011_AHH2Y2DRXY/Results/*_RNA*/*.vcf ./
```

I risultati sono in `/data/novaseq_results/research/validation/vcf/results/rna`

Scritp utilizzato:

`/data/novaseq_results/research/validation/vcf/results/rna/rna_annotation.sh`





