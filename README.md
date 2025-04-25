# Neutrophil RNA-Seq Pipeline: GM-CSF & TNF-α Response

This project analyzes transcriptomic responses of **human neutrophils** stimulated with **GM-CSF**, **TNF-α**, or left untreated, using a **Nextflow-based RNA-seq pipeline**.

---

## Biological Context

- **GM-CSF** reprograms neutrophils into a transcriptionally distinct, pro-migratory, neuro-immune-activated state
- Enriched GO terms include: `axonogenesis`, `neuron projection development`, `dendritic morphogenesis`
- **TNF-α** shows a milder transcriptional effect

---

## Pipeline Overview

Built using **Nextflow DSL2** + **Conda**, the pipeline includes:

| Step            | Tool            |
|-----------------|-----------------|
| Quality Control | FastQC          |
| Alignment       | STAR            |
| Quantification  | featureCounts   |
| DE Analysis     | DESeq2 (R)      |
| Enrichment      | clusterProfiler |

---

```bash
nextflow run main.nf -profile standard

