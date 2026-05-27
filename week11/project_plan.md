# Project Plan: Hippocampal Estradiol, Chromatin State, and Stress-Induced Memory Vulnerability

_Short project proposal for re-analysis of public CUT&RUN data from [Hokenson et al., Neuron 2026.](https://www.cell.com/neuron/fulltext/S0896-6273(25)00993-6)_

## Topic

This project will re-analyze public CUT&RUN data from Hokenson et al., *Neuron* 2026, to evaluate whether hippocampal estradiol state and estrogen receptor (ER$\beta$) activation are associated with chromatin states that predict vulnerability to multiple acute concurrent stresses (MACS). The focus is on H3K4me3 and H3K27me3 profiles in mouse hippocampus across sex, estrous-cycle phase, MACS exposure, and ER$\beta$ agonism.

These are (imo) the crucial results of the publication:

- MACSs induced memory disturbances in male and proestrus female mice whereas estrus females were protected from memory disturbances (Figure 1 A-E)
- High estradiol levels in hippocampus were required for MACSs induced memory problems (Figure 3 B-D)
- Stimulating ER$\beta$ in estrus females led to MACS-induced memory problems and changes in chromatin (Figure 7 A-B)

## Background

Post-traumatic stress disorder is more common in women than in men, even after accounting for trauma type. Hokenson et al. use a MACS mouse model to study sex-specific vulnerability to stress-induced memory disturbances. In their model, MACS disrupts memory in males and proestrus females, whereas estrus females are relatively protected.

Estradiol differs across both biological sex and estrous-cycle phase. In mice, the estrous cycle lasts approximately 4-5 days and includes proestrus, estrus, metestrus, and diestrus. Circulating estradiol is high during proestrus and low in males, but hippocampal estradiol shows a different pattern: it is reported to be high in males and proestrus females and lower in estrus females.

Estrogen receptors are ligand-activated transcription factors. The two main nuclear receptor subtypes are ER$\alpha$ and ER$\beta$. After ligand binding, they dimerize and can bind estrogen response elements, but many estrogen-responsive genes lack canonical estrogen response elements and may instead be regulated through interactions with DNA-bound transcription factors such as AP-1 and Sp1.

## Data

Primary dataset: [Hokenson et al., DOI: 10.1016/j.neuron.2025.12.037](https://doi.org/10.1016/j.neuron.2025.12.037). Public data are available through [GEO GSE313359](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313359) and [BioProject PRJNA1379317](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1379317).

The assay is CUT&RUN-seq on *Mus musculus* hippocampus, profiling H3K4me3 and H3K27me3.
The analysis can either start from processed bigWig tracks or raw reads in FASTQ format.

| Condition | Experimental block | H3K4me3 libraries | H3K27me3 libraries |
|---|---|---|---|
| Naive male | baseline sex comparison | n = 3 | n = 3 |
| Naive proestrus female | high-estradiol estrous phase | n = 3 | n = 3 |
| Naive estrus female | low-estradiol estrous phase | n = 3 | n = 3 |
| Estrus female control-vehicle | perturbation control | n = 4 | n = 4 |
| Estrus female MACS-vehicle | stress without ER$\beta$ agonist | n = 3 | n = 2 |
| Estrus female control-DPN | ER$\beta$ agonist without stress | n = 4 | n = 4 |
| Estrus female MACS-DPN | stress with ER$\beta$ agonist | n = 3 | n = 3 |
| Total | 7 biological conditions across two marks | 23 libraries | 22 libraries |

## Proposed Analysis

The authors claim that high estradiol levels generate permissive chromatin states which then allow MACSs to induce memory problems.

I want to examine this statement with the following analyses:

1. Compare H3K4me3 and H3K27me3 peak patterns across naive males, naive proestrus females, and naive estrus females to assess differences in chromatin "permissivity".

2. Test whether DPN treatment moves estrus females toward a proestrus-like chromatin state, relative to estrus control-vehicle samples.

3. Test whether MACS changes chromatin differently in estrus females with versus without ERβ activation (DPN treatment).

4. Identify genes associated with permissive chromatin in males and proestrus females, and do some enrichment analysis to see if they are associated to plasticity, memory, or stress-relevant pathways.

5. Perform motif enrichment on H3K4me3 peaks to test whether ER, AP-1, or Sp1 motifs are enriched in permissive states.

## Links

- Article: [https://doi.org/10.1016/j.neuron.2025.12.037](https://doi.org/10.1016/j.neuron.2025.12.037)
- GEO series: [GSE313359](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313359)
- BioProject / SRA: [PRJNA1379317](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1379317)

