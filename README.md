# SEATRAC TB Hackday

Wednesday, December 13, 2023
9 AM to 3 PM (Pacific Time)

UW SLU F106/F107 and Zoom (University of Washington South Lake Union Campus)
([Zoom link, registration required](https://urldefense.com/v3/__https://washington.zoom.us/meeting/register/tJ0sf-GtrjMtHtAekQf7Ew-m4_oYOvm78vVB__;!!GuAItXPztq0!nEXklIjP7jjBX1StCKX1p08_DtEWjHU6skNOBaO0ntWCECjyQD9WB8glSFwDTxCp98VbKIU_x7RSTnc2p_c$))

Goal: Collaborating in small groups to develop new analytical insights
from 4 published bulk and single-cell RNAseq datasets that profile NHP
host gene expression during Mtb infection and BCG vaccination.

Datasets: PMIDs: 37097292, 35483355, 37267955, 37390827

Program: Agenda here [link to agenda]. Short talks to introduce datasets, workshop tutorials on
analysis techniques, and open-ended hack time to plan and execute
analyses.

Breakfast included starting 8:30am, lunch included at noon, happy
hour snacks included starting 3pm.

Registration:
Follow link: https://washington.zoom.us/meeting/register/tJ0sf-GtrjMtHtAekQf7Ew-m4_oYOvm78vVB


DATASET(S) 

Foreman et al., 2023: bulk RNAseq from sorted cells in NHP Mtb challenge model 

Gideon et al., 2022: scRNAseq from granulomas in NHP Mtb challenge model 

Darrah et al., 2023: scRNAseq from BAL in NHP Mtb challenge (BCG route, correlates of protection) 

Liu et al., 2023: whole-blood bulk RNAseq from NHP Mtb challenge (BCG route and IV BCG dose, should match Darrah et al. study) 

 

============================================================================= 

Foreman TW, Nelson CE, Sallin MA, Kauffman KD, Sakai S, Otaizo-Carrasquero F, Myers TG, Barber DL. CD30 co-stimulation drives differentiation of protective T cells during Mycobacterium tuberculosis infection. J Exp Med. 2023 Aug 7;220(8):e20222090. doi: 10.1084/jem.20222090. Epub 2023 Apr 25. PMID: 37097292; PMCID: PMC10130742. 

PDF: https://rupress.org/jem/article-pdf/220/8/e20222090/1451373/jem_20222090.pdf 

MS: https://rupress.org/jem/article/220/8/e20222090/214054/CD30-co-stimulation-drives-differentiation-of 

DATA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227653 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228114 for the superset that contains both the macaque and mouse data for that study 

DATA SUMMARY: n=4 Rhesus macaques, Mtb challenge (40–80 CFU of Mtb-Erdman-mCherry and euthanized 6–7 wk after infection); bulk RNAseq from FACS sorted T cells 

 

Figure 1   

 

============================================================================= 

Gideon HP, Hughes TK, Tzouanas CN, Wadsworth MH 2nd, Tu AA, Gierahn TM, Peters JM, Hopkins FF, Wei JR, Kummerlowe C, Grant NL, Nargan K, Phuah JY, Borish HJ, Maiello P, White AG, Winchell CG, Nyquist SK, Ganchua SKC, Myers A, Patel KV, Ameel CL, Cochran CT, Ibrahim S, Tomko JA, Frye LJ, Rosenberg JM, Shih A, Chao M, Klein E, Scanga CA, Ordovas-Montanes J, Berger B, Mattila JT, Madansein R, Love JC, Lin PL, Leslie A, Behar SM, Bryson B, Flynn JL, Fortune SM, Shalek AK. Multimodal profiling of lung granulomas in macaques reveals cellular correlates of tuberculosis control. Immunity. 2022 May 10;55(5):827-846.e10. doi: 10.1016/j.immuni.2022.04.004. Epub 2022 Apr 27. PMID: 35483355; PMCID: PMC9122264. 

PAPER: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9122264/ 

We have performed high-throughput single-cell mRNA sequencing on 6 granulomas from 2 non-human primates at 4 weeks post infection with M. tuberculosis to understand cellular and molecular factors associated with bacterial control (n = 10,006 cells). 

We have performed high-throughput single-cell mRNA sequencing on 26 granulomas from 4 non-human primates at 10 weeks post infection with M. tuberculosis to understand cellular and molecular factors associated with bacterial control (n = 109,584 cells) 

 

DATA 4 WKS: https://singlecell.broadinstitute.org/single_cell/study/SCP1749/cellular-ecology-of-m-tuberculosis-granulomas-4-week-dataset#study-summary 

DATA 10 WKS: https://singlecell.broadinstitute.org/single_cell/study/SCP257/cellular-ecology-of-m-tuberculosis-granulomas-10-week-dataset#study-summary 

DATA SUMMARY: single-cell RNAseq (SeqWell); 6 granulomas from 2 NHP at 4 wks and 26 granulomas from 4 NHP at 10 weeks post-infection (includes CFU per granuloma) 

============================================================================= 

Darrah PA, Zeppa JJ, Wang C, Irvine EB, Bucsan AN, Rodgers MA, Pokkali S, Hackney JA, Kamath M, White AG, Borish HJ, Frye LJ, Tomko J, Kracinovsky K, Lin PL, Klein E, Scanga CA, Alter G, Fortune SM, Lauffenburger DA, Flynn JL, Seder RA, Maiello P, Roederer M. Airway T cells are a correlate of i.v. Bacille Calmette-Guerin-mediated protection against tuberculosis in rhesus macaques. Cell Host Microbe. 2023 Jun 14;31(6):962-977.e8. doi: 10.1016/j.chom.2023.05.006. Epub 2023 Jun 1. PMID: 37267955; PMCID: PMC10355173. 

https://www.nature.com/articles/s41586-019-1817-8.pdf 

https://www.nature.com/articles/s41586-019-1817-8 

DATA: https://singlecell.broadinstitute.org/single_cell/study/SCP796/prevention-of-mycobacterium-tuberculosis-infection-and-disease-in-nonhuman-primates-following-intravenous-bcg-vaccination?scpbr=the-alexandria-project 

 

DATA SUMMARY: single-cell; n= 15 Rhesus macaques with BCG vaccination and Mtb challenge; 3 individuals per group for AE, IDhigh, IDlow, IV, Naïve-controls; BAL collected at Weeks 13 and 25 prior to challenge at Week 26; stimulated and unstimulated conditions; 60 samples, 1000 – 5000 cells each. 

Mycobacterium tuberculosis (Mtb) is the leading cause of death from infection worldwide. Intradermal (ID) vaccination with BCG has variable efficacy against pulmonary tuberculosis, the major cause of mortality and disease transmission. Here we show that the route and dose of BCG vaccination alters circulating and lung resident T cells and subsequent protection against Mtb challenge in nonhuman primates (NHP). NHP immunized with BCG by the intravenous (IV) route induced substantially higher antigen-specific CD4 (Th1 or Th17) and CD8 responses in blood, spleen, bronchoalveolar lavage (BAL), and lung lymph nodes compared to the same BCG dose administered by ID or aerosol (AE) routes. Moreover, IV immunization was the only route that induced a high frequency of antigen-specific tissue resident T cells in lung parenchyma. Six months after BCG vaccination, NHP were challenged with virulent Mtb. Strikingly, 9 of 10 NHP that received BCG IV were highly protected, with 6 NHP showing no detectable infection as determined by PET CT imaging, mycobacterial growth, pathology, granuloma formation, or de novo immune responses to Mtb-specific antigens. The finding that BCG IV prevents or significantly limits Mtb infection in NHP has important implications for vaccine development and provides a model for determining immune correlates and mechanisms of protection against TB. 

 

============================================================================= 

Liu YE, Darrah PA, Zeppa JJ, Kamath M, Laboune F, Douek DC, Maiello P, Roederer M, Flynn JL, Seder RA, Khatri P. Blood transcriptional correlates of BCG-induced protection against tuberculosis in rhesus macaques. Cell Rep Med. 2023 Jul 18;4(7):101096. doi: 10.1016/j.xcrm.2023.101096. Epub 2023 Jun 29. PMID: 37390827; PMCID: PMC10394165. 

https://www.sciencedirect.com/science/article/pii/S266637912300215X?via%3Dihub#sec4.1 

 

DATA IV BCG DOSE STUDY (167 samples, 34 NHP): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218270 

DATA BCG ROUTE STUDY (144 samples, 36 NHP):  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218157 

 

DATA SUMMARY: Darah et al.  IV BCG dose and BCG route study matching; Pre BCG and Day 2, Wk 4, Wk 12 for IV BCG dose study;  Pre, Day 2, Wk2 and Wk12 for route study; bulk RNAseq from whole blood; includes protection data for each individual NHP 