# Mapping protein targets' Uniprot, TTD, UMLS, Ensembl, DisGeNET ids:

This project includes a simple Java program that produces a [TSV file](https://github.com/iit-Demokritos/target-id-mapping/blob/main/target-mappings_latest.tsv?raw=true), listing thousand of known target protein available ids that could be found in target databases.

In particular, we start from retrieving the drug information included in the latest Therapeutic Target Database [1] (VERSION 7.1.01, RELEASED ON 2019.07.14) in a file. We then enrich the drug fields by querying the UMLS Metathesaurus vocabulary Database[2], using the MetamorphoSys tool  and by parsing the DisGenet_uniprot_mapping file, provided by DisGeNET [3].


## Licence & Required Citation
For any use of the drug-mappings.tsv file in your work, **a citation to the following paper is expected:**

*Aisopos, F., Paliouras, G. Comparing methods for drug–gene interaction prediction on the biomedical literature knowledge graph: performance versus explainability. BMC Bioinformatics 24, 272 (2023), [DOI](https://doi.org/10.1186/s12859-023-05373-2).*

target_id_mapping - NCSR Demokritos module Copyright 2021 Fotis Aisopos
The Java code and TSV file are provided **only for academic/research use and are licensed under the Apache License, Version 2.0 (the "License")**; you may not use this file except in compliance with the License. You may obtain a copy of the License at: https://www.apache.org/licenses/LICENSE-2.0 .

## Mapping TSV file data format

The resulting file ([target-mappings_latest.tsv](https://github.com/iit-Demokritos/target-id-mapping/blob/main/target-mappings_latest.tsv?raw=true)) includes a tab-separated entry for each target, including multiple ids that could be found and crossed-checked from the aforementioned sources.
For ids not found in none of the above sources, 'null' string is added. Multiple TTD ids/gene names for a specific target are separated with a comma separator(,).
An example of the format of the TSV data file is as follows:

```sh
CUI	Uniprot_Name	Uniprot_id	Gene_name	EnsembleGeneId	TTD_id	disgenet_gene_id
C1420086	NTCP_HUMAN	Q14973	SLC10A1	ENSG00000100652	T99189,T71907 6554
C1442506	MTR1B_HUMAN	P49286	MTNR1B	ENSG00000134640	T48268  4544
C1417177	MKRN1_HUMAN	Q9UHC7	MKRN1	ENSG00000133606	T86877  23608
C1416588	KCNJ6_HUMAN	P48051	KCNJ6	ENSG00000157542	T17721  3763
C1413043	CAH2_HUMAN	P00918	CA2	ENSG00000104267	T20401  760
...
```

## Java Project File Structure & running

The code includes the CreateTargetMappings class (main class), and the auxiliary TargetEntry class, representing a target object with the various properties.
To run the aforementioned Java project, it is obvious that we need to have access to the following sources:
- TTD (to download the targets' information file in raw format)
- Entrez Programming Utilities (E-utilities) API (query PUG for PubChem ids and obtain a token to query for a TGT and an API key)
and also include needed jar libraries in the CLASSPATH.

## References
[1]:  Y. X. Wang, S. Zhang, F. C. Li, Y. Zhou, Y. Zhang, R. Y. Zhang, J. Zhu, Y. X. Ren, Y. Tan, C. Qin, Y. H. Li, X. X. Li, Y. Z. Chen* and F. Zhu*. Therapeutic Target Database 2020: enriched resource for facilitating research and early development of targeted therapeutics. Nucleic Acids Research. 48(D1): D1031-D1041 (2020). PubMed ID: 31691823

[2]: Bodenreider, O. (2004). The unified medical language system (UMLS): integrating biomedical terminology. Nucleic acids research, 32(suppl_1), D267-D270.

[3]: Janet Piñero, Juan Manuel Ramírez-Anguita, Josep Saüch-Pitarch, Francesco Ronzano, Emilio Centeno, Ferran Sanz, Laura I Furlong.
The DisGeNET knowledge platform for disease genomics: 2019 update. Nucl. Acids Res. (2019) doi:10.1093/nar/gkz1021
