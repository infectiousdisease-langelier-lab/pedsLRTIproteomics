This is the repository for the metadata, proteomic data, and code necessary for data analysis and figure generation for the manuscript titled "Tracheal aspirate and plasma proteomics reveals the local and systemic host immune response to severe pediatric lower respiratory viral infections." Link to manuscript will be added once accepted/published. 

## File descriptions

1) "ViralLRTIProteomics_metadata.csv" includes limited de-idenfitied metadata for the cohort.
2) "TAdf.csv" includes the normalized, log2-transformed relative protein aptamer concentrations measured in the tracheal aspirate from all subjects, as provided by SomaScan®. The rows are subjects and the columns are proteins.
3) "plasmadf.csv" includes the normalized, log2-transformed relative protein aptamer concentrations measured in the plasma from all subjects, as provided by SomaScan®. The rows are subjects and the columns are proteins.
4) "Uniprotein key.csv" is the key linking the protein name in files #2 and #3 with the UniProtein ID, which is needed for the enrichment analyses performed in the paper. This is provided by SomaScan®
5) "ViralLRTIProteomics_code.R" is the code used for all of the analysis and figures. The version of R used is R v4.3.2.
