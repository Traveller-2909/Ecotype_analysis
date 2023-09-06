# Ecotype_analysis
Code and data of the paper titled "Ecotype diversification in northern pike (Esox lucius) as a mechanism for coexistence in an extreme brackish habitat"
The RMarkdown file "Ecotype_analysis_main" accessed by the RProject contains all relevant steps to reproduce the results presented in the paper. It consists of 4 main sections (headers)
      Life history phenotypes: Code to analyze elemental transects in sequentially formed aragonitic otoliths of pike to describe lifelong behavioral phenotypes
      Genotyping individual pike: Code to convert output from a STRUCTURE analysis into genotype assignment probabilities and assign discrete genotypes
      Growth analysis of life history phenotypes and genotype: Code to perform age-specific analyses of otolith growth increments via GLMM models and lifelong analyses of otolith radius at age via von Bertalanffy
                                                               growth curves
      Phenotype - genotype matching: Code to match observed lifelong behavioral phenotypes to underlying genotypes in order to distinguish ecotypes
      
The folder "Data" contains all data relevant to the analysis, with subfolders denoting which general part of the analysis the data belongs to

The RMarkdown file is currently not intended to be knitted as HTML, but rather serve as a series of code chunks that can be run independently to generate the necessary data tables, figures and tests to
reproduce the analysis of the paper. This I might adress in future versions, but as of now, I don´t recommend trying to knit this (if you do, you´ve been warned... ;))
