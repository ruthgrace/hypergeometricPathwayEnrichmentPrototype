hypergeometricPathwayEnrichmentPrototype
========================================

This is a NetBeans project. There is no GUI at this time. It simply processes and outputs files.

The hypergeometricWithWikiPathways method generates genesets from WikiPathways .gpml files (the path for which must be provided in the WIKIPATHWAYSFOLDER variable), and performs the hypergeometric test on the annotated VCF file (given in the VCFFILE variable). The output is a table of pathway names and p-values in output/enriched_pathways_and_Pvalues.txt, a list of genes found in the VCF file but not in the genesets (output/genes_not_in_genesets.txt), and a folder (output/gpmlfiles/) that contains all the .gpml pathways with the rectangle outlining mutated genes bolded.

The .gpml files can be opened in Cytoscape with a plugin (http://apps.cytoscape.org/apps/gpmlplugin) and unfortunately only works with Cytoscape 2.7 and 2.8.

Here is the Cytoscape visualization of a .gpml file before running the hypergeometric test.

![alt tag](https://github.com/ruthgrace/hypergeometricPathwayEnrichmentPrototype/blob/master/readmefiles/before.png)

Here is the Cytoscape visualization of the same .gpml file output after running the hypergeometric test.

![alt tag](https://github.com/ruthgrace/hypergeometricPathwayEnrichmentPrototype/blob/master/readmefiles/after.png)

The program can be run using GeneSets without WikiPathways by specifying the GENESETFILE variable and by using the hypergeometricWithGeneSets method instead of the hypergeometricWithWikiPathways method.


Dependencies include jsc.jar for the hypergeometric test.
