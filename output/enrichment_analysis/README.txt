This file describes the contents of the output/enrichment_analysis folder.

expandedInteractomes/

	cytoscapeFormat_[bait].txt: tab-delimited file containing the corum complex data for the corum interactions related to the bait. Can be combined with original interactome bait-prey file for the full corum-expanded interactome.
	
	[bait]_complexMembers.txt: text file containing gene names of corum-expanded interactions

	[bait]_complexNames.txt: text file containing the names of corum complexes that were included to expand the interactomes.


#####################
commonComplexes/

	commonComplexes_list.txt : this file contains a list of unique complex identifiers that are shared between multiple interactomes. 

	commonComplexes_cytoscapeFormat.txt - a tab-delimited file containing the protein interaction data from the complexes listed in commonComplexes_list. Contains the following columns:
		Interactor1	
		Interactor2	
		Complex.name	
		Complex.id	
		geneName1	
		geneName2	
		GATA4	
		TBX5	
		NKX25

	original_interactome_overlaps.txt: from the originally-provided interactomes, which proteins are shared across networks? A "1" in a bait column indicates that the protein is present in this bait's interactome.
	
	expanded_interactome_overlaps.txt: from the forum-expanded interactomes, which proteins are shared? A "1" in the bait column indicates that the protein is present in this bait's corum-expanded interactome. 

#####################
enriched4interactors/

	[bait]_complexes.txt: a list of corum complexes that showed significant overrepresentation of original interactome members (by Fisher's p-value, with a background set of all possible corum interactors.)

#####################
enriched4PCGC/

	interactome_complex_enrichment.csv: comma-delimited file of all corum complexes that contain at least one protein in the three bait interactomes, the Fisher's p-value for number of DNV/LoF/rec/combined mutations in that corum complex, and which bait interactome(s) potentially interact with this corum complex. Columns include the following:
		ComplexID	
		ComplexName	
		DNV_pval	
		LoF_pval	
		rec_pval	
		combined_pval	
		GATA4	(0 indicates bait is not known to interact with this complex)
		TBX5	
		NKX25	
		Genes	(List of genes that make up the complex)
