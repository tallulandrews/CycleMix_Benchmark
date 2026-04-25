source("~/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Stuff.R")

convert_geneset <- function(gene_tab, in.org, out.org, in.name, out.name) {
	new_names <- General_Map(gene_tab[,1], in.org=in.org, out.org=out.org, in.name=in.name, out.name=out.name)
	remove <- new_names == ""
	gene_tab <- gene_tab[!remove,]
	gene_tab[,1] <- new_names[!remove]
	gene_tab <- unique(gene_tab)
	return(gene_tab)
}
