#### 
# LEANR+ORA Function:
####
# Author: Ashley Morrison
# Date: 20230821
# Purpose: Demonstration of the LEANR_ORA pipeline and the function for ORA from the LEANR output. 
####

## Workflow:
# Differential Expression analysis --> LEANR --> LEANR_ORA
# Subnetwork analysis (LEANR) and ORA from subnetworks (LEANR_ORA) are run in two seaprate steps for ease of use. 
# An example workflow and the LEANR_ORA function are below. 

##################################################################
###  LEANR:
# The input for LEANR is the list of adjusted p-values from differential expression (list values) and protein IDs (list names). 
# IMPORTANT!! ALL protein IDs MUST be consistent between the PPI and p-value list for LEANR, as well as the proteins in the signature lists. 
# It is easiest to just convert all names to STRING_IDs as the PPI is made using the data from the STRING_db
# A full list of STRING aliases can be downloaded from the STRING db, or are in "../../leanr_cluster_proteomics/source_data/9606.protein.aliases.v11.5.txt")
# When converting names, most proteins are present in the ____ and ____ data sources, but the conversion needs to happen using one source at the time as there are repeated proteins

### Example PPI creation code with generation of proteins_in_graph list (used for LEANR_ORA analysis):
#library(igraph)
#ppi <- read.table("../../STRING_db/9606.protein.links.v11.5.txt",header=T) <- full list of STRING interaction data 
#ppi_exp <- ppi 

## filter to JUST proteins from experimental data:
#prot <- limma_output$STRING_id <- the STRING_ID for proteins/genes in limma analysis
#ppi_filt <- ppi_exp %>% filter(protein1 %in% prot & protein2 %in% prot) <- remove any interactions that aren't between experimental proteins

## filter out edges below a user specified condifence threshold (in this case, 600 (the same as 0.6))
#ppi_filt_600 <- ppi_filt %>% filter(combined_score >= 600)

#nrow(ppi_filt) <- 313538 edges after filtering to just experimental proteins
#nrow(ppi_filt_600) <- 10600 edges after filtering to 600 CS

## make graph using igraph R package:
#ppi_graph_600 <- graph_from_data_frame(ppi_filt_600, directed=FALSE, vertices=NULL) <- make graph
#E(ppi_graph_600)$weight <- ppi_filt_600$combined_score <- LEANR required a weighted network, but the weights are NOT used in calculations

## Make list of proteins for later use in ORA
#proteins_in_graph_600 <- V(ppi_graph_600)$name <- pull network protein names for use in LEANR_ORA function

### Example p-value list creation:
#LEANR_pval <- (as.vector(limma_results$adj.P.Val)) <- adjusted p-values from LIMMA
#names(LEANR_pval) <- (limma_results$STRING_id) <- STRING_IDs for the genes/proteins in limma analysis

### Running LEANR:
#library(LEANR)
#LEAN_results <- run.lean(LEANR_pval,ppi_graph_600, n_reps = 10000)


##################################################################

### LEANR_ORA:

## The input for the LEANR_ORA function is:
# AD_signatures: a signature dataframe with at least 2 columns: 'STRING_id' (the protein's STRING ID) and signature_key (the signatures that the protein/gene is assigned to)
# leanr_res: the entire output of the LEANR function
# proteins_in_graph: The list of proteins in the experimental PPI (the final list of proteins used in the analysis, NOT the full human PPI)
# overlap_cutoff: the minimum overlap required between a subnetwork and signature for it to be considered significantly enriched, if not set then will be 10
# min_sig_size: the minimum signature size for analysis, if not set then will be 10
# max_sig_size: the maximum signature size for analysis, if not set then will be 5,000

## LEANR_ORA Function:
LEANR_ORA <- function(AD_signatures,leanr_res,proteins_in_graph,overlap_cutoff,min_sig_size,max_sig_size){
  
  ## Set the optional argument values:
  if(missing(overlap_cutoff)) {
    overlap_cutoff <- 10
  } 
  
  if(missing(min_sig_size)) {
    min_sig_size <- 10
  } 
  
  if(missing(max_sig_size)) {
    max_sig_size <- 5000
  } 
  
  ## remove non-coding genes and the proteins that didn't match the STRING database
  AD_sig_no_na <- AD_sig[!is.na(AD_sig$STRING_id),]
  GS_list <- split(x = AD_sig_no_na$STRING_id, f = AD_sig_no_na$signature_key) #final gene set list
  
  ## Pull significant subnetworks from the LEANR output:
  leanr_res_full <- as.data.frame(leanr_res$restab) #pull results tab as data frame
  leanr_res_full_signf <- leanr_res_full %>% filter(PLEAN <= .05) #pull central proteins to pull subnetworks for network generation
  
  #end the function if no significant subnetworks by PLEAN:
  if (nrow(leanr_res_full_signf)==0) {
    return(print("no significant subnetworks"))
    break
  }
  
  leanr_res_signf_list <- rownames(leanr_res_full_signf) # make list of just the significant central proteins to filter subnetworks by
  
  leanr_res_nhs <- leanr_res$nhs #make list of all subnetworks from the LEANR output
  leanr_res_signf_subnetworks <- leanr_res_nhs[names(leanr_res_nhs) %in% leanr_res_signf_list] # pull JUST the subnetworks that were significantly enriched
  
  ## Run ORA on each subnetwork individually:
  fgsea_ora_res <- lapply(seq_along(leanr_res_signf_subnetworks), function(i){
    fora(pathways=GS_list, #list of AD gene sets
         genes = leanr_res_signf_subnetworks[[i]], #run through proteins in each significant subnetwork
         universe = proteins_in_graph, #all proteins in graph to compare to as background
         minSize = min_sig_size, # cutoffs for signature size
         maxSize = max_sig_size
    ) %>% 
      mutate(subnetwork = names(leanr_res_signf_subnetworks)[i]) #add in subnetwork column 
  }) %>% 
    data.table::rbindlist() %>% #combine all previously generated tables
    filter(padj < 0.05) %>% #filter output to only include significantly enriched signatures by ORA
    filter(overlap>= overlap_cutoff) %>% #cut off for number of overlapped proteins between signature and subnetwork
    arrange(subnetwork, padj) #arrange by adjusted pvalue
  
  ## add subnetwork size for each enriched subnetwork:
  fgsea_ora_res$network_size <- "NA"
  for (i in 1:nrow(fgsea_ora_res)) {
    fgsea_ora_res$network_size[i] <- length(leanr_res_signf_subnetworks[[fgsea_ora_res$subnetwork[i]]])
  }
  
  ## add column headers to differentiate the different statistical values:
  colnames(fgsea_ora_res) <- c("signature","ORA_pval","ORA_padj","subnetwork_overlap_w_signature","signature_size","subnetwork_proteins_in_signature","subnetwork_center_protein","subnetwork_size")
  
  ## add in original subnetwork PLEAN values from LEANR analysis:
  leanr_res_full_signf_plean <- leanr_res_full_signf %>% dplyr::select(PLEAN)
  leanr_res_full_signf_plean$protein <- rownames(leanr_res_full_signf_plean)
  
  fgsea_ora_res <- fgsea_ora_res %>% left_join(leanr_res_full_signf_plean, by=c("subnetwork_center_protein"="protein"))
  colnames(fgsea_ora_res)[colnames(fgsea_ora_res)=="PLEAN"] <- "subnetwork_PLEAN_score"
  
  ## add -ln(ORA_padj) for plotting visualization:
  fgsea_ora_res$ORA_ln_padj <- -log(fgsea_ora_res$ORA_padj)
  
  ## re-order columns of results table for better visualization:
  fgsea_ora_res <- fgsea_ora_res[, c("signature", "ORA_pval", "ORA_padj", "ORA_ln_padj","subnetwork_PLEAN_score","subnetwork_overlap_w_signature","subnetwork_size","signature_size",  "subnetwork_center_protein", "subnetwork_proteins_in_signature")]
  
  ## plot the most significant (aka lowest padj) adjust p-value for each significant signature
  fgsea_ora_res_max_score <- aggregate(fgsea_ora_res$ORA_ln_padj, by = list(fgsea_ora_res$signature), max)
  colnames(fgsea_ora_res_max_score) <- c("signature","-ln_padj")
  
  p <- ggplot(fgsea_ora_res_max_score, aes(x=`-ln_padj`,y=signature, fill=`-ln_padj`))+
    geom_bar(stat="identity", width = 0.40)+
    labs(fill = "-ln(padj)")+
    ggtitle("Enriched Signatures")+
    theme(axis.title.y=element_blank())+
    xlab("-ln(padj)")+ 
    #guides(fill=guide_legend(title="-ln(padj)"))+
    scale_fill_gradient(low = "#56B4E9", high = "#0072B2")
  
  ## final result list:
  res <- list(fgsea_ora_res,p)
  names(res) <- c("LEANR_ORA_results","enrichment_plot")
  
  
  return(res)
}

## example LEANR_ORA
#ORA_600 <- LEANR_ORA(LEAN_results,proteins_in_graph_600,overlap_cutoff=10)

##################################################################
### Function for running JUST ORA on limma output:

DEP_ORA <- function(AD_signatures,limma_res,overlap_cutoff,min_sig_size,max_sig_size) {
  # limma_res_tab is full limma output that will be filtered in the analysis
  #temp var to test function:
  #limma_res <- C1_all
  #overlap_cutoff <- 10
  
  ## Set the optional argument values:
  if(missing(overlap_cutoff)) {
    overlap_cutoff <- 10
  } 
  
  if(missing(min_sig_size)) {
    min_sig_size <- 10
  } 
  
  if(missing(max_sig_size)) {
    max_sig_size <- 5000
  } 
  
  ## remove non-coding genes and the proteins that didn't match the STRING database
  AD_sig_no_na <- AD_sig[!is.na(AD_sig$STRING_id),]
  GS_list <- split(x = AD_sig_no_na$STRING_id, f = AD_sig_no_na$signature_key) #final gene set list
  
  
  limma_res_signif <- limma_res %>% filter(adj.P.Val <= .05 & !is.na(STRING_id))
  
  #pull list of DE proteins:
  DE_list <- limma_res_signif$STRING_id
  
  #pull full list of proteins for universal set:
  universe_list <- limma_res[!is.na(limma_res$STRING_id),"STRING_id"]
  
  #Run ORA
  just_ora <-fora(pathways=GS_list, #list of AD gene sets
                  genes=DE_list, #list of proteins to do ORA on
                  universe = universe_list, #all proteins in graph to compare to as background
                  minSize = 10, #set to 10 bc this is the minimum overlap
                  maxSize = 5000 #set to be larger than the biggest gene set
  ) %>% 
    filter(overlap >= overlap_cutoff)  #cut off for number of overlapped proteins
  
  return(just_ora)
}



