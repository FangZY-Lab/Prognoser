#' Prognoser: Integrative Prognostic Gene Signature Discovery
#'
#' This function integrates multiple expression datasets and corresponding
#' survival information to discover prognostic gene signatures. It performs
#' univariate Cox regression in each dataset, meta-analyzes the results,
#' screens genes based on effect direction and significance, and finally
#' clusters gene sets to remove redundancy.
#'
#' @param expression_accession_vector A character vector of dataset names.
#'        Each dataset (e.g., "GSE12345") must exist as a data frame in the
#'        global environment, with rows as genes and columns as samples.
#'        For each dataset, a corresponding survival data frame named
#'        \code{paste0(dataset, "_survival")} must also exist, containing
#'        columns \code{time} (numeric) and \code{status} (numeric, 0=alive,
#'        1=dead).
#' @param denovo_threshold Numeric. P-value threshold for univariate Cox
#'        regression to classify a gene as a "protection" (HR < 1) or "risk"
#'        (HR >= 1) gene within each dataset. Default is 0.05.
#' @param knowledge_genesets A list of pre-defined gene sets (character
#'        vectors) to be incorporated into the candidate gene set pool.
#' @param na_ratio Numeric. Maximum allowed proportion of missing values
#'        (NA) for a gene across datasets when building the effect-size
#'        matrix. Genes with NA proportion exceeding this value are removed.
#'        Default is 0.3.
#' @param gene_threshold Numeric. P-value threshold (from robust rank
#'        aggregation) for selecting genes into the final gene pools.
#'        Default is 0.001.
#' @param min_genes Integer. Minimum number of genes required in a final
#'        clustered gene set. Sets with fewer genes are discarded.
#'        Default is 5.
#' @param similarity_threshold Numeric. Jaccard similarity threshold for
#'        clustering gene sets. Gene sets with similarity >= this value are
#'        merged. Default is 0.5.
#' @param linkage Character. Linkage method for the hierarchical-clustering
#'        step. Must be one of \code{"ward.D"}, \code{"ward.D2"},
#'        \code{"single"}, \code{"complete"}, \code{"average"},
#'        \code{"mcquitty"}, \code{"median"} or \code{"centroid"}.
#'        Ignored when \code{de_redundant_basis = "igraph_component"}.
#'        Default is \code{"ward.D"}.
#' @param de_redundant_basis Character. How biologically redundant gene sets
#'        are consolidated. \code{"agglomeration"} (default) cuts the
#'        hierarchical clustering of the Jaccard distances at
#'        \code{1 - similarity_threshold} and unions each cluster.
#'        \code{"igraph_component"} builds an undirected graph whose edges
#'        are the gene-set pairs with Jaccard similarity at or above
#'        \code{similarity_threshold}, then unions each connected component.
#'        The two are related: with \code{linkage = "single"} the
#'        agglomeration partition is identical to the connected components of
#'        that graph, so \code{"igraph_component"} names the same
#'        consolidation explicitly, is faster on the sparse graph and returns
#'        the network itself.
#'
#' @return A list with three elements, plus two more when
#' \code{de_redundant_basis = "igraph_component"}:
#' \itemize{
#'   \item \code{genesets}: A list of final non-redundant gene signatures,
#'         named with \code{"LRGS"} (low-risk / protection) and \code{"HRGS"}
#'         (high-risk) prefixes.
#'   \item \code{combined_matrix_protection}: A matrix of combined effect
#'         scores (sign * -log10(P)) for protection-direction genes, including
#'         a column \code{p_protection} with robust rank aggregation scores.
#'   \item \code{combined_matrix_risk}: A matrix of combined effect scores
#'         for risk-direction genes, including a column \code{p_risk} with
#'   \item \code{network_protection}: The Jaccard edge list
#'         (\code{from}, \code{to}, \code{jaccard}) consolidated into the
#'         protection signatures. Only present with
#'         \code{de_redundant_basis = "igraph_component"}.
#'   \item \code{network_risk}: The same edge list for the risk signatures.
#'         robust rank aggregation scores.
#' }
#'
#' @importFrom igraph components graph_from_data_frame
#' @import survival
#' @import RobustRankAggreg
#' @import meta
#'
#' @examples
#' \dontrun{
#' # Assuming datasets "GSE10000" and "GSE20000" exist in the environment,
#' # along with their survival data frames "GSE10000_survival" and
#' # "GSE20000_survival".
#' my_known_genesets <- list(immune = c("CD8A", "GZMB", "PRF1"))
#' result <- Prognoser(
#'   expression_accession_vector = c("GSE10000", "GSE20000"),
#'   knowledge_genesets = my_known_genesets,
#'   denovo_threshold = 0.05,
#'   gene_threshold = 0.001,
#'   similarity_threshold = 0.6
#' )
#' final_signatures <- result$genesets
#' }
#'
#' @export


Prognoser=function(expression_accession_vector,
                   denovo_threshold=0.05,
                   knowledge_genesets,
                   na_ratio=0.3,
                   gene_threshold=0.001,
                   min_genes=5,
                   similarity_threshold=0.5,
                   linkage="ward.D",
                   de_redundant_basis=c("agglomeration", "igraph_component")){
  linkage <- match.arg(linkage, c("ward.D", "ward.D2", "single", "complete",
                                  "average", "mcquitty", "median", "centroid"))
  de_redundant_basis <- match.arg(de_redundant_basis)
  if (!length(expression_accession_vector)) {
    stop("`expression_accession_vector` must name at least one dataset.", call. = FALSE)
  }
  needed = c(expression_accession_vector,
             paste0(expression_accession_vector, "_survival"))
  missing_objects = needed[!vapply(needed, exists, logical(1))]
  if (length(missing_objects) > 0) {
    stop("Objects not found in the calling environment: ",
         paste(missing_objects, collapse = ", "),
         ". Each dataset needs `<name>` and `<name>_survival`.", call. = FALSE)
  }
  vector=expression_accession_vector
  list_data=list()
  list_survival_data=list()
  print("Start loading expression data and prognosis information.")
  for (q in seq_along(vector)) {
    list_data[[q]]=get(vector[q])
    list_survival_data[[q]]=get(paste0(vector[q],"_survival"))
    names(list_data)[q]=paste0(vector[q])
    names(list_survival_data)[q]=paste0(vector[q],"_survival")
  }
  modified_vector=lapply(vector, function(x) {
    ifelse(grepl("_",x),substr(x,1,regexpr("_", x)-1),x)
  })
  modified_vector=unique(as.character(modified_vector))
  large_list=c()
  for (item in seq_along(modified_vector)) {
    assign(paste0("list_", modified_vector[item]),list_data[grep(paste0("^",modified_vector[item]), names(list_data))])
    large_list=c(large_list,paste0("list_", modified_vector[item]))
  }
  print("Cox univariate analysis were started.")
  genesets_protection=list()
  genesets_risk=list()
  for (j in seq_along(large_list)) {
    expList=get(large_list[j])
    for(i in seq_along(expList)){
      print("Begin processing expression data and prognostic data.")
      data=expList[[i]]
      dimnames=list(rownames(data),colnames(data))
      data=matrix(as.numeric(as.matrix(data)),nrow=nrow(data),dimnames=dimnames)
      data=as.data.frame(data)
      gene_var=apply(data, 1, var)
      data=data[gene_var > 0, ]
      prognosis=list_survival_data[[paste0(names(expList)[i],"_survival")]]
      prognosis$time=as.numeric(prognosis$time)
      prognosis$status=as.numeric(prognosis$status)
      data=as.data.frame(t(data))
      prognosis_data=merge(prognosis,data,by="row.names")
      rownames(prognosis_data)=prognosis_data[,1]
      prognosis_data=prognosis_data[,-1]
      outTab=data.frame()
      for(sur in colnames(prognosis_data[,3:ncol(prognosis_data)])){
        cox = coxph(Surv(time,status) ~ prognosis_data[,sur], 
                    data = prognosis_data)
        coxSummary = summary(cox)
        outTab=rbind(outTab,
                     cbind(gene_id=sur,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           P.Value=coxSummary$coefficients[,"Pr(>|z|)"]))
      }
      outTab$HR=as.numeric(outTab$HR)
      outTab$P.Value=as.numeric(outTab$P.Value)
      outTab=na.omit(outTab)
      assign(paste0(names(expList)[i],"_result"),outTab)
    }
    if (length(expList)>1){
      for (mp in seq_along(expList)){
        if (mp == 1) { 
          result_single=get(paste0(names(expList)[mp],"_result"))
          result_p=result_single[,c("gene_id","P.Value")]
          colnames(result_p)[2]=c("P.Value1")
          result_HR=result_single[,c("gene_id","HR")]
          colnames(result_HR)[2]=c("HR1")
          result_p_all=result_p
          result_HR_all=result_HR
        } else {
          result_single=get(paste0(names(expList)[mp],"_result"))
          result_p=result_single[,c("gene_id","P.Value")]
          colnames(result_p)[2]=paste0("P.Value",mp)
          result_HR=result_single[,c("gene_id","HR")]
          colnames(result_HR)[2]=paste0("HR",mp)
          result_p_all=merge(result_p_all,result_p,by="gene_id")
          result_HR_all=merge(result_HR_all,result_HR,by="gene_id")
        }
      }
      rownames(result_p_all)=result_p_all[,1]
      result_p_all=result_p_all[,-1]
      rownames(result_HR_all)=result_HR_all[,1]
      result_HR_all=result_HR_all[,-1]
      result_p_all = as.data.frame(
        lapply(result_p_all, function(x) {
          non_zero = x[x != 0]
          if(length(non_zero) > 0) {
            min_val = min(non_zero)
            x[x == 0] = min_val
          }
          return(x)
        }),
        row.names = rownames(result_p_all)
      )
      result_HR_all = as.data.frame(apply(result_HR_all, 2, function(col) {
        finite_vals = col[is.finite(col) & col != 0]
        if (length(finite_vals) > 0) {
          col_min = min(finite_vals, na.rm = TRUE)
          col_max = max(finite_vals, na.rm = TRUE)
        } else {
          col_min = NA
          col_max = NA
        }
        col[col == 0] = col_min
        col[is.infinite(col) & col > 0] = col_max
        col[is.infinite(col) & col < 0] = col_min
        return(col)
      }))
      common_genes = intersect(rownames(result_HR_all), rownames(result_p_all))
      result_HR_all = result_HR_all[common_genes, ]
      result_p_all = result_p_all[common_genes, ]
      meta_results = data.frame(
        gene_id = common_genes,
        combined_HR = NA,
        combined_HR_lower = NA,
        combined_HR_upper = NA,
        combined_P = NA,
        Q = NA,
        I2 = NA,
        p_heterogeneity = NA,
        n_studies = NA,
        stringsAsFactors = FALSE
      )
      for(mt in seq_along(common_genes)) {
        gene = common_genes[mt]
        hr_values = as.numeric(result_HR_all[gene, ])
        p_values = as.numeric(result_p_all[gene, ])
        valid_idx = !is.na(hr_values) & !is.na(p_values) & is.finite(hr_values)
        hr_values = hr_values[valid_idx]
        p_values = p_values[valid_idx]
        if(length(hr_values) < 1) next
        z_scores = sign(log(hr_values)) * qnorm(1 - p_values/2)
        log_hr = log(hr_values)
        se_loghr = log_hr / z_scores
        meta_data = data.frame(
          study = colnames(result_HR_all)[valid_idx],
          TE = log_hr,
          seTE = se_loghr,
          stringsAsFactors = FALSE
        )
        tryCatch({
          meta_analysis = metagen(TE = TE,
                                   seTE = seTE,
                                   data = meta_data,
                                   studlab = study,
                                   comb.fixed = FALSE,
                                   comb.random = TRUE,
                                   method.tau = "DL",
                                   hakn = FALSE,
                                   prediction = FALSE,
                                   sm = "HR")
          meta_results[mt, "combined_HR"] = exp(meta_analysis$TE.random)
          meta_results[mt, "combined_HR_lower"] = exp(meta_analysis$lower.random)
          meta_results[mt, "combined_HR_upper"] = exp(meta_analysis$upper.random)
          meta_results[mt, "combined_P"] = meta_analysis$pval.random
          meta_results[mt, "Q"] = meta_analysis$Q
          meta_results[mt, "I2"] = meta_analysis$I2
          meta_results[mt, "p_heterogeneity"] = meta_analysis$pval.Q
          meta_results[mt, "n_studies"] = length(hr_values)
          
        }, error = function(e) {
          cat(paste("Gene", gene, "meta-analysis error:", e$message, "\n"))
        })
      }
      result=data.frame(gene_id=meta_results$gene_id,HR=meta_results$combined_HR,P.Value=meta_results$combined_P)
      assign(paste0(gsub("_.*$","",names(expList)[1]),"_result"),result)
    }else{
      print("Genesets are not processed for merging.")
      result_single=get(paste0(names(expList)[1],"_result"))
      result=data.frame(gene_id=result_single$gene_id,HR=result_single$HR,P.Value=result_single$P.Value)
      assign(paste0(gsub("_.*$","",names(expList)[1]),"_result"),result)
    }
  }
  print("Start screening genes.")
  for (mq in seq_along(modified_vector)) {
    result=get(paste0(modified_vector[mq],"_result"))
    result=result[!is.na(result$HR) & !is.na(result$P.Value), , drop=FALSE]
    result_protection=result[result$HR < 1,]
    result_risk=result[result$HR >= 1,]
    result_protection=result_protection[order(result_protection$P.Value), ]
    result_risk=result_risk[order(result_risk$P.Value), ]
    result_protection=result_protection[result_protection$P.Value < denovo_threshold, ]
    result_risk=result_risk[result_risk$P.Value < denovo_threshold, ]
    result_protection=result_protection$gene_id
    assign(paste0(modified_vector[mq],"_protection"),result_protection)
    result_risk=result_risk$gene_id
    assign(paste0(modified_vector[mq],"_risk"),result_risk)
  }
  combined_protection=list()
  combined_risk=list()
  for (mn in seq_along(modified_vector)) {
    combined_protection[[mn]]=get(paste0(modified_vector[mn],"_protection"))
    names(combined_protection)[mn]=paste0(modified_vector[mn], "_protection")
    combined_risk[[mn]]=get(paste0(modified_vector[mn],"_risk"))
    names(combined_risk)[mn]=paste0(modified_vector[mn], "_risk")
  }
  padded_list_protection=c(combined_protection,knowledge_genesets)
  padded_list_risk=c(combined_risk,knowledge_genesets)
  names(padded_list_protection)=sprintf("LRGS%d", seq_along(padded_list_protection))
  names(padded_list_risk)=sprintf("HRGS%d", seq_along(padded_list_risk))
  for (i in seq_along(padded_list_protection)) {
    padded_list_protection[[i]]=unique(padded_list_protection[[i]])
  }
  print("Check for duplicates within each low risk gene set.")
  sapply(padded_list_protection, function(x) any(duplicated(x)))
  for (i in seq_along(padded_list_risk)) {
    padded_list_risk[[i]]=unique(padded_list_risk[[i]])
  }
  print("Check for duplicates within each high risk gene set.")
  sapply(padded_list_risk, function(x) any(duplicated(x)))
  print("Start building the effect size matrix.")
  for (mn in seq_along(modified_vector)) {
    if (mn == 1){
      pmatrix=get(paste0(modified_vector[mn],"_result"))[,c("gene_id","P.Value")]
      colnames(pmatrix)[mn+1]=paste0(modified_vector[mn])
      hrmatrix=get(paste0(modified_vector[mn],"_result"))[,c("gene_id","HR")]
      colnames(hrmatrix)[mn+1]=paste0(modified_vector[mn])
    }else{
      pmatrix=merge(pmatrix,get(paste0(modified_vector[mn],"_result"))[,c("gene_id","P.Value")],by="gene_id",all=T)
      colnames(pmatrix)[mn+1]=paste0(modified_vector[mn])
      hrmatrix=merge(hrmatrix,get(paste0(modified_vector[mn],"_result"))[,c("gene_id","HR")],by="gene_id",all=T)
      colnames(hrmatrix)[mn+1]=paste0(modified_vector[mn])
    }
  }
  rownames(pmatrix)=pmatrix[,1]
  pmatrix=pmatrix[,-1,drop=FALSE]
  pmatrix$na_count=rowSums(is.na(pmatrix))/ncol(pmatrix)
  pmatrix=pmatrix[pmatrix$na_count <= na_ratio,]
  pmatrix=pmatrix[,!colnames(pmatrix) %in% "na_count",drop=FALSE]
  rownames(hrmatrix)=hrmatrix[,1]
  hrmatrix=hrmatrix[,-1,drop=FALSE]
  hrmatrix$na_count=rowSums(is.na(hrmatrix))/ncol(hrmatrix)
  hrmatrix=hrmatrix[hrmatrix$na_count <= na_ratio,]
  hrmatrix=hrmatrix[,!colnames(hrmatrix) %in% "na_count",drop=FALSE]
  hrmatrix=ifelse(hrmatrix>1,1,-1)
  pmatrix=-log10(pmatrix)
  combined_matrix=hrmatrix*pmatrix
  combined_matrix_protection=combined_matrix
  combined_matrix_risk=combined_matrix
  rank_list_protection = lapply(combined_matrix_protection, function(x) {
    if(is.numeric(x)) {
      ranks = rank(x, na.last = "keep", ties.method = "min")
      rownames(combined_matrix_protection)[order(ranks, na.last = NA)]
    }
  })
  aggregated_ranks = aggregateRanks(rank_list_protection)
  combined_matrix_protection$p_protection = aggregated_ranks$Score[match(rownames(combined_matrix_protection), aggregated_ranks$Name)]
  rank_list_risk = lapply(combined_matrix_risk, function(x) {
    if(is.numeric(x)) {
      ranks = rank(-x, na.last = "keep", ties.method = "min")
      rownames(combined_matrix_risk)[order(ranks, na.last = NA)]
    }
  })
  aggregated_ranks = aggregateRanks(rank_list_risk)
  combined_matrix_risk$p_risk = aggregated_ranks$Score[match(rownames(combined_matrix_risk), aggregated_ranks$Name)]
  protection_gene_pool=combined_matrix_protection[combined_matrix_protection$p_protection < gene_threshold,]
  risk_gene_pool=combined_matrix_risk[combined_matrix_risk$p_risk < gene_threshold,]
  update_protection_geneset=list()
  protection_geneset=padded_list_protection
  for (mi in seq_along(names(protection_geneset))){
    coldn=protection_geneset[[mi]]
    coldn=coldn[coldn != ""]
    col_intersect_dn=intersect(coldn,rownames(protection_gene_pool))
    update_protection_geneset[[mi]]=col_intersect_dn
    names(update_protection_geneset)[mi]=paste0(names(protection_geneset)[mi],"_purification")
  }
  update_risk_geneset=list()
  risk_geneset=padded_list_risk
  for (mi in seq_along(names(risk_geneset))){
    coldn=risk_geneset[[mi]]
    coldn=coldn[coldn != ""]
    col_intersect_dn=intersect(coldn,rownames(risk_gene_pool))
    update_risk_geneset[[mi]]=col_intersect_dn
    names(update_risk_geneset)[mi]=paste0(names(risk_geneset)[mi],"_purification")
  }
  GENE_SETS_protection=update_protection_geneset[
    sapply(update_protection_geneset, function(x) {
      sum(x == "" | is.na(x) | x == " ") < length(x)
    })
  ]
  update_GENE_SETS_protection=list()
  for (i in seq_along(names(GENE_SETS_protection))){
    coldn=GENE_SETS_protection[[i]]
    coldn=coldn[coldn != ""]
    update_GENE_SETS_protection[[i]]=coldn
    names(update_GENE_SETS_protection)[i]=paste0(names(GENE_SETS_protection)[i])
  }
  
  GENE_SETS_risk=update_risk_geneset[
    sapply(update_risk_geneset, function(x) {
      sum(x == "" | is.na(x) | x == " ") < length(x)
    })
  ]
  update_GENE_SETS_risk=list()
  for (i in seq_along(names(GENE_SETS_risk))){
    colup=GENE_SETS_risk[[i]]
    colup=colup[colup != ""]
    update_GENE_SETS_risk[[i]]=colup
    names(update_GENE_SETS_risk)[i]=paste0(names(GENE_SETS_risk)[i])
  }
  integration_ratio=function(set1, set2) {
    intersection=length(intersect(set1, set2))
    union=length(union(set1, set2))
    ratio=intersection/union
    return(ratio)
  }
  print("Remove redundancy by clustering.")
  aux.merge = function(gsets, minjac, linkage_method, basis){
    if (length(gsets) < 2L) {
      empty_net = data.frame(from = character(0), to = character(0),
                             jaccard = numeric(0), stringsAsFactors = FALSE)
      return(list(sets = lapply(gsets, function(x) sort(unique(x))),
                  network = if (basis == "igraph_component") empty_net else NULL))
    }
    integration_matrix = outer(seq_along(gsets), seq_along(gsets),
                               Vectorize(function(x, y) integration_ratio(gsets[[x]],gsets[[y]])))
    dimnames(integration_matrix) = list(names(gsets), names(gsets))
    integration_matrix[is.na(integration_matrix)] = 0
    integration_matrix[is.infinite(integration_matrix)] = 0
    integration_matrix = pmax(pmin(integration_matrix, 1), 0)
    if (basis == "igraph_component") {
      upper = upper.tri(integration_matrix)
      keep = upper & integration_matrix >= minjac
      network = data.frame(from = rownames(integration_matrix)[row(integration_matrix)[keep]],
                           to = colnames(integration_matrix)[col(integration_matrix)[keep]],
                           jaccard = integration_matrix[keep],
                           stringsAsFactors = FALSE)
      if (nrow(network) == 0) {
        print("No gene-set pair reaches the similarity threshold; nothing to merge.")
        merged = lapply(gsets, function(x) sort(unique(x)))
      } else {
        g = igraph::graph_from_data_frame(network[, c("from", "to", "jaccard")],
                                          directed = FALSE)
        comp = igraph::components(g)
        groups = split(names(comp$membership), comp$membership)
        merged = lapply(groups, function(nm) sort(unique(unlist(gsets[nm]))))
        isolated = gsets[setdiff(names(gsets), names(comp$membership))]
        if (length(isolated) > 0) {
          merged = c(lapply(isolated, function(x) sort(unique(x))), merged)
        }
      }
      return(list(sets = merged, network = network))
    }
    dist_matrix = as.dist(1 - integration_matrix)
    hc = hclust(dist_matrix, method = linkage_method)
    ## centroid / median linkage can produce inversion ("crossing")
    ## heights. cutree() refuses an unsorted tree, but reordering
    ## hc$merge would break the parent references, so make the heights
    ## monotone and leave the merge structure untouched.
    if (is.unsorted(hc$height)) hc$height = cummax(hc$height)
    clust = cutree(hc, h = 1 - minjac)
    list(sets = as.list(tapply(seq_along(gsets), clust, function(i){sort(unique(unlist(gsets[i])))})),
         network = NULL)
  }
  merged_protection = aux.merge(gsets = update_GENE_SETS_protection, minjac = similarity_threshold, linkage_method = linkage, basis = de_redundant_basis)
  merged_risk = aux.merge(gsets = update_GENE_SETS_risk, minjac = similarity_threshold, linkage_method = linkage, basis = de_redundant_basis)
  update_GENE_SETS_disjunction_protection = merged_protection$sets
  update_GENE_SETS_disjunction_protection=update_GENE_SETS_disjunction_protection[!duplicated(lapply(update_GENE_SETS_disjunction_protection, sort))]
  update_GENE_SETS_disjunction_risk = merged_risk$sets
  update_GENE_SETS_disjunction_risk=update_GENE_SETS_disjunction_risk[!duplicated(lapply(update_GENE_SETS_disjunction_risk, sort))]
  update_GENE_SETS_disjunction_protection=update_GENE_SETS_disjunction_protection[sapply(update_GENE_SETS_disjunction_protection,length) >= min_genes]
  update_GENE_SETS_disjunction_risk=update_GENE_SETS_disjunction_risk[sapply(update_GENE_SETS_disjunction_risk,length) >= min_genes]
  names(update_GENE_SETS_disjunction_protection)=sprintf("LRGS%d", seq_along(update_GENE_SETS_disjunction_protection))
  names(update_GENE_SETS_disjunction_risk)=sprintf("HRGS%d", seq_along(update_GENE_SETS_disjunction_risk))
  update_GENE_SETS_disjunction=c(update_GENE_SETS_disjunction_risk,update_GENE_SETS_disjunction_protection)
  output=list(update_GENE_SETS_disjunction,combined_matrix_protection,combined_matrix_risk)
  names(output)=c("genesets","combined_matrix_protection","combined_matrix_risk")
  if (de_redundant_basis == "igraph_component") {
    output$network_protection = merged_protection$network
    output$network_risk = merged_risk$network
  }
  return(output)
}














