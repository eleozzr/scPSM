#' @include utils.R
#'
NULL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Main Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' scPSM integration
#'
#' Correct for batch effects, impute dropouts, and denoise data
#' using the propensity score matching method for single-cell RNA-sequencing data.
#'
#' @param batches One or more log-expression matrices where genes correspond to rows and cells correspond to columns.
#' Each matrix should contain the same number of rows, corresponding to the same genes in the same order.
#' Each matrix represents a batch.
#' @param markers A vector specifying which features used as marker genes to compute propensity scores.
#' @param hvg A vector specifying which features used as HVGs for identifying MNN group.
#' @param k.self An integer scalar specifying the number of nearest neighbors in searching KNNs.
#' @param k.mnn An integer scalar specifying the number of nearest neighbors in matching MNN pairs.
#' @param correct.all A logical scalar specifying whether correction should be applied to all genes, even if only a subset is used for the MNN group identification.
#' @param merge.order An integer vector containing the linear merge order of batches.
#'
#' @importFrom Matrix t
#' @importFrom stats predict glm formula binomial
#' @importFrom BiocNeighbors queryKNN KmknnParam
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' #Please see the tutorial for the preprocessing and then run the following code
#' #psm.data <- psm_integrate(batches = batches, markers = markers, hvg = hvg, merge.order = 1:4)


psm_integrate <- function(batches, markers, hvg, k.self=10, k.mnn=10,
                          correct.all=TRUE, merge.order=1:4)
{
  prep.out <- .prepare_input_data(batches=batches, cos.norm.in=T, cos.norm.out=T,
                                  subset.row=hvg, correct.all=T)
  in.batches <- prep.out$In # matrices of hvgs
  out.batches <- prep.out$Out # matrices of all genes
  subset.row <- prep.out$Subset # position of hvgs
  same.set <- prep.out$Same

  in.batches <- lapply(in.batches, t)
  if (!same.set) {
    out.batches <- lapply(out.batches, t)
  }

  restrict <- NULL
  merge.tree <- .create_tree_predefined(in.batches, restrict, merge.order)
  NEXT <- .get_next_merge
  UPDATE <- .update_tree

  # Putting the out.batches information in the 'extras'.
  merge.tree <- .add_out_batches_to_tree(merge.tree, if (same.set) NULL else out.batches)

  nbatches <- length(unlist(merge.tree))
  nmerges <- nbatches - 1L
  mnn.pairings <- left.set <- right.set <- vector("list", nmerges)
  for (mdx in seq_len(nmerges)) {
    # Traversing the merge tree to find the next two batches to merge.
    next.step <- NEXT(merge.tree) # list of nmerges
    left <- next.step$left # Large MNN_treenode
    right <- next.step$right # Large MNN_treenode

    left.data <- .get_node_data(left) # dgCMatrix of cell * hvgs
    left.restrict <- .get_node_restrict(left) # NULL
    left.index <- .get_node_index(left) # 1L
    left.origin <- .get_node_origin(left) # seq of 1
    left.extras <- .get_node_extras(left)[[1]] # dgCMatrix of cell * genes

    right.data <- .get_node_data(right)
    right.restrict <- .get_node_restrict(right)
    right.index <- .get_node_index(right)
    right.origin <- .get_node_origin(right)
    right.extras <- .get_node_extras(right)[[1]]

    # Really no point being too cute here with other matrix representations,
    # there are coercions for the NN search _and_ for the dense per-gene output.
    left.data <- as.matrix(left.data)
    right.data <- as.matrix(right.data)

    # Find the k-nearest neighbors in data for each point in query
    W_right <- queryKNN(right.data, query=right.data, k=k.self,
                        BNPARAM=KmknnParam(), BPPARAM=SerialParam(), get.distance=FALSE)

    # Find mnn pairs using hvgs
    mnn.sets <- .restricted_mnn(left.data, left.restrict, right.data, right.restrict,
                                k=k.mnn, prop.k=NULL, BNPARAM=KmknnParam(), BPPARAM=SerialParam())

    s1 <- mnn.sets$first
    s2 <- mnn.sets$second

    mnn.pairings[[mdx]] <- DataFrame(left=s1, right=s2) # DF of mnn pairs
    left.set[[mdx]] <- left.index
    right.set[[mdx]] <- right.index

    # calculate corrected right. data and extras
    # to project right to left, for each cell in right, find mnn pairs in left to cells of knn pairs in right
    n_left <- nrow(left.data)
    n_right <- nrow(right.data)
    parings <- lapply(1:n_right, function(x){
      subset(as.data.frame(mnn.pairings[[mdx]]), right %in% W_right[["index"]][x,])})

    # calculate propensity scores and weights
    # left - comparison, right - treatment, project treat to comparison
    batch.id <- rep(c(0,1), c(n_left, n_right))
    if (correct.all) {
      raw <- rbind(left.extras, right.extras)
    } else {
      raw <- rbind(left.data, right.data)
    }
    psm.data <- as.matrix(raw)[, markers]
    psm.data <- cbind(as.vector(batch.id), psm.data)
    colnames(psm.data)[1] <- "treat"

    # for all cells in right batch
    psm.datas <- lapply(1:n_right, function(x){
      tmp0=as.data.frame(psm.data[c(parings[[x]]$left, parings[[x]]$right+n_left),])
      #tmp0$treat=factor(tmp0$treat)
      return(tmp0)
    })

    psm.logit <- lapply(1:n_right, function(x){
      if (length(parings[[x]]$left) > 1) {
        glm(formula(psm.datas[[x]]),
            data = psm.datas[[x]],
            family = binomial)
      } else NULL
    })

    w.logit <- lapply(1:n_right, function(x){
      if (length(parings[[x]]$left) > 1) {
        #tmp0=psm.datas[[x]]
        exp(predict(psm.logit[[x]],subset(psm.datas[[x]],treat==0)))
      } else NULL
    })
    w.logit <- lapply(w.logit, function(x) replace(x, !is.finite(x), 1))
    #save(psm.logit,psm.datas,n_right,w.logit,parings,left.extras,left.data,right.data,file="~/Desktop/test_00.RData")
    # Applying the correction and expanding the reference batch.
    right.extras.cor <- sapply(1:n_right, function(x){
      if (length(parings[[x]]$left) > 1) {
        Matrix::colSums(w.logit[[x]]*left.extras[parings[[x]]$left, ])/sum(w.logit[[x]])
      } else {
        right.extras[x,]
      }
    })
    right.data.cor <- sapply(1:n_right, function(x){
      if (length(parings[[x]]$left) > 1) {
        Matrix::colSums(w.logit[[x]]*left.data[parings[[x]]$left, ])/sum(w.logit[[x]])
      } else {
        right.data[x,]
      }
    })

    right.extras.cor[is.na(right.extras.cor)] <- 0 # matrix of hvgs*cells
    right.data.cor[is.na(right.data.cor)] <- 0 # matrix of genes*cells
    colnames(right.extras.cor) <- row.names(right.extras)
    colnames(right.data.cor) <- row.names(right.data)

    merge.tree <- UPDATE(merge.tree, next.step$chosen,
                         data=rbind(left.data, t(right.data.cor)),
                         index=c(left.index, right.index),
                         restrict=.combine_restrict(left.data, left.restrict, t(right.data.cor), right.restrict),
                         origin=c(left.origin, right.origin),
                         extras=list(rbind(left.extras, t(right.extras.cor))))
  }
  if (!correct.all) {
    full.data <- .get_node_data(merge.tree)
  } else {
    full.data <- .get_node_extras(merge.tree)[[1]]
  }
  full.data
}

