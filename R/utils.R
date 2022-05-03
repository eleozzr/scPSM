.restricted_mnn <- function(left.data, left.restrict, right.data, right.restrict, k, prop.k=NULL, ...) {
  if (!is.null(left.restrict)) {
    left.data <- left.data[left.restrict,,drop=FALSE]
  } else {
    left.data <- left.data
  }

  if (!is.null(right.restrict)) {
    right.data <- right.data[right.restrict,,drop=FALSE]
  } else {
    right.data <- right.data
  }

  k1 <- .choose_k(k, prop.k, nrow(left.data))
  k2 <- .choose_k(k, prop.k, nrow(right.data))

  pairs <- BiocNeighbors::findMutualNN(left.data, right.data, k1=k1, k2=k2, ...)
  pairs$first <- .unrestrict_indices(pairs$first, left.restrict)
  pairs$second <- .unrestrict_indices(pairs$second, right.restrict)
  pairs
}

.choose_k <- function(k, prop.k, N) {
  if (is.null(prop.k)) {
    k
  } else {
    min(N, max(k, round(prop.k * N)))
  }
}

.unrestrict_indices <- function(index, restrict) {
  if (!is.null(restrict)) index <- restrict[index]
  index
}

# Methods for creating a tree with a pre-defined structure.
#' @importFrom utils tail relist
.create_tree_predefined <- function(batches, restrict, merge.order) {
  if (is.null(merge.order)) {
    merge.order <- seq_along(batches)
  }

  if (!is.list(merge.order) && length(merge.order) > 1L) {
    merge.tree <- list(merge.order[1], merge.order[2])
    for (i in tail(merge.order, -2L)) {
      merge.tree <- list(merge.tree, i)
    }
  } else {
    merge.tree <- merge.order
  }

  merge.tree <- .binarize_tree(merge.tree)

  # Checking validity of leaf identities.
  leaves <- unlist(merge.tree)
  if (!is.numeric(leaves)) {
    leaves <- match(leaves, names(batches))
  } else {
    leaves <- as.integer(leaves)
  }
  if (any(is.na(leaves)) || anyDuplicated(leaves) || any(leaves < 1) || any(leaves > length(batches))) {
    stop("invalid leaf nodes specified in 'merge.order'")
  }

  merge.tree <- relist(leaves, merge.tree)
  .fill_tree(merge.tree, batches, restrict)
}

.binarize_tree <- function(merge.tree) {
  if (!is.list(merge.tree) && length(merge.tree)==1L) {
    return(merge.tree)
  }

  N <- length(merge.tree)
  if (N==1L) {
    # Get rid of useless internal nodes with only one child.
    .binarize_tree(merge.tree[[1]])
  } else if (N==2L) {
    # Keep going.
    list(.binarize_tree(merge.tree[[1]]), .binarize_tree(merge.tree[[2]]))
  } else if (N > 2L) {
    # Progressive merge.
    current <- list(.binarize_tree(merge.tree[[1]]), .binarize_tree(merge.tree[[2]]))
    for (i in 3:N) {
      current <- list(current, .binarize_tree(merge.tree[[i]]))
    }
    current
  } else {
    stop("merge tree contains a node with no children")
  }
}

.fill_tree <- function(merge.tree, batches, restrict) {
  if (!is.list(merge.tree)) {
    val <- batches[[merge.tree]]
    return(MNN_treenode(index=merge.tree, data=val, restrict=restrict[[merge.tree]]))
  }
  if (length(merge.tree)!=2L) {
    stop("merge tree structure should contain two children per node")
  }
  merge.tree[[1]] <- .fill_tree(merge.tree[[1]], batches, restrict)
  merge.tree[[2]] <- .fill_tree(merge.tree[[2]], batches, restrict)
  merge.tree
}

#' @importFrom methods new
MNN_treenode <- function(index, data, restrict, origin=rep(index, nrow(data)), extras=list()) {
  new("MNN_treenode", index=index, data=data, restrict=restrict, origin=origin, extras=extras)
}

.get_node_index <- function(node) node@index

.get_node_data <- function(node) node@data

.get_node_origin <- function(node) node@origin

.get_node_restrict <- function(node) node@restrict

.get_node_extras <- function(node) node@extras

.add_out_batches_to_tree <- function(merge.tree, out.batches) {
  if (!is.list(merge.tree)) {
    merge.tree@extras <- list(out.batches[[.get_node_index(merge.tree)]])
    return(merge.tree)
  }
  merge.tree[[1]] <- .add_out_batches_to_tree(merge.tree[[1]], out.batches)
  merge.tree[[2]] <- .add_out_batches_to_tree(merge.tree[[2]], out.batches)
  merge.tree
}

.get_next_merge <- function(merge.tree, path=NULL) {
  if (!is.list(merge.tree[[1]]) && !is.list(merge.tree[[2]])) {
    list(left=merge.tree[[1]], right=merge.tree[[2]], chosen=path)
  } else if (is.list(merge.tree[[2]])) {
    .get_next_merge(merge.tree[[2]], c(path, 2L))
  } else {
    .get_next_merge(merge.tree[[1]], c(path, 1L))
  }
}

.update_tree <- function(merge.tree, path, ...) {
  if (length(path)==0L) {
    return(MNN_treenode(...))
  }
  merge.tree[[path[1]]] <- .update_tree(merge.tree[[path[[1]]]], path[-1], ...)
  merge.tree
}

.prepare_input_data <- function(batches, cos.norm.in, cos.norm.out, subset.row, correct.all) {
  nbatches <- length(batches)
  in.batches <- out.batches <- batches
  same.set <- TRUE

  # Subsetting to the desired subset of genes.
  if (!is.null(subset.row)) {
    subset.row <- .row_subset_to_index(batches[[1]], subset.row)
    if (identical(subset.row, seq_len(nrow(batches[[1]])))) {
      subset.row <- NULL
    } else {
      in.batches <- lapply(in.batches, "[", i=subset.row, , drop=FALSE) # Need the extra comma!
      if (correct.all) {
        same.set <- FALSE
      } else {
        out.batches <- in.batches
      }
    }
  }

  # Applying cosine normalization for MNN identification.
  # We use the L2 norm for the subsetted input to adjust the output,
  # to ensure that results are consistent regardless of the manner of subsetting.
  if (cos.norm.in) {
    norm.scaling <- vector("list", nbatches)
    for (b in seq_len(nbatches)) {
      current.in <- in.batches[[b]]
      cos.out <- batchelor::cosineNorm(current.in, mode="all")
      in.batches[[b]] <- cos.out$matrix
      norm.scaling[[b]] <- cos.out$l2norm
    }
  }
  if (cos.norm.out) {
    if (!cos.norm.in) {
      norm.scaling <- lapply(in.batches, batchelor::cosineNorm, mode="l2norm")
    }
    out.batches <- mapply(FUN=.apply_cosine_norm, out.batches, norm.scaling, SIMPLIFY=FALSE)
  }

  if (cos.norm.out!=cos.norm.in) {
    same.set <- FALSE
  }

  list(In=in.batches, Out=out.batches, Subset=subset.row, Same=same.set)
}

#' @importFrom S4Vectors normalizeSingleBracketSubscript
.row_subset_to_index <- function(x, index) {
  if (is.null(index)) {
    seq_len(nrow(x))
  } else {
    S4Vectors::normalizeSingleBracketSubscript(index, x)
  }
}

#' @importFrom scater normalizeCounts
.apply_cosine_norm <- function(x, l2) {
  l2 <- pmax(1e-8, l2) # protect against zero-L2.
  scater::normalizeCounts(x, size_factors=l2, center_size_factors=FALSE, log=FALSE)
}

.combine_restrict <- function(left.data, left.restrict, right.data, right.restrict) {
  if (is.null(left.restrict) && is.null(right.restrict)) {
    NULL
  } else {
    if (is.null(left.restrict)) {
      left.restrict <- seq_len(nrow(left.data))
    }
    if (is.null(right.restrict)) {
      right.restrict <- seq_len(nrow(right.data))
    }
    c(left.restrict, right.restrict + nrow(left.data))
  }
}
