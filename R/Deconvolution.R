#######################################
####### DECONVOLUTION FUNCTIONS #######
#######################################

############################################
#' Basis Matrix
#' @description Basis matrix construction
#' @name SCDC_basis
#' @param x ExpressionSet object for single cells
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @return a list of basis matrix, sum of cell-type-specific library size, sample variance matrix, basis matrix by mvw, mvw matrix.
#' @export
SCDC_basis <- function(x, ct.sub = NULL, ct.varname, sample){
  # select only the subset of cell types of interest
  if (is.null(ct.sub)){
    ct.sub <- unique(x@phenoData@data[,ct.varname])
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]
  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0,]
  # calculate sample mean & sample variance matrix: genes by cell types
  countmat <- exprs(x.sub)
  ct.id <- droplevels(as.factor(x.sub@phenoData@data[,ct.varname]))
  sample.id <- as.character(x.sub@phenoData@data[,sample])
  ct_sample.id <- paste(ct.id,sample.id, sep = '%')

  mean.mat <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y,1,sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))

  sigma <- sapply(unique(mean.id[,1]), function(id){
    y = mean.mat[,mean.id[,1] %in% id]

    #with there are not all samples for each cluster
    if(class(y)=='matrix')
      apply(y,1,var, na.rm = TRUE)
    else{
      print('Fixing problem with cluster and sample.')
      #I add zero to the first column in order to create a simulated data that should create a proportion of zero.
      y <-numeric(length(y))
      apply(cbind(y,numeric(length(y))),1,var, na.rm = TRUE)
    }
  })

  sum.mat2 <- sapply(unique(sample.id), function(sid){
    sapply(unique(ct.id), function(id){
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% sid])
      sum(y)/ncol(y)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  sum.mat <- rowMeans(sum.mat2, na.rm = T)

  basis <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # weighted basis matrix
  my.max <- function(x,...){
    y <- apply(x,1,max, na.rm = TRUE)
    y / median(y, na.rm = T)
  }

  # MATCH DONOR, CELLTYPE, GENES!!!!!!!!!!!!!!!!
  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid,
                   drop = FALSE]
      apply(y,1,var, na.rm=T)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)

  q15 <- apply(var.adj,2,quantile, probs = 0.15, na.rm =T)
  q85 <- apply(var.adj,2,quantile, probs = 0.85, na.rm =T)

  var.adj.q <- t(apply(var.adj, 1,
                       function(y){y[y<q15] <- q15[y<q15]
                       y[y>q85] <- q85[y>q85]
                       return(y)})) + 1e-4

  message("Creating Basis Matrix adjusted for maximal variance weight")
  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id){
    sid = unlist(strsplit(id,'%'))[2]
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q[,sid]), '/')
    apply(yy,1,sum, na.rm = TRUE)/sum(yy)
  })

  basis.mvw <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat.mvw)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # reorder columns
  basis.mvw <- basis.mvw[,ct.sub]
  sigma <- sigma[, ct.sub]
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]

  return(list(basis = basis, sum.mat = sum.mat,
              sigma = sigma, basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

#############################################
#' Basis matrix for single cells from one subject
#' @description Basis matrix construction for single cells from one subject
#' @name SCDC_basis_ONE
#' @param x ExpressionSet object for single cells
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @return a list of basis matrix, sum of cell-type-specific library size, sample variance matrix, basis matrix by mvw, mvw matrix.
#' @export
SCDC_basis_ONE <- function(x , ct.sub = NULL, ct.varname, sample){
  # select only the subset of cell types of interest
  if (is.null(ct.sub)){
    ct.sub <- unique(x@phenoData@data[,ct.varname])[!is.na(unique(x@phenoData@data[,ct.varname]))]
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]
  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(exprs(x.sub)) > 0,]
  # calculate sample mean & sample variance matrix: genes by cell types
  countmat <- exprs(x.sub)
  # ct.id <- droplevels(as.factor(x.sub@phenoData@data[,ct.varname]))
  ct.id <- x.sub@phenoData@data[,ct.varname]
  sample.id <- x.sub@phenoData@data[,sample]
  ct_sample.id <- paste(ct.id,sample.id, sep = '%')

  mean.mat <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y,1,sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))

  # by subj, then take avg????
  sum.mat2 <- sapply(unique(sample.id), function(sid){
    sapply(unique(ct.id), function(id){
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% sid])
      sum(y)/ncol(y)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)
  sum.mat <- rowMeans(sum.mat2, na.rm = T)

  basis <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat)*z)
    # id = unique(mean.id[,1])[1]
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # weighted basis matrix
  my.max <- function(x,...){
    y <- apply(x,1,max, na.rm = TRUE)
    y / median(y, na.rm = T)
  }

  # MATCH DONOR, CELLTYPE, GENES!!!!!!!!!!!!!!!!
  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid,
                   drop = FALSE]
      apply(y,1,var, na.rm=T)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)

  q15 <- apply(var.adj,2,quantile, probs = 0.15, na.rm =T)
  q85 <- apply(var.adj,2,quantile, probs = 0.85, na.rm =T)

  var.adj.q <- as.matrix(apply(var.adj, 1,
                               function(y){y[y<q15] <- q15[y<q15]
                               y[y>q85] <- q85[y>q85]
                               return(y)}) + 1e-4)

  message("Creating Basis Matrix adjusted for maximal variance weight")
  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q), '/')
    apply(yy,1,sum, na.rm = TRUE)/sum(yy)
  })

  basis.mvw <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat.mvw)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # reorder columns
  basis.mvw <- basis.mvw[,ct.sub]
  sigma <- NULL # in the one subject case, no variance is calculated.
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]

  return(list(basis = basis, sum.mat = sum.mat,
              sigma = sigma, basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

#################################
#' Clustering QC
#' @description Single cells Clustering QC
#' @name SCDC_qc
#' @import pheatmap
#' @param sc.eset ExpressionSet object for single cells
#' @param ct.varname variable name for 'cell type'
#' @param sample variable name for subject/sample
#' @param scsetname the name for the single cell dataset
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param arow annotation of rows for pheatmap
#' @param qcthreshold the probability threshold used to filter out questionable cells
#' @param generate.figure logical. If generate the heatmap by pheatmap or not. default is TRUE.
#' @return a list including: 1) a probability matrix for each single cell input; 2) a clustering QCed ExpressionSet object; 3) a heatmap of QC result.
#' @export
SCDC_qc <- function (sc.eset, ct.varname, sample, scsetname = "Single Cell",
                   ct.sub, iter.max = 1000, nu = 1e-04, epsilon = 0.01, arow =NULL,
                   qcthreshold = 0.7, generate.figure = F,
                   cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                   parallelize= F, core_number = NULL,
                   ...) {
  sc.basis = SCDC_basis(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  M.S <- sc.basis$sum.mat[ct.sub]
  xsc <- getCPM0(exprs(sc.eset)[rownames(sc.basis$basis.mvw),])
  N.sc <- ncol(xsc)

  m.basis <- sc.basis$basis.mvw[, ct.sub]
  sigma <- sc.basis$sigma[, ct.sub]
  valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(m.basis)) ==
                                                0) & (!is.na(M.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  m.basis <- m.basis[, valid.ct]
  M.S <- M.S[valid.ct]
  sigma <- sigma[, valid.ct]

  prop.qc <- NULL

  if(!parallelize){

    for (i in 1:N.sc) {
      #message("Begin iterative weighted estimation...")
      basis.temp <- m.basis
      xsc.temp <- xsc[, i]
      sigma.temp <- sigma
      ### weighting scheme:
      lm.qc <- nnls::nnls(A=basis.temp,b=xsc.temp)
      delta <- lm.qc$residuals
      wt.gene <- 1/(nu + delta^2 + colSums((lm.qc$x)^2*t(sigma.temp)))
      x.wt <- xsc.temp*sqrt(wt.gene)
      b.wt <- sweep(basis.temp,1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      prop.wt <- lm.wt$x/sum(lm.wt$x)
      delta <- lm.wt$residuals
      for (iter in 1:iter.max){
        wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x)^2*t(sigma.temp)))
        x.wt <- xsc.temp*sqrt(wt.gene)
        b.wt <- sweep(basis.temp,1,sqrt(wt.gene),"*")
        lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
        delta.new <- lm.wt$residuals
        prop.wt.new <- lm.wt$x/sum(lm.wt$x)
        if (sum(abs(prop.wt - prop.wt.new) < epsilon )){
          prop.wt <- prop.wt.new
          delta <- delta.new
          #message("Converged at iteration ", iter)
          break
        }
        prop.wt <- prop.wt.new
        delta <- delta.new
      }

      prop.qc <- rbind(prop.qc, prop.wt)

      if((i%%100) == 0.0){
        print(paste0('Progress: ', i, ' cells of ', N.sc, ' ', round(((i/N.sc)*100),2), '%'))
      }
    }

  }else{
    prop.qc <- qc_iteration_parallelized(m.basis, xsc, sigma, N.sc, nu, epsilon, iter.max, core_number)
  }


  # name col and row
  colnames(prop.qc) <- colnames(m.basis)
  rownames(prop.qc) <- colnames(xsc)
  if (generate.figure){
    library(pheatmap)
    heat.anno <- pheatmap(prop.qc, annotation_row = arow,
                          annotation_names_row=FALSE, show_rownames = F,
                          annotation_names_col=FALSE, cutree_rows = length(ct.sub),
                          color = cbPalette[2:4],
                          cluster_rows = T, cluster_cols = F)
  } else {
    heat.anno <- NULL
  }

  prop.qc.keep <- rowSums(prop.qc > qcthreshold) ==1 # truncated values -> F or T
  sc.eset.qc <- sc.eset[,prop.qc.keep]
  return(list(prop.qc = prop.qc, sc.eset.qc = sc.eset.qc, heatfig = heat.anno))
}

############################################
#' Paralellization of the for procedure for QC process
#' @description
#' @name iteration_over_clusters_parallelized
#' @param m.basis
#' @param xsc
#' @param N.sc
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param iter.max the maximum number of iteration in WNNLS
#' @param core_number Number of cores that the processors uses, for default uses # processors -1
#' @return
#' @export
qc_iteration_parallelized <- function(m.basis, xsc, sigma, N.sc, nu, epsilon, iter.max, core_number = NULL, ...){

  gc()

  print('Executing parallelized function...')

  library(foreach)
  library(doParallel)

  if(is.null(core_number)){
    # Calculate the number of cores
    no_cores <- detectCores() - 1
  }else{
    no_cores <- core_number
  }

  print(paste0('Number of cores to use: ', no_cores))

  # Initiate cluster
  temporal_directory <- paste0(getwd(), '/', 'temporal_log/')
  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  #create the directory
  dir.create(temporal_directory)

  registerDoParallel(no_cores)

  prop.qc <- NULL

  #definition of the parallel function
  prop.qc <- foreach(i = 1:N.sc,
                            .combine = rbind)  %dopar%
    {
        gc()

        message1 <- paste0("Begin iterative weighted estimation...")

        basis.temp <- m.basis
        xsc.temp <- xsc[, i]
        sigma.temp <- sigma
        ### weighting scheme:
        lm.qc <- nnls::nnls(A=basis.temp,b=xsc.temp)
        delta <- lm.qc$residuals
        wt.gene <- 1/(nu + delta^2 + colSums((lm.qc$x)^2*t(sigma.temp)))
        x.wt <- xsc.temp*sqrt(wt.gene)
        b.wt <- sweep(basis.temp,1,sqrt(wt.gene),"*")
        lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
        prop.wt <- lm.wt$x/sum(lm.wt$x)
        delta <- lm.wt$residuals
        for (iter in 1:iter.max){
          wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x)^2*t(sigma.temp)))
          x.wt <- xsc.temp*sqrt(wt.gene)
          b.wt <- sweep(basis.temp,1,sqrt(wt.gene),"*")
          lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
          delta.new <- lm.wt$residuals
          prop.wt.new <- lm.wt$x/sum(lm.wt$x)
          if (sum(abs(prop.wt - prop.wt.new) < epsilon )){
            prop.wt <- prop.wt.new
            delta <- delta.new
            #message("Converged at iteration ", iter)
            break
          }
          prop.wt <- prop.wt.new
          delta <- delta.new
        }

        if((i%%1000) == 0.0){
          progress <- paste0('Progress: ', i, ' cells of ', N.sc, ' ', round(((i/N.sc)*100),2), '%')
          message2 <- paste0("Begin iterative weighted estimation...", '\n', progress)
          save_log_file(temporal_directory, paste0('status_',i,'.log'), message1, message2)
        }

        prop.wt
    }
  stopImplicitCluster()

  #recover log files from temp directory
  list_log_files <- list.files(temporal_directory)
  print(paste0('Processing ', length(list_log_files), ' log files.'))

  for(counter in 1:length(list_log_files)){
    file_name_path <- paste0(temporal_directory, list_log_files[counter])
    print(paste0('Recovering log file: ', list_log_files[counter]))
    recovery_file_content <- readChar(file_name_path, file.info(file_name_path)$size)

    split_file <- unlist(strsplit(recovery_file_content,'\n'))
    for(counter in 1:length(split_file)){
      print(split_file[counter])
    }
  }

  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  print(paste0('parallelized process finished :', nrow(prop.qc)))

  gc()
  prop.qc
}


#################################
#' Clustering QC for single cells from one subject
#' @description Clustering QC for single cells from one subject
#' @name SCDC_qc_ONE
#' @param sc.eset ExpressionSet object for single cells
#' @param ct.varname variable name for 'cell type'
#' @param sample variable name for subject/sample
#' @param scsetname the name for the single cell dataset
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param arow annotation of rows for pheatmap
#' @param qcthreshold the probability threshold used to filter out questionable cells
#' @param generate.figure logical. If generate the heatmap by pheatmap or not. default is TRUE.
#' @return a list including: 1) a probability matrix for each single cell input; 2) a clustering QCed ExpressionSet object; 3) a heatmap of QC result.
#' @export
SCDC_qc_ONE <- function(sc.eset, ct.varname, sample, scsetname = "Single Cell",
                  ct.sub, iter.max = 1000, nu = 1e-04, epsilon = 0.01,
                    arow = NULL, weight.basis = F, qcthreshold = 0.7,
                  generate.figure = T,
                  cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                  ...){
  sc.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  if (weight.basis){
    basis.mvw <- sc.basis$basis.mvw[, ct.sub]
  } else {
    basis.mvw <- sc.basis$basis[, ct.sub]
  }

  xsc <- getCPM0(exprs(sc.eset))
  N.sc <- ncol(xsc)

  ALS.S <- sc.basis$sum.mat[ct.sub]
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]
  prop.est.mvw <- NULL

  # prop estimation for each sc sample:
  for (i in 1:N.sc) {
    xsc.i <- xsc[, i]*100 # why times 100 if not normalize???
    gene.use <- intersect(rownames(basis.mvw), names(xsc.i))
    basis.mvw.temp <- basis.mvw[gene.use,]
    xsc.temp <- xsc.i[gene.use]
    message(paste(colnames(xsc)[i], "has common genes", sum(xsc[, i] != 0), "..."))
    # first NNLS:

    lm <- nnls::nnls(A=basis.mvw.temp,b=xsc.temp)
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xsc.temp*sqrt(wt.gene)
    b.wt <- sweep(basis.mvw.temp,1,sqrt(wt.gene),"*")

    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals

    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xsc.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.temp,1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.new <- lm.wt$x/sum(lm.wt$x)

      if (sum(abs(prop.wt.new - prop.wt)) < epsilon){
        prop.wt <- prop.wt.new
        delta <- delta.new
        message("WNNLS Converged at iteration ", iter)
        break
      }
      prop.wt <- prop.wt.new
      delta <- delta.new
    }

    prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
  }
  colnames(prop.est.mvw) <- colnames(basis.mvw)
  rownames(prop.est.mvw) <- colnames(xsc)

  ### plot steps:
 if (generate.figure){
   heat.anno <- pheatmap(prop.est.mvw, annotation_row = arow,
                         annotation_names_row=FALSE, show_rownames = F,
                         annotation_names_col=FALSE, cutree_rows = length(ct.sub),
                         color = cbPalette[2:4],
                         cluster_rows = T, cluster_cols = F) #, main = scsetname
 } else {
   heat.anno <- NULL
 }

  prop.qc.keep <- rowSums(prop.est.mvw > qcthreshold) ==1 # truncated values -> F or T
  sc.eset.qc <- sc.eset[,prop.qc.keep]
  return(list(prop.qc = prop.est.mvw, sc.eset.qc = sc.eset.qc, heatfig = heat.anno))
}

######################################
#' Proportion estimation
#' @description Proportion estimation function for multi-subject case
#' @name SCDC_prop
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param truep true cell-type proportions for bulk samples if known
#' @param Transform_bisque The bulk sample transformation from bisqueRNA. Aiming to reduce the systematic difference between single cells and bulk samples.
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
SCDC_prop <- function (bulk.eset, sc.eset, ct.varname, sample, ct.sub, iter.max = 1000,
                       nu = 1e-04, epsilon = 0.01, truep = NULL, weight.basis = T,
                       Transform_bisque = F, ...)
{
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset)) > 0, , drop = FALSE]
  ct.sub <- intersect(ct.sub, unique(sc.eset@phenoData@data[,
                                                            ct.varname]))
  sc.basis <- SCDC_basis(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname,
                         sample = sample)
  commongenes <- intersect(rownames(sc.basis$basis.mvw), rownames(bulk.eset))
  if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
    stop("Too few common genes!")
  }
  message(paste("Used", length(commongenes), "common genes..."))
  if (weight.basis) {
    basis.mvw <- sc.basis$basis.mvw[commongenes, ct.sub]
  }
  else {
    basis.mvw <- sc.basis$basis[commongenes, ct.sub]
  }
  # link to bisqueRNA, bulk transformation method. https://github.com/cozygene/bisque
  if (Transform_bisque) {
    GenerateSCReference <- function(sc.eset, ct.sub) {
      cell.labels <- base::factor(sc.eset[[ct.sub]])
      all.cell.types <- base::levels(cell.labels)
      aggr.fn <- function(ct.sub) {
        base::rowMeans(Biobase::exprs(sc.eset)[,cell.labels == ct.sub, drop=F])
      }
      template <- base::numeric(base::nrow(sc.eset))
      sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
      return(sc.ref)
    }
    sc.ref <- GenerateSCReference(sc.eset, cell.types)[genes, , drop = F]
    ncount <- table(sc.eset@phenoData@data[, sample], sc.eset@phenoData@data[, ct.varname])
    true.prop <- ncount/rowSums(ncount, na.rm = T)
    sc.props <- round(true.prop[complete.cases(true.prop), ], 2)
    Y.train <- sc.ref %*% t(sc.props[, colnames(sc.ref)])
    dim(Y.train)
    X.pred <- exprs(bulk.eset)[commongenes, ]
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    template <- base::numeric(base::length(sample.names))
    base::names(template) <- sample.names
    SemisupervisedTransformBulk <- function(gene, Y.train, X.pred) {
      Y.train.scaled <- base::scale(Y.train[gene, , drop = T])
      Y.center <- base::attr(Y.train.scaled, "scaled:center")
      Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
      n <- base::length(Y.train.scaled)
      shrink.scale <- base::sqrt(base::sum((Y.train[gene, , drop = T] - Y.center)^2)/n + 1)
      X.pred.scaled <- base::scale(X.pred[gene, , drop = T])
      Y.pred <- base::matrix((X.pred.scaled * shrink.scale) +
                               Y.center, dimnames = base::list(base::colnames(X.pred),
                                                               gene))
      return(Y.pred)
    }
    Y.pred <- base::matrix(base::vapply(X = commongenes,
                                        FUN = SemisupervisedTransformBulk, FUN.VALUE = template,
                                        Y.train, X.pred, USE.NAMES = TRUE), nrow = base::length(sample.names))
    indices <- base::apply(Y.pred, MARGIN = 2, FUN = function(column) {
      base::anyNA(column)
    })
    if (base::any(indices)) {
      if (sum(!indices) == 0) {
        base::stop("Zero genes left for decomposition.")
      }
      Y.pred <- Y.pred[, !indices, drop = F]
      sc.ref <- sc.ref[!indices, , drop = F]
    }

    results <- base::as.matrix(base::apply(Y.pred, 1, function(b) {
      sol <- lsei::pnnls(sc.ref, b, sum = 1)
      return(sol$x)
    }))

    prop.est.mvw <- t(results)
    colnames(prop.est.mvw) <- colnames(sc.ref)
    rownames(prop.est.mvw) <- colnames(bulk.eset)
    yhat <- sc.ref %*% results
    colnames(yhat) <- colnames(bulk.eset)
    yobs <- exprs(bulk.eset)
    yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
    peval <- NULL
    if (!is.null(truep)) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw,
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
  } else {
    xbulk <- getCPM0(exprs(bulk.eset)[commongenes, ])
    sigma <- sc.basis$sigma[commongenes, ct.sub]
    ALS.S <- sc.basis$sum.mat[ct.sub]
    N.bulk <- ncol(bulk.eset)
    valid.ct <- (colSums(is.na(sigma)) == 0) & (colSums(is.na(basis.mvw)) ==
                                                  0) & (!is.na(ALS.S))
    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }
    message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    basis.mvw <- basis.mvw[, valid.ct]
    ALS.S <- ALS.S[valid.ct]
    sigma <- sigma[, valid.ct]
    prop.est.mvw <- NULL
    yhat <- NULL
    yhatgene.temp <- rownames(basis.mvw)
    for (i in 1:N.bulk) {
      basis.mvw.temp <- basis.mvw
      xbulk.temp <- xbulk[, i]*100
      sigma.temp <- sigma
      message(paste(colnames(xbulk)[i], "has common genes",
                    sum(xbulk[, i] != 0), "..."))
      lm <- nnls::nnls(A = basis.mvw.temp, b = xbulk.temp)
      delta <- lm$residuals
      wt.gene <- 1/(nu + delta^2 + colSums((lm$x * ALS.S)^2 *
                                             t(sigma.temp)))
      x.wt <- xbulk.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
      lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
      prop.wt <- lm.wt$x/sum(lm.wt$x)
      delta <- lm.wt$residuals
      for (iter in 1:iter.max) {
        wt.gene <- 1/(nu + delta^2 + colSums((lm.wt$x * ALS.S)^2 *
                                               t(sigma.temp)))
        x.wt <- xbulk.temp * sqrt(wt.gene)
        b.wt <- sweep(basis.mvw.temp, 1, sqrt(wt.gene), "*")
        lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
        delta.new <- lm.wt$residuals
        prop.wt.new <- lm.wt$x/sum(lm.wt$x)
        if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
          prop.wt <- prop.wt.new
          delta <- delta.new
          R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*%
                          as.matrix(lm.wt$x))/var(xbulk.temp)
          message("WNNLS Converged at iteration ",
                  iter)
          break
        }
        prop.wt <- prop.wt.new
        delta <- delta.new
      }
      R2 <- 1 - var(xbulk.temp - basis.mvw.temp %*% as.matrix(lm.wt$x))/var(xbulk.temp)
      prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
      yhat.temp <- basis.mvw.temp %*% as.matrix(lm.wt$x)
      yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
      yhat <- cbind(yhat[yhatgene.temp, ], yhat.temp[yhatgene.temp,
                                                     ])
    }
    colnames(prop.est.mvw) <- colnames(basis.mvw)
    rownames(prop.est.mvw) <- colnames(xbulk)
    colnames(yhat) <- colnames(xbulk)
    yobs <- exprs(bulk.eset)
    yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
    peval <- NULL
    if (!is.null(truep)) {
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw,
                          pest.names = c("SCDC"), select.ct = ct.sub)
    }
  }

  return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw,
              yhat = yhat, yeval = yeval, peval = peval))
}

############################################
#' Proportion estimation function for one-subject case
#' @description Proportion estimation function for one-subject case
#' @name SCDC_prop_ONE
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param sample variable name for subject/samples
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param truep true cell-type proportions for bulk samples if known
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
SCDC_prop_ONE <- function (bulk.eset, sc.eset, ct.varname, sample, truep = NULL,
                           ct.sub, iter.max = 2000, nu = 1e-10, epsilon = 0.01, weight.basis = T,
                           ...) {
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset)) > 0, , drop = FALSE]
  sc.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.sub,
                             ct.varname = ct.varname, sample = sample)
  if (weight.basis) {
    basis <- sc.basis$basis.mvw
  }
  else {
    basis <- sc.basis$basis
  }
  commongenes <- intersect(rownames(basis), rownames(bulk.eset))
  if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])) {
    stop("Too few common genes!")
  }
  message(paste("Used", length(commongenes), "common genes..."))
  basis.mvw <- basis[commongenes, ct.sub]
  xbulk <- getCPM0(exprs(bulk.eset)[commongenes, ])
  ALS.S <- sc.basis$sum.mat[ct.sub]
  N.bulk <- ncol(bulk.eset)
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))
  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]
  prop.est.mvw <- NULL
  yhat <- NULL
  yhatgene.temp <- rownames(basis.mvw)
  for (i in 1:N.bulk) {
    xbulk.temp <- xbulk[, i]
    message(paste(colnames(xbulk)[i], "has common genes",
                  sum(xbulk[, i] != 0), "..."))
    lm <- nnls::nnls(A = basis.mvw, b = xbulk.temp)
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xbulk.temp * sqrt(wt.gene)
    b.wt <- sweep(basis.mvw, 1, sqrt(wt.gene), "*")
    lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
    prop.wt <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals
    for (iter in 1:iter.max) {
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xbulk.temp * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw, 1, sqrt(wt.gene), "*")
      lm.wt <- nnls::nnls(A = b.wt, b = x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.new <- lm.wt$x/sum(lm.wt$x)
      if (sum(abs(prop.wt.new - prop.wt)) < epsilon) {
        prop.wt <- prop.wt.new
        delta <- delta.new
        message("WNNLS Converged at iteration ",
                iter)
        break
      }
      prop.wt <- prop.wt.new
      delta <- delta.new
    }
    prop.est.mvw <- rbind(prop.est.mvw, prop.wt)
    yhat.temp <- basis.mvw %*% as.matrix(lm.wt$x)
    yhatgene.temp <- intersect(rownames(yhat.temp), yhatgene.temp)
    yhat <- cbind(yhat[yhatgene.temp, ], yhat.temp[yhatgene.temp,
                                                   ])
  }
  colnames(prop.est.mvw) <- colnames(basis.mvw)
  rownames(prop.est.mvw) <- colnames(bulk.eset)
  colnames(yhat) <- colnames(bulk.eset)
  yobs <- exprs(bulk.eset)
  yeval <- SCDC_yeval(y = yobs, yest = yhat, yest.names = c("SCDC"))
  peval <- NULL
  if (!is.null(truep)) {
    if (all(rownames(truep) == rownames(prop.est.mvw))){
      peval <- SCDC_peval(ptrue = truep, pest = prop.est.mvw,
                          pest.names = c("SCDC"), select.ct = ct.sub)
    } else {
      message("Your input sample names for proportion matrix and bulk.eset do not match! Please make sure sample names match.")
    }

  }
  return(list(prop.est.mvw = prop.est.mvw, basis.mvw = basis.mvw,
              yhat = yhat, yeval = yeval, peval = peval))
}

############################################
#' Tree-guided proportion estimation
#' @description Proportion estimation function for multi-subject case, and apply tree-guided deconvolution
#' @name SCDC_prop_subcl_marker
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param fl.varname variable name for first-level 'meta-clusters'
#' @param sample variable name for subject/samples
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.fl.sub 'cell types' for first-level 'meta-clusters'
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param weight.basis logical, use basis matrix adjusted by MVW, default is T.
#' @param select.marker logical, select marker genes to perform deconvolution in tree-guided steps. Default is T.
#' @param markers A set of marker gene that input manully to be used in deconvolution. If NULL, then
#' @param marker.varname variable name of cluster groups when selecting marker genes. If NULL, then use ct.varname.
#' @param allgenes.fl logical, use all genes in the first-level deconvolution
#' @param pseudocount.use a constant number used when selecting marker genes, default is 1.
#' @param LFC.lim a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @param truep true cell-type proportions for bulk samples if known
#' @param iteration.use_final_foldchange TRUE/FALSE. If at the end the cluster has zero genes if this parameter is true, the boostraping is going to be calculated over the foldchange with <0.05, not with zero.
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
SCDC_prop_subcl_marker <- function(bulk.eset, sc.eset, ct.varname, fl.varname, sample,
                                   ct.sub = NULL, ct.fl.sub, iter.max = 3000, nu = 1e-04, epsilon = 0.001,
                                   weight.basis = T, truep = NULL,
                                   select.marker = T, markers = NULL, marker.varname = NULL, allgenes.fl = F,
                                   pseudocount.use = 1, LFC.lim = 0.5, parallelize= F, core_number = NULL, fix_number_genes = NULL, marker_gene_strategy = 'boostrap_outliers',
                                   iteration.minimun_number_markers = 28, iteration.use_maximum = FALSE, iteration.maximo_genes = 35, iteration.use_final_foldchange = FALSE,
                                   bootstrap.sample_size = NULL,...) {

  sc.eset.orig <- sc.eset

  if (is.null(ct.sub)){
    ct.sub <- unique(sc.eset@phenoData@data[,ct.varname])[!is.na(unique(sc.eset@phenoData@data[,ct.varname]))]
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  ct.fl.sub <- ct.fl.sub[!is.na(ct.fl.sub)]
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset))>0, , drop = FALSE]
  sc.eset <- sc.eset[,sc.eset@phenoData@data[,ct.varname] %in% ct.sub]
  message('SCDC basis for main cluster...')
  gc()
  sc.basis <- SCDC_basis(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  message('SCDC basis for secondary cluster (metacluster)...')
  sc.fl.basis <- SCDC_basis(x = sc.eset, ct.sub = ct.fl.sub[!is.na(ct.fl.sub)],
                            ct.varname = fl.varname, sample = sample)
  if (select.marker){
    if (is.null(marker.varname)){
      marker.varname <- ct.varname
    }
    # wilcox test on two groups of cells for marker gene selection... (refer to seurat::FindMarkers)
    countmat <- exprs(sc.eset)
    ct.group <- sc.eset@phenoData@data[,marker.varname]
    markers.wilcox <- NULL

    global.min.LFC <- 0
    global.genes.diff.LFC_max.top.cluster <- 0
    #this should be parametric and corresponds to the minimun number of markers on each cluster
    top_genes <- 20

    # u=1
    if(!parallelize){

      for(u in 1:length(unique(ct.group))){
        ct.group.temp <- ct.group == unique(ct.group)[u]

        #for checking progress for each cluster
        message(paste0('Selecting markers (1 cluster vs others ): ',u, '->',unique(ct.group)[u]))

        #Using custom_apply for sparse matrix
        group.1 <- custom_apply(X = countmat[,ct.group.temp],
                                MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                                    pseudocount.use))
        group.2 <- custom_apply(X = countmat[,! ct.group.temp],
                                MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                                    pseudocount.use))
        genes.diff <- rownames(sc.eset)[(group.1 - group.2) > LFC.lim]
        genes.diff.LFC_max <- max(group.1 - group.2)

        #We got the top N values (top_genes)
        genes.diff.LFC_max.top.cluster <- min(tail(sort(group.1 - group.2),top_genes))


        count.use <- countmat[rownames(sc.eset) %in% genes.diff,]

        #in order to check the global maximum value of LFC that it should be configured
        if(genes.diff.LFC_max < global.min.LFC | global.min.LFC==0){
          global.min.LFC <- genes.diff.LFC_max
        }

        #in order to check the global maximum value of LFC that it should be configured. With top N markers
        if(genes.diff.LFC_max.top.cluster < global.genes.diff.LFC_max.top.cluster | global.genes.diff.LFC_max.top.cluster==0){
          global.genes.diff.LFC_max.top.cluster <- genes.diff.LFC_max.top.cluster
        }

        #it's because there is just 1 gene
        if(class(count.use)=='numeric'){
          count.use <- t(count.use)
          rownames(count.use) <- c(genes.diff)
        }

        #check if there is none gene for the cluster
        if(nrow(count.use)==0){
          message(paste0("Zero Markers selected..."), 'Maximun LFC: ', genes.diff.LFC_max)
        }else{
          ##
          p_val <- sapply(1:nrow(count.use), function(x){
            wilcox.test(count.use[x,] ~ ct.group.temp)$p.value
          })

          p_val_adj <- p.adjust(p = p_val, method = "bonferroni",
                                n = nrow(count.use))
          markers.temp <- rownames(count.use)[p_val_adj < 0.05]

          #before adding, I will have the top n of the markers genes. THIS HAS TO BE DELETED FOR NORMAL BEHAVIOUR
          limit <- 5
          if(length(markers.temp) <= limit){
            top_markers.temp <- markers.temp
          }else{
            top_markers.temp <- markers.temp[limit,]
          }
          message("Markers MODIFIED selected(",length(top_markers.temp), "): ", paste(shQuote(top_markers.temp), collapse=", "))


          markers.wilcox <- c(markers.wilcox, markers.temp)

          #For checking the genes selected as markers for an specific cluster
          message("Markers selected(",length(markers.temp), "): ", paste(shQuote(markers.temp), collapse=", "))
        }
      }

    }else{

      message('Marker gene strategy: ', marker_gene_strategy)

      #There are three different implementation of marker gene function including the original(default)
      if(marker_gene_strategy=='default'){
        #call the original function with parallelization and some additional messages.
        markers.wilcox <- iteration_over_clusters_default_parallelized(ct.group, LFC.lim, sc.eset, pseudocount.use, core_number)
      }else if(marker_gene_strategy=='default_improved'){
        #call to the parallelized function
        markers.wilcox <- iteration_over_clusters_parallelized(ct.group, LFC.lim, sc.eset, top_genes, pseudocount.use, core_number, fix_number_genes)
      }else if(marker_gene_strategy=='wilcox_outlier'){
        #improve of library
        markers.wilcox <- iteration_over_clusters_parallelized_wilcox_outlier(sc.eset =sc.eset, ct.group = ct.group, core_number = core_number)
      }else if(marker_gene_strategy=='boostrap_outliers'){
        #improve with boostrap
        markers.wilcox <- iteration_over_clusters_parallelized_wilcox_boostrap(sc.eset = sc.eset.orig, ct.group = ct.group, core_number = core_number,
                                                                               iteration.minimun_number_markers = iteration.minimun_number_markers, iteration.use_maximum = iteration.use_maximum, iteration.maximo_genes = iteration.maximo_genes,
                                                                               iteration.use_final_foldchange = iteration.use_final_foldchange, bootstrap.sample_size = bootstrap.sample_size)
      }
    }

    markers <- unique(markers.wilcox)
    message("Global minimun LFC (1 marker): ", global.min.LFC)
    message("Global minimun LFC (Minimum ", top_genes, ' genes marker): ', global.genes.diff.LFC_max.top.cluster)
    message("Selected ",length(markers), " marker genes by Wilcoxon test...")
    message("Genes:", paste(shQuote(markers), collapse=", "))

  } # else need input of marker genes for clustering

  # match genes / cells first
  if (weight.basis){
    basis <- sc.basis$basis.mvw
    basis.fl <- sc.fl.basis$basis.mvw
  } else {
    basis <- sc.basis$basis
    basis.fl <- sc.fl.basis$basis
  }
  if (!is.null(markers)){
    commongenes <- Reduce(intersect, list(rownames(basis), rownames(bulk.eset), markers))
    commongenes.fl <- Reduce(intersect, list(rownames(basis.fl), rownames(bulk.eset), markers))
    non_commongenes  <- Reduce(setdiff, commongenes, markers)
    non_commongenes.fl <- Reduce(setdiff, commongenes.fl, markers)
  } else {
    commongenes <- intersect(rownames(basis), rownames(bulk.eset))
    commongenes.fl <- intersect(rownames(basis.fl), rownames(bulk.eset))
    # stop when few common genes exist...
    if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])){
      stop('Too few common genes!')
    }
  }

  message(paste("Used", length(commongenes), "common genes for all cell types, \n",
                "Used", length(commongenes.fl), "common genes for first level cell types..."))

  message('Genes that are not shared between SC and bulk:')
  message('*Cluster(',length(non_commongenes), '):',paste(shQuote(non_commongenes), collapse=", "))
  message('**Metacluster(', length(non_commongenes.fl), '):', paste(shQuote(non_commongenes.fl), collapse=", "))

  basis.mvw <- basis[commongenes, ct.sub]
  basis.mvw.fl <- basis.fl[commongenes.fl, ct.fl.sub]

  xbulk0 <- getCPM0(exprs(bulk.eset)[commongenes,])
  xbulk <- as.matrix(xbulk0) ## whether to normalize all /common genes
  colnames(xbulk) <- colnames(bulk.eset)
  xbulk1 <- getCPM0(exprs(bulk.eset)[commongenes.fl,])
  xbulk.fl <- as.matrix(xbulk1)

  ALS.S <- sc.basis$sum.mat[ct.sub]
  N.bulk <- ncol(bulk.eset)
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))
  ALS.S.fl <- sc.fl.basis$sum.mat[ct.fl.sub]
  valid.ct.fl <- (colSums(is.na(basis.mvw.fl)) == 0) & (!is.na(ALS.S.fl))

  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution...\n",
                "Used", sum(valid.ct.fl),"first level cell types ..."))

  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]
  basis.mvw.fl <- basis.mvw.fl[, valid.ct.fl]
  ALS.S.fl <- ALS.S[valid.ct.fl]

  prop.est <- NULL
  rsquared <- NULL

  # prop estimation for each bulk sample:
  for (i in 1:N.bulk) {
    # i=1
    xbulk.temp <- xbulk[, i] ## *1e3  will affect a little bit
    message(paste(colnames(xbulk)[i], "has common genes", sum(xbulk[, i] != 0), "..."))
    if (allgenes.fl){
      markers.fl <- names(xbulk.temp)
    } else {
      markers.fl <- Reduce(intersect, list(markers, names(xbulk.temp)))
    }

    # first level NNLS:
    lm <- nnls::nnls(A=basis.mvw.fl[markers.fl,],b=xbulk.temp[markers.fl])
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xbulk.temp[markers.fl] *sqrt(wt.gene)
    b.wt <- sweep(basis.mvw.fl[markers.fl,],1,sqrt(wt.gene),"*")

    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt.fl <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals

    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xbulk.temp[markers.fl] * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.fl[markers.fl,],1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.fl.new <- lm.wt$x/sum(lm.wt$x)

      if (sum(abs(prop.wt.fl.new - prop.wt.fl)) < epsilon){
        prop.wt.fl <- prop.wt.fl.new
        delta <- delta.new
        message("WNNLS for First level clusters Converged at iteration ", iter)
        break
      }
      prop.wt.fl <- prop.wt.fl.new
      delta <- delta.new
    }
    names(prop.wt.fl) <- colnames(basis.mvw.fl)

    # relationship between first level and overall
    rt <- table(sc.eset@phenoData@data[,ct.varname], sc.eset@phenoData@data[,fl.varname])
    rt <- rt[,ct.fl.sub]
    rt.list <- list()
    prop.wt <- NULL

    # prop.wt
    for (j in 1:ncol(rt)){ # for each first level cluster
      # j=1
      rt.list[[j]] <- rownames(rt)[rt[,j] >0]
      names(rt.list)[j] <- colnames(rt)[j]
      sub.cl <- rownames(rt)[rt[,j] >0]
      if (length(sub.cl) > 1 & prop.wt.fl[colnames(rt)[j]] > 0) {
        if (is.null(dim(prop.wt.fl))){
          # specify genes in xbulk.j ... first level genes ...
          xbulk.j <- basis.mvw.fl[,j]*prop.wt.fl[j] + (xbulk.temp - basis.mvw.fl %*% lm.wt$x)*prop.wt.fl[j]
        } else {
          xbulk.j <- basis.mvw.fl[,j]*prop.wt.fl[,j] + (xbulk.temp - basis.mvw.fl %*% lm.wt$x)*prop.wt.fl[,j]
        }

        markers.sl <- Reduce(intersect, list(markers, rownames(xbulk.j)))

        ##############################################################################
        # make markers.sub as a list, for each of the first-level intra clusters.
        ##############################################################################

        basis.sl <- basis.mvw[markers.sl,rownames(rt)[rt[,j] >0]]
        lm.sl <- nnls::nnls(A=basis.sl,b=xbulk.j[markers.sl,])
        delta.sl <- lm.sl$residuals
        wt.gene.sl <- 1/(nu + delta.sl^2)
        x.wt.sl <- xbulk.j[markers.sl,]*sqrt(wt.gene.sl)
        b.wt.sl <- sweep(basis.sl,1,sqrt(wt.gene.sl),"*")

        lm.wt.sl <- nnls::nnls(A=b.wt.sl, b=x.wt.sl)
        prop.wt.sl <- lm.wt.sl$x/sum(lm.wt.sl$x)
        delta.sl <- lm.wt.sl$residuals

        for (iter in 1:iter.max){
          wt.gene.sl <- 1/(nu + delta.sl^2)
          x.wt.sl <- xbulk.j[markers.sl,] * sqrt(wt.gene.sl)
          b.wt.sl <- sweep(basis.sl,1,sqrt(wt.gene.sl),"*")
          lm.wt.sl <- nnls::nnls(A=b.wt.sl, b=x.wt.sl)
          delta.sl.new <- lm.wt.sl$residuals
          prop.wt.sl.new <- lm.wt.sl$x/sum(lm.wt.sl$x)

          if (sum(abs(prop.wt.sl.new - prop.wt.sl)) < epsilon){
            prop.wt.sl <- prop.wt.sl.new
            delta.sl <- delta.sl.new
            message("WNNLS for Second level clusters",paste(shQuote(rownames(rt)[rt[,j] >0]), collapse=", ")," Converged at iteration ", iter)
            break
          }
          prop.wt.sl <- prop.wt.sl.new
          delta.sl <- delta.sl.new
        }
        names(prop.wt.sl) <- sub.cl
        prop.wt <- c(prop.wt, prop.wt.sl*prop.wt.fl[colnames(rt)[j]])
      } else if (length(sub.cl) == 1){
        # j=2
        prop.wt <- c(prop.wt, prop.wt.fl[colnames(rt)[j]])
      } else if (length(sub.cl) > 1 & prop.wt.fl[colnames(rt)[j]] == 0){
        prop.wt.sl <- rep(0, length(sub.cl))
        names(prop.wt.sl) <- sub.cl
        prop.wt <- c(prop.wt, prop.wt.sl)
      }

    }
    prop.est <- rbind(prop.est, prop.wt)
  }
  rownames(prop.est) <- colnames(bulk.eset)

  peval <- NULL
  if (!is.null(truep)){
    peval <- SCDC_peval(ptrue= truep, pest = prop.est, pest.names = c('SCDC'),
                       select.ct = ct.sub)
  }

  return(list(prop.est = prop.est, prop.wt.fl = prop.wt.fl, basis.mvw = basis.mvw, peval = peval,
              sc.basis = sc.basis, sc.fl.basis = sc.fl.basis))
}

#function to save message in a log file for each loop
save_log_file <- function(file_path, file_name, message1, message2){

  file_name_path <- paste0(file_path, '/', file_name)
  write(paste0(message1, '\n', message2), file=file_name_path,append=TRUE)
}

############################################
#' parallelization of the for the original procedure with some additional lines for:
#' @description  parallelization of the for the original procedure with some additional lines for: 1. Show the selected genes. 2. Add parallelization. 3. Add logging function in a temporal file.
#' @name iteration_over_clusters_default_parallelized
#' @param ct.group  List of clusters that will be analyzed
#' @param LFC.lim Foldchange limit for the comparison between each cluster amoung the others
#' @param sc.eset ExpressionSet object for single cells
#' @param pseudocount.use
#' @param core_number Number of cores that will be used for the process.
#' @return List with the marker genes for all selected clusters
#' @export
iteration_over_clusters_default_parallelized <- function(ct.group, LFC.lim, sc.eset, pseudocount.use, core_number = NULL,  ...){

  gc()
  countmat <- exprs(sc.eset)

  print('Executing parallelized function...')

  library(foreach)
  library(doParallel)

  if(is.null(core_number)){
    # Calculate the number of cores
    no_cores <- detectCores() - 1
  }else{
    no_cores <- core_number
  }

  print(paste0('Number of cores to use: ', no_cores))

  # Initiate cluster
  temporal_directory <- paste0(getwd(), '/', 'temporal_log/')
  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  #create the directory
  dir.create(temporal_directory)

  registerDoParallel(no_cores)

  iterations <- length(unique(ct.group))

  #definition of the parallel function
  markers.wilcox <- foreach(u = 1:iterations,
                            .combine = c)  %dopar%
    {
      gc()

      ct.group.temp <- ct.group == unique(ct.group)[u]

      #for checking progress for each cluster
      message1 <- paste0('Selecting markers (1 cluster vs others ): ',u, '->',unique(ct.group)[u])

      #Using custom_apply for sparse matrix
      group.1 <- custom_apply(X = countmat[,ct.group.temp],
                              MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                                  pseudocount.use))
      group.2 <- custom_apply(X = countmat[,! ct.group.temp],
                              MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                                  pseudocount.use))
      genes.diff <- rownames(sc.eset)[(group.1 - group.2) > LFC.lim]


      count.use <- countmat[rownames(sc.eset) %in% genes.diff,]

      #it's because there is just 1 gene
      if(class(count.use)=='numeric'){
        count.use <- t(count.use)
        rownames(count.use) <- c(genes.diff)
      }

      markers.temp = NULL

      #check if there is none gene for the cluster
      if(nrow(count.use) == 0){
        message2 <- paste0("Zero Markers selected...")
        save_log_file(temporal_directory, paste0(unique(ct.group)[u],'.log'), message1, message2)
      }else{
        ##
          p_val <- sapply(1:nrow(count.use), function(x){
            wilcox.test(count.use[x,] ~ ct.group.temp)$p.value
          })

          p_val_adj <- p.adjust(p = p_val, method = "bonferroni",
                                n = nrow(count.use))

          markers.temp <- rownames(count.use)[p_val_adj < 0.05]
        }

      #For checking the genes selected as markers for an specific cluster
      message2 <- paste0("Markers selected(",length(markers.temp), "): ", paste(shQuote(markers.temp), collapse=", "))
      save_log_file(temporal_directory, paste0(unique(ct.group)[u],'.log'), message1, message2)

      markers.temp
    }
  stopImplicitCluster()

  #recover log files from temp directory
  list_log_files <- list.files(temporal_directory)
  print(paste0('Processing ', length(list_log_files), ' log files.'))

  for(counter in 1:length(list_log_files)){
    file_name_path <- paste0(temporal_directory, list_log_files[counter])
    print(paste0('Recovering log file: ', list_log_files[counter]))
    recovery_file_content <- readChar(file_name_path, file.info(file_name_path)$size)

    split_file <- unlist(strsplit(recovery_file_content,'\n'))
    for(counter in 1:length(split_file)){
      print(split_file[counter])
    }
  }

  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  print(paste0('Number unique marker genes:', length(unique(markers.wilcox))))
  gc()
  markers.wilcox
}

############################################
#' Procedure with multiples modifications created through all the experimental process.
#' @description It contains some experiments but is not the final version. It is just saved for documentation purposes.
#' @name iteration_over_clusters_parallelized
#' @param ct.group  List of clusters that will be analyzed
#' @param LFC.lim a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @param sc.eset ExpressionSet object for single cells
#' @param top_genes Number of top genes that we want to consider to check the maximum LFC for each cluster. This is just for informative purposes.
#' @param pseudocount.use
#' @param core_number Number of cores that will be used for the process.
#' @param fix_number_genes  Maximum number of marker genes that has to be returned by cluster
#' @return List with the marker genes for all selected clusters
#' @export
iteration_over_clusters_parallelized <- function(ct.group, LFC.lim, sc.eset, top_genes, pseudocount.use, core_number = NULL, fix_number_genes = NULL,  ...){

  gc()
  countmat <- exprs(sc.eset)
  global.min.LFC <- 0
  genes.diff.LFC_max <- 0
  genes.diff.LFC_max.top.cluster <- 0
  global.genes.diff.LFC_max.top.cluster <- 0

  print('Executing parallelized function...')

  library(foreach)
  library(doParallel)

  if(is.null(core_number)){
    # Calculate the number of cores
    no_cores <- detectCores() - 1
  }else{
    no_cores <- core_number
  }

  print(paste0('Number of cores to use: ', no_cores))

  # Initiate cluster
  temporal_directory <- paste0(getwd(), '/', 'temporal_log/')
  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  #create the directory
  dir.create(temporal_directory)

  registerDoParallel(no_cores)

  iterations <- length(unique(ct.group))
  #iterations <- 1

  #definition of the parallel function
  markers.wilcox <- foreach(u = 1:iterations,
          .combine = c)  %dopar%
  {
    gc()

    global.genes.diff.LFC_max.top.cluster <- 0

    ct.group.temp <- ct.group == unique(ct.group)[u]

    #for checking progress for each cluster
    message1 <- paste0('Selecting markers (1 cluster vs others ): ',u, '->',unique(ct.group)[u])

    #Using custom_apply for sparse matrix
    group.1 <- custom_apply(X = countmat[,ct.group.temp],
                            MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                                pseudocount.use))
    group.2 <- custom_apply(X = countmat[,! ct.group.temp],
                            MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                                pseudocount.use))
    genes.diff <- rownames(sc.eset)[(group.1 - group.2) > LFC.lim]

    genes.diff.LFC_max <- max(group.1 - group.2)

    ###########################################
    print(paste0('LFC.lim: ', LFC.lim))
    genes.diff.LFC <- cbind(names=rownames(sc.eset), LFC=as.numeric((group.1 - group.2)))
    genes.diff.LFC <- genes.diff.LFC[genes.diff.LFC[,2] > LFC.lim,]
    genes.diff <- genes.diff.LFC[,1]


    ########## doing the foldchange for all genes with wilcox and bonferroni analysis
      save_wilcox_analysis = FALSE
    if(save_wilcox_analysis){
      #save the foldchange values for the genes
      foldChange <- (group.1 - group.2)

      p_val.aux <- sapply(1:nrow(countmat), function(x){
        wilcox.test(countmat[x,] ~ ct.group.temp)$p.value
      })

      p_val_adj.aux <- p.adjust(p = p_val.aux, method = "bonferroni",
                                n = nrow(countmat))

      foldchange_pvalue_padjust <- data.frame(gene_name = rownames(countmat), foldchange = foldChange, p_val = p_val.aux, p_val_adj = p_val_adj.aux)
      rownames(foldchange_pvalue_padjust) <- rownames(countmat)

      save(foldchange_pvalue_padjust, file = paste0(getwd(), '/', 'cluster_foldchange/', unique(ct.group)[u]))
    }
    ##########    ##########    ##########



    #We got the top N values (top_genes)
    genes.diff.LFC_max.top.cluster <- min(tail(sort(group.1 - group.2),top_genes))

    count.use <- countmat[rownames(sc.eset) %in% genes.diff,]


    #in order to check the global maximum value of LFC that it should be configured
    if(genes.diff.LFC_max < global.min.LFC | global.min.LFC==0){
      global.min.LFC <- genes.diff.LFC_max
    }

    #in order to check the global maximum value of LFC that it should be configured. With top N markers
    if(genes.diff.LFC_max.top.cluster < global.genes.diff.LFC_max.top.cluster | global.genes.diff.LFC_max.top.cluster==0){
      global.genes.diff.LFC_max.top.cluster <- genes.diff.LFC_max.top.cluster
    }

    #it's because there is just 1 gene
    if(class(count.use)=='numeric'){
      count.use <- t(count.use)
      rownames(count.use) <- c(genes.diff)
    }

    #check if there is none gene for the cluster
    if(nrow(count.use)==0){
      message2 <- paste0("Zero Markers selected...", 'Maximun LFC: ', genes.diff.LFC_max)
      save_log_file(temporal_directory, paste0(unique(ct.group)[u],'.log'), message1, message2)
    }else{

      p_val <- sapply(1:nrow(count.use), function(x){
        wilcox.test(count.use[x,] ~ ct.group.temp)$p.value
      })

      p_val_adj <- p.adjust(p = p_val, method = "bonferroni",
                            n = nrow(count.use))

      genes_p_val_adj <- cbind(gene_name = rownames(count.use), p_val_adj=p_val_adj)


      markers.temp <- rownames(count.use)[p_val_adj < 0.05]


      ################################## Part for have the more significant genes but get those in the top base on the fold change

      #joining
      genes_LFC = data.frame(gene_name = genes.diff.LFC[,1], LFC = as.numeric(genes.diff.LFC[,2]))
      # data frame 2
      genes_p_value_adj = data.frame(gene_name = as.character(genes_p_val_adj[,1]), p_val_adj = as.numeric(genes_p_val_adj[,2]))
      genes_p_value_adv_LFC <-merge(x=genes_LFC,y=genes_p_value_adj,by="gene_name")

      genes_p_value_adv_LFC.filtered <- genes_p_value_adv_LFC[genes_p_value_adv_LFC$p_val_adj < 0.05, ]
      #order by LFC or P VAL ADJ
      genes_p_value_adv_LFC.filtered.ordered <- genes_p_value_adv_LFC.filtered[order(-genes_p_value_adv_LFC.filtered$LFC),]

      print(paste0('markers: ', length(markers.temp), ', new markers: ', length(genes_p_value_adv_LFC.filtered.ordered)))

      markers.temp <- as.character(genes_p_value_adv_LFC.filtered.ordered$gene_name)

      #print('marker with all the data....')
      print(genes_p_value_adv_LFC.filtered.ordered)

      ##################################


      ###########################checking how many genes filtering p-value-adjusted=0
      genes_filtered <- genes_p_value_adv_LFC.filtered.ordered[genes_p_value_adv_LFC.filtered.ordered$p_val_adj==0,]
      message3 <- paste0("Markers GENES with P-VALUE-ADJ=0:(",nrow(genes_filtered), ", min-FC:" , min(genes_filtered$LFC) , ", max-FC:" ,max(genes_filtered$LFC) ,  "): ", paste(shQuote(genes_filtered$gene_name), collapse=", "))
      message1 <- paste0(message1, '\n', message3)
      ###########################checking how many genes filtering p-value-adjusted=0


      top_markers.temp <- NULL

      #before adding, I will have the top n of the markers genes. THIS HAS TO BE DELETED FOR NORMAL BEHAVIOUR
      if(!is.null(fix_number_genes)){
        print(paste0('lenght markers: ', length(markers.temp), ' fix number genes:', fix_number_genes))
        if(length(markers.temp) <  fix_number_genes){
          top_markers.temp <- markers.temp
        }else{

          cluster_more_genes <- paste0('cluster_',c())
          cluster_more_genes_2 <- paste0('cluster_',c())

          #change some cluster with one more gene
          if(unique(ct.group)[u] %in% cluster_more_genes){
            top_markers.temp <- markers.temp[1:(fix_number_genes+10)]
          }else if(unique(ct.group)[u] %in% cluster_more_genes_2)
          {
            top_markers.temp <- markers.temp[1:(fix_number_genes-20)]
          }else if(unique(ct.group)[u]=='cluster_xxx'){
            top_markers.temp <- markers.temp[1:(fix_number_genes+50)]
          }else{
            top_markers.temp <- markers.temp[1:(fix_number_genes)]
          }
        }

        #replace the temp por the top n
        markers.temp <- top_markers.temp
      }

      #For checking the genes selected as markers for an specific cluster
      message2 <- paste0('Maximun LFC intra-cluster: ', genes.diff.LFC_max, '\n', "Markers selected(",length(markers.temp), "): ", paste(shQuote(markers.temp), collapse=", "))
      save_log_file(temporal_directory, paste0(unique(ct.group)[u],'.log'), message1, message2)


      markers.temp
    }
  }
  stopImplicitCluster()

  #recover log files from temp directory
  list_log_files <- list.files(temporal_directory)
  print(paste0('Processing ', length(list_log_files), ' log files.'))

  for(counter in 1:length(list_log_files)){
    file_name_path <- paste0(temporal_directory, list_log_files[counter])
    print(paste0('Recovering log file: ', list_log_files[counter]))
    recovery_file_content <- readChar(file_name_path, file.info(file_name_path)$size)

    split_file <- unlist(strsplit(recovery_file_content,'\n'))
    for(counter in 1:length(split_file)){
      print(split_file[counter])
    }
  }

  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  print(paste0('Number unique marker genes:', length(unique(markers.wilcox))))

  gc()
  markers.wilcox
}


############################################
#' Marker gene selection with parallelization using wilcox method from the Seurat library
#' @description  Marker gene selection with parallelization using wilcox method from the Seurat library
#' @name iteration_over_clusters_parallelized_wilcox
#' @param sc.eset ExpressionSet object for single cells
#' @param ct.group  List of clusters that will be analyzed
#' @param core_number Number of cores that will be used for the process.
#' @param LFC.lim a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @param minimun_number_markers  If the variable minimum_genes_pipeline is active and also if in the first step with the wilcox text there are at least this genes, the process is not going to use outlier detection with dbscan
#' @param minimum_genes_pipeline  TRUE/FALSE. Parameter that enable the minimum gene process using outlier detection. If is FALSE the process is not going to run and the genes would be the ones found for wilcox test
#' @param add_zero_p_val_adj  TRUE/FALSE. If after the dbscan we want to add the genes that were found with the normal wilcox detection using the threshold (normally 0 or 0.05)
#' @param minimum_p_val_adj p-value-adjusted for the wilcox process. Could be zero or 0.05
#' @param maximo_genes  If we want to limitate the maximum number of marker genes at the end of the process.
#' @return List with the marker genes for all selected clusters
#' @export
iteration_over_clusters_parallelized_wilcox_outlier <- function(sc.eset, ct.group, core_number = NULL, LFC.lim=0.20,  minimun_number_markers = 20, minimum_genes_pipeline = TRUE,
                                                                add_zero_p_val_adj = TRUE, minimum_p_val_adj=0, maximo_genes=20,...){

  gc()

  pbmc <- CreateSeuratObject(counts = sc.eset@assayData$exprs,
                                          project = "Deconvolution_bulk_brain_data",
                                          assay = "RNA")

  Idents(pbmc) <- sc.eset$cluster_normalized
  cells_by_cluster <- as.matrix(table(Idents(pbmc)))

  print('Executing parallelized function...')

  library(foreach)
  library(doParallel)
  library(Seurat)
  library("fpc")

  if(is.null(core_number)){
    # Calculate the number of cores
    no_cores <- detectCores() - 1
  }else{
    no_cores <- core_number
  }

  print(paste0('Number of cores to use: ', no_cores))

  # Initiate cluster
  temporal_directory <- paste0(getwd(), '/', 'temporal_log/')
  print(paste0("Temporal folder: ", temporal_directory))
  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  #create the directory
  dir.create(temporal_directory)

  registerDoParallel(no_cores)

  ct.group.final <- unique(ct.group)

  iterations <- length(ct.group.final)

  #definition of the parallel function
  markers.wilcox <- foreach(u = 1:iterations,
                            .combine = c)  %dopar%
  {
    gc()
    cluster_name <- ct.group.final[u]

    #calculate number of cell for the cluster
    number_cells <- cells_by_cluster[rownames(cells_by_cluster)==cluster_name,]

    print(paste0('Cluster: ', cluster_name, ' Number of cells: ', number_cells))

    markers.temp = NULL
    message1 <- ''
    message2 <- ''

    if(cluster_name %in% c('cluster_to_be_filtered')){
      message1 <- paste0('The cluster: ', cluster_name, ' has too many problems')
    }else if(number_cells <= 3){
      message1 <- paste0("The cluster ", cluster_name ," doesn't have enough cells: ", number_cells)
    }else{

      markers.wilcox <- FindMarkers(pbmc, ident.1 = cluster_name, ident.2 = NULL,only.pos=TRUE, logfc.threshold = LFC.lim, verbose = F)

      #I added the verification of NA because the cluster 53 doesn't have any pvalues. This is going to apply for other cluster with the same behaviour
      filtered.wilcox <- rownames(markers.wilcox[markers.wilcox$p_val_adj <= minimum_p_val_adj & !is.na(markers.wilcox$p_val_adj),])

      #if the cluster doesn't have any marker
      if((length(filtered.wilcox)==0 | length(filtered.wilcox) < minimun_number_markers) & minimum_genes_pipeline){

        message1 <- paste0(cluster_name, "(ORIGINAL:", length(filtered.wilcox), ') -> Number of cells: ', number_cells)
        message1 <- paste0(message1, ', Genes: ',paste(shQuote(filtered.wilcox), collapse=", "))

        # Compute DBSCAN using fpc package
        set.seed(123)
        if(minimum_p_val_adj==0){
          db <- fpc::dbscan(markers.wilcox$avg_logFC, eps = 0.02, MinPts = 4)
        }else{
          markers.wilcox <- markers.wilcox[markers.wilcox$p_val_adj <= minimum_p_val_adj,]
          db <- fpc::dbscan(markers.wilcox$avg_logFC, eps = 0.02, MinPts = 4)
        }

        #select genes in the cluster 0
        marker_genes.temp <- as.character(rownames(markers.wilcox)[db$cluster==0])

        #if is true, the zero_p_val_adj will be added to the dbscan genes, false ONLY dbscan genes.
        if(add_zero_p_val_adj){
          #add the genes to the current ones and deduplique
          filtered.wilcox <- c(filtered.wilcox, marker_genes.temp)
          filtered.wilcox <- unique(filtered.wilcox)
        }else{
          filtered.wilcox <- marker_genes.temp
        }

        message1 <- paste0(message1, cluster_name, " (FIXED:", length(marker_genes.temp), ')')
        message1 <- paste0(message1, ', Genes: ',paste(shQuote(marker_genes.temp), collapse=", "))

      }else{
        #for checking progress for each cluster
        message1 <- paste0(cluster_name, "(", length(filtered.wilcox), ') -> Number of cells: ', number_cells)
      }

      message2 <- paste0('Final Genes', " (", length(filtered.wilcox), '):',paste(shQuote(filtered.wilcox), collapse=", "))


      #have the top N of genes
      if(length(filtered.wilcox) < maximo_genes & maximo_genes != 0){
        markers.temp <- filtered.wilcox[1:maximo_genes]
      }else{
        markers.temp <- filtered.wilcox
      }

    }

    #For checking the genes selected as markers for an specific cluster
    save_log_file(temporal_directory, paste0(cluster_name,'.log'), message1, message2)
    print(paste0(message1, ' ', message2))

    markers.temp
  }

  stopImplicitCluster()

  #recover log files from temp directory
  list_log_files <- list.files(temporal_directory)
  print(paste0('Processing ', length(list_log_files), ' log files.'))

  for(counter in 1:length(list_log_files)){
    file_name_path <- paste0(temporal_directory, list_log_files[counter])
    print(paste0('Recovering log file: ', list_log_files[counter]))
    recovery_file_content <- readChar(file_name_path, file.info(file_name_path)$size)

    split_file <- unlist(strsplit(recovery_file_content,'\n'))
    for(counter in 1:length(split_file)){
      print(split_file[counter])
    }
  }

  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  print(paste0('Number unique marker genes:', length(unique(markers.wilcox))))

  gc()
  markers.wilcox
}

############################################
#' Marker gene selection with parallelization using wilcox method from the Seurat library
#' @description  Marker gene selection with parallelization using wilcox method from the Seurat library
#' @name iteration_over_clusters_parallelized_wilcox_boostrap
#' @param sc.eset ExpressionSet object for single cells
#' @param ct.group  List of clusters that will be analyzed
#' @param core_number Number of cores that will be used for the process.
#' @param LFC.lim Fa threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @param minimun_number_markers  If in the first step with the wilcox text there are at least this genes, the process is not going to use outlier detection with dbscan
#' @param use_maximum  TRUE/FALSE. If the process selects just a fixed number of genes.
#' @param iteration.use_final_foldchange TRUE/FALSE. If at the end the cluster has zero genes if this parameter is true, the boostraping is going to be calculated over the foldchange with <0.05, not with zero.
#' @param maximo_genes  If we want to limitate the maximum number of marker genes at the end of the process.
#' @return List with the marker genes for all selected clusters
#' @export
iteration_over_clusters_parallelized_wilcox_boostrap <- function(sc.eset, ct.group, core_number = NULL, LFC.lim = 0.20,
                                                                 iteration.minimun_number_markers = 28, iteration.use_maximum = TRUE, iteration.maximo_genes = 35, iteration.use_final_foldchange = FALSE,
                                                                 bootstrap.sample_size = NULL,...){

  print(paste0('Parameters for iteration: ', "iteration.minimun_number_markers:", iteration.minimun_number_markers, ", iteration.use_maximum:", iteration.use_maximum, " iteration.maximo_genes:", iteration.maximo_genes))

  gc()
  library(Seurat)
  seurat_object <- CreateSeuratObject(counts = sc.eset@assayData$exprs,
                             project = "Deconvolution_bulk_brain_data",
                             assay = "RNA")

  Idents(seurat_object) <- sc.eset$cluster_normalized
  cells_by_cluster <- as.matrix(table(Idents(seurat_object)))

  print('Executing parallelized function...')

  library(foreach)
  library(doParallel)
  library(Seurat)
  library(fpc)

  if(is.null(core_number)){
    # Calculate the number of cores
    no_cores <- detectCores() - 1
  }else{
    no_cores <- core_number
  }

  print(paste0('Number of cores to use: ', no_cores))

  # Initiate cluster
  temporal_directory <- paste0(getwd(), '/', 'temporal_log/')
  print(paste0("Temporal folder: ", temporal_directory))
  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  #create the directory
  dir.create(temporal_directory)

  registerDoParallel(no_cores)

  ct.group.final <- unique(ct.group)
  iterations <- length(ct.group.final)

  print(paste0("Cluster order: ", paste(shQuote(ct.group.final), collapse=", ")))

  markers.wilcox.final <- NULL

  #definition of the parallel function
  markers.wilcox.final <- foreach(u = 1:iterations,
                            .combine = c)  %dopar%
    {
      rm(.Random.seed, envir=globalenv())
      set.seed(Sys.time())

      gc()
      library("stringr")
      cluster_name <- ct.group.final[u]
      cluster_number <- as.numeric(str_remove(cluster_name, 'cluster_'))

      #calculate number of cell for the cluster
      number_cells <- as.numeric(cells_by_cluster[rownames(cells_by_cluster)==cluster_name,])

      markers.temp = NULL
      message1 <- ''
      message2 <- ''

      message1 <- (paste0('Cluster: ', cluster_name, ' Number of cells: ', number_cells))
      print(message1)

      markers.wilcox <- NULL
      check_boostrap = F

      if(number_cells <= 3){
        message1 <- paste0(message1, '\n', "The cluster ", cluster_name ," doesn't have enough cells: ", number_cells)
        print(message1)
        check_boostrap = FALSE
      }else{

        its_well_calculated <- TRUE

        result = tryCatch({
          markers.wilcox <- FindMarkers(seurat_object, ident.1 = cluster_name, ident.2 = NULL, only.pos = TRUE, logfc.threshold = LFC.lim, verbose = F)
        }, warning = function(w) {
          print(paste0("warning", w))
        }, error = function(e) {
          print(paste0("error", e))
          its_well_calculated <- FALSE
        }, finally = {
        })

        if(its_well_calculated & !is.null(markers.wilcox)){

          markers.wilcox.filtered <- markers.wilcox[markers.wilcox$p_val_adj == 0 & !is.na(markers.wilcox$p_val_adj) & !is.na(rownames(markers.wilcox)) & rownames(markers.wilcox) != 'NA',]
          markers.temp <- rownames(markers.wilcox.filtered)

          #for checking progress for each cluster
          message1 <- paste0(message1, '\n', cluster_name, "(", length(markers.temp), ') -> Number of cells: ', number_cells, ' *****Genes: ', paste(shQuote(markers.temp), collapse=", "))

          #add other markers with boostrap
          if(nrow(markers.wilcox.filtered) <= iteration.minimun_number_markers){
            check_boostrap = T
          }else{
            check_boostrap = F
          }
        }else{
          check_boostrap = F
        }
      }

      if(check_boostrap){
        #call the analizer
        bootstrap_genes = bootstrap_gene_finder(cluster_number, seurat_object, bootstrap.sample_size = bootstrap.sample_size)
        message1 <- paste0(message1, '\n', 'Boostrap generation: ', "(", length(bootstrap_genes), ') -> Number of cells: ', number_cells, '. *****Genes: ', paste(shQuote(bootstrap_genes), collapse=", "))
        markers.temp <- c(markers.temp, bootstrap_genes)
        markers.temp <- unique(markers.temp)
      }

      #I add normal behaviour if there are zero genes with boostraping strategy (< iteration.minimun_number_markers)
      if(length(markers.temp) < iteration.minimun_number_markers & iteration.use_final_foldchange){

        #call the analizer with  min.p.adj_value = 0.05 that corresponds to the original version of the algorithm
        foldchange_genes = bootstrap_gene_finder(cluster_number, seurat_object, bootstrap.sample_size = bootstrap.sample_size, min.p.adj_value = 0.05)
        message1 <- paste0(message1, '\n', 'Foldchange generation: ', "(", length(foldchange_genes), ') -> Number of cells: ', number_cells, '. *****Genes: ', paste(shQuote(foldchange_genes), collapse=", "))
        markers.temp <- c(markers.temp, foldchange_genes)
        markers.temp <- unique(markers.temp)
      }

      #have the top N of genes
      if(length(markers.temp) >= iteration.maximo_genes & iteration.use_maximum){
        markers.temp <- markers.temp[1:iteration.maximo_genes]
      }

      message2 <- paste0('\n','Final Genes', " (", length(markers.temp), '):',paste(shQuote(markers.temp), collapse=", "), "\n\n")

      #For checking the genes selected as markers for an specific cluster
      save_log_file(temporal_directory, paste0(cluster_name,'.log'), message1, message2)
      print(paste0(message1, ' ', message2))

      markers.temp
    }

  stopImplicitCluster()

  #recover log files from temp directory
  list_log_files <- list.files(temporal_directory)
  print(paste0('Processing ', length(list_log_files), ' log files.'))

  for(counter in 1:length(list_log_files)){
    file_name_path <- paste0(temporal_directory, list_log_files[counter])
    print(paste0('Recovering log file: ', list_log_files[counter]))
    recovery_file_content <- readChar(file_name_path, file.info(file_name_path)$size)

    split_file <- unlist(strsplit(recovery_file_content,'\n'))
    for(counter in 1:length(split_file)){
      print(split_file[counter])
    }
  }

  if(file.exists(temporal_directory)){
    #delete the directory
    unlink(temporal_directory, recursive = TRUE)
  }

  print(paste0('Number unique marker genes:', length(unique(markers.wilcox.final))))

  gc()
  markers.wilcox.final
}


############################################
#' Function that perform a boostrap process from one cluster over a set of other clusters applying at the end outlier analysis with dbscan.
#' @description  Function that perform a boostrap process from one cluster over a set of other clusters applying at the end outlier analysis with dbscan.
#' @name bootstrap_gene_finder
#' @param LFC.lim a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @param cluster_number  Cluster number that the function has to process (TODO create a generalization for dataset with cluster with different names)
#' @param seurat_object  Seural object for single cells that enable us to apply the function FindMarkers
#' @param use_limit TRUE/FALSE. If the process selects just a fixed number of genes.
#' @param maximo_genes Maximum number of genes that have to be selected.
#' @param bootstrap_number  How many boostraping loops the algorithm is going to execute
#' @param min.p.adj_value   Value that should be greater or equal to zero. Normally it is 0.05
#' @return List with the marker genes
#' @export
bootstrap_gene_finder <- function(cluster_number, seurat_object, LFC.lim = 0.20, use_limit = FALSE, maximo_genes = 20, bootstrap_number = 100, bootstrap.sample_size = NULL, min.p.adj_value = 0,...){

  #restarting the seed to allow different results
  rm(.Random.seed, envir=globalenv())
  set.seed(Sys.time())

  #I have to calculate dynamically the number of clusters in the cluster_normalized that at the end is set to the Identity in the Seurat object.
  number_of_clusters <- length(unique(Idents(seurat_object)))

  #Calculate the sample size having a 15% of the total of clusters. For example 73*.15=10.95=10 or 29*.15 = 4.35=4
  if(is.null(bootstrap.sample_size)){
    bootstrap.sample_size <- as.integer(number_of_clusters*0.15)
  }

  final_result <- NULL
  for(counter in 1:bootstrap_number){

    set.seed(Sys.time())

    #TODO create a generalization for other datasets without the requeriment of having a cluster_normalized vector

    #Create a list with the number of clusters
    all_clusters <- c(1:number_of_clusters)
    all_clusters <- all_clusters[all_clusters!=cluster_number]

    other_clusters <- paste0('cluster_', sample(all_clusters, bootstrap.sample_size, replace = F))

    its_well_calculated <- TRUE

    result = tryCatch({
      markers.wilcox <- FindMarkers(seurat_object, ident.1 = paste0("cluster_", cluster_number), ident.2 = other_clusters , only.pos = TRUE, logfc.threshold = LFC.lim)
    }, warning = function(w) {
      print(paste0("warning", w))
    }, error = function(e) {
      print(paste0("error", e))
      its_well_calculated <- FALSE
    }, finally = {
    })

    if(its_well_calculated & !is.null(markers.wilcox)){
      #the p_val_adj could be 0 or less than 0.05 normally
      filtered.wilcox <- markers.wilcox[markers.wilcox$p_val_adj <= min.p.adj_value & !is.na(markers.wilcox$p_val_adj) & markers.wilcox$avg_logFC >= 0.0 & !is.na(rownames(markers.wilcox)) & rownames(markers.wilcox) != 'NA',]
      filtered.wilcox.genes <- rownames(markers.wilcox[markers.wilcox$p_val_adj < 0.05 & !is.na(rownames(markers.wilcox)) & rownames(markers.wilcox) != 'NA',])

      #I get all new results
      temporal_result <- rownames(filtered.wilcox)

      #have the top N of genes
      if(length(temporal_result) >= maximo_genes & use_limit){
        temporal_result <- temporal_result[1:maximo_genes]
      }else{
        temporal_result <- calculate_genes_using_dbscan(filtered.wilcox = filtered.wilcox)
      }

      final_result <- c(final_result, temporal_result)
      final_result <- unique(final_result)
    }
  }

  final_result
}

############################################
#' Function that selects the markers genes given a wilcox object by using dbscan algorithm and selecting just the genes that are considered like outliers (Without any cluster)
#' @description  Function that selects the markers genes given a wilcox object by using dbscan algorithm and selecting just the genes that are considered like outliers (Without any cluster)
#' @name calculate_genes_using_dbscan
#' @param filtered.wilcox Wilcox object with the Foldchange analysis base on the comparison between each cluster amoung the others
#' @return List with the marker genes
#' @export
calculate_genes_using_dbscan <- function(filtered.wilcox, plot.dbscan.results = FALSE){
  library("fpc")
  library("dbscan")
  final_result <- NULL

  print(paste0('markers with 0: ', nrow(markers.wilcox)))

  if(nrow(filtered.wilcox)>0){
    # Compute DBSCAN using fpc package
    db <- fpc::dbscan(filtered.wilcox[,c(2)], eps = 0.5, MinPts = 5)
    length(db$cluster[db$cluster==0])

    # Plot DBSCAN results
    if(plot.dbscan.results){
      plot(db, filtered.wilcox$avg_logFC, main = "DBSCAN", frame = FALSE)
    }

    #select genes in the cluster 0
    rownames(filtered.wilcox)[db$cluster==0]
    marker_genes <- NULL
    marker_genes.temp <- as.character(rownames(filtered.wilcox)[db$cluster==0])

    final_result <- c(final_result, marker_genes.temp)

    print(paste0('Markers after dbscan: ', length(final_result)))
  }else{
    print('Not markers with 0.0')
  }

  final_result
}

############################################
#' Function that calculates the foldchange between a given cluster over the rest
#' @description Function that calculates the foldchange between a given cluster over the rest (It doesn't use an external library to apply the wilcox algorithm). At the end it uses bonferroni correction
#' @name calculate_foldchange
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param LFC.lim a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
calculate_foldchange <- function(sc.eset, ct.varname, ct.group.temp, LFC.lim){
  gc()

  countmat <- exprs(sc.eset)

  #for checking progress for each cluster
  message1 <- paste0('Selecting markers (1 cluster vs others ): ',u, '->',unique(ct.group)[u])

  #Using custom_apply for sparse matrix
  group.1 <- custom_apply(X = countmat[,ct.group.temp],
                          MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                              pseudocount.use))
  group.2 <- custom_apply(X = countmat[,! ct.group.temp],
                          MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                              pseudocount.use))
  genes.diff <- rownames(sc.eset)[(group.1 - group.2) > LFC.lim]
  genes.diff.LFC_max <- max(group.1 - group.2)

  ###########################################
  print(paste0('LFC.lim: ', LFC.lim))
  genes.diff.LFC <- cbind(names=rownames(sc.eset), LFC=as.numeric((group.1 - group.2)))
  genes.diff.LFC <- genes.diff.LFC[genes.diff.LFC[,2] > LFC.lim,]
  genes.diff <- genes.diff.LFC[,1]

  count.use <- countmat[rownames(sc.eset) %in% genes.diff,]

  p_val <- sapply(1:nrow(count.use), function(x){
    wilcox.test(count.use[x,] ~ ct.group.temp)$p.value
  })

  p_val_adj <- p.adjust(p = p_val, method = "bonferroni",
                        n = nrow(count.use))

  genes_p_val_adj <- cbind(gene_name = rownames(count.use), p_val_adj=p_val_adj)


  ################################## Part for have the more significant genes but get those in the top base on the fold change

  #joining
  genes_LFC = data.frame(gene_name = genes.diff.LFC[,1], LFC = as.numeric(genes.diff.LFC[,2]))

  # data frame 2
  genes_p_value_adj = data.frame(gene_name = as.character(genes_p_val_adj[,1]), p_val_adj = as.numeric(genes_p_val_adj[,2]))
  genes_p_value_adv_LFC <-merge(x=genes_LFC,y=genes_p_value_adj,by="gene_name")

  genes_p_value_adv_LFC.filtered <- genes_p_value_adv_LFC[genes_p_value_adv_LFC$p_val_adj < 0.05, ]

  #order by LFC or P VAL ADJ
  genes_p_value_adv_LFC.filtered.ordered <- genes_p_value_adv_LFC.filtered[order(-genes_p_value_adv_LFC.filtered$LFC),]

  print(paste0('markers: ', length(markers.temp), ', new markers: ', length(genes_p_value_adv_LFC.filtered.ordered)))

  markers.temp <- as.character(genes_p_value_adv_LFC.filtered.ordered$gene_name)

  print(genes_p_value_adv_LFC.filtered.ordered)

  genes_p_value_adv_LFC.filtered.ordered
}

############################################
#' Tree-guided proportion estimation for ONE subject
#' @description Proportion estimation function for ONE-subject case, and apply tree-guided deconvolution
#' @name SCDC_prop_ONE_subcl_marker
#' @param bulk.eset ExpressionSet object for bulk samples
#' @param sc.eset ExpressionSet object for single cell samples
#' @param ct.varname variable name for 'cell types'
#' @param fl.varname variable name for first-level 'meta-clusters'
#' @param sample variable name for subject/samples
#' @param ct.sub a subset of cell types that are selected to construct basis matrix
#' @param ct.fl.sub 'cell types' for first-level 'meta-clusters'
#' @param iter.max the maximum number of iteration in WNNLS
#' @param nu a small constant to facilitate the calculation of variance
#' @param epsilon a small constant number used for convergence criteria
#' @param weight.basis logical, use basis matrix adjusted by MVW, default is T.
#' @param select.marker logical, select marker genes to perform deconvolution in tree-guided steps. Default is T.
#' @param markers A set of marker gene that input manully to be used in deconvolution. If NULL, then
#' @param marker.varname variable name of cluster groups when selecting marker genes. If NULL, then use ct.varname.
#' @param allgenes.fl logical, use all genes in the first-level deconvolution
#' @param pseudocount.use a constant number used when selecting marker genes, default is 1.
#' @param LFC.lim a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.
#' @param truep true cell-type proportions for bulk samples if known
#' @return Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
#' @export
SCDC_prop_ONE_subcl_marker <- function(bulk.eset, sc.eset, ct.varname, fl.varname, sample, truep = NULL,
                                       ct.sub = NULL, ct.fl.sub, iter.max = 3000, nu = 1e-04, epsilon = 0.001,
                                       weight.basis = F, bulk_disease = NULL, select.marker = T, markers = NULL, marker.varname = NULL,
                                       pseudocount.use = 1, LFC.lim = 0.5, allgenes.fl = F,
                                       ...)
{

  if (is.null(ct.sub)){
    ct.sub <- unique(sc.eset@phenoData@data[,ct.varname])[!is.na(unique(sc.eset@phenoData@data[,ct.varname]))]
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  ct.fl.sub <- ct.fl.sub[!is.na(ct.fl.sub)]
  bulk.eset <- bulk.eset[rowSums(exprs(bulk.eset))>0, , drop = FALSE]
  message('SCDC basis for main cluster...')
  sc.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.sub, ct.varname = ct.varname, sample = sample)
  message('SCDC basis for secondary cluster (metacluster)...')
  sc.fl.basis <- SCDC_basis_ONE(x = sc.eset, ct.sub = ct.fl.sub[!is.na(ct.fl.sub)],
                                ct.varname = fl.varname, sample = sample)
  if (select.marker){
    if (is.null(marker.varname)){
      marker.varname <- ct.varname
    }
    # wilcox test on two groups of cells for marker gene selection... (refer to seurat::FindMarkers)
    countmat <- exprs(sc.eset)
    ct.group <- sc.eset@phenoData@data[,marker.varname]
    markers.wilcox <- NULL

    global.min.LFC <- 0
    global.genes.diff.LFC_max.top.cluster <- 0
    #this should be parametric and corresponds to the minimun number of markers on each cluster
    top_genes <- 10

    #u=1
    for(u in 1:length(unique(ct.group))){
      ct.group.temp <- ct.group == unique(ct.group)[u]

      #for checking progress for each cluster
      message(paste0('Selecting markers (1 cluster vs others ): ',u, '->',unique(ct.group)[u]))

      #Using custom_apply for sparse matrix
      group.1 <- custom_apply(X = countmat[,ct.group.temp],
                       MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                           pseudocount.use))
      group.2 <- custom_apply(X = countmat[,! ct.group.temp],
                       MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
                                                           pseudocount.use))
      genes.diff <- rownames(sc.eset)[(group.1 - group.2) > LFC.lim]

      genes.diff.LFC_max <- max(group.1 - group.2)
      #We got the top N values (top_genes)
      genes.diff.LFC_max.top.cluster <- min(tail(sort(group.1 - group.2),top_genes))

      count.use <- countmat[rownames(sc.eset) %in% genes.diff,]

      #in order to check the global maximum value of LFC that it should be configured
      if(genes.diff.LFC_max < global.min.LFC | global.min.LFC==0){
        global.min.LFC <- genes.diff.LFC_max
      }

      #in order to check the global maximum value of LFC that it should be configured. With top N markers
      if(genes.diff.LFC_max.top.cluster < global.genes.diff.LFC_max.top.cluster | global.genes.diff.LFC_max.top.cluster==0){
        global.genes.diff.LFC_max.top.cluster <- genes.diff.LFC_max.top.cluster
      }

      #it's because there is just 1 gene
      if(class(count.use)=='numeric'){
        count.use <- t(count.use)
        rownames(count.use) <- c(genes.diff)
      }

      #check if there is none gene for the cluster
      if(nrow(count.use)==0){
        message(paste0("Zero Markers selected..."), 'Maximun LFC: ', genes.diff.LFC_max)
      }else{
        ##
        p_val <- sapply(1:nrow(count.use), function(x){
          wilcox.test(count.use[x,] ~ ct.group.temp)$p.value
        })
        p_val_adj <- p.adjust(p = p_val, method = "bonferroni",
                              n = nrow(count.use))
        markers.temp <- rownames(count.use)[p_val_adj < 0.05]
        markers.wilcox <- c(markers.wilcox, markers.temp)

        #For checking the genes selected as markers for an specific cluster
        message("Markers selected(",length(markers.temp), "): ", paste(shQuote(markers.temp), collapse=", "))
      }
    }
    markers <- unique(markers.wilcox)
    message("Global minimun LFC (1 marker): ", global.min.LFC)
    message("Global minimun LFC (Minimum ", top_genes, ' genes marker): ', global.genes.diff.LFC_max.top.cluster)
    message("Selected ",length(markers), " marker genes by Wilcoxon test...")
  } # else need input of marker genes for clustering

  # match genes / cells first
  if (weight.basis){
    basis <- sc.basis$basis.mvw
    basis.fl <- sc.fl.basis$basis.mvw
  } else {
    basis <- sc.basis$basis
    basis.fl <- sc.fl.basis$basis
  }
  if (!is.null(markers)){
    commongenes <- Reduce(intersect, list(rownames(basis), rownames(bulk.eset), markers))
    commongenes.fl <- Reduce(intersect, list(rownames(basis.fl), rownames(bulk.eset), markers))
    non_commongenes  <- Reduce(setdiff, commongenes, markers)
    non_commongenes.fl <- Reduce(setdiff, commongenes.fl, markers)
  } else {
    commongenes <- intersect(rownames(basis), rownames(bulk.eset))
    commongenes.fl <- intersect(rownames(basis.fl), rownames(bulk.eset))
    # stop when few common genes exist...
    if (length(commongenes) < 0.2 * min(dim(sc.eset)[1], dim(bulk.eset)[1])){
      stop('Too few common genes!')
    }
  }

  message(paste("Used", length(commongenes), "common genes for all cell types, \n",
                "Used", length(commongenes.fl), "common genes for first level cell types..."))

  message('Genes that are not shared between SC and bulk:')
  message('*Cluster(',length(non_commongenes), '):',paste(shQuote(non_commongenes), collapse=", "))
  message('**Metacluster(', length(non_commongenes.fl), '):', paste(shQuote(non_commongenes.fl), collapse=", "))

  basis.mvw <- basis[commongenes, ct.sub]
  basis.mvw.fl <- basis.fl[commongenes.fl, ct.fl.sub]

  xbulk0 <- getCPM0(exprs(bulk.eset)[commongenes,])
  xbulk <- as.matrix(xbulk0) ## whether to normalize all /common genes
  colnames(xbulk) <- colnames(bulk.eset)
  xbulk1 <- getCPM0(exprs(bulk.eset)[commongenes.fl,])
  xbulk.fl <- as.matrix(xbulk1)

  ALS.S <- sc.basis$sum.mat[ct.sub]
  N.bulk <- ncol(bulk.eset)
  valid.ct <- (colSums(is.na(basis.mvw)) == 0) & (!is.na(ALS.S))

  ALS.S.fl <- sc.fl.basis$sum.mat[ct.fl.sub]
  valid.ct.fl <- (colSums(is.na(basis.mvw.fl)) == 0) & (!is.na(ALS.S.fl))

  if (sum(valid.ct) <= 1) {
    stop("Not enough valid cell type!")
  }
  message(paste("Used", sum(valid.ct), "cell types in deconvolution...\n",
                "Used", sum(valid.ct.fl),"first level cell types ..."))

  basis.mvw <- basis.mvw[, valid.ct]
  ALS.S <- ALS.S[valid.ct]

  basis.mvw.fl <- basis.mvw.fl[, valid.ct.fl]
  ALS.S.fl <- ALS.S[valid.ct.fl]

  prop.est <- NULL
  rsquared <- NULL

  # prop estimation for each bulk sample:
  for (i in 1:N.bulk) {
    # i=1
    xbulk.temp <- xbulk[, i] *1e3 ## will affect a little bit
    message(paste(colnames(xbulk)[i], "has common genes", sum(xbulk[, i] != 0), "..."))
    if (allgenes.fl){
      markers.fl <- names(xbulk.temp)
    } else {
      markers.fl <- Reduce(intersect, list(markers, names(xbulk.temp)))
    }

    # first level NNLS:
    lm <- nnls::nnls(A=basis.mvw.fl[markers.fl,],b=xbulk.temp[markers.fl])
    delta <- lm$residuals
    wt.gene <- 1/(nu + delta^2)
    x.wt <- xbulk.temp[markers.fl] *sqrt(wt.gene)
    b.wt <- sweep(basis.mvw.fl[markers.fl,],1,sqrt(wt.gene),"*")

    lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
    prop.wt.fl <- lm.wt$x/sum(lm.wt$x)
    delta <- lm.wt$residuals

    for (iter in 1:iter.max){
      wt.gene <- 1/(nu + delta^2)
      x.wt <- xbulk.temp[markers.fl] * sqrt(wt.gene)
      b.wt <- sweep(basis.mvw.fl[markers.fl,],1,sqrt(wt.gene),"*")
      lm.wt <- nnls::nnls(A=b.wt, b=x.wt)
      delta.new <- lm.wt$residuals
      prop.wt.fl.new <- lm.wt$x/sum(lm.wt$x)

      if (sum(abs(prop.wt.fl.new - prop.wt.fl)) < epsilon){
        prop.wt.fl <- prop.wt.fl.new
        delta <- delta.new
        message("WNNLS for First level clusters Converged at iteration ", iter)
        break
      }
      prop.wt.fl <- prop.wt.fl.new
      delta <- delta.new
    }
    names(prop.wt.fl) <- colnames(basis.mvw.fl)

    # relationship between first level and overall
    rt <- table(sc.eset@phenoData@data[,ct.varname], sc.eset@phenoData@data[,fl.varname])
    rt <- rt[,ct.fl.sub]
    rt.list <- list()
    prop.wt <- NULL

    # prop.wt
    for (j in 1:ncol(rt)){ # for each first level cluster
      # j=1
      rt.list[[j]] <- rownames(rt)[rt[,j] >0]
      names(rt.list)[j] <- colnames(rt)[j]
      sub.cl <- rownames(rt)[rt[,j] >0]
      if (length(sub.cl) > 1 & prop.wt.fl[colnames(rt)[j]] > 0) {
        if (is.null(dim(prop.wt.fl))){
          # specify genes in xbulk.j??? first level genes?
          xbulk.j <- basis.mvw.fl[,j]*prop.wt.fl[j] + (xbulk.temp - basis.mvw.fl %*% lm.wt$x)*prop.wt.fl[j]
        } else {
          xbulk.j <- basis.mvw.fl[,j]*prop.wt.fl[,j] + (xbulk.temp - basis.mvw.fl %*% lm.wt$x)*prop.wt.fl[,j]
        }

        markers.sl <- Reduce(intersect, list(markers, rownames(xbulk.j)))

        ##############################################################################
        # make markers.sub as a list, for each of the first-level intra clusters.
        ##############################################################################

        basis.sl <- basis.mvw[markers.sl,rownames(rt)[rt[,j] >0]]
        lm.sl <- nnls::nnls(A=basis.sl,b=xbulk.j[markers.sl,])
        delta.sl <- lm.sl$residuals
        wt.gene.sl <- 1/(nu + delta.sl^2)
        x.wt.sl <- xbulk.j[markers.sl,]*sqrt(wt.gene.sl)
        b.wt.sl <- sweep(basis.sl,1,sqrt(wt.gene.sl),"*")

        lm.wt.sl <- nnls::nnls(A=b.wt.sl, b=x.wt.sl)
        prop.wt.sl <- lm.wt.sl$x/sum(lm.wt.sl$x)
        delta.sl <- lm.wt.sl$residuals

        for (iter in 1:iter.max){
          wt.gene.sl <- 1/(nu + delta.sl^2)
          x.wt.sl <- xbulk.j[markers.sl,] * sqrt(wt.gene.sl)
          b.wt.sl <- sweep(basis.sl,1,sqrt(wt.gene.sl),"*")
          lm.wt.sl <- nnls::nnls(A=b.wt.sl, b=x.wt.sl)
          delta.sl.new <- lm.wt.sl$residuals
          prop.wt.sl.new <- lm.wt.sl$x/sum(lm.wt.sl$x)

          if (sum(abs(prop.wt.sl.new - prop.wt.sl)) < epsilon){
            prop.wt.sl <- prop.wt.sl.new
            delta.sl <- delta.sl.new
            message("WNNLS for Second level clusters",paste(shQuote(rownames(rt)[rt[,j] >0]), collapse=", ")," Converged at iteration ", iter)
            break
          }
          prop.wt.sl <- prop.wt.sl.new
          delta.sl <- delta.sl.new
        }
        names(prop.wt.sl) <- sub.cl
        prop.wt <- c(prop.wt, prop.wt.sl*prop.wt.fl[colnames(rt)[j]])
      } else if (length(sub.cl) == 1){
        # j=2
        prop.wt <- c(prop.wt, prop.wt.fl[colnames(rt)[j]])
      } else if (length(sub.cl) > 1 & prop.wt.fl[colnames(rt)[j]] == 0){
        prop.wt.sl <- rep(0, length(sub.cl))
        names(prop.wt.sl) <- sub.cl
        prop.wt <- c(prop.wt, prop.wt.sl)
      }

    }
    prop.est <- rbind(prop.est, prop.wt)
  }
  rownames(prop.est) <- colnames(bulk.eset)

  peval <- NULL
  if (!is.null(truep)){
    peval <- SCDC_eval(ptrue= truep, pest = prop.est, pest.names = c('SCDC'),
                       dtname = 'Perou', select.ct = ct.sub, bulk_obj = bulk.eset,
                       bulk_disease = bulk_disease)
  }

  return(list(prop.est = prop.est, prop.wt.fl = prop.wt.fl, basis.mvw = basis.mvw, peval = peval,
              sc.basis = sc.basis, sc.fl.basis = sc.fl.basis))
}
