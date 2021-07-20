#' Visualize a pedigree with cancer affection statuses and other information
#'
#' This function relies on the kinship2 package
#'
#' @param ped a pedigree that contains \code{ID, Sex, MotherID, FatherID,
#' isProband, CurAge}
#' It is assumed that in this data.frame, Sex = 1 if male, Sex = 0 if female
#' and Sex = NA if unkonwn
#' Affliction is coded in "isAffXX" columns, where XX are short names of cancers
#' The age of diagnosis is correspondingly coded in "AgeXX" columns
#' @param annot.cancers the cancers shortnames to display. When set to
#' default /code{NULL}, the top four cancers will be displayed.
#' @param annot.feature ONE feature that we would get annotation from pedigree,
#'  one of the choices from \code{c("Ancestry","Twins","CurAge","race")}
#' @param title a string for the title on the plot
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics arrows frame legend lines par points polygon segments strheight strwidth text title
#' @examples
#' visPed(test_fam_1)
#' @export
visPed <- function(ped, annot.cancers = "all", annot.feature = "CurAge",
                   title = "Your Pedigree") {
  # Translate columns of pedigree into kinship2 standard
  # In kinship2, Sex = 1 when male, Sex = 2 when female and Sex = 3 when unknown
  # We assume that ped has -999 for missing MotherID or FatherID
  ks_df <- with(ped, {
    Sex <- ifelse(Sex == 0, 2, Sex)
    Sex <- ifelse(is.na(Sex), 3, Sex)

    MotherID <- ifelse(MotherID == -999, NA, MotherID)
    FatherID <- ifelse(FatherID == -999, NA, FatherID)

    cbind(ID, FatherID, MotherID, Sex, isDead)
  })

  # Translate column names
  colnames(ks_df) <- c("id", "dadid", "momid", "sex", "censor")

  # Create family object(s) in kinship2
  fam_vec <- kinship2::makefamid(
    id = ks_df[, "id"],
    father.id = ks_df[, "dadid"],
    mother.id = ks_df[, "momid"]
  )

  # Exit if there are disconnected families
  if (length(unique(fam_vec)) > 1) {
    rlang::abort(paste0("Pedigree contains disconnected families, ",
                 "please first double check with kinship2::makefamid"),
      level = "DisconnectedFamily"
    )
  }

  # Currently supported cancers: BC OC BRA COL ENDO GAS KID MELA PANC PROS SI (SMA renamed)
  supported_cancers <- c("BRA", "BC", "CER", "COL", "ENDO",
                         "GAS", "KID", "LEUK", "MELA", "OC",
                         "OST", "PANC", "PROS", "SI",
                         "STS", "THY", "UB", "HEP")

  # Get the cancers to plot
  cancers_in_pedigree <- substring(colnames(ped[, grepl("isAff", names(ped)),
                                                drop = FALSE]), 6)

  if (annot.cancers == "all") {
    cancer_status <- ped[, paste0("isAff", cancers_in_pedigree), drop = FALSE]
    non_appearing_cancers <- paste0("isAff", setdiff(supported_cancers,
                                                     cancers_in_pedigree))
    if (length(non_appearing_cancers) == 1 && non_appearing_cancers == "isAff") {
      # if all cancers present, don't need to do anything
    } else {
      # when there are some missing cancers, set them all to 0
      cancer_status[non_appearing_cancers] <- 0
    }
    cancers_to_plot <- cancers_in_pedigree
  }

  # Reorder cancer columns so that they are always consistent across plots
  cancer_status <- cancer_status[, order(colnames(cancer_status))]

  # Set up the pedigree object in kinship2
  ks_ped <- kinship2::pedigree(
    id = ks_df[, "id"],
    dadid = ks_df[, "dadid"],
    momid = ks_df[, "momid"],
    sex = ks_df[, "sex"],
    affected = as.matrix(cancer_status),
    status = ks_df[, "censor"]
  )


  # Set up annotations
  annotations <- NULL
  if (!is.null(annot.cancers)) {
    cancer_ages <- ped[, paste0("Age", cancers_to_plot), drop = FALSE]

    age_list <- list()
    for (cancer in cancers_to_plot) {
      # Get the affection status age - it could be NA
      # Avoid writing anything if there is an NA
      affection_of_cancer <- ped[[paste0("isAff", cancer)]]
      age_of_cancer <- ped[[paste0("Age", cancer)]]
      # Remove incorrectly placed ages - make sure that person is affected
      ages_to_annotate <- ifelse(affection_of_cancer,
        age_of_cancer,
        NA
      )
      # If they are affected, but missing age, put "NA" as the age
      ages_to_annotate[is.na(ages_to_annotate) & affection_of_cancer] <- "NA"
      ages_to_annotate[is.na(ages_to_annotate)] <- ""
      annotation_string <- ifelse(ages_to_annotate == "",
        "",
        paste(cancer, ages_to_annotate, sep = ": ")
      )
      age_list[[cancer]] <- annotation_string
    }

    annotations <- do.call(paste, c(age_list, list(sep = "\n")))

    # Define a recursive function which cuts of excess \ns from the end
    cut_n <- function(string) {
      if (regmatches(string, regexpr(".{2}$", string)) == "\n\n") {
        output <- substr(string, 1, nchar(string) - 1)
        return(cut_n(output))
      } else {
        return(string)
      }
    }

    # Loop over the annotations to fix them
    for (i in 1:length(annotations)) {
      if (annotations[i] == "\n" ||
          all(unlist(strsplit(annotations[i], split = "\n")) == "")) {
        # If the entire string is \n or \n repeated, get rid of it
        annotations[i] <- ""
      } else if (unlist(strsplit(annotations[i], split = "\n"))[1] == "") {
        # If there is a \n at the start, get rid of it
        annotations[i] <- substring(annotations[i], 2)
      }
      if (annotations[i] != "" && annotations[i] != "\n") {
        # Put a \n at the end if there is actually a notation
        # But don't if there is already one
        annotations[i] <- cut_n(paste0(annotations[i], "\n"))
      }

      # There might still be "\n" at the start... remove them
      # We might need to do this recursively...
      while (!identical(strsplit(annotations[i], split = "\n")[[1]], character(0)) &
        substr(annotations[i], 1, 1) == "\n") {
        annotations[i] <- substring(annotations[i], 2)
      }

      # If there are multiple \n between cancers, make sure there is only one
      double_ns <- strsplit(annotations[i], split = "\n")[[1]] == ""
      if (!identical(which(double_ns), integer(0))) {
        # Get rid of the indexes
        annotations[i] <- paste(paste(strsplit(annotations[i], split = "\n")[[1]][!double_ns],
          collapse = "\n"
        ),
        "\n",
        collapse = ""
        )
      }
    }
  }

  if (!is.null(annot.feature)) {
    annot.feature <- intersect(
      c("Twins", "Ancestry", "CurAge", "race", "ID"),
      annot.feature
    )
    features_cols <- ped[c(annot.feature, "isDead")]
    feature_annotations <- lapply(1:nrow(ped), function(i) {
      mark <- substr(features_cols[i, 1], 1, 3)
      if (annot.feature == "CurAge") {
        annot.feature <- "Age"
      }
      annot.feature <- substr(annot.feature, 1, 3)
      # If the mark is na, write that
      mark[is.na(mark)] <- "NA"
      # mark <- mark[!is.na(mark)]
      # Do not plot Age if isDead==1
      # if (!(annot.feature == "Age" && features_cols[i, "isDead"] == 1)) {
      #   paste(annot.feature, mark, sep = ": ")
      # }
      paste(annot.feature, mark, sep = ": ")
    })
    annotations <- paste0(annotations, sapply(feature_annotations,
      paste0,
      collapse = "\n"
    ))
  }

  # Identify proband
  probands <- which(ped$isProband == 1)

  # Plot using kinship2 engine
  visEngine(ks_ped,
    annot = annotations,
    feature.name = annot.feature,
    which.proband = probands,
    main_title = title
  )

  invisible(ks_ped)
}

#' Visualization engine inspired by kinship2
#'
#' @param x a pedigree data.frame
#' @param annot the annotation string
#' @param feature.name the names of the features
#' @param which.proband which ID is the proband
#' @param add.ids which additional IDs to add
#' @param id ID column from x pedigree
#' @param status alive or dead status
#' @param affected affected or unaffected column
#' @param cex spacing
#' @param col column number, either 1 or the number of individuals
#' @param symbolsize symbol size
#' @param branch branch size
#' @param packed boolean to pack the pedigree plot
#' @param align align parameters
#' @param width width parameter
#' @param density density parameters
#' @param main_title title on plot
#' @param mar legend margins
#' @param angle angle of labels
#' @param keep.par boolean to keep or discard old graphics parameters
#' @param subregion boolean whether or not to split plot into subregions
#' @param pconnect controls parent connections
#' @import kinship2
visEngine <- function(x, annot, feature.name = NULL,
                      which.proband = NULL, add.ids = x$id[x$add == 1],
                      id = x$id, status = x$status, affected = x$affected,
                      cex = 1, col = 1, symbolsize = 1, branch = 0.6, packed = TRUE,
                      align = c(1.5, 2), width = 10, density = c(-1, 35, 65, 20),
                      main_title = "Your Pedigree",
                      mar = c(4.1, 1, 4.1, 1), angle = c(90, 65, 40, 0), keep.par = FALSE,
                      subregion, pconnect = 0.5) {
  getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

  Call <- match.call()
  affected[is.na(affected)] <- -1

  colmat <- getPalette(ncol(affected))
  density <- rep(-1, ncol(affected))
  angle <- rep(90, ncol(affected))

  n <- length(x$id)
  if (is.null(status)) {
    status <- rep(0, n)
  } else {
    if (!all(status == 0 | status == 1)) {
      stop("Invalid status code")
    }
    if (length(status) != n) {
      stop("Wrong length for status")
    }
  }
  if (!missing(id)) {
    if (length(id) != n) {
      stop("Wrong length for id")
    }
  }
  if (is.null(affected)) {
    affected <- matrix(0, nrow = n)
  }
  else {
    if (is.matrix(affected)) {
      if (nrow(affected) != n) {
        stop("Wrong number of rows in affected")
      }
      if (is.logical(affected)) {
        affected <- 1 * affected
      }
      # the maximum number of cancers a individual can have
      affect_indicator <- affected != 0
      maxncancer <- max(apply(affect_indicator, 1, sum))
      # if (maxncancer > length(angle) || maxncancer > length(density))
      #     stop("More columns in the affected matrix than angle/density values")
    }
    else {
      if (length(affected) != n) {
        stop("Wrong length for affected")
      }
      if (is.logical(affected)) {
        affected <- as.numeric(affected)
      }
      if (is.factor(affected)) {
        affected <- as.numeric(affected) - 1
      }
    }
    # if (max(affected, na.rm = TRUE) > min(affected, na.rm = TRUE)) {
    #     affected <- matrix(affected - min(affected, na.rm = TRUE),
    #                        nrow = n)
    # }
    # else {
    #     affected <- matrix(affected, nrow = n)
    # }
    # if (!all(affected == 0 | affected == 1 | affected == -1))
    #     print(affected)
    #     print("adfasd")
    #     stop("Invalid code for affected status")
  }
  if (length(col) == 1) {
    col <- rep(col, n)
  } else if (length(col) != n) {
    stop("Col argument must be of length 1 or n")
  }

  # sub-region -------

  # plist

  # subreg

  subregion2 <- function(plist, subreg) {
    if (subreg[3] < 1 || subreg[4] > length(plist$n)) {
      stop("Invalid depth indices in subreg")
    }
    lkeep <- subreg[3]:subreg[4]
    for (i in lkeep) {
      if (!any(plist$pos[i, ] >= subreg[1] & plist$pos[i, ] <= subreg[2])) {
        stop(paste("No subjects retained on level", i))
      }
    }
    nid2 <- plist$nid[lkeep, ]
    n2 <- plist$n[lkeep]
    pos2 <- plist$pos[lkeep, ]
    spouse2 <- plist$spouse[lkeep, ]
    fam2 <- plist$fam[lkeep, ]
    if (!is.null(plist$twins)) {
      twin2 <- plist$twins[lkeep, ]
    }
    for (i in 1:nrow(nid2)) {
      keep <- which(pos2[i, ] >= subreg[1] & pos2[i, ] <=
        subreg[2])
      nkeep <- length(keep)
      n2[i] <- nkeep
      nid2[i, 1:nkeep] <- nid2[i, keep]
      pos2[i, 1:nkeep] <- pos2[i, keep]
      spouse2[i, 1:nkeep] <- spouse2[i, keep]
      fam2[i, 1:nkeep] <- fam2[i, keep]
      if (!is.null(plist$twins)) {
        twin2[i, 1:nkeep] <- twin2[i, keep]
      }
      if (i < nrow(nid2)) {
        tfam <- match(fam2[i + 1, ], keep, nomatch = 0)
        fam2[i + 1, ] <- tfam
        if (any(spouse2[i, tfam] == 0)) {
          stop("A subregion cannot separate parents")
        }
      }
    }
    n <- max(n2)
    out <- list(
      n = n2[1:n], nid = nid2[, 1:n, drop = F],
      pos = pos2[, 1:n, drop = F], spouse = spouse2[, 1:n, drop = F], fam = fam2[, 1:n, drop = F]
    )
    if (!is.null(plist$twins)) {
      out$twins <- twin2[, 1:n, drop = F]
    }
    out
  }

  plist <- kinship2::align.pedigree(x, packed = packed, width = width, align = align)

  op <- options(repr.plot.height = nrow(plist$pos) * 10)
  # on.exit(options(op))
  if (!missing(subregion)) {
    plist <- subregion2(plist, subregion)
  }
  xrange <- range(plist$pos[plist$nid > 0])
  maxlev <- nrow(plist$pos)

  sinf <- dim(plist$nid)

  ### size ###

  # dev.new(width = 0.2*sinf[2], height = 20*sinf[1])
  # if (dev.cur() <= 3) dev.off(dev.cur())
  # dev.new(height=30, units = "units")
  frame()

  # allow additional space for cancer legend
  oldpar <- par(xpd = T, mar = mar + c(0, 0, 4, 0))
  # oldpar <- par(mar = mar, xpd = TRUE)
  psize <- par("pin")
  stemp1 <- strwidth("ABC", units = "inches", cex = cex) * 2.5 / 3
  stemp2 <- strheight("1g", units = "inches", cex = cex)
  stemp3 <- max(strheight(id, units = "inches", cex = cex))
  ht1 <- psize[2] / maxlev - (stemp3 + 1.5 * stemp2)
  if (ht1 <= 0) {
    stop("Labels leave no room for the graph, reduce cex")
  }
  ht2 <- psize[2] / (maxlev + (maxlev - 1) / 2)
  wd2 <- 0.8 * psize[1] / (0.8 + diff(xrange))
  boxsize <- symbolsize * min(ht1, ht2, stemp1, wd2)
  hscale <- (psize[1] - boxsize) / diff(xrange)
  vscale <- (psize[2] - (stemp3 + stemp2 / 2 + boxsize)) / max(1, maxlev - 1)
  boxw <- boxsize / hscale
  boxh <- boxsize / vscale
  labh <- stemp2 / vscale
  legh <- min(1 / 4, boxh * 1.5)
  par(usr = c(
    xrange[1] - boxw / 2, xrange[2] + boxw / 2,
    maxlev + boxh + stemp3 + stemp2 / 2, 1
  ))
  circfun <- function(nslice, n = 50) {
    nseg <- ceiling(n / nslice)
    theta <- -pi / 2 - seq(0, 2 * pi, length = nslice + 1)
    out <- vector("list", nslice)
    for (i in 1:nslice) {
      theta2 <- seq(theta[i], theta[i + 1], length = nseg)
      out[[i]] <- list(
        x = c(0, cos(theta2) / 2),
        y = c(0, sin(theta2) / 2) + 0.5
      )
    }
    out
  }
  polyfun <- function(nslice, object) {
    zmat <- matrix(0, ncol = 4, nrow = length(object$x))
    zmat[, 1] <- object$x
    zmat[, 2] <- c(object$x[-1], object$x[1]) - object$x
    zmat[, 3] <- object$y
    zmat[, 4] <- c(object$y[-1], object$y[1]) - object$y
    ns1 <- nslice + 1
    theta <- -pi / 2 - seq(0, 2 * pi, length = ns1)
    x <- y <- double(ns1)
    for (i in 1:ns1) {
      z <- (tan(theta[i]) * zmat[, 1] - zmat[, 3]) / (zmat[, 4] - tan(theta[i]) * zmat[, 2])
      tx <- zmat[, 1] + z * zmat[, 2]
      ty <- zmat[, 3] + z * zmat[, 4]
      inner <- tx * cos(theta[i]) + ty * sin(theta[i])
      indx <- which(is.finite(z) & z >= 0 & z <= 1 & inner > 0)
      x[i] <- tx[indx]
      y[i] <- ty[indx]
    }
    nvertex <- length(object$x)
    temp <- data.frame(
      indx = c(1:ns1, rep(0, nvertex)),
      theta = c(theta, object$theta), x = c(x, object$x),
      y = c(y, object$y)
    )
    temp <- temp[order(-temp$theta), ]
    out <- vector("list", nslice)
    for (i in 1:nslice) {
      rows <- which(temp$indx == i):which(temp$indx == (i + 1))
      out[[i]] <- list(x = c(0, temp$x[rows]), y = c(0, temp$y[rows]) + 0.5)
    }
    out
  }
  if (ncol(affected) == 1) {
    polylist <- list(
      square = list(list(x = c(-1, -1, 1, 1) / 2, y = c(0, 1, 1, 0))),
      circle = list(list(x = 0.5 * cos(seq(0, 2 * pi, length = 50)), y = 0.5 * sin(seq(0, 2 * pi, length = 50)) + 0.5)),
      diamond = list(list(x = c(0, -0.5, 0, 0.5), y = c(0, 0.5, 1, 0.5))),
      triangle = list(list(x = c(0, -0.56, 0.56), y = c(0, 1, 1)))
    )
  }
  else {
    nc <- maxncancer
    square <- polyfun(nc, list(
      x = c(-0.5, -0.5, 0.5, 0.5),
      y = c(-0.5, 0.5, 0.5, -0.5), theta = -c(3, 5, 7, 9) * pi / 4
    ))
    circle <- circfun(nc)
    diamond <- polyfun(nc, list(x = c(0, -0.5, 0, 0.5), y = c(-0.5, 0, 0.5, 0), theta = -(1:4) * pi / 2))
    triangle <- polyfun(nc, list(x = c(-0.56, 0, 0.56), y = c(-0.5, 0.5, -0.5), theta = c(-2, -4, -6) * pi / 3))
    polylist <- list(square = square, circle = circle, diamond = diamond, triangle = triangle)
  }


  # drawbox ----- draw one person ---------
  drawbox <- function(x, y, sex, affected, status, colmat, polylist,
                      density, angle, boxw, boxh) {
    total_fractions <- 1:maxncancer
    aff_cids <- which(affected == 1)
    unknown_cids <- which(affected == -1)


    if (length(aff_cids) != 0) {
      unknown_fractions <- total_fractions[-(1:length(aff_cids))]
      for (i in seq_along(aff_cids)) {
        polygon(x + (polylist[[sex]])[[i]]$x * boxw,
          y + (polylist[[sex]])[[i]]$y * boxh,
          col = colmat[aff_cids[i]],
          border = 1, density = density[i], angle = angle[i]
        )
      }
    }
    else {
      unknown_fractions <- total_fractions
    }

    if (length(unknown_cids) != 0) {
      for (i in seq_along(unknown_cids)) {
        unf <- unknown_fractions[i]
        lab <- unknown_cids[i]
        # polygon(x + (polylist[[sex]])[[unf]]$x * boxw,
        #         y + (polylist[[sex]])[[unf]]$y * boxh, col = NA,
        #         border = col)
        midx <- x + mean(range(polylist[[sex]][[unf]]$x * boxw))
        midy <- y + mean(range(polylist[[sex]][[unf]]$y * boxh))
        # print(c(midx, midy))
        points(midx, midy,
          pch = "?",
          cex = min(1, cex * 2 / length(affected)),
          col = colmat[lab]
        )
      }
    }

    for (i in 1:maxncancer) {
      polygon(x + (polylist[[sex]])[[i]]$x * boxw,
        y + (polylist[[sex]])[[i]]$y * boxh,
        col = NA,
        border = 1
      )
    }
    if (status == 1) {
      segments(x - 0.6 * boxw, y + 1.1 * boxh, x + 0.6 *
        boxw, y - 0.1 * boxh)
    }
  }

  sex <- ifelse(x$sex == "female", 2, 1)
  for (i in 1:maxlev) {
    for (j in seq_len(plist$n[i])) {
      k <- plist$nid[i, j]
      drawbox(
        plist$pos[i, j], i, sex[k], affected[k, ],
        status[k], colmat, polylist, density, angle,
        boxw, boxh
      )
      # Omit the ID number for simplicity
      text_target <- paste(annot[k], sep = "\n")
      # text_target <- paste(id[k], annot[k], sep = "\n")
      tts <- strsplit(text_target, split = "\n")[[1]]
      for (tti in seq_along(tts)) {
        tt <- tts[tti]
        text(plist$pos[i, j], i + 0.8 * boxh + labh * 1.2 * tti, tt,
          # cex = ifelse(tti == 1, cex, cex*0.7), adj = c(0.5, 1))
          cex = ifelse(tti == 1, cex * 0.7, cex * 0.7), adj = c(0.5, 1)
        )
      }

      if (k %in% which.proband) {
        x0 <- plist$pos[i, j] - 0.3 * boxw - 0.2
        y0 <- i + 0.8 * boxh + labh
        x1 <- x0 + 0.1
        y1 <- y0 - 0.1
        length <- 0.05
        arrows(x0, y0, x1, y1, length = length)
      }
    }
  }

  maxcol <- ncol(plist$nid)
  for (i in 1:maxlev) {
    tempy <- i + boxh / 2
    if (any(plist$spouse[i, ] > 0)) {
      temp <- (1:maxcol)[plist$spouse[i, ] > 0]
      segments(
        plist$pos[i, temp] + boxw / 2, rep(
          tempy,
          length(temp)
        ), plist$pos[i, temp + 1] - boxw / 2,
        rep(tempy, length(temp))
      )
      temp <- (1:maxcol)[plist$spouse[i, ] == 2]
      if (length(temp)) {
        tempy <- tempy + boxh / 10
        segments(
          plist$pos[i, temp] + boxw / 2, rep(
            tempy,
            length(temp)
          ), plist$pos[i, temp + 1] - boxw / 2,
          rep(tempy, length(temp))
        )
      }
    }
  }
  for (i in 2:maxlev) {
    zed <- unique(plist$fam[i, ])
    zed <- zed[zed > 0]
    for (fam in zed) {
      xx <- plist$pos[i - 1, fam + 0:1]
      parentx <- mean(xx)
      who <- (plist$fam[i, ] == fam)
      if (is.null(plist$twins)) {
        target <- plist$pos[i, who]
      } else {
        twin.to.left <- (c(0, plist$twins[i, who])[1:sum(who)])
        temp <- cumsum(twin.to.left == 0)
        tcount <- table(temp)
        target <- rep(tapply(
          plist$pos[i, who], temp,
          mean
        ), tcount)
      }
      yy <- rep(i, sum(who))
      segments(plist$pos[i, who], yy, target, yy - legh)
      if (any(plist$twins[i, who] == 1)) {
        who2 <- which(plist$twins[i, who] == 1)
        temp1 <- (plist$pos[i, who][who2] + target[who2]) / 2
        temp2 <- (plist$pos[i, who][who2 + 1] + target[who2]) / 2
        yy <- rep(i, length(who2)) - legh / 2
        segments(temp1, yy, temp2, yy)
      }
      if (any(plist$twins[i, who] == 3)) {
        who2 <- which(plist$twins[i, who] == 3)
        temp1 <- (plist$pos[i, who][who2] + target[who2]) / 2
        temp2 <- (plist$pos[i, who][who2 + 1] + target[who2]) / 2
        yy <- rep(i, length(who2)) - legh / 2
        text((temp1 + temp2) / 2, yy, "?")
      }
      segments(min(target), i - legh, max(target), i - legh)
      if (diff(range(target)) < 2 * pconnect) {
        x1 <- mean(range(target))
      } else {
        x1 <- pmax(
          min(target) + pconnect,
          pmin(max(target) - pconnect, parentx)
        )
      }
      y1 <- i - legh
      if (branch == 0) {
        segments(x1, y1, parentx, (i - 1) + boxh / 2)
      } else {
        y2 <- (i - 1) + boxh / 2
        x2 <- parentx
        ydelta <- ((y2 - y1) * branch) / 2
        segments(
          c(x1, x1, x2), c(y1, y1 + ydelta, y2 - ydelta),
          c(x1, x2, x2), c(y1 + ydelta, y2 - ydelta, y2)
        )
      }
    }
  }
  arcconnect <- function(x, y) {
    xx <- seq(x[1], x[2], length = 15)
    yy <- seq(y[1], y[2], length = 15) + (seq(-7, 7))^2 / 98 -
      0.5
    lines(xx, yy, lty = 2)
  }
  uid <- unique(plist$nid)
  for (id in unique(uid[uid > 0])) {
    indx <- which(plist$nid == id)
    if (length(indx) > 1) {
      tx <- plist$pos[indx]
      ty <- ((row(plist$pos))[indx])[order(tx)]
      tx <- sort(tx)
      for (j in 1:(length(indx) - 1)) {
        arcconnect(tx[j +
          0:1], ty[j + 0:1])
      }
    }
  }
  ckall <- x$id[is.na(match(x$id, x$id[plist$nid]))]

  add_legend <- function(...) {
    opar <- par(
      fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
      mar = c(0, 2, 0, 0), new = TRUE
    )
    on.exit(par(opar))
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend(...)
  }

  cancer_legend <- add_legend("topleft",
    legend = substring(dimnames(x$affected)[[2]], 6),
    border = NA,
    col = colmat,
    horiz = TRUE,
    cex = 0.8,
    bty = "n",
    # pch = c(rep(".", length(colmat)), "[C]"),
    fill = colmat
  )

  if (!is.null(feature.name)) {
    fl_x <- cancer_legend$rect$left + 0.05 * boxh
    fl_y <- cancer_legend$rect$top - 0.7 * boxh #- cancer_legend$rect$h
    if (feature.name != "CurAge") {
      fl <- add_legend(
        x = fl_x - 0.09, y = fl_y,
        legend = paste0(substr(feature.name, 1, 3), "  ", feature.name),
        col = 1,
        border = NULL,
        horiz = TRUE,
        # title = "Feature Annotation",
        # pch = substr(feature.name, 1, 3),
        pt.cex = 0.8,
        # fill = NA,
        bty = "n",
        cex = 0.8
      )
    }
  }
  if (!is.null(which.proband)) {
    pl_x <- cancer_legend$rect$left
    pl_y <- cancer_legend$rect$top - boxh
    pl <- add_legend(
      x = pl_x*0.97, y = pl_y,
      legend = "Proband",
      cex = 0.8,
      bty = "n"
    )
    arrows(pl_x + 0.02, pl_y - 0.08, pl_x + 0.05, pl_y - 0.035, length = 0.03)
  }

  title(main_title)
  if (length(ckall > 0)) {
    cat("Did not plot the following people:", ckall, "\n")
  }
  if (!keep.par) {
    par(oldpar)
  }
  tmp <- match(1:length(x$id), plist$nid)
  invisible(list(
    plist = plist, x = plist$pos[tmp], y = row(plist$pos)[tmp],
    boxw = boxw, boxh = boxh, call = Call
  ))
}
