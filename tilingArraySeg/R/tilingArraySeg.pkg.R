DATE <-
"Fri Nov 15 22:11:49 2013"
VERSION <-
"0.0.5"
.onLoad <-
function( libname, pkgname ) { ##.onAttach
    cat( "Loading ", pkgname, " version ", VERSION, " (", DATE, ")\n", sep="" )
    cat( "Copyright (C) David J Reiss, Institute for Systems Biology.\n" )
    cat( "Please email dreiss.isb@gmail.com.org if you run into any issues.\n" )
  }
`%betw%` <-
function (a, b) 
if (b[2] >= b[1]) a >= b[1] & a <= b[2] else a >= b[2] & a <= 
    b[1]
classify.start.vs.stop.sites <-
function (brk.sites, brk.deltas, return.data = F, use = c("refs", 
    "rats", "pred", "cor"), near.gene = NA, window = NA, nearest = 1, 
    do.average = T) 
{
    gene.coords <- get.gene.coords()
    gene.coo <- as.data.frame(gene.coords)
    deltas <- breaks <- chrs <- dirs <- NULL
    for (chr in names(brk.sites)) {
        for (dir in names(brk.sites[[1]])) {
            deltas <- rbind(deltas, brk.deltas[[chr]][[dir]])
            breaks <- c(breaks, brk.sites[[chr]][[dir]])
            chrs <- c(chrs, rep(chr, nrow(brk.deltas[[chr]][[dir]])))
            dirs <- c(dirs, rep(dir, nrow(brk.deltas[[chr]][[dir]])))
        }
    }
    deltas[is.na(deltas)] <- 0
    rat.cols <- grep("rats", colnames(deltas))
    if (length(rat.cols) > 0) 
        deltas[, rat.cols] <- abs(deltas[, rat.cols])
    chip.cols <- grep("chip", colnames(deltas))
    if (length(chip.cols) > 0) 
        deltas[, chip.cols] <- abs(deltas[, chip.cols])
    deltas <- deltas[, grepl(paste(use, collapse = "|"), colnames(deltas)), 
        drop = F]
    halo.starts <- as.integer(as.character(gene.coo$Start.new))
    halo.starts[!is.na(gene.coo$Start.new)] <- as.integer(as.character(gene.coo$Start.new))[!is.na(gene.coo$Start.new)]
    halo.genes <- as.character(gene.coo$canonical_Name)
    halo.stops <- as.integer(as.character(gene.coo$Stop.new))
    halo.stops[!is.na(gene.coo$Stop.new)] <- as.integer(as.character(gene.coo$Stop.new))[!is.na(gene.coo$Stop.new)]
    ddirs <- ifelse(dirs == "FORWARD", "Rev", "For")
    if (!exists("probs.expressed")) 
        probs.expressed <- get.orf.prob.expressed()
    probs.expressed[probs.expressed == -Inf] <- -20
    gene.coo[!is.na(gene.coo$Start.new), "Start"] <- as.integer(as.character(gene.coo[!is.na(gene.coo$Start.new), 
        "Start.new"]))
    gene.coo[!is.na(gene.coo$Stop.new), "Stop"] <- as.integer(as.character(gene.coo[!is.na(gene.coo$Stop.new), 
        "Stop.new"]))
    if (!require(multicore, quietly = T, warn.conflicts = F)) 
        mclapply <- lapply
    tmp <- mclapply(1:length(breaks), function(i) {
        br <- breaks[i]
        chr <- chrs[i]
        dir <- dirs[i]
        inds <- gene.coo$where == chr & gene.coo$Orientation == 
            ddirs[i]
        start.dists <- abs(br - halo.starts[inds])
        start.dist <- min(start.dists, na.rm = T)
        gene.start <- halo.genes[inds][which(start.dists == start.dist)][1]
        start.prob.exp <- probs.expressed[gene.start]
        if (is.na(start.prob.exp)) 
            start.prob.exp <- 0
        stop.dists <- abs(br - halo.stops[inds])
        stop.dist <- min(stop.dists, na.rm = T)
        gene.stop <- halo.genes[inds][which(stop.dists == stop.dist)][1]
        stop.prob.exp <- probs.expressed[gene.stop]
        if (is.na(stop.prob.exp)) 
            stop.prob.exp <- 0
        is.start <- start.dist < stop.dist
        start.stop.prob.exp <- if (is.start) 
            start.prob.exp
        else stop.prob.exp
        names(start.stop.prob.exp) <- NULL
        gene.hit <- if (is.start) 
            gene.start
        else gene.stop
        if (i%%200 == 0) 
            cat(chr, dir, br, is.start, gene.hit, min(c(start.dist, 
                stop.dist)), start.stop.prob.exp, "\n")
        list(c(brk = br, dist = min(c(start.dist, stop.dist)), 
            is.start = as.numeric(is.start), prob.exp = start.stop.prob.exp, 
            deltas[i, ]), c(chr = chr, dir = dir))
    })
    fit.data <- fit.dat <- list()
    for (i in 1:length(tmp)) {
        fit.dat[[i]] <- tmp[[i]][[1]]
        fit.data[[i]] <- tmp[[i]][[2]]
    }
    rm(tmp)
    all.data <- cbind(data.frame(do.call(rbind, fit.data)), data.frame(do.call(rbind, 
        fit.dat)))
    if (!is.na(near.gene)) 
        return(all.data)
    fit.data <- all.data[, -(1:3)]
    if (do.average) {
        for (i in 0:(length(nearest) - 1)) {
            ref.regex <- paste("^refs\\.\\d+\\.", i, "$", sep = "")
            rat.regex <- paste("^rats\\.\\d+\\.", i, "$", sep = "")
            cor.regex <- paste("^cor\\..*\\.v\\.\\d+\\.", i, 
                "$", sep = "")
            if (i == 0) {
                ref.regex <- "^refs\\.\\d+$"
                rat.regex <- "^rats\\.\\d+$"
                cor.regex <- "^cor\\..*\\.v\\.\\d+$"
            }
            if (length(grep(ref.regex, colnames(fit.data), perl = T)) > 
                0) 
                fit.data <- cbind(ref.mean = apply(fit.data[, 
                  grep(ref.regex, colnames(fit.data), perl = T)], 
                  1, mean), fit.data)
            if (length(grep(rat.regex, colnames(fit.data), perl = T)) > 
                0) 
                fit.data <- cbind(rat.mean = apply(fit.data[, 
                  grep(rat.regex, colnames(fit.data), perl = T)], 
                  1, mean), fit.data)
            if (length(grep(cor.regex, colnames(fit.data), perl = T)) > 
                0) 
                fit.data <- cbind(cor.mean = apply(fit.data[, 
                  grep(cor.regex, colnames(fit.data), perl = T)], 
                  1, mean), fit.data)
        }
        rownames(fit.data) <- rownames(all.data)
        fit.data <- fit.data[, -c(grep("^refs", colnames(fit.data), 
            perl = T), grep("^rats", colnames(fit.data), perl = T), 
            grep("^cor\\..*\\.v", colnames(fit.data), perl = T))]
    }
    tmp.fit.data <- subset(fit.data, dist < 200)
    for (i in 0:20) {
        if (i > 0) {
            tmp.fit.data$dist[abs(pred - 0.5) < 0.1 & tmp.fit.data$dist < 
                3 * probe.spacing] <- tmp.fit.data$dist[abs(pred - 
                0.5) < 0.1 & tmp.fit.data$dist < 3 * probe.spacing] * 
                2
        }
        wts <- log10((1/(tmp.fit.data$dist/300 + 1)^2 + tmp.fit.data$prob.exp^2)/100) + 
            4
        wts[wts < 0] <- 0
        fit <- glm(is.start ~ . - dist - prob.exp, data = tmp.fit.data, 
            family = binomial(), weights = wts)
        pred <- predict(fit, newdata = tmp.fit.data, type = "response")
        cat(i, ":", extractAIC(fit)[2], sum(abs(pred - 0.5) < 
            0.1), cor(pred[fit.data$prob.exp <= -1.3], fit.data$is.start[fit.data$prob.exp <= 
            -1.3]), "\n")
    }
    print(summary(fit))
    pred <- predict(fit, newdata = fit.data, type = "response")
    prob.start.stop <- list()
    for (chr in names(brk.sites)) {
        prob.start.stop[[chr]] <- list()
        for (dir in names(brk.sites[[1]])) {
            prob.start.stop[[chr]][[dir]] <- pred[all.data$chr == 
                chr & all.data$dir == dir]
            names(prob.start.stop[[chr]][[dir]]) <- all.data$br[all.data$chr == 
                chr & all.data$dir == dir]
        }
    }
    cat("Total sites:", sum(sapply(prob.start.stop, sapply, length)), 
        "\n")
    cat("Probable start sites:", sum(sapply(prob.start.stop, 
        sapply, function(i) sum(i >= 0.5))), "\n")
    cat("Probable stop sites:", sum(sapply(prob.start.stop, sapply, 
        function(i) sum(i <= 0.5))), "\n")
    if (return.data) 
        prob.start.stop$all.data <- all.data
    prob.start.stop
}
export.data <-
function (org) 
{
    breaks <- break.sites$break.points
    breaks.start <- break.sites.start$break.points
    breaks.stop <- break.sites.stop$break.points
    break.densities <- break.densities.start <- break.densities.stop <- probe.probs <- rats.cors <- list()
    for (chr in names(chr.map)) {
        for (dir in c("FORWARD", "REVERSE")) {
            probes <- get.probes.in.window(1, 9e+09, chr = chr, 
                dir = dir, order = T)
            cat(chr, dir, length(probes), "\n")
            if (!exists("probe.spacing")) {
                tmp.rats <- index.by(rats, probes)
                probe.spacing <<- as.numeric(names(which.max(table(abs(diff(tmp.rats$POSITION))))))
                cat("Probe spacing:", probe.spacing, "\n")
                rm(tmp.rats)
            }
            probe.probs[[chr]][[dir]] <- index.by(reg.fit$pred, 
                probes)
            dd <- density(break.sites$breaks[[chr]][[dir]], bw = probe.spacing, 
                n = 2^18, from = min(break.sites$breaks[[chr]][[dir]], 
                  na.rm = T) - 100, to = max(break.sites$breaks[[chr]][[dir]], 
                  na.rm = T) + 100, na.rm = T)
            break.densities[[chr]][[dir]] <- cbind(dd$x, dd$y/max(dd$y, 
                na.rm = T))
            dd <- density(break.sites.start[[chr]][[dir]], bw = probe.spacing, 
                n = 2^18, from = min(break.sites.start[[chr]][[dir]], 
                  na.rm = T) - 100, to = max(break.sites.start[[chr]][[dir]], 
                  na.rm = T) + 100, na.rm = T)
            break.densities.start[[chr]][[dir]] <- cbind(dd$x, 
                dd$y/max(dd$y, na.rm = T))
            dd <- density(break.sites.stop[[chr]][[dir]], bw = probe.spacing, 
                n = 2^18, from = min(break.sites.stop[[chr]][[dir]], 
                  na.rm = T) - 100, to = max(break.sites.stop[[chr]][[dir]], 
                  na.rm = T) + 100, na.rm = T)
            break.densities.stop[[chr]][[dir]] <- cbind(dd$x, 
                dd$y/max(dd$y, na.rm = T))
            tmp.rats <- as.matrix(index.by(rats[, 9:ncol(rats)], 
                probes))
            cor.window <- (-3):3
            if (dir == "REVERSE") 
                wind <- cor.window
            else if (dir == "FORWARD") 
                wind <- cor.window * -1
            cor.inds <- as.matrix(expand.grid(wind[wind <= 0], 
                wind[wind >= 0]))
            cor.inds <- cor.inds[!(cor.inds[, 1] == 0 & cor.inds[, 
                2] == 0), ]
            nct <- nrow(tmp.rats)
            cor.tmp <- t(mcsapply(1:nrow(tmp.rats), function(i) {
                inds <- cor.inds + i
                apply(inds, 1, function(j) if (j[1] <= 0 || j[2] > 
                  nct) 
                  NA
                else cor(tmp.rats[j[1], ], tmp.rats[j[2], ], 
                  use = "pairwise", method = "pearson"))
            }))
            if (nrow(cor.tmp) == 1) 
                cor.tmp <- t(cor.tmp)
            rownames(cor.tmp) <- rownames(tmp.rats)
            colnames(cor.tmp) <- gsub("-", "m", paste("cor", 
                cor.inds[, 1], "v", cor.inds[, 2], sep = "."))
            rats.cors[[chr]][[dir]] <- cor.tmp
            rm(cor.tmp, tmp.rats, nct)
        }
    }
    save(org, probe.spacing, refs, rats, rats.cors, breaks, break.densities, 
        probs.expressed, probe.probs, gene.coords, breaks.start, 
        break.densities.start, breaks.stop, break.densities.stop, 
        file = paste("tilingArraySeg_", org, "_export.RData", 
            sep = ""))
}
find.peaks.in.density <-
function (brks, bw = 20, cutoff = 0.1) 
{
    if (length(brks) <= 1) 
        return(numeric())
    dens <- density(brks, n = diff(range(brks))/2, bw = bw, na.rm = T)
    dens$y <- dens$y/max(dens$y)
    lower <- c(diff(dens$y), NA)
    upper <- c(NA, diff(dens$y))
    deriv.zero <- abs(lower + upper) < sd(lower + upper, na.rm = T)
    deriv.2nd.neg <- lower < 0 & upper > 0
    dens$x[deriv.zero & deriv.2nd.neg & dens$y > cutoff & !is.na(lower) & 
        !is.na(upper)]
}
get.coding.rgns <-
function (slop = 0, chr = NA, dirs = "BOTH") 
{
    coords <- get.gene.coords()
    if (is.na(chr)) 
        chr <- names(chr.map)
    out <- list()
    for (w in chr) {
        cc <- coords[as.character(coords$where) == w, ]
        if (nrow(cc) <= 0) {
            tmp <- range(subset(rats, SEQ_ID == chr.map[w])$POSITION)
            coo <- rep(FALSE, tmp[2])
            out[[w]] <- coo
            next
        }
        coo <- rep(FALSE, max(c(cc$Start, cc$Stop)))
        if (dirs == "FORWARD") 
            cc <- cc[as.character(cc$Orientation) == "For", ]
        else if (dirs == "REVERSE") 
            cc <- cc[as.character(cc$Orientation) == "Rev", ]
        for (i in 1:nrow(cc)) coo[cc$Start[i]:cc$Stop[i]] <- TRUE
        if (slop != 0) {
            for (i in 1:slop) {
                c.up <- c(FALSE, coo[1:(length(coo) - 1)])
                c.down <- c(coo[2:length(coo)], FALSE)
                if (slop > 0) 
                  coo <- coo & c.up & c.down
                else if (slop < 0) 
                  coo <- coo | c.up | c.down
            }
        }
        out[[w]] <- coo
    }
    out
}
get.gene.coords <-
function (get.rnas = T) 
{
    gene.coords
}
get.orf.prob.expressed <-
function (genes = NA) 
{
    if (!exists("gene.coords")) 
        gene.coords <- get.gene.coords()
    if (is.na(genes)) 
        genes <- unique(as.character(gene.coords$canonical_Name))
    p.vals <- numeric()
    tmp.sel <- list()
    if (!require(multicore, quietly = T, warn.conflicts = F)) 
        mclapply <- lapply
    p.vals <- mclapply(as.list(genes), function(gene) {
        tmp <- subset(gene.coords, canonical_Name == gene)
        if (nrow(tmp) == 0) 
            tmp <- subset(gene.coords, Gene_Name == gene)
        if (nrow(tmp) == 0) {
            warning("Invalid gene: ", gene)
            return(NA)
        }
        chr <- as.character(tmp$where[1])
        start <- min(c(tmp$Start[1], tmp$Stop[1]))
        end <- max(c(tmp$Start[1], tmp$Stop[1]))
        dir <- if (tmp$Orientation[1] == "For") 
            "REVERSE"
        else "FORWARD"
        probes <- get.probes.in.window(low = start, upp = end, 
            dir = dir, chr = chr)
        if (is.null(probes)) 
            return(NA)
        pred <- reg.fit$pred[probes]
        pred <- pred[!is.na(pred)]
        p.val <- sum(log10(1 - pred))/length(probes)
        p.val[p.val == -Inf] <- -20
        if (which(genes == gene)%%25 == 0) 
            cat(gene, ":", start, end, chr, dir, "\t", p.val, 
                "\n")
        p.val
    })
    p.vals <- unlist(p.vals)
    names(p.vals) <- genes
    p.vals
}
get.probes.in.window <-
function (low, upp, chr = "Chr", dir = "both", order = T) 
{
    if (low > upp) {
        tmp <- low
        low <- upp
        upp <- tmp
    }
    tmp <- subset(rats, POSITION %betw% c(low, upp))
    if (dir != "both") {
        tmp <- subset(tmp, GENE_EXPR_OPTION == dir)
        if (chr != "all") {
            tmp <- subset(tmp, SEQ_ID == chr.map[chr])
        }
    }
    else {
        if (chr != "all") {
            tmp <- subset(tmp, SEQ_ID == chr.map[chr])
        }
    }
    out <- as.character(tmp$PROBE_ID)
    if (order) 
        out <- out[order(tmp$POSITION)]
    out
}
get.rgn.prob.expressed <-
function (low, upp, dir, chr) 
{
    probes <- get.probes.in.window(low = low, upp = upp, dir = dir, 
        chr = chr)
    pred <- reg.fit$pred[probes]
    pred <- pred[!is.na(pred)]
    p.val <- sum(log10(1 - pred))/length(probes)
    cat(gene, ":", start, end, chr, dir, "\t", p.val, "\n")
    p.vals[gene] <- p.val
}
grey.image <-
function (mat, n.gray = 32, x = 1:nrow(mat), y = 1:ncol(mat), 
    col = gray((0:n.gray)/n.gray), ...) 
{
    image(x, y, mat, col = col, ...)
}
index.by <-
function (x, y) 
{
    ttmp1 <- rank(y)
    if (is.vector(x)) {
        ttmp2 <- x[names(x) %in% y]
        ttmp2 <- ttmp2[order(names(ttmp2))]
        return(ttmp2[ttmp1])
    }
    else {
        ttmp2 <- x[rownames(x) %in% y, ]
        ttmp2 <- ttmp2[order(rownames(ttmp2)), ]
        return(ttmp2[ttmp1, ])
    }
}
joint.breakpoints <-
function (chrs = NA, dirs = NA, lower = NA, upper = NA, cor.window = (-3):3, 
    weights = c(ref = 1, rat = 1, cor = 0.05, pred = 0), weights.add = NA, 
    window = 10000, overlap = 2, n.boot = 1, cp = 1e-09, return.deltas = F, 
    xv = "1se", xval = 100, xvmult = 10, minbucket = 5, boot.vary.pos = probe.spacing, 
    log.refs = F, ...) 
{
    in.args <- c(mget(names(formals()), env = as.environment(-1)), 
        sapply(as.list(substitute({
            ...
        })[-1]), deparse))
    in.args$function.called <- joint.breakpoints
    if (!require(multicore, quietly = T, warn.conflicts = F)) 
        mclapply <- lapply
    if (is.na(window) && !is.na(lower) && !is.na(upper)) 
        window <- upper - lower + 1000
    if (is.na(chrs)) 
        chrs <- names(chr.map)
    if (is.na(dirs)) 
        dirs <- c("REVERSE", "FORWARD")
    if (!is.na(weights.add)) 
        weights[names(weights.add)] <- weights.add
    print(weights)
    all.breaks <- all.deltas <- list()
    low <- lower
    upp <- upper
    if (is.na(low)) 
        low <- -99999999999
    if (is.na(upp)) 
        upp <- 99999999999
    for (chr in chrs) {
        all.breaks[[chr]] <- all.deltas[[chr]] <- list()
        for (dir in dirs) {
            posns <- rats$POSITION
            selected <- rats$SEQ_ID == chr.map[chr] & posns %betw% 
                c(low, upp)
            selected <- which(selected & rats$GENE_EXPR_OPTION == 
                dir)
            cat(chr, dir, length(selected), "\n")
            selected <- selected[order(posns[selected])]
            posns <- posns[selected]
            col.wts <- numeric()
            data <- data.frame(posns, rownames = rownames(rats)[selected])
            if (weights["ref"] > 0) {
                tmp.refs <- as.matrix(refs[selected, 9:ncol(refs), 
                  drop = F])
                if (exists("refs.seq")) {
                  tmp.refs.seq <- as.matrix(refs.seq[selected, 
                    9:ncol(refs.seq), drop = F])
                  tmp.refs <- tmp.refs - tmp.refs.seq
                }
                colnames(tmp.refs) <- paste("refs", 1:ncol(tmp.refs), 
                  sep = ".")
                if (log.refs) 
                  tmp.refs <- log(tmp.refs)
                data <- cbind(data, tmp.refs)
                col.wts <- c(col.wts, rep(weights["ref"]/ncol(tmp.refs), 
                  ncol(tmp.refs)))
            }
            tmp.rats <- NULL
            if (weights["rat"] > 0 || weights["cor"] > 0) {
                tmp.rats <- as.matrix(rats[selected, 9:ncol(rats), 
                  drop = F])
                colnames(tmp.rats) <- paste("rats", 1:ncol(tmp.rats), 
                  sep = ".")
                data <- cbind(data, tmp.rats)
                col.wts <- c(col.wts, rep(weights["rat"]/ncol(tmp.rats), 
                  ncol(tmp.rats)))
            }
            pred.tmp <- NULL
            if (exists("reg.fit") && weights["pred"] > 0) {
                pred.tmp <- reg.fit$pred[selected]
                data <- cbind(data, pred = pred.tmp)
                col.wts <- c(col.wts, weights["pred"])
            }
            cor.tmp <- NULL
            if (weights["cor"] > 0) {
                nct <- nrow(tmp.rats)
                if (dir == "REVERSE") 
                  wind <- cor.window
                else if (dir == "FORWARD") 
                  wind <- cor.window * -1
                cor.inds <- as.matrix(expand.grid(wind[wind <= 
                  0], wind[wind >= 0]))
                cor.inds <- cor.inds[!(cor.inds[, 1] == 0 & cor.inds[, 
                  2] == 0), ]
                cor.tmp <- t(mcsapply(1:nrow(tmp.rats), function(i) {
                  inds <- cor.inds + i
                  apply(inds, 1, function(j) if (j[1] <= 0 || 
                    j[2] > nct) 
                    NA
                  else cor(tmp.rats[j[1], ], tmp.rats[j[2], ], 
                    use = "pairwise", method = "pearson"))
                }))
                if (nrow(cor.tmp) == 1) 
                  cor.tmp <- t(cor.tmp)
                rownames(cor.tmp) <- rownames(tmp.rats)
                colnames(cor.tmp) <- gsub("-", "m", paste("cor", 
                  cor.inds[, 1], "v", cor.inds[, 2], sep = "."))
                data <- cbind(data, cor.tmp)
                col.wts <- c(col.wts, rep(weights["cor"]/ncol(cor.tmp), 
                  ncol(cor.tmp)))
            }
            tmp <- mclapply(as.list(seq(min(posns), max(posns), 
                by = window)), function(w) {
                breaks <- numeric()
                out.deltas <- NULL
                wh <- which(posns %betw% c(w, w + (window * overlap)))
                dat <- data[wh, ]
                dat[, -(1:2)] <- apply(dat[, -(1:2), drop = F], 
                  2, function(i) (i - mean(i, na.rm = T))/sd(i, 
                    na.rm = T))
                if (ncol(dat) > 3) 
                  dat[, -(1:2)] <- t(apply(dat[, -(1:2), drop = F], 
                    1, "*", col.wts))
                else dat[, -(1:2)] <- apply(dat[, -(1:2), drop = F], 
                  1, "*", col.wts)
                low <- lower
                upp <- upper
                if (is.na(low)) 
                  low <- min(posns)
                if (is.na(upp)) 
                  upp <- max(posns)
                if (!any(posns[wh] %betw% c(low, upp))) {
                  cat(chr, dir, w, "\n")
                  next
                }
                sample.vec <- c(1, -1)
                sample.prob <- c(0.5, 0.5)
                for (boot in 1:n.boot) {
                  if (boot > 1) {
                    res <- resid
                    for (i in 1:ncol(res)) res[, i] <- res[, 
                      i] * sample(sample.vec, nrow(res), prob = sample.prob, 
                      replace = T)
                    dat[, -(1:2)] <- the.fit + res
                    if (boot.vary.pos > 0) 
                      dat[, "posns"] <- dat[, "posns"] + runif(nrow(dat), 
                        -boot.vary.pos, boot.vary.pos)
                  }
                  tr <- NULL
                  capture.output(tr <- my.mvpart(data.matrix(dat[, 
                    -(1:2)]) ~ posns, dat, keep.y = F, plot = F, 
                    na.action = na.exclude, cp = cp, xv = xv, 
                    xval = xval, xvmult = xvmult, minbucket = minbucket, 
                    ...))
                  new.breaks <- tr$splits[, "index"]
                  breaks <- c(breaks, new.breaks)
                  if (return.deltas) {
                    pred <- predict(tr)
                    out.deltas <- rbind(out.deltas, t(sapply(new.breaks, 
                      function(br) {
                        inds <- order(abs(dat[, "posns"] - br))[1:2]
                        if (dir == "FORWARD") 
                          inds <- rev(inds)
                        pred[inds[1], ] - pred[inds[2], ]
                      })))
                    colnames(out.deltas) <- colnames(dat[, -(1:2)])
                  }
                  cat(chr, dir, w, posns[length(posns)], length(breaks), 
                    range(breaks), boot, "\n")
                  if (boot == 1) {
                    if (n.boot > 1) {
                      the.fit <- predict(tr)
                      resid <- dat[, -(1:2)] - the.fit
                    }
                  }
                }
                return(list(breaks = breaks, deltas = out.deltas))
            })
            breaks <- numeric()
            out.deltas <- NULL
            for (i in 1:length(tmp)) {
                breaks <- c(breaks, tmp[[i]]$breaks)
                if (return.deltas) 
                  out.deltas <- rbind(out.deltas, tmp[[i]]$deltas)
            }
            rm(tmp)
            all.breaks[[chr]][[dir]] <- if (n.boot == 1) 
                unique(breaks)
            else breaks
            if (return.deltas) 
                all.deltas[[chr]][[dir]] <- out.deltas
        }
    }
    invisible(list(breaks = all.breaks, args = in.args, deltas = all.deltas))
}
mcsapply <-
function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, mc.preschedule = TRUE, 
    mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("cores")) 
{
    FUN <- match.fun(FUN)
    require(multicore, quietly = T, warn.conflicts = F)
    answer <- mclapply(X, FUN, ..., mc.preschedule = mc.preschedule, 
        mc.set.seed = mc.set.seed, mc.silent = mc.silent, mc.cores = mc.cores)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (simplify && length(answer) && length(common.len <- unique(unlist(lapply(answer, 
        length)))) == 1L) {
        if (common.len == 1L) 
            unlist(answer, recursive = FALSE)
        else if (common.len > 1L) 
            array(unlist(answer, recursive = FALSE), dim = c(common.len, 
                length(X)), dimnames = if (!(is.null(n1 <- names(answer[[1L]])) & 
                is.null(n2 <- names(answer)))) 
                list(n1, n2))
        else answer
    }
    else answer
}
my.mvpart <-
function (form, data, minauto = TRUE, size, xv = c("1se", "min", 
    "pick", "none"), xval = 10, xvmult = 0, xvse = 1, snip = FALSE, 
    plot.add = TRUE, text.add = TRUE, digits = 3, margin = 0, 
    uniform = FALSE, which = 4, pretty = TRUE, use.n = TRUE, 
    all.leaves = FALSE, bars = TRUE, legend, bord = FALSE, xadj = 1, 
    yadj = 1, prn = FALSE, branch = 1, rsq = FALSE, big.pts = FALSE, 
    pca = FALSE, interact.pca = FALSE, wgt.ave.pca = FALSE, keep.y = TRUE, 
    ...) 
{
    call <- match.call()
    number <- function(x) {
        match(x, sort(unique(x)))
    }
    cv.var <- function(x, cv = 10) {
        x <- match(x, sort(unique(x)))
        luni <- length(unique(x))
        if (luni >= cv) {
            grps <- ceiling((cv * cumsum(table(x)))/length(x))
            x <- number(grps[x])
        }
        x
    }
    if (length(xval) > 1) {
        if (xvmult > 1) 
            xvalvar <- xval
        xval <- cv.var(xval)
    }
    choice <- c("1se", "min", "pick", "none")
    xv <- choice[pmatch(xv[1], choice)]
    if (!missing(size) || xv == "none") 
        xval <- 0
    if (minauto) {
        n <- nrow(data)
        minsplit <- ceiling(log2(n))
        minbucket <- ceiling(minsplit/3)
    }
    z <- rpart(form, data = data, ...)
    if (all(z$where == 1)) {
        cat("No splits possible -- try decreasing cp\n")
        return(z)
    }
    if (!is.null(z)) {
        xval <- z$control$xval
        if (xvmult > 1) {
            zresse <- zres <- matrix(NA, nrow = nrow(z$cptable), 
                ncol = xvmult)
            zres[, 1] <- z$cptable[, 4]
            zresse[, 1] <- z$cptable[, 5]
            cat("X-Val rep : 1")
            for (i in 2:xvmult) {
                if (length(xval) == nrow(data)) 
                  xval <- cv.var(xvalvar)
                ztemp <- rpart(form, data = data, ...)$cptable[, 
                  4:5]
                zres[, i] <- ztemp[, 1]
                zresse[, i] <- ztemp[, 2]
                cat(" ", i)
                NULL
            }
            cat("\n")
            z$cptable[, 4] <- apply(zres, 1, mean)
            z$cptable[, 5] <- apply(zresse, 1, mean)
            tabmins <- apply(zres, 2, function(x, nc, sizes) {
                sizes[x == min(x)][1]
            }, nc = nrow(zres), sizes = z$cptable[, 2] + 1)
            cat("Minimum tree sizes\n")
            print(table(tabmins))
        }
        if (missing(size)) {
            if (xv == "pick") {
                if (xvmult <= 1) 
                  plotcp(z, xvse, pch = 16, col = 2)
                else plotcp(z, xvse, pch = 16, col = 2, tab = table(tabmins))
                size.loc <- locator(1)
                if (!is.null(size.loc)) {
                  splt <- round(size.loc$x)
                  if (splt < 2) 
                    splt <- 2
                  else if (splt > length(z$cptable[, 1])) 
                    splt <- length(z$cptable[, 1])
                  cpp <- z$cptable[, 1][splt]
                  z <- prune.rpart(z, cp = cpp)
                }
            }
            else if ((xv == "1se" | xv == "min") && (xval[1] != 
                0)) {
                xerror <- z$cptable[, 4]
                xstd <- z$cptable[, 5]
                if (xv == "min") 
                  splt <- min(seq(along = xerror)[xerror == min(xerror)])
                else splt <- min(seq(along = xerror)[xerror <= 
                  min(xerror) + xvse * xstd])
                if (!is.na(splt)) {
                  if (splt == 1) 
                    splt <- 2
                  cpp <- z$cptable[, 1][splt]
                  z <- prune.rpart(z, cp = cpp)
                }
                else {
                  (cat("No pruning possible : size 2 tree produced ?? \n"))
                  use.size <- TRUE
                  size <- 2
                }
            }
        }
        else {
            if (size <= 2) 
                cpp <- z$cptable[2, 1]
            else if (size >= max(z$cptable[, 2] + 1)) 
                cpp <- z$cptable[dim(z$cptable)[1], 1]
            else cpp <- z$cptable[, 1][min(abs(size - z$cptable[, 
                2] - 1)) == abs(size - z$cptable[, 2] - 1)][1]
            z <- prune.rpart(z, cp = cpp)
        }
        if (snip) {
            z <- snip.rpart(z)
        }
        if (rsq && xval != 0 && z$method != "class") {
            rsq.rpart(z)
            locator(1)
        }
        if (pca) {
            locator(1)
            rpart.pca(z, interact = interact.pca, wgt.ave = wgt.ave.pca)
        }
    }
    else {
        cat("No splits could be formed\n")
    }
    if (!is.null(z)) {
        if (!keep.y) 
            z$y <- NULL
        z$call <- call
        invisible(z)
    }
}
plot.domains.in.window <-
function (center, where, window, svglinks = F) 
{
    data.file <- "data/Pfam_UCSC.txt"
    if (!exists("domain.coords")) {
        if (!file.exists(data.file)) 
            return()
        domain.coords <- read.table(data.file, head = T)
    }
    window <- round(c(center - window/2, center + window/2))
    w.expand <- round(c(window[1] - diff(window)/5, window[2] + 
        diff(window)/5))
    in.window <- function(x, w) {
        x >= w[1] & x <= w[2]
    }
    plot(0, 0, xlim = window, ylim = c(0, 2.5), ann = F, xaxt = "n", 
        yaxt = "n", bty = "n")
    axis(side = 1, pos = 0, mgp = c(0, 0.5, 0))
    genes.in <- domain.coords[gsub("PLASMID_", "", toupper(domain.coords$chrom)) == 
        toupper(where) & (in.window(domain.coords[, "chromStart"], 
        w.expand) | in.window(domain.coords[, "chromEnd"], w.expand)), 
        ]
    if (nrow(genes.in) <= 0) 
        return()
    genes.x <- apply(genes.in[, c("chromStart", "chromEnd")], 
        1, mean)
    genes.is.fwd <- genes.in$strand == "+"
    for (i in 1:nrow(genes.in)) {
        start <- genes.in[i, "chromStart"]
        end <- genes.in[i, "chromEnd"]
        name <- genes.in[i, "name"]
        if (svglinks) {
            setSVGShapeToolTip(title = name, desc = paste("pfamAC =", 
                genes.in[i, "pfamAC"]))
            setSVGShapeURL(paste("http://pfam.sanger.ac.uk/search/keyword?query=", 
                genes.in[i, "pfamAC"], sep = ""))
        }
        if (genes.in[i, "strand"] == "+") {
            rect(start, 1, end, 1.5, col = "indianred", border = "black", 
                lwd = 1)
            text(genes.x[i], 1.25, labels = name, adj = c(0.5, 
                0.5), col = "black", cex = 0.7)
        }
        else if (genes.in[i, "strand"] == "-") {
            rect(start, 0.25, end, 0.75, col = "indianred2", 
                border = "black", lwd = 1)
            text(genes.x[i], 0.5, labels = name, adj = c(0.5, 
                0.5), col = "black", cex = 0.7)
        }
    }
}
plot.everything <-
function (gene = NA, window = NA, dirs = c("REVERSE", "FORWARD"), 
    chr = NA, lower = 10000, upper = 20000, plot.cor.matrix = F, 
    window.cor = 3, col.func = topo.colors, main = NULL, segment = F, 
    method = "pearson", mindev = 0.0075, new.dev = F, dir = "pdfs", 
    ...) 
{
    if (!is.na(gene) && length(gene) > 1) {
        if (dir != "") 
            dir.create(dir, recursive = T)
        lapply(as.list(gene), function(g, ...) {
            cairo_pdf(paste(dir, "/", g, ".pdf", sep = ""))
            plot.everything(g, window = window, dirs = dirs, 
                main = g)
            dev.off()
        })
        return()
    }
    else if (length(lower) > 1) {
        if (dir != "") 
            dir.create(dir, recursive = T)
        lapply(as.list(1:length(lower)), function(i, ...) {
            cat(lower[i], upper[i], paste(dir, "/", lower[i], 
                "-", upper[i], ".pdf", sep = ""), "\n")
            cairo_pdf(paste(dir, "/", lower[i], "-", upper[i], 
                ".pdf", sep = ""))
            plot.everything(lower = lower[i], upper = upper[i], 
                window = window, dirs = dirs, main = paste(lower[i], 
                  "-", upper[i]))
            dev.off()
        })
        return()
    }
    plot.genes.in.window <- function(center, where, window, svglinks = F) {
        window <- round(c(center - window/2, center + window/2))
        w.expand <- round(c(window[1] - diff(window)/5, window[2] + 
            diff(window)/5))
        in.window <- function(x, w) {
            x >= w[1] & x <= w[2]
        }
        plot(0, 0, xlim = window, ylim = c(0, 2.5), ann = F, 
            xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
        axis(side = 1, pos = 1, mgp = c(0, 0.5, 0), cex.axis = 0.7)
        axis(side = 3, pos = 1, mgp = c(0, 0.5, 0), cex.axis = 0.7)
        genes.in <- gene.coords[gene.coords[, "where"] == where & 
            (in.window(gene.coords[, "Start"], w.expand) | in.window(gene.coords[, 
                "Stop"], w.expand)), ]
        if (nrow(genes.in) <= 0) 
            return()
        genes.x <- apply(genes.in[, c("Start", "Stop")], 1, mean)
        genes.is.fwd <- genes.in[, "Orientation"] == "For"
        for (i in 1:nrow(genes.in)) {
            start <- genes.in[i, "Start"]
            end <- genes.in[i, "Stop"]
            name <- genes.in[i, "Gene_Name"]
            if (svglinks) {
                setSVGShapeToolTip(title = name, desc = genes.in[i, 
                  "canonical_Name"])
            }
            if (genes.in[i, "Orientation"] == "For") {
                rect(start, 1, end, 1.5, col = "yellow", border = "black", 
                  lwd = 1)
                text(genes.x[i], 1.25, labels = name, adj = c(0.5, 
                  0.5), col = "black", cex = 0.7)
            }
            else if (genes.in[i, "Orientation"] == "Rev") {
                rect(start, 0.25, end, 0.75, col = "orange", 
                  border = "black", lwd = 1)
                text(genes.x[i], 0.5, labels = name, adj = c(0.5, 
                  0.5), col = "black", cex = 0.7)
            }
        }
    }
    if (!is.na(gene)) {
        if (!exists("gene.coords")) 
            gene.coords <- get.gene.coords()
        tmp <- gene.coords[gene.coords$canonical_Name == gene, 
            , drop = F]
        if (nrow(tmp) == 0) 
            tmp <- gene.coords[gene.coords$Gene_Name == gene, 
                , drop = F]
        if (nrow(tmp) == 0) {
            warning("Invalid gene:", gene)
            return()
        }
        if (is.na(window)) 
            window <- 5000
        start <- min(c(tmp$Start[1], tmp$Stop[1]))
        end <- max(c(tmp$Start[1], tmp$Stop[1]))
        lower <- start - window
        upper <- end + window
        if (is.na(chr)) 
            chr <- as.character(tmp$where[1])
        cat(gene, ":", lower, upper, chr, "\n")
    }
    if (is.na(chr)) 
        chr <- names(chr.map)[1]
    if (!exists("probe.spacing")) {
        tmp.rats <- rats[get.probes.in.window(chr = chr, dir = "FORWARD", 
            low = 0, upp = 9e+09, order = T), ]
        probe.spacing <<- as.numeric(names(which.max(table(abs(diff(tmp.rats$POSITION))))))
        cat("Probe spacing:", probe.spacing, "\n")
        rm(tmp.rats)
    }
    if (segment) {
        seg.wind <- 10000
        break.sites <- joint.breakpoints(lower = lower, upper = upper, 
            window = seg.wind, chr = chr, ...)
        break.sites$break.points <- lapply(break.sites$breaks, 
            lapply, find.peaks.in.density, bw = probe.spacing, 
            cutoff = 0.2)
    }
    probes.dir <- list()
    for (dir in dirs) {
        probes.dir[[dir]] <- get.probes.in.window(chr = chr, 
            dir = dir, low = lower, upp = upper, order = T)
    }
    tmp <- cors <- list()
    for (dir in dirs) {
        tmp[[dir]] <- index.by(rats, probes.dir[[dir]])
        qtmp <- as.matrix(tmp[[dir]][, 9:ncol(tmp[[dir]]), drop = F])
        if (nrow(qtmp) > 1) {
            if (method != "cor.test") {
                cors[[dir]] <- cor(t(qtmp), method = method, 
                  use = "pairwise")
            }
            else {
                ttmp <- matrix(nrow = nrow(qtmp), ncol = nrow(qtmp))
                for (i in 2:(nrow(ttmp) - 1)) {
                  for (j in (i + 1):(nrow(ttmp) - 1)) ttmp[i, 
                    j] <- ttmp[j, i] <- cor.test(qtmp[i, ], qtmp[j, 
                    ])$p.value
                }
                cors[[dir]] <- -log(ttmp + 1e-09)
            }
        }
    }
    dcors <- ccors <- list()
    if (length(cors) >= 1) {
        ccors <- cors[[1]]
        diag(ccors) <- NA
        if (length(dirs) > 1) {
            if (nrow(cors[[2]]) > nrow(ccors)) 
                cors[[2]] <- cors[[2]][-nrow(cors[[2]]), -nrow(cors[[2]])]
            else if (nrow(ccors) > nrow(cors[[2]])) 
                ccors <- ccors[-nrow(ccors), -nrow(ccors)]
            ccors[lower.tri(ccors)] <- cors[[2]][lower.tri(cors[[2]])]
        }
        for (dir in dirs) {
            cl <- cors
            if (dir == "REVERSE") {
                dcor <- sapply(1:window.cor, function(i) c(cl[[dir]][cbind(1:(ncol(cl[[dir]]) - 
                  i), (i + 1):ncol(cl[[dir]]))], rep(NA, i - 
                  1)))
                dcors[[dir]] <- apply(dcor, 1, mean, na.rm = T)
            }
            else {
                dcor <- sapply(1:window.cor, function(i) c(cl[[dir]][cbind(rev(1:(ncol(cl[[dir]]) - 
                  i)), rev((i + 1):ncol(cl[[dir]])))], rep(NA, 
                  i - 1)))
                dcors[[dir]] <- rev(apply(dcor, 1, mean, na.rm = T))
            }
        }
    }
    new.dev <- new.dev || (dev.cur() == 1)
    if (new.dev) 
        dev.new(height = 10)
    if (plot.cor.matrix) {
        layout(t(matrix(c(1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
            10, 11, 12), 1, 15)))
    }
    else {
        layout <- 1:10
        if (exists("tfbd.raw")) 
            layout <- c(layout, max(layout) + 1, max(layout) + 
                2)
        layout(t(matrix(layout, 1, length(layout))))
    }
    old.par <- par(oma = c(1, 2, 1, 1), mar = c(0.3, 1.5, 0.3, 
        0.3), mgp = c(3, 1, 0) * 0.1, lwd = 1, xaxs = "i")
    on.exit(par(old.par))
    if (plot.cor.matrix) {
        grey.image(exp(ccors), x = tmp[[1]]$POSITION[1:nrow(ccors)], 
            y = tmp[[1]]$POSITION[1:nrow(ccors)], col = col.func(256))
    }
    plot.genes.in.window(center = mean(c(lower, upper)), where = chr, 
        window = upper - lower, svglinks = (names(dev.cur()) == 
            "devSVG"))
    plot.domains.in.window(center = mean(c(lower, upper)), where = chr, 
        window = upper - lower, svglinks = (names(dev.cur()) == 
            "devSVG"))
    if (!is.null(main)) 
        text(mean(c(lower, upper)), 0, main, adj = c(0.5, 0), 
            font = 2)
    zlim <- NA
    for (dir in dirs) {
        plot.ratios.matrix(tmp[[dir]], NULL, NULL, lower, upper, 
            col.func, zlim = zlim)
    }
    if (length(dcors) >= 1) {
        for (dir in dirs) {
            plot(tmp[[dir]]$POSITION[1:length(dcors[[dir]])], 
                exp(dcors[[dir]]), pch = 20, cex = 0.5, xaxt = "n", 
                ylab = "", ylim = exp(c(-1, 1)))
            points(tmp[[dir]]$POSITION[1:length(dcors[[dir]])][dcors[[dir]] > 
                0.9], exp(dcors[[dir]][dcors[[dir]] > 0.9]), 
                pch = 20, cex = 0.5, col = "red")
            if (FALSE) {
                require(tree)
                grid <- data.frame(x.in = tmp[[dir]]$POSITION[1:length(dcors[[dir]])], 
                  y.in = dcors[[dir]])
                tr <- tree(y.in ~ x.in, grid, minsize = 5, mindev = mindev * 
                  20000/abs(upper - lower))
                fitted <- predict(tr, grid)
                lines(grid$x.in, exp(fitted), col = "red")
            }
        }
    }
    if (exists("refs")) {
        yrng <- quantile(as.matrix(refs[, 9:ncol(refs)]), probs = c(0.001, 
            0.999), na.rm = T)
        for (dir in dirs) {
            tmp <- index.by(refs, probes.dir[[dir]])
            ttmp <- as.matrix(tmp[, 9:ncol(tmp), drop = F])
            y <- apply(ttmp, 1, mean, na.rm = T)
            yerr <- apply(ttmp, 1, sd, na.rm = T)
            if (exists("refs.seq")) {
                tmp.q <- index.by(refs.seq, probes.dir[[dir]])[, 
                  9:ncol(refs.seq)]
                tmp.q2 <- as.matrix(tmp.q)
                ttmp2 <- ttmp - tmp.q2
                y2 <- apply(ttmp2, 1, mean, na.rm = T)
                yerr2 <- apply(ttmp2, 1, sd, na.rm = T)
                plot(tmp$POSITION, y, pch = 20, cex = 0.5, xaxt = "n", 
                  ylim = yrng, col = "gray", xlab = "", ylab = "")
                arrows(tmp$POSITION, y + yerr, tmp$POSITION, 
                  y - yerr, length = 0, col = "gray")
                points(tmp$POSITION, y2, pch = 20, cex = 0.5)
                arrows(tmp$POSITION, y2 + yerr2, tmp$POSITION, 
                  y2 - yerr2, length = 0)
                rm(tmp.q, tmp.q2, ttmp2, y2, yerr2)
                gc()
            }
            else {
                plot(tmp$POSITION, y, pch = 20, cex = 0.5, xaxt = "n", 
                  ylim = yrng, xlab = "", ylab = "")
                arrows(tmp$POSITION, y + yerr, tmp$POSITION, 
                  y - yerr, length = 0)
            }
            if (exists("refs.seg")) {
                tmp.q <- index.by(refs.seg, probes.dir[[dir]])[, 
                  9:ncol(refs.seg)]
                ttmp2 <- as.matrix(tmp.q)
                y3 <- apply(ttmp2, 1, mean, na.rm = T)
                lines(tmp$POSITION, y3, col = "red")
                rm(tmp.q, ttmp2, y3)
                gc()
            }
        }
    }
    for (dir in dirs) {
        if (exists("reg.fit") && class(reg.fit) != "try-error") {
            tmp <- index.by(rats, probes.dir[[dir]])$POSITION
            pred <- index.by(reg.fit$pred, probes.dir[[dir]])
            if (is.data.frame(pred)) 
                pred <- pred[, 1]
            plot(tmp, pred, typ = "n", pch = 20, cex = 0.5, xaxt = "n", 
                ylim = c(-0.1, 1.1), xlab = "", ylab = "")
            points(tmp, pred, pch = 20, cex = 0.5, typ = "b")
            points(tmp[pred >= 0.8], pred[pred >= 0.8], pch = 20, 
                cex = 0.5, col = "palevioletred")
            points(tmp[pred >= 0.9], pred[pred >= 0.9], pch = 20, 
                cex = 0.5, col = "violetred")
            posn.range <- range(tmp)
            brk.cols <- c(break.sites = "green", break.sites.stop = "red", 
                break.sites.start = "blue")
            for (brk in names(brk.cols)) {
                if (exists(brk)) {
                  brks <- get(brk)
                  if (!is.null(brks$breaks)) 
                    pks <- brks$breaks[[chr]][[dir]]
                  else if (!is.null(brks$Chr)) 
                    pks <- brks[[chr]][[dir]]
                  if (is.data.frame(pks)) 
                    pks <- pks[, 1]
                  pks <- pks[pks %betw% posn.range]
                  if (length(pks) > 0) {
                    dd <- density(pks, bw = probe.spacing, n = 16384, 
                      from = posn.range[1], to = posn.range[2], 
                      na.rm = T)
                    lines(dd$x, dd$y/max(dd$y), col = brk.cols[brk])
                  }
                }
            }
        }
        else {
            tmp <- index.by(rats, probes.dir[[dir]])$POSITION
            plot(tmp, rep(1, length(tmp)), typ = "n", pch = 20, 
                cex = 0.5, xaxt = "n", ylim = c(-0.1, 1.1), xlab = "", 
                ylab = "")
        }
        if (exists("break.sites") || exists("break.sites.start") || 
            exists("break.sites.stop") || exists("hand.sites")) {
            posn.range <- c(lower, upper)
            brk.cols <- c(break.sites = "green", break.sites.stop = "red", 
                break.sites.start = "blue")
            for (brk in names(brk.cols)) {
                if (exists(brk)) {
                  brks <- get(brk)
                  if (!is.null(brks$break.points)) {
                    pks <- brks$break.points[[chr]][[dir]]
                    if (is.data.frame(pks)) 
                      pks <- pks[, 1]
                    pks <- pks[pks %betw% posn.range]
                    if (length(pks) > 0) {
                      arrows(pks, -999, pks, 999, col = brk.cols[brk], 
                        length = 0)
                      arrows(pks, -999, pks, 999, col = brk.cols[brk], 
                        lty = if (dir == "REVERSE") 
                          3
                        else 4, length = 0, xpd = NA)
                    }
                  }
                }
            }
        }
        if (exists("hand.sites")) {
            chr2 <- if (chr == "Chr") 
                "chromosome"
            else chr
            dir2 <- if (dir == "REVERSE") 
                "For"
            else "Rev"
            hand.brks <- hand.sites[hand.sites$chr == chr2 & 
                hand.sites$dir == dir2, ]
            hand.brks <- hand.brks[hand.brks$start %betw% posn.range | 
                hand.brks$end %betw% posn.range, ]
            hand.brks <- hand.brks[!is.na(hand.brks$i), ]
            lty <- ifelse(hand.brks$annotator == "dave", 2, 3)
            lwd <- ifelse(hand.brks$flag == TRUE, 1, 3)
            if (length(hand.brks$start) > 0) 
                arrows(hand.brks$start, -999, hand.brks$start, 
                  999, col = "blue", lwd = lwd, lty = lty, length = 0, 
                  xpd = NA)
            if (length(hand.brks$end) > 0) 
                arrows(hand.brks$end, -999, hand.brks$end, 999, 
                  col = "red", lwd = lwd, lty = lty, length = 0, 
                  xpd = NA)
        }
    }
    invisible(list(cors = ccors, dcors = dcors))
}
plot.ratios.matrix <-
function (rats, dir = NULL, chr = NULL, lower = 10000, upper = 20000, 
    col.func = topo.colors, zlim = NA) 
{
    if (is.null(dir)) 
        tmp <- rats
    else {
        probes <- get.probes.in.window(lower, upper, chr, dir)
        tmp <- index.by(rats, probes)
    }
    posns <- as.numeric(tmp$POSITION)
    if (any(diff(posns) == 0)) {
        tmp <- tmp[diff(posns) > 0, ]
        posns <- posns[diff(posns) > 0]
    }
    if (is.na(zlim[1])) 
        zlim <- quantile(as.matrix(rats[, 9:ncol(rats)]), probs = c(0.001, 
            0.999), na.rm = T)
    grey.image(as.matrix(tmp[, 9:ncol(tmp)]), x = posns, col = col.func(256), 
        xaxt = "n", yaxt = "n", xlab = "", ylab = "", zlim = zlim)
    if (exists("rats.seg")) {
        ttmp <- rats.seg[rownames(tmp), 9:ncol(rats.seg)]
        breaks <- apply(ttmp, 2, function(i) diff(i) != 0)
        sapply(1:ncol(breaks), function(i) {
            brk <- which(breaks[, i])
            if (length(brk) > 0) {
                arrows(brk, i - 0.5, brk, i + 0.5, length = 0, 
                  col = "red")
            }
        })
    }
}
regress.probe.expression.probs <-
function (window.cors = (-3):3, window.refs = (-3):3, slop = -50, 
    use.reference.option = "raw", regression.opt = "glm", interactions = F, 
    add.to.formula = "", include.log.od = F, include.rev = F, 
    n.iters = 25, in.fit = NULL, ...) 
{
    in.args <- c(mget(names(formals()), env = as.environment(-1)), 
        sapply(as.list(substitute({
            ...
        })[-1]), deparse))
    in.args$function.called <- regress.probe.expression.probs
    if (is.null(in.fit)) {
        coding.frac <- numeric()
        params <- NULL
        window.cors <- window.cors[window.cors != 0]
        for (chr in rev(names(chr.map))) {
            for (dir in c("REVERSE", "FORWARD")) {
                tmp.par <- NULL
                selected <- which(rats$GENE_EXPR_OPTION == dir & 
                  rats$SEQ_ID == chr.map[chr])
                cat(chr, dir, length(selected), "\n")
                tmp <- rats[selected, ]
                tmp <- tmp[order(tmp$POSITION), ]
                if (!is.na(window.cors) && length(window.cors) > 
                  0) {
                  tmp2 <- t(as.matrix(tmp[, 9:ncol(tmp)]))
                  nct <- ncol(tmp2)
                  if (dir == "REVERSE") 
                    window <- window.cors
                  else if (dir == "FORWARD") 
                    window <- window.cors * -1
                  cor.tmp <- t(mcsapply(1:ncol(tmp2), function(i) {
                    wind <- window + i
                    wind[wind < 1 | wind > nct] <- NA
                    cor(tmp2[, i], tmp2[, wind], use = "pairwise", 
                      method = "pearson")
                  }))
                  rownames(cor.tmp) <- rownames(tmp)
                  colnames(cor.tmp) <- gsub("-", "m", paste("cor.v.", 
                    window.cors, sep = ""))
                  tmp.par <- cbind(tmp.par, cors = cor.tmp)
                  if (!include.log.od) 
                    rm(cor.tmp, tmp2)
                }
                if (include.log.od) {
                  rownames(tmp2) <- sapply(strsplit(rownames(tmp2), 
                    ".", fixed = T), function(i) paste(i[9:11], 
                    collapse = "."))
                  tmp.par <- cbind(tmp.par, t(abs(tmp2)))
                  rm(cor.tmp, tmp2)
                }
                if (use.reference.option == "segments") 
                  tmp <- refs.seg[selected, ]
                else tmp <- refs[selected, ]
                tmp <- tmp[order(tmp$POSITION), ]
                tmp2 <- t(as.matrix(tmp[, 9:ncol(tmp)]))
                if (use.reference.option == "raw.m.seq") {
                  tmp.q <- refs.seq[selected, ]
                  tmp.q <- tmp.q[order(tmp.q$POSITION), ]
                  tmp.q2 <- t(as.matrix(tmp.q[, 9:ncol(tmp.q)]))
                  tmp2 <- tmp2 - tmp.q2
                  rm(tmp.q, tmp.q2)
                }
                i.tmp <- t(apply(tmp2, 2, mean, na.rm = T))
                if (dir == "REVERSE") 
                  window <- window.refs
                else if (dir == "FORWARD") 
                  window <- window.refs * -1
                nct <- ncol(tmp2)
                ref.tmp <- t(mcsapply(1:ncol(tmp2), function(i) {
                  wind <- window + i
                  wind[wind < 1 | wind > nct] <- NA
                  i.tmp[, wind]
                }))
                if (length(window) == 1) 
                  ref.tmp <- t(ref.tmp)
                rownames(ref.tmp) <- rownames(tmp)
                colnames(ref.tmp) <- gsub("-", "m", paste("ref.", 
                  window.refs, sep = ""))
                tmp.par <- cbind(tmp.par, refs = ref.tmp)
                rm(ref.tmp, tmp2, i.tmp)
                real.dir <- if (dir == "REVERSE") 
                  "FORWARD"
                else "REVERSE"
                coding <- get.coding.rgns(slop = slop, chr = chr, 
                  dir = real.dir)[[1]]
                if (length(coding) <= 0) {
                  cat("Skipping.\n")
                  next
                }
                params <- rbind(params, tmp.par)
                rm(tmp.par)
                gc()
                tmp3 <- as.matrix(tmp[, 4:6])
                coding.tmp <- mcsapply(1:nrow(tmp3), function(i) sum(coding[tmp3[i, 
                  1]:tmp3[i, 2]], na.rm = T))
                names(coding.tmp) <- rownames(tmp3)
                coding.tmp[coding.tmp > as.integer(as.character(tmp$PROBE_LENGTH))] <- 0
                coding.frac <- c(coding.frac, coding.tmp)
                rm(coding.tmp)
                gc()
            }
        }
        if (include.rev) {
            rev.pars <- grep("cor\\.|ref\\.", colnames(params), 
                perl = T)
            new.rev.pars <- paste("rev", colnames(params)[rev.pars], 
                sep = ".")
            new.par <- NULL
            for (chr in rev(names(chr.map))) {
                selected.rev <- which(rats$GENE_EXPR_OPTION == 
                  "REVERSE" & rats$SEQ_ID == chr.map[chr])
                selected.for <- which(rats$GENE_EXPR_OPTION == 
                  "FORWARD" & rats$SEQ_ID == chr.map[chr])
                par.rev <- params[selected.rev, ][order(rats$POSITION[selected.rev]), 
                  ]
                par.for <- params[selected.for, ][order(rats$POSITION[selected.for]), 
                  ]
                nr.for <- nrow(par.for)
                nr.rev <- nrow(par.rev)
                new.par <- if (nr.rev < nr.for) {
                  rbind(new.par, cbind(par.for, rbind(par.rev[, 
                    rev.pars], rep(NA, length(rev.pars)))))
                }
                else {
                  rbind(new.par, cbind(par.for, par.rev[1:nrow(par.for), 
                    rev.pars]))
                }
                new.par <- if (nr.for < nr.rev) {
                  rbind(new.par, cbind(par.rev, rbind(par.for[, 
                    rev.pars], rep(NA, length(rev.pars)))))
                }
                else {
                  rbind(new.par, cbind(par.rev, par.for[1:nrow(par.rev), 
                    rev.pars]))
                }
            }
            colnames(new.par)[(ncol(params) + 1):ncol(new.par)] <- new.rev.pars
            new.par <- new.par[rownames(params), ]
            params <- new.par
        }
        if (any(names(coding.frac) != rownames(rats))) 
            coding.frac <- coding.frac[rownames(rats)]
        if (any(rownames(params) != rownames(rats))) 
            params <- params[rownames(rats), ]
        rownames(params) <- names(coding.frac) <- rownames(rats)
    }
    else {
        params <- in.fit$params
        coding.frac <- in.fit$coding.frac
    }
    cat("PARAMS:", dim(params), colnames(params), "\n")
    out.pred <- NULL
    coding.frac.cutoff <- 50
    prob.cutoff <- 0.8
    data <- data.frame(coding.frac = coding.frac/max(coding.frac, 
        na.rm = T), params)
    form <- paste("coding.frac ~", paste(colnames(data)[-1], 
        collapse = "+"))
    if (add.to.formula != "") 
        form <- paste(form, "+", add.to.formula)
    cat("Using formula:", form, "\n")
    if (!interactions) {
        formula <- as.formula(form)
    }
    else {
        formula <- as.formula(paste("coding.frac ~ (", paste(colnames(data)[-1], 
            collapse = "+"), ")^2"))
    }
    for (iter in 1:n.iters) {
        weights <- rep(1, length(coding.frac))
        names(weights) <- names(coding.frac)
        new.weight <- (n.iters - iter + 1)/n.iters
        if (iter > 1) {
            weights[coding.frac < coding.frac.cutoff & pred >= 
                prob.cutoff] <- new.weight
        }
        if (regression.opt == "glm") {
            fit <- glm(formula, data = data, weights = weights, 
                family = binomial(), trace = F, na.action = na.exclude)
            pred <- predict(fit, type = "response", na.action = na.exclude)
        }
        if (!all(names(pred) == rownames(params))) 
            pred <- index.by(pred, rownames(params))
        names(pred) <- rownames(params)
        out.pred <- cbind(out.pred, pred)
        reg.fit <- list(params = params, coding.frac = coding.frac, 
            fit = fit, pred = pred, pred.iter = out.pred, args = in.args)
        if (regression.opt == "glm") 
            cat(iter, new.weight, slop, range(window.cors), extractAIC(reg.fit$fit), 
                log(sum(reg.fit$pred[reg.fit$coding.frac >= coding.frac.cutoff], 
                  na.rm = T) + sum(1 - reg.fit$pred[reg.fit$coding.frac < 
                  coding.frac.cutoff], na.rm = T)), "\n")
    }
    print(summary(fit))
    invisible(reg.fit)
}
remove.duplicate.break.sites <-
function (break.sites, dist = 40) 
{
    for (chr in names(break.sites)) {
        for (dir in names(break.sites[[chr]])) {
            brks <- break.sites[[chr]][[dir]]
            newbrk <- numeric()
            while (length(brks) > 1) {
                hits <- which(abs(brks - brks[1]) < dist)
                if (length(hits) > 1) {
                  cat("Merging hits at", chr, dir, brks[hits], 
                    "\n")
                  posn <- mean(brks[hits])
                  newbrk <- c(newbrk, posn)
                }
                else {
                  newbrk <- c(newbrk, brks[1])
                }
                brks <- brks[-hits]
            }
            break.sites[[chr]][[dir]] <- newbrk
        }
    }
    break.sites
}
runnit <-
function (org = "mm", combine.dye.swap = T, plot = T, export = T, 
    ...) 
{
    try(X11.options(type = "Xlib"))
    if (!exists("rats")) {
        if (org == "halo") {
            load("all_data_halo.RData")
            load("halo_coords.RData")
            gene.coords <- halo.coords
            chr.map <- c(Chr = "NC_002607", pNRC100 = "NC_001869", 
                pNRC200 = "NC_002608")
        }
        else {
            rdata.file <- list.files(path = sprintf("./data/%s/", 
                org), patt = glob2rx("*all_intensity.norm*.Rdata"), 
                full = T)
            load(rdata.file)
            if (exists("norm")) {
                rats <- norm
            }
            else {
                if (org == "mm") {
                  rats <- mm_all_intensity.norm
                  rm(mm_all_intensity.norm)
                }
                else if (org == "ss") {
                  rats <- ss_all_intensity.norm
                  rm(ss_all_intensity.norm)
                }
                else if (org == "pf") {
                  if (exists("pf_new_all_intensity.norm")) {
                    pf_all_intensity.norm <- pf_new_all_intensity.norm
                    rm(pf_new_all_intensity.norm)
                  }
                  rats <- pf_all_intensity.norm
                  rm(pf_all_intensity.norm)
                }
            }
            print(dim(rats))
            tmp <- as.data.frame(rats)
            tmp <- cbind(PROBE = rownames(rats), tmp)
            probe.info.file <- list.files(path = sprintf("./data/%s/", 
                org), patt = glob2rx("*probe*info.txt.gz"), full = T)
            if (length(probe.info.file) <= 0) {
                if (org == "mm") 
                  probe.info <- read.delim(gzfile("data/mm/mm_probe_023308.info.txt.gz"), 
                    head = T)
                else if (org == "ss") 
                  probe.info <- read.delim(gzfile("data/ss/ss_probe_023310.info.txt.gz"), 
                    head = T)
                else if (org == "pf") 
                  probe.info <- read.delim(gzfile("data/pf/pf_probe_023309.info.txt.gz"), 
                    head = T)
            }
            else {
                probe.info <- read.delim(gzfile(probe.info.file), 
                  head = T)
            }
            rats <- merge(probe.info, tmp, by = "PROBE")
            rm(tmp, probe.info)
            colnames(rats)[colnames(rats) == "GENOME_ACCN"] <- "SEQ_ID"
            colnames(rats)[colnames(rats) == "ORIENTATION"] <- "GENE_EXPR_OPTION"
            colnames(rats)[colnames(rats) == "LENGTH"] <- "PROBE_LENGTH"
            colnames(rats)[colnames(rats) == "PROBE"] <- "PROBE_ID"
            tmp <- as.character(rats$GENE_EXPR_OPTION)
            tmp[tmp == "FORWARD"] <- "reverse"
            tmp[tmp == "REVERSE"] <- "forward"
            tmp <- toupper(tmp)
            rats$GENE_EXPR_OPTION <- as.factor(tmp)
            rm(tmp)
            rats$POSITION <- as.integer(as.character(rats$POSITION))
            rats$START <- as.integer(as.character(rats$START))
            rats$END <- as.integer(as.character(rats$END))
            if (length(grep(".ref", colnames(rats), fixed = T)) > 
                0) {
                refs <- rats[, c(1:8, grep(".ref", colnames(rats), 
                  fixed = T))]
                rats <- rats[, grep(".ref", colnames(rats), fixed = T, 
                  invert = T)]
                chr.map <- c(Chr = unique(as.character(rats$SEQ_ID)))
            }
            else {
                if (org == "mm") {
                  refs <- rats[, c(1:8, grep("^mma_R", colnames(rats)))]
                  rats <- rats[, grep("^mma_R", colnames(rats), 
                    invert = T)]
                  chr.map <- c(Chr = unique(as.character(rats$SEQ_ID)))
                }
                else if (org == "ss") {
                  refs <- rats[, c(1:8, grep("^Solfataricus_R", 
                    colnames(rats)))]
                  rats <- rats[, grep("^Solfataricus_R", colnames(rats), 
                    invert = T)]
                  chr.map <- unique(as.character(rats$SEQ_ID))
                  names(chr.map) <- c("Chr", "SSV1")
                }
                else if (org == "pf") {
                  refs <- rats[, c(1:8, grep("^p.fu_R", colnames(rats)))]
                  rats <- rats[, grep("^p.fu_R", colnames(rats), 
                    invert = T)]
                  chr.map <- c(Chr = unique(as.character(rats$SEQ_ID)))
                }
            }
            rownames(refs) <- rownames(rats) <- as.character(refs$PROBE_ID)
            gene.info.file <- list.files(path = sprintf("./data/%s/", 
                org), patt = glob2rx("*gene.info.txt.gz"), full = T)
            gene.coords <- read.delim(gzfile(gene.info.file), 
                head = T)
            gc()
            if (combine.dye.swap) {
                if (length(grep(".ref", colnames(refs), fixed = T)) > 
                  0) {
                  rat.conds <- colnames(rats)[9:ncol(rats)]
                  ref.conds <- colnames(refs)[9:ncol(refs)]
                  new.rats <- rats[, 1:8]
                  new.refs <- refs[, 1:8]
                  new.rat.conds <- new.ref.conds <- NULL
                  done.rat.cond <- rep(FALSE, length(rat.conds))
                  names(done.rat.cond) <- rat.conds
                  done.ref.cond <- rep(FALSE, length(ref.conds))
                  names(done.ref.cond) <- ref.conds
                  for (i in rat.conds) {
                    if (done.rat.cond[i]) 
                      next
                    conds <- grep(i, rat.conds, fixed = T, val = T)
                    if (length(conds) > 1) {
                      cat(conds, "\n")
                      tmp <- apply(rats[, conds], 1, sum, na.rm = T)
                      new.rats <- cbind(new.rats, tmp)
                      new.rat.conds <- c(new.rat.conds, i)
                      done.rat.cond[conds] <- TRUE
                      conds <- grep(i, ref.conds, fixed = T, 
                        val = T)
                      if (length(conds) > 1) {
                        cat("\t ==>", conds, "\n")
                        tmp <- apply(refs[, conds], 1, sum, na.rm = T)
                        new.refs <- cbind(new.refs, tmp)
                        new.ref.conds <- c(new.ref.conds, i)
                        done.ref.cond[conds] <- TRUE
                      }
                    }
                  }
                  colnames(new.rats)[9:ncol(new.rats)] <- new.rat.conds
                  rats <- new.rats
                  colnames(new.refs)[9:ncol(new.rats)] <- new.ref.conds
                  refs <- new.refs
                }
                else {
                  inds <- seq(9, ncol(rats), by = 2)
                  tmp <- (as.matrix(rats[, inds]) + as.matrix(rats[, 
                    (inds + 1)]))
                  colnames(tmp) <- colnames(rats)[inds]
                  rats <- cbind(rats[, 1:8], as.data.frame(tmp))
                  tmp <- (as.matrix(refs[, inds]) + as.matrix(refs[, 
                    (inds + 1)]))
                  colnames(tmp) <- colnames(refs)[inds]
                  refs <- cbind(refs[, 1:8], as.data.frame(tmp))
                  rm(tmp, inds)
                }
            }
            rats[, 9:ncol(rats)] <- log2(rats[, 9:ncol(rats)]/refs[, 
                9:ncol(refs)])
            refs[, 9:ncol(rats)] <- log2(refs[, 9:ncol(rats)])
        }
        rats <<- rats
        refs <<- refs
        gene.coords <<- gene.coords
        chr.map <<- chr.map
        cat("Data size:", dim(rats), dim(refs), "\n")
        tmp.rats <- rats[get.probes.in.window(chr = "Chr", dir = "FORWARD", 
            low = 0, upp = 9e+09, order = T), ]
        probe.spacing <<- as.numeric(names(which.max(table(abs(diff(tmp.rats$POSITION))))))
        cat("Probe spacing:", probe.spacing, "\n")
        rm(tmp.rats)
        if (plot) 
            plot.everything()
        gc()
    }
    if (!exists("reg.fit")) {
        reg.fit <<- regress.probe.expression.probs(slop = 0, 
            ...)
        if (min(reg.fit$pred, na.rm = T) > 0) {
            cat("Warning - regression did not set probs to zero (", 
                min(reg.fit$pred, na.rm = T), ") -- rescaling...\n")
            tmp <- reg.fit$pred - min(reg.fit$pred, na.rm = T)
            tmp <- tmp/max(tmp, na.rm = T)
            reg.fit$pred <<- tmp
            rm(tmp)
        }
        plot.everything()
    }
    if (!exists("break.sites")) {
        break.sites <<- joint.breakpoints(return.deltas = T, 
            n.boot = 20, weights.add = c(ref = 2, rat = 1, cor = 2, 
                pred = 0.25), ...)
        plot.everything()
    }
    if (!exists("mclapply")) 
        mclapply <- lapply
    cat(sum(sapply(break.sites$breaks, sapply, length)), "break points\n")
    if (is.null(break.sites$break.points)) {
        break.sites$break.points <<- mclapply(break.sites$breaks, 
            lapply, find.peaks.in.density, bw = probe.spacing, 
            cutoff = 0.2)
        plot.everything()
    }
    cat(sum(sapply(break.sites$break.points, sapply, length)), 
        "break sites\n")
    if (is.null(break.sites$prob.start.stop)) {
        if (!exists("probs.expressed")) 
            probs.expressed <<- get.orf.prob.expressed()
        use <- c("pred", "refs")
        if (org %in% c("halo", "mm")) 
            use <- c(use, "cors")
        break.sites$prob.start.stop <<- classify.start.vs.stop.sites(break.sites$breaks, 
            break.sites$deltas, use = use)
    }
    dd <- density(unlist(break.sites$prob.start.stop), n = 512)
    pdf(paste("dd_", org, ".pdf", sep = ""))
    plot(dd)
    dev.off()
    cutoff <- dd$x[128:356][which(abs(dd$y[128:356] - dd$y[129:357]) < 
        0.001)]
    if (length(cutoff) > 1) 
        cutoff <- cutoff[which.min(abs(cutoff - 0.5))]
    else if (length(cutoff) <= 0) 
        cutoff <- dd$x[128:356][which.min(abs(dd$y[128:356] - 
            dd$y[129:357]))]
    cat("Prob start/stop cutoff:", cutoff, "\n")
    if (!exists("break.sites.start")) {
        break.sites.start <<- mclapply(break.sites$prob.start.stop, 
            lapply, function(i) as.integer(names(i[i >= cutoff])))
        break.sites.start$break.points <<- mclapply(break.sites.start, 
            lapply, find.peaks.in.density, bw = probe.spacing, 
            cutoff = 0.2)
    }
    cat(sum(sapply(break.sites.start$break.points, sapply, length)), 
        "start break sites\n")
    if (!exists("break.sites.stop")) {
        break.sites.stop <<- mclapply(break.sites$prob.start.stop, 
            lapply, function(i) as.integer(names(i[i < cutoff])))
        break.sites.stop$break.points <<- mclapply(break.sites.stop, 
            lapply, find.peaks.in.density, bw = probe.spacing, 
            cutoff = 0.2)
    }
    cat(sum(sapply(break.sites.stop$break.points, sapply, length)), 
        "stop break sites\n")
    plot.everything()
    save.image(paste("tilingArraySeg_", org, ".RData", sep = ""))
    if (export) 
        export.data(org)
}
