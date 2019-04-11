 plot.logregtree2=
function (x, nms, full = TRUE, and.or.cx = 1, leaf.sz = 1, leaf.txt.cx = 1, 
    coef.cx = 1, indents = rep(0, 4), coef = TRUE, coef.rd = 4, nny=0.25,
    ...) 
{
    ltree <- x
    if (class(ltree) != "logregtree") 
        stop("ltree is not of class logregtree")
    ntchx1 <- indents[2]
    ntchy1 <- indents[1]
    ntchx2 <- indents[4]
    ntchy2 <- indents[3]
    logregtmp <- hist(1:10, plot = FALSE)
    if (length(logregtmp$mids) > 0) 
        logregwhite <- 1
    if (length(logregtmp$mids) == 0) 
        logregwhite <- 1
    data <- ltree$trees
    cf <- ltree$coef
    if (is.na(cf)) 
        coef <- FALSE
    names(data) <- c("number", "conc", "knot", "neg", "pick")
    if (data$pick[1] == 0) {
        cat("empty tree", "\n")
        plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", 
            ylab = "")
        text(0.5, 0.5, "Empty tree", cex = 1.5 * coef.cx)
    }
    else {
        if (!missing(nms)) {
            which <- (1:dim(data)[1])[data$conc == 3]
            data$knot[which] <- nms[data$knot[which]]
        }
        level <- ceiling(log(1 + max(data$number[data$pick >= 
            1]))/log(2))
        max.p <- 2^level - 1
        data <- data[1:max.p, ]
        data$x <- rep(0, max.p)
        data$y <- rep(0, max.p)
        left.b <- 2^(0:(level - 1))
        right.b <- 2^(1:level) - 1
        l.max <- max(data$number[left.b][data$pick[left.b] >= 
            1])
        r.max <- max(data$number[right.b][data$pick[right.b] >= 
            1])
        if (full) {
            l.shift <- 1/(l.max + 1)
            r.shift <- 1/(r.max + 1)
            w.shift <- min(l.shift, r.shift)
            plot(c(0.85 * (l.shift - ntchx1), 0.95 * (1 - r.shift + 
                ntchx2)), c(-1 + ntchy1, 1.05 * (-level - ntchy2)), 
                type = "n", xlab = "", ylab = "", xaxt = "n", 
                yaxt = "n", col = 1)
            if (coef) {
                text(0.85 * (l.shift - ntchx1), -1 + ntchy1, 
                  paste("Parameter =", round(cf, coef.rd)), cex = 1.5 * 
                    coef.cx, adj = 0)
            }
        }
        else {
            plot(c(0, 1), c(-1, -level), type = "n", xlab = "", 
                ylab = "", xaxt = "n", yaxt = "n", col = 1)
            text(0.02, -1.2, paste("Parameter =", round(cf, coef.rd)), 
                cex = 1.5 * coef.cx, adj = 0)
        }
        for (k in level:1) {
            y.val <- (-k)
            min.val <- 2^(-y.val - 1)
            max.val <- 2^(-y.val) - 1
            diff.x <- 1/(min.val + 1)
            x.val <- seq(diff.x, 1 - diff.x, diff.x)
            if (level == 1) 
                x.val <- x.val * (0.95 * (1 - r.shift + ntchx2) - 
                  0.85 * (l.shift - ntchx1)) + 0.85 * (l.shift - 
                  ntchx1)
            for (j in 1:length(x.val)) {
                which <- 2^(k - 1) + j - 1
                if (data$pick[which] == 1) {
                  data$x[which] <- x.val[j]
                  data$y[which] <- y.val
                  if (data$conc[which] == 1) {
                    text(x.val[j], y.val, "and", cex = 2 * and.or.cx)
                  }
                  if (data$conc[which] == 2) {
                    text(x.val[j], y.val, "or", cex = 2 * and.or.cx)
                  }
                  if (data$conc[which] == 3) {
                    if (data$neg[which] == 0) {
                      points(x.val[j], y.val, pch = 0, cex = 7 * 
                        leaf.sz)
                      text(x.val[j], y.val-nny, as.character(data$knot[which]), 
                        cex = 0.5 * leaf.txt.cx)
                    }
                    if (data$neg[which] == 1) {
                      points(x.val[j], y.val, pch = 15, cex = 7 * 
                        leaf.sz)
                      text(x.val[j], y.val-nny, as.character(data$knot[which]), 
                        col = logregwhite, cex = 0.5 * leaf.txt.cx)
                    }
                  }
                  if (data$conc[which] == 4) {
                    if (data$neg[which] == 0) {
                      points(x.val[j], y.val, pch = 1, cex = 7 * 
                        leaf.sz)
                      text(x.val[j], y.val-nny, as.character(data$knot[which]), 
                        cex = 0.5)
                    }
                    if (data$neg[which] == 1) {
                      points(x.val[j], y.val, pch = 16, cex = 7 * 
                        leaf.sz)
                      text(x.val[j], y.val-nny, as.character(data$knot[which]), 
                        cex = 0.5, col = logregwhite)
                    }
                  }
                }
            }
        }
        for (k in 1:max.p) {
            if ((ceiling(log(k)/log(2)) < level) & (2 * k < max.p)) {
                if ((data$pick[k] == 1) & (data$pick[2 * k] == 
                  1)) {
                  lines(data$x[c(k, 2 * k)], c(data$y[k] - 0.2, 
                    data$y[2 * k] + 0.2))
                  lines(data$x[c(k, 2 * k + 1)], c(data$y[k] - 
                    0.2, data$y[2 * k + 1] + 0.2))
                }
            }
        }
        invisible()
    }
}


plot.logregmodel2=
function (x, pscript = FALSE, title = TRUE, nms, ...) 
{
    if (class(x) != "logregmodel") 
        stop("x not of class logregmodel")
    logregtmp <- hist(1:10, plot = FALSE)
    if (length(logregtmp$mids) > 0) 
        logregwhite <- 1
    if (length(logregtmp$mids) == 0) 
        logregwhite <- 1
    ntr <- x$ntrees[1]
    msz <- x$nleaves[1]
    if (msz == -1) 
        msz <- 0
    scores <- x$coef
    lscores <- length(scores)
    scores <- scores[(lscores - ntr + 1):lscores]
    for (j in 1:ntr) {
        if (pscript) 
            postscript(paste("tree", ntr, ".", msz, ".", j, ".ps", 
                sep = ""), print.it = FALSE, horizontal = TRUE, 
                ...)
        if (missing(nms)) 
            plot.logregtree2(x$trees[[j]], indents = c(0.2, 0, 
                0, 0), coef = scores[j], ...)
        else plot.logregtree2(x$trees[[j]], indents = c(0.2, 0, 
            0, 0), coef = scores[j], nms = nms, ...)
        if (title && msz > 0) 
 #           title(main = paste("tree", j, "out of", ntr, "total size is", 
   #             msz))
        if (title && msz <= 0) 
#            title(main = paste("tree", j, "out of", ntr))
        if (pscript) 
            dev.off()
    }
}