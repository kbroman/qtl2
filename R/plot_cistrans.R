#' cis-trans plot for eQTL results
#'
#' Scatterplot of gene location vs QTL location to display eQTL/pQTL results
#'
#' @param peaks Data frame of QTL results, as output by [find_peaks()]
#'     Should contain columns `lodcolumn`, `chr`, and `pos`.
#'
#' @param map Marker map, as a list of chromosomes, each being a
#'     vector of positions; used to get chromosome names and start/end
#'
#' @param pheno_pos Data frame with chromosome, position, and name of
#'     the phenotypes in `peaks`. Should contain columns `pheno`,
#'     `chr`, and `pos`.
#'
#' @param gap Gap between chromosomes in the plot
#'
#' @param pty Plot type; the default "s" forces a square plot; use "m"
#'     to use the maximal plotting region.
#'
#' @param bgcolor Background color for chromosome rectangles
#'
#' @param altbgcolor Alternate background color for chromosome rectangles
#'
#' @param pch Point type for "trans" points
#'
#' @param col Point color for "trans" points
#'
#' @param bg Background color for "trans" points
#'
#' @param cex Character expansion size for "trans" points
#'
#' @param pch_cis Point type for "cis" points
#'
#' @param col_cis Point color for "cis" points
#'
#' @param bg_cis Background color for "cis" points
#'
#' @param cex_cis Character expansion size for "cis" points
#'
#' @param cis_window Window size that defines cis-QTL (on same
#'     chromosome, and position of gene within this distance of QTL
#'     position)
#'
#' @param ... Additional graphics arguments passed to [graphics::points()]
#'
#' @return Invisibly returns a data frame with `peaks` merged with
#'     `pheno_pos`, with an additional column `cis` that indicates
#'     which QTL appear to be cis (vs trans).
#'
#' @seealso [plot_scan1()], [find_peaks()], [plot_lodpeaks()]
#'
#' @importFrom graphics rect par axis title points
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example pQTL results (from Keele et al. 2026, https://doi.org/10.1016/j.xgen.2025.101069)
#' # contains qtl, map_endpoints, and pheno_pos
#' url <- "https://kbroman.org/qtl2/assets/sampledata/pqtl_data.RData"
#' tempfile <- file.path(tempdir(), basename(url))
#' download.file(url, tempfile)
#' load(tempfile)
#' unlink(tempfile)
#'
#' plot_cistrans(qtl, map, pheno_pos, cis_window=5)
#' }

plot_cistrans <-
    function(peaks, map, pheno_pos, gap=0, pty="s",
             bgcolor="gray90", altbgcolor="gray80",
             pch=21, col="slateblue", bg="slateblue", cex=0.8,
             pch_cis=21, col_cis="violetred", bg_cis="violetred", cex_cis=0.8,
             cis_window=1, ...)
{
    # check that peaks conforms to expectations
    if(!is.data.frame(peaks) || !all(c("lodcolumn", "chr", "pos") %in% colnames(peaks)))
        stop("peaks should be a data frame with columns lodcolumn, chr, and pos")
    # check that map conforms to expectations
    stopifnot(is.list(map))
    # check that pheno_pos conforms to expectations
    if(!is.data.frame(pheno_pos) || !all(c("pheno", "chr", "pos") %in% colnames(pheno_pos)))
        stop("peaks should be a data frame with columns pheno, chr, and pos")
    # check that the phenotypes in peaks are also in pheno_pos
    if(!all(peaks$lodcolumn %in% pheno_pos$pheno))
        stop("Not all lodcolumns in peaks are found in the pheno column of pheno_pos")

    # merge peaks and pheno_pos
    peaks <- base::merge(peaks, pheno_pos, by.x="lodcolumn", by.y="pheno")
    # drop loci that don't map to the provided map
    peaks <- peaks[peaks$chr.x %in% names(map) & peaks$chr.y %in% names(map), ]
    # should give xpos_scan1 a "strict" argument, and make it NA if outside of chr's range
    peaks$xpos <- xpos_scan1(map, gap=gap, thechr=peaks$chr.x, thepos=peaks$pos.x)
    peaks$ypos <- xpos_scan1(map, gap=gap, thechr=peaks$chr.y, thepos=peaks$pos.y)

    plot_cistrans_internal <-
        function(peaks, map, pheno_pos, gap=0, pty="s", xlab="QTL position",
                 ylab="Gene position", xlim=NULL, ylim=NULL,
                 mgp=NULL, mgp.x=NULL, mgp.y=NULL,
                 bgcolor="gray90", altbgcolor="gray80",
                 col=col, bg = bg, pch=pch, las=1, cex=cex,
                 type="p", ...)
    {
        # get x- and y-axis range
        if(is.null(xlim) || is.null(ylim)) {
            start <- xpos_scan1(map, chr=names(map), gap=gap,
                                names(map)[1], min(map[[1]]))
            end <- xpos_scan1(map, chr=names(map), gap=gap,
                              names(map)[length(map)], max(map[[length(map)]]))

            if(is.null(xlim)) xlim <- c(start, end)
            if(is.null(ylim)) ylim <- c(start, end)
        }

        if(is.null(mgp.x)) {
            if(!is.null(mgp)) mgp.x <- mgp
            else mgp.x <- c(1.7, 0.3, 0)
        }
        if(is.null(mgp.y)) {
            if(!is.null(mgp)) mgp.y <- mgp
            else mgp.y <- c(2.0, 0.4, 0)
        }


        par(pty=pty)
        plot(0,0, type="n", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i",
             xlab="", ylab="", xaxt="n", yaxt="n")
        chr <- names(map)
        chr_midpt <- sapply(seq_along(map), function(i) xpos_scan1(map, names(map), gap=gap,
                                                                   names(map)[i], mean(range(map[[i]], na.rm=TRUE))))

        for(i in seq_along(chr_midpt)) {
            graphics::axis(side=1, at=chr_midpt[i], chr[i], mgp=mgp.x, tick=FALSE, las=las)
            graphics::axis(side=2, at=chr_midpt[i], chr[i], mgp=mgp.y, tick=FALSE, las=las)
        }
        graphics::title(xlab=xlab, mgp=mgp.x)
        graphics::title(ylab=ylab, mgp=mgp.y)

        # add background rectangles
        chrstart <- xpos_scan1(map, gap=gap, thechr=names(map), thepos=sapply(map, min))
        chrend <- xpos_scan1(map, gap=gap, thechr=names(map), thepos=sapply(map, max))

        rect(min(chrstart), min(chrstart), max(chrend), max(chrend), col=bgcolor, border=NA)
        for(chri in seq_along(map)) {
            for(chrj in seq_along(map)) {
                if((chri+chrj) %% 2 == 1) {
                    rect(chrstart[chri], chrstart[chrj], chrend[chri], chrend[chrj],
                         col=altbgcolor, border=NA)
                }
            }
        }

        # plot points
        cis <- (!is.na(peaks$chr.x) & !is.na(peaks$chr.y) & !is.na(peaks$pos.x) & !is.na(peaks$pos.y) &
                peaks$chr.x == peaks$chr.y & abs(peaks$pos.x - peaks$pos.y) < cis_window)

        if(type != "n") {
            if(any(!cis)) points(peaks$xpos[!cis], peaks$ypos[!cis], pch=pch, bg=bg, col=col, cex=cex, ...)
            if(any(cis)) points(peaks$xpos[cis], peaks$ypos[cis], pch=pch_cis, bg=bg_cis, col=col_cis, cex=cex_cis, ...)
        }

        peaks$cis <- cis
        invisible(peaks)
    }

    plot_cistrans_internal(peaks, map, pheno_pos, gap=0, pty=pty,
                           bgcolor=bgcolor, altbgcolor=altbgcolor,
                           pch=pch, cex=cex, bg=bg, col=col,
                           ...)
}
