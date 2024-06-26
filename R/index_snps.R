#' Create index of equivalent SNPs
#'
#' For a set of SNPs and a map of marker/pseudomarkers, partition the
#' SNPs into groups that are contained within common intervals and
#' have the same strain distribution pattern, and then create an index
#' to a set of distinct SNPs, one per partition.
#'
#' @param map Physical map of markers and pseudomarkers; generally
#'     created from [insert_pseudomarkers()] and
#'     used for a set of genotype probabilities (calculated with
#'     [calc_genoprob()]) that are to be used to
#'     interpolate SNP genotype probabilities (with
#'     [genoprob_to_snpprob()]).
#' @param snpinfo Data frame with SNP information with the following columns:
#' * `chr` - Character string or factor with chromosome
#' * `pos` - Position (in same units as in the `"map"`).
#' * `sdp` - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' * `snp` - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' @param tol Tolerance for determining whether a SNP is exactly at a
#' position at which genotype probabilities were already calculated.
#'
#' @return A data frame containing the input `snpinfo` with three
#' added columns: `"index"` (which indicates the groups of
#' equivalent SNPs), `"interval"` (which indicates the map
#' interval containing the SNP, with values starting at 0), and
#' `on_map` (which indicates that the SNP is within
#' `tol` of a position on the map). The rows get reordered,
#' so that they are ordered by chromosome and position, and the
#' values in the `"index"` column are _by chromosome_.
#'
#' @details We split the SNPs by chromosome and identify the intervals
#' in the `map` that contain each. For SNPs within `tol`
#' of a position at which the genotype probabilities were
#' calculated, we take the SNP to be at that position. For each
#' marker position or interval, we then partition the SNPs into
#' groups that have distinct strain distribution patterns, and
#' choose a single index SNP for each partition.
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DO_Recla/recla.zip")
#' recla <- read_cross2(file)
#'
#' # founder genotypes for a set of SNPs
#' snpgeno <- rbind(m1=c(3,1,1,3,1,1,1,1),
#'                  m2=c(1,3,1,3,1,3,1,3),
#'                  m3=c(1,1,1,1,3,3,3,3),
#'                  m4=c(1,3,1,3,1,3,1,3))
#' sdp <- calc_sdp(snpgeno)
#' snpinfo <- data.frame(chr=c("19", "19", "X", "X"),
#'                       pos=c(40.36, 40.53, 110.91, 111.21),
#'                       sdp=sdp,
#'                       snp=c("m1", "m2", "m3", "m4"), stringsAsFactors=FALSE)
#'
#' # update snp info by adding the SNP index column
#' snpinfo <- index_snps(recla$pmap, snpinfo)
#' }
#'
#' @seealso [genoprob_to_snpprob()], [scan1snps()], [find_index_snp()]
#' @export
index_snps <-
    function(map, snpinfo, tol=1e-8)
{
    if(is.null(map)) stop("map is NULL")
    if(is.null(snpinfo)) stop("snpinfo is NULL")
    if(!is_nonneg_number(tol)) stop("tol should be a single non-negative number")

    if(nrow(snpinfo)==0) {
        # no rows; return data frame with no rows
        return( cbind(snpinfo, index=numeric(0), interval=numeric(0), on_map=logical(0)) )
    }

    uchr <- unique(snpinfo$chr)
    chrID <- names(map)
    if(!all(uchr %in% chrID)) {
        mischr <- uchr[!(uchr %in% chrID)]
        stop("Not all chr found in genoprobs: ", paste(mischr, collapse=","))
    }
    # reorder
    uchr <- factor(factor(uchr, levels=chrID))

    # reorder snp info by chromosome and position
    snpinfo <- snpinfo[order(factor(snpinfo$chr, levels=uchr), snpinfo$pos),,drop=FALSE]

    # if more than one chromosome:
    if(length(uchr) > 1) {
        uchr <- as.character(uchr) # convert back to character

        ### split snpinfo by chr, ordered as in genoprobs
        chr <- factor(snpinfo$chr, levels=uchr)
        snpinfo_spl <- split(snpinfo, chr)

        ### loop over chromosomes, using recursion
        for(i in seq_along(uchr))
            snpinfo_spl[[i]] <- index_snps(map, snpinfo_spl[[i]])

        # combine the results
        snpinfo <- snpinfo_spl[[1]]
        for(i in 2:length(uchr))
            snpinfo <- rbind(snpinfo, snpinfo_spl[[i]])

        return(snpinfo)
    }

    # one chromosome

    # make chromosome a character string again
    uchr <- as.character(uchr)

    ### find snps in map
    this_map <- map[[uchr]]
    if(any(is.na(this_map))) stop("Missing values in map on chr ", uchr)
    snploc <- find_intervals(snpinfo$pos, this_map, tol)
    interval <- snploc[,1]
    on_map <- (snploc[,2]==1)

    # drop snps outside of range
    snps2drop <- (interval < 0 | (interval >= length(this_map)-1 & !on_map))
    if(any(snps2drop)) {
        snpinfo <- snpinfo[!snps2drop,,drop=FALSE]
        interval <- interval[!snps2drop]
        on_map <- on_map[!snps2drop]
    }
    if(nrow(snpinfo) == 0)
        stop("No SNPs within range")

    ### find unique (interval, on_map, sdp) patterns
    pat <- paste(interval, on_map, snpinfo$sdp, sep=":")
    upat <- unique(pat)
    ### order unique patterns by position
    pos <- snpinfo$pos[match(upat, pat)]
    upat <- upat[order(pos)]
    snpinfo$index <- match(pat, upat)
    snpinfo$interval <- interval
    snpinfo$on_map <- on_map
    ### first row with each unique pattern
    urow <- seq_len(nrow(snpinfo))[match(upat, pat)]
    ### update index as linking to row number
    snpinfo$index <- urow[snpinfo$index]

    snpinfo
}
