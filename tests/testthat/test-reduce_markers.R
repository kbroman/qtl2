context("reduce markers")

test_that("reduce_markers matches qtl::pickMarkerSubset", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    map <- grav2$gmap

    wts <- lapply(vapply(map, length, 1),
                  function(n) runif(n, 1, 5))

    submap <- reduce_markers(map, 1, wts)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    for(i in seq(along=map)) {
        submapi <- reduce_markers(map[[i]], 1, wts[[i]])
        expect_equal(submapi, submap[[i]])

        expect_equal(names(submapi), qtl::pickMarkerSubset(map[[i]], 1, wts[[i]]))
    }

    # repeat with 5 cM
    submap <- reduce_markers(map, 5, wts)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    for(i in seq(along=map)) {
        submapi <- reduce_markers(map[[i]], 5, wts[[i]])
        expect_equal(submapi, submap[[i]])

        expect_equal(names(submapi), qtl::pickMarkerSubset(map[[i]], 5, wts[[i]]))
    }

    # without weights, get same result
    set.seed(8764972)
    submap <- reduce_markers(map, 1)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    set.seed(8764972)
    for(i in seq(along=map)) {
        submapi <- reduce_markers(map[[i]], 1)
        expect_equal(submapi, submap[[i]])
    }

    # repeat at 5 cM
    set.seed(8764972)
    submap <- reduce_markers(map, 5)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    set.seed(8764972)
    for(i in seq(along=map)) {
        submapi <- reduce_markers(map[[i]], 5)
        expect_equal(submapi, submap[[i]])
    }

})


test_that("reduce_markers works", {

    map <- list(c1=c(marker=5),
                c2=c(marker1=5, marker2=7))
    wts <- list(1, c(1,2))
    expect_equal(reduce_markers(map), map)
    expect_equal(reduce_markers(map, 3, wts),
                 list(c1=c(marker=5), c2=c(marker2=7)))

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))

    RNGkind("Mersenne-Twister")
    set.seed(20151205)
    # use proportion typed as weights, but randomize a bit so there are no ties
    prtyped <- lapply(grav2$geno, function(a) colMeans(!is.na(a) & a>0) + runif(ncol(a), 0, 0.1))

    # expected results at min_dist=0.5
    expected <-  structure(list(`1` = c(PVV4 = 0, `AXR-1` = 6.250674, `HH.335C-Col/PhyA` = 9.303868,
                                        EC.480C = 12.577629, EC.66C = 18.39283, GD.86L = 19.648542, `CH.160L-Col` = 25.23229,
                                        `CC.98L-Col/101C` = 25.882095, AD.121C = 35.669554, `AD.106L-Col` = 36.631747,
                                        GB.112L = 44.511901, GD.97L = 50.602339, `EG.113L/115C` = 53.631519,
                                        CD.89C = 57.405778, `BF.206L-Col` = 62.323561, CH.200C = 67.24062,
                                        EC.88C = 72.296762, GD.160C = 76.255226, HH.375L = 81.273196,
                                        CH.215L = 81.900811, BF.116C = 88.639045, `GH.157L-Col` = 90.914789,
                                        CC.318C = 97.121385, `CD.173L/175C-Col` = 99.711975, `GH.127L-Col/ADH` = 104.80612,
                                        `HH.360L-Col` = 109.519692),
                                `2` = c(AD.156C = 0, BF.325L = 6.241796, GH.580L = 11.678237, DF.225L = 12.929907, BH.145C = 13.543347,
                                        FD.226C = 14.160665, `EC.495C-Col` = 14.777988, FD.81L = 16.350906,
                                        BF.105C = 17.282463, CH.284C = 19.850147, `FD.222L-Col` = 23.453866,
                                        CD.245L = 25.685576, CH.65C = 27.258244, CH.1500C = 30.511191,
                                        BF.221L = 31.128921, FD.85C = 35.085099, `GB.150L-Col` = 37.651404,
                                        FD.150C = 39.551472, `GD.460L-Col` = 41.451601, Erecta = 42.701883,
                                        `BH.195L-Col` = 47.01718, GD.298C = 47.638238, GH.247L = 49.221621,
                                        `BH.120L-Col` = 52.14781, DF.140C = 55.429739, `EG.357C/359L-Col` = 62.530946),
                                `3` = c(DF.77C = 0, `GB.120C-Col` = 2.613844, `GD.248C-Col/249L` = 7.032127,
                                        CH.322C = 8.291564, `FD.111L-Col/136C` = 11.480925, CC.266L = 15.051063,
                                        `BF.270L-Col/271C` = 16.346984, `CC.110L/127C` = 18.094189, EC.58C = 18.776988,
                                        `GH.321L/323C-Col` = 21.422012, `GH.226C/227L-Col` = 22.0854,
                                        `EC.83C/84L` = 25.665317, `GD.318C/320L` = 28.194604, AD.92L = 30.995068,
                                        HH.410C = 34.67091, CD.800C = 38.086009, `BF.134C-Col` = 39.704548,
                                        EG.122C = 43.418084, `GD.113C-Col` = 44.350416, `BH.323C-Col` = 45.169461,
                                        DF.303C = 48.653206, DF.76L = 49.640633, BF.128C = 50.30031,
                                        HH.117C = 52.688081, `GD.296C-Col` = 55.438502, DF.328C = 58.346106,
                                        FD.98C = 59.596912, AD.182C = 66.588704, GD.106C = 68.488354,
                                        `HH.171C-Col/173L` = 70.695509, `AD.495L-Col` = 71.31314, `AD.112L-Col` = 72.897197,
                                        `BH.109L-Col` = 75.144036, `HH.90L-Col` = 76.392762),
                                `4` = c(ANL2 = 0, GH.250C = 6.388021, GH.165L = 17.889887, BF.151L = 18.5072, EC.306L = 22.755474,
                                        `BH.92L-Col` = 27.797576, FD.154L = 30.442727, CH.238C = 34.758259,
                                        `CD.84C-Col/85L` = 35.37558, SC5 = 39.942822, `DF.108L-Col` = 40.635774,
                                        AD.307C = 45.254346, `AD.115L-Col` = 46.537827, `FD.167L-Col` = 48.152558,
                                        `CH.70L/71C-Col` = 53.206339, `GH.433L-Col` = 56.459811, GB.490C = 60.771852,
                                        GB.750C = 67.771693, `BH.342C/347L-Col` = 73.668364),
                                `5` = c(FD.207L = 0, CH.690C = 2.94061, AD.292L = 5.263751, BH.144L = 6.232987, `EC.198L-Col` = 11.517411,
                                        BH.180C = 12.450901, BH.325L = 13.120403, `BH.107L-Col` = 17.456274,
                                        BF.269C = 18.082939, `AD.114C-Col` = 21.390456, DF.231C = 23.320809,
                                        `DF.184L-Col` = 28.452585, GH.473C = 33.209409, GH.117C = 34.154297,
                                        `GH.121L-Col` = 38.533295, DF.154C = 41.782964, HH.480C = 45.506147,
                                        `BH.96L-Col` = 51.937046, CC.188L = 53.180697, EC.395C = 56.869536,
                                        DF.300C = 57.82089, `GB.235C-Col` = 59.761174, CD.179L = 61.70335,
                                        CH.88L = 63.311003, CC.216C = 65.25314, EC.151L = 69.302894,
                                        DFR = 70.2914, AD.254C = 76.061598, `AD.75C-Col` = 76.680571,
                                        `GD.350L-Col` = 77.933341, CD.160L = 78.551316, GB.223C = 79.168573,
                                        DF.450C = 80.418873, `HH.445L-Col` = 81.991829, GD.118C = 82.92342,
                                        CC.262C = 85.155018, `GB.102L-Col/105C` = 87.390432, HH.143C = 93.991816,
                                        `BF.168L-Col` = 95.563679, DF.119L = 101.176465, `CH.331L-Col` = 102.460924,
                                        `GD.222C-Col` = 103.094847, g2368 = 103.932196, EG.205L = 107.488997)),
                           is_x_chr = c(`1` = FALSE, `2` = FALSE, `3` = FALSE, `4` = FALSE, `5` = FALSE))

    attr(expected, "is_x_chr") <- attr(grav2$gmap, "is_x_chr")

    result <- reduce_markers(grav2$gmap, 0.5, prtyped)

    expect_equal(result, expected)

})
