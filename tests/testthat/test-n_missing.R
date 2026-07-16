context("n_missing and n_typed")

test_that("n_missing and n_typed work for iron", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    expect_equal(n_missing(iron[1:4,]), c("1"=0,  "2"=36, "3"=0,  "4"=0))
    expect_equal(n_typed(iron[1:4,]),   c("1"=66, "2"=30, "3"=66, "4"=66))
    expect_equal(n_missing(iron[1:4,], sum="prop"), c("1"=0, "2"=18/33, "3"=0, "4"=0))
    expect_equal(n_typed(iron[1:4,], sum="prop"),   c("1"=1, "2"=15/33, "3"=1, "4"=1))

    expect_equal(n_missing(iron[,3], "mar"),         c("D3Mit22"=129,"D3Mit18"=129))
    expect_equal(n_typed(iron[,3], "mar"),           c("D3Mit22"=155,"D3Mit18"=155))
    expect_equal(n_missing(iron[,3], "mar", "prop"), c("D3Mit22"=129/284,"D3Mit18"=129/284))
    expect_equal(n_typed(iron[,3], "mar", "prop"),   c("D3Mit22"=155/284,"D3Mit18"=155/284))

})

test_that("n_missing and n_typed work for grav2", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))

    expect_equal(n_missing(grav2[125:128,]), c("125"=2,  "126"=3, "127"=4,  "128"=1))
    expect_equal(n_typed(grav2[125:128,]),   c("125"=232,  "126"=231, "127"=230,  "128"=233))
    expect_equal(n_missing(grav2[125:128,], sum="prop"), c("125"=2/234, "126"=3/234, "127"=4/234, "128"=1/234))
    expect_equal(n_typed(grav2[125:128,], sum="prop"), c("125"=232/234, "126"=231/234, "127"=230/234, "128"=233/234))

    expected <- c(AD.156C = 1, BF.325L = 2, GH.580L = 2, DF.225L = 0, AD.77L = 1,
                  CH.266C = 1, CH.610C = 0, HH.258C = 0, BH.145C = 0, `BF.226C/BH.58L` = 0,
                  FD.226C = 0, GD.145C = 0, GH.94L = 1, BF.82C = 0, GD.465C = 0,
                  FD.306L = 0, `EC.495C-Col` = 0, BH.460L = 3, FD.81L = 0, BF.105C = 0,
                  CH.284C = 1, `FD.222L-Col` = 1, CD.245L = 0, EG.66L = 1, CH.65C = 0,
                  CH.1500C = 4, BF.221L = 0, FD.85C = 0, `GB.150L-Col` = 0, FD.150C = 1,
                  `GD.460L-Col` = 0, CC.332C = 3, Erecta = 0, `CH.145L-Col/150C` = 2,
                  `AD.191L-Col` = 3, `BH.195L-Col` = 0, GD.298C = 1, GH.247L = 2,
                  `BH.120L-Col` = 1, DF.140C = 4, `EG.357C/359L-Col` = 2, `EC.235L-Col/247C` = 3)
    expect_equal(n_missing(grav2[,2], "mar"), expected)

})
