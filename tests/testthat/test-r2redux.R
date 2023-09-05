test_that("r2_diff works",
          {
            dat=dat1
            expected_output=c(0.03896678,0.03836254,0.0001473686,0.0001436128,2.321425e-06,0.4278911,0.4366883,0.2183442,0.0006042383,0.00488788,-0.0005576171)
            nv=length(dat$V1)
            v1=c(1,2)
            v2=c(1)
            expect_equal(as.numeric(r2_diff(dat1,v1,v2,nv)),expected_output,tolerance=1e-7)
          })


test_that("r2_beta_var works",
          {
            dat=dat2
            expected_output=c(0.01118301,0.004980285,7.072931e-05,3.161929e-05,0.000162113,-2.988221e-05,0.03037793,-0.00123582,0.02490076,-0.005127546)
            nv=length(dat$V1)
            v1=c(1)
            v2=c(2)
            expect_equal(as.numeric(r2_beta_var(dat,v1,v2,nv)),expected_output,tolerance=1e-7)
          })



test_that("r2_enrich_beta works",
          {
            dat=dat2
            expected_output=c(0.01118301,0.004980285,0.4392572,0.1956205,0.08042288,0.0431134,0.9950922,-0.1165778,0.6025904,-0.2113493,0.1591692,0.07958459,0.000232035,0.0001160175)
            nv=length(dat$V1)
            v1=c(1)
            v2=c(2)
            expected_ratio=0.04
            expect_equal(as.numeric(r2_enrich_beta(dat,v1,v2,nv,expected_ratio)),expected_output,tolerance=1e-7)
          })



test_that("r2_var works",
          {
            dat=dat1
            expected_output=c(0.0001436128,3.98997e-10 , 1.188162e-10,0.03836254,0.06433782,0.01764252)
            nv=length(dat$V1)
            v1=c(1)
            expect_equal(as.numeric(r2_var(dat,v1,nv)),expected_output,tolerance=1e-7)
          })


test_that("r_diff works",
          {
            dat=dat1
            expected_output=c(0.1958636,0.197006,0.0009247466,0.0009238836,3.65286e-06,0.5500319,0.2750159,-0.001142375,0.002603666,-0.004888417)
            nv=length(dat$V1)
            v1=c(1)
            v2=c(2)
            expect_equal(as.numeric(r_diff(dat1,v1,v2,nv)),expected_output,tolerance=1e-7)
          })


