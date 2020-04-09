test_that("`drain_dem()` works", {
  library(magrittr)
  library(demdrainer)
  dem <- FedData::get_ned(demdrainer::mcphee, label = "McPhee")
  drained_dem <- drain_dem(dem, label = "McPhee")
  expect_true(!identical(dem, drained_dem))

})
