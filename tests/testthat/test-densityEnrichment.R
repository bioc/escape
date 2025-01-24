# test script for densityEnrichment.R - testcases are NOT comprehensive!

test_that("densityEnrichment works", {
  
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
             Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
  set.seed(42)
  expect_doppelganger(
    "denistyEnrichment_default_plot",
    densityEnrichment(
      seuratObj, 
      gene.set.use = "Tcells",
      gene.sets = GS)
  )
  
  set.seed(42)
  expect_doppelganger(
    "denistyEnrichment_group.by_plot",
    densityEnrichment(
      seuratObj, 
      gene.set.use = "Bcells",
      gene.sets = GS,
      group.by = "groups")
  )
  
})
