# test script for geyserEnrichment.R - testcases are NOT comprehensive!

test_that("geyserEnrichment works", {
  
  seuratObj <- getdata("runEscape", "pbmc_small_ssGSEA")
  
  set.seed(42)
  expect_doppelganger(
    "geyserEnrichment_default_plot",
    geyserEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Tcells"
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "geyserEnrichment_scale_plot",
    geyserEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Tcells",
      scale = TRUE
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "geyserEnrichment_facet_plot",
    geyserEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Tcells",
      facet.by = "groups"
    )
  )

  set.seed(42)
  expect_doppelganger(
    "geyserEnrichment_order_plot",
    geyserEnrichment(
      seuratObj, 
      order.by = "mean",
      assay = "escape",
      gene.set = "Tcells"
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "geyserEnrichment_gradient_plot",
    geyserEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Tcells",
      color.by = "Tcells"
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "geyserEnrichment_gradient_reorder_plot",
    geyserEnrichment(
      seuratObj, 
      assay = "escape",
      order.by = "mean",
      gene.set = "Tcells",
      color.by = "Tcells"
    )
  )
  
  set.seed(42)
  expect_doppelganger(
    "geyserEnrichment_gradient_facet_plot",
    geyserEnrichment(
      seuratObj, 
      assay = "escape",
      gene.set = "Tcells",
      color.by = "Tcells",
      facet.by = "groups"
    )
  )
  
})
