# test script for performNormalization.R - testcases are NOT comprehensive!

test_that("performNormalization works", {
  
  seuratObj <- getdata("performPCA", "pbmc_hallmarks")
  GS.hallmark <- getdata("performNormalization", "GS.Hallmark")
  
  
  seuratObj.p <- performNormalization(seuratObj, 
                                      assay = "escape.ssGSEA", 
                                      gene.sets = GS.hallmark,
                                      make.positive = TRUE)
  seuratObj.pg <- performNormalization(seuratObj, 
                                      assay = "escape.ssGSEA", 
                                      gene.sets = GS.hallmark,
                                      make.positive = TRUE, groups=20)
  
  expect_equal(seuratObj.p@assays$escape.ssGSEA, 
               seuratObj.pg@assays$escape.ssGSEA)
  
  seuratObj.n <- performNormalization(seuratObj, 
                                      assay = "escape.ssGSEA", 
                                      gene.sets = GS.hallmark,
                                      make.positive = FALSE)
  seuratObj.ng <- performNormalization(seuratObj, 
                                      assay = "escape.ssGSEA", 
                                      gene.sets = GS.hallmark,
                                      make.positive = FALSE, groups=20)
  
  expect_equal(seuratObj.n@assays$escape.ssGSEA_normalized, 
               getdata("performNormalization", "performNormalization_nonpositive"))
  expect_equal(seuratObj.n@assays$escape.ssGSEA_normalized, 
               seuratObj.ng@assays$escape.ssGSEA_normalized)
  
})
