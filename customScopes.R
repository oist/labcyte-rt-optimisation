msScope_qcSI <- function(libs) {
  libs$Tag_dust     <- libs$extracted   - libs$rdna - libs$spikes - libs$cleaned
  libs$rDNA         <- libs$rdna
  libs$Spikes       <- libs$spikes
  libs$Unmapped     <- libs$cleaned     - libs$mapped
  libs$Non_proper   <- libs$mapped      - libs$properpairs
  libs$Duplicates   <- libs$properpairs - libs$librarySizes - libs$strandInvaders
  libs$Invaders     <- libs$strandInvaders
  libs$Counts       <- libs$librarySizes
  list( libs    = libs
      , columns = c( "Tag_dust", "rDNA", "Spikes", "Unmapped"
                   , "Non_proper", "Duplicates", "Invaders", "Counts")
      , total   = libs$extracted)
}

msScope_counts <- function(libs) {
  libs$Promoter   <- libs$promoter
  libs$Exon       <- libs$exon
  libs$Intron     <- libs$intron
  libs$Intergenic <- libs$librarySizes - libs$promoter - libs$intron - libs$exon
  libs$Invaders   <- libs$strandInvaders
  list( libs    = libs
      , columns = c("Promoter", "Exon", "Intron", "Intergenic", "Invaders")
      , total   = libs$librarySizes + libs$strandInvaders)
}

msScope_libSizeNormByBarcode <- function(libs) {
  libs$Yield   <- libs$libSizeNormByBarcode
  list( libs    = libs
      , columns = c("Yield")
      , total   = libs$Yield)
}

msScope_libSizeNormByIndex <- function(libs) {
  libs$Yield   <- libs$libSizeNormByIndex
  list( libs    = libs
      , columns = c("Yield")
      , total   = libs$Yield)
}

msScope_libSize <- function(libs) {
  libs$Yield   <- libs$librarySizes
  list( libs    = libs
      , columns = c("Yield")
      , total   = libs$Yield)
}
