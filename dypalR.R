#################################################################################
######################### Dypal Console Grid Simulations ########################
#################################################################################
# We empty any previous data
rm(list = ls())

require(rgl)
require(emdbook)

########################
kappa.F <- function(ref, comp){
  
  non.na <- !is.na(comp==ref)
  ref    <- as.numeric( ref[non.na])
  comp   <- as.numeric(comp[non.na])
  
  lev.r  <-  unique(ref)
  r      <- length(lev.r)
  k.mat  <- matrix(NA, nr=r, nc=r, dimnames = list(lev.r, lev.r))
  
  for (i in 1:r) {
    for (j in 1:r) {
      k.mat[i, j] <- sum(ref==lev.r[i] & comp==lev.r[j])}}
  
  global <- as.numeric(sum(k.mat))
  marg   <- sum(as.numeric(apply(k.mat, 1, sum)) * 
    as.numeric(apply(k.mat, 2, sum)))
  tr.mat <- sum(diag(k.mat))
  
  return((global * tr.mat - marg) / (global ^ 2  - marg))
}


# The functions below aim at launching a set of simulation using Dypal in console mode
# Dypal is called in dypal.compute() that writes a new batch through batch.write()
# upon lambda{i} values previously defined in a grid through make.grid().
# For each row of the grid made, we thus write the batch, launch Dypal 
# and then wait for the end of the simulation.
# The kappa value calculated by comparing simulated landscapes and a reference shape is then read
# and stored in a results_temp file. At the end of the simulation, the final results
# file is written in a dedicated folder and temporary files are removed.

# Manipulate and write the .xml batch
batch.write <- function(batch, new.lambdas=c(1, 1)) {
  # we first read the batch
  batch.rw  <- readLines(batch, warn=FALSE)
  # we split it using the '<OPERATION>' balises
  batch.rw  <- strsplit(batch.rw, "<OPERATION>")[[1]]
  # we check that the number of new lambdas provided
  # correspond to the number of operations
  if (length(new.lambdas) != (length(batch.rw)-1)){
    stop("The number of new lamdbas provided is not equal to
          the number of lambdas in the batch file")}
  # the first element initiates the building of the new batch
  batch.new <- batch.rw[1]
  for (i in 2:length(batch.rw)) {
    # testing below: if a raster-based operation is detected
    # the structure of values and their enclosing balises must be adapted accordingly
    if(grepl("<MATRIX_EQUATION>", batch.rw[i])) {
      batch.new <- paste(batch.new,
                         sub("</OPERATOR><VALUE>.*</VALUE>",
                             paste("</OPERATOR><VALUE>", new.lambdas[i-1], "</VALUE>", sep=""),
                             batch.rw[i]), sep="<OPERATION>")
    } else {  #ie homogen dilation
      batch.new <- paste(batch.new,
                         sub("<DILATION_SIZE>.*</DILATION_SIZE>",
                             paste("<DILATION_SIZE>", new.lambdas[i-1], "</DILATION_SIZE>", sep=""),
                             batch.rw[i]), sep="<OPERATION>")}}
  # we finally rewrite the batch
  write.table(batch.new, file=batch, quote=FALSE, col.names=FALSE, row.names=FALSE)}

# Prepare the order sheet
make.grid <- function(seqs=list(lambda1=0:2, lambda2=0:2)) {
  grid <- expand.grid(seqs)
	grid <- cbind(grid, kappa.v=numeric(nrow(grid)), kappa.a=numeric(nrow(grid)), duration=numeric(nrow(grid)))
  return(grid)}

# Launch Dypal in console mode from R
# wd = working directory path. A "dist" folder in the wd must contains dypal.jar and libs
# name = name of the folder to be created
# grid = a grid of lambda values that will also store kappa values
# shape = the shapefile to compare the simulation to
dypal.compute <- function(wd     = "C:/dypal",
                          name   = "",
                          batch  = "",
                          grid   = make.grid(),
                          shape  = "shape/c2000.shp",
                          kappa  = "kappa_temp.txt"){
  setwd(wd)
  # we create an R folder in the working directory if it doesnt exist yet
  dir.create("R", showWarnings=FALSE)
  # we create another folder with the time and name of simulation in the R folder
  logdir <- paste(format(Sys.time(), "%m-%d-%y_%H%M"),name,sep="_")
  dir.create(paste("R", logdir, sep="/"), showWarnings=FALSE)                         
  # we convert paths so that Windows (...) can understand R
  wd.r    <- gsub("/", "\\\\", getwd())
  batch.r <- paste(wd.r, batch, sep="\\")
  shape.r <- paste(wd.r, shape, sep="\\")
  kappa.r <- paste(wd.r, "R", logdir, kappa,  sep="\\")
  dypal.r <- paste(wd.r, "dist\\dypal.jar",   sep="\\")
  java.r  <- "java -jar"
  # we create the kappa-temp file
  write.table(NA, kappa.r, col.n=FALSE, row.n=FALSE)
  # we copy the batch file
  file.copy(batch.r,  paste(wd.r, "R", logdir, batch,  sep="\\"))
  nr  <- nrow(grid)
  results <- grid
  results.path <- paste(wd, "R", logdir,"results-temp.txt", sep="/")
  launch.time <- Sys.time()
  # here we loop to launch Dypal
  for (i in 1:nr) {
    iter.time <- Sys.time()
    batch.write(batch.r, grid[i, 1:2])
    # cosmetics
    cat("*********************************************************************\n")
    cat("R: Performing iteration ", i, " over ", nr, "\n")
    cat("R: Launch time ",format(launch.time),"\n")
    # we launch Dypal
    system(paste(java.r, dypal.r, "-run", batch.r,"-kappa", shape.r,  kappa.r), wait=TRUE)
    end.time <- Sys.time()
    results$duration[i] <- as.numeric(end.time-iter.time)
    # we read the kappa-temp file
    results$kappa.v[i] <- scan(kappa.r) # vectorial kappa
    
    #kappa.ascii
    #ref  <- as.matrix(read.table("C:/dypal/Ascii_2000.txt", h=F, na.strings=-9999))
    #comp <- as.matrix(read.table("C:/dypal/Ascii.txt", h=F, na.strings=-9999))
    #results$kappa.a[i] <- kappa.F(ref, comp)
    
    # cosmetics
    cat("R: Kappa.vect     = ", results$kappa.v[i], "\n")
#     cat("R: Kappa.ascii    = ", results$kappa.a[i], "\n")
#     cat("R: Run duration   = ", results$duration[i], "\n")
#     cat("R: Ellapsed time  = ", format(end.time - launch.time, dig=2), "\n")
#       end <- launch.time+(mean(results$duration[1:i])*(nr/i))
#     cat("R: End time       = ", format(end, dig=2), "\n")
    cat("*********************************************************************\n")
    # we write temporary results
    res.path <- paste(wd,"R",logdir,"results.txt", sep="/")
    write.table(results, res.path, quote=FALSE, row.names=FALSE)}
  # at the end we write the final results file
  write.table(results, res.path, quote=FALSE, row.names=FALSE)
  file.remove(kappa.r)}

# # # # # ##### Simulation example####################################
dypal.compute(name   = "voro_paper_validate3",
              batch  = "batch-voro_paper.xml",
              grid   =  make.grid(list(lambda1=seq(1, 2, 0.1),
                                       lambda2=seq(2, 3, 0.1))),
              shape  = "shapes\\voro_after.shp")




######################################################################
require(MultiScaleR)
rmap <- function(path) {return(as.matrix(read.table(path, skip=6, na.s=-9999)))} 
hsi        <- rmap("C:/dypal/layer-HSI.txt")
flammap.ff <- rmap("C:/dypal/layer-flammap_fullforest.txt")

flammap.fs <- rmap("C:/dypal/layer-flammap_fullsavanna.txt")
flammap.fs <- flammap.fs[round(seq(1, 496, length=497)), round(seq(1, 732, length=733))]
flammap.fs[is.na(hsi)] <- NA                         

curvature  <- rmap("C:/dypal/layer-curvature.txt")
curvature  <- curvature[, round(seq(1, 734, length=733))]
curvature[is.na(hsi)] <- NA                         

                                                  
slope      <- rmap("C:/dypal/layer-slope.txt")
roads      <- rmap("C:/dypal/layer_distlog200_m.txt")

mhm.plot(hsi)
mhm.plot(flammap.ff)
mhm.plot(flammap.fs)
mhm.plot(curvature)
mhm.plot(slope)
mhm.plot(roads)

hsi.m <- mhm(hsi, NEI=crc, wdw.r=seq(5, 50, 5))
mhm.plot(apply(hsi.m, 1:2, mean))
write.table(apply(hsi.m, 1:2, mean), file="C:/dypal/layer-hsiM.txt", quote=F, col.n=F, row.n=F)

flammap.ff.m <- mhm(flammap.ff, NEI=crc, wdw.r=seq(5, 50, 5))
mhm.plot(apply(flammap.ff.m, 1:2, mean))

flammap.fs.m <- mhm(flammap.fs, NEI=crc, wdw.r=seq(5, 50, 5))
mhm.plot(apply(flammap.fs.m, 1:2, mean))

curvature.m <- mhm(curvature, NEI=crc, wdw.r=seq(5, 50, 5))
mhm.plot(apply(curvature.m, 1:2, mean))

slope.m <- mhm(slope, NEI=crc, wdw.r=seq(5, 50, 5))
mhm.plot(apply(slope.m, 1:2, mean))

