#################################################################################
######################### Dypal Console Grid Simulations ########################
#################################################################################

# The functions below aim at launching a set of simulation using Dypal in console mode
# Dypal is called in dypal.compute() that writes a new batch through batch.write()
# upon lambda{i} values previously defined in a grid through make.grid().
# For each row of the grid made, we thus write the batch, launch Dypal 
# and then wait for the end of the simulation.
# The kappa value calculated by comparing simulated landscapes and a reference shape is then read
# and stored in a results_temp file. At the end of the simulation, the final results
# file is written in a dedicated folder and temporary files are removed.

# Manipulate and write the .xml batch
batch.write <- function(batch, new.lambdas=c(1,2)) {
  batch.rw  <- readLines(batch)
  # we split the batch according to the '<OPERATION>' balises
  batch.rw  <- strsplit(batch.rw,"<OPERATION>")[[1]]
  if (length(new.lambdas) != (length(batch.rw)-1)){
    stop("The number of new lamdbas provided is not equal to
          the number of lambdas in the batch file")}
  batch.new <- batch.rw[1]
	for (i in 2:length(batch.rw)) {
    # Not really clear but this is where we manipulate the .xml file
	  batch.new <- paste(batch.new,
	               sub("<DILATION_SIZE>.*</DILATION_SIZE>",
	                   paste("<DILATION_SIZE>", new.lambdas[i-1], "</DILATION_SIZE>", sep=""),
	                   batch.rw[i]), sep="<OPERATION>")
	}
  # we finally rewrite the batch
	write.table(batch.new, file=batch, quote=FALSE, col.names=FALSE, row.names=FALSE)
}


# Prepare the order sheet
make.grid <- function(seqs=list(0:2,0:2)) {
  names(seqs) <- paste("lambda",1:length(seqs), sep="")
  grid <- expand.grid(seqs)
  n    <- nrow(grid)
	grid <- cbind(grid, kappa=numeric(n), kappa.ascii=numeric(n),
                duration=numeric(n), poly=numeric(n), errors=numeric(n))
  return(grid)
}

# wd = working directory path. A "dist" folder in the wd must contains dypal.jar and libs
# name = name of the folder to be created
# Launch Dypal in console mode from R
dypal.compute <- function(wd     = "C:/dypal",
                          name   = "",
                          batch  = "",
                          grid   = make.grid(),
                          shape  = "2000_AspectOS_SingleTol10.shp",
                          kappa  = "kappa_temp.txt",
                          shape.ascii = "2000_AspectOS_SingleTol10Ascii.txt",
                          kappa.ascii = "kappa_temp_ascii.txt",
                          poly   = "poly_temp.txt",
                          errors = "errors_temp.txt",
                          ascii  = FALSE){

# a quick test for arguments classes
if(!is.character(wd)    |
   !is.character(name)  |
   !is.character(batch) |
   !is.character(shape) |
   !is.character(kappa) |
   !is.character(poly)  |
   !is.character(errors)){
   stop("All arguments except 'grid' must be character strings")}                 
if(!is.data.frame(grid)) {
  stop("'grid' must be a data.frame")}
    
setwd(wd)

# we create an R folder in the working directory if it doesnt exist yet
dir.create("R", showWarnings=FALSE)

# we create another folder with the time and name of simulation in the R folder
logdir <- paste(format(Sys.time(), "%m-%d-%y_%H%M"),name,sep="_")
dir.create(paste("R", logdir, sep="/"), showWarnings=FALSE)                         

# we convert paths so that Windows (...) can understand R
# NB : dypal have to be in its 'dist' folder, within the wd folder
# NB : change '-Xmx4G' for '-Xmx2G' if you have 4Go of memory, or let it empty
wd.r      <- gsub("/","\\\\",getwd())
batch.r   <- paste(wd.r, batch,               sep="\\")
shape.r   <- paste(wd.r, shape,               sep="\\")
kappa.r   <- paste(wd.r, "R", logdir, kappa,  sep="\\")
poly.r    <- paste(wd.r, "R", logdir, poly, sep="\\")
errors.r  <- paste(wd.r, "R", logdir, errors, sep="\\")
if (ascii) {
  kappa.comp  <- paste(wd.r, "dist\\kappa.jar",  sep="\\")
  kappa.ascii.r <-  paste(wd.r, "R", logdir, kappa.ascii,  sep="\\")
  shape.ascii.r <-  paste(wd.r, "R", logdir, shape.ascii,  sep="\\")
  dypal.r     <- "C:\\dypal\\dist\\dypal.jar -f"
} else {
  dypal.r   <- paste(wd.r, "dist\\dypal.jar",   sep="\\")}
java.r    <- "java -jar -Xmx4G"

nb.lamb <- ncol(grid)-5
nb.run  <- nrow(grid)
results <- grid
results.path      <- paste(wd,"R",logdir,"results-temp.txt", sep="/")

launch.time <- Sys.time()
for (i in 1:nb.run) {
  iter.time <- Sys.time()
  cat("*********************************************************************\n")
  cat("R: Iteration ", i, " over ", nb.run, "\n")
  cat("R: Launch time ",format(launch.time),"\n")
  cat("*********************************************************************\n")  
  batch.write(batch.r, grid[i, 1:nb.lamb])
  system(paste(java.r, dypal.r, batch.r, shape.r, kappa.r, errors.r, poly.r), wait=TRUE)
  kappa  <- scan(kappa.r)
  if (ascii) {
    system(paste(java.r, kappa.comp,
               paste(sub(".shp","",shape),"_",sub(".xml", "",batch),"Ascii.txt", sep=""),
               shape.ascii,
               kappa.ascii.r,2), wait=TRUE)
    results$kappa.ascii[i] <- scan(kappa.ascii.r)
  }
  end.time <- Sys.time()
  #errors <- scan(errors.r)
  #poly   <- scan(poly.r)
  results$kappa[i]  <- kappa
  results$duration[i] <- as.numeric(end.time - iter.time)
  cat("*********************************************************************\n")
  for (j in 1:length(grid[i, 1:nb.lamb])) {
    cat("R: Lambda", j, "            = ", as.numeric(grid[i,j]),"\n")}
    cat("R: Kappa              = ", kappa, "\n")
    cat("R: Ellapsed time      = ", end.time - launch.time,"\n")
    cat("R: Run duration       = ", as.numeric(end.time-launch.time)," sec \n")
    cat("R: Estimated duration = ", launch.time+mean(results$duration[1:i])*(nb.run-i),"sec \n")
    cat("*********************************************************************\n")
  write.table(results, paste(wd,"R",logdir,"results-temp.txt", sep="/") , quote=FALSE, row.names=FALSE)
  }
  write.table(results, paste(wd,"R",logdir,"results.txt", sep="/"), quote=FALSE, row.names=FALSE)
  file.remove(paste(wd,"R",logdir,"results-temp.txt", sep="/"), errors.r, poly.r, kappa.r)
}

##### Temp
grid1 <- make.grid(list(seq(0,20,1),seq(0,3,0.25)))
dypal.compute(wd     = "C:/dypal",
              name   = "kappa-ascii",
              batch  = "batchinter.xml",
              grid   = grid1,
              shape  = "2000_AspectOS_SingleTol10.shp",
              ascii=F)



#################################################################################
######################### Universal grid plotter ################################
#################################################################################
require(rgl)
col.summer <-colorRampPalette(c("deepskyblue","gold","firebrick1"))
surf.res <- function(res) {
  if (missing(res)) {res <- read.table(choose.files(default="C:/dypal/R"), h=T)}
  if (!is.data.frame(res)) {
    stop("A data.frame must be provided")}
  if(sum(names(res) %in% c("lambda1","lambda2","kappa"))!=3) {
    stop("The data.frame must contain 'lambda1','lambda2' and kappa columns")}
  x <- as.numeric(levels(as.factor(res$lambda1)))
  y <- as.numeric(levels(as.factor(res$lambda2))) 
  z <- res$kappa
  if (length(x) * length(y) != length(z)) {
    stop("The nb of kappa values is not equal to nb(lambda1) x nb(lambda2)")}
  m <- matrix(res$kappa, nr=length(x), nc=length(y),
              dimnames=list(x,y), byrow=FALSE)
  plot3d(res$lambda1, res$lambda2, res$kappa, box=F,
    xlab="Lambda1", ylab="Lambda2", zlab="Kappa",
    size=5, col=col.summer(25)[cut(res$kappa, 25, labels=F)])
  surface3d(x, y, m, col=col.summer(25)[cut(m,25,labels=F)])
  rgl.quads(c(0,max(x),max(x),0),
            c(0,0,max(y),max(y)), rep(res$kappa[1],4), col="#FFFFFF")
  q <- which(m>m[1,1], arr=T)
#   text3d(x[q[,"row"]], y[q[,"col"]],
#          m[q[,"row"],
#          q[,"col"]]+par3d("bbox")[6]*0.01,
#          round(m[q[,"row"], q[,"col"]],3),
#          col="black",
#          cex=0.7)
  return(m)
}
surf.res()  
#################################################################################
#################################################################################

