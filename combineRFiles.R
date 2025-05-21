combineRFiles <- function(folderName="R/", outputFileName="all_funcs_combined_new.R", modelRunFileName="runModelRCS.R") {
  sink(outputFileName)


  filePrefix <- function(fileName="") {
    return(paste(folderName, fileName, sep=""))
  }

  i <- 1
  files <- filePrefix(list.files(filePrefix()))
  files <- c(files, modelRunFileName)

  for (file in files) {
    current_file = readLines(file)
    cat("\n\n#### Current file:",file,"\n\n")
    cat(current_file, sep ="\n")
  }

  sink()

}

outputFileName <- commandArgs(trailingOnly = TRUE)[1]
modelRunFileName <- commandArgs(trailingOnly = TRUE)[2]

if(outputFileName != "") {
  if(modelRunFileName != "") {
    combineRFiles(outputFileName=outputFileName, modelRunFileName=modelRunFileName)
  } else {
    combineRFiles(outputFileName=outputFileName)
  }
} else if (modelRunFileName != "") {
  combineRFiles(modelRunFileName=modelRunFileName)
} else {
  combineRFiles()
}
