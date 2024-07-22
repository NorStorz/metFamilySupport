exportData <- function(obj){

 
  ## First three lines have to be added
  write.table(obj, file = "data/test.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
}