args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
e <- new.env()
load(args[1], envir=e)
x <- e$CellChatDB.human$interaction
stopifnot(nrow(x)==3233, all(c('ligand.symbol','receptor.symbol','annotation') %in% names(x)))
write.csv(x, args[2])
cat('Exported', nrow(x), 'interaction rows\n')
