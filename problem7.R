library("igraph")

#setwd("Downloads/gplus")

cir_files <- list.files(pattern = "*.circles",recursive = FALSE)

multi_cir <- list()
count = 0
#lapply(cir_files, function(x) {
for (i in 1: length(cir_files)){
  t <- tryCatch(read.table(cir_files[i], fill = TRUE, numerals = "no.loss"), error=function(e) NULL)
  if (!is.null(t)){
    if (dim(t)[1] == 5 ){
      edge_file <- gsub("circles", "edges", cir_files[i])
      count = count + 1
      g = read.graph(file = edge_file,format = "ncol", directed=TRUE)
      summary(g)
      #add the ego node to the graph
      g <- add.vertices(g, 1)
      for (m in 1: (length(V(g))-1))
      {
        g<-add.edges(g, c(length(V(g)), m))
      }
      summary(g)
      m_community = infomap.community(g)
      for (j in 1:dim(t)[1])
      {
        member = numeric(dim(t)[2])
        for (k in 2: dim(t)[2]){
          if(!is.null(t[j,k])){
            index = which(V(g)$name==as.character(t[j,k]))
            if (length(index) == 1) member[k-1] = m_community$membership[index]
            else member[k-1] = 0
          }
        }
        title=paste("Communicty Membership in Circle", as.character(j), "of Person", as.character(i))
        h1 <- hist(member, main = title, breaks=seq(-0.5, by=1, length.out = max(member)+2))
      }
      
    }
  }
}