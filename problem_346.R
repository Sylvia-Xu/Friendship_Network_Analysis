library("igraph")

g = read.graph(file = "facebook_combined.txt", format = "edgelist", directed = F)
#1
# c = is.connected(g)
# d = diameter(g)
# 
# degreesVector <- degree(g)
# h <- hist(degreesVector, breaks=seq(-0.5, by=1 , length.out=max(degreesVector)+2))
# pl <- data.frame(x=h$mids[2:1046], y=h$density[2:1046])
# plot(pl , type="o")
# 
# dat = data.frame(x=h$mids, y=h$density)
# 
# model = lm(y ~ poly(x, n=8, raw=TRUE), data=dat)
# model2 = nls(y ~ exp(a*x+b), data=dat, start = list(a=0, b=0))
# 
# newdat = data.frame(x=seq(0,1046,0.5))
# newdat$y = predict(model, newdat)
# newdat$yy = predict(model2, newdat)
# plot(y ~ x, dat)
# lines(y ~ x, newdat, col = 2, lwd = 2)
# lines(yy ~ x, newdat, col = 3, lwd = 2)

#2
nei = neighbors(g,1)
g1 = setdiff(1:vcount(g),nei)
g2 = setdiff(g1,1)
nei1 <- delete.vertices(g, g2)
nei2 <- delete.vertices(g, g1)

#3
d = degree(g)
core_index = which(d > 200)
core_deg = d[core_index]
degree_average = mean(core_deg)

struct1 = fastgreedy.community(nei1)
struct2 = edge.betweenness.community(nei1)
struct3 = infomap.community(nei1)

ver_size = numeric(vcount(nei1)) + 3
ver_size[1] = 8
ver_col = struct1$membership
plot(nei1, edge.width = 0.1, vertex.size = ver_size, vertex.label = NA, vertex.frame.color = NA, vertex.color = ver_col)

ver_size = numeric(vcount(nei1)) + 3
ver_size[1] = 8
ver_col = struct2$membership
plot(nei1, edge.width = 0.1, vertex.size = ver_size, vertex.label = NA, vertex.frame.color = NA, vertex.color = ver_col)

ver_size = numeric(vcount(nei1)) + 3
ver_size[1] = 8
ver_col = struct3$membership
plot(nei1, edge.width = 0.1, vertex.size = ver_size, vertex.label = NA, vertex.frame.color = NA, vertex.color = ver_col)

#4
struct4 = fastgreedy.community(nei2)
struct5 = edge.betweenness.community(nei2)
struct6 = infomap.community(nei2)

ver_size = numeric(vcount(nei2)) + 3
ver_col = struct4$membership
plot(nei2, edge.width = 0.1, vertex.size = ver_size, vertex.label = NA, vertex.frame.color = NA, vertex.color = ver_col)

ver_size = numeric(vcount(nei2)) + 3
ver_col = struct5$membership
plot(nei2, edge.width = 0.1, vertex.size = ver_size, vertex.label = NA, vertex.frame.color = NA, vertex.color = ver_col)

ver_size = numeric(vcount(nei2)) + 3
ver_col = struct6$membership
plot(nei2, edge.width = 0.1, vertex.size = ver_size, vertex.label = NA, vertex.frame.color = NA, vertex.color = ver_col)

cl <- clusters(nei2)
gccIndex = which.max(cl$csize)
nonGccNodes <- (1:vcount(nei2))[cl$membership != gccIndex]
nei3 <- delete.vertices(nei2, nonGccNodes)

struct7 = fastgreedy.community(nei3)
struct8 = edge.betweenness.community(nei3)
struct9 = infomap.community(nei3)

#5
embed = numeric(0)
disp = numeric(0)
for (i in 1:length(core_index))
{
  nei_index = neighbors(g,core_index[i])
  g1 = setdiff(1:vcount(g),nei_index)
  core_nei <- delete.vertices(g, g1)
  embed = c(embed, degree(core_nei))
  temp = numeric(vcount(core_nei))
  for (j in 1:vcount(core_nei))
  {
    core_nei_push = core_nei
    nei_sub = neighbors(core_nei, j)
    ed = incident(core_nei, j)
    delete.edges(core_nei, ed)
    if(length(nei_sub) <= 1)
    {
      temp[j] = 0
      j = j + 1
    }
    else
    {
      dis = 0
      for (p in 1:(length(nei_sub)-1))
      {
        for(q in (p+1):length(nei_sub))
        {
          dis = dis + distances(core_nei, nei_sub[p], nei_sub[q])
        }
      }
      temp[j] = dis
      core_nei = core_nei_push
    }
  }
  disp = c(disp, temp)
}
h = hist(embed, breaks=seq(-0.5, by=1, length.out=max(embed)+2))
pl <- data.frame(x=h$mids, y=h$density)
plot(pl, type="o")

h = hist(disp, breaks=seq(-0.5, by=1, length.out=max(disp)+2))
pl <- data.frame(x=h$mids, y=h$density)
plot(pl, type="o")

#6
feature = matrix(0,0,4)
for (i in 1:length(core_index))
{
  nei = neighbors(g,core_index[i])
  g1 = setdiff(1:vcount(g),nei)
  g2 = setdiff(g1,1)
  sub_g <- delete.vertices(g, g2)
  
  struct = fastgreedy.community(sub_g)
  for (j in 1:length(struct))
  {
    if(sizes(struct)[j] >= 10)
    {
      nonNodes <- (1:vcount(sub_g))[struct$membership != j]
      com_g <- delete.vertices(sub_g, nonNodes)
      attr_line = numeric(4)
      attr_line[1] = vcount(com_g)
      attr_line[2] = 2 * ecount(com_g) / (attr_line[1] * attr_line[1] - attr_line[1])
      attr_line[3] = transitivity(com_g)
      attr_line[4] = i
      feature = rbind(feature, attr_line)
    }
  }
}


