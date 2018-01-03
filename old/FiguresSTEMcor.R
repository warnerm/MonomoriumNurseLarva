load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")

d = dist(fpkm)
mypar()
hc <- hclust(d)
