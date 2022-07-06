#Yirun 'Aaron' Wang
#aaron.wang@wisc.edu

rm(list=ls())
install.packages("FITSio")
require("FITSio")
cB58 = readFrameFromFITS("cB58_Lyman_break.fit") 
n=dim(cB58)[1]
files = list.files('data')
n_files = length(files)
distances = vector(mode = "numeric", length = n_files)
best_shift = vector(mode = "numeric", length = n_files)
new_cB58 = (cB58$FLUX-mean(cB58$FLUX))/sd(cB58$FLUX)
for (i in 1:n_files) {
  path = paste(sep="","data/", files[i])
  data = readFrameFromFITS(path)
  spec = data$flux
  mind = Inf
  besti = 1
  for (point in 1:(length(spec)-n+1)) {
    new_spec=(spec[point:(point+n-1)]-mean(spec[point:(point+n-1)]))/sd(spec[point:(point+n-1)])
    d = sqrt(sum((new_cB58-new_spec)^2))
    if (d < mind) {
      mind = d
      besti = point
    }
  }
  distances[i] = mind
  best_shift[i] = besti
}

distance = sort(distances, decreasing = FALSE)
spectrumID = files[order(distances)]
i = best_shift[order(distances)]
df = data.frame(distance, spectrumID, i)
write.csv(df, file = "hw2.csv")