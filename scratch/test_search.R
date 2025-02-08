

setwd( R'(C:\Users\James.Thorson\Desktop\Git\VAST\R)' )
files = list.files()
search_terms = c("ThorsonUtilities", "FishStatsUtils", "TMBhelper")

for(j in seq_along(search_terms)){
for(i in seq_along(files)){
  read = scan( files[i], what="character", quiet = TRUE )
  if(length(grep( search_terms[j], read ))){
    stop( "Check ", files[i], " for ", search_terms[j] )
  }
}}

setwd( R'(C:\Users\James.Thorson\Desktop\Git\VAST)' )
devtools::document()
