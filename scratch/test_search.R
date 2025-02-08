

setwd( R'(C:\Users\James.Thorson\Desktop\Git\VAST\R)' )
files = list.files()

for(i in seq_along(files)){
  read = scan( files[i], what="character" )
  if(length(grep( "jpeg", read ))){
    stop( "Check ", files[i] )
  }
}

setwd( R'(C:\Users\James.Thorson\Desktop\Git\VAST)' )
devtools::document()
