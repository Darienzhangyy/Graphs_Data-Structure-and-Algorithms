
is_connected <- function(g,v1,v2){
  if (is_valid(g)==FALSE){
    stop("error");break
  }
  if(v1 %in% names(g) && v2 %in% names(g)){
    path=c(g[[v1]]$edges)
    node=NULL
    for (i in 1:length(g)){
      node=c(node,names(g[path[i]]))
      path=c(path,g[[node[i]]]$edges)
    }
  }
  if (v2 %in% node){
    return(TRUE)
  }else{
    return(FALSE)
  }
}




