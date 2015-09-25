---
title: "graph.R"
output: pdf_document
---

        ###1. Function - is_valid
        is_valid <- function(g){
                count <- 0
                ##Check if there are names for the primary list that they are all unique
                primary <- names(g)
                if (sum(duplicated(primary)) != 0){count=count+1}
                if (length(primary)==0){count=count+1}
                for (i in 1 : length(g)){
                        ##Check that object is a list of lists
                        if (typeof(g[[i]])!="list") {count=count+1;break}
                        ##Check that each secondary list contains only edges and weights vectors 
                        if (length(g[[i]])!=2) {count=count+1;break}
                        if (!names(g[[i]])[1] %in% c("edges","weights")) {count=count+1;break}
                        if (!names(g[[i]])[2] %in% c("edges","weights")) {count=count+1;break}
                        if (names(g[[i]])[1]==names(g[[i]])[2]) {count=count+1;break}
                        ##Check duplicated edges
                        if (sum(duplicated(g[[i]][[1]])) != 0){count=count+1;break}
                        ##Check that edges and weights vectors that are of the appropriate type
                        if (!is.null(g[[i]][[1]])){
                                if (!is.vector(g[[i]][[1]])) {count=count+1;break}
                                if (min(g[[i]][[1]]) <= 0 | NA %in% g[[i]][[1]]) {count=count+1;break}
                                if (!is.numeric(g[[i]][[2]])) {count=count+1;break}
                                if (!is.integer(g[[i]]$edges)) {count=count+1;break}
                                ##Check that there are not any edges to non-existent vertices
                                if (max(g[[i]][[1]])>length(g))  {count=count+1;break}
                                ##Check that all weights are not less than or equal to 0
                                if (min(g[[i]][[2]])<=0 | NA %in% g[[i]][[2]]) {count=count+1;break}
                                ##Check that every edge has a weight
                                if (length(g[[i]][[1]]) != length(g[[i]][[2]])){count=count+1;break}
                        }
                }
                ##As long as "count" is not equal to zero, there must be a mistake somewhere
                if (count == 0) {print (TRUE)}
                else {print(FALSE)}
        }


is_undirected<- function(g){ 
        #check if the graph is valid; if it's not return false
        if (is_valid(g)==F) {
                stop ("error"); break
        } else { 
                if (length(g)<=1) {
                        return(TRUE); break} else {
                                
                                #create a matrix of 0's with nrow=ncol=length(g)         
                                m0<-matrix(0, nrow=length(g), ncol =length(g))   
                                #give initial value n=0 
                                n=0
                                for (j in 1:length(g)) {
                                        for (i in 1:length(g)) {
                                                #store each entry of m0[i,j] with a weight that correponds to its edge for each vertex
                                                #vertex i is transformed to m[i,]
                                                #if vertex i directly connects to vertex j(j is a edge value), store its corresponding weight value into m0[i,j]
                                                m0[i,g[[i]]$edges]<-1/g[[i]]$weights    
                                                #if i and j have the same weight and i and j are not equal and they're not zero's
                                                #record n=n+1
                                                if (m0[i,j]==m0[j,i]& i!=j & m0[i,j]!=0  ) {
                                                        n=n+1
                                                }  else { if (length(g)==2){
                                                        n=n+1
                                                } else {
                                                        n=n
                                                } 
                                                }
                                        }
                                }
                                
                                #if n>0 and m0 is symetric that means vertex i and vertex j are undirected 
                                if (n>0 & (isSymmetric(m0))){
                                        return (TRUE); break 
                                }else{
                                        return (FALSE)
                                }
                        }
        }
}



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
