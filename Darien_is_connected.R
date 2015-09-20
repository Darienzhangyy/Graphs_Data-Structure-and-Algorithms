gt = list(A = list(edges   = c(2L),
                   weights = c(1 )),
          B = list(edges   = c(3L),
                   weights = c(1 )),
          C = list(edges   = c(5L),
                   weights = c(1 )),
          D = list(edges   = c(2L),
                   weights = c(1 )),
          E = list(edges   = c(4L,6L),
                   weights = c(1,1  )),
          F = list(edges   = c(),
                   weights = c())
)


####################################################
is_connected_child <- function(g,v1,v2){
  #past a vector to record the visited nodes
  past <- c(NULL)
  #childnode is a vector of neighbor nodes
  childnode <- g[[v1]]$edges
  #get a vector called childnode_bool with TRUE when its values corresponds to the values in childnode 
  #and FALSE elsewhere
  childnode_bool <- rep(FALSE,length(names(g)))
  for (i in 1:length(childnode)){
    for (j in 1:length(childnode_bool)){
      if (j==childnode[i]) {
        childnode_bool[j]=TRUE
      }
    }
  }
  #"names(g)[childnode_bool" creates a vector containing the labels of the neighbor nodes    
  if (v2 %in% names(g)[childnode_bool]){
    return(TRUE);stop
  }else{
    #use the fucntion recursively on all the neighbornodes
    for (k in 1:length(childnode)){
      #past keeps track of all the visited nodes,renew past
      past <- c(past,v1)
      #get the label of one of the new neighbor node in this iteration of k
      new=names(g)[g[[v1]]$edges[k]]
      #If the new node is not visited before
      if(!new %in% past){
#Recursion: Implement is_connected_child to the new node to check whether it is conected to the destination
        is_connected_child(g,new,v2)
      }
    }
  }
}

is_connected_child (gt,"A","B")
is_connected_child (gt,"A","C")