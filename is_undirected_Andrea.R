##Test Data
graph1 = list(A = list(edges   = c(2L),
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

graph2 = list(A = list(edges   = c(2L),
                       weights = c(14)),
              B = list(edges   = c(3L,4L),
                       weights = c(23,13)),
              D = list(edges   = c(1L),
                       weights = c(5) ),
              F = list(edges   = c(1L,5L),
                       weights = c(43,33)),
              N = list(edges   = c(1L,2L,4L),
                       weights = c(33,22,11))
)

graph3 =list(A = list(edges = c(2L, 3L),
                      weights = c(1, 1)),
             B = list(edges = c(1L, 4L),
                      weights = c(1, 1)),
             C = list(edges = c(1L),
                      weights = c(1)),
             D = list(edges = c(2L),
                      weights =c(1))
)
                 

#graph4 is slightly different from graph3 in weight of D to B
graph4 =list(A = list(edges = c(2L, 3L),
                      weights = c(1, 1)),
             B = list(edges = c(1L, 4L),
                      weights = c(1, 1)),
             C = list(edges = c(1L),
                      weights = c(1)),
             D = list(edges = c(2L),
                      weights =c(3))
)




##Function- is_undirected: Check if the graph object is undirected,
##                         this is true if all directed edges have a complementary directed 
##                         edge with the same weight in the opposite direction.
##
##Input: g, graph object
##Output: TRUE if g is valid, FALSE ifnot
is_undirected <- function(g){
##initialize flag, which helps us to distinguish if g is undirected   
  flag_edg <- 1
  flag_wgt <- 1
  flag_is_undrt <-TRUE
  
#loop for every vertex in the graph     
  for (i in 1 : length(g)){    
#loop for all vertices connected to the vertex i    
    for (j in 1 : length(g[[i]][[1]])) {    
      vertex_num <- g[[i]][[1]][j]
#check if the directed edge g[[i]][[1]][j] have a complementary directed edge
      flag_edg <- match(i, g[[vertex_num]][[1]], 
                    nomatch = 0)
#if there is no complementary edge, break directly and the graph is not undirected
      if ( flag_edg == 0 ){
        break
      }
#if there is a complementary edge, check if the edge g[[i]][[1]][j] has the same 
#weight with the complementary edge.
      flag_wgt <- match(g[[i]][[2]][j], g[[j]][[2]][flag_edg],
                        nomatch = 0)
#if they do not have the same weight, break directly and the graph is not undirected
      if ( flag_wgt == 0){
        break
      }
    }
    if (flag_edg == 0 || flag_wgt == 0){
      break
    }
  }
#after finishing the loops, if flag_edg and flag_wgt are still nonzero, it means the
#graph is undirected
  if (flag_edg !=0 && flag_wgt !=0){
    flag_is_undrt = TRUE
  }
  else flag_is_undrt = FALSE
  return (flag_is_undrt)
}

is_undirected(graph1)
is_undirected(graph2)
is_undirected(graph3)
is_undirected(graph4)
