source('graph.R')

# minPos()
######################################################################################################
# Input - col, a column in a matrix.
#
# Output - the minimum value greater than zero.
#
# Description - Returns the minimum positive value in a column of integers.

minPos = function(col) { 
  if(sum(col)==0) {
    return(0) 
  } else {
    return(min(col[(col>0)]))
  }
}


# min_span_tree()
######################################################################################################
# Input - g, a graph object.
#
# Output - g, a minimum spanning tree (MST) of the original graph object.
#
# Description - Implements Boruvka's algorithm, as described in the relevant Wikipedia article.

min_span_tree = function(g) {
  if(!suppressWarnings(is_valid(g))) { 
    return(stop('Invalid graph.', call.=F)) 
  }
  if(length(g)==1) { 
    g[[1]] = list(edges=NULL, weights=NULL); return(g) 
  }
# Save the original vertex names and rename as character-coerced integers.
  origNames = names(g)
  names(g) = seq(length(g))
# Construct adjacency matrix, a, and equally-dimensioned matrix of zeroes, b.
  a = b = adjacencyMatrix(g); b[] = 0
# Initialize t, a forest, as a listed set of one-vertex trees.
  t = sapply(seq(nrow(a)), function(x) { c(x) } , simplify=F)
  setLength = Inf
# 
  while(setLength>1) {
    c = list(t[[1]])
#   Run through each set in the forest.
    for(i in seq(length(t))) {
#     Identify the minimally-weighted out edge for each vertex in the set.
      currentRows = a[unlist(t[[i]]),-unlist(t[[i]]), drop=F]
      currentRow = apply(currentRows, 2, minPos)
      minLogical = currentRow==min(currentRow[(currentRow>0)])
      minValue = currentRow[minLogical][1]
      minVertex = as.integer(names(minValue))
      altVertex = as.integer(rownames(currentRows)[currentRows[,as.character(minVertex)]==minPos(currentRows[,as.character(minVertex)])])
      b[altVertex,minVertex] = b[minVertex,altVertex] = minValue
      d = lapply(c, function(x) { c(unlist(t[[i]]), minVertex) %in% unlist(x) } )
      d = unlist(lapply(d, function(x) { any(unlist(x))==T } ))
#     Update the temporary forest, c, based upon the current set.
      if(sum(d)>1) { 
        c[[seq(length(c))[d][1]]] = unique(c(unlist(c[d]), c(unlist(t[[i]]), minVertex)))
        c = c[-(seq(length(c))[d][-1])]
      } else if(sum(d)==1) {
        c[[seq(length(c))[d]]] = unique(c(c[[seq(length(c))[d]]], c(unlist(t[[i]]), minVertex)))
      } else {
        c[[length(c)+1]] = c(unlist(t[[i]]), minVertex)
      }
    }
#   Update the permanent forest, t,  based upon the temporary forest, c, and the size of the forest.
    t = c
    setLength = length(t)
  }
# Rename the vertices.
  g = listGraph(b)
  names(g) = origNames
  return(g)
}
