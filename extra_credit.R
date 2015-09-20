# adjacencyMatrix()
######################################################################################################
# Input - g, a graph object.
#
# Output - a, an adjacency matrix labeled with vertex names.
#
# Description - Constructs adjacency matrix from graph object, with vertices as dimnames.

adjacencyMatrix = function(g) {
  a = matrix(0, nrow=length(g), ncol=length(g))
  adjDF = as.data.frame(a)
  rownames(adjDF) = colnames(adjDF) = names(g)
  childNodes = lapply(g, function(x) { names(g)[x$edges] } )
  edgeWeights = lapply(g, function(x) { x$weights } )
  for(i in seq(length(g))) {
    a[i, match(unlist(childNodes[i]), colnames(adjDF))] = unlist(edgeWeights[match(names(g)[i], colnames(adjDF))])
  }
  dimnames(a) = list(names(g), names(g))
  return(a)
}

# listGraph()
######################################################################################################
# Input - a, an adjacency matrix labeled with vertex names.
#
# Output - g, a graph object.
#
# Description - Constructs list-based graph object from adjacency matrix, 
#               with dimnames as vertex labels.

listGraph = function(a) {
  rowToList = function(row) { list(edges=seq(ncol(a))[(row>0)], weights=unname(row[(row>0)])) }
  g = apply(a, 1, rowToList)
  makeNull = seq(ncol(a))[(rowSums(a)==0 & colSums(a)>0)]
  if(length(makeNull)>0) { g[makeNull] = list(list(edges=NULL, weights=NULL)) }
  if(is.null(dimnames(a)[[1]])) { names(g) = NULL }
  return(g)
}

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
  if(!suppressWarnings(is_valid(g))) { return(stop('Invalid graph.', call.=F)) }
  if(length(g)==1) { g[[1]] = list(edges=NULL, weights=NULL); return(g) }
  origNames = names(g)
  names(g) = seq(length(g))
  a = b = adjacencyMatrix(g); b[] = 0
  t = sapply(seq(nrow(a)), function(x) {c(x)} , simplify=F)
  setLength = Inf
  while(setLength>1) {
    c = list(t[[1]])
    for(i in seq(length(t))) {
      currentRows = a[unlist(t[[i]]),-unlist(t[[i]]), drop=F]
      currentRow = apply(currentRows, 2, minPos)
      minLogical = currentRow==min(currentRow[(currentRow>0)])
      minValue = currentRow[minLogical][1]
      minVertex = as.integer(names(minValue))
      altVertex = as.integer(rownames(currentRows)[currentRows[,as.character(minVertex)]==minPos(currentRows[,as.character(minVertex)])])
      b[altVertex,minVertex] = b[minVertex,altVertex] = minValue
      d = lapply(c, function(x) { c(unlist(t[[i]]), minVertex) %in% unlist(x) } )
      d = unlist(lapply(d, function(x) { any(unlist(x))==T } ))
      if(sum(d)>1) { 
        c[[seq(length(c))[d][1]]] = unique(c(unlist(c[d]), c(unlist(t[[i]]), minVertex)))
        c = c[-(seq(length(c))[d][-1])]
      } else if(sum(d)==1) {
        c[[seq(length(c))[d]]] = unique(c(c[[seq(length(c))[d]]], c(unlist(t[[i]]), minVertex)))
      } else {
        c[[length(c)+1]] = c(unlist(t[[i]]), minVertex)
      }
    }
    t = c
    setLength = length(t)
  }
  g = listGraph(b)
  names(g) = origNames
  return(g)
}