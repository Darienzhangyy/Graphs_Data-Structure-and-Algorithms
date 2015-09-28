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
  if(length(makeNull)>0) { 
    g[makeNull] = list(list(edges=NULL, weights=NULL)) 
  }
  if(is.null(dimnames(a)[[1]])) { 
    names(g) = NULL 
  }
  return(g)
}


# is_valid()
######################################################################################################
# Input - g, a graph object.
# 
# Output - TRUE if g is valid, FALSE if not.
# 
# Description - Validate the graph object to ensure that it meets all requirements: 
#               (1) Check that object is a list of non-null lists. 
#               (2) Check that there are names for the primary list that they are all unique. 
#               (3) Check that each secondary list contains only edges and weights vectors.
#               (4) Check that there are no NAs.
#               (5) Check that edges vectors are of the appropriate type. 
#               (6) Check that weights vectors are of the appropriate type.
#               (7) Check that there are no edges to non-existent vertices. 
#               (8) Check that all weights are strictly greater than 0. 
#               (9) Check that every edge has a weight.
#              (10) Check that no edges are duplicated.

is_valid = function(g) {
  warnMe = function(index) {
    warningSet = c('Invalid list structure to graph object.',
                   'Invalid vertex label(s) detected.',
                   'Invalid vertex attribute structure(s) detected; only "edges" and "weights" accepted.',
                   'NA(s) detected.',
                   'Non-integer edge(s) detected.',
                   'Non-numeric weight(s) detected.',
                   'Nonexistent edge(s) detected.',
                   'Nonpositive weight(s) detected.',
                   'Edge-weight mismatch(es) detected.',
                   'Duplicate edge(s) detected.')
    return(warning(warningSet[index], call.=F))
  }
# (1) Check that object is a list of non-null lists. 
  if((!is.list(g)) | (is.null(unlist(g))) | (!is.list(unlist(g, recursive=F))) | (is.list(unlist(unlist(g, recursive=F), recursive=F)))) { 
    warnMe(1); return(F) 
  }
# (2) Check that there are names for the primary list that they are all unique. 
  if((is.null(names(g))) | (length(unique(names(g)))!=length(names(g)))) { 
    warnMe(2); return(F) 
  }
# (3) Check that each secondary list contains only edges and weights vectors.
  if(!all(unlist(lapply(g, function(x) { ifelse(sort(names(x))!=c('edges', 'weights'), F, T) } )))) { 
    warnMe(3); return(F) 
  }
# (4) Check that there are no NAs.
  h = g;  names(h) = NULL;  k = unlist(h, recursive=F)
  if(any(is.na(unlist(k[names(k)=='edges']))) | any(is.na(unlist(k[names(k)=='weights'])))) { 
    warnMe(4); return(F) 
  }
# (5) Check that edges vectors are of the appropriate type. 
  if(typeof(unlist(k[names(k)=='edges']))!='integer') { 
    warnMe(5); return(F) 
  } 
# (6) Check that weights vectors are of the appropriate type.
  if(typeof(unlist(k[names(k)=='weights']))!='double') { 
    warnMe(6); return(F) 
  }
# (7) Check that there are no edges to non-existent vertices. 
  if(!all(unique(unlist(k[names(k)=='edges'])) %in% seq(length(g)))) { 
    warnMe(7); return(F) 
  }
# (8) Check that all weights are strictly greater than 0. 
  if(!all(unique(unlist(k[names(k)=='weights']))>0)) { 
    warnMe(8); return(F) 
  }
# (9) Check that every edge has a weight.
  if(length(unlist(k[names(k)=='edges']))!=length(unlist(k[names(k)=='weights']))) { 
    warnMe(9); return(F) 
  }
# (10) Check that no edges are duplicated.
  m = unlist(g, recursive=F)
  if(!all(unlist(lapply(m[names(k)=='edges'], function(x) { length(unique(x))==length(x) } )))) { 
    warnMe(10); return(F) 
  }
  return(T)
}


# is_undirected()
######################################################################################################
# Input - g, a graph object.
# 
# Output - TRUE if undirected, FALSE if not.
# 
# Description - Check if the graph object is undirected, this is true if all directed edges have a 
#               complementary directed edge with the same weight in the opposite direction.
#
# Idea - Create an adjacency matrix and check if it's symmetric.

is_undirected = function(g) {
  if(!suppressWarnings(is_valid(g))) { 
    return(stop('Invalid graph.', call.=F)) 
  }
  a = adjacencyMatrix(g)
  return(ifelse(identical(a, t(a)), T, F))
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

