rm(list=ls())
graph1 = list(A = list(edges=c(2L), weights=c(1)),
              B = list(edges=c(3L), weights=c(1)),
              C = list(edges=c(5L), weights=c(1)),
              D = list(edges=c(2L), weights=c(1)),
              E = list(edges=c(4L,6L), weights=c(1,1)),
              F = list(edges=c(), weights=c())
             )
graph2 = list(A = list(edges=c(2L), weights=c(14)),
              B = list(edges=c(3L,4L), weights=c(23,13)),
              D = list(edges=c(1L), weights=c(5)),
              F = list(edges=c(1L,5L), weights=c(43,33)),
              N = list(edges=c(1L,2L,4L), weights=c(33,22,11))
             )
graph3 = list(list(edges=c(2L), weights=c(1)),
              list(edges=c(3L), weights=c(1)),
              list(edges=c(5L), weights=c(1)),
              list(edges=c(2L), weights=c(1)),
              list(edges=c(4L,6L), weights=c(1,1)),
              list(edges=c(), weights=c())
             )
graph4 = list(list(edges=c(2L), weights=c(14)),
              list(edges=c(3L,4L), weights=c(23,13)),
              list(edges=c(1L), weights=c(5)),
              list(edges=c(1L,5L), weights=c(43,33)),
              list(edges=c(1L,2L,4L), weights=c(33,22,11))
             )
graph5 = list(A = list(edges=c(2L), weights=c(1)),
              E = list(edges=c(4L,6L), weights=c(1,1)),
              B = list(edges=c(3L), weights=c(1)),
              D = list(edges=c(2L), weights=c(1)),
              C = list(edges=c(5L), weights=c(1)),
              F = list(edges=c(), weights=c())
             )         
graph6 = list(list(edges=c(2L), weights=c(1)),
              list(edges=c(4L,6L), weights=c(1,1)),
              list(edges=c(3L), weights=c(1)),
              list(edges=c(2L), weights=c(1)),
              list(edges=c(5L), weights=c(1)),
              list(edges=c(), weights=c())
             )        
graph7 = list(list(edges=c(2L), weights=c(1)),
              list(edges=c(4L,6L), weights=c(1,5)),
              list(edges=c(3L), weights=c(1)),
              list(edges=c(2L), weights=c(1)),
              list(edges=c(5L), weights=c(1)),
              list(edges=c(), weights=c())
             )        
graph8 = list(A = list(edges=c(2L,3L,4L,5L), weights=c(14,5,43,33)),
              B = list(edges=c(1L,3L,4L,5L), weights=c(14,23,13,22)),
              D = list(edges=c(1L,2L), weights=c(5,23)),
              F = list(edges=c(1L,2L,5L), weights=c(43,13,11)),
              N = list(edges=c(1L,2L,4L), weights=c(33,22,11))
             )
graph9 = list(A = list(edges=c(2L,3L,4L,5L), weights=c(14,5,43,33)),
              N = list(edges=c(1L,2L,4L), weights=c(33,22,11)),
              D = list(edges=c(1L,2L), weights=c(5,23)),
              F = list(edges=c(1L,2L,5L), weights=c(43,13,11)),
              B = list(edges=c(1L,3L,4L,5L), weights=c(14,23,13,22))
             )
graph10 = list(list(A = list(edges=c(2L), weights=c(1)),
                    B = list(edges=c(3L), weights=c(1)),
                    C = list(edges=c(5L), weights=c(1)),
                    D = list(edges=c(2L), weights=c(1)),
                    E = list(edges=c(4L,6L), weights=c(1,1)),
                    F = list(edges=c(), weights=c())
                   )) # invalid list structure
graph11 = list(A = list(edges=c(2L), weights=c(1)),
               B = list(edges=c(3L), weights=c(1)),
               A = list(edges=c(5L), weights=c(1)),
               D = list(edges=c(2L), weights=c(1)),
               E = list(edges=c(4L,6L), weights=c(1,1)),
               F = list(edges=c(), weights=c())
              ) # non-unique vertex labels
graph12 = list(A = list(edges=c(2L), weights=c(1)),
               B = list(edges=c(3L), weights=c(1)),
               C = list(edges=c(5L), hello=c(1)),
               D = list(edges=c(2L), weights=c(1)),
               E = list(edges=c(4L,6L), weights=c(1,1)),
               F = list(edges=c(), weights=c())
              ) # invalid vertex attribute structure
graph13 = list(A = list(edges=c(2L), weights=c(1)),
               B = list(edges=c(3L), weights=c(1)),
               C = list(edges=c(5L), weights=c(1)),
               D = list(edges=c(2), weights=c(1)),
               E = list(edges=c(4L,6L), weights=c(1,1)),
               F = list(edges=c(), weights=c())
              ) # invalid edge type
graph14 = list(A = list(edges=c(2L), weights=c(1)),
               B = list(edges=c(3L), weights=c(1)),
               C = list(edges=c(5L), weights=c(1)),
               D = list(edges=c(2L), weights=c('1')),
               E = list(edges=c(4L,6L), weights=c(1,1)),
               F = list(edges=c(), weights=c())
              ) # invalid weight type
graph15 = list(A = list(edges=c(2L), weights=c(1)),
               B = list(edges=c(3L), weights=c(1)),
               C = list(edges=c(5L), weights=c(1)),
               D = list(edges=c(8L), weights=c(1)),
               E = list(edges=c(4L,6L), weights=c(1,1)),
               F = list(edges=c(), weights=c())
              ) # edge to nonexistent vertex
graph16 = list(A = list(edges=c(2L), weights=c(1)),
               B = list(edges=c(3L), weights=c(1)),
               C = list(edges=c(5L), weights=c(-1)),
               D = list(edges=c(2L), weights=c(1)),
               E = list(edges=c(4L,6L), weights=c(1,1)),
               F = list(edges=c(), weights=c())
              ) # nonpositive weight
graph17 = list(A = list(edges=c(2L), weights=c(1)),
               B = list(edges=c(3L), weights=c(1)),
               C = list(edges=c(5L), weights=c(1)),
               D = list(edges=c(2L), weights=c(1)),
               E = list(edges=c(4L,6L), weights=c(1)),
               F = list(edges=c(), weights=c())
              ) # missing weight
graph18 = list(list(edges=c(2L,3L,4L,5L), weights=c(14,5,43,33)),
               list(edges=c(1L,3L,4L,5L), weights=c(14,23,13,22)),
               list(edges=c(1L,2L), weights=c(5,23)),
               list(edges=c(1L,2L,5L), weights=c(43,13,11)),
               list(edges=c(1L,2L,4L), weights=c(33,22,11))
              )



# randomGraph()
######################################################################################################
# Input - length, the number of desired vertices; 
#         maxWeight, the maximum desired edge weight; 
#         seed, optional seed number for reproducability.
#
# Output - graph, a graph object;
#          seed, the seed used to construct graph.
#
# Description - Constructs a valid graph object with (length) vertices and edge weights no 
#               larger than (maxWeight).

randomGraph = function(length, maxWeight=10, seed=NULL) {
  seed = ifelse(is.null(seed), round(runif(1, 1, 1e5)), seed)
  set.seed(seed)
  g = list()
  for(i in seq(length)) {
    g[i] = NA
  }
  filler = function(x) {
    n = ceiling(runif(1, 1, sqrt(length)))
    uniqueEdges = unique(as.integer(ceiling(runif(n, 1, length))))
    n = length(uniqueEdges)
    filled = list(edges=uniqueEdges, 
                  weights=as.numeric(ceiling(runif(n, 1e-10, maxWeight))))
    return(filled)
  }
  g = lapply(g, filler)
  names(g) = as.character(seq(length))
  return(list(graph=g, seed=seed))
}


# adjacencyMatrix()
######################################################################################################
# Input - g, a graph object.
#
# Output - adjacency, an adjacency matrix labeled with vertex names.
#
# Description - Constructs adjacency matrix from graph object, with vertices as dimnames.

adjacencyMatrix = function(g) {
  adjacency = matrix(0, nrow=length(g), ncol=length(g))
  adjDF = as.data.frame(adjacency)
  rownames(adjDF) = colnames(adjDF) = names(g)
  childNodes = lapply(g, function(x) { names(g)[x$edges] } )
  edgeWeights = lapply(g, function(x) { x$weights } )
  for(i in seq(length(g))) {
    adjacency[i, match(unlist(childNodes[i]), colnames(adjDF))] = unlist(edgeWeights[match(names(g)[i], colnames(adjDF))])
  }
  dimnames(adjacency) = list(names(g), names(g))
  return(adjacency)
}


# tidyLabels()
######################################################################################################
# Input - g, a graph object; 
#         v1, a vertex label in g; 
#         v2, a vertex label in g.
#
# Output - G, a graph object; 
#          V1, a vertex label in G; 
#          V2, a vertex label in G; 
#          Warn, the content of any warning; 
#          Exit, whether to return from function.
#
# Description - Labels graph and ensures vertices are valid matches.

tidyLabels = function(g, v1=1, v2=2) {
  warn = NULL
  exit = F
  charNodes = any(c(is.character(v1), is.character(v2)))
  charGraph = !is.null(names(g))
  if(charNodes==T) {
    if(charGraph==F) {
      names(g) = as.character(seq(length(g)))
    }
    if(anyNA(suppressWarnings(as.integer(c(v1, v2))))) {
      if(charGraph==F & length(g)<=26) {
        names(g) = LETTERS[seq(length(g))]
        unnamed = F
        warn = 'Unlabeled graph vertices now labeled A-Z.'
        warning(warn, call.=F)
      } else {
        v1 = toupper(as.character(v1))
        v2 = toupper(as.character(v2))
      }
    } else {
      if(max(as.integer(c(v1, v2)))>length(g)) {
        warn = 'Invalid vertices; please enter valid vertices.'
        exit = T
        return(list(G=g, V1=v1, V2=v2, Warn=warn, Exit=exit))
      }
      v1 = names(g)[as.integer(v1)]
      v2 = names(g)[as.integer(v2)]
      warn = 'Vertex inputs interpreted as indexes.'
      warning(warn, call.=F)
    }
  } else {
    v1 = as.integer(v1)
    v2 = as.integer(v2)
    if(charGraph==F) {
      names(g) = as.character(seq(length(g)))
    }
    v1 = names(g)[v1]
    v2 = names(g)[v2]
  }
  if(!all(c(v1, v2) %in% names(g))) {
    warn = 'Invalid vertices; please enter valid vertices.'
    exit = T
    return(list(G=g, V1=v1, V2=v2, Warn=warn, Exit=exit))
  }
  return(list(G=g, V1=v1, V2=v2, Warn=warn, Exit=exit))
}


# connectionChecker()
######################################################################################################
# Input - g, a graph object.
#
# Output - a three-column data frame containing every possible vertex pair and a logical indicating
#          whether the pair is connected.
#
# Description - Creates sorted matrix of vertex pairs and applies is_connected().

connectionChecker = function(g) {
  tidied = suppressWarnings(tidyLabels(g, v1=names(g)[1], v2=names(g)[2]))
  g = tidied$G
  possibilities = t(cbind(combn(names(g), 2), 
                          combn(names(g), 2)[c(2,1),], 
                          matrix(rep(names(g), 2), nrow=2, byrow=T)))
  possibilities = possibilities[order(possibilities[,1], possibilities[,2]),]
  out = cbind(as.data.frame(possibilities, stringsAsFactors=F), 
              connected=apply(possibilities, 1, function(x) { suppressWarnings(is_connected(g, x[1], x[2])) } ))
  percentage = sum(out$connected)/nrow(out)
  out = list(out, percentagePairsConnected=percentage)
  return(out)
}


# followVertexInner()
######################################################################################################
# Input - path, a path object (a primary list of paths, each containing a secondary list with two 
#                              vectors: $vertices and $weight);
#         A, an adjacency matrix
#
# Output - a path object with an extra (top) layer of hierarchical listing reflecting parent vertices.
#
# Description - Takes one step out from the terminal vertex of every path, returning an updated path
#               list which takes into account edge weights. An extra layer of listing requires removal.

followVertexInner = function(path, A=a) {
  tail = path$vertices[length(path$vertices)]
  out = as.data.frame(cbind(unname(A[tail,!match(A[tail,], 0, nomatch=F)]), 
                            colnames(A)[!match(A[tail,], 0, nomatch=F)]))
  dimnames(out) = list(seq(nrow(out)), c('weight', 'vertices'))
  out = unname(apply(out, 1, function(x) { list(vertices=c(path$vertices, unname(x[2])), 
                                                weight=sum(path$weight, as.integer(unname(x[1])))) } ))
  return(out)
}


# followVertex()
######################################################################################################
# Input - pathList, a path object with an extra (top) layer of hierarchical listing;
#         a, an adjacency matrix.
#
# Output - a path object.
#
# Description - Removes the extra layer of listing.

followVertex = function(pathList, a) { 
  out = unlist(lapply(pathList, followVertexInner, A=a), recursive=F)
  if(length(out)==0) { return(NULL) } 
  return(out)
}

# traceVertexInner()
######################################################################################################
# Input - path, a path object (a primary list of paths, each containing a secondary list with two 
#                              vectors: $vertices and $weight);
#         A, an adjacency matrix
#
# Output - a path object with an extra (top) layer of hierarchical listing reflecting child vertices.
#
# Description - Takes one step back from the initial vertex of every path, returning an updated path
#               list which takes into account edge weights. An extra layer of listing requires removal.

traceVertexInner = function(path, A=a) {
  head = path$vertices[1]
  out = as.data.frame(cbind(unname(A[!match(A[,head], 0, nomatch=F),head]),
                            rownames(A)[!match(A[,head], 0, nomatch=F)]))
  dimnames(out) = list(seq(nrow(out)), c('weight', 'vertices'))
  out = unname(apply(out, 1, function(x) { list(vertices=c(unname(x[2]), path$vertices),
                                                weight=sum(as.integer(unname(x[1])), path$weight)) } ))
  return(out)
}


# traceVertex()
######################################################################################################
# Input - pathList, a path object with an extra (top) layer of hierarchical listing;
#         a, an adjacency matrix.
#
# Output - a path object.
#
# Description - Removes the extra layer of listing.

traceVertex = function(pathList, a) { 
  out = unlist(lapply(pathList, traceVertexInner, A=a), recursive=F)
  if(length(out)==0) { return(NULL) } 
  return(out)
}


# getTails()
######################################################################################################
# Input - pathList, a path object.
#
# Output - a vector of terminal vertices for all paths in pathList.
#
# Description - Subsets on list elements and converts to a vector.

getTails = function(pathList) { unlist(lapply(pathList, function(x) { x$vertices[length(x$vertices)] } )) }


# getHeads()
######################################################################################################
# Input - pathList, a path object.
#
# Output - a vector of initial vertices for all paths in pathList.
#
# Description - Subsets on list elements and converts to a vector.

getHeads = function(pathList) { unlist(lapply(pathList, function(x) { x$vertices[1] } )) }


# weights()
######################################################################################################
# Input - pathList, a path object.
#
# Output - a vector of weights for all paths in pathList.
#
# Description - Subsets on list elements and converts to a vector.

weights = function(pathList) { unlist(lapply(pathList, function(x) { x$weight } )) }


# is_linked()
######################################################################################################
# Input - pathDn, a descending path object;
#         pathUp, an ascending path object;
#         heads, a vector of initial vertices for descending paths;
#         tails, a vector of terminal vertices for ascending paths.
#
# Output - FALSE if no ascending path links to a descending path;
#          list(linked, a logical TRUE;
#               best, a scalar indicating the weight of the weighted shortest path through (heads) and
#                     (tails);
#               bestPath, a character vector containing the vertices of the weighted shortest path.)
#
# Description - Checks whether the intersection of terminal and initial vertex sets is empty. If not, 
#               calculates the path weights for each path through the intersection, reporting the path 
#               and weight of the smallest.

is_linked = function(pathDn, pathUp, heads, tails) {
  bridges = NULL
  if(length(intersect(heads, tails))>0) {
    for(i in seq(length(intersect(heads, tails)))) {
      tmp = data.frame(bridge=intersect(heads, tails)[i],
                       minDn=min(weights(pathDn)[match(intersect(heads, tails)[i], tails)]),
                       minUp=min(weights(pathUp)[match(intersect(heads, tails)[i], heads)]),
                       stringsAsFactors=F)
      bridges = rbind(bridges, tmp)
    }
    best = min(rowSums(bridges[,2:3]))
    bestRow = bridges[(rowSums(bridges[,2:3])==best),]
    topHalf = pathDn[which((weights(pathDn)==bestRow[1,2]) & (tails==bestRow[1,1]))]
    botHalf = pathUp[which((weights(pathUp)==bestRow[1,3]) & (heads==bestRow[1,1]))]
    bestPath = c(unlist(topHalf, recursive=F)$vertices, unlist(botHalf, recursive=F)$vertices[-1])
    return(list(linked=T, best=best, bestPath=bestPath))
  } else {
    return(list(linked=F))
  }
}


# is_valid()
######################################################################################################
# Input - g, a graph object.
# 
# Output - TRUE if g is valid, FALSE if not.
# 
# Description - Validate the graph object to ensure that it meets all requirements: 
#               (1) Check that object is a list of lists. 
#               (2) Check if there are names for the primary list that they are all unique. 
#               (3) Check that each secondary list contains only edges and weights vectors.
#               (4) Check edges and weights vectors are of the appropriate type. 
#               (5) Check that there are no edges to non-existent vertices. 
#               (6) Check that all weights are strictly greater than 0. 
#               (7) Check that every edge has a weight.

is_valid = function(g, requireLabels=T) {
  warnMe = function(index) {
    warningSet = c('Invalid list structure to graph object.',
                   'Invalid vertex label(s) detected.',
                   'Invalid vertex attribute structure(s) detected.',
                   'Invalid edge type(s) detected.',
                   'Invalid weight type(s) detected.',
                   'Nonexistent edge(s) detected.',
                   'Nonpositive weight(s) detected.',
                   'Edge-weight mismatch(es) detected.')
    return(warning(warningSet[index], call.=F))
  }
  if(!is.list(g)) { warnMe(1); return(F) }
  if(!is.list(unlist(g, recursive=F))) { warnMe(1); return(F) }
  if(is.list(unlist(unlist(g, recursive=F), recursive=F))) { warnMe(1); return(F) }
  if(requireLabels) { if(is.null(names(g))) { warnMe(2); return(F) } }
  if(suppressWarnings(any(unique(names(g))!=names(g)))) { warnMe(2); return(F) }
  if(!all(unlist(lapply(g, function(x) { ifelse(names(x)!=c('edges', 'weights'), F, T) } )))) { warnMe(3); return(F) }
  secondaryTypes = unlist(lapply(unlist(g, recursive=F), function(x) { typeof(x) } ))
  edges = as.logical(seq(length(secondaryTypes)) %% 2)
  if(length(setdiff(secondaryTypes[edges], c('integer', 'NULL')))>0) { warnMe(4); return(F) }
  if(length(setdiff(secondaryTypes[!edges], c('double', 'NULL')))>0) { warnMe(5); return(F) }
  h = g;  names(h) = NULL;  k = unlist(h, recursive=F)
  if(!all(unique(unlist(k[names(k)=='edges'])) %in% seq(length(g)))) { warnMe(6); return(F) }
  if(!all(unique(unlist(k[names(k)=='weights']))>0)) { warnMe(7); return(F) }
  if(length(unlist(k[names(k)=='edges']))!=length(unlist(k[names(k)=='weights']))) { warnMe(8); return(F) }
  return(T)
}


# is_undirected()
######################################################################################################
# Input - g, a graph object;
#         requireLabels, a logical (default is TRUE) indicating whether graph must be labeled.
# 
# Output - TRUE if undirected, FALSE if not.
# 
# Description - Check if the graph object is undirected, this is true if all directed edges have a 
#               complementary directed edge with the same weight in the opposite direction.
#
# Idea - Create an adjacency matrix and check if it's symmetric.

is_undirected = function(g, requireLabels=T) {
  if(!is_valid(g, requireLabels)) { return(invisible(NULL)) }
  if(is.null(names(g))) { names(g) = as.character(seq(length(g))) }
  a = adjacencyMatrix(g)
  return(ifelse(identical(a, t(a)), T, F))
}


# is_isomorphic
######################################################################################################
# Input - g1, a graph object; 
#         g2, a graph object;
#         requireLabels, a logical (default is TRUE) indicating whether graph must be labeled.
#
# Output - TRUE if g1 and g2 are isomorphic, FALSE if not.
#
# Description - Check if the graph objects are isomorphic, meaning all vertices, edges, and weights 
#               are identical. Comparison of vertices should be based on names not indexes, indexes 
#               should only be used if vertex labels are not defined.
#
# Note - Likely slow for graphs with more than 150 vertices due to for loop implementation, with a
#        maximum of (n^2-n)/2 iterations when vertices of g1 and g2 are perfectly reversed.

is_isomorphic = function(g1, g2, requireLabels=T) {
  if(!(is_valid(g1, requireLabels) & is_valid(g2, requireLabels))) { return(invisible(NULL)) }
  if(length(g1)!=length(g2)) { return(F) }
  if(suppressWarnings(any(sort(names(g1))!=sort(names(g2))))) { return(F) }
  if((is.null(names(g1)) | is.null(names(g2)))) { names(g1) = names(g2) = NULL }
  if(all(is.null(names(g1)) & is.null(names(g2)))) {
    g2nodes = seq(length(g2))
    for(i in seq(length(g1))) {
      matches = vector(mode='logical', length=length(g1))
      for(j in g2nodes) {
        matches[j] = identical(g1[i], g2[j])
      }
      if(!any(matches)) { return(F) } 
      g2nodes = setdiff(g2nodes, which(matches==T)[1])
    }
  } else {
    a1 = adjacencyMatrix(g1); a2 = adjacencyMatrix(g2)
    if(!identical(a1, a2)) { return(F) }
  }
  return(T)
}


# is_connected()
######################################################################################################
# Input - g, a graph object; 
#         v1, a vertex label in g; 
#         v2, a vertex label in g;
#         requireLabels, a logical (default is TRUE) indicating whether graph must be labeled;
#         unweighted, a logical (default is FALSE) indicating whether edge weights should be used.
#
# Output - TRUE if there exists a path from v1 to v2 in g, FALSE if not.
#
# Description - Determine if there is any path between vertex v1 and vertex v2 in graph g. 
#               If v1 or v2 are not in g then throw an error.
#
# Idea - Construct an adjacency matrix, g, and take it to the 1st, 2nd, ..., nth power, where n is the
#        number of vertices. If the [v1,v2] entry of the ith product is ever positive, then the 
#        unweighted shortest path is of length i. If the [v1,v2] entry of all products is zero, the
#        vertices are unconnected.

is_connected = function(g, v1, v2, requireLabels=T, unweighted=F) {
  if(!is_valid(g, requireLabels)) { return(invisible(NULL)) }
  if(!any(is.character(c(v1, v2)) | mode(c(v1, v2))=='numeric')) { 
    warning('Invalid vertices; please enter valid vertices.', call.=F)
    return(F)
  }
  if(!is.logical(unweighted)) { warning('Unweighted must be logical; defaulting to FALSE.', call.=F) }
  tidied = suppressWarnings(tidyLabels(g, v1, v2))
  g = tidied$G; v1 = tidied$V1; v2 = tidied$V2
  if(!is.null(tidied$Warn)) { warning(tidied$Warn, call.=F) }
  if(tidied$Exit==T) { return(F) }
  g = h = adjacencyMatrix(g)
  for(i in seq(sqrt(length(g)))) { 
    if(h[v1, v2]>0) { 
      if(unweighted==T) {
        return(list(connected=T, length_of_shortest_unweighted_path=i))
      } else {
        return(T) 
      }
    }
    h = h %*% g
  } 
  return(F)
}


# is.connected()                                                an alternate version of is_connected()
######################################################################################################
# Input - g, a graph object; 
#         v1, a vertex label in g; 
#         v2, a vertex label in g.
#
# Output - TRUE if there is a path from v1 to v2 in g, FALSE if not.
#
# Description - Determine if there is any path between vertex v1 and vertex v2 in graph g. If v1 or v2
#               are not in g then throw an error.
#
# Idea - keep tracing out paths until 
#        (1) a path from v1 to v2 is found, or 
#        (2) the set of terminal nodes contains nothing new (i.e. all active paths are in cycle).

is.connected = function(g, v1, v2) {
  if(!any(is.character(c(v1, v2)) | mode(c(v1, v2))=='numeric')) { 
    warning('Invalid vertices; please enter valid vertices.', call.=F)
    return(F)
  }
  tidied = suppressWarnings(tidyLabels(g, v1, v2))
  g = tidied$G; v1 = tidied$V1; v2 = tidied$V2
  if(!is.null(tidied$Warn)) { warning(tidied$Warn, call.=F) }
  if(tidied$Exit==T) { return(F) }
  childNodes = lapply(g, function(x) { names(g)[x$edges] } )
  pathsForward = list(v1)
  terminalNodes = NULL
  uniqueCheck = 1
  while((!(v2 %in% terminalNodes) | is.null(terminalNodes)) & (!(uniqueCheck==0) | is.null(terminalNodes))) {
    toAdd = toDelete = children = NULL
    for(i in seq(length(pathsForward))) {
      children = c(children, unname(unlist(childNodes[pathsForward[[i]][length(pathsForward[[i]])]])))
    }  
    uniqueCheck = terminalNodes = length(setdiff(children, unlist(pathsForward)))
    if(!((!(v2 %in% terminalNodes) | is.null(terminalNodes)) & (!(uniqueCheck==0) | is.null(terminalNodes)))) {
      out = F
      break
    }
    for(i in seq(length(pathsForward))) {
      children = unname(unlist(childNodes[pathsForward[[i]][length(pathsForward[[i]])]]))
      if(length(children)==0) {
        toDelete = c(toDelete, i)
      } else if(length(children)==1) {
        pathsForward[[i]] = c(pathsForward[[i]], children)
      } else {
        for(j in seq(2, length(children), 1)) {
          toAdd = append(toAdd, list(c(pathsForward[[i]], children[j])))
        }
        pathsForward[[i]] = c(pathsForward[[i]], children[1])
      }
    }
    pathsForward = ifelse(is.null(toDelete), pathsForward, pathsForward[-toDelete])
    if(!is.null(toAdd)) { 
      pathsForward = append(pathsForward, toAdd) 
    }
    terminalNodes = unlist(lapply(pathsForward, function(x) { x[length(x)] } ))
  }
  if(v2 %in% terminalNodes) {
    return(T)
  } else if(uniqueCheck==0) {
    return(F)
  } else {
    return(warning('Unknown error encountered.', call.=F))
  }
}


# shortest_path()
######################################################################################################
# Input - g, graph object; 
#         v1, a vertex label in g; 
#         v2, a vertex label in g;
#         requireLabels, a logical (default is TRUE) indicating whether graph must be labeled;
#         weight, a logical (default is FALSE) indicating whether the edge weight of the shortest path 
#         should be reported.
#
# Output - a vector of the names (or indexes if unlabeled) of vertices that make up the shortest 
#          path, in order. If (weight) is TRUE, a list is returned, with the path vector ($shortestPath)
#          as the first element and the edge weight ($pathWeight) as the second element. If no path
#          exists [is_connected(g) = FALSE], returns an empty vector.
#
# Description - Find the shortest path from vertex v1 to vertex v2 using the edges of graph g. Note
#               that there may not be a unique solution for any given graph, you are only required to 
#               return one path.
#
# Idea - (1) Initialize a shortest path length as infinity. 
#        While upward and downward paths are active (non-NULL):
#        | (2) Take one step down from terminal vertices of the downward path.
#        | (3) Eliminate initial vertices of the upward path if the edge weight to reach them 
#        |     exceeds (shortest path - min(edge weight to reach terminal vertices).
#        | (4) If the graph is directed, eliminate previously traversed terminal vertices.
#        | (5) Check for complete paths, choose shortest, compare to shortest path, and update.
#        | (6) Take one step up from initial vertices of the upward path.
#        | (7) Eliminate terminal vertices of the downward path if the edge weight to reach  
#        |     them exceeds (shortest path - min(edge weight to reach initial vertices).
#        | (8) If the graph is directed, eliminate previously traversed initial vertices.
#        | (9) Check for complete paths, choose shortest, compare to shortest path, and update.

shortest_path = function(g, v1, v2, requireLabels=T, weight=F) {
  if(!is_valid(g, requireLabels)) { return(invisible(NULL)) }
  if(!any(is.character(c(v1, v2)) | mode(c(v1, v2))=='numeric')) { 
    return(warning('Invalid vertices; please enter valid vertices.', call.=F))
  }
  if(is.logical(weight)==F) { warning('Invalid weight; defaulting to FALSE.', call.=F) }
  tidied = suppressWarnings(tidyLabels(g, v1, v2))
  g = tidied$G; v1 = tidied$V1; v2 = tidied$V2
  if(!is.null(tidied$Warn)) { warning(tidied$Warn, call.=F) }
  if(tidied$Exit==T) { return(invisible(NULL)) }
  if(!suppressWarnings(is_connected(g, v1, v2))) { 
    warning('Unconnected vertices; please enter connected vertices.', call.=F)
    return(character(0)) 
  } 
  a = adjacencyMatrix(g)
  best = Inf
  pathDn = list(list(vertices=v1, weight=0))
  pathUp = list(list(vertices=v2, weight=0))
  tails = v1
  heads = v2
  oldHeads = oldTails = NULL
  while(!(is.null(pathDn) | is.null(pathUp))) {
    pathDn = followVertex(pathDn, a)
    if(is.null(pathDn)) { break }
    keepers = (best - min(weights(pathDn))) > weights(pathUp)
    pathUp = pathUp[keepers]
    oldTails = unique(c(oldTails, tails))
    heads = getHeads(pathUp)
    tails = getTails(pathDn)
    if(!is_undirected(g) & !all(is.na(match(tails, oldTails)))) {
        pathDn = pathDn[-which(tails %in% oldTails)]
        tails = tails[-which(tails %in% oldTails)]
    }
    linkCheck = is_linked(pathDn, pathUp, heads, tails)
    if(linkCheck$linked==T) { 
      possBest = linkCheck$best
      possBestPath = linkCheck$bestPath 
      if(possBest<best) { best=possBest; bestPath=possBestPath }
    }
    pathUp = traceVertex(pathUp, a)
    if(is.null(pathUp)) { break }
    keepers = (best - min(weights(pathUp))) > weights(pathDn)
    pathDn = pathDn[keepers]
    oldHeads = unique(c(oldHeads, heads))
    heads = getHeads(pathUp)
    tails = getTails(pathDn)
    if(!is_undirected(g) & !all(is.na(match(heads, oldHeads)))) {
      pathUp = pathUp[-which(heads %in% oldHeads)]
      heads = heads[-which(heads %in% oldHeads)]
    }
    linkCheck = is_linked(pathDn, pathUp, heads, tails)
    if(linkCheck$linked==T) { 
      possBest = linkCheck$best
      possBestPath = linkCheck$bestPath 
      if(possBest<best) { best=possBest; bestPath=possBestPath }
    }
  }
  if(weight==T) {
    return(list(shortestPath=bestPath, pathWeight=best))
  } else {
    return(bestPath)
  }
}