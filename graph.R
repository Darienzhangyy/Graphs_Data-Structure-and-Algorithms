# graphs
######################################################################################################
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
graph8 = list(A = list(edges=c(2L,3L,4L,5L), weights=c(14,5,43,33)),
              B = list(edges=c(1L,3L,4L,5L), weights=c(14,23,13,22)),
              D = list(edges=c(1L,2L), weights=c(5,23)),
              F = list(edges=c(1L,2L,5L), weights=c(43,13,11)),
              N = list(edges=c(1L,2L,4L), weights=c(33,22,11))
) # undirected


# randomGraph()
######################################################################################################
# Input - length, the number of desired vertices; 
#         maxWeight, the maximum desired edge weight; 
#         seed, optional seed number for reproducability;
#         undirected, logical (default: FALSE) indicating whether graph should be undirected;
#         adjacency, logical (default: FALSE returns list-based graph) indicating whether adjacency 
#                    matrix should be returned.
#
# Output - graph, a graph object;
#          seed, the seed used to construct graph.
#
# Description - Constructs a valid graph object with (length) vertices and edge weights no 
#               larger than (maxWeight).

randomGraph = function(length, maxWeight=10, seed=NULL, undirected=F, adjacency=F) {
# Check if random seed is provided and sets seed accordingly.
  seed = ifelse(is.null(seed), round(runif(1, 1, 1e5)), seed)
  set.seed(seed)
# Create a list-based graph object with NA vertices.
  g = list()
  for(i in seq(length)) {
    g[i] = NA
  }
# Define a function to generate random edges and weights, given a graph of specified size.
  filler = function(x) {
    n = ceiling(runif(1, 1, sqrt(length)))
    uniqueEdges = unique(as.integer(ceiling(runif(n, 1, length))))
    n = length(uniqueEdges)
    filled = list(edges=uniqueEdges, 
                  weights=as.numeric(ceiling(runif(n, 1e-10, maxWeight))))
    return(filled)
  }
# Lapply function to fill in list-based graph object, providing character-coerced integers as vertex names.
  g = lapply(g, filler)
  names(g) = as.character(seq(length))
# If instructed, create symmetric adjacency matrix from the upper triangle, and return graph as specified by the user.
  if(undirected) {
    g = adjacencyMatrix(g)
    g[lower.tri(g)] = 0
    g[lower.tri(g)] = g[lower.tri(g)] + t(g)[lower.tri(g)]
    if(!adjacency) { 
      g = listGraph(g) 
    }
  } else {
    if(adjacency) { 
      g = adjacencyMatrix(g) 
    }
  }
# Return the random graph AND the seed used to generate it.
  return(list(graph=g, seed=seed))
}


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
  if(length(out)==0) { 
    return(NULL) 
  } 
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
  if(length(out)==0) { 
    return(NULL) 
  } 
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
# Check for connected v1-v2 path.
  if(length(intersect(heads, tails))>0) {
#   Construct a data frame which contains, for each path:
#     bridge: the node connecting the downward and upward paths, 
#     minDn:  the minimally weighted downward path to the bridge node, and
#     minUp:  the minimally weighted upward path to the bridge node.
    for(i in seq(length(intersect(heads, tails)))) {
      tmp = data.frame(bridge=intersect(heads, tails)[i],
                       minDn=min(weights(pathDn)[match(intersect(heads, tails)[i], tails)]),
                       minUp=min(weights(pathUp)[match(intersect(heads, tails)[i], heads)]),
                       stringsAsFactors=F)
      bridges = rbind(bridges, tmp)
    }
#   Identify the shortest of these v1-v2 paths, and extract the weight and the nodes from which it's composed.
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

is_undirected<- function(g){ 
        #check if the graph is valid; if it's not return false
 if (invisible(is_valid(g))==F) {
   stop ("error")} 
 if (length(g)<=1) {
   return(TRUE)} else {
#create a matrix of 0's with nrow=ncol=length(g)         
   m0<-matrix(0, nrow=length(g), ncol =length(g))   
#give initial value n=0 
   n=0
 for (j in 1:length(g)) {
 for (i in 1:length(g)) {
#store each entry of m0[i,j] with a weight that correponds to its edge for each vertex
#vertex i is transformed to m[i,]
#if vertex i directly connects to vertex j(j is a edge value), store its corresponding weight value into m0[i,j]
m0[i,g[[i]]$edges]<-g[[i]]$weights    
#if i and j have the same weight and i and j are not equal and they're not zero's
#record n=n+1
if (m0[i,j]==m0[j,i]& i!=j & m0[i,j]!=0  ) {n=n+1
}  else { if (length(g)==2){ n=n+1
} else {n=n
}}}}
#if n>0 and m0 is symetric that means vertex i and vertex j are undirected 
if (n>0 & (isSymmetric(m0))){
return (TRUE)} else{
return (FALSE)}}}


# is_isomorphic
######################################################################################################
# Input - g1, a graph object; 
#         g2, a graph object.
#
# Output - TRUE if g1 and g2 are isomorphic, FALSE if not.
#
# Description - Check if the graph objects are isomorphic, meaning all vertices, edges, and weights 
#               are identical. Comparison of vertices should be based on names not indexes, indexes 
#               should only be used if vertex labels are not defined.

is_isomorphic = function(g1, g2) {
  if(!(suppressWarnings(is_valid(g1)) & suppressWarnings(is_valid(g2)))) { 
    return(stop('Invalid graph.', call.=F)) 
  }
# Check that the number and names of vertices are the same for each graph.
  if((length(g1)!=length(g2)) | (suppressWarnings(any(sort(names(g1))!=sort(names(g2)))))) { 
    return(F) 
  }
# Define a function to sort the vertices of a list-based graph object.
  sortAndLabel = function(vertex, g) {
    df = data.frame(edges=unlist(names(g)[vertex$edges]), weights=unlist(vertex$weights))
    df = df[order(df$edges), , drop=F]
    return(list(edges=df$edges, weights=df$weights))
  }
  g1 = lapply(g1, sortAndLabel, g=g1); g2 = lapply(g2, sortAndLabel, g=g2)
  g1 = g1[sort(names(g1))];            g2 = g2[sort(names(g2))]
# Convert list-based graph objects to adjacency matrices and check that they are identical.
  a1 = adjacencyMatrix(g1);            a2 = adjacencyMatrix(g2)
  if(!identical(a1, a2)) { 
    return(F) 
  }
  return(T)
}


# is_connected()
######################################################################################################
# Input - g, a graph object; 
#         v1, a vertex label in g; 
#         v2, a vertex label in g.
#
# Output - TRUE if there exists a path from v1 to v2 in g, FALSE if not.
#
# Description - Determine if there is any path between vertex v1 and vertex v2 in graph g. 
#               If v1 or v2 are not in g then throw an error.
#
# Idea - Construct an unweighted adjacency matrix, g, and take it to the 1st, 2nd, ..., nth power, 
#        where n is the number of vertices. If the [v1,v2] entry of the ith product is ever positive,  
#        then the unweighted shortest path is of length i. If the [v1,v2] entry of all products is zero, 
#        the vertices are unconnected. SOURCE: Strang (2009)

is_connected = function(g, v1, v2) {
  if(!suppressWarnings(is_valid(g))) { 
    return(stop('Invalid graph.', call.=F)) 
  }
  names(g) = toupper(names(g)); v1 = toupper(v1); v2 = toupper(v2)
# Check the validity of vertices: character vectors of length one whose elements are named vertices in g.
  if(!(is.character(v1) & is.character(v2) & all(c(v1, v2) %in% names(g)) & length(v1)==1 & length(v2)==1)) { 
    return(stop('Invalid vertices; please enter valid vertices.', call.=F))
  }
# Convert g to an unweighted adjacency matrix, a, and duplicate as h.
  a = adjacencyMatrix(g)
  a[(a>0)] = 1; h = a
# Iteratively raise a to higher powers, checking the [v1,v2] entry for positive values.
  for(i in seq(length(g))) { 
    if(h[v1, v2]>0) { 
      return(T) 
    }
    h = h %*% a
  } 
  return(F)
}

# shortest_path()
######################################################################################################
# Input - g, graph object; 
#         v1, a vertex label in g; 
#         v2, a vertex label in g.
#
# Output - a vector of the names (or indexes if unlabeled) of vertices that make up the shortest 
#          path, in order. If no path exists, returns an empty vector.
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

shortest_path = function(g, v1, v2) {
  if(!suppressWarnings(is_valid(g))) { 
    return(stop('Invalid graph.', call.=F)) 
  }
  names(g) = toupper(names(g)); v1 = toupper(v1); v2 = toupper(v2)
# Check the validity of vertices: character vectors of length one whose elements are named vertices in g.
  if(!(is.character(v1) & is.character(v2) & all(c(v1, v2) %in% names(g)) & length(v1)==1 & length(v2)==1)) { 
    return(stop('Invalid vertices; please enter valid vertices.', call.=F))
  }
# Check that v1 and v2 are in fact connected.
  if(!suppressWarnings(is_connected(g, v1, v2))) { 
    warning('Unconnected vertices; please enter connected vertices.', call.=F)
    return(character(0)) 
  } 
# Convert the list-based graph object to an adjacency matrix.
  a = adjacencyMatrix(g)
# Initialize the shortest path, downward/upward paths, terminal/initial nodes, and visited terminal/initial nodes.
  best = Inf
  pathDn = list(list(vertices=v1, weight=0))
  pathUp = list(list(vertices=v2, weight=0))
  tails = v1
  heads = v2
  oldHeads = oldTails = NULL
# Follow the downward path and trace the upward path until one has been pruned to NULL.
  while(!(is.null(pathDn) | is.null(pathUp))) {
#   Update the downward path and break if NULL.
    pathDn = followVertex(pathDn, a)
    if(is.null(pathDn)) { 
      break 
    }
#   Prune any upward path that is too long, given the updated downward paths, to form a new shortest path.
    keepers = (best - min(weights(pathDn))) > weights(pathUp)
    pathUp = pathUp[keepers]
#   Update the visited terminal node set and update the terminal and initial nodes.
    oldTails = unique(c(oldTails, tails))
    heads = getHeads(pathUp)
    tails = getTails(pathDn)
#   After at least one iteration on a directed graph, prune the downward path and terminal nodes, unless v1=v2.
    if(!is_undirected(g) & !all(is.na(match(tails, oldTails)))) {
      if(pathUp[[1]]$weight==0) {
        return(pathDn[[1]]$vertices)
      } else {
        pathDn = pathDn[-which(tails %in% oldTails)]
        tails = tails[-which(tails %in% oldTails)]
      }
    }
#   Check for v1-v2 paths via is_linked() and update current shortest path accordingly.
    linkCheck = is_linked(pathDn, pathUp, heads, tails)
    if(linkCheck$linked==T) { 
      possBest = linkCheck$best
      possBestPath = linkCheck$bestPath 
      if(possBest<best) { 
        best=possBest; bestPath=possBestPath 
      }
    }
#   Update the upward path and break if NULL.
    pathUp = traceVertex(pathUp, a)
    if(is.null(pathUp)) { 
      break 
    }
#   Prune any downward path that is too long, given the updated upward paths, to form a new shortest path.
    keepers = (best - min(weights(pathUp))) > weights(pathDn)
    pathDn = pathDn[keepers]
#   Update the visited initial node set and update the terminal and initial nodes.
    oldHeads = unique(c(oldHeads, heads))
    heads = getHeads(pathUp)
    tails = getTails(pathDn)
#   After at least one iteration on a directed graph, prune the upward path and initial nodes.
    if(!is_undirected(g) & !all(is.na(match(heads, oldHeads)))) {
      pathUp = pathUp[-which(heads %in% oldHeads)]
      heads = heads[-which(heads %in% oldHeads)]
    }
#   Check for v1-v2 paths via is_linked() and update current shortest path accordingly.
    linkCheck = is_linked(pathDn, pathUp, heads, tails)
    if(linkCheck$linked==T) { 
      possBest = linkCheck$best
      possBestPath = linkCheck$bestPath 
      if(possBest<best) { 
        best=possBest; bestPath=possBestPath 
      }
    }
  }
  return(bestPath)
}

