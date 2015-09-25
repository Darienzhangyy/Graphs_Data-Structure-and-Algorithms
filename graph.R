---
title: "graph.R"
output: pdf_document
---
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
graph3 = list(list(edges=c(2L), weights=c(1)),
              list(edges=c(3L), weights=c(1)),
              list(edges=c(5L), weights=c(1)),
              list(edges=c(2L), weights=c(1)),
              list(edges=c(4L,6L), weights=c(1,1)),
              list(edges=c(), weights=c())
) # unlabeled graph1
graph4 = list(list(edges=c(2L), weights=c(14)),
              list(edges=c(3L,4L), weights=c(23,13)),
              list(edges=c(1L), weights=c(5)),
              list(edges=c(1L,5L), weights=c(43,33)),
              list(edges=c(1L,2L,4L), weights=c(33,22,11))
) # unlabeled graph2
graph5 = list(A = list(edges=c(2L), weights=c(1)),
              E = list(edges=c(4L,6L), weights=c(1,1)),
              B = list(edges=c(3L), weights=c(1)),
              D = list(edges=c(2L), weights=c(1)),
              C = list(edges=c(5L), weights=c(1)),
              F = list(edges=c(), weights=c())
) # reordered graph1         
graph6 = list(list(edges=c(2L), weights=c(1)),
              list(edges=c(4L,6L), weights=c(1,1)),
              list(edges=c(3L), weights=c(1)),
              list(edges=c(2L), weights=c(1)),
              list(edges=c(5L), weights=c(1)),
              list(edges=c(), weights=c())
) # reordered graph3       
graph7 = list(list(edges=c(2L), weights=c(1)),
              list(edges=c(4L,6L), weights=c(1,5)),
              list(edges=c(3L), weights=c(1)),
              list(edges=c(2L), weights=c(1)),
              list(edges=c(5L), weights=c(1)),
              list(edges=c(), weights=c())
) # reweighted graph6        
graph8 = list(A = list(edges=c(2L,3L,4L,5L), weights=c(14,5,43,33)),
              B = list(edges=c(1L,3L,4L,5L), weights=c(14,23,13,22)),
              D = list(edges=c(1L,2L), weights=c(5,23)),
              F = list(edges=c(1L,2L,5L), weights=c(43,13,11)),
              N = list(edges=c(1L,2L,4L), weights=c(33,22,11))
) # undirected
graph9 = list(A = list(edges=c(2L,3L,4L,5L), weights=c(14,5,43,33)),
              N = list(edges=c(1L,2L,4L), weights=c(33,22,11)),
              D = list(edges=c(1L,2L), weights=c(5,23)),
              F = list(edges=c(1L,2L,5L), weights=c(43,13,11)),
              B = list(edges=c(1L,3L,4L,5L), weights=c(14,23,13,22))
) # reordered graph8
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
graphX = list(A = list(edges=c(2L), weights=c(1)),
              B = list(edges=c(3L), weights=c(1)),
              C = list(edges=c(5L), weights=c(1)),
              D = list(edges=c(2L,7L), weights=c(1,1)),
              E = list(edges=c(4L,6L), weights=c(1,1)),
              F = list(edges=c(), weights=c()),
              X = list(edges=c(), weights=c())
) # two null vertices, both connected


# randomGraph()
######################################################################################################
# Input - length, the number of desired vertices; 
#         maxWeight, the maximum desired edge weight; 
#         seed, optional seed number for reproducability;
#         undirected, logical (default: FALSE) indicating whether graph should be undirected;
#         adjacency, logical (default: FALSE returns list-based graph) indicating whether adjacency 
#         matrix should be returned.
#
# Output - graph, a graph object;
#          seed, the seed used to construct graph.
#
# Description - Constructs a valid graph object with (length) vertices and edge weights no 
#               larger than (maxWeight).

randomGraph = function(length, maxWeight=10, seed=NULL, undirected=F, adjacency=F) {
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
  if(undirected==T) {
    g = adjacencyMatrix(g)
    g[lower.tri(g)] = 0
    g[lower.tri(g)] = g[lower.tri(g)] + t(g)[lower.tri(g)]
    if(!adjacency) { g = listGraph(g) }
  } else {
    if(adjacency) { g = adjacencyMatrix(g) }
  }
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
  if(length(makeNull)>0) { g[makeNull] = list(list(edges=NULL, weights=NULL)) }
  if(is.null(dimnames(a)[[1]])) { names(g) = NULL }
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
  if(length(g1)!=length(g2)) { return(F) }
  if(suppressWarnings(any(sort(names(g1))!=sort(names(g2))))) { return(F) }
  sortAndLabel = function(vertex, g) {
    df = data.frame(edges=unlist(names(g)[vertex$edges]), weights=unlist(vertex$weights))
    df = df[order(df$edges), , drop=F]
    return(list(edges=df$edges, weights=df$weights))
  }
  g1 = lapply(g1, sortAndLabel, g=g1); g2 = lapply(g2, sortAndLabel, g=g2)
  g1 = g1[sort(names(g1))];            g2 = g2[sort(names(g2))]
  a1 = adjacencyMatrix(g1);            a2 = adjacencyMatrix(g2)
  if(!identical(a1, a2)) { return(F) }
  return(T)
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
  if(!suppressWarnings(is_valid(g))) { return(stop('Invalid graph.', call.=F)) }
  names(g) = toupper(names(g)); v1 = toupper(v1); v2 = toupper(v2)
  if(!(is.character(v1) & is.character(v2) & all(c(v1, v2) %in% names(g)))) { 
    return(stop('Invalid vertices; please enter valid vertices.', call.=F))
  }
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
      if(pathUp[[1]]$weight==0) {
        return(pathDn[[1]]$vertices)
      } else {
        pathDn = pathDn[-which(tails %in% oldTails)]
        tails = tails[-which(tails %in% oldTails)]
      }
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
  return(bestPath)
}
