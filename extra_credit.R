# adjacencyMatrix()
######################################################################################################
# Input - g, a graph object.
#
# Output - a, an adjacency matrix labeled with vertex names.
#
# Description - Constructs adjacency matrix from graph object, with vertices as dimnames.

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

