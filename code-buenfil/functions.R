#Please reference Meila and Bao JMLR paper for definitions of various statistics within the model

#maps permutation to a matrix. The matrix has dimensions maxitem x tmax, 
#which can be specified
#or deduced from the single permutation passed in
pi_to_Pi <- function(pi1,maxitem=NULL,tmax = NULL) {
  if (is.null(maxitem)) {
    maxitem = max(unlist(pi1))
  }
  if (is.null(tmax)) {
    tmax = length(pi1)
  }
  Pi = matrix(data= 0, nrow=maxitem,ncol = tmax)
  for (j  in 1:length(pi1)) {
    Pi[pi1[[j]],j] = 1
  }
  a = 1:length(pi1)
  return(Pi)
}

#maps the list of all permutations to a list of their matrix representations. 
#The matrix representations have dimensions
#maxitem x tmax, where maxitem is the largest item index that appears across 
#all permutations in S_N, and tmax is
#the length of the largest permutation
#we can also specify maxitem to be larger than the max item that appears in S_N
#we can also specify tmax to be larger than the largest permutation length that 
#appears in S_N
SN_to_Pis <- function(S_N,maxitem=NULL,tmax = NULL) {
  if (is.null(maxitem)) {
    maxitem = max(unlist(S_N))
  }
  if (is.null(tmax)) {
    lengths = lapply(S_N,length)
    tmax = max(unlist(lengths))
  }
  
  Pis = list()
  for (i in 1:length(S_N)) {
    Pis[[i]] = pi_to_Pi(S_N[[i]],maxitem,tmax)
  }
  return(Pis)
}

#maps the matrix representation of a permutation to its list representation
Pi_to_pi <- function(Pi) {
  indices = which(Pi==1,arr.ind=TRUE)
  return(as.list(indices[,1]))
}

#computes the N_j statistics, contained in a vector N of size tmax. N_j, the jth
#entry of N, is the number of permutations in
#S_N that have length greater than or equal to j.
#tmax can also be specified.
computeN <- function(S_N,tmax = NULL) {
  if (is.null(tmax)) {
    lengths = unlist(lapply(S_N,length))
    tmax = max(unlist(lengths))
  }
  M = length(S_N)
  
  mat1 = matrix(data = rep(lengths,tmax),nrow = M)
  a = 1:tmax
  mat2 = matrix(data = rep(1:tmax,M), byrow=TRUE,ncol=tmax)
  comparematrix = mat1 >= mat2
  
  N = apply(comparematrix,2,sum)
  return(N)
}

#compute sum of entries of the lower triangular portion of a matrix A(below the 
#diagonal, not including it)
L <- function(A) {
  return(sum(A[lower.tri(A, diag=FALSE)]))
}

#computes the q statistic from the sample in matrix form Pis. q is stored as a 
#matrix of size maxitem x tmax, where the (i,j) entry equals
#the number of times item i is observed in rank j of the data.
#to compute Pis automatically from S_N, you can specify S_N, which will override
#what is specified for Pis
compute_q <- function(Pis,S_N = NULL) {
  if (!is.null(S_N)) {
    Pis = SN_to_Pis(S_N)
  }
  
  #sum matrices in Pis
  return(Reduce('+', Pis))
}

#compute the Q statistic from the sample S_N. Q is stored as an array of size 
#tmax x maxitem x maxitem, where Q[j,,] = Q_j
#we can also specifiy maxitem or tmax to be larger.
#also returns Pis for convenience
compute_Q <- function(S_N,maxitem=NULL,tmax = NULL) {
  if (is.null(maxitem)) {
    maxitem = max(unlist(S_N))
  }
  if (is.null(tmax)) {
    lengths = unlist(lapply(S_N,length))
    tmax = max(unlist(lengths))
  }
  M = length(S_N)
  Pis = SN_to_Pis(S_N,maxitem,tmax)
  
  #first we compute Pitilde from Pi, where in each row, we replace all 0s before
  #a 1 with ones. If there are no 1s, then
  #the entire row is replaced with 1s. We do this for each Pi to obtain Ptildes,
  #a list of the Ptilde matrices.
  #we also convert Pis into an array here
  Ptildes = array(data = 0, dim = c(M,maxitem,tmax))
  Parray = array(data = 0, dim = c(M,maxitem,tmax))
  for (m in 1:M) {
    Pi = Pis[[m]]
    Ptildes[m,,] = Pi
    Parray[m,,] = Pi
    for (k in 1:maxitem) {
      if (length(which(Pi[k,]==1)) == 0) {
        Ptildes[m,k,] = 1
      }
      else {
        index = which(Pi[k,]==1)
        Ptildes[m,k,1:index] = 1
      }
    }
  }
  #computation of Q
  Q = array(data =0, dim = c(tmax,maxitem,maxitem))
  for (i in 1:maxitem) {
    for (k in 1:maxitem) {
      temp = list()
      for (m in 1:M) {
        temp[[m]] = (1-Ptildes[m,k,])&Parray[m,i,]
      }
      Q[,i,k] = Reduce('+', temp)
    }
  }
  
  return(list(Q=Q,Pis=Pis))
}

#computes R using q and Q. R is a list of length tmax, where each entry is a 
#matrix of size maxitem x maxitem
compute_R <- function(q,Q) {
  maxitem = dim(q)[1]
  tmax = dim(q)[2]
  R = list()
  for (j in 1:tmax) {
    matj = matrix(data = rep(q[,j],maxitem),ncol = maxitem)
    R[[j]] = matj - Q[j,,]
  }
  
  return(R)
}

#computes the representations of a permutation pi in terms of s_j. sig is 
#permutation in list form of
#sigma_0. sig should have at least all the elements present in pi.
pi_to_sj <- function(pi1,sig) {
  
  t = length(pi1)
  maxitem = max(unlist(sig))
  Pi = pi_to_Pi(pi1,maxitem,t)
  
  
  #compute the matrix representation of sig
  Sig = pi_to_Pi(sig)
  
  s = rep(0,length(pi1))
  
  #compute R from Q and q
  #(following Proposition 2)
  q = compute_q(list(Pi))
  Q = compute_Q(list(pi1),maxitem,t)$Q
  R = compute_R(q,Q)
  for (j in 1:t) {
    s[j] = L(t(Sig)%*%R[[j]]%*%Sig)
  }
  
  return(s)
}

#constructs pi using s and sigma. pi is a list, while s is a vector.
#sigma needs to be at least length N, where N = length(s) + max(s).
#if sigma is not specified, then sigma = 1,2, ... n is automatically used
sj_to_pi <- function(s,sigma0 = NULL) {
  t = length(s)
  maxs = max(s)
  n = t + maxs
  if (is.null(sigma0)) {
    sigma0 = as.list(1:n)
  }
  if (length(sigma0) < n) {
    print("error: sigma0 wasn't long enough")
    return(0)
  }
  
  pi = list()
  for (j in 1:t) {
    pi[[j]] = sigma0[[s[j]+1]]
    sigma0 = sigma0[-(s[j]+1)]
  }
  
  return(pi)
}

#EstimateSigmaTheta function from Meila 2010 paper. Finds maximum likelihood 
#sigma and theta, given the R_j (list of maxitem x maxitem matrices), 
#N (vector of length tmax containing the N_j), and theta_init, the initial 
#theta_j parameter values for j =1, ... tmax.
#if theta_init is not specified, then theta_init is set to be all 1s.
#if max_iterations not specified, then max_iterations set to be 10
#TODO: implement convergence criterion?
EstimateSigmaTheta <- function(R,N,theta_init = NULL,max_iterations = NULL) {
  tmax = length(N)
  if  (is.null(theta_init)) {
    theta_init = rep(1,tmax)
  }
  if  (is.null(max_iterations)) {
    max_iterations = 10
  }
  theta = theta_init
  print("Iteration:")
  for (i in 1:max_iterations) {
    print(i)
    Rtheta = lapply(1:tmax,function(j) theta[j]*R[[j]])
    Rthetasum = Reduce('+', Rtheta)
    
    #compute sigma1
    sigma1 = greedyR(Rthetasum)
    Sigma1 = pi_to_Pi(sigma1)
    
    #compute theta
    #first compute L(t(Sigma1)%*%R[[j]]%*%Sigma1) for every, we need to check 
    #that some aren't 0
    temp = rep(0,tmax)
    for (j in 1:tmax) {
      temp[j] = L(t(Sigma1)%*%R[[j]]%*%Sigma1)
      if (temp[j] == 0) {
        #if L(t(Sigma1)%*%R[[j]]%*%Sigma1) is 0, then we set this to be 
        #something small
        temp[j] = .001
      }
    }
    #using R element-wise operations here
    theta = log(1 + N/temp)
  }
  
  return(list(sigma1 = sigma1,theta = theta))
}

#greedyR algorithm from Meila 2010. Estimates sigma in a greedy way, 
#input is sufficient statistic Rthetasum 
#(typically calculated within function EstimateSigmaTheta) 
#(this is a maxitem x maxitem matrix)
greedyR <- function(Rthetasum) {
  maxitem = dim(Rthetasum)[1]
  V = as.list(1:maxitem)
  sigma1 = list()
  for (j in 1:(maxitem-1)) {
    #obtain submatrix
    temp = Rthetasum[unlist(V),unlist(V)]
    
    #find column sums of submatrix
    columnSums = colSums(temp)
    
    minindex = which.min(columnSums)
    
    #set selected item to be next item in sigma1
    sigma1[[j]] = V[[minindex]]
    
    #remove the selected item from V
    V = V[-minindex]
  }
  
  #after for loop, V has only 1 item left, so set this to be the last element of
  #sigma1
  sigma1[[maxitem]] = V[[1]]
  
  return(sigma1)
}

#full maximum likelihood estimation from S_N
#if theta_init or max_iterations not specified, they are set to the default 
#values in EstimateSigmaTheta
#(these values are described in the description of that function).
ML_estimation <- function(S_N,theta_init = NULL,max_iterations = NULL) {
  #compute Q (and also Pis)
  output = compute_Q(S_N)
  Q = output$Q
  Pis = output$Pis
  
  #compute q
  q = compute_q(Pis)
  
  #compute N
  N = computeN(S_N)
  
  #compute R
  R = compute_R(q,Q)
  
  #maximum likelihood estimation
  output1 = EstimateSigmaTheta(R,N,theta_init,max_iterations)
  sigma1 = output1$sigma1
  theta = output1$theta
  
  return(list(sigma1 = sigma1, theta=theta))
}

#sample from GMMs (generalized mallows S model)
#t is vector with the lengths we want the sampled pis to have, and its length 
#is the number of
#pis we want to sample
#theta is the vector of parameters. It must have length larger than max(t).
sample_GMM <- function(t,theta) {
  num_samples = length(t)
  maxrank = length(theta)
  
  #we don't need to sample maxrank every time, but sampling is fast so it's ok
  S = matrix(data =0,nrow = num_samples,ncol = maxrank)
  for (j in 1:maxrank) {
    S[,j] = rgeom(num_samples,prob = 1- exp(-theta[j]))
  }
  pis = list()
  for (i in 1:num_samples) {
    #print(S[i,1:t[i]])
    pis[[i]] = sj_to_pi(S[i,1:t[i]])
  }
  return(pis)
  
}

#given S_N with arbitrary items, we want to decide a way to map them to 1:n where n is the number of items
#returns list of 3 elements, first element is the renamed S_N (S_Nnew)
#second element is reordermap, a vector that can map the old names of items to the new ones
#third element is unordermap, a vector that can map the new names of items to their old ones
rename_SN <- function(S_N) {
  #first find the unique elements in S_N
  unique_items = unique(unlist(S_N))
  maxitem = max(unique_items)
  
  #build the mapping. The mapping we build is, put the items in the order they occurred in S_N.
  #Label these in order from 1 ... N.
  #mapping is a vector where, at index j, contains the new item name for item j.
  mapping = rep(0,maxitem)
  for (j in 1:maxitem) {
    #check where item j occurred in unique_items, if it occurred
    index = which(unique_items==j)
    
    #if the item did occur, then put in the index where it occurred.
    if (!identical(index, integer(0))) {
      mapping[j] = index
    }
  }
  
  #now that we have the mapping, we can build the new S_N
  S_Nnew = list()
  for (i in 1:length(S_N)) {
    S_Nnew[[i]] = list()
    for (j in 1:length(S_N[[i]])) {
      #print(S_N[[i]][[j]])
      S_Nnew[[i]][[j]] = mapping[S_N[[i]][[j]]]
    }
  }
  return(list(S_Nnew = S_Nnew, reordermap = mapping,unordermap = unique_items))
}

#using items that have been renamed using rename_SN, we can go back to the
#original items using unordermap
original_names_SN <- function(S_Nnew,unordermap) {
  S_Norg = list()
  for (i in 1:length(S_Nnew)) {
    S_Norg[[i]] = list()
    for (j in 1:length(S_Nnew[[i]])) {
      #print(S_N[[i]][[j]])
      S_Norg[[i]][[j]] = unordermap[S_Nnew[[i]][[j]]]
    }
  }
  return(S_Norg)
}

#this functions completes two rankings relative to one another. It returns the completed rankings, which
#have the same length and contain the same items.
complete_rankings <-function(pi1,pi2) {
  pi1vec = unlist(pi1)
  pi2vec = unlist(pi2)
  #first we find the elements in common between the two rankings
  intersection = intersect(pi1vec,pi2vec)
  
  #we will remove the elements in the intersection item by item
  #(this can be made faster)
  #start with the full rankings
  pi1vec_minus_intersection = pi1vec
  pi2vec_minus_intersection = pi2vec
  for (i in 1:length(intersection)) {
    item_to_remove = intersection[i]
    
    #modify pi1
    indices_to_remove = which(pi1vec_minus_intersection == item_to_remove)
    pi1vec_minus_intersection = pi1vec_minus_intersection[-indices_to_remove]
    
    #modify pi2
    indices_to_remove = which(pi2vec_minus_intersection == item_to_remove)
    pi2vec_minus_intersection = pi2vec_minus_intersection[-indices_to_remove]
  }
  
  #now we append pi1vec_minus_intersection onto the end of pi2vec, and
  #append pi2vec_minus_intersection onto the end of pi1vec, completing both rankings
  pi2vec_complete = c(pi2vec,pi1vec_minus_intersection)
  pi1vec_complete = c(pi1vec,pi2vec_minus_intersection)
  return(list(pi1_complete = as.list(pi1vec_complete),pi2_complete = as.list(pi2vec_complete)))
}

#computes the deltas (number of inversions at each rank) between two rankings of possibly different lengths
compute_deltas <-function(pi1,pi2) {
  #In order to compute the number of inversions at each rank between two rankings,
  #we first complete both lists so that they are the same length and
  #contain the same items.
  output = complete_rankings(pi1,pi2)
  pi1_complete = output$pi1_complete
  pi2_complete = output$pi2_complete
  
  #number of items in each completed list
  t = length(pi1_complete)
  #initialize C1_0 C2_0
  C1= list()
  C2 = list()
  #at each rank, we compute the number of inversions added, stored in deltas
  deltas = rep(0,t)
  for (j in 0:(t-1)) {
    print(j)
    item1 = pi1_complete[[j+1]]
    item2 = pi2_complete[[j+1]]
    
    #there are 5 different cases which determine how we compute the number of
    #inversions delta_{j+1} as well as how C1_j and C2_j are modified to form
    #C1_{j+1} and C2_{j+1}
    #which case we're in is determined by new1 and new2, binary indicators of
    #whether item1 in already in C2 and whether item2 is already in C1,
    #respectively.
    if (length(C1) == 0) {
      new1 = 1
    } else {
      #max in R acts as an OR statement
      new1 = 1-max((C2 == item1))
    }
    
    if (length(C2) ==0) {
      new2 = 1
    } else {
      #max in R acts as an OR statement
      new2 = 1-max((C1 == item2))
    }
    
    if (item1 == item2) {#Case 1: item1 = item2
      deltas[j+1] = length(C1) + length(C2)
      #C1 and C2 remain unchanged in case 1
      
      #In the remaining cases, item1 != item2 since we've ruled that out
    } else if ((new1==TRUE) & (new2==TRUE)) {#Case 2: item1 and item2 are both new
      deltas[j+1] = length(C1) + length(C2) + 1
      C1 = append(C1,item1)
      C2 = append(C2,item2)
    } else if ((new1==TRUE) & (new2==FALSE)) {#Case 3: item1 is new, item2 is old
      index = which(C1==item2)
      deltas[j+1] = index - 1 + length(C2)
      #remove item2 from C1, add item1
      C1 = append(C1[-index],item1)
      #C2 doesn't change
    } else if ((new1==FALSE) & (new2==TRUE)) {#Case 4: item2 is new, item1 is old
      index = which(C2==item1)
      deltas[j+1] = index - 1 + length(C1)
      #remove item1 from C2, add item2
      C2 = append(C2[-index],item2)
      #C1 doesn't change
    } else if ((new1==FALSE) & (new2==FALSE)) { #Case 5: item1 and item2 are both old
      index1 = which(C2==item1)
      index2 = which(C1 ==item2)
      deltas[j+1] = index1 + index2 - 2
      #remove item2 from C1, and remove item1 from C2
      C1 = C1[-index2]
      C2 = C2[-index1]
    }
  }
  return(deltas)
}