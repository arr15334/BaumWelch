
#' Create sequence from HMM
#'
#' Create a sequence from a simulated HMM given transition and emission matrices
#' @param transition.mat Transition Matrix specifying the probability of changing or remaining in the same state
#' @param emission.mat Emission matrix specifying probability of emitting certain symbol from a state
#' @param n The length of the sequence to be generated
#' @param start.prob Vector of initial probabilities, probaility of starting from a state
#' @return list containing generated sequence and the vector of states under which each symbol was generated
#' @export
create_sequence <- function(transition.mat, emission.mat, n, start.prob) {
  states <- c('H', 'L')
  observation.set <- c('A','C','G','T')
  sequence <- vector()
  h.states <- vector()
  h.states[1] <- sample(states, 1, prob = start.prob)
  sequence[1] <- sample(observation.set, 1, prob = emission.mat[h.states[1],])

  for (i in 2:n) {
    h.states[i] <- sample(states, 1, prob = transition.mat[h.states[i-1],])
    sequence[i] <- sample(observation.set, 1, prob = emission.mat[h.states[i],])
  }
  return(list('sequence'=sequence, 'hidden_states'=h.states))
}

#' Forward algorithm
#'
#' Calculates the probability of the next emission, using the computed probability of the current emission
#' @param obs Vector of visible observations
#' @param transition Transition matrix, columns must sum 1
#' @param E Emission matrix, rows must sum 1
#' @return alpha matrix of probabilities of being in a certain state for each step
#' @export
forward <- function(obs, transition, E) {
  x <- c(obs[1],obs)
  A <- transition
  #print(dim(A))
  ftable <- matrix(0, nrow=length(x), ncol=dim(A)[1])
  ftable[1,] <- 0
  ftable[1,1] <- 1
  for (i in 2:length(x)) {
    for (l in 1:dim(A)[1]) {
      K <- ftable[i-1,] * A[,l]
      K <- sum(K)
      ftable[i,l] <- E[l, x[i]] * K
    }
  }
  p <- sum(ftable[length(x),])
  alpha <- ftable[2:length(x),]
  return(alpha)
}


#' Backward Algorithm
#'
#' Calculates the probability of the sequence ending given state i at time t
#' @param obs Vector of visible observations
#' @param transition Transition Matrix
#' @param emission Emission matrix
#' @return beta matrix of probabilities of model being in state s at time t
#' @export
backwd <- function (obs, transition, emission) {
  Tn = length(obs)
  M = nrow(transition)
  beta = matrix(1, Tn, M)

  for(t in (Tn-1):1){
    tmp = as.matrix(beta[t+1, ] * emission[, obs[t+1]])
    beta[t, ] = t(transition %*% tmp)
  }
  return(beta)
}

#' Baum-Welch algorithm
#'
#' Updates transition and emission matrices using forward-backward algorithm
#' @param visible Vector of visible observations
#' @param t.mat Transition Matrix
#' @param e.mat Emission matrix
#' @param states The names of the hidden states
#' @param symbols The array of symbols that can be produced by the model
#' @return list containing a: transition matrix, b: emission matrix, initial_distribution vector
#' @export
BaumWelch = function(visible, t.mat, e.mat, initial_distribution,
                     states,symbols,n.iter = 100){
  A.diff = vector()
  B.diff = vector()
  for(i in 1:n.iter){
    A.old = t.mat
    B.old = e.mat
    N = length(visible)
    M = nrow(t.mat)
    K = ncol(e.mat)
    alpha = forward(visible, t.mat, e.mat)
    beta = backwd(visible,t.mat, e.mat)
    xi = array(0, dim=c(M, M, N-1)) #3D matrix
    for(t in 1:N-1){
      denominator = ((alpha[t,] %*% t.mat) * e.mat[,visible[t+1]]) %*% matrix(beta[t+1,])
      for(s in 1:M){
        numerator = alpha[t,s] * t.mat[s,] * e.mat[,visible[t+1]] * beta[t+1,]
        xi[s,,t]=numerator/as.vector(denominator)
      }
    }

    xi.all.t = rowSums(xi, dims = 2)
    t.mat = xi.all.t/rowSums(xi.all.t)
    gamma = apply(xi, c(1, 3), sum)
    gamma = cbind(gamma, colSums(xi[, , N-1]))
    for(l in 1:K){
      if(is.null(colnames(e.mat))) {
        e.mat[, l] = rowSums(gamma[, which(visible==l)])
      } else {
        e.mat[, colnames(e.mat)[l]] = rowSums(gamma[, which(visible == colnames(e.mat)[l])])
      }
    }
    e.mat = e.mat/rowSums(e.mat)
    A.diff[i] = sum(abs(A.old - t.mat))
    B.diff[i] = sum(abs(B.old - e.mat))
  }
  row.names(t.mat) <- colnames(t.mat) <- row.names(e.mat) <- states
  colnames(e.mat) <- symbols
  return(list(a = t.mat, b = e.mat,
              initial_distribution = initial_distribution,
              A.diff=A.diff, B.diff=B.diff))
}


