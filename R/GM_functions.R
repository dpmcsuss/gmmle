#' @importFrom ggplot2 aes as_labeller facet_grid geom_errorbar geom_line geom_point ggplot scale_color_manual scale_linetype_manual scale_shape_manual scale_x_log10 theme_bw xlab ylab
#' @importFrom magrittr `%>%`
#' @importFrom igraph get.adjacency gorder sample_gnm sample_gnp sample_pa  sample_smallworld

are_you_sure <- function(ask){
  if (ask) {
    answer <- utils::askYesNo(
      paste0("WARNING: This simulation can take several hours. ",
      "Are you sure you want to run it?"), default = FALSE)
    if (is.na(answer) || !answer) {
      message("Cancelling simulation.")
      return(FALSE)
    }
  }
  return(TRUE)
}


#' @title Statistics for errors induced by shuffling
#'
#' @description Calculate edge disagreements and probability bounds
#'   after k vertices are randomly shuffled in a graph G
#' @param k number of shuffled vertices
#' @param G graph object (igraph)
#' @param m number of random permutations
#' 
#'
#' @return A vector of length 4. The first three elements
#'  are the mean, sd, and minimum for the number of edge
#'  disagreements. The 4th, pmax is a theoretical bound
#'  on maximum error for which graph matching will be
#'  possible.
#' 
#' @examples
#'  kshuffle_mean_sd(10, karate, m = 10)
#' @export
kshuffle_mean_sd <- function(k, G, m = 1000) {
  n <- gorder(G)
  rand_perm_error <- function(G, k) {
    flag <- FALSE
    i <- 0
    while(!flag) {
      permutation <- 1:n
      permute <- sample(n, k, replace = F)
      permutation[permute] <-
        permute[sample(1:k, k, replace = F)]
      flag <- sum(permutation != 1:n) == k
      i <- i + 1
    }
    A <- get.adjacency(G)
    sum(abs(A - A[permutation, permutation]))/2
  }
  res <- sapply(1:m, function(i) rand_perm_error(G, k))
  pmax <- 0.5 - 0.5*sqrt(1 - exp( -6*k*log(n) / min(res)))
  tibble::tibble(K = k, mean = mean(res), sd = stats::sd(res),
    min = min(res), pmax = pmax)
}

#' @title Shuffling Simulation of different graph models
#'
#' @description Run a simulation instance of different graph models
#'    by shuffling vertices and computing disagreement statistics.
#'    Models include Erdos Renyi (ER), Watts-Strogatz small world
#'    graphs with rewiring probabilities of 0.05 and 0.7, and
#'    preferential attachment (PA) models with powers of 1, 1, and 2,
#'    and zero.appeal of 0, 500, and 0. All graphs are undirected with
#'    n=500 vertices. Each graph is sampled once and the shuffles are
#'    repeated m times.
#'
#' @param k_grid sequence of values for number of shuffled vertices
#' @param d average degree of the random graph
#' seed = random seed
#' @param seed Seed for random number generator
#' @param n number of nodes
#' @param m number of Monte Carlo replicates
#' @param ask Ask about running big simulations.
#'
#' @return Tibble with rows for each graph-K pair and
#'  7 columns, one for the model, parameters, and K.
#'  Columns 4:6 give the mean, sd, and min for the
#'  number of edge disagreements. Column 7 gives
#'  a theoretical estimate for the maximum probability
#'  of edge flips for which graph matching will still
#'  be possible.
#' 
#' @examples
#'  simulation_instance(k_grid = c(2, 10),
#'    m = 10, n = 30, ask = FALSE)
#' 
#' @export
simulation_instance <- function(k_grid = c(2, 10, 50, 100, 200, 400),
  d = 1, seed = 0, m = 1000, n = 500, ask = TRUE) {

  set.seed(seed)
  sim_df <- tibble::tribble(
    ~model, ~G, ~param,
    "ER", sample_gnm(n, m = d*n, directed = F),
      list(m = d*n, n=n),
    "WS", sample_smallworld(1, n, d, 0.05),
      list(d = d, p = 0.05, dim = 1),
    "WS", sample_smallworld(1, n, d, 0.7),
      list(d = d, p = 0.7, dim = 1),
    "PA0", sample_pa(n, power = 1, m = d,
        directed = F, zero.appeal = 0),
      list(power = 1, m = d, zero.appeal = 0),
    "PA1", sample_pa(n, power = 1, m = d,
        directed = F, zero.appeal = 500),
      list(power = 1, m = d, zero.appeal = 500),
    "PA2", sample_pa(n, power = 2, m = d,
        directed = F, zero.appeal = 0),
      list(power = 2, m = d, zero.appeal = 0)
    )

  res <- sim_df %>% dplyr::mutate(res = purrr::map(G,
                                                   ~purrr::map_df(k_grid, kshuffle_mean_sd, G = .x, m = m)
  )) %>%
    dplyr::select(-G) %>%
    tidyr::unnest(res)
  df <- res[, c(4, 5, 6, 7, 1, 1, 3)]
  names(df) <- c("mean", "sd", "min", "p.max", "model", "param", "K")
  df$modparam <- interaction( df$param,df$model)
  df$param <- as.factor(df$param)
  df

}


#' @title Shuffling Simulation for ER graphs
#'
#' @description Run a simulation instance of different graph models
#'    by shuffling vertices and computing disagreement statistics.
#'    Models include Erdos Renyi (ER), Watts-Strogatz small world
#'    graphs with rewiring probabilities of 0.05 and 0.7, and
#'    preferential attachment (PA) models with powers of 1, 1, and 2,
#'    and zero.appeal of 0, 500, and 0. All graphs are undirected with
#'    n=500 vertices. Each graph is sampled once and the shuffles are
#'    repeated m times.
#'
#' @param k_grid sequence of values for number of shuffled vertices
#' @param seed Seed for random number generator
#' @param m number of Monte Carlo replicates
#' @param n number of nodes
#' @param p_vec vectors of edge probabilities
#' @param ask Ask about running big simulations.
#'
#' @return Tibble with rows for each graph-K pair and
#'  7 columns, one for the model, parameters (list column), and K.
#'  Columns 4:6 give the mean, sd, and min for the
#'  number of edge disagreements. Column 7 gives
#'  a theoretical estimate for the maximum probability
#'  of edge flips for which graph matching will still
#'  be possible.
#' simulation_instance_ER(k_grid = c(2, 10),
#'   p_vec = (1:2)/10, n = 30, m = 10, seed = 0, ask = FALSE)
#' 
#' @export
simulation_instance_ER <- function(
  k_grid = c(2, 10, 50, 100, 200, 400),
  p_vec = (1:5)/10, n = 500, m = 1000,
  seed = 0, ask = TRUE) {
  set.seed(seed)
  res <- tibble::tibble(p = p_vec) %>%
    dplyr::mutate(G = purrr::map(p, ~sample_gnp(n, .x, directed = F))) %>%
    dplyr::mutate(res = purrr::map(G,
                                   ~purrr::map_df(k_grid, kshuffle_mean_sd, G = .x, m = m)
    )) %>%
    dplyr::select(-G) %>%
    tidyr::unnest(res)
  df <- res[, c(3, 4, 5, 6)]
  df$model <- res$p
  df$param <- "ER"
  df$K <- res$K
  df$modparam <- interaction( df$param,df$model)
  df$param <- as.factor(df$param)
  names(df)[1:4] <- c("mean", "sd", "min", "p.max")
  return(df)
}

#' @title Shuffling Simulation for for Watts-Strogatz graphs
#'
#' @description Run a simulation instance of different graph models
#'    by shuffling vertices and computing disagreement statistics.
#'    Models include Erdos Renyi (ER), Watts-Strogatz small world
#'    graphs with rewiring probabilities of 0.05 and 0.7, and
#'    preferential attachment (PA) models with powers of 1, 1, and 2,
#'    and zero.appeal of 0, 500, and 0. All graphs are undirected with
#'    n=500 vertices. Each graph is sampled once and the shuffles are
#'    repeated m times.
#'
#' @param k_grid sequence of values for number of shuffled vertices
#' @param seed random seed
#' @param seed Seed for random number generator
#' @param m number of Monte Carlo replicates
#' @param n number of nodes
#' @param d average degree
#' @param p_vec vectors of edge probabilities
#' @param ask Ask about running big simulations.
#'
#' @return Tibble with rows for each graph-K pair and
#'  7 columns, one for the model, parameters, and K.
#'  Columns 4:6 give the mean, sd, and min for the
#'  number of edge disagreements. Column 7 gives
#'  a theoretical estimate for the maximum probability
#'  of edge flips for which graph matching will still
#'  be possible.
#' 
#' @examples
#' simulation_instance_WS(k_grid = c(2, 10),
#'  p_vec = c(0, .05), d = 10, n = 50, m = 10, ask = FALSE)
#' 
#' @export
simulation_instance_WS <- function(
  k_grid = c(2, 10, 50, 100, 200, 400),
  p_vec = c(0, .05, .1, .25, 1), d = 50, n = 500, m = 1000,
  seed = 0, ask = TRUE) {
  set.seed(seed)
  res <- tibble::tibble(p = p_vec) %>%
    dplyr::mutate(G = purrr::map(p, 
                                 ~sample_smallworld(
                                   dim = 1, size = n, nei = d, p = .x)
    )
    ) %>%
    dplyr::mutate(res = purrr::map(G,
                                   ~purrr::map_df(k_grid, kshuffle_mean_sd, G = .x, m = m)
    )) %>%
    dplyr::select(-G) %>%
    tidyr::unnest(res)
  df <- res[, c(3, 4, 5, 6)]
  df$model <- res$p
  df$param <- "WS"
  df$K <- res$K
  df$modparam <- interaction( df$param,df$model)
  df$param <- as.factor(df$param)
  names(df)[1:4] <- c("mean", "sd", "min", "p.max")
  return(df)
}

#' Generate corrupted graph from A with probability p
#' 
#' Each edge/non-edge is flipped with probability p.
#' 
#' @param A adjacency matrix
#' @param p corrupting probability
#' 
#' @return Corrupted graph
#' 
#' @export
corrupting_channel <- function(A, p) {
  n <- ncol(A)
  X <- get.adjacency(sample_gnp(n, p, directed = F))
  B <- A * (1-X) + (1-A) * X
  diag(B) = 0
  sigma <- sample(1:n, n)
  P <- Matrix::Diagonal(n)[sigma,]
  return(list(B = B, P = P))
}

#' Vertex matching errors for a corrupted graph
#' 
#' The graph A is corrupted with probability p and then
#'  the corrupted graph and the original graphs are matched
#'  with initialization at the identity.
#' @param A adjacency
#' @param p corrupting probability
#' @param seed random seed
#' 
#' @return list of diagonal of match and the proportion of errors
#' 
#' @examples
#' local_gm_errs(karate[], .05)
#' 
#' @export
local_gm_errs <- function(A, p, seed = 0) {
  set.seed(seed)
  n <- ncol(A)
  Bc <- corrupting_channel(A, p)
  B <- Bc$B
  if(n > 10000){
    # for arxiv network, set max_iter = 3
    match1 <- iGraphMatch::graph_match_FW(A, B,
      start = Matrix::Diagonal(n),
      max_iter = 3)
  } else {
    match1 <- iGraphMatch::graph_match_FW(A, B,
      start = Matrix::Diagonal(n))
  }

  res <- list(Pdiag = Matrix::diag(match1$P),
    errors =  sum(Matrix::diag(match1$P))/n)
  return(res)
}

#' Cumulative matching errors by degree
#' 
#' Calculate cumulative vertex matching errors ordered by degree (decreasing)
#' 
#' @param A adjacency matrix
#' @param errors_list list with number of errors, each entry is the result of a simulation
#' @param p_grid Grid of error probabilities
#' 
#' @return data frame with three columns for corrupting probability,
#'  cumulative errors, and degree rank.
highdegree_errors <-  function(A, errors_list, p_grid) {
  degree <- Matrix::colSums(A)
  n <- ncol(A)
  ord <- order(degree, decreasing = T)
  pdiag_Table <- lapply(1:length(p_grid), function(i_p) {
    sapply(errors_list, function(err_sim) err_sim[[i_p]]$Pdiag)
  })
  highdegree_tab <- lapply(pdiag_Table, function(pdiag_t) {
    apply(pdiag_t, 2, function(pdiag_t_i) {
      sapply(1:n, function(i)  {
        sum(pdiag_t_i[ord[1:i]])/i
      })
    })
  })
  highdegree_df <-
    Reduce(rbind, lapply(1:length(p_grid),
      function(p_i)
        data.frame(p.channel = p_grid[p_i],
                   errors.by.degree = c(highdegree_tab[[p_i]]),
                   degree_rank = rep((1:n)/n, length(errors_list) )
                 ))[1:6])

  return(highdegree_df)
}
