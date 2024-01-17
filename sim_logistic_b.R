#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

setwd("~/stan_halcyon")
set.seed(1066)

NSIM = 30 # 30 # simulations to fit
DNAL = 120 # alignment length

Q <- matrix(c(-1.075162, 0.186970, 0.696268, 0.191923, 0.181082, -0.873473, 0.255492, 0.436899, 0.674340, 0.255493, -1.164645, 0.234813, 0.191924, 0.451108, 0.242449, -0.885481),
            nrow=4, ncol=4, byrow=T)

PI <- c(0.246000, 0.254000, 0.254000, 0.246000)

null_tree = ape::read.tree("T128.txt")
null_tree = reorder(null_tree, "postorder")

E = nrow(null_tree$edge)
L=length(null_tree$tip.label)
N=2*L-1
R=null_tree$edge[E,1] # If postorder, then last node is root and going up the list is preorder.

# Simulate from range of scaling.
symetric_sim_r = function(n) { exp(rnorm(n, mean=0, sd=1)) }

cp_scale_tree = function(null_tree, r) {
  scale_tree = null_tree
  scale_tree$edge.length <- r * scale_tree$edge.length
  scale_tree
}

sig_fn = function(h, a, b, log_r) {
  h / (1 + exp(-(a + b*log_r)))
}

halcyon_sim = function() {
  # Scale the input tree.
  sim_r = symetric_sim_r(E)
  cne_tree = cp_scale_tree(null_tree, sim_r)
  # Simulate alignment with phangorn.
  sim_X = phangorn::simSeq(cne_tree, l=DNAL, Q=Q, bf=PI, type="DNA")
  sim_X = do.call(rbind, sim_X)
  # Simulate trait with ape.
  sim_h = runif(1, 0.1, 4)
  sim_a = rnorm(1, 0, 2)
  sim_b = 0 # No effect, intercept only.
  sim_yR= runif(1,-1, 1)
  y_tree = cp_scale_tree(null_tree, sig_fn(sim_h, sim_a, sim_b, log(sim_r)))
  sim_y = ape::rTraitCont(y_tree, sigma=1)
  sim_y = sim_y + sim_yR
  list(sim_h=sim_h, sim_a=sim_a, sim_b=sim_b, sim_r=sim_r, sim_X=sim_X, sim_y=sim_y, sim_yR=sim_yR)
}

the_sims = lapply(1:NSIM, function(i) { halcyon_sim() })

save.image(file="sim_logistic_b.RData")

# Make a path to leaf matrix. It is N nodes x E edges with N,E set if edge E is on path to node N.

make_ne_mat = function(tree, N, E) {
  ne_mat = matrix(0, nrow=N, ncol=E)
  for (e in seq(E, 1, -1)) {
    i = tree$edge[e,1]
    j = tree$edge[e,2]
    ne_mat[j,]  = ne_mat[i,]
    ne_mat[j,e] = 1
  }
  ne_mat
}

path_to_leaf=make_ne_mat(null_tree, N, E)

get_ith_sim = function(ith) {
  the_sim = the_sims[[ith]]
  sim_X = the_sim$sim_X
  sim_y = the_sim$sim_y
  the_dat <- list(L=L,N=N,E=E,R=R,S=ncol(sim_X),
                  X=sim_X,
                  Q=Q, pi_distro=PI,
                  edge=null_tree$edge,
                  edge_length=null_tree$edge.length,
                  path_to_leaf=path_to_leaf,
                  y=sim_y)
}

for (i in 1:NSIM) {
  sim_i = get_ith_sim(i)
  cmdstanr::write_stan_json(sim_i, sprintf("sim_logistic_b/sim_logistic_b_%d.json", i))
}
