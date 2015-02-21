
## modified on 9 Feb, 2015
## copied and modified on 9 Dec, 2014
## written on 14 May, 2014
###
### To generate small synthetic network with hidden common cause

#
# to randomly generate (synthetic) gene network in the form of matrices,
# where each link a_ij has a number representing the
# effect of gene i on gene j: +ve for activation, -ve for inhibition.
# Associated with each a_ij =/= 0 is t_ij > 0, which represents the time delay
# of the effect of gene i on gene j.

var.transition.mat <- function(links, delays) {
  # the matrix A of the Vector Autoregressive Model, if the time series
  # x(t) = a_1*x(t-1) + a_2*x(t-2) + ... + a_d*x(t-d)
  # written as X(t) = A X(t-1), where X(t) = [x'(t) x'(t-1) ... x'(t-d)]',
  #         |a_1 a_2 ... a_d|
  #         | I   0  ...  0 |
  # and A = | 0   I  ...  0 |
  #         |        ...    |
  #         |          I  0 |

  n <- nrow(links);
  md <- max(delays);
  m <- matrix(0, nrow=md*n, ncol=md*n);

  # the a's
  for(i in 1:n) {
    for(j in 1:n) {
      if(delays[i,j] > 0) {
        m[j, i + (delays[i,j]-1)*n] <- links[i,j]; # links[i,j] is effect of gene i on gene j
      }
    }
  }
  # the I's
  for(i in 1:((md-1)*n)) {
    m[n+i,i] <- 1;
  }
  #
  m
}

max.abs.eigenvalue <- function(links, delays) {
  x <- eigen(var.transition.mat(links,delays));
  max(abs(x$values))
}

###
random.fill <- function(L, es) {
  ## fill in es for length L as much as possible, the rest are randomly picked
  if(L < length(es)) {
    sample(es,L, replace=FALSE)
  } else {
    r <- rep(es, L %/% length(es));
    n <- L %% length(es);
    if(n > 0) {
      r <- c(r, sample(es, n, replace=FALSE));
    }
    sample(r,L, replace=FALSE)
  }
}

gen.misc.hidden <- function(ng, nh, max.parents, max.delay) {
  # ng observed nodes (the first ng), nh hidden nodes.
  # ng should be large enough.
  # each hidden node will have 0,1,2 or 3 parents, and 2,3,4 or 5 children, all distinct.
  # the hidden node is placed as the last
  # Returns a list of two matrices of (ng+nh) by (ng+nh), one is the links, the other
  # is the delays.
  n <- ng + nh;
  rlink <- matrix(0, nrow=n, ncol=n);
  rdelay <- matrix(0, nrow=n, ncol=n);
  ps <- c(0,1,2,3);
  cs <- c(2,3,4,5);
  np <- random.fill(nh, ps); ## number of parents for each hidden node
  nc <- random.fill(nh, cs); ## number of children for each hidden node
  nps <- sum(np);
  ncs <- sum(nc);
  n.other <- ng - nps - ncs;
  mp <- sample(0:max.parents, ng - ncs, replace=TRUE); ## number of parents
  for(i in 1:(ng - ncs)) { ## fill in the columns one by one, (i,j) is i --> j
    if(i <= n.other) { ## others, can use any observed gene as parents
      idx <- sample(1:ng, mp[i]);
    } else { ## parents of hidden nodes, parents cannot be children of hidden nodes
      idx <- sample(1:(ng-ncs), mp[i]);
    }
    rlink[idx,i] <- (1-2*rbinom(mp[i],1,0.5))*runif(mp[i],0.5,1.5);
    rdelay[idx,i] <- sample(1:max.delay, mp[i], replace=TRUE);
  }
  ## the hidden nodes
  pj <- n.other;
  cj <- ng - ncs;
  for(i in 1:nh) {
    ## parents
    if(np[i] > 0) {
      idx <- pj + (1:np[i]);
      rlink[idx,ng+i] <- (1-2*rbinom(np[i],1,0.5))*runif(np[i],0.5,1.5);
      rdelay[idx,ng+i] <- sample(1:max.delay, np[i], replace=TRUE);
    }
    ## children, all nc > 0
    idx <- cj + (1:nc[i]);
    rlink[ng+i,idx] <- (1-2*rbinom(nc[i],1,0.5))*runif(nc[i],0.5,1.5);
    rdelay[ng+i,idx] <- sample(1:max.delay, nc[i], replace=TRUE);
    ##
    pj <- pj + np[i];
    cj <- cj + nc[i];
  }
  ##
  list(links=rlink, delays=rdelay)
}

permute_grn <- function(links, delays, ng=nrow(links)) {
  # to randomly permute the genes so that the position has no special advantage
  x <- 1:nrow(links);
  x[1:ng] <- sample(1:ng);
  list(links=links[x,x], delays=delays[x,x])
}

### example of The synthetic network: 
# a_ij is effect of gene i on gene j.
# mdl.links <- matrix(c(0    , 0.7 , 0    , 0.31, 0.8 ,
#                       0    , 0   , -0.85, 0   , 0   ,
#                       -0.55, 0.64, 0    , 0   , 0   ,
#                       0    , 0   , 0    , 0   , -0.6,
#                       0    , 0   , 0.4  , 0   , 0   ),
#                     nrow=5, ncol=5, byrow=TRUE);
# mdl.delay <- matrix(c(0, 1, 0, 4, 2,
#                       0, 0, 1, 0, 0,
#                       2, 1, 0, 0, 0,
#                       0, 0, 0, 0, 1,
#                       0, 0, 1, 0, 0),
#                     nrow=5, ncol=5, byrow=TRUE);
# 
#

sim_grn <- function(links, delays, N, init.mean=1, init.sigma2=1, alpha=1) {
  # links is a matrix encoding the links of a gene network (direction and magnitude)
  # delays contains the time delay (in steps) for each link
  # Returns a matrix a matrix of N by n, where n is the number of genes in the gene network.
  # The expression should be interpreted in log scale, since its value may be negative.
  # The initial expression is from N(init.mean, init.sigma2)
  # alpha controls the gaussianity of the noise, the noise is sign(e)*(abs(e)^alpha) where e is N(0, err.sigma2)
  # and the err.sigma2 of each noise will be heuristically determined to be delta*variance for that variable.
  T <- delays;
  maxd <- max(delays);
  ng <- nrow(links);
  r <- matrix(0,nrow=(N+maxd),ncol=ng);
  tmp.e <- rnorm((N+maxd)*ng, mean=0, sd=1);
  init.sd <- sqrt(init.sigma2);
  err <- matrix(sign(tmp.e)*(abs(tmp.e)^alpha), nrow=N+maxd, ncol=ng, byrow=TRUE);
  # re-scale each column such that its variance is init.sigma2
  for(i in 1:ng) {err[,i] <- init.sd*err[,i]/sd(err[,i]);}
  # generate the initial expression
  for(i in 1:maxd) {
    r[i,] <- init.mean + err[i,];
  }
  #
  for(i in (maxd+1):(maxd+N)) {
    for(j in 1:ng) {
      x <- err[i,j];
      for(k in 1:ng) {
        if(T[k,j] != 0) {
          x <- x + r[i-T[k,j],k]*links[k,j];
        }
      }
      r[i,j] <- x;
    }
  }
  r[(maxd+1):(maxd+N),]
}

#
plot_exp <- function(r) {
  # r is a n by g matrix, where n is the number of time points,
  # g is the number of genes
  # The values are the expression of the genes at the different time points
  n <- nrow(r);
  g <- ncol(r);
  legend.r <- if(is.null(colnames(r))) rep("",g) else colnames(r);
  for(i in 1:g) {legend.r[i] <- paste("i=",i,legend.r[i]);}
  plot(x=c(1,n),y=range(r), type="n", main="Expressions",xlab="Time",ylab="Expression");
  for(i in 1:g) {
    lines(x=1:n, y=r[,i], type="b",col=i,pch=i);
  }
  legend(x="topright",legend=legend.r,pch=1:g,col=1:g);
}

##
print.true.grn <- function(links, delays) {
  # print in a format acceptable to cmp_grn.py
  n <- nrow(links);
  # links[i,j] is effect of gene i on gene j
  # 0-based indices
  for(j in 1:n) {
    for(i in 1:n) {
      if(delays[i,j] > 0) {
        cat("To: ",j-1, " From: ",i-1, " Delay: ",delays[i,j],
            " Coef: ",links[i,j], "\n", sep="");
      }
    }
  }
}

print.exp <- function(e) {
  # print expression data in tab separated file
  for(i in 1:nrow(e)) {
    cat(e[i,], sep="\t");
    cat("\n");
  }
}

##
gen.one.misc.case <- function(ng,nh,mp,md, e2, alpha, nps, outprefix) {
  ## ng observed genes, nh hidden nodes
  tmp.tn.o <- gen.misc.hidden(ng,nh, mp,md);
  # scale it so that the eigen value of VAR model has modulus
  # will not be too large
  e <- max.abs.eigenvalue(tmp.tn.o$links, tmp.tn.o$delays);
  if(e >= 1) {
    tmp.tn.o$links <- tmp.tn.o$links * 0.5 / e;
  }
  
  tmp.tn <- permute_grn(tmp.tn.o$links, tmp.tn.o$delays, ng);
  ## leave out the last column, which is the hidden node
  tmp.sr <- sim_grn(tmp.tn$links, tmp.tn$delays, max(nps), 1,e2, alpha);

  ##
  sink(paste(outprefix,"_grn.txt", sep=""));
  print.true.grn(tmp.tn$links, tmp.tn$delays);
  sink();
  
  for(np in nps) {
    sink(paste(outprefix,"_nps",np,".txt", sep=""));
    print.exp(tmp.sr[1:np,1:ng]);
    sink();
    ##
    sink(paste(outprefix,"_nps",np,"_complete.txt", sep=""));
    print.exp(tmp.sr[1:np,]);
    sink();
  }
}

gen.misc.cases <- function() {
  out <- "misc_cases/";
  md <- 4;
  nps <- 800; # need only generate the large number that we will use, then take fragments if we want to test shorter series
  dir.create(out, mode="0755");
  for(ng in c(50,100)) { # number of observed genes
    ## with 10% hidden nodes
    nh <- ng/10;
    dir.create(paste(out,"n",ng, sep=""), recursive=TRUE, mode="0755");
    for(e2 in c(0.5,1,2,3,4)) {
      for(alpha in c(0.5, 1, 2, 3)) {
        for(r in 1:40) { # replicates
          name <- paste("n",ng,"e",e2,"a",10*alpha,"r",r, sep="");
          cat("Processing", name, "\n");
          gen.one.misc.case(ng,nh,4,md, e2,alpha, c(20,50,100,200,400,800), paste(out,"n",ng,"/",name, sep=""));
        }
      }
    }
  }
}

###
read.one.grn <- function(ng, name) {
  ## since the GRN may not have connection for all genes, so the number of genes is given.
  ## read GRN in this format: (0-based)
  ## To: 0 From: 2 Delay: 1 Coef: -0.908147 Score: 75.699439
  ## To: 1 From: 2 Delay: 1 Coef: 0.848796 Score: 55.699844
  ## To: 2 From: 1 Delay: 4 Coef: -0.184466 Score: 2.015631
  tmp <- try(read.table(name));
  rL <- matrix(0, nrow=ng,ncol=ng);
  rD <- matrix(0, nrow=ng,ncol=ng);
  if(! inherits(tmp, 'try-error')) {
    for(i in 1:nrow(tmp)) {
      from <- 1 + tmp[i,4];
      to <- 1 + tmp[i,2];
      if(tmp[i,6] > 0) {
        rL[from,to] <- tmp[i,8];
        rD[from,to] <- tmp[i,6];
      }
    }
  }
  list(links=rL, delays=rD)
}

###
## generate segment data for the misc cases, using the same GRN
gen.one.seg.misc.case <- function(ng, e2, alpha, segs, links,delays, outprefix) {
  ## the first ng are observed genes, the rest are hidden

  ##
  sink(paste(outprefix,"_grn.txt", sep=""));
  print.true.grn(links, delays);
  sink();

  for(r in 1:length(segs)) {
    tmp.sr <- sim_grn(links, delays, segs[r], 1,e2, alpha);

    sink(paste(outprefix,"_s",r,".txt", sep=""));
    print.exp(tmp.sr[,1:ng]);
    sink();

    sink(paste(outprefix,"_s",r,"_c.txt", sep=""));
    print.exp(tmp.sr);
    sink();
  }
}

gen.seg.misc.cases <- function() {
  out <- "seg_misc_cases/";
  indir <- "misc_cases";
  md <- 4;

  dir.create(out, mode="0755");
  for(ng in c(50,100)) { # number of observed genes
    ## with 10% hidden nodes
    nh <- ng/10;
    dir.create(paste(out,"n",ng, sep=""), recursive=TRUE, mode="0755");
    ## simulate the segments, 32 segments, each with length uniformly
    ## chosen in [20 .. 30] to get a total of about 640 to 960 points
    segs <- sample(20:30, 32, replace=TRUE);
    for(e2 in c(0.5,1,2,3,4)) {
      for(alpha in c(0.5, 1, 2, 3)) {
        dir.create(paste(out,"n",ng,"/e",e2,"a",10*alpha, sep=""), recursive=TRUE, mode="0755");
        for(r in 1:40) { # replicates
          name <- paste("n",ng,"e",e2,"a",10*alpha,"r",r, sep="");
          cat("Processing", name, "\n");
          ## use the existing GRN from misc_cases
          tmp.grn <- read.one.grn(ng+nh, paste(indir, "/n",ng,"/",name, "_grn.txt", sep=""));
          if(sum(tmp.grn$delays > 0) == 0) {
            cat("Read empty GRN for ", name, "\n");
          } else {
            gen.one.seg.misc.case(ng, e2,alpha, segs, tmp.grn$links,tmp.grn$delays, paste(out,"n",ng,"/e",e2,"a",10*alpha,"/r",r, sep=""));
          }
        }
      }
    }
  }
}


###
