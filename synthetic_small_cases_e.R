
## modified on 9 Feb, 2015
## copied and modified on 20 Sep, 2014
## written on 14 May, 2014
###
### To generate small synthetic network with hidden common cause

#
# to randomly generate (synthetic) gene network in the form of matrices,
# where each link a_ij has a number representing the
# effect of gene i on gene j: +ve for activation, -ve for inhibition.
# Associated with each a_ij =/= 0 is t_ij > 0, which represents the time delay
# of the effect of gene i on gene j.

gen.small.hidden <- function(m,n, max.delay) {
  # m nodes as parent, of a hidden node, then n children of the hidden node
  # the hidden node is placed as the last
  # Returns a list of two matrices of (m+n+1) by (m+n+1), one is the links, the other
  # is the delays.
  ng <- m+n+1;
  rlink <- matrix(0, nrow=ng, ncol=ng);
  rdelay <- matrix(0, nrow=ng, ncol=ng);
  if(m > 0) {
    rlink[1:m,ng] <- (1-2*rbinom(m,1,0.5))*runif(m,0.5,1.5);
    rdelay[1:m,ng] <- sample(1:max.delay, m, replace=TRUE);
  }
  rlink[ng,m+1:n] <- (1-2*rbinom(n,1,0.5))*runif(n,0.5,1.5);
  rdelay[ng,m+1:n] <- sample(1:max.delay, n, replace=TRUE);
  
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

# use relative error level.
cal.sigma2 <- function(A, delta) {
  # A[i,j] is effect of x_i -> x_j.
  # Let H_i = var(x_i), E_i = var(e_i), where x_i = \sum_k {a_{ki}*x_k + e_i}
  # then H_i ~= \sum_k {a_{ki}^2 * H_k + E_i} since e_i's are independent of others, and ignore the covariance among x_k
  # want E_i = \delta * H_i if x_i has no parents, and E_i = 1 if x_i has parents
  # Let B_{ij} = A_{ji}^2, H = [H_1 H_2 ... H_n]', E = [E_1 E_2 ... E_n]
  #  then H = BH + E, and E = CH + V,
  #  where C_{ij} = 0 if i \= j, C_{ii} = \delta if x_i has parents, C_{ii} = 0 if x_i has no parents;
  #  and V_i = 0 if x_i has parents, V_i = 1 if x_i has no parents.
  # So we have (I - B - C)H = V
  # To solve for H and then get E, as an estimate of the initial variance of errors to use.
  n <- nrow(A);
  # D is to be I - B - C.
  D <- t(-A^2);
  V <- ifelse(rowSums(D) >= 0, 1,0);
  for(i in 1:n) {D[i,i] <- D[i,i] + 1 - delta*(1-V[i]);}
  H <- solve(D,V);
  E <- delta*(1-V)*H + V;
  E
}

sim_grn <- function(links, delays, N, init.mean=1, init.sigma2=1, alpha=1) {
  # links is a matrix encoding the links of a gene network (direction and magnitude)
  # delays contains the time delay (in steps) for each link
  # Returns a matrix a matrix of N by n, where n is the number of genes in the gene network.
  # The expression should be interpreted in log scale, since its value may be negative.
  # The initial expression is from N(init.mean, init.sigma2)
  # alpha controls the gaussianity of the noise, the noise is sign(e)*(abs(e)^alpha) where e is N(0, err.sigma2)
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

gen.one.hidden.case <- function(m,n,md, e2, alpha, nps, outprefix) {
  ## m parents, n children, with one hidden node between them
  ng <- m + n; # number of observed genes
  tmp.tn.o <- gen.small.hidden(m,n, md);
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

gen.hidden.cases <- function() {
  out <- "small_cases_e/";
  md <- 4;
  nps <- 200; # need only generate the large number that we will use, then take fragments if we want to test shorter series
  dir.create(out, mode="0755");
  for(m in c(0,1,2,3)) { # number of parents
    for(n in c(2,3,4,5)) { # number of children
      ## with one hidden node between them
      ng <- m + n; ## number of observed genes
      dir.create(paste(out,"p",m,"c",n, sep=""), recursive=TRUE, mode="0755");
      for(e2 in c(0.5,1,2,3,4)) {
        for(alpha in c(0.5, 1, 2, 3)) {
          for(r in 1:20) { # replicates
            name <- paste("p",m,"c",n,"e",e2,"a",10*alpha,"r",r, sep="");
            cat("Processing", name, "\n");
            gen.one.hidden.case(m,n,md, e2,alpha, c(20,50,100,200), paste(out,"p",m,"c",n,"/",name, sep=""));
          }
        }
      }
    }
  }
}
##
## to generate expression data without hidden nodes, using output GRN from clinde on data with hidden nodes
gen.one.non.hidden.case <- function(md, e2, alpha, nps, links,delays, outprefix) {
  ng <- nrow(links); # number of observed genes
  tmp.sr <- sim_grn(links, delays, max(nps), 1,e2, alpha);

  ##
  sink(paste(outprefix,"_grn.txt", sep=""));
  print.true.grn(links, delays);
  sink();
  
  for(np in nps) {
    sink(paste(outprefix,"_nps",np,".txt", sep=""));
    print.exp(tmp.sr[1:np,1:ng]);
    sink();
  }
}

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

gen.non.hidden.cases <- function(indir) {
  ## organized the same way as for the hidden node data
  ## indir is the name of directory containing the GRN by clinde that we want
  ## we use the GRN from clinde on st2 and nps200
  out <- "small_cases_e_nh/";
  md <- 4;
  nps <- 200; # need only generate the large number that we will use, then take fragments if we want to test shorter series
  dir.create(out, mode="0755");
  for(m in c(0,1,2,3)) { # number of parents
    for(n in c(2,3,4,5)) { # number of children
      ng <- m + n; ## number of observed genes
      dir.create(paste(out,"p",m,"c",n, sep=""), recursive=TRUE, mode="0755");
      for(e2 in c(0.5,1,2,3,4)) {
        for(alpha in c(0.5, 1, 2, 3)) {
          for(r in 1:20) { # replicates
            name <- paste("p",m,"c",n,"e",e2,"a",10*alpha,"r",r, sep="");
            cat("Processing", name, "\n");
            tmp.grn <- read.one.grn(ng, paste(indir,"/p",m,"c",n,"/",name, "_nps200_st2_hiddenCL_outgrn_s2.txt", sep=""));
            ## try not to have empty GRN
            if(sum(tmp.grn$delays > 0) == 0) {
              tmp.grn <- gen.small.hidden(m,n-1,md);
            }
            gen.one.non.hidden.case(md, e2,alpha, c(20,50,100,200), tmp.grn$links,tmp.grn$delays, paste(out,"p",m,"c",n,"/",name, sep=""));
          }
        }
      }
    }
  }
}

###
## to generate segment data for the same GRNs for hidden and non-hidden cases.

gen.one.seg.hidden.case <- function(ng, e2, alpha, segs, links,delays, outprefix) {
  ## the first ng are observed genes, the rest are hidden
  ##
  sink(paste(outprefix,"_grn.txt", sep=""));
  print.true.grn(links, delays);
  sink();

  for(r in 1:length(segs)) {
    tmp.sr <- sim_grn(links, delays, segs[r], 1,e2, alpha);

    sink(paste(outprefix,"_seg",r,".txt", sep=""));
    print.exp(tmp.sr[,1:ng]);
    sink();

    sink(paste(outprefix,"_seg",r,"_complete.txt", sep=""));
    print.exp(tmp.sr);
    sink();
  }
}

gen.seg.hidden.cases <- function() {
  out <- "seg_small_cases_e/";
  indir <- "small_cases_e";
  md <- 4;

  dir.create(out, mode="0755");
  for(m in c(0,1,2,3)) { # number of parents
    for(n in c(2,3,4,5)) { # number of children
      ## with one hidden node between them
      ng <- m + n; ## number of observed genes
      dir.create(paste(out,"p",m,"c",n, sep=""), recursive=TRUE, mode="0755");

      ## simulate the segments, 8 segments, each with length uniformly
      ## chosen in [20 .. 30] to get a total of about 160 to 240 points
      segs <- sample(20:30, 8, replace=TRUE);
      for(e2 in c(0.5,1,2,3,4)) {
        for(alpha in c(0.5, 1, 2, 3)) {
          for(r in 1:20) { # replicates
            name <- paste("p",m,"c",n,"e",e2,"a",10*alpha,"r",r, sep="");
            cat("Processing", name, "\n");
            ## use the existing GRN from small_cases_e
            tmp.grn <- read.one.grn(ng+1, paste(indir, "/p",m,"c",n,"/",name, "_grn.txt", sep=""));
            if(sum(tmp.grn$delays > 0) == 0) {
              cat("Read empty GRN for ", name, "\n");
            } else {
              gen.one.seg.hidden.case(ng,e2,alpha, segs, tmp.grn$links,tmp.grn$delays, paste(out,"p",m,"c",n,"/",name, sep=""));
            }
          }
        }
      }
    }
  }
}

### for the non-hidden cases, also generate segment data
gen.one.seg.non.hidden.case <- function(ng, e2, alpha, segs, links,delays, outprefix) {
  ## the first ng are observed genes, the rest are hidden
  ##
  sink(paste(outprefix,"_grn.txt", sep=""));
  print.true.grn(links, delays);
  sink();

  for(r in 1:length(segs)) {
    tmp.sr <- sim_grn(links, delays, segs[r], 1,e2, alpha);

    sink(paste(outprefix,"_seg",r,".txt", sep=""));
    print.exp(tmp.sr[,1:ng]);
    sink();
  }
}

gen.seg.non.hidden.cases <- function() {
  out <- "seg_small_cases_e_nh/";
  indir <- "small_cases_e_nh";
  md <- 4;

  dir.create(out, mode="0755");
  for(m in c(0,1,2,3)) { # number of parents
    for(n in c(2,3,4,5)) { # number of children
      ## with one hidden node between them
      ng <- m + n; ## number of observed genes
      dir.create(paste(out,"p",m,"c",n, sep=""), recursive=TRUE, mode="0755");

      ## simulate the segments, 8 segments, each with length uniformly
      ## chosen in [20 .. 30] to get a total of about 160 to 240 points
      segs <- sample(20:30, 8, replace=TRUE);
      for(e2 in c(0.5,1,2,3,4)) {
        for(alpha in c(0.5, 1, 2, 3)) {
          for(r in 1:20) { # replicates
            name <- paste("p",m,"c",n,"e",e2,"a",10*alpha,"r",r, sep="");
            cat("Processing", name, "\n");
            ## use the existing GRN from small_cases_e
            tmp.grn <- read.one.grn(ng, paste(indir, "/p",m,"c",n,"/",name, "_grn.txt", sep=""));
            if(sum(tmp.grn$delays > 0) == 0) {
              cat("Read empty GRN for ", name, "\n");
            } else {
              gen.one.seg.non.hidden.case(ng,e2,alpha, segs, tmp.grn$links,tmp.grn$delays, paste(out,"p",m,"c",n,"/",name, sep=""));
            }
          }
        }
      }
    }
  }
}


###
