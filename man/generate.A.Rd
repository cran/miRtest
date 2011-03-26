\name{generate.A}
\alias{generate.A}
\title{Turn a data...}
\usage{generate.A(df, X, Y, verbose=TRUE)}
\description{Turn a data.frame indicating gene sets into the allocation matrix.}
\value{Allocation matrix A necessary for "miR.test" function.}
\author{Stephan Artmann}
\arguments{\item{df}{data.frame with mRNAs in its first and miRNAs in its second column.}
\item{X}{Expression matrix of miRNAs whose row names will be used to generate the list of miRNAs.}
\item{Y}{Expression matrix of mRNAs whose row names will be used to generate the list of mRNAs.}
\item{verbose}{Logical. Shall progress be printed?}}
\examples{#######################################
### Generate random expression data ###
#######################################
# Generate random miRNA expression data of 3 miRNAs
# with 8 replicates
set.seed(1)
X = rnorm(24);
dim(X) = c(3,8);
rownames(X) = 1:3;
# Generate random mRNA expression data with 20 mRNAs
# and 10 replicates
Y = rnorm(200);
dim(Y) = c(20,10);
rownames(Y) = 1:20;
# Let's assume that we want to compare 2 miRNA groups, each of 4 replicates:
group.miRNA = factor(c(1,1,1,1,2,2,2,2));
# ... and that the corresponding mRNA experiments had 5 replicates in each group
group.mRNA = factor(c(1,1,1,1,1,2,2,2,2,2));
####################
### Perform Test ###
####################
library(miRtest)
#Let miRNA 1 attack mRNAs 1 to 9 and miRNA 2 attack mRNAs 10 to 17.
# mRNAs 18 to 20 are not attacked. miRNA 3 has no gene set.
miR = c(rep(1,9),c(rep(2,8)));
mRNAs = 1:17;
A = data.frame(mRNAs,miR); # Note that the miRNAs MUST be in the second column!
A
set.seed(1)
P = miR.test(X,Y,A,group.miRNA,group.mRNA)
P


#####################################################
### For a faster result: use other gene set tests ###
#####################################################
# Wilcoxon two-sample test is recommended for fast results
# Note that results may vary depending on how much genes correlate

P.gsWilcox = miR.test(X,Y,A,group.miRNA,group.mRNA,gene.set.tests="W")
P.gsWilcox
############################################
### We can use an allocation matrix as A ###
############################################
A = generate.A(A,X=X,Y=Y,verbose=FALSE);
A
# Now we can test as before
set.seed(1)
P = miR.test(X,Y,A,group.miRNA,group.mRNA,allocation.matrix=TRUE)
P


#####################
### Other Designs ###
#####################

# Some more complicated designs are implemented, check the vignette "miRtest" for details.
group.miRNA = 1:8
group.mRNA = 1:10
covariable.miRNA = factor(c(1,2,3,4,1,2,3,4))    ### A covariable in miRNAs.
covariable.mRNA = factor(c(1,2,3,4,5,1,2,3,4,5)) ### A covariable in mRNAs.

library(limma)
design.miRNA = model.matrix(~group.miRNA + covariable.miRNA)
design.mRNA =  model.matrix(~group.mRNA + covariable.mRNA)

P = miR.test(X,Y,A,design.miRNA=design.miRNA,design.mRNA=design.mRNA,allocation.matrix=TRUE)
P}
