\name{miR.test}
\alias{miR.test}
\title{Main Function of miRtest package.}
\usage{miR.test(X, Y, A, group.miRNA, group.mRNA, gene.set.tests="romer",
    design.miRNA, design.mRNA, adjust="none", permutation=FALSE,
    nrot=1000, allocation.matrix=FALSE, verbose=FALSE, errors=TRUE)}
\description{Main Function of miRtest package.}
\author{Stephan Artmann}
\value{Matrix with testing results for every miRNA in its rows and the applied gene set test in its columns. Note that result will depend on whether multiple hypothesis testing correction was applied or not.}
\references{Artmann, Stephan and Jung, Klaus and Bleckmann, Annalen and Beissbarth, Tim (2012).
Detection of simultaneous group effects in microRNA expression and 
related functional gene sets. PLoS ONE 7(6):e38365, PMID: 22723856.

Brunner, E. (2009) Repeated measures under non-sphericity.
Proceedings of the 6th St. Petersburg Workshop on Simulation,
605-609.

Jelle J. Goeman, Sara A. van de Geer, Floor de Kort, Hans C. van
Houwelingen (2004) A global test for groups of genes: testing
association with a clinical outcome. Bioinformatics 20, 93-99.

Jung, Klaus and Becker, Benjamin and Brunner, Edgar and Beissbarth, Tim (2011).
Comparison of Global Tests for Functinoal Gene Sets in
Two-Group Designs and Selection of Potentially
Effect-causing Genes. Bioinformatics, 27: 1377-1383.

Majewski, IJ, Ritchie, ME, Phipson, B, Corbin, J, Pakusch, M,
Ebert, A, Busslinger, M, Koseki, H, Hu, Y, Smyth, GK, Alexander,
WS, Hilton, DJ, and Blewitt, ME (2010). Opposing roles of polycomb
repressive complexes in hematopoietic stem and progenitor cells.
_Blood_, published online 5 May 2010.

Mansmann, U. and Meister, R., 2005, Testing differential gene
expression in functional groups, _Methods Inf Med_ 44 (3).

Smyth, G. K. (2004). Linear models and empirical Bayes methods for
assessing differential expression in microarray experiments.
_Statistical Applications in Genetics and Molecular Biology_,
Volume *3*, Article 3.

Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader,
JE, and Smyth, GK (2010). ROAST: rotation gene set tests for
complex microarray experiments. _Bioinformatics_, published online
7 July 2010.}
\arguments{\item{X}{miRNA expression matrix with genes in rows and replicates in columns}
\item{Y}{mRNA expression matrix with genes in rows and replicates in columns}
\item{A}{Allocation data.frame or Allocation matrix. An allocation data.frame contains the mRNAs in its first column and the miRNAs in its second column. See vignette `miRtest' for information on Allocation matrices.}
\item{group.miRNA}{Vector of miRNA group membership, being either numeric or a factor (**this makes a difference**). E. g. if you have four replicates in a control group and three replicates in a treated group, you may choose c(1,1,1,1,2,2,2)}
\item{design.miRNA}{If specified, group.miRNA will be ignored. Here you can specify a design matrix as it is returned from the model.matrix `limma' function.}
\item{design.mRNA}{If specified, group.mRNA will be ignored. Here you can specify a design matrix as it is returned from the model.matrix `limma' function.}
\item{group.mRNA}{Vector of mRNA group membership, being either numeric or a factor (**this makes a difference**).E. g. if you have four replicates in a control group and three replicates in a treated group, you may choose c(1,1,1,1,2,2,2)}
\item{gene.set.tests}{Test to be applied for gene set testing. Can be one or more of the following: `globaltest', `GA', `RHD', `KS', `W', `Fisher', `roast', `romer', or `all' if you want to do all tests.}
\item{adjust}{Muliple hypothesis testing adjustment. Same options as in "p.adjust" function.}
\item{permutation}{Number of permutations for `globaltest' or `GlobalAncova' gene set tests. Put to "FALSE" to use the approximate p-values instead of permutation ones.}
\item{nrot}{Number of rotations for rotation tests `ROAST' and `romer'}
\item{allocation.matrix}{Logical, is A an allocation matrix with mRNAs in its columns and miRNAs in its rows, or is it an allocation data.frame?}
\item{verbose}{Defaults to FALSE. If TRUE, output on progress is printed.}
\item{errors}{Defaults to TRUE. If set to FALSE, some errors checking correct sizes of matrices are turned into warning messages.}}
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
