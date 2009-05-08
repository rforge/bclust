bclust<-function(x,transformed.par,x.id=1:nrow(x),effect.family="gaussian",
var.select=TRUE,labels=NULL)
{
if (missing(x)){stop("data are not specified")}
if (!is.matrix(x)){stop("data must be a matrix")}
if (missing(transformed.par)){stop("transformed.par is missing")}
if (!is.vector(transformed.par)){stop("transformed.par must be a vector")}
if (is.null(transformed.par)){stop("model parameters are not specified")}
if (!(length(x.id)==nrow(x))){stop("x.id mismatches")}
if ((sum(is.na(x))>0) | (sum(is.na(x.id))>0) | 
(sum(is.na(transformed.par))>0)) stop("NA is not allowed")
if(!(effect.family%in%c("gaussian","alaplace"))) {stop(
"just gaussian and alaplace families are allowed")}
id.order<-order(x.id)
repno<-as.vector(table(x.id[id.order]))
if (is.null(labels)){if (is.null(rownames(x))){labels<-paste(1:length(repno))} else  
{if (length(repno)==nrow(x)){labels<-rownames(x)}} } else 
{if (!(length(labels)==length(repno))){warning("label length do not match the data")}}

y<-x[id.order,]
labels<-labels[id.order]

bclustG<-function(y,repno,hyperparameters)
{
out<-.C("RfastbclustG",PACKAGE="bclust",y=as.double(y),nrowy=as.integer(nrow(y)),ncoly=as.integer(ncol(y)),repno=as.double(repno),
nrepno=as.integer(length(repno)),theta=as.double(hyperparameters),merge=as.double(rep(0,2*(length(repno)-1))),
height=as.double(rep(0,(length(repno)-1))))
mytree<-list(merge=matrix(out$merge,ncol=2,byrow=TRUE),height=out$height,order=nocrossing.order(matrix(out$merge,ncol=2,byrow=TRUE))
)
return(mytree)
}

bclustvsG<-function(y,repno,hyperparameters)
{
out<-.C("RfastbclustvsG",PACKAGE="bclust",y=as.double(y),nrowy=as.integer(nrow(y)),ncoly=as.integer(ncol(y)),repno=as.double(repno),
nrepno=as.integer(length(repno)),theta=as.double(hyperparameters),merge=as.double(rep(0,2*(length(repno)-1))),
height=as.double(rep(0,(length(repno)-1))))
mytree<-list(merge=matrix(out$merge,ncol=2,byrow=TRUE),height=out$height,order=nocrossing.order(matrix(out$merge,ncol=2,byrow=TRUE))
)
return(mytree)
}

bclustAL<-function(y,repno,hyperparameters)
{
out<-.C("RfastbclustAL",PACKAGE="bclust",y=as.double(y),nrowy=as.integer(nrow(y)),ncoly=as.integer(ncol(y)),repno=as.double(repno),
nrepno=as.integer(length(repno)),theta=as.double(hyperparameters),merge=as.double(rep(0,2*(length(repno)-1))),
height=as.double(rep(0,(length(repno)-1))))
mytree<-list(merge=matrix(out$merge,ncol=2,byrow=TRUE),height=out$height,order=nocrossing.order(matrix(out$merge,ncol=2,byrow=TRUE))
)
return(mytree)
}


bclustvsAL<-function(y,repno,hyperparameters)
{
out<-.C("RfastbclustvsAL",PACKAGE="bclust",y=as.double(y),nrowy=as.integer(nrow(y)),ncoly=as.integer(ncol(y)),repno=as.double(repno),
nrepno=as.integer(length(repno)),theta=as.double(hyperparameters),merge=as.double(rep(0,2*(length(repno)-1))),
height=as.double(rep(0,(length(repno)-1))))
mytree<-list(merge=matrix(out$merge,ncol=2,byrow=TRUE),
height=out$height,order=nocrossing.order(matrix(out$merge,ncol=2,byrow=TRUE)))
return(mytree)
}

nocrossing.order<-function(mergemat)
{
	myorder<-matrix(NA,nrow(mergemat),nrow(mergemat)+1)
			na.add <- function(vec)
	{
		return(c(vec,rep(NA,nrow(mergemat)+1-length(vec))))
	}

	myorder[1,]<-na.add(mergemat[1,])
			i<-2
			while (i<=nrow(mergemat))
	{
		helpvec1<-mergemat[i,][1]
				helpvec2<-mergemat[i,][2]
				if (helpvec1>0)
		{
			helpvec1<-na.exclude(myorder[helpvec1,])
		}
		if (helpvec2>0)
		{
			helpvec2<-na.exclude(myorder[helpvec2,])
		}
		myorder[i,]<-na.add(c(helpvec1,helpvec2))
				i<-i+1
	}
	return(-myorder[i-1,])
}

if (effect.family=="gaussian")
{
if (var.select)
{
if (!(length(transformed.par)==6)) {stop("transformed.par is of a wrong size")}
	{
	btree<-bclustvsG(y,repno,transformed.par)
	}
}
if (!var.select)
{
if (!(length(transformed.par)==5)) {stop("transformed.par is of a wrong size")}
	{
	btree<-bclustG(y,repno,transformed.par)
	}
}
}


if (effect.family=="alaplace")
{
if (var.select)
{
if (!(length(transformed.par)==7)) {stop("transformed.par is of a wrong size")}
	{
	btree<-bclustvsAL(y,repno,transformed.par)
	}
}
if (!var.select)
{
if (!(length(transformed.par)==6)) {stop("transformed.par is of a wrong size")}
	{
	btree<-bclustAL(y,repno,transformed.par)
	}
}
}
if (sum(is.na(btree$height))>0) {stop("dendrogram heights includes NAs, transformed.par values may be inadequately adjusted")}
bclust.tree<-btree
bclust.tree$logposterior<- -btree$height
minheight<-min(btree$height)
minindex<-which(btree$height==minheight)
increment<-diff(btree$height,lag=1)
bclust.tree$height<-diffinv(abs(increment),lag=1)
bclust.tree$clust.number<-((length(repno)-1):1)
bclust.tree$cut<-bclust.tree$height[minindex]
bclust.tree$data<-y
bclust.tree$repno<-repno
bclust.tree$transformed.par<-transformed.par
bclust.tree$labels<-labels
oldClass(bclust.tree)<-c("bclustvs","hclust")
return(bclust.tree)
}



loglikelihood<-function(x.mean,x.css,
repno,transformed.par,effect.family="gaussian",var.select=TRUE)
{
n<-length(as.vector(x.mean))
if (effect.family=="gaussian")
{
	if(var.select)
	{
	return(-.C("RlogmargtmvsG",PACKAGE="bclust",
 	theta=as.double(transformed.par), 
 	mean= as.double(x.mean), 
  	css=as.double(x.css),
  	nrow=as.integer(nrow(x.css)),
   	ncol=as.integer(ncol(x.css)),r=as.double(repno),
    	nr=as.integer(length(repno)), result=as.double(0))$result)
	} else
	{
	return(-.C("RlogmargtmG",PACKAGE="bclust",
 	theta=as.double(transformed.par), 
   	mean= as.double(x.mean),
   	css=as.double(x.css), n=as.integer(n),
   	r=as.double(repno), nr=as.integer(length(repno)),
   	result=as.double(0))$result)
	}
}
if (effect.family=="alaplace")
{
	if(var.select)
	{
	return(-.C("RlogmargtmvsAL",PACKAGE="bclust",
 	theta=as.double(transformed.par), 
  	mean= as.double(as.vector(x.mean)),
  	css=as.double(as.vector(x.css)), 
	nrow=as.integer(nrow(x.css)), ncol=as.integer(ncol(x.css)),
  	r=as.double(repno), nr=as.integer(length(repno)),
  	result=as.double(0))$result)
	} else
	{
	return(-.C("RlogmargtmAL",PACKAGE="bclust",
 	theta=as.double(transformed.par), 
 	mean= as.double(x.mean),
 	css=as.double(x.css), n=as.integer(n), r=as.double(repno),
  	nr=as.integer(length(repno)), result=as.double(0))$result)
	}
 stop("function arguments mismatches")	
}
}



meancss<-function(x,x.id=1:nrow(x))
{
id.order<-order(x.id)
repno<-as.vector(table(x.id[id.order]))
y<-x[id.order,]
gendiag<-function(x.vec)
{
y<-rep(sum(x.vec),2*length(x.vec)-1)
y[2*1:length(x.vec)-1]<-x.vec
matrix(rep( rep(1:0,length=2*length(x.vec)-1),y),sum(x.vec),length(x.vec))
}
matsum<-function(x.mat,rowsby) {t(t(x.mat)%*%gendiag(rowsby))}
matrep<-function(x.mat,eachrow)
{
result<-c()
for (i in 1:nrow(x.mat))
result<-c(result,rep(x.mat[i,],eachrow[i]))
return(matrix(result,ncol=ncol(x.mat),byrow=TRUE))
}

R<-matrix(rep(repno,ncol(y)),ncol=ncol(y))
x.mean<-matsum(y,repno)/R
x.css<-matsum((y-matrep(x.mean,repno))^2,repno) 
return(list(mean=x.mean,css=x.css,repno=repno))
}



