##################
#Author: bramach2#
##################

#clear workspace
rm(list = ls(all=T))

#load required libraries
library(MASS)
library(ellipse)        

#set random seed
set.seed(1)

#set alpha value
alpha <- 0.0001
CV <- 1.8692

#simulate 2-d gaussian data with 4 centers
data<-mvrnorm(n=100, mu=c(5,5), Sigma = matrix(c(10,3,3,2),2,2))
data<-rbind(data,mvrnorm(n=200, mu=c(20,5), Sigma = matrix(c(5,2,2,8),2,2)))            
data<-rbind(data,mvrnorm(n=50, mu=c(5,20), Sigma = matrix(c(3,6,6,20),2,2)))
data<-rbind(data,mvrnorm(n=150, mu=c(20,20), Sigma = matrix(c(15,6,6,15),2,2)))
data<-data.frame(data)
names(data)<-c("X1","X2")

#load dataset provided
data<-read.table("hw5-3d-data.csv",sep = ",",header = T)

#generalize for n-dimensional data
p<-ncol(data)+1

#initialize one cluster membership for all data points
data[,p]<-rep(1,nrow(data))
names(data)[p]<-"cluster"

#initalize first center
first<-kmeans(data[,-p],centers = 1)
centers<-first$centers

#L2 norm calculation
norm_vec <- function(x) sqrt(sum(x^2))

#keep track of number of clusters
numclusters = 1

repeat{
        numcenters <- nrow(centers)
        done = 0
        for(j in 1:numcenters)
        {
                cj <- centers[j,]
                X <- subset(data, cluster == j)
                if(nrow(X) <= 7)
                {
                        next
                }
                E <- eigen(cov(as.matrix(X[,-p])))
                m <- E$vectors[,1]*(sqrt(2*E$values[1]/pi))
                newc1 <- centers[j,] + m
                newc2 <- centers[j,] - m
                km <- kmeans(X[,-p], rbind(newc1,newc2))
                centersdash <- km$centers
                newclusters = km$cluster
                v = centersdash[1,] - centersdash[2,]
                Xdash <- as.matrix(X[,-p])%*%as.matrix(v) / norm_vec(v)
                Xdash = (Xdash - mean(Xdash)) / sd(Xdash)
                sorted <- sort(Xdash)
                z <- pnorm(sorted)
                n = length(z)
                sum = 0
                for(i in 1:n)
                {
                        zi = z[i]
                        zni = z[n+1-i]
                        sum = sum + (2*i-1)*(log(zi) + log(1-zni))
                }
                Asq = -n - sum/n
                Asq_st = Asq*(1+ 4/n - 25/(n^2))
                if((Asq_st < -CV) | (Asq_st > CV))
                {
                        #reject the null hypothesis and accept new centers
                        done = 1 #to indicate that we are not yet done
                        numclusters = numclusters + 1
                        centers[j,] = centersdash[1,]
                        centers <- rbind(centers, centersdash[2,])
                        #got to be doubly careful about ordering of centers and cluster numbers
                        d1 = dist(rbind(X[which(newclusters == 1)[1],-p],centersdash[1,]))
                        d2 = dist(rbind(X[which(newclusters == 1)[1],-p],centersdash[2,]))
                        if(d1 < d2)
                        {
                                #means that cluster "1" has center centersdash[1]
                                newclusters[newclusters==2] = numclusters
                                newclusters[newclusters==1] = j
                        }else{
                                #means that cluster "1" has center centersdash[2]
                                newclusters[newclusters==1] = numclusters
                                newclusters[newclusters==2] = j
                        }
                        data[which(data$cluster == j),p] = newclusters
                        for(u in 1:(p-2))
                        {
                                plot(data[,c(u,u+1)], col = data$cluster, pch = '.', cex = 5)
                                #draw ellipses
                                for(l in 1:numclusters)
                                {
                                        x <- cov(data[which(data$cluster == l),c(u,u+1)])
                                        el <- ellipse(x, centre = centers[l,c(u,u+1)])
                                        lines(el, col = l, lwd = 2)
                                        points(x = centers[l,u], y = centers[l,u+1], pch = '+', cex = 2)
                                }
                                
                        }
                        
                }
        }
        if(done == 0)
        {
                break;
        }
}


