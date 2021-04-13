#data is the data table stemming from a bismarck coverage file and should be in the form: 
#start | end | seqnames | spin 
#with the meanings chromosome, start and end position, methylation state
#annot a data table  the annotation file for a genomic region and  should be in the form seqnames(chromosome) start end
library(data.table)
library(magrittr)

#compute.overlaps returns the overlap of data with and annotation data table in bed format
compute.overlaps <- function(data, annot){
  
  # compute the overalps between the coverage and the annotation
  start0 <- data[,region := foverlaps(data,setkey(annot,seqnames,start,end), which=T, mult="first", nomatch=NA)] %>%
    
    # remove base pairs that are not in the annotation
    .[!is.na(region)] %>%
    
    # define the central position for each sequence of the annotation
    .[, center := (annot[.BY$region, end]-annot[.BY$region, start])/2,by=region] %>%
    .[, .SD[which.min((start-center)^2),], by=region] %>%
    .[, .(seqnames, start, end)]
  
  return(start0)
}

#avg.corr.fun.annot computes correlation functions for a set of genomic annotations in bed format

avg.corr.fun.annot <- function(data, annot){
  start0 <- data[,region := foverlaps(data,setkey(annot,seqnames,start,end), which=T, mult="first", nomatch=NA)] %>%
    .[!is.na(region)] %>%
    .[, center := (annot[.BY$region, end]-annot[.BY$region, start])/2,by=region] %>%
    .[, .SD[which.min((start-center)^2),], by=region] %>%
    .[, .(seqnames, start, end)]
  
  return(avg.corr.fun(data, setkey(start0,seqnames,start, end)))
}


#######################################
# Helper functions
#######################################
#avg.corr.fun is the first step for computing correlation functions, start is the result of the overlaps
#between the selected annotation and the coverage.

avg.corr.fun <- function(data, start){
  #Create a matrix large enough to contain all the two points spatial correlation functions (mat)
  mat <- matrix(data = -2,1e9,4)
  return(get_avg_corr(mat, data[order(seqnames,start)], start[order(seqnames,start)]))
  #GO to get_avg_corr
}


raw.corr.fun <- function(data, start){
  mat <- matrix(data = -2,1e9,4)
  return(corr_wrapper2(mat, data,start))
}


# C function for computing two point spatial correlations

cppFunction('List spatial_correlation(NumericMatrix correlation,NumericVector x, IntegerVector position, IntegerVector position0){
           
           
           /* insize it is the the length in base pairs of the chromosome of the coverage */
            unsigned int insize = x.size();
            
            /* position0sizeit is the the length in base pairs of the chromosome of the annotation */
            unsigned int position0size = position0.size();
            
            /* max_length is the maximum distance you want to compute correlation functions up to */
            unsigned int max_length = 10000;        
            
            /* set initial counters */
            
            unsigned int center = 0;
            unsigned int length = 0;
            unsigned int i = 0;
            unsigned int j= 0;
            unsigned int counter =0;
            unsigned int iposition0=0;
            counter = center;
            
            
            /* do a loop over all the base pairs*/
            while(counter<insize){
            
            /* compute two point correlations up to max_length */
            if(length<max_length){
            
            /* control to break the loops in case the input empty matrix was too small */
            if(j>correlation.size()-1){
            break;
            }
            
            /* compute the value (methylation at position center)* (methylation at position center +i) */
            correlation(j,0) = x[center]*x[center+i];
            
            /* find the  distance in base pairs between the two sites*/
            correlation(j,1) = position[center+i] - position[center];
            
            /* single sites values*/
            correlation(j,2) = x[center];
            correlation(j,3) = x[center+i];
            
            /* assing the real distance to length*/
            length = position[center+i] - position[center];
            i += 1;
            j += 1;
            counter += 1;}
            
            /* After the loop exceed max_length you find the next position of the annotation*/
            else{
              iposition0++;
              
              /* Control to quit the function if you exceed the last positionn of the annotation*/
              if(iposition0>position0size){
              printf("size methylation0 exceed");
              break;
               }
                  
              while(position[center] < position0[iposition0]){
                center++;
                if(center>insize){
              printf("size methylation exceed");
              break;
              }  
              }
            counter = center;
            i = 0;
            length=0;
            } 
            
            
            }
            
            /* return the filled matrix with correlations, last is the nrows of the matrix */
            
            List ret;
            ret["corr"] = correlation;
            ret["last"] = j;
            return ret;
            }
            ')


#corr_wrapper_2 computes the correlation functions for each annotationn and each chromosomes
#dt is a data table containing CpG methylation information as above

corr_wrapper2 <- function(mat, dt, start0){
  print("computing correlation...")
  
  #corr.tmp is a list containing correlation functionns for each chromosome, go to corr_wrapper
  corr.tmp <- dt[,list(list(corr=corr_wrapper(mat, spin,start,start0[seqnames==.BY$seqnames,start]))) ,by=seqnames]
  
  # tmp is a list of data tablesc created from the filled matrix of correlation functions 
  tmp <- lapply(seq_along(corr.tmp$seqnames), 
                function(x) return(as.data.table(corr.tmp$V1[[x]])[,
                                                                   seqnames:=corr.tmp[x,seqnames]]))
  return(rbindlist(tmp) %>% setnames(.,c("si.sj","distance","i","j","seqnames")))
}

#corr_wrapper is the last step. It takes as an input the matrix creates, spin(methylation) the start of each base pair a
# and the start0, which is the central position of each annotation

corr_wrapper <- function(mat, spin, start, start0){
  #go to the C function spatial_correlation
  tmp <- spatial_correlation(mat, spin,start,start0)
  
  return(tmp$corr[1:tmp$last,])
}


# GET_AVG_CORR returns tmp, which is the result of the function corr_wrapper_2
get_avg_corr <- function(mat, dt,start0){
  tmp <- corr_wrapper2(mat, dt,start0)
  return(tmp)
}

