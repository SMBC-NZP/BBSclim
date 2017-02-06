#' write_one_matrix
#'
#' Function to write .dm file

write_one_matrix <- function(dmnum,dm) {
  mat1=NULL
  if (is.null(dm)) {
    mat1=sprintf('%d 0 0',dmnum)
  } else {
    nr=dim(dm)[1]; nc=dim(dm)[2]; if (is.null(colnames(dm))) nc=0
    mat1=paste(dmnum,nr+1,nc+1,sep=' ')   #  write design matrix no. rows,cols
    mat1=paste(mat1,paste(c('-',colnames(dm)),collapse=','),sep='\n')
    for (i in 1:dim(dm)[1])
      mat1=paste(mat1,paste(c(rownames(dm)[i],dm[i,]),collapse=','),sep='\n')
  }
  return(mat1)
}

#' parse_output
#'
#' Function to parse Presence .out file

parse_output <- function(outname) {
  output=readLines(outname)
  if (length(output)<4) {
    cat('\nerror - nothing in output file\n\n')
    return(NULL)
  } else {
    i=grep('^[ABCDEDG][0123456789]+ ',output)
    v=t(matrix(unlist(strsplit(gsub(' +',',',output[i]),',')),nrow=5))
    betas=as.numeric(v[,4]); names(betas)=v[,2]; beta.se=as.numeric(v[,5])
    i=grep('^-2log.like',output);  negll2=as.numeric(strsplit(output[i],'=')[[1]][2])
    i=grep('^AIC ',output);  aic=as.numeric(strsplit(output[i],'=')[[1]][2])
    i=grep('^Number of param',output); npar=as.numeric(strsplit(output[i],'=')[[1]][2])
    i=grep('^==>name=',output); modname=gsub('.+name=','',output[i])
    dmx=list()
    for (i in 1:6) {
      j=grep(paste('^Matrix',i),output);
      nr=as.numeric(gsub('.+=','',gsub(',.+','',output[j])))
      nc=as.numeric(gsub('.+cols=','',output[j]))
      dm=NULL;
      if (nr>0 && nc>1) {
        dm=matrix(unlist(strsplit(gsub(' +',',',output[j+2:nr]),',')),nr-1,nc,byrow=T)
        rnames=dm[,1]; cnames=paste0(c('a','b','c','d','e','f')[i],1:(nc-1));
        dm=dm[,-1]; dim(dm)=c(nr-1,nc-1);
        rownames(dm)=rnames; colnames(dm)=cnames;
      }
      dmx[[i]]=dm
    }
    j=grep('of <',output); l=length(j); real.est=real.se=real.name=NULL
    if (l>0) {
      j1=which(nchar(output[j[l]:length(output)])==0)[1]+j[l]-1
      j=c(j,j1);
      for (i in 1:l) {
        i1=j[i]; i2=j[i+1]
        k=grep(':',output[i1:i2])
        for (m in 1:length(k)) {
          v=strsplit(gsub(' +',',',gsub('.+: ','',output[i1+k[m]-1])),',')[[1]][2:3]
          real.est=c(real.est,as.numeric(v[1])); real.se=c(real.se,as.numeric(v[2]))
          real.name=c(real.name,gsub(' +:.+','',output[i1+k[m]-1]))
        }
      }
      names(real.est)=real.name; names(real.se)=real.name
    }
    return(list(outfile=outname,modname=modname,betas=betas,beta.se=beta.se,
                aic=aic,negll2=negll2,npar=npar,dm=dmx,
                real.est=real.est,real.se=real.se))
  }
}


#' write_dm_and_run2
#'
#' Write design matrices to file and run correlated detection model in Presence. Note that this function calls Presence in a slightly different way than the `write_dm_and_run()` function in the current version of RPresence
#' @param pao .pao file for analysis
#' @param cov_list List containing the covariates for each parameter
#' @param out Which output folder should results be written to?
#' @param dm_list List containing the design matrices for each parameter
#' @param modname Optional name for .out file with model results
#' @param fixed vector of fixed values (eg., fixed(1)=0, fixed(4)=.1 -> fixed=c(0,.1); names(fixed)=c(1,4))
#' @param initvals vector of beta initial values
#' @param maxfn maximum number if iterations for optimization routine
#' @return A .out file containing the fitted model and parameters estimates (see ?occ.mod for details)
#' @export

write_dm_and_run2 <- function(pao, cov_list, is.het,
                              dm_list = NULL, modname,
                              fixed = FALSE, inits = FALSE,
                              maxfn = 32000, alpha = NULL, parse = FALSE) {

  if(inits){
    total.betas	<- sum(1 + length(cov_list$psi.cov),
                       1 + length(cov_list$th0.cov),
                       1 + length(cov_list$th1.cov),
                       1 + length(cov_list$gam.cov),
                       1 + length(cov_list$eps.cov),
                       (1 + is.het) * (1 + length(cov_list$p1.cov)),
                       is.het)  # last two are for p2 and pi
    initvals	<- rep(0, total.betas)
  }else{
    initvals <- NULL
  }


  if(fixed) {
    fixedpars <- rep("eq", pao$nseasons);
    r1 <- dim(dm_list$dm1)[1] + dim(dm_list$dm2)[1] + dim(dm_list$dm3)[1] + dim(dm_list$dm4)[1]
    names(fixedpars) <- as.character((r1 + 1):(r1 + pao$nseasons))
  }else{
    fixedpars = NULL
  }

  outname <- paste('inst/output/', alpha, "/pres/", modname,'.out',sep='')

  if (file.exists(outname)) cat('\n**** output file exists - model not run\n***********\n') else {
   dmname = paste0('inst/output/', alpha,"/pres/", modname, '_dm.dm')
    if (!is.null(dm_list$dm1)) {
      if (file.exists(dmname)) file.remove(dmname)
      dmmat = write_one_matrix(0, dm_list$dm1)
      dmmat = paste(dmmat, write_one_matrix(1,dm_list$dm2), sep='\n')
      dmmat = paste(dmmat, write_one_matrix(2,dm_list$dm3), sep='\n')
      dmmat = paste(dmmat, write_one_matrix(3,dm_list$dm4), sep='\n')
      dmmat = paste(dmmat, write_one_matrix(4,dm_list$dm5), sep='\n')
      dmmat = paste(dmmat, write_one_matrix(5,dm_list$dm6), sep='\n')
      if (!is.null(fixedpars)) {
        sfx=paste('fixed(',names(fixedpars),')=',fixedpars,sep='')
        dmmat=paste(c(dmmat,sfx),collapse='\n')
      }
      if (!is.null(initvals)) {
        si=paste('init(',1:length(initvals),')=',initvals,sep='',collapse='\n')
        dmmat=paste(dmmat,si,sep='\n')
      }
      cat(dmmat,'\n',sep='',file=dmname)
    }
    if (length(find("PRES_FLDR")) <= 0) {
      if (file.exists("C://Progra~2/Presence")) {
        PRES_FLDR = "C://Progra~2/Presence"
      }
      else {
        PRES_FLDR = "C://Program Files/Presence"
      }
    }
   s = c(paste0("i=", pao$paoname), paste0("j=", dmname), paste0("l=",
       outname), "VC", "quiet", paste0("maxfn=", maxfn), paste0("name=",
       modname))

    t1 = Sys.time(); cat(c(s,'\n'));
    i = .C("_Z8presencePiPPc", as.integer(length(s)), as.character(s), PACKAGE = "RPresence")
    t2 = Sys.time(); cat(c('\nCPU time:',t2-t1,'\n'));
    #if (rand.dm) unlink(dmname)
  }
  v=1; if (parse) v=parse_output(outname)
  return(v)
}
