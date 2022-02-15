
# reading data from phases extracted by the Python script
# compute HMM score and hidden phases
# let Python paste detected and hidden phases into a CSV file 
# write HMM results file

library(reticulate)
library(hmm.discnp)

py_available(initialize=TRUE)
pyphase <- import_from_path("extractphasesfromjsonR", ".")
pygreat <- import_from_path("greatlogmakerR", ".")


trainhmm <- function(y) {
  m <- diag(5)  # create a "dirty" matrix as a start point for the learning. No zeros or ones.
  m[,]<-0.01
  diag(m) = 0.96
  yval <- c(1,2,3,4,5)  
  h <- hmm(y,yval, par0 = list(tpm=m, Rho=m), K=5, stationary=FALSE, cis=FALSE, newstyle=FALSE)
  return(h)
}

seqprob <- function(hmm,y) {
  # compute (kind of) probability of generating sequence y by hmm
  # first apply viterbi to attain most probable states
  s  <- viterbi(y,Rho=hmm$Rho, tpm=hmm$tpm, ispd = c(0.96,0.01,0.01,0.01,0.01))
  # then compute probability of output given these states
  prs <- pr(s,y,Rho=hmm$Rho, tpm=hmm$tpm, ispd = c(0.96,0.01,0.01,0.01,0.01))
  return(list(s=s, p=prs))
}

score.hmm <- function(hmm,seq) {
  result <- list(log.like = hmm$log.like)
  result$tpm.score = (sum(diag(hmm$tpm)) + sum(diag(hmm$tpm[,-1])))/5  # our measure for how close to the ideal (= all 1 diag)
  result$Rho.score = sum(diag(hmm$Rho))/5  # our measure for how close to the ideal (= all 1 diag)
  result$seq.prob = seqprob(hmm, seq)
  result$score = (result$tpm.score + result$Rho.score) / 2
  class(result) = "hmm.score"
  return(result)
} 


computehmm <- function(seq)  {
  h <- trainhmm(seq)
  sc <- score.hmm(h,seq)
  r <- list(seq = seq,
            tpmsc = sc$tpm.score,
            rhosc = sc$Rho.score,
            score = sc$score,
            pr = sc$seq.prob$p,
            vseq = sc$seq.prob$s,
            loglike = sc$log.like,
            len = length(h$y)
  )
  return(r)  
}

writeresult <- function(myfile, res) {
  options(digits=16)
  fileConn<-file(paste(myfile,"_hmmres.txt",sep=""))
  writeLines(c( paste("HMM results for ",myfile, "_003.json",sep=""),    
                paste("detect seq = ",paste(res$seq,collapse=","),sep=""),
                paste("hidden seq = ",paste(res$vseq,collapse=","),sep=""),
                paste("tmpsc = ",res$tpmsc,sep=""),
                paste("rhosc = ",res$rhosc,sep=""),
                paste("score = ",res$score,sep=""),
                paste("pr = ",res$pr,sep=""),
                paste("loglike = ",res$loglike,sep=""),
                paste("len = ",res$len,sep="")
                ),fileConn)
  close(fileConn)
}
 

 
myfile <- "test"

seqtext = pyphase$detectphases(myfile)
seq <- as.factor(as.numeric(strsplit(seqtext,",")[[1]]))

res = computehmm(seq)

hmmres = paste(res$vseq,collapse=",")
pyphase$insertphases(myfile,hmmres)

writeresult(myfile,res)

pygreat$makegreatlog("test_phase.csv","DP_test.csv","test_great.csv")
