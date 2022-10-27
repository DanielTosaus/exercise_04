library("Biostrings")
#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param seq1 DNAString object representing NT or AA sequence to align
#' @param seq2 DNAString object representing NT or AA sequence to align
#' @param align A list of DNAString objects with alignment of input sequences
#' @param match An integer value of a score for matching bases
#' @param mismatch An integer value of score for mismatching bases
#' @param gap An integer value of penalty for gap insertion
hirschberg_template <- function(seq1, seq2, align, match, mismatch, gap){
    
    first_align_row <- align[[1]] # initialize the first row of alignment
    second_align_row <- align[[2]] # initialize the second row of alignment
  
    
    if(length(seq1)==0) # length of seq1 is equal to zero
    {
        for(i in (1:length(seq2))) # for each character in seq2
        {
            first_align_row <- c(first_align_row,'-') # add gap
            second_align_row <- c(second_align_row, seq2[i]) # add character from seq2
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if(length(seq2)==0) # length of seq2 is equal to zero
    {
        for(i in (1:length(seq1))) # for each character in seq1
        {
            first_align_row <- c(first_align_row, seq1[i]) # add character from seq1
            second_align_row <- c(second_align_row,'-') # add gap
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if(length(seq1)==1 & length(seq2)==1) # length of seq1 and seq2 is equal to 1
    {
        first_align_row <- seq1[1]# add character from seq1
        second_align_row <- seq2[1]# add character from seq2
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else  #NW algorithm
    {
        
        x_len <- length(seq1)   # length of seq1
        x_mid <- x_len/2        # half of the length of seq1
        y_len <- length(seq2)   # length of seq2
        
        left_score <- NWScore(seq1[1:x_mid], seq2)# NW score for the first half of seq1 and the whole seq2
        right_score <-NWScore(reverse(seq1[(x_mid+1):x_len]), reverse(seq2)) # NW score for the second half of seq1 and the whole seq2 (both are reversed)
        print(right_score)
        print(right_score)
        sum <- left_score + rev(right_score)
        y_mid <- which.max(sum)-1  # index of division for seq2
        print(y_mid)
        # The first half
        if(ymid==0) # index of division for seq2 is equal to 0
        {
            align <- hirschberg_template(seq1[1:x_mid], DNAString(""), align, match, mismatch, gap)# call hirschberg function for the first half of seq1 and for an empty DNAString object
        }
        else
        {
            align <- hirschberg_template(seq1[1:x_mid], seq2[1:y_mid], align, match, mismatch, gap)# call hirschberg function for the first half of seq1 and for the first part of seq2
        }
        
        # The second half
        if ((x_mid + 1) > x_len) # seq1 cannot be further divided
        {
            align <- hirschberg_template(DNAString(""), seq2[(y_mid+1), y_len], align, match, mismatch, gap)# call hirschberg function for an empty DNAString object and the second half of seq2
        }
        else if ((y_mid + 1) > y_len) # seq2 cannot be further divided
        {
            align <- hirschberg_template(seq1[(x_mid+1):x_len], DNAString(""), align, match, mismatch, gap)# call hirschberg function for the second half of seq1 and for an empty DNAString object
        }
        else 
        {
            align <- hirschberg_template(seq1[(x_mid+1):x_len], seq2[(y_mid+1), y_len], align, match, mismatch, gap)# call hirschberg function for the second half of seq1 and the second part of seq2
        }
    }
    return(align)
}

NWScore <- function(X,Y){
  gap = -2
  match = 2
  mismatch = -1
  S <- matrix(nrow=length(X)+1, ncol=length(Y)+1)
  S[1,1] <- 0
  
  for(j in (2:length(Y))){
    S[1,j] <- gap*(j-1)
  }
  
  for(i in (2:length(X))){
    S[i,1] <- gap*(i-1)
    for(j in (2:length(Y))){
      if(X[i]==Y[j]){
        scoreSub = S[i-1, j - 1] + match
      }
      else{
        scoreSub = S[i-1, j - 1] + mismatch
      }
      
      scoreDel = S[i-1, j] + gap
      scoreIns = S[i, j - 1] + gap
      S[1, j] = max(scoreSub, scoreDel, scoreIns)
    }
  }
  LastLine <- c(1:length(Y))
  for(j in (1:length(Y))){
    LastLine[j] = S[1, j]
  }
  LastLine
}
seq1 <- DNAString("AGTACGCA")
seq2 <- DNAString("TATGC")
match <- 2
mismatch <- -1
gap <- -2

hirschberg_template(seq1, seq2, c(), match, mismatch, gap)





X <- seq2
Y <- seq1
S <- matrix(nrow=length(X)+1, ncol=length(Y)+1)
S[1,1] <- 0

for(j in (2:(length(Y)+1))){
  S[1,j] <- gap*(j-1)
}
for(i in (2:(length(X)+1))){
  S[i,1] <- gap*(i-1)
}
for(i in (2:length(X))){
  
  for(j in (2:length(Y))){
    if(X[i]==Y[j]){
      scoreSub = S[i-1, j - 1] + match
    }
    else{
      scoreSub = S[i-1, j - 1] + mismatch
    }
    
    scoreDel = S[i-1, j] + gap
    scoreIns = S[i, j - 1] + gap
    S[i, j] = max(scoreSub, scoreDel, scoreIns)
  }
}
LastLine <- c(1:length(Y))
for(j in (1:length(Y))){
  LastLine[j] = S[1, j]
}
  