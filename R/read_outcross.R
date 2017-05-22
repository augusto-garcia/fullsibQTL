## FUNCTION: MTMIM_read.outcross
## Reads in markers and traits data
## Luciano s function changed for me to receive 2 or 3 elements in the first line and phenotype with comma (,) to separate phenotypes, as done with genotypes.

read.outcross.pheno <- function (dir, file) {
  if (missing(file))
    stop("missing file")
  if (!missing(dir) && dir != "")
    file <- file.path(dir, file)
  n.lines <- length(scan(file, what = character(), skip = 0, 
                         nlines = 0, blank.lines.skip = FALSE, quiet = TRUE, sep = "\n"))
  cur.mar <- 0
  cur.pheno <- 0
  n.phen <- 0
  flag <- 0
  for (i in 1:n.lines) {
    a <- scan(file, what = character(), skip = i - 1, nlines = 1,
              blank.lines.skip = TRUE, quiet = TRUE)
    if (length(a) == 0)
      next
    if (length(grep("#", a[1])) != 0)
    next
    if (flag == 0) { #reading first line
      flag <- 1
      #if (length(a) != 3) #Rod
      #  stop("The first line of the input file must have the following information: 'number of individuals', 'number of markers', and 'number of traits'. These numbers must be separated with an empty space. For instance, 10 5 0.", call.= TRUE)
      n.ind <- as.numeric(a[1])
      n.mar <- as.numeric(a[2])
      if(length(a) == 3)#Rod
        n.phen <- as.numeric(a[3]) #
      else
        n.phen <- 0          
      cat(" Working...\n\n")
      marnames <- rep("", n.mar)
      geno <- matrix(0, ncol = n.mar, nrow = n.ind)
      segr.type <- character(n.mar)
      if (n.phen == 0) {
        pheno <- numeric(0) #matrix(1:n.ind, ncol = 1)
        phenonames <- character(0) #c("number")
      }
      else {
        pheno <- matrix(0, ncol = n.phen, nrow = n.ind)
        phenonames <- rep("", n.phen)
      }
    } #finishes reading first line in the file (flag==0)
    else {#reading lines of markers and traits
      if (substring(a[1], 1, 1) == "*") {#reading lines of markers and traits that start with "*"
        cur.mar <- cur.mar + 1
        cur.row <- 1
        if (cur.mar > n.mar) {#reading lines of traits that start with "*"
          cur.pheno <- cur.pheno + 1
          if (cur.pheno > n.phen) 
            next
          phenonames[cur.pheno] <- substring(a[1], 2)
          if (length(a) > 1) {
            p <- a[-1]
            p[p == "-"] <- NA
            p <- paste(p, collapse = "")#Rod
            p <- unlist(strsplit(p, ","))#Rod
            n <- length(p)
            oldna <- is.na(p)
            numerp <- suppressWarnings(as.numeric(p))
            newna <- is.na(numerp)
            wh <- !oldna & newna
            if (any(wh)) {
              droppedasmissing <- unique(p[wh])
              if (length(droppedasmissing) > 1) {
                themessage <- paste("The values", paste("\"", 
                                                        droppedasmissing, "\"", sep = "", collapse = " "))
                themessage <- paste(themessage, " for pheno \"", 
                                    phenonames[cur.pheno], "\" were", sep = "")
              }
              else {
                themessage <- paste("The value \"", droppedasmissing, 
                                    "\" ", sep = "")
                themessage <- paste(themessage, " for pheno \"", 
                                    phenonames[cur.pheno], "\" was", sep = "")
              }
              themessage <- paste(themessage, "interpreted as missing.")
              warning(themessage, call.= FALSE)
            }
            pheno[cur.row + (0:(n - 1)), cur.pheno] <- numerp
          }
          else n <- 0
          cur.row <- cur.row + n
        }#finishes reading lines of traits that start with "*" (cur.mar > n.mar)
        else {#reading lines of markers that start with "*" (cur.mar <= n.mar)
          marnames[cur.mar] <- substring(a[1], 2)
          if (length(a) < 2) {
            stop("the segregation type of marker ", marnames[cur.mar],
                 " should be placed next to its name (on the same line)")
          }
          segr.type[cur.mar] <- a[2]
          if (length(a) > 2) {
            g <- paste(a[c(-1, -2)], collapse = "")
            g <- unlist(strsplit(g, ","))
            n <- length(g)
            geno[cur.row + (0:(n - 1)), cur.mar] <- as.character(g)
          }
          else n <- 0
          cur.row <- cur.row + n
        }#finishes reading lines of markers that start with "*" (cur.mar <= n.mar)
      }#finishes reading lines of markers and traits that start with "*"
      else {#reading lines of markers and traits that do not start with "*"
        if (cur.mar > n.mar) {
          a <- paste(a, collapse = "")#Rod
          a <- unlist(strsplit(a, ","))#Rod
          a[a == "-" | a == "NA"] <- NA

          n <- length(a)
          pheno[cur.row + (0:(n - 1)), cur.pheno] <- as.numeric(a)
          cur.row <- cur.row + n
        }#finishes reading lines of traits that do not start with "*"
        else {
          g <- paste(a, collapse = "")
          g <- unlist(strsplit(g, ","))
          n <- length(g)
          geno[cur.row + (0:(n - 1)), cur.mar] <- as.character(g)
          cur.row <- cur.row + n
        }#finishes reading lines of markers that do not start with "*"
      }#finishes reading lines of markers and traits that do not start with "*"
    }#finishes reading lines of markers and traits 
  } #close loop for lines (i)
  colnames(geno) <- marnames
  geno[!is.na(geno) & geno == "-"] <- NA
  temp.data <- codif.data(geno, segr.type)
  geno <- temp.data[[1]]
  segr.type.num <- temp.data[[2]]
  rm(temp.data)
  if(n.phen != 0) {
    colnames(pheno) <- phenonames
    pheno[!is.na(pheno) & pheno == "-"] <- NA
  }
  cat(" --Read the following data:\n")
  cat("\tNumber of individuals: ", n.ind, "\n")
  cat("\tNumber of markers:     ", n.mar, "\n")
  cat("\tNumber of traits:      ", n.phen, "\n")
  if(n.phen != 0) {
    miss.value.pheno <- apply((apply(pheno, 2,is.na)),2,sum)
    cat("\tMissing trait values:      ", "\n")
    for(i in 1:n.phen) {
      cat("\t",formatC(paste(colnames(pheno)[i],":",sep=""),width=max(nchar(paste(colnames(pheno),":",sep="")))), miss.value.pheno[i], "\n")
    }
  }

  structure(list(geno = geno, n.ind = n.ind, n.mar = n.mar,
                 segr.type = segr.type, segr.type.num = segr.type.num,
                 n.phen = n.phen, pheno = pheno, 
                 input = file),
            class = "outcross")
}

print.outcross <- function (x, ...) {
    if (!any(class(x) == "outcross"))
        stop("this is not an object of class 'outcross'")
    cat("  This is an object of class 'outcross'\n")
    cat("    No. individuals:   ", x$n.ind, "\n")
    cat("    No. markers:       ", x$n.mar, "\n")
    cat("    Segregation types:\n")
    quant <- cbind(table(x$segr.type))
    for (i in 1:length(quant)) {
        cat(paste("       ", rownames(quant)[i], ":\t", quant[i], 
                  "\n", sep = ""))
    }
    cat("    No. traits:        ", x$n.phen, "\n")
    if(x$n.phen > 0) {
        miss.value <- apply((apply(x$pheno, 2,is.na)),2,sum)
        cat("    Missing trait values:", "\n")
        for (i in 1:x$n.phen) {
                                        #cat(paste("       ", colnames(x$pheno)[i], "\n", sep = ""))
            cat("\t",formatC(paste(colnames(x$pheno)[i],":", sep=""), width= max(nchar(paste(colnames(x$pheno),":",sep="")))), miss.value[i], "\n")
        }
    }
}

