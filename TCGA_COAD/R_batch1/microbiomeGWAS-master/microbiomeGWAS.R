funED <- function(distMat, nPerm = 1e6, nPermEach = 1e4, soFile1 = "dExp1.so", soFile2 = "dExp2.so"){
  dyn.load(soFile1)
  dyn.load(soFile2)
  nSample <- nrow(distMat)
  out <- .C("dExp1", nSample = as.integer(nSample), distVec = as.vector(distMat), nD = 23L, eD = rep(0.1, 23L))
  eD1 <- out$eD
  nPermGroup <- ceiling(nPerm / nPermEach)
  eD2 <- rep(0, 23)
  for(i in 1:nPermGroup){
    idxMat <- matrix(NA, nSample, nPermEach)
    for(s in 1:nPermEach)idxMat[, s] <- sample(nSample)
    out <- .C("dExp2", nSample = as.integer(nSample), distVec = as.vector(distMat), nPerm = as.integer(nPermEach), idxMat = as.vector(idxMat - 1L), nD = 23L, eD = rep(0.1, 23L))
    eD2 <- eD2 + out$eD
  }
  eD2 <- eD2 / nPermGroup
  eD <- eD2;
  eD[c(1:6, 9:12)] <- eD1[c(1:6, 9:12)]
  names(eD) <- c("eDij", "eDij2", "eDijDik", "eDij3", "eDij2Dik", "eDijDjkDik", "eDijDjkDkl", "eDijDikDil", "eDij4", "eDij3Dik", "eDij2Dik2", "eDij2DjkDik", "eDij2DjkDkl", "eDijDjk2Dkl", "eDij2DikDil", "eDijDjkDikDil", "eDijDjkDklDil", "eDij2Dkl2", "eDijDjkDklDlm", "eDijDikDilDim", "eDijDikDilDlm", "eDijDikDlm2", "eDijDikDlmDln")
  return(eD)
}

funEG <- function(f){
  if(any(is.na(f)))stop("No NA is allowed in MAF!")
  if(is.vector(f)){
    tempLevel <- unique(f)
    p0 <- (1 - tempLevel)^2
    p1 <- 2 * tempLevel * (1 - tempLevel)
    tempIdx <- as.integer(factor(f, levels = tempLevel))
  }else if(is.matrix(f) && ncol(f) == 2){
    tempPaste <- paste(f[, 1], f[, 2])
    tempIdx <- which(!duplicated(tempPaste))
    tempLevel <- tempPaste[tempIdx]
    p0 <- f[tempIdx, 1]
    p1 <- f[tempIdx, 2]
    tempIdx <- as.integer(factor(tempPaste, levels = tempLevel))
  }else stop("f must be a vector or a matrix with 2 columns!")
  p2 <- 1 - (p0 + p1)
  pMat <- cbind(p0, p1, p2)
  pMat[pMat < 0 & pMat > -1e-10] <- 0
  pMat[pMat > 1 & pMat - 1 < 1e-10] <- 1
  if(any(pMat < 0) || any(pMat > 1))stop("Invalid value of f")
  eGij <- eGij2 <- eGijGik <- eGij3 <- eGij2Gik <- eGijGjkGik <- eGijGjkGkl <- eGijGikGil <- eGij4 <- eGij3Gik <- eGij2Gik2 <- eGij2GjkGik <- eGij2GjkGkl <- eGijGjk2Gkl <- eGij2GikGil <- eGijGjkGikGil <- eGijGjkGklGil <- eGijGjkGklGlm <- eGijGikGilGim <- eGijGikGilGlm <- 0
  for(i in 0:2){
    for(j in 0:2){
      Pij <- pMat[, i + 1] * pMat[, j + 1]
      eGij <- eGij + Pij * abs(i - j)
      eGij2 <- eGij2 + Pij * abs(i - j)^2
      eGij3 <- eGij3 + Pij * abs(i - j)^3
      eGij4 <- eGij4 + Pij * abs(i - j)^4
      for(k in 0:2){
        Pijk <- Pij * pMat[, k + 1]
        eGijGik <- eGijGik + Pijk * abs(i - j) * abs(i - k)
        eGij2Gik <- eGij2Gik + Pijk * abs(i - j)^2 * abs(i - k)
        eGijGjkGik <- eGijGjkGik + Pijk * abs(i - j) * abs(j - k) * abs(i - k)
        eGij3Gik <- eGij3Gik + Pijk * abs(i - j)^3 * abs(i - k)
        eGij2Gik2 <- eGij2Gik2 + Pijk * abs(i - j)^2 * abs(i - k)^2
        eGij2GjkGik <- eGij2GjkGik + Pijk * abs(i - j)^2 * abs(j - k) * abs(i - k)
        for(l in 0:2){
          Pijkl <- Pijk * pMat[, l + 1]
          eGijGjkGkl <- eGijGjkGkl + Pijkl * abs(i - j) * abs(j - k) * abs(k - l)
          eGijGikGil <- eGijGikGil + Pijkl * abs(i - j) * abs(i - k) * abs(i - l)
          eGij2GjkGkl <- eGij2GjkGkl + Pijkl * abs(i - j)^2 * abs(j - k) * abs(k - l)
          eGijGjk2Gkl <- eGijGjk2Gkl + Pijkl * abs(i - j) * abs(j - k)^2 * abs(k - l)
          eGij2GikGil <- eGij2GikGil + Pijkl * abs(i - j)^2 * abs(i - k) * abs(i - l)
          eGijGjkGikGil <- eGijGjkGikGil + Pijkl * abs(i - j) * abs(j - k) * abs(i - k) * abs(i - l)
          eGijGjkGklGil <- eGijGjkGklGil + Pijkl * abs(i - j) * abs(j - k) * abs(k - l) * abs(i - l)
          for(m in 0:2){
            Pijklm <- Pijkl * pMat[, m + 1]
            eGijGjkGklGlm <- eGijGjkGklGlm + Pijklm * abs(i - j) * abs(j - k) * abs(k - l) * abs(l - m)
            eGijGikGilGim <- eGijGikGilGim + Pijklm * abs(i - j) * abs(i - k) * abs(i - l) * abs(i - m)
            eGijGikGilGlm <- eGijGikGilGlm + Pijklm * abs(i - j) * abs(i - k) * abs(i - l) * abs(l - m)
          }
        }
      }
    }
  }
  eG <- cbind(eGij = eGij, 
              eGij2 = eGij2, 
              eGijGik = eGijGik, 
              eGij3 = eGij3, 
              eGij2Gik = eGij2Gik, 
              eGijGjkGik = eGijGjkGik, 
              eGijGjkGkl = eGijGjkGkl, 
              eGijGikGil = eGijGikGil, 
              eGij4 = eGij4, 
              eGij3Gik = eGij3Gik, 
              eGij2Gik2 = eGij2Gik2,
              eGij2GjkGik = eGij2GjkGik, 
              eGij2GjkGkl = eGij2GjkGkl, 
              eGijGjk2Gkl = eGijGjk2Gkl, 
              eGij2GikGil = eGij2GikGil, 
              eGijGjkGikGil = eGijGjkGikGil, 
              eGijGjkGklGil = eGijGjkGklGil, 
              eGij2Gkl2 = eGij2^2, 
              eGijGjkGklGlm = eGijGjkGklGlm, 
              eGijGikGilGim = eGijGikGilGim, 
              eGijGikGilGlm = eGijGikGilGlm, 
              eGijGikGlm2 = eGijGik * eGij2, 
              eGijGikGlmGln = eGijGik^2)
  eG <- eG[tempIdx, ]
  if(nrow(eG) == 1)eG <- eG[1, ]
  return(eG)
}

funEstimate <- function(N, eD, eG){
  eDij <- eD["eDij"]
  eDij2 <- eD["eDij2"]
  eDijDik <- eD["eDijDik"]
  eDij3 <- eD["eDij3"]
  eDij2Dik <- eD["eDij2Dik"]
  eDijDjkDik <- eD["eDijDjkDik"]
  eDijDjkDkl <- eD["eDijDjkDkl"]
  eDijDikDil <- eD["eDijDikDil"]
  eDij4 <- eD["eDij4"]
  eDij3Dik <- eD["eDij3Dik"]
  eDij2Dik2 <- eD["eDij2Dik2"]
  eDij2DjkDik <- eD["eDij2DjkDik"]
  eDij2DjkDkl <- eD["eDij2DjkDkl"]
  eDijDjk2Dkl <- eD["eDijDjk2Dkl"]
  eDij2DikDil <- eD["eDij2DikDil"]
  eDijDjkDikDil <- eD["eDijDjkDikDil"]
  eDijDjkDklDil <- eD["eDijDjkDklDil"]
  eDij2Dkl2 <- eD["eDij2Dkl2"]
  eDijDjkDklDlm <- eD["eDijDjkDklDlm"]
  eDijDikDilDim <- eD["eDijDikDilDim"]
  eDijDikDilDlm <- eD["eDijDikDilDlm"]
  eDijDikDlm2 <- eD["eDijDikDlm2"]
  eDijDikDlmDln <- eD["eDijDikDlmDln"]
  if(is.matrix(eG)){
    eGij <- eG[, "eGij"]
    eGij2 <- eG[, "eGij2"]
    eGijGik <- eG[, "eGijGik"]
    eGij3 <- eG[, "eGij3"]
    eGij2Gik <- eG[, "eGij2Gik"]
    eGijGjkGik <- eG[, "eGijGjkGik"]
    eGijGjkGkl <- eG[, "eGijGjkGkl"]
    eGijGikGil <- eG[, "eGijGikGil"]
    eGij4 <- eG[, "eGij4"]
    eGij3Gik <- eG[, "eGij3Gik"]
    eGij2Gik2 <- eG[, "eGij2Gik2"]
    eGij2GjkGik <- eG[, "eGij2GjkGik"]
    eGij2GjkGkl <- eG[, "eGij2GjkGkl"]
    eGijGjk2Gkl <- eG[, "eGijGjk2Gkl"]
    eGij2GikGil <- eG[, "eGij2GikGil"]
    eGijGjkGikGil <- eG[, "eGijGjkGikGil"]
    eGijGjkGklGil <- eG[, "eGijGjkGklGil"]
    eGijGjkGklGlm <- eG[, "eGijGjkGklGlm"]
    eGijGikGilGim <- eG[, "eGijGikGilGim"]
    eGijGikGilGlm <- eG[, "eGijGikGilGlm"]
  }else if(is.vector(eG)){
    eGij <- eG["eGij"]
    eGij2 <- eG["eGij2"]
    eGijGik <- eG["eGijGik"]
    eGij3 <- eG["eGij3"]
    eGij2Gik <- eG["eGij2Gik"]
    eGijGjkGik <- eG["eGijGjkGik"]
    eGijGjkGkl <- eG["eGijGjkGkl"]
    eGijGikGil <- eG["eGijGikGil"]
    eGij4 <- eG["eGij4"]
    eGij3Gik <- eG["eGij3Gik"]
    eGij2Gik2 <- eG["eGij2Gik2"]
    eGij2GjkGik <- eG["eGij2GjkGik"]
    eGij2GjkGkl <- eG["eGij2GjkGkl"]
    eGijGjk2Gkl <- eG["eGijGjk2Gkl"]
    eGij2GikGil <- eG["eGij2GikGil"]
    eGijGjkGikGil <- eG["eGijGjkGikGil"]
    eGijGjkGklGil <- eG["eGijGjkGklGil"]
    eGijGjkGklGlm <- eG["eGijGjkGklGlm"]
    eGijGikGilGim <- eG["eGijGikGilGim"]
    eGijGikGilGlm <- eG["eGijGikGilGlm"]
  }else stop("Something wrong!")
  varSM1 <- N * (N - 1) / 2 * (eGij2 - eGij^2) * eDij2
  varSM2 <- N * (N - 1) * (N - 2) * (eGijGik - eGij^2) * eDijDik
  varSM <- varSM1 + varSM2
  eSM31 <- N * (N - 1) / 2 * (eGij3 - 3 * eGij2 * eGij + 2 * eGij^3) * eDij3
  eSM32 <- 3 * N * (N - 1) * (N - 2) * (eGij2Gik - eGij2 * eGij - 2 * eGijGik * eGij + 2 * eGij^3) * eDij2Dik
  eSM33 <- N * (N - 1) * (N - 2) * (eGijGjkGik - 3 * eGijGik * eGij + 2 * eGij^3) * eDijDjkDik
  eSM34 <- 3 * N * (N - 1) * (N - 2) * (N - 3) * (eGijGjkGkl - 2 * eGijGik * eGij + eGij^3) * eDijDjkDkl
  eSM35 <- N * (N - 1) * (N - 2) * (N - 3) * (eGijGikGil - 3 * eGijGik * eGij + 2 * eGij^3) * eDijDikDil
  eSM3 <- eSM31 + eSM32 + eSM33 + eSM34 + eSM35
  eSM41 <- N * (N - 1) / 2 * (eGij4 - 4 * eGij3 * eGij + 6 * eGij2 * eGij^2 - 3 * eGij^4) * eDij4
  eSM42 <- 4 * N * (N - 1) * (N - 2) * (eGij3Gik - (3 * eGij2Gik + eGij3) * eGij + 3 * (eGijGik + eGij2) * eGij^2 - 3 * eGij^4) * eDij3Dik
  eSM43 <- 3 * N * (N - 1) * (N - 2) * (eGij2Gik2 - 4 * eGij2Gik * eGij + (4 * eGijGik + 2 * eGij2) * eGij^2 - 3 * eGij^4) * eDij2Dik2
  eSM44 <- 6 * N * (N - 1) * (N - 2) * (eGij2GjkGik - 2 * (eGijGjkGik + eGij2Gik) * eGij + (5 * eGijGik + eGij2) * eGij^2 - 3 * eGij^4) * eDij2DjkDik
  eSM45 <- 12 * N * (N - 1) * (N - 2) * (N - 3) * (eGij2GjkGkl - (2 * eGijGjkGkl + eGij2Gik) * eGij + 3 * eGijGik * eGij^2 - eGij^4) * eDij2DjkDkl
  eSM46 <- 6 * N * (N - 1) * (N - 2) * (N - 3) * (eGijGjk2Gkl - 2 * (eGijGjkGkl + eGij2Gik) * eGij + (4 * eGijGik + eGij2) * eGij^2 - 2 * eGij^4) * eDijDjk2Dkl
  eSM47 <- 6 * N * (N - 1) * (N - 2) * (N - 3) * (eGij2GikGil - 2 * (eGijGikGil + eGij2Gik) * eGij + (5 * eGijGik + eGij2) * eGij^2 - 3 * eGij^4) * eDij2DikDil
  eSM48 <- 12 * N * (N - 1) * (N - 2) * (N - 3) * (eGijGjkGikGil - (eGijGikGil + eGijGjkGik + 2 * eGijGjkGkl) * eGij + 5 * eGijGik * eGij^2 - 2 * eGij^4) * eDijDjkDikDil
  eSM49 <- 3 * N * (N - 1) * (N - 2) * (N - 3) * (eGijGjkGklGil - 4 * eGijGjkGkl * eGij + 4 * eGijGik * eGij^2 - eGij^4) * eDijDjkDklDil
  eSM410 <- 3 / 4 * N * (N - 1) * (N - 2) * (N - 3) * (eGij2 - eGij^2)^2 * eDij2Dkl2
  eSM411 <- 12 * N * (N - 1) * (N - 2) * (N - 3) * (N - 4) * (eGijGjkGklGlm - 2 * eGijGjkGkl * eGij + eGijGik * eGij^2) * eDijDjkDklDlm
  eSM412 <- N * (N - 1) * (N - 2) * (N - 3) * (N - 4) * (eGijGikGilGim - 4 * eGijGikGil * eGij + 6 * eGijGik * eGij^2 - 3 * eGij^4) * eDijDikDilDim
  eSM413 <- 12 * N * (N - 1) * (N - 2) * (N - 3) * (N - 4) * (eGijGikGilGlm - (2 * eGijGjkGkl + eGijGikGil) * eGij + 3 * eGijGik * eGij^2 - eGij^4) * eDijDikDilDlm
  eSM414 <- 3 * N * (N - 1) * (N - 2) * (N - 3) * (N - 4) * (eGijGik - eGij^2) * (eGij2 - eGij^2) * eDijDikDlm2
  eSM415 <- 3 * N * (N - 1) * (N - 2) * (N - 3) * (N - 4) * (N - 5) * (eGijGik - eGij^2)^2 * eDijDikDlmDln
  eSM4 <- eSM41 + eSM42 + eSM43 + eSM44 + eSM45 + eSM46 + eSM47 + eSM48 + eSM49 + eSM410 + eSM411 + eSM412 + eSM413 + eSM414 + eSM415
  estimate <- cbind(varSM = varSM, eSM3 = eSM3, skewSM = eSM3 / (varSM)^1.5, eSM4 = eSM4, kurtSM = eSM4 / varSM^2 - 3)
  if(nrow(estimate) == 1)estimate <- estimate[1, ]
  return(estimate)
}

###########################################################################

#funAdjust

funAdjust <- function(Z, skew, kurt){
  xi <- (sqrt(1 + 2 * skew * Z) - 1) / skew
  sigma <- sqrt(1 + skew * xi)
  pSkew <- exp(xi^2 / 2 + skew / 6 * xi^3) * pnorm(sigma * xi, lower.tail = FALSE) / exp(xi * (Z - sigma^2 * xi / 2))
  alpha <- kurt / 6
  beta <- skew / 2
  eta <- beta / 6 / alpha^2 - beta^3 / 27 / alpha^3 + Z / 2 / alpha
  delta <- eta^2 + (1 / 3 / alpha - beta^2 / 9 / alpha^2)^3
  deltaPositive <- pmax(0, delta)
  deltaNegative <- pmin(0, delta)
  temp1 <- eta + sqrt(deltaPositive)
  temp2 <- eta - sqrt(deltaPositive)
  xi <- -beta / 3 / alpha + sign(temp1) * (sign(temp1) * temp1)^(1 / 3) + sign(temp2) * (sign(temp2) * temp2)^(1 / 3)
  tempIdx <- which(delta < 0)
  if(length(tempIdx) > 0){
    complexDelta <- complex(real = deltaNegative)
    temp1 <- (eta + sqrt(complexDelta))^(1 / 3)
    temp2 <- (eta - sqrt(complexDelta))^(1 / 3)
    xi1 <- Re(-beta / 3 / alpha + temp1 + temp2)
    xi2 <- Re(-beta / 3 / alpha + temp1 * complex(real = -1 / 2, imaginary = sqrt(3) / 2) + temp2 * complex(real = -1 / 2, imaginary = -sqrt(3) / 2))
    xi3 <- Re(-beta / 3 / alpha + temp1 * complex(real = -1 / 2, imaginary = -sqrt(3) / 2) + temp2 * complex(real = -1 / 2, imaginary = sqrt(3) / 2))
    xiMat <- cbind(xi1, xi2, xi3)[tempIdx, ]
    xiMat[xiMat > Z[tempIdx] | xiMat < 0] <- 0
    xi[tempIdx] <- apply(xiMat, 1, max)
  }
  sigma <- sqrt(1 + skew * xi + kurt / 3 * xi^2)
  pKurt <- exp(xi^2 / 2 + skew / 6 * xi^3 + kurt / 24 * xi^4) * pnorm(sigma * xi, lower.tail = FALSE) / exp(xi * (Z - sigma^2 * xi / 2))
  return(cbind(p_adj_skew = pSkew, p_adj_skew_and_kurtosis = pKurt))
}


sm <- function(snp1, distmat){
  print(snp1[1:3])
  Gij = abs(outer(snp1, snp1, '-'))
  dij = distmat - mean(c(distmat[lower.tri(distmat)]))
  return(sum(Gij[lower.tri(distmat)] * dij[lower.tri(distmat)]))
}

distRes = function(distMat, dataCovariate){
  tempList <- lapply(dataCovariate, function(x)outer(x, x, function(a, b)abs(a - b)))
  tempList <- lapply(tempList, function(x)x[lower.tri(x)])
  dataModel <- data.frame(y = distMat[lower.tri(distMat)], tempList)
  tempFormula <- as.formula(paste0("y ~ 1 + ", paste(colnames(dataCovariate), collapse = " + ")))
  tempModel <- lm(tempFormula, data = dataModel)
  distMatRes <- distMat; distMatRes[, ] <- 0
  distMatRes[lower.tri(distMatRes)] <- tempModel$residuals
  distMatRes <- distMatRes + t(distMatRes)
  return(distMatRes)
}

funMain <- function(SM, distMat, eD, eG){
  estimate <- funEstimate(nrow(distMat), eD, eG)
  varSM <- estimate[, "varSM"]
  skewSM <- estimate[, "skewSM"]
  kurtSM <- estimate[, "kurtSM"]
  ZM <- SM / sqrt(varSM)
  pM <- pMS <- pMK <- pnorm(ZM, lower.tail = FALSE)
  tempIdx <- which(!is.na(ZM) & ZM > 0)
  pMSK <- funAdjust(ZM[tempIdx], skewSM[tempIdx], kurtSM[tempIdx])
  pMS[tempIdx] <- pMSK[, 1]
  pMK[tempIdx] <- pMSK[, 2]
  return(cbind(ZM, pM , pMS, pMK ))
}

# test results 
# 
# packageDir = ('/fh/fast/sun_w/licai/R_batch5/microbiomeGWAS-master/')
# setwd(packageDir)
# system(paste0('cd ', packageDir, '; sh compile.src.sh'))
# 
# soFile1 <- paste0(packageDir, "/lib/dExp1.so")
# soFile2 <- paste0(packageDir, "/lib/dExp2.so")
# soFile3 <- paste0(packageDir, "/lib/parsePlink.so")
# soFile4 <- paste0(packageDir, "/lib/parsePlink2.so")
# 
# snp1 = system("cut -f7 myfavpedfile.ped", intern = T)
# snp1
# snp1[which(snp1 == 'G G')] = 2
# snp1[which(snp1 == 'G A')] = 1
# snp1[which(snp1 == 'A A')] = 0
# snp1 = as.numeric(snp1)
# 
# snp2 = system("cut -f11 myfavpedfile.ped", intern = T)
# snp2
# snp2[which(snp2 == 'G G')] = 2
# snp2[which(snp2 == 'G A')] = 1
# snp2[which(snp2 == 'A A')] = 0
# snp2 = as.numeric(snp2)
# 
# snps = rbind(snp1, snp2)
# colnames(snps) = paste0("S", 1:379)
# 
# distMat <- as.matrix(read.table("data/distMat379.txt"))
# eD <- funED(distMat, soFile1 = soFile1, soFile2 = soFile2)
# eD
# 
# dataCovariate <- read.table('data/dataCovariate379.txt', header = TRUE, check.names = FALSE)
# 
# tempList <- lapply(dataCovariate, function(x)outer(x, x, function(a, b)abs(a - b)))
# tempList <- lapply(tempList, function(x)x[lower.tri(x)])
# dataModel <- data.frame(y = distMat[lower.tri(distMat)], tempList)
# tempFormula <- as.formula(paste0("y ~ 1 + ", paste(colnames(dataCovariate), collapse = " + ")))
# tempModel <- lm(tempFormula, data = dataModel)
# distMatRes <- distMat; distMatRes[, ] <- 0
# distMatRes[lower.tri(distMatRes)] <- tempModel$residuals
# distMatRes <- distMatRes + t(distMatRes)
# 
# SM = apply(snps, 1,sm, distMatRes)
# mafs = apply(snps, 1, mean)/2
# eG = funEG(mafs)
# estimate <- funEstimate(nrow(distMat), eD, eG)
# 
# varSM <- estimate[, "varSM"]
# skewSM <- estimate[, "skewSM"]
# kurtSM <- estimate[, "kurtSM"]
# ZM <- SM / sqrt(varSM)
# pM <- pMS <- pMK <- pnorm(ZM, lower.tail = FALSE)
# tempIdx <- which(!is.na(ZM) & ZM > 0)
# pMSK <- funAdjust(ZM[tempIdx], skewSM[tempIdx], kurtSM[tempIdx])
# pMS[tempIdx] <- pMSK[, 1]
# pMK[tempIdx] <- pMSK[, 2]
# 


