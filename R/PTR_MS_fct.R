##############
### Start ####
##############

# Library ####
# library(calibrate)
# library(carData)
# library(car)
# library(cluster)
# library(DescTools)
# library(dygraphs)
# library(ggplot2)
# library(hrbrthemes)
# library(ica)
# library(kableExtra)
# library(knitr)
# library(magrittr)
# library(MALDIquant)
# library(MASS)
# library(rhdf5)
# library(rnirs)
# library(scales)
# library(stringr)
# library(tidyverse)
# library(tinytex)
# library(viridis)
# library(zoo)

# Fonction d'importation de sp : ####

#' import_h5
#' import tous les fichiers .h5 présents dans le dossier h5 du wd
#' @param wd
#'
#' @return sp
#' @export
#'
#' @examples
import_h5 <- function(){

  if(("Figures" %in% dir())==FALSE){
    dir.create("Figures")
  }

  # Importation des datas et reduction ####
  f_h5 <- dir("h5") %>% grep(".h5",.)         # recherche des fichiers h5
  ls_h5 <- vector("list", length(f_h5))   # prepare une liste des acquisitions
  names(ls_h5) <- nm_ls(f_h5)             # avec les bons noms

  for (i in 1:length(ls_h5)){
    ls_h5[[i]] <- paste0("h5/",dir("h5")[f_h5[i]]) %>% H5Fopen() # import des fichiers
  }

  mt_h5 <- meta_h5(ls_h5)                 # prepare le fichier meta

  sp <- sp_red(ls_h5,mt_h5)              # reduction du spectre

  sp$h5 <- list(f_h5 = f_h5,
                ls_h5 = ls_h5,
                mt_h5 = mt_h5,
                wd = getwd())

  h5closeAll()

  print.h("Calcul des AUC en cours")

  length(citation_list) %>% sample(1) %>% citation_list[[.]] %>% cat()
  sp$MS.AUC <- peak_AUC(sp)

  print.h("Export des figures d'evalutation de la calibrations")
  # print_calib(42.0344, "01_Acetonitrile",     sp$xMS, sp$MS, addProton = FALSE)
  # print_calib(45.0340, "02_Acetaldehyde",     sp$xMS, sp$MS, addProton = FALSE)
  # print_calib(54.0344, "03_Acrylonitrile",    sp$xMS, sp$MS, addProton = FALSE)
  print_calib(58.0418, "04_Acetone",          sp$xMS, sp$MS)
  # print_calib(69.0704, "05_Isopropene",       sp$xMS, sp$MS, addProton = FALSE)
  # print_calib(73.0653, "06_Butanone",         sp$xMS, sp$MS, addProton = FALSE)
  # print_calib(79.0548, "07_Benzene",          sp$xMS, sp$MS, addProton = FALSE)
  # print_calib(93.0704, "08_Toluene",          sp$xMS, sp$MS, addProton = FALSE)
  # print_calib(107.0861,"09_O_xylene",         sp$xMS, sp$MS, addProton = FALSE)
  # print_calib(121.1017,"10_Trimethylbenzene", sp$xMS, sp$MS, addProton = FALSE)
  print_calib(136.1252,"11_Terpene",          sp$xMS, sp$MS)
  print_calib(204.1878,"12_Sesquiterpenes",   sp$xMS, sp$MS)
  # print_calib(191.2695,"DEET", sp$xMS, sp$MS)

  print.h("Sauvegarde de l'environnement (.RData)")
  str_split(wd,"/")[[1]] %>% length() %>% str_split(wd,"/")[[1]][.] %>% paste0(".RData") %>% save.image()

  return(sp)
}

# Fonctions d'importations et reductions : ####

nm_ls <- function(f_h5){
  nm_h5 <- str_split(dir("h5")[f_h5],"_20......_.......h5",  simplify = TRUE)[,1]
  if(length(nm_h5) != length(unique(nm_h5))){
    unm <- unique(nm_h5)
    for (i in 1:length(unm)){
      inm <- which(nm_h5 == unm[i])
      if(length(inm) > 1){
        eph <- log(length(inm),10) %>% floor() %>% add(1)
        nm_h5[inm] <- paste0("000", 1:length(inm)) %>% str_sub(-eph) %>% paste(nm_h5[inm], ., sep = "_")
      }
    }
  }
  return(nm_h5)
}
# retourne le nom des acquisitions en supprimant la date et le .h5 tout en numerotant
# les acquisitons eponymes.

meta_h5 <- function(ls_h5){
  mt_h5 <- matrix(NA, nrow = length(ls_h5), ncol = 3, dimnames = list(names(ls_h5), c("nbr_MS","start","end")))
  mt_h5[1,c(1,3)] <- dim(ls_h5[[1]]$PeakData$PeakData)[4]
  mt_h5[1, 2] <- 1
  if(length(ls_h5)>1){
    for(i in seq(2, length(ls_h5), by = 1)){
      mt_h5[i,1] <- dim(ls_h5[[i]]$PeakData$PeakData)[4]
      mt_h5[i,2] <- mt_h5[i-1,3] + 1
      mt_h5[i,3] <- mt_h5[i,2] + mt_h5[i,1] -1
    }
  }
  return(mt_h5)
}
# retourne un premier fichier meta a partir de la liste des acquisitions.
# Ce fichier comportement uniquement le nombre d'acquisitions et leur ordonnancement.

reduction2 <- function(ls_h5){
  # Import des datas
  mat_h5 <- red_ls(ls_h5)
        # import la matrices des intensites
  xMS_r <- colnames(mat_h5) %>% as.numeric()
        # defini l'abscisse
  seuil <- apply(mat_h5, 2, sd) %>% median() %>% multiply_by(20)
        # seuil = 20 fois l'ecart-type du bruit de fond

  # Suppressions des datas non informatives
  x_coef <- apply(mat_h5, 1, which.sup, threshold = seuil) %>% unlist() %>% unique() %>% sort()

    # ajoute des colonnes frontieres (pour l'alignement)
    x_zero <- c(1,x_coef[1]-4,x_coef[1]-3,x_coef[1]-2,x_coef[1]-1)
    if((x_coef[1] - x_coef[2]) != -1){
      x_zero <- c(x_zero, x_coef[1]+1,x_coef[1]+2,x_coef[1]+3,x_coef[1]+4)
    }

    for(i in 2:(length(x_coef)-1)){
      if( (x_coef[i] - x_coef[i-1]) != 1){
        x_zero <- c(x_zero, x_coef[i]-4,x_coef[i]-3,x_coef[i]-2,x_coef[i]-1)
      }
      if( (x_coef[i] - x_coef[i+1]) != -1){
        x_zero <- c(x_zero, x_coef[i]+4,x_coef[i]+3,x_coef[i]+2,x_coef[i]+1)
      }
    }
    x_zero <- c(x_zero, x_coef[i+1]+4,x_coef[i+1]+3,x_coef[i+1]+2,x_coef[i+1]+1)

    x_zero <- unique(x_zero) %>% sort()
    x_coef <- union(x_zero, x_coef) %>% unique() %>% sort()

  MS <- matrix(NA,nrow(mat_h5),length(x_coef),dimnames = list(rownames(mat_h5),xMS_r[x_coef]))
    # une matrice vide de bonne dimension
  eph <- mat_h5
  eph[,x_zero] <- runif(length(x_zero)*nrow(eph)) %>% round(2) %>% matrix(ncol=length(x_zero))
    # genere un bruit dans les colonnes a 0
  MS <- eph[,x_coef]
  xMS <- xMS_r[x_coef]
  return(list(MS,xMS,seuil))
}
# retourne une matrice et une abscisse allegees des zones sans information.

red_ls <- function(ls_h5){
  cit_rdm <- length(ls_h5) %>% divide_by(2) %>% ceiling() %>% sample(1)
    # pour passer le temps

  for(j in 1:length(ls_h5)){
    # extraction des abscisses
    xMS <- ls_h5[[j]]$FullSpectra$MassAxis
                    # extraction de l'abscisse [~160 000 pts]
    lMS_r <- coef_xMS_r(xMS, bn_x_m = 30, bn_x_M = 500, nbr_dec = 3)
                    # les coordonnes correspondantes aux abscisses
                    # les bornes sont limitees entre 30 et 500 m/z
                    # la precsion est au millieme
    xMS_r <- names(lMS_r) %>% as.numeric()
                    # les valeurs de l'abscisse

    # extraction des intensites
    all_MS <- t(ls_h5[[j]]$FullSpectra$TofData[,1,1,])
                    # extraction des intensites [ nbr d'acquisition x 160 000 pts]
    MS_r <- matrix(NA,nrow(all_MS),length(xMS_r),dimnames = list(1:nrow(all_MS),xMS_r))
    for (i in 1:length(lMS_r)){
      if(length(lMS_r[[i]])>1){ # si plusieurs coordonnees pour une meme abscisse
        MS_r[,i] <- rowMeans(all_MS[,lMS_r[[i]]]) # moyenne des intensites
      }else{                    # sinon
        MS_r[,i] <- all_MS[,lMS_r[[i]]] # intensite exacte
      }
    }
    row.names(MS_r) <- paste(names(ls_h5)[[j]],row.names(MS_r), sep = "_")

    # concatenation de la matrice intensite
    if(j == 1){
      mat_h5 <- MS_r
      # la premiere matrice
    }else if(ncol(MS_r) == ncol(mat_h5)){
      mat_h5 <- rbind(mat_h5, MS_r)
      # si tout se passe bien
    }else if(ncol(MS_r) > ncol(mat_h5)){
      mat_h5 <- rbind(mat_h5, MS_r[,1:ncol(mat_h5)])
      eph <- ncol(MS_r) - ncol(mat_h5)
      eph1 <- paste(str_flatten(head(colnames(mat_h5), n = 3), " "), "...", str_flatten(tail(colnames(mat_h5), n = eph), " "), "length", ncol(mat_h5))
      eph2 <- paste(str_flatten(head(xMS_r, n = 3), " "), "...", str_flatten(tail(xMS_r, n = eph), " "), "length", ncol(MS_r))
      dl <- paste("WARNING !", "xMS general :", eph1, paste("xMS",names(ls_h5)[j],":"), eph2, sep = "\n")
      cat(dl)
      # pour les matrice trop grandes
    }else if(ncol(MS_r) < ncol(mat_h5)){
      mat_h5 <- rbind(mat_h5[,1:ncol(MS_r)], MS_r)
      eph <- ncol(mat_h5) - ncol(MS_r)
      eph1 <- paste(str_flatten(head(colnames(mat_h5), n = 3), " "), "...", str_flatten(tail(colnames(mat_h5), n = eph), " "), "length", ncol(mat_h5))
      eph2 <- paste(str_flatten(head(xMS_r, n = 3), " "), "...", str_flatten(tail(xMS_r, n = eph), " "), "length", ncol(MS_r))
      dl <- paste("WARNING !", "xMS general :", eph1, paste("xMS",names(ls_h5)[j],":"), eph2, sep = "\n")
      cat(dl)
      # pour les matrices trop petites
    }

    print.h(paste0(names(ls_h5)[[j]]," # ",j,"/",length(ls_h5)))
      # donne une indication sur l'etat d'avancement

    if(j == cit_rdm){ # s'il est temps d'etre cite
      length(citation_list) %>% sample(1) %>% citation_list[[.]] %>% print()
    }
  }
  return(mat_h5)
}
# retourne une matrice composee des spectres de masses reduits (xMS arrondi au centieme).

coef_xMS_r <- function(xMS, bn_x_m = 30, bn_x_M = 500, nbr_dec = 3){
  xMS_round <- round(xMS, nbr_dec)

  bn_r_m <- det_c(bn_x_m,xMS_round)
  bn_r_M <- det_c((bn_x_M+0.01),xMS_round)
  xMS_r <- unique(xMS_round[bn_r_m:bn_r_M])
  lMS_r <- vector("list", length(xMS_r))
  names(lMS_r) <- xMS_r

  eph <- bn_r_m
  for(i in 1:length(lMS_r)){
    lMS_r[[i]] <- which(xMS_r[i] == xMS_round[(eph):(eph+10)]) + eph - 1
    eph <- eph + length(lMS_r[[i]])
  }
  return(lMS_r)
}
# retourne l'index de chaque dixieme de masse.

sp_red <- function(ls_h5,mt_h5){
  sp <- reduction2(ls_h5)

  names(sp) <- c("MS", "xMS", "threshold")
  sp <- c(sp,list("indX" = bucketing_xMS(sp$xMS, 31, 500)))

  print.h("Alignement")
  MS_unalign <- sp$MS
  xMS_unalign <- sp$xMS
  indX_n_al <- bucketing_xMS(xMS_unalign, 31, 500)
  seuil <- sp$threshold
  ls_align <- alignement_sp(sp)
  MS_align <- ls_align$MS_align
  peaks <- ls_align$peaks
  xMS <- colnames(MS_align) %>% as.numeric()
  indX <- bucketing_xMS(xMS, 31, 500)

  rm(sp)

  sp <- list("MS" = MS_align,
             "xMS" = xMS,
             "indX" = indX,
             "MS_n_al" = MS_unalign,
             "xMS_n_al" = xMS_unalign,
             "indX_n_al" = indX_n_al,
             "peaks" = peaks,
             "threshold" = seuil)

  print.h("Mise en forme")
  sp <- c(sp,
          list("MS_m" = t(apply(mt_h5,1,mean_mat, MS = sp$MS)),
               "MS_max" = t(apply(mt_h5,1,max_mat, MS = sp$MS)),
               "xT" = extract_timing(ls_h5, sp$MS),
               "Tacq" = unlist(lapply(ls_h5,acq_time)),
               "Dacq" = unlist(lapply(ls_h5,acq_day)),
               "xT_init" = extract_timing(ls_h5, sp$MS),
               "Tacq_init" = unlist(lapply(ls_h5,acq_time)),
               "Dacq_init" = unlist(lapply(ls_h5,acq_day))))

  print.h("Detection des pics")
  sp <- c(sp,
          list("pk" = sp_peaks(sp)))
  return(sp)
}
# retourne la liste sp comportant les spectres et l'axe des absisces

# Fonctions d'alignement : ####
align_print_ctr <- function(mat, sp, suffixe = "raw"){
  xMS <- colnames(mat) %>% as.numeric()
  indX <- bucketing_xMS(xMS, 31, 500)
  sup_threshold <- sapply(indX, sp_max, data = mat)
  sp_indX_th <- which(sup_threshold > 1000) %>% as.numeric() %>% indX[.]

  for(f in seq(1,length(sp_indX_th),by = 9)){
    title.tiff <- names(sp_indX_th)[c(f,f+8)]
    if(length(which.na(title.tiff)) >0){
      title.tiff[which.na(title.tiff)] <- names(sp_indX_th)[length(sp_indX_th)]
    }
    title.tiff <- str_flatten(title.tiff,"to") %>%
      paste0("Figures/Alignement/mass_at",.,"_",suffixe,".tiff")
    tiff(title.tiff, width = 800, height = 800,units = "px")
    layout(matrix(1:9, 3, 3, byrow = TRUE))
    par(mar = c(3,2,.5,.5), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0), xpd = FALSE, yaxt = "n")
    for(i in f:(f+8)){
      if(i <= length(sp_indX_th)){
        p.max <- apply(mat[,sp_indX_th[[i]]],1,which.max)
        y.max.all <- apply(mat[,sp_indX_th[[i]]],1,max)
        p.shift <- median(p.max) %>% subtract(p.max,.)
        x.max <- xMS[sp_indX_th[[i]]][p.max]
        y.max <- max(y.max.all)
        y.max.rdm <- length(x.max) %>% runif(y.max*1.1, y.max*1.3)

        matplot(xMS[sp_indX_th[[i]]], t(mat[,sp_indX_th[[i]]]), type = "l", ylim = c(0,y.max*1.3),
                ylab = " ", xlab = " ", col = alpha(viridis(length(p.max)),0.3), lty = 1)
        abline(v = median(x.max), col = "midnightblue", lwd = 2)
        points(x.max, y.max.rdm, col = alpha("chartreuse3",0.2), pch = 16)
        points(x.max, y.max.all, col = alpha("blue3",0.2), pch = 16)
        legend("bottomright", legend = round(y.max,0), bty = "n", cex = 2)
      }
    }
    dev.off()
  }
}
# export des graphes des pics alignes

re.align <- function(sp){
  # Question : Quelle masse ?
  erreur_saisie<-TRUE
  while(erreur_saisie == TRUE){
    print("Quel masse doit etre re-aligne ? (un entier entre 30 et 500)")
    mass.align <- scan(file="", what="",sep = " ", nmax=1) %>% as.numeric()
    i <- names(sp$indX) %>% w.equal(mass.align)
    if(length(i) == 1){
      erreur_saisie<-FALSE
    }else{
      print("Il n'y a aucun spectre qui peut etre aligne sur cette masse.")
      erreur_saisie<-TRUE
    }
  }
  p.max <- apply(sp$MS_n.al[,sp$indX[[i]]],1,which.max)
  y.max.all <- apply(sp$MS_n.al[,sp$indX[[i]]],1,max)
  x.max <- sp$xMS[sp$indX[[i]]][p.max]
  y.max <- max(y.max.all)
  y.max.rdm <- length(x.max) %>% runif(y.max*1.1, y.max*1.3)

  par(mar = c(3,2,3,.5), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0), xpd = FALSE, yaxt = "n")
  matplot(sp$xMS[sp$indX[[i]]], t(sp$MS_n.al[,sp$indX[[i]]]), type = "l",
          ylab = " ", xlab = " ", main = paste("re-alignement de la masse",mass.align),
          col = alpha(viridis(length(p.max)),0.3), lty = 1)
  points(x.max, y.max.all, col = "blue3", pch = 16)

  # QUestion : combien de peak ?
  erreur_saisie <- TRUE
  while(erreur_saisie == TRUE){
    print("Combien de pics estimes ? (2, 3 ou 4)")
    print("[voir graphe associe]")
    nb.pk <- scan(file="", what="",sep = " ", nmax=1) %>% as.numeric()
    if(length(nb.pk) == 1){
      erreur_saisie <- FALSE
    }else{
      print("Il faut un chiffre...")
      erreur_saisie <- TRUE
    }
  }

  rdm.pk <- TRUE
  while(rdm.pk == TRUE){
    p.k <- cbind(x.max,y.max.all) %>% scale() %>% t() %>% multiply_by(c(3,1)) %>% t() %>% kmeans(nb.pk*3)
    p.col <- rep(NA,length(p.max))
    res.k <- sapply(1:(nb.pk*3),list)
    x.k.mean <- rep(NA,nb.pk*3)
    x.k.med <- rep(NA,nb.pk*3)
    y.k.mean <- rep(NA,nb.pk*3)
    for(j in 1:(nb.pk*3)){
      p.col[which(p.k$cluster==j)] <- alpha(plasma(nb.pk*3),0.7)[j]
      eph <- which(p.k$cluster==j)
      x.k.mean[j] <- x.max[eph] %>% mean()
      x.k.med[j] <- x.max[eph] %>% median()
      y.k.mean[j] <- y.max.all[eph] %>% mean()
      res.k[[j]] <- data.frame("X" = x.max[eph], "Y" = y.max.all[eph]) %>% as.matrix()
      eph <- unique(res.k[[j]][,1]) %>% length()
      if(eph == 1){
        res.k[[j]][,1] <- nrow(res.k[[j]]) %>% runif(0,0.001) %>% add(res.k[[j]][,1])
      }
    }
    ellipse <- lapply(res.k, ellipsoidhull) %>% lapply(predict)

    yplot <- t(sp$MS_n.al[,sp$indX[[i]]])

    par(mar = c(3,2,3,.5), cex.main=2, cex.lab = 2, cex.axis = 2, bg = "ivory", mgp = c(3.5,1.5,0), xpd = FALSE, yaxt = "n")
    matplot(sp$xMS[sp$indX[[i]]], yplot, type = "l", ylim = c(0,max(yplot))*1.3,
            ylab = " ", xlab = " ", main = paste("re-alignement de la masse",mass.align),
            col = alpha(viridis(length(p.max)),0.3), lty = 1)
    aaa <- lapply(ellipse, polygon, col = alpha("chocolate3",.3), border = "chocolate3", lwd = 1.5)
    points(x.max, y.max.all, bg = p.col, pch = 21, cex = 1.2, col = alpha("black",0.2))
    points(x.k.mean, y.k.mean, col = "black", pch = LETTERS[1:j], cex = 1.8)
    points(x.k.mean, y.k.mean, col = "white", pch = LETTERS[1:j], cex = 1.6)
    legend("topleft", legend = paste0(LETTERS[1:j], " (size : ", p.k$size,")"), pch = 16,
           col = alpha(plasma(j),0.7), ncol = ceiling(nb.pk*.3), bty = "n", cex = 1)


    lst.pk <- list()
    let.pk <- rep(NA,nb.pk)
    # Question : /!\ faire Grrr
    g <- 0
    while(g < nb.pk){
      g <- add(g,1)
      cat(" Quelles groupes de points correspondent a l'alignement",g,"?\n (une ou plusieurs lettre(s) : ABC ... [voir graphe associe])")
      if(g==1){
        cat("Si vous voulez un nouveau tirage aleatoire pour mieux repartir les groupes,\n montrez les dents et dites \"Grrr\" ")
      }
      gr <- scan(file="", what="",sep = "#", nmax=1)
      if(gr == "Grrr"){
        g <- nb.pk
      }else{
        gr <- str_remove_all(gr, " ") %>% str_remove_all(",") %>% str_to_upper() %>% str_split("")
        let.pk[g] <- gr[[1]] %>% str_flatten()
        lst.pk <- c(lst.pk, list(sapply(gr, match, table = LETTERS, simplify = TRUE)))
        names(lst.pk)[g] <- paste0("Gr",g)
        lst.pk[[g]] <- sapply(lst.pk[[g]],w.equal,p.k$cluster) %>% unlist() %>% sort()
        rdm.pk <- FALSE
      }
    }
  }

  # Mass et decalage
  x.max.u <- unique(x.max) %>% sort()
  y.max.u <- max(yplot) %>% rep(length(x.max.u))
  segments(x.max.u, 0, y1 = y.max.u*1.1, col = alpha("black", 0.2), lty = 3, lwd = 1.5)
  points(x.max.u, y.max.u*1.1, pch = 16, cex = 1.2, col = alpha("chocolate",0.4))
  text(x.max.u, y.max.u*1.2, col = "black", labels = seq(1,length(x.max.u),1), cex = 0.5)

  x.ref <- rep(NA,nb.pk*2) %>% matrix(2, nb.pk, dimnames = list(c("index","x.pos"),ch(seq(1,nb.pk,1))))
  # Question : quelle reference ?
  for(g in 1:nb.pk){
    erreur_saisie<-TRUE
    while(erreur_saisie == TRUE){
      if(g == 1){
        print(paste0("Quelle est la position du Groupe1 (",let.pk[g],") ? (un chiffre entre 1 et ",length(x.max.u),")"))
      }else{
        print(paste0("Par rapport a la position du Groupe1 (en vert) ou doit se touver le Groupe",g," (",let.pk[g],") ?"))
      }
      eph <- scan(file="", what="",sep = " ", nmax=1) %>% as.numeric()

      if(length(eph) == 1){
        x.ref[1,g] <- eph
        x.ref[2,g] <- x.max.u[x.ref[1,g]]
        if(g == 1){
          segments(x.ref[2,g], 0, y1 = y.max.u*1.1, col = alpha("chartreuse3", 0.6), lty = 3, lwd = 1.5)
          points(x.ref[2,g], max(yplot)*1.1, pch = 21, cex = 1.2, col = "chocolate", bg = alpha("chartreuse3",0.6))
        }else{
          segments(x.ref[2,g], 0, y1 = y.max.u*1.1, col = alpha("darkblue", 0.6), lty = 3, lwd = 1.5)
          points(x.ref[2,g], max(yplot)*1.1, pch = 21, cex = 1.2, col = "chocolate", bg = alpha("darkblue",0.6))
        }
        erreur_saisie<-FALSE
      }else{
        print("Il faut un chiffre...")
        erreur_saisie<-TRUE
      }
    }
  }

  # alignement des groupes :
  mat.al <- sp$MS_n.al[, sp$indX[[i]]]
  for(g in 1:nb.pk){
    p.max <- apply(sp$MS_n.al[lst.pk[[g]], sp$indX[[i]]], 1, which.max)
    shift.med <- det_c(x.ref[2,g],sp$xMS[sp$indX[[i]]])
    p.shift <- subtract(p.max,shift.med)
    shift.wdth <- range(p.shift) %>% diff()
    shift.start <- subtract(p.shift, max(p.shift)) %>% abs()
    shift.end <- subtract(p.shift, min(p.shift))

    MS.ctr <- matrix(NA, nrow = length(p.max), ncol = add(length(sp$indX[[i]]), shift.wdth))
    for(j in 1:nrow(MS.ctr)){
      a <- runif(shift.start[j],0,1) %>% round(2)
      b <- runif(shift.end[j],0,1) %>% round(2)
      MS.ctr[j,] <- c(a, mat.al[lst.pk[[g]][j],], b)
    }

    a <- max(p.shift) %>% add(1)
    b <- length(sp$indX[[i]]) %>% add(a) %>% subtract(1)

    mat.al[lst.pk[[g]],] <- MS.ctr[,a:b]
  }

  par(mar = c(3,2,3,.5), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0), xpd = FALSE, yaxt = "n")
  matplot(sp$xMS[sp$indX[[i]]], t(mat.al), type = "l",
          ylab = " ", xlab = " ", main = paste("re-alignement de la masse",mass.align),
          col = alpha(viridis(nrow(mat.al)),0.3), lty = 1)

  # Graphe de control
  p.max <- apply(mat.al, 1, which.max)
  x.max <- sp$xMS[sp$indX[[i]]][p.max]
  y.max.all <- apply(mat.al,1,max)
  y.max <- max(y.max.all)
  y.max.rdm <- length(p.max) %>% runif(y.max*1.1, y.max*1.3)

  shift.med <- median(p.max)
  p.shift <- subtract(p.max,shift.med)
  shift.wdth <- range(p.shift) %>% diff()
  shift.start <- subtract(p.shift, max(p.shift)) %>% abs()
  shift.end <- subtract(p.shift, min(p.shift))

  title.tiff <- paste0("Figures/Alignement/re_align_at_mass_",mass.align,".tiff")
  tiff(title.tiff, width = 600, height = 400,units = "px")
  par(mar = c(5,5,3,0.1), cex.main=2, cex.lab = 2, cex.axis = 2,mgp = c(3.5,1.5,0),xpd = FALSE)
  matplot(sp$xMS[sp$indX[[i]]], t(mat.al), type = "l", ylim = c(0,y.max*1.3),
          xlab = "m/z", ylab = "Relative intensity (u.a.)",
          main = paste("re-alignement de la masse",mass.align),
          col = alpha(viridis(length(p.max)),0.3), lty = 1)
  abline(v = median(x.max), col = "midnightblue", lwd = 2)
  points(x.max, y.max.rdm, col = alpha(p.col,0.2), pch = 16)
  points(x.max, y.max.all, col = alpha(p.col,0.2), pch = 16)
  legend("bottomright", legend = round(y.max,0), bty = "n", cex = 2)
  dev.off()


  print("Mise en forme")
  sp$MS[,sp$indX[[i]]] <- mat.al

  sp <- c(sp,
          list("MS_m" = t(apply(mt_h5,1,mean_mat, MS = sp$MS)),
               "MS_max" = t(apply(mt_h5,1,max_mat, MS = sp$MS))))

  sp <- c(sp,
          "pk" = list("MS" = pk.list(sp$MS, sp),
                      "MS_m" = pk.list(sp$MS_m, sp),
                      "MS_max" = pk.list(sp$MS_max, sp)))
  print("Ok !")
}
# permet de re-aligner a la main les spectres a une masse donnee.
# /!\ a utiliser avec rigueur et intelligence.

create_local_MS <- function(MS, xMS){
  createMassSpectrum(mass = xMS, intensity = MS)
}
# creer une lsite d'objet 'spectra'

round_xMS <- function(sp){return(round(sp@mass,3))}
# arrondi au millieme

list_to_mat <- function(sp = spectra[[1]], vec = DA_mill){
  vec_mass <- sp@mass %>% round(3)
  vec_inte <- sp@intensity
  dup <- 1
  while(dup != 0){
    eph <- duplicated(vec_mass)
    eph <- which(eph == TRUE)
    if(length(eph)>0){
      eph <- which(vec_mass== vec_mass[eph[1]])
      vec_inte[eph] <- mean(vec_inte[eph])
      vec_inte <- vec_inte[-eph[-1]]
      vec_mass <- vec_mass[-eph[-1]]

      eph <- duplicated(vec_mass)
      eph <- which(eph == TRUE)
    }
    dup <- length(eph)
  }
  coord_l_DA <- match(vec_mass,vec)
  MS_vec <- rep(NA,length(vec))
  MS_vec[coord_l_DA] <- vec_inte
  MS_vec <- zoo(MS_vec, vec)
  MS_vec <- na.spline(MS_vec)
  # plot(DA_mill, MS_vec, type = "l", xlim = c(59,59.2))
  # points(DA_mill, MS_vec, pch = 16, col = "black")
  # points(vec_mass, vec_inte, pch = 16, col = "blue")

  return(MS_vec)
}
# passe d'un objet 'spectra' a une matrice d'intensite

plotSpectra <- function(xMS_unalign, MS_unalign, xMS_align, MS_align, range=c(137, 137.2), nbr_sp = 100){
  par(mfrow=c(1, 2))
  if(nbr_sp > nrow(MS_align)){nbr_sp <- nrow(MS_align)}

  eph <- seq(1, nrow(MS_unalign), length.out = nbr_sp) %>% round()
  color <- viridis(nbr_sp, 0.5)
  brn_a <- det_c(range[1],xMS_unalign)
  brn_b <- det_c(range[2],xMS_unalign)

  matplot(xMS_unalign[brn_a:brn_b], t(MS_unalign[eph, brn_a:brn_b]),
          main = paste("unalign spectra (mass",range[1],":",range[2],"Da)"),
          xlab = "MS unalign", ylab = "intensity (a.u.)",
          type = "l", col = color)

  brn_a <- det_c(range[1],xMS_align)
  brn_b <- det_c(range[2],xMS_align)

  matplot(xMS_align[brn_a:brn_b], t(MS_align[eph, brn_a:brn_b]),
          main = paste("align spectra (mass",range[1],":",range[2],"Da)"),
          xlab = "MS align", ylab = "intensity (a.u.)",
          type = "l", col = color)

  par(mfrow=c(1, 1))
}
# plot 2 spectre pour comparer les alignements

alignement_sp <- function(sp){
  spectra <- apply(sp$MS,1, create_local_MS, xMS = sp$xMS)
  spectra <- alignSpectra(spectra, halfWindowSize = 20, SNR = 2, tolerance = 0.002, warpingMethod="lowess")

  avgSpectra <- averageMassSpectra(spectra, labels= names(spectra), method="mean")
  peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=5, SNR=10)
  peaks <- binPeaks(peaks, tolerance=0.002)
  peaks <- filterPeaks(peaks, minFrequency=0.25)
  featureMatrix <- intensityMatrix(peaks, avgSpectra)

  DA_mill <- lapply(spectra, round_xMS) %>% unlist() %>% unique() %>% sort()
  xMS <- DA_mill %>% as.numeric()

  MS_align <- sapply(spectra, list_to_mat, vec = DA_mill)
  MS_align <- t(MS_align)
  colnames(MS_align) <- DA_mill
  return(list("MS_align" = MS_align, "peaks" = featureMatrix))
}
# aligne les spectres
  # zoom <- c(142.9,143.3) # c(59,59.1) # c(137,137.2) # c(141.9, 142.3)
  #
  # plotSpectra(xMS_unalign = sp$xMS, MS_unalign = sp$MS,
  #             xMS_align = xMS, MS_align = MS_align,
  #             range = zoom, nbr_sp = 50)

# Fonctions de mise en formes : ####

bucketing_xMS <- function(xMS, bn_x_m = 31, bn_x_M = 500){
  k = 0
  for (i in (bn_x_m + 1):(bn_x_M - 1)){
    if(length(which((xMS >= i-.5)&( xMS < i+.5))) > 0){
      k = k + 1
    }
  }
  indX <- vector("list", k)

  k = 0
  for (i in (bn_x_m + 1):(bn_x_M - 1)){
    if(length(which((xMS >= i-.5)&( xMS < i+.5))) > 0){
      k = k + 1
      indX[[k]] <- which((xMS >= i-.5)&( xMS < i+.5))
      names(indX)[k] <- i
    }
  }
  return(indX)
}
# retourne les indices de xMS pour chaque masse unitaire
# ou il y a un pic detecte.

extract_timing <- function(ls_h5 = ls_h5, MS = sp$MS){
  xT <- ls_h5[[1]]$TimingData$BufTimes[1,]
  if(length(ls_h5)>1){
    for(i in 2:length(ls_h5)){
      xT <- c(xT, ls_h5[[i]]$TimingData$BufTimes[1,])
    }
  }
  names(xT) <- row.names(MS)
  return(xT)
}
# retourne le moment d'acquisition de chaque ms a partir de T0.

# Fonction de detection des pics : ####
pk.list2 <- function(sp){
  pks <- sp$peaks
  val <- list()
  for(i in 1:nrow(pks)){
    val[[i]] <- boxplot.stats(pks[i,])[[4]] %>% sort(decreasing = TRUE)
  }
  eph <- lapply(val,length) %>% unlist() %>% min()
  wei <- matrix(NA, 4,eph, dimnames = list(c("W", "int", "pos", "coor"),1:eph))

  pkl <- list()
  for(i in 1:nrow(pks)){
    pk_mat <- boxplot.stats(pks[i,])[[4]] %>% sort(decreasing = TRUE)
    wei[1,] <- names(pk_mat[1:eph]) %>% as.numeric() %>% round(0)
    wei[2,] <- pk_mat[1:eph] %>% unname() %>% round(0)
    wei[3,] <- sapply(wei[1,], match, table = names(sp$indX))
    wei[4,] <- names(pk_mat[1:eph]) %>% as.numeric() %>% round(3)
    wei[4,] <- sapply(wei[4,], match, table = sp$xMS)
    pkl[[i]] <- wei
    names(pkl)[[i]] <- row.names(sp$MS)[i]
  }

 return(pkl)
}
# retourne une liste de peak pour chaque spectre de mat (adapte a MS)

pk.list3 <- function(mat, sp){
  spectra <- apply(mat, 1, create_local_MS, xMS = sp$xMS)
  avgSpectra <- averageMassSpectra(spectra, labels= names(spectra), method="mean")
  peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=5, SNR=10)
  peaks <- binPeaks(peaks, tolerance=0.002)
  peaks <- filterPeaks(peaks, minFrequency=0.25)
  pks <- intensityMatrix(peaks, avgSpectra)

  val <- list()
  for(i in 1:nrow(pks)){
    val[[i]] <- boxplot.stats(pks[i,])[[4]] %>% sort(decreasing = TRUE)
  }
  eph <- lapply(val,length) %>% unlist() %>% min()
  wei <- matrix(NA, 4,eph, dimnames = list(c("W", "int", "pos", "coor"),1:eph))

  pkl <- list()
  for(i in 1:nrow(pks)){
    pk_mat <- boxplot.stats(pks[i,])[[4]] %>% sort(decreasing = TRUE)
    wei[1,] <- names(pk_mat[1:eph]) %>% as.numeric() %>% round(0)
    wei[2,] <- pk_mat[1:eph] %>% unname() %>% round(0)
    wei[3,] <- sapply(wei[1,], match, table = names(sp$indX))
    wei[4,] <- names(pk_mat[1:eph]) %>% as.numeric() %>% round(3)
    wei[4,] <- sapply(wei[4,], match, table = sp$xMS)
    pkl[[i]] <- wei
    names(pkl)[[i]] <- row.names(sp$MS)[i]
  }
  return(pkl)
}
# retourne une liste de peak pour chaque spectre de mat (adapte a MS_max, MS_m, loading,...)

pk.list4 <- function(mat, sp){
  mat_p <- mat
  mat_p[which(mat<0)] <- 0

  mat_m <- mat %>% abs()
  mat_m[which(mat>0)] <- 0

  spectra <- apply(mat, 1, create_local_MS, xMS = sp$xMS)
  avgSpectra <- averageMassSpectra(spectra, labels= names(spectra), method="mean")
  peaks <- detectPeaks(avgSpectra, method="MAD", halfWindowSize=5, SNR=10)
  peaks <- binPeaks(peaks, tolerance=0.002)
  peaks <- filterPeaks(peaks, minFrequency=0.25)
  pks <- intensityMatrix(peaks, avgSpectra)

  val <- list()
  for(i in 1:nrow(pks)){
    val[[i]] <- boxplot.stats(pks[i,])[[4]] %>% sort(decreasing = TRUE)
  }
  eph <- lapply(val,length) %>% unlist() %>% min()
  wei <- matrix(NA, 4,eph, dimnames = list(c("W", "int", "pos", "coor"),1:eph))

  pkl <- list()
  for(i in 1:nrow(pks)){
    pk_mat <- boxplot.stats(pks[i,])[[4]] %>% sort(decreasing = TRUE)
    wei[1,] <- names(pk_mat[1:eph]) %>% as.numeric() %>% round(0)
    wei[2,] <- pk_mat[1:eph] %>% unname() %>% round(0)
    wei[3,] <- sapply(wei[1,], match, table = names(sp$indX))
    wei[4,] <- names(pk_mat[1:eph]) %>% as.numeric() %>% round(3)
    wei[4,] <- sapply(wei[4,], match, table = sp$xMS)
    pkl[[i]] <- wei
    names(pkl)[[i]] <- row.names(sp$MS)[i]
  }
  return(pkl)
}
# retourne une liste de peak pour chaque spectre de mat (adapte a MS_max, MS_m, loading,...)

pk.red <- function(pk_x = pk_mat[1,1], mat = pk_mat, w.sub = 4){
  eph <- subtract(mat[1,],pk_x) %>% abs() %>% multiply_by(-1) %>% which.sup(-(w.sub+1))
  return(eph[which.max(mat[2,eph])])
}
# trouve les pics proches

pk.short <- function(pk_mat){
  eph <- sapply(pk_mat[1,], pk.red, mat = pk_mat) %>% unique() %>% sort()
  pk_mat <- pk_mat[,eph]
  colnames(pk_mat) <- pk_mat[1,]
  return(pk_mat)
}
# reduit la liste des pics

sp_peaks <- function(sp){
  options(warn=-1)

  peaks_all <- pk.list2(sp) %>% lapply(pk.short)
  peaks_max <- pk.list3(sp$MS_max, sp) %>% lapply(pk.short)
  peaks_mean <- pk.list3(sp$MS_m, sp) %>% lapply(pk.short)

  options(warn=0)
  return(list("MS" = peaks_all,
              "MS_m" = peaks_mean,
              "MS_max" = peaks_max))
}
# fait les listes des peaks pour MS, MS_m et MS_max lors de l'import

# Fonctions de gestion des meta-data : ####
empty.meta <- function(mt_h5){
  eph <- c("names","acq_number", colnames(mt_h5), "use4analysis", "blank(acq_number)", "color",
           "concentration","unit","acq_T0", "delta_T", "grp1", "grp2", "...")
  eph <- matrix("", nrow = nrow(mt_h5), ncol = length(eph)-6) %>%
    cbind(row.names(mt_h5), 1:nrow(mt_h5), mt_h5, rep(TRUE,nrow(mt_h5)),.) %>%
    rbind(eph,.)
  write.table(eph, file = "meta_empty.csv", sep = ";", dec = ",", row.names = FALSE, col.names = FALSE)
}

import.meta <- function(nm){
  read.table(paste0(nm,".csv"), sep = ";", dec = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
}

export.para <- function(mt){
  # eph <- sapply(mt,length) %>% w.equal(1)
  eph <- 1:length(mt)
  eph <- eph[-match("h5",names(mt))]

  mat <- matrix(NA,length(eph),2)
  for(i in 1:length(eph)){
    if(length(mt[[eph[i]]]) == 1){
      mat[i,] <- c(names(mt)[eph[i]],mt[[eph[i]]])
    }else if(length(mt[[eph[i]]]) > 1){
      arg <- mt[[eph[i]]] %>% ch() %>% str_flatten(" ")
      mat[i,] <- c(names(mt)[eph[i]],arg)
    }
  }
  nm <- dir() %>% grep("PTR_MS_para",.) %>% length() %>% add(1) %>% paste0("PTR_MS_para_",.,".csv")
  write.table(mat, file = nm, sep = ";", dec = ".", row.names = FALSE, col.names = FALSE)
}

import.para <- function(nm, para = mt){
  eph <- read.table(paste0(nm,".csv"), sep = ";", row.names = 1, stringsAsFactors = FALSE)
  for(i in 1:nrow(eph)){
    j <- match(rownames(eph)[i], names(para))

    if(str_detect(eph[i,1]," ")==FALSE){
      para[[j]] <- eph[i,1]
    }else{
      para[[j]] <- str_split(eph[i,1]," ", simplify = TRUE) %>% as.numeric()
    }
  }
  return(para)
}

# Fonctions de gestions du temps d'acquisition : ####

acq_time <- function(ls = ls_h5[[1]]){
  oldw <- getOption("warn")
  options(warn = -1)
  eph <- ls$AcquisitionLog$Log$timestring[1]
  Th <- str_sub(eph,12,13) %>% as.numeric()
  Tm <- str_sub(eph,15,16) %>% as.numeric() %>% divide_by(60)
  Ts <- str_sub(eph,18,19) %>% as.numeric() %>% divide_by(3600)
  T_acq <- Th + Tm + Ts

  options(warn = oldw)
  return(T_acq)
}
# retourne le T0 de chaque acquisition.

acq_day <- function(ls = ls_h5[[1]]){
  oldw <- getOption("warn")
  options(warn = -1)
  eph <- ls$AcquisitionLog$Log$timestring[1]
  eph <- str_sub(eph,1,10) %>% str_split("-", simplify = TRUE)
  day_by_month <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  eph <- sum(day_by_month[1:(as.numeric(eph[2])-1)]) %>% add(as.numeric(eph[3]))
  options(warn = oldw)
  return(eph)
}
# retourne le D0 de chaque acquisition.

re.init.T.para <- function(sp){
  sp <- c(sp,
          list("xT" = sp$xT_init,
               "Tacq" = sp$Tacq_init))
  return(sp)
}
# reinitialise le temps

re.calc.T.para <- function(sp, mt){

  #Tzero <- mt$h5[,"acq_T0"] %>% unlist() %>% paste0("_1") %>% sp$xT_init[.]

  eph <- which.na(mt$h5[,"acq_T0"])
  if(length(eph) > 0){  mt$h5[eph,"acq_T0"] <- rownames(mt$h5)[eph] }

  eph <- which(mt$h5[,"acq_T0"] == "")
  if(length(eph) > 0){  mt$h5[eph,"acq_T0"] <- rownames(mt$h5)[eph] }

  eph <- which.na(mt$h5[,"delta_T"])
  if(length(eph) > 0){  mt$h5[eph,"delta_T"] <- 0 }

  eph <- which(mt$h5[,"delta_T"] == "")
  if(length(eph) > 0){  mt$h5[eph,"delta_T"] <- 0 }

  ref_T0 <- sapply(mt$h5[,"acq_T0"], grep, x=rownames(mt$h5))

  for(i in mt$acq){
    eph <- seq(mt$h5$start[i],mt$h5$end[i],by = 1)
    d.T <- (sp$Tacq[i] - sp$Tacq[ref_T0[i]]) + (sp$Dacq[i] - sp$Dacq[ref_T0[i]])*24
    d.T <- d.T*60*60

    sp$xT[eph] <- add(sp$xT_init[eph], d.T) %>% add(mt$h5[,"delta_T"][i])
  }
  return(sp)

}
#recalcule le temps a partir des references de la colonne "mt$h5$XXX"

# Fonctions d'analyses (visualisation) : ####
ctrl_plot_fct <- function(ls_h5 = ls_h5, sp = sp, mt = mt){

  if(("Control" %in% dir("Figures"))==FALSE){
    dir.create("Figures/Control")
  }

  tiff(filename = paste0("Figures/Control/Ctrl_plot_red_",mt$ctrl_plot_xmin,"_",mt$ctrl_plot_xmax,".tif"),
       width = 1000, height = 580)
  matplot(sp$xMS,t(sp$MS),type = "l", col = alpha("red",0.1),
          xlim = c(mt$ctrl_plot_xmin, mt$ctrl_plot_xmax),
          ylim = c(mt$ctrl_plot_ymin, mt$ctrl_plot_ymax),
          xlab = "m/z", ylab = "Relative intensity (u.a.)")
  dev.off()

  tiff(filename = paste0("Figures/Control/Ctrl_plot_raw_",mt$ctrl_plot_xmin,"_",mt$ctrl_plot_xmax,".tif"),
       width = 1000, height = 580)
  matplot(0,0, type = "l",
          xlim = c(mt$ctrl_plot_xmin, mt$ctrl_plot_xmax),
          ylim = c(mt$ctrl_plot_ymin, mt$ctrl_plot_ymax),
          xlab = "m/z", ylab = "Raw intensity (cps)")
  for(j in 1:length(ls_h5)){
    matplot(ls_h5[[j]]$FullSpectra$MassAxis, ls_h5[[j]]$FullSpectra$TofData[,1,1,],
            type = "l", col = alpha("black",0.1), add = TRUE)
    print(paste(j,"/",length(ls_h5)))
  }
  dev.off()
}
# permet de verifier si le spectre reduit correspond au spectre brut.

view_plot_fct <- function(ls_h5 = ls_h5, sp = sp, mt = mt){
  if(("Produits pures" %in% dir("Figures"))==FALSE){
    dir.create("Figures/Produits pures")
  }

  if((mt$view_plot_peak[1] == FALSE)&(mt$view_each.group == FALSE)){
    xmin <- mt$view_plot_xmin
    xmax <- mt$view_plot_xmax
    cmin <- det_c(xmin,sp$xMS)
    cmax <- det_c(xmax,sp$xMS)
    for(j in 1:length(ls_h5)){
      pk.MS <- sp$pk$MS_max[[j]]
      titre <- paste0("Figures/Produits pures/",row.names(mt$h5)[j],"_",xmin,"_",xmax,".tif")
      ymax <- max(pk.MS[2,])*1.1

      tiff(filename = titre, width = 1000, height = 580)
        par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

            matplot(sp$xMS[cmin:cmax], sp$MS_max[j,cmin:cmax],
                    type = "l", col = "grey95",
                    ylim = c(0,ymax), lwd = 5,
                    xlab = "m/z", ylab = "Relative intensity (u.a.)",
                    main = row.names(mt$h5)[j], xaxt="n")
            rect(0,-ymax*2, xmax*2, ymax*2, col="grey90")

            matplot(sp$xMS[cmin:cmax], sp$MS_max[j,cmin:cmax], add = TRUE,
                    type = "l", col = alpha("turquoise2",0.5), lwd = 5)

            matplot(sp$xMS[cmin:cmax], t(sp$MS[mt$h5[j,"start"]:mt$h5[j,"end"],cmin:cmax]),
                    type = "l", col = alpha("grey50",0.3), add = TRUE)

            axis(side=1,0:600, tcl=0,labels=FALSE)
            axis(side=2,(-ymax:ymax)*2,tcl=0,labels=FALSE)
            axis(side=3,0:600, tcl=0,labels=FALSE)
            axis(side=4,(-ymax:ymax)*2,tcl=0,labels=FALSE)

            axis(1, at = seq(dizaine(xmin), dizaine(xmax), 10), lwd.ticks = 2, tck = -0.03)
            axis(1, at = seq(dizaine(xmin), dizaine(xmax) +10, 5), labels = FALSE, tck = -0.03)
            axis(1, at = seq(dizaine(xmin), dizaine(xmax) +10, 1), labels = FALSE, tck = -0.01)
            text(pk.MS[1,], pk.MS[2,], labels = pk.MS[1,], cex = 0.8, pos = 3, offset = 0.5)
      dev.off()
    }
  }

  if((mt$view_plot_peak[1] == FALSE)&(mt$view_each.group != FALSE)){

    xmin <- mt$view_plot_xmin
    xmax <- mt$view_plot_xmax
    cmin <- det_c(xmin,sp$xMS)
    cmax <- det_c(xmax,sp$xMS)

    grp <- mt$view_each.group %>% mt$h5[,.]
    grp[grp == ""] <- NA
    grp_u <- unique(grp)
    if(length(which.na(grp_u))>0){grp_u <- grp_u[-which.na(grp_u)]}

    for(g in grp_u){
      ind_G <- which(grp == g)
      ind_mt <- mt$h5[ind_G[1],"start"]:mt$h5[ind_G[1],"end"]
      pk.MS <- sp$pk$MS_max[[ind_G[1]]]

      for(f in ind_G[-1]){
        ind_mt <- c(ind_mt, mt$h5[f,"start"]:mt$h5[f,"end"])
        pk.MS <- cbind(pk.MS, sp$pk$MS_max[[f]])
      }

      pk.MS <- pk.short(pk.MS)

      titre <- paste0("Figures/Produits pures/",row.names(mt$h5)[ind_G[1]],"_",g,"_",xmin,"_",xmax,".tif")
      ymax <- max(pk.MS[2,])*1.1

      tiff(filename = titre, width = 1000, height = 580)
      par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

      if(length(ind_G)==1){ind_G <- c(ind_G,ind_G)}

      matplot(sp$xMS[cmin:cmax], t(sp$MS_max[ind_G,cmin:cmax]),
              type = "l", col = "grey95",
              ylim = c(0,ymax), lwd = 5,
              xlab = "m/z", ylab = "Relative intensity (u.a.)",
              main = paste0(row.names(mt$h5)[ind_G[1]],"_",g), xaxt="n")
      rect(0,-ymax*2, xmax*2, ymax*2, col="grey90")

      matplot(sp$xMS[cmin:cmax], t(sp$MS_max[ind_G,cmin:cmax]), add = TRUE,
              type = "l", col = alpha("turquoise2",0.5), lwd = 5)

      matplot(sp$xMS[cmin:cmax], t(sp$MS[ind_mt,cmin:cmax]),
              type = "l", col = alpha("grey50",0.3), add = TRUE)

      axis(side=1,0:600, tcl=0,labels=FALSE)
      axis(side=2,(-ymax:ymax)*2,tcl=0,labels=FALSE)
      axis(side=3,0:600, tcl=0,labels=FALSE)
      axis(side=4,(-ymax:ymax)*2,tcl=0,labels=FALSE)

      axis(1, at = seq(dizaine(xmin), dizaine(xmax), 10), lwd.ticks = 2, tck = -0.03)
      axis(1, at = seq(dizaine(xmin), dizaine(xmax) +10, 5), labels = FALSE, tck = -0.03)
      axis(1, at = seq(dizaine(xmin), dizaine(xmax) +10, 1), labels = FALSE, tck = -0.01)
      text(pk.MS[1,], pk.MS[2,], labels = pk.MS[1,], cex = 0.8, pos = 3, offset = 0.5)
      dev.off()
    }
  }

  if((mt$view_plot_peak[1] != FALSE)&(mt$view_each.group != FALSE)){
    grp <- mt$view_each.group %>% mt$h5[,.]
    grp[grp == ""] <- NA
    grp_u <- unique(grp)
    if(length(which.na(grp_u))>0){grp_u <- grp_u[-which.na(grp_u)]}

    for(g in grp_u){
      ind_G <- which(grp == g)
      ind_mt <- mt$h5[ind_G[1],"start"]:mt$h5[ind_G[1],"end"]

      for(f in ind_G[-1]){
      ind_mt <- c(ind_mt, mt$h5[f,"start"]:mt$h5[f,"end"])
    }

      for(i in mt$view_plot_peak){
      pp <- list()
      pp$Mu <- i
      pp$xmin <- pp$Mu - 0.5
      pp$xmax <- pp$Mu + 0.5
      pp$ind_Mu <- which(names(sp$indX)==ch(pp$Mu)) %>% sp$indX[.] %>% unlist() %>% unname()
      pp$mat_Mu <- sp$MS[ind_mt,pp$ind_Mu] %>% t()
      pp$mat_X <- sp$xMS[pp$ind_Mu]
      pp$titre <- paste0("Figures/Produits pures/",g,"_MS_",pp$Mu,".tif")
      pp$main <- g
      pp$col <- ncol(pp$mat_Mu) %>% viridis(0.3)

      plot_peak(pp)
    }
    }
  }

  if((mt$view_plot_peak[1] != FALSE)&(mt$view_each.group == FALSE)){
    for(j in mt$acq){
      for(i in mt$view_plot_peak){
        pp <- list()
        pp$Mu <- i
        pp$xmin <- pp$Mu - 0.5
        pp$xmax <- pp$Mu + 0.5
        pp$ind_Mu <- which(names(sp$indX)==ch(pp$Mu)) %>% sp$indX[.] %>% unlist() %>% unname()
        pp$mat_Mu <- sp$MS[mt$h5[j,"start"]:mt$h5[j,"end"],pp$ind_Mu] %>% t()
        pp$mat_X <- sp$xMS[pp$ind_Mu]
        pp$titre <- paste0("Figures/Produits pures/",row.names(mt$h5)[j],"_MS_",pp$Mu,".tif")
        pp$main <- row.names(mt$h5)[j]
        pp$col <- ncol(pp$mat_Mu) %>% viridis(0.3)

        plot_peak(pp)
      }
    }
  }

}
# creer des graphes de visualisation des acquistions

view_one_plot_fct <- function(n_acq = mt$view_plot, ls_h5, sp, mt){
  if(("Produits pures" %in% dir("Figures"))==FALSE){
    dir.create("Figures/Produits pures")
  }

  pk.MS <- sp$pk.MS[[n_acq]]
  if(is.null(ncol(pk.MS)) == TRUE){
    pk.MS <- rbind(pk.MS, pk.MS) %>% t()
  }

  if(mt$view_plot_peak[1] == FALSE){
    n_ex <- which((mt$h5$start <= n_acq)&(mt$h5$end >= n_acq))
    xmin <- mt$view_plot_xmin
    xmax <- mt$view_plot_xmax
    eph <- det_c(xmin,sp$xMS):det_c(xmax,sp$xMS)
    ymax <- max(sp$MS[n_acq,eph])*1.1
    titre <- paste0("Figures/Produits pures/",row.names(mt$h5)[n_ex],"_acq_",n_acq,"_from",
                    xmin,"_to",xmax,".tif")

    tiff(filename = titre, width = 1000, height = 580)
      par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))
        matplot(sp$xMS, sp$MS[n_acq,], xaxt="n",
                type = "l", col = mt$h5$color[n_ex], lwd = 2,
                xlim = c(xmin, xmax), ylim = c(0, ymax),
                xlab = "m/z", ylab = "Relative intensity (u.a.)",
                main = row.names(mt$h5)[n_ex])
        axis(1, at = seq(dizaine(xmin), dizaine(xmax), 10), lwd.ticks = 2, tck = -0.03)
        axis(1, at = seq(dizaine(xmin), dizaine(xmax)+10, 5), labels = FALSE, tck = -0.03)
        axis(1, at = seq(dizaine(xmin), dizaine(xmax)+10, 1), labels = FALSE, tck = -0.01)
        text(pk.MS[1,], pk.MS[2,], labels = pk.MS[1,], cex = 0.8, pos = 3, offset = 0.5)
    dev.off()
  }else{
    for(i in 1:length(mt$view_plot_peak)){
      n_ex <- which((mt$h5$start <= n_acq)&(mt$h5$end >= n_acq))
      xmin <- mt$view_plot_peak[i] - 0.5
      xmax <- mt$view_plot_peak[i] + 0.5
      eph <- det_c(xmin,sp$xMS):det_c(xmax,sp$xMS)
      ymax <- max(sp$MS[n_acq,eph])*1.1

      titre <- paste0("Figures/Produits pures/",row.names(mt$h5)[n_ex],"_acq_",n_acq,"_MS_",mt$view_plot_peak[i],".tif")

      tiff(filename = titre, width = 1000, height = 580)
        par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))
          matplot(sp$xMS, sp$MS[n_acq,], xaxt="n",
                  type = "l", col = mt$h5$color[n_ex], lwd = 2,
                  xlim = c(xmin, xmax), ylim = c(0, ymax),
                  xlab = "m/z", ylab = "Relative intensity (u.a.)",
                  main = row.names(mt$h5)[n_ex])
          axis(1, at = mt$view_plot_peak, labels = FALSE, lwd.ticks = 2, tck = -0.03)
          axis(1, at = seq(xmin, xmax, 0.1), tck = -0.03)
          axis(1, at = seq(xmin, xmax, 0.01), labels = FALSE, tck = -0.01)
      dev.off()
    }
  }
}
# creer des graphes de visualisation des acquistions

plot_peak <- function(pp){
  tiff(filename = pp$titre, width = 1000, height = 580)
  par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))
  matplot(pp$mat_X, pp$mat_Mu,
          type = "l", col = pp$col,
          xlim = c(pp$xmin, pp$xmax),
          xlab = "m/z", ylab = "Relative intensity (u.a.)",
          main = pp$main, xaxt="n")
  axis(1, at = pp$Mu, labels = FALSE, lwd.ticks = 2, tck = -0.03)
  axis(1, at = seq(pp$xmin, pp$xmax, 0.1), tck = -0.03)
  axis(1, at = seq(pp$xmin, pp$xmax, 0.01), labels = FALSE, tck = -0.01)
  dev.off()
}
# creer un graphe centre sur une masse unitaire

print_calib <- function(mass_mol, name_mol, x_ax = xMS, y_mat = all_MS, addProton = TRUE){

  if(addProton == TRUE){
    mass_mol <- add(mass_mol,1.00782)
  }
  mass_mol <- round(mass_mol,4)
  x_m <- mass_mol - 0.35
  x_M <- mass_mol + 0.35
  xbr_m <- det_c(x_m,x_ax)
  xbr_M <- det_c(x_M,x_ax)

  if(("Control" %in% dir("Figures"))==FALSE){
    dir.create("Figures/Control")
  }

  ybr_M <- apply(y_mat[,xbr_m:xbr_M],1,which.max) %>% add(xbr_m) %>% subtract(1)
  xbr_m_red <- min(ybr_M) - 5
  xbr_M_red <- max(ybr_M) + 5

  eph <- x_ax[ybr_M] %>% median()
  if(eph <= mass_mol){
    par_fig <- c(0.65,0.99, 0.5, 1)
  }else{
    par_fig <- c(0.07,0.4, 0.5, 1)
  }
  err <- subtract(eph, mass_mol) %>% abs() %>% round(4)

  tiff(file = paste0("Figures/Control/Calibration_control_",name_mol,".tiff"), width = 700, height = 450,units = "px")
    # graphe general :
    par(mar = c(5,5,3,0.1), cex.main=2, cex.lab = 2, cex.axis = 2,mgp = c(3.5,1.5,0),xpd = FALSE)

      matplot(c(x_m,x_M),c(0,1), type = "l", col = "white",
              xlab = "m/z", ylab = "Relative intensity (u.a.)",
              main = paste(name_mol,"+ H+ (Mass : ", str_replace(ch(mass_mol),"\\.",","),
                           "; error : ", str_replace(ch(err),"\\.",","),")"))

      for(i in 1:nrow(y_mat)){
        matplot(x_ax[xbr_m:xbr_M], y_mat[i,xbr_m:xbr_M]/max(y_mat[i,xbr_m:xbr_M]),
                type = "l", col = alpha(viridis(nrow(y_mat))[i],0.3), add = TRUE)
      }
      abline(v = mass_mol, lty = 2, lwd = 3)
      abline(v = eph, lty = 3, lwd = 2)

      expression()

    # insert :
    par(fig = par_fig,mar = c(0.1,3,3.5,0.1),
        new = T, bty = "o", yaxt = "n", cex = 0.8, mgp = c(2,1,0))

      matplot(c(x_ax[xbr_m_red],x_ax[xbr_M_red]),c(0.7,1), type = "l", col = "white",
              xlab = " ", ylab = " ")
      for(i in 1:nrow(y_mat)){
        matplot(x_ax[xbr_m_red:xbr_M_red], y_mat[i,xbr_m_red:xbr_M_red]/max(y_mat[i,xbr_m_red:xbr_M_red]),
                type = "l", col = alpha(viridis(nrow(y_mat))[i],0.3), add = TRUE)
      }
      abline(v = mass_mol, lty = 2, lwd = 3)
      abline(v = eph, lty = 3, lwd = 2)
  dev.off()
}
# returne un graphe centre sur mass_mol

print_spectre_circle <- function(num_acq, sp, mode = "threshold", thr = 1000){
  data <- colnames(sp$MS.AUC) %>% as.numeric() %>% t() %>% as.data.frame()
  dimnames(data) <- dimnames(sp$MS.AUC[num_acq,])

  data <- rbind(data, sp$MS.AUC[num_acq,]) %>% t() %>% as.data.frame()
  colnames(data) <- c("mass", "Value")
  data$Value <- abs(data$Value)

  titre <- rownames(sp$MS.AUC)[num_acq]
  # Order data
  tmp <- data %>% arrange(desc(Value)) %>% mutate(mass = factor(mass, mass))

  if(mode == "threshold"){
    eph <- which(tmp$Value >= thr)
    tmp <- tmp[eph,]
    number_of_bar <- nrow(tmp)
  }else if(mt$view_circular_mode == "bar"){
    tmp <- tmp[1:thr,]
    number_of_bar <- nrow(tmp)
  }

  tmp$Value <- round(tmp$Value,0)
  tmp$id <- seq(1, nrow(tmp))
  tmp$Lab <- paste0(tmp$mass," (",tmp$Value,")")
  angle <- 90 - 360 * (tmp$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  tmp$hjust <- ifelse( angle < -90, 1, 0)
  tmp$angle <- ifelse(angle < -90, angle+180, angle)
  tmp$Value_raw <- tmp$Value
  eph <- which(tmp$Value>9999)
  if(length(eph) >=1 ){
    tmp$Value[eph] <- 9999
  }

  eph <- max(tmp$Value_raw) %>% viridis(0.8, direction = -1)
  tmp$color <- eph[tmp$Value_raw]

  eph <- divide_by(number_of_bar,2) %>% round(0)
  tmp_lab <- tmp[eph,]
  tmp_lab$Lab <- titre

  # Make the plot
  ggplot(tmp, aes(x=as.factor(id), y=Value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity", fill= tmp$color) +
    ylim(-4000,18000) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")
    ) +
    coord_polar(start = 0) +
    geom_text(data = tmp, aes(x = id, y = Value+600, label = Lab), color = "black", fontface = "bold",
              alpha = 0.6, size = 4, angle = tmp$angle, hjust = tmp$hjust, inherit.aes = FALSE ) +
    geom_text(data = tmp_lab, aes(x = id, y = Value+12000, label = Lab), color = "black", fontface = "bold",
              alpha = 0.6, size = 4,  inherit.aes = FALSE )

  paste0("Figures/Circular/Sp_circle_",titre,".png") %>% ggsave(width = 3.5, height = 4)
}
# creer des graphes circualires pour visualiser les data les plus importantes

# Fonctions d'analyses (graphe dynamique) : ####

dy.spectra <- function(sp, mode = "align", select.sp = NULL){

  if(mode == "max"){
    sp_sel <- sp$MS_max
    xMS <- sp$xMS
    titre <- "sp_max"
  }else if(mode == "mean"){
    sp_sel <- sp$MS_m
    xMS <- sp$xMS
    titre <- "sp_mean"
  }else if(mode == "raw"){
    sp_sel <- sp$MS_n_al
    xMS <- sp$xMS_n_al
    titre <- "sp_not_align"
  }else{
    sp_sel <- sp$MS
    xMS <- sp$xMS
    titre <- "sp_align"
  }

  if(is.null(select.sp) == TRUE){
    sel <- demande.print.spectra(sp_sel)
  }else if(select.sp == "all"){
    sel <- 1:nrow(sp_sel)
  }else{
    sel <- match(select.sp, rownames(sp_sel))
  }

  eph <- paste0(getwd(),"/Figures/") %>% dir() %>% grep(titre, .) %>% length() %>% add(1)
  eph <- paste0(getwd(),"/Figures/dy_",titre,"_",eph)

  if(length(sel)==1){
    dysp <- sp_sel[sel,] %>% cbind(xMS,.) %>% as.data.frame()
    colnames(dysp)[2] <- rownames(sp_sel)[sel]
  }else{
    dysp <- t(sp_sel[sel,]) %>% cbind(xMS,.) %>% as.data.frame()
  }
  rownames(dysp) <- xMS

  rmarkdown::render(paste0(fct_dir,"/print_dyspectre.Rmd"), output_file = eph)

  return(colnames(dysp)[-1])
}

demande.print.spectra <- function(sp_sel){
  erreur_demande <- TRUE
  while(erreur_demande == TRUE){
    erreur_saisie <- TRUE
    hd <- head(rownames(sp_sel)) %>% str_flatten(", ")
    tl <- tail(rownames(sp_sel)) %>% str_flatten(", ")

    while (erreur_saisie == TRUE){
      cat("Il y a plusieurs MS :\n\n", hd,"\n\n [...] \n\n", tl,
          "\n\nCombien en analyse-t-on ? \nSaisi un chiffre entre 1 et", nrow(sp_sel),", puis 'Enter' :")
      bl <- scan(file = "", what = "", sep = " ", nmax = 1) %>% as.numeric()

      if((max(bl) <= nrow(sp_sel))&(min(bl)>0)){ # verification que le nombre de fichiers selectionnes ne soit pas trop grand
        erreur_saisie <- FALSE
      }else{ # sinon on se fout un peu de la gueule du user et on recommence
        cat("   O.0' Combien ?\n")
        erreur_saisie <- TRUE
      }
    }

    # a partir de le on va recuperer 'cl' alias les coordonnees des data a analyser
    if(bl < nrow(sp_sel)){ #plusieurs fichiers... on sent que ca va gaver
      erreur_saisie <- TRUE
      while(erreur_saisie == TRUE){
        if(bl == 1){
          cat("Saisi le numero de l'acquisition que tu veux, puis 'Enter'.\n")
        }else{
          cat(paste("Saisi les", bl, "numeros des acquisitions que tu veux. 'Enter' entre chaque numero.\n"))
        }

        paste0(1:nrow(sp_sel), ": ", rownames(sp_sel)," \\\n") %>% cat()
        cl <- scan( file = "", what = "", sep = " ", nmax = bl) %>% as.numeric()
        erreur_saisie <- FALSE
      }
    }else{cl <- 1:bl}

    cl <- sort(cl) %>% unique()
    eph <- cl %in% 1:nrow(sp_sel) %>% grep(FALSE,.)
    if(length(eph) == 0){
      cat("Ok. C'est parti pour :\n")
      cat(rownames(sp_sel)[cl])

      erreur_saisie <- TRUE
      while(erreur_saisie == TRUE){
        cat("\nC'est bien ca ? (ok/nop)")
        your.ok <- scan( file = "", what = "", sep = " ", nmax = 1)
        if(your.ok == "ok"){
          erreur_saisie <- FALSE
          erreur_demande <- FALSE
        }else if(your.ok == "nop"){
          erreur_saisie <- FALSE
          cat("Hum... du coup on reprend.\n")
        }
      }
    }else if(length(eph) > 0){ # la, c'est que ca merde...
      cat(" Il y a un probleme ! -_-'\n On reprend ...")
    }
  }
  return(cl)
}

# Fonctions d'analyses (blanc) : ####

blank_subtraction <- function(sp, mt){
  index_raw <- which.not.na(mt$h5$blank.acq_number.)
  index_blk <- mt$h5$blank.acq_number.[index_raw]

  if(length(sp$blk_sub) == 0){ sp$blk_sub <- FALSE}

  if(sp$blk_sub == FALSE){
    sp <- c(sp, list("MS_raw" = sp$MS, "MS_m_raw" = sp$MS_m))
    sp$blk_sub <- TRUE
  }

  for(i in 1:length(index_raw)){
    ir <- index_raw[i]
    ib <- index_blk[i]
    sp$MS_m[ir,] <- sp$MS_m[ir,] - sp$MS_m[ib,]
    cir <- coor_mt(ir)
    sp$MS[cir,] <- sp$MS[cir,] - sp$MS_m[ib,]
  }
  cir <- coor_mt(index_raw)
  rownames(sp$MS_m[index_raw,]) <- paste0(rownames(sp$MS_m[index_raw,]),"_bs")
  rownames(sp$MS[cir,]) <- paste0(rownames(sp$MS[cir,]),"_bs")

  sp$pk.MS = pk.list(sp$MS, sp)
  sp$pk.MS_m = pk.list(sp$MS_m, sp)
  sp$pk.MS_max = pk.list(sp$MS_max, sp)

  if(("Blank" %in% dir("Figures"))==FALSE){
    dir.create("Figures/Blank")
  }

  for(ir in index_raw){
    tiff(filename = paste0("Figures/Blank/", row.names(sp$MS_m)[ir], "_blank.tif"), width = 1000, height = 580)
    matplot(sp$xMS, sp$MS_m[ir,],
            type = "l", col = alpha("grey50",0.3),
            xlim = c(mt$blanck_plot_xmin, mt$blanck_plot_xmax),
            xlab = "m/z", ylab = "Relative intensity (u.a.)",
            main = row.names(sp$MS_m)[ir])
    text(sp$pk.MS_m[[ir]][1,], sp$pk.MS_m[[ir]][2,], labels = sp$pk.MS_m[[ir]][1,], cex = 0.8, pos = 3, offset = 0.5)
    dev.off()
  }
  return(sp)
}
# supprime le blanc moyen

reverse_blank_subtraction <- function(sp){
  if(sp$blk_sub == TRUE){
    sp$MS <- sp$MS_raw
    sp$MS_m <- sp$MS_m_raw
    sp$blk_sub <- FALSE
  }
  return(sp)
}
# reviens aux spectres d'origines (avant soustraction des blancs)

# Fonctions d'analyses (AUC) : ####

peak_AUC <- function(sp, mt){
  MS.AUC <- matrix(NA, nrow = nrow(sp$MS), ncol = length(sp$indX),
                   dimnames = list(rownames(sp$MS),names(sp$indX)))
  for (i in 1:nrow(sp$MS)){
    for (j in 1:length(sp$indX)){
      MS.AUC[i,j] <- AUC(sp$xMS[sp$indX[[j]]],sp$MS[i,sp$indX[[j]]]) %>% round(2)
    }
  }
  return(as.data.frame(MS.AUC))
}
# retourne une matrice mesurant l'AUC de chaque masse detectee.

plot_AUC <- function(sp, mt, acq = mt$acq, M.num = mt$M.num){
  Tbn <- c(0,0)
  Ibn <- c(0,0)
  for(i in acq){
    Tbn[2] <- max(Tbn[2],sp$xT[indS(row.names(mt$h5)[i])])
    for(j in M.num){
      Ibn[2] <- max(Ibn[2],sp$MS.AUC[indS(row.names(mt$h5)[i]),ch(j)])
    }
  }

  unit = "s"
  Tdiv = 1
  if(Tbn[2]>600){
    unit = "min"
    Tdiv = 60
  }
  if(Tbn[2]>36000){
    unit = "h"
    Tdiv = 3600
  }

  AUC_col <- indColor(acq, mt$h5)

  tiff(file = mt$titre, width = 1200, height = 600,units = "px")
    par(mar = c(5,5,2,16), cex.main=2, cex.lab = 2, cex.axis = 2,mgp = c(3.5,1.5,0),xpd = NA)

      matplot(Tbn/Tdiv, Ibn, type = "l", col = "white", xlab = paste0("time (",unit,")"), ylab = "AUC (a.u.)", log = mt$AUC_plot_exp)
      nq <- 0
      for(i in acq){
        cl <- 15
        nq <- nq + 1
        for(j in M.num){
          coor <- indS(row.names(mt$h5)[i])
          Mnam <- ch(j)
          if(length(coor)>1){
            matplot(sp$xT[coor]/Tdiv,sp$MS.AUC[coor,Mnam], type = "l", lwd = 2,
                    col = AUC_col[nq], add = TRUE)
          }else{
            matplot(sp$xT[coor]/Tdiv,sp$MS.AUC[coor,Mnam], pch = cl, col = AUC_col[nq], add = TRUE)
          }

          # points(sp$xT[coor]/Tdiv,sp$MS.AUC[coor,Mnam], pch = cl, cex = 0.7, col = AUC_col[nq])
          cl <- cl + 1
        }
      }

      if(length(acq)<= 10){
        l.acq <- acq
      }else{
        eph <- length(acq) %>% subtract(4)
        l.acq <- acq[c(1:5,eph:length(acq))]
      }

      if(length(M.num)<= 10){
        l.num <- M.num
      }else{
        eph <- length(M.num) %>% subtract(4)
        l.num <- M.num[c(1:5,eph:length(M.num))]
      }

      legend("topright", bty = "n", cex = 1.5, xpd = NA, inset = c(-0.26,0),
             legend = c("Sample(s) :",row.names(mt$h5)[l.acq]," ","Masse(s) :", l.num),
             lty = c(NA,rep(1,length(l.acq)), NA, NA, rep(NA,length(l.num))), lwd = 2,
             pch = c(NA,rep(NA,length(l.acq)), NA, NA, 14 + seq(1,length(l.num))),
             col = c(NA,AUC_col, NA, NA, rep("black", length(l.num))))

  dev.off()
}
# creer un graph tracant le suivi de l'AUC de chaque acquisiton(acq) a chaque masse(mt$M.num).
# acq et mt$M.num peuvent etre egal a 1 ou plus.

monitor_plot_AUC <- function(sp, mt){
  acq = mt$acq
  if(("AUC" %in% dir("Figures"))==FALSE){
    dir.create("Figures/AUC")
  }

  if(mt$AUC_plot_exp == TRUE){
    mt$AUC_plot_exp <- "y"
  }else{
    mt$AUC_plot_exp <- ""
  }

  if(mt$AUC_each.Mass == TRUE){

   for(ma in mt$M.num){

     if(mt$AUC_each.group != FALSE){

       grp <- mt$h5[,mt$AUC_each.group][mt$acq] %>% ch()
       for(u in unique(grp)){
         ind.AUC <- mt$acq[which(u == grp)]
         mt$titre <- c("Figures/AUC/AUC_at", ma,"of",mt$AUC_each.group,u) %>% str_flatten(" ") %>% paste0(".tiff")
         plot_AUC(sp, mt, acq = ind.AUC, M.num = ma)
       }
     }else{
       mt$titre <- c("Figures/AUC/AUC_at", ma,"of",head(row.names(mt$h5)[acq])) %>% str_flatten(" ") %>% paste0(".tiff")
       plot_AUC(sp, mt, M.num = ma)
     }
   }
  }else{

    if(mt$AUC_each.group != FALSE){

      grp <- mt$h5[,mt$AUC_each.group][mt$acq] %>% ch()
      for(u in unique(grp)){
        ind.AUC <- mt$acq[which(u == grp)]
        mt$titre <- c("Figures/AUC/AUC_at", mt$M.num,"of",mt$AUC_each.group,u) %>% str_flatten(" ") %>% paste0(".tiff")
        plot_AUC(sp, mt, acq = ind.AUC, M.num = mt$M.num)
      }
    }else{
      mt$titre <- c("Figures/AUC/AUC_at", mt$M.num,"of",head(row.names(mt$h5)[acq])) %>% str_flatten(" ") %>% paste0(".tiff")
      plot_AUC(sp, mt)
    }
  }

}

export_AUC <- function(sp, mt){
  if(str_ends(mt$export.AUC, ".csv") == FALSE){
    mt$export.AUC <- paste0(mt$export.AUC, ".csv")
  }

  eph <- colnames(sp$MS.AUC) %>% as.numeric() %>% rbind(sp$MS.AUC)
  row.names(eph)[1] <- "Mass (a.u.)"

  write.table(x = eph, file = mt$export.AUC, sep = ";", append = FALSE,
              dec = ".", row.names = TRUE, col.names = FALSE)

}

# Fonctions d'analyses (AUC dynamiques) : ####
colNtoH <- function(nm_col){
  eph <- col2rgb(nm_col)
  eph <- rgb(eph[1], eph[2],eph[3], maxColorValue = 255)
  return(eph)
}

dy.mat.AUC <- function(clmat = row.names(mt$h5)[1], Mnam = mt$M.num, sp_dy = sp, Tdiv_dy = Tdiv){
  coor_MA <- match(ch(Mnam), colnames(sp_dy$MS.AUC))
  coor <- indS(clmat)

  eph <- cbind(sp_dy$xT[coor]/Tdiv_dy, sp_dy$MS.AUC[coor,coor_MA]) %>% as.matrix()
  colnames(eph) <- c("xT", paste(clmat,"mass", Mnam, sep = "_"))
  return(eph)
}

dy.AUC <- function(sp, mt, tit = titre, ACQ = mt$acq, MA = mt$M.num){

  mt$titre <- paste0(mt$wd,"/Figures/AUC/", tit,".html")

  Tbn <- c(0,0)
  Ibn <- c(0,0)
  for(i in ACQ){
    Tbn[2] <- max(Tbn[2],sp$xT[indS(row.names(mt$h5)[i])])
    for(j in MA){
      Ibn[2] <- max(Ibn[2],sp$MS.AUC[indS(row.names(mt$h5)[i]),ch(j)])
    }
  }

  unit = "s"
  Tdiv = 1
  if(Tbn[2]>600){
    unit = "min"
    Tdiv = 60
  }
  if(Tbn[2]>36000){
    unit = "h"
    Tdiv = 3600
  }

  AUC_col <- indColor(ACQ, mt$h5) %>% sapply(colNtoH) %>% rep(each = length(MA))

  dy_mat <- row.names(mt$h5)[ACQ] %>% sapply(dy.mat.AUC, Mnam = MA)
  ldy <- length(dy_mat)

  if(ldy == 1){
    dysp <- dy_mat %>% as.matrix()
  }else{
    dysp <- merge(dy_mat[[1]], dy_mat[[2]], all=TRUE)
    if(ldy >= 3){
      for(i in 3:ldy){
        dysp <- merge(dysp, dy_mat[[i]], all=TRUE)
      }
    }
  }

  rmarkdown::render(paste0(fct_dir,"/print_dyAUC.Rmd"), output_file = mt$titre)
}

monitor_plot2_AUC <- function(sp, mt){
  # Mise en forme :

  if(("AUC" %in% dir("Figures"))==FALSE){
    dir.create("Figures/AUC")
  }

  if(mt$AUC_plot_exp == TRUE){
    mt$AUC_plot_exp <- "y"
  }else{
    mt$AUC_plot_exp <- ""
  } # Exponentiel ?

  # Graphe :

  # Etape 1 : chaque mass ?
  if(mt$AUC_each.Mass == TRUE){
    # OUI

    for(ma in mt$M.num){

      # Etape 2 : chaque groupe ?
      if(mt$AUC_each.group != FALSE){
        # OUI
        grp <- mt$h5[,mt$AUC_each.group][mt$acq] %>% ch() # on repere les groupes
        for(u in unique(grp)){
          ind.AUC <- mt$acq[which(u == grp)]              # on repere les indices de chaques groupes

          # Etape 3 : plot statique ou dynamique ?
          if(mt$AUC_plot_dy == FALSE){
            # plot statique

            mt$titre <- c("Figures/AUC/AUC_at", ma,"of",u) %>% str_flatten(" ") %>% paste0(".tiff")
            plot_AUC(sp, mt, acq = ind.AUC, M.num = ma)

          }else if(mt$AUC_plot_dy == TRUE){
            # plot dynamique
            titre <- c("AUC_at", ma,"of",u) %>% str_flatten(" ")
            dy.AUC(sp, mt, tit = titre, ACQ = ind.AUC, MA = ma)
          }
        }

      }else{
        # NON

        # Etape 3 : plot dynamique ou statique ?
        if(mt$AUC_plot_dy == FALSE){
          # plot statique

          mt$titre <- c("Figures/AUC/AUC_at", ma,"of",head(row.names(mt$h5)[acq])) %>% str_flatten(" ") %>% paste0(".tiff")
          plot_AUC(sp, mt, M.num = ma)
        }else if(mt$AUC_plot_dy == TRUE){
          # plot dynamique

          titre <- c("AUC_at", ma,"of all selected") %>% str_flatten(" ")
          dy.AUC(sp, mt, tit = titre, ACQ = mt$acq, MA = ma)
        }
      }
    }
  }else{
    # NON

    # Etape 2 : chaque groupe ?
    if(mt$AUC_each.group != FALSE){
      # OUI

      grp <- mt$h5[,mt$AUC_each.group][mt$acq] %>% ch()
      for(u in unique(grp)){

        ind.AUC <- mt$acq[which(u == grp)]

        # Etape 3 : plot statique ou dynamique ?
        if(mt$AUC_plot_dy == FALSE){
          # plot statique

          mt$titre <- c("Figures/AUC/AUC_at", mt$M.num,"of",u) %>% str_flatten(" ") %>% paste0(".tiff")
          plot_AUC(sp, mt, acq = ind.AUC, M.num = mt$M.num)
        }else if(mt$AUC_plot_dy == TRUE){
          # plot dynamique

          titre <- c("AUC_at", mt$M.num,"of",u) %>% str_flatten(" ")
          dy.AUC(sp, mt, tit = titre, ACQ = ind.AUC, MA = mt$M.num)
        }
      }
    }else{
      # NON

      # Etape 3 : plot statique ou dynamique ?
      if(mt$AUC_plot_dy == FALSE){
        # plot statique

        mt$titre <- c("Figures/AUC/AUC_at", mt$M.num,"of",head(row.names(mt$h5)[acq])) %>% str_flatten(" ") %>% paste0(".tiff")
        plot_AUC(sp, mt)
      }else if(mt$AUC_plot_dy == TRUE){
        # plot dynamique

        titre <- c("AUC_at", mt$M.num,"of all selected") %>% str_flatten(" ")
        dy.AUC(sp, mt, titre, ACQ = mt$acq , MA = mt$M.num)
      }
    }
  }
}

# Fonctions d'analyses (PCA) : ####

PCA_plot_scores <- function(mt, MSpca){
  ax.pca <- c(mt$PCx, mt$PCy)

  if(length(row.names(mt$h5)[mt$acq]) > 5){
    l_ech <- length(row.names(mt$h5)[mt$acq])
    titre <- c("PCA/PC",ax.pca[1],ax.pca[2],"of",row.names(mt$h5)[mt$acq[1]],"to",row.names(mt$h5)[mt$acq[l_ech]]) %>% str_flatten("_") %>% paste0(".tiff")
    legende <- c(row.names(mt$h5)[mt$acq[1:3]], row.names(mt$h5)[mt$acq[(l_ech-2):l_ech]])
    p_ch <- c(NA, rep(16,6))
    col_l <- c(NA, alpha(mt$h5$color[mt$acq[c(1:3,((l_ech-2):l_ech))]],0.5))
  }else{
    titre <- c("PCA/PC",ax.pca[1],ax.pca[2],"of",row.names(mt$h5)[mt$acq]) %>% str_flatten("_") %>% paste0(".tiff")
    legende <- c(row.names(mt$h5)[mt$acq])
    p_ch <- c(NA, rep(16,length(mt$acq)))
    col_l <- c(NA, alpha(mt$h5$color[mt$acq],0.5))
  }

  tiff(file = titre, width = 680, height = 450,units = "px")
    par(mar = c(5,5,1,15), cex.main=2, cex.lab = 2, cex.axis = 2,mgp = c(3.5,1.5,0),xpd = FALSE)
      plot(MSpca$Tr[,ax.pca[1]], MSpca$Tr[,ax.pca[2]], pch = 16, col = alpha(MSpca$mt_col,0.5),
           xlab = paste("PC",ax.pca[1],"(",MSpca$EV[ax.pca[1]],"%)"),
           ylab = paste("PC",ax.pca[2],"(",MSpca$EV[ax.pca[2]],"%)"))
    abline(h =0, v = 0, lty = 2)

    legend("topright", bty = "n", cex = 1.5, xpd = NA, inset = c(-0.5,0),
           legend = c("Sample(s) :",legende) , pch = p_ch, col = col_l, pt.cex = 1)
  dev.off()
}
# plot des scores (mt$PCx et mt$PCy)

PCA_plot_loadings <- function(sp, mt, MSpca){
  pc <- mt$PCA_plot_loading_PC

  if(length(row.names(mt$h5)[mt$acq]) > 5){
    l_ech <- length(row.names(mt$h5)[mt$acq])
    titre <- c("PCA/Loading_PC",pc,"of",row.names(mt$h5)[mt$acq[1]],"to",row.names(mt$h5)[mt$acq[l_ech]]) %>% str_flatten("_") %>% paste0(".tiff")
    legende <- c(row.names(mt$h5)[mt$acq[1:3]], row.names(mt$h5)[mt$acq[(l_ech-2):l_ech]])
    p_ch <- c(NA, rep(16,6))
    col_l <- c(NA, alpha(mt$h5$color[mt$acq[c(1:3,((l_ech-2):l_ech))]],0.5))
  }else{
    titre <- c("PCA/Loading_PC",pc,"of",row.names(mt$h5)[mt$acq]) %>% str_flatten("_") %>% paste0(".tiff")
    legende <- c(row.names(mt$h5)[mt$acq])
    p_ch <- c(NA, rep(16,length(mt$acq)))
    col_l <- c(NA, alpha(mt$h5$color[mt$acq],0.5))
  }

  tiff(file = titre, width = 800, height = 350,units = "px")
  par(mar = c(5,5,2.5,0.2), cex.main=2, cex.lab = 2, cex.axis = 2,mgp = c(3.5,1.5,0),xpd = FALSE)

  matplot(sp$xMS_PCA, MSpca$P[,pc], type = "l",
          xlim = c(mt$pca_plot_xmin, mt$pca_plot_xmax),
          ylim = c(min(MSpca$P[,pc]), max(MSpca$P[,pc]))*1.2,
          xlab = "m/z", ylab = "Relative intensity (u.a.)",
          main = paste0("loadings of PC", pc, " (", MSpca$EV[pc], " %)"))

  bitext(MSpca$pk.ld[[pc]][1,], MSpca$pk.ld[[pc]][2,], text = MSpca$pk.ld[[pc]][1,],
         bicex = 1, off7 = 0.3)
  dev.off()
}
# plot des loadings (mt$PCA_plot_loading_PC)

# Fonctions d'analyses (autres) : ####
coor_mt <- function(acq, mta = mt$h5){
  cr <- row.names(mta)[acq] %>% sapply(indS) %>% unlist()
  if(length(acq)==1){cr <- as.integer(cr)}

  nm_eph <- "plop"
  for(i in acq){
    eph <- log(mta[i,"nbr_MS"],10) %>% floor() %>% add(1)
    nm_eph <- str_sub(paste0("000",1:mta[i,"nbr_MS"]), -eph) %>% paste(row.names(mta)[i],., sep="_") %>% c(nm_eph,.)
  }
  names(cr) <- nm_eph[-1]
  return(cr)
}
# retourne l'index de chaque spectre de chaque acquisition donnee en entree.

indS <- function(nam, mta = mt$h5){
  eph <- which(nam == rownames(mta))
  seq(mta[eph,"start"],mta[eph,"end"]) %>% return()
}
# retourne les index de tous les MS de l'acquistion.

# Fonctions utilitaires : ####

bitext <- function(x, y, text, off7, bicex = 1, colo = NULL){
  ip <- which(y >0)
  im <- which(y <0)
  if(length(ip)>0){text(x[ip], y[ip], text[ip], pos = 3, offset = off7, cex = bicex, col = colo)}
  if(length(im)>0){text(x[im], y[im], text[im], pos = 1, offset = off7, cex = bicex, col = colo)}
}
# comme text() mais repere les intensites negatives et positives,
# puis ajuste en fonction la position du texte.

ch <- function(num){as.character(num)}
# simplifie as.character() ^^ .

ctrl_color <- function(mt){
  mt_col <- mt$h5$color[mt$acq]
  eph <- which.na(mt_col)
  if(length(eph) > 0){
    mt_col[eph] <- viridis(length(eph)) %>% alpha(0.5)
  }

  eph <- which(mt_col == "")
  if(length(eph) > 0){
    mt_col[eph] <- viridis(length(eph)) %>% alpha(0.5)
  }

  return(mt_col)
}
# control de la colour

det_c <- function(brn,vec){
  subtract(vec,brn) %>% sapply(abs) %>% which.min()
}
# retourne l'index valable d'une borne le long d'un vecteur.

dizaine <- function(x){
  eph <- log(x,10) %>% floor() %>% multiply_by(10)
  divide_by(x,eph) %>% floor() %>% multiply_by(eph)
}
# arrondi a la dizaine inferieure (utilise dans rapport_pca.Rmd).

export_for_fit <- function(x_range = c(136.5,137.5), mat = sp$MS_m, out.num = mt$acq){
  br_m <- det_c(x_range[1], sp$xMS)
  br_M <- det_c(x_range[2], sp$xMS)

  for(i in out.num) {
    eph <- paste(rownames(mat)[i], "to", x_range[1], "at", x_range[2], sep = "_") %>%
      str_replace_all("\\.","u") %>% paste0("_fityk.txt")
    cbind(sp$xMS[br_m:br_M], mat[br_m:br_M,out.num]) %>%
      write.table(file = eph, row.names = FALSE, col.names = FALSE)
  }
}
# export un courbe pour fityk

heure <- function(){
  eph <- Sys.time()
  eph <- str_split(eph,pattern = " ")
  return(eph[[1]][2])
}
# donne l'heure

indColor <- function(acq, mt = mt_h5){
  if(length(mt$color) >0){
    ind_col <- mt$color[acq]
    ind_col[match("", ind_col)] <- "black"
  }else{
    ind_col <- rainbow(length(acq))
  }
  return(ind_col)
}
# retourne les couleurs de chaque acquistion.

max_mat <- function(meta, MS){
  apply(MS[meta[2]:meta[3],],2,max)
}
# retourne le spectre max de tous les spectres d'une acquisition.

mean_mat <- function(meta, MS){
  colMeans(MS[meta[2]:meta[3],])
}
# retourne le spectre moyen de tous les spectres d'une acquisition.

meta_col <- function(acq, mt = mt_h5){
  eph <- indColor(acq, mt) %>% which.na(.)
  mt$color[acq[eph]] <- rainbow(length(eph))
  return(mt)
}
# ajoute des couleurs dans le fichier meta mt_h5 pour chaque acquistion
# ou la couleur n'est pas precisee.

nS <- function(Poisson_Rouge = "?"){row.names(mt$h5)}
# une fonction rapide pour se souvenirs du numero des acquisitons.

print.h <- function(txt = "hello there"){
  heure() %>% paste0(txt,", ",.) %>% print()
}
# print le texte suivi de l'heure

rep_mtm <- function(col.nam, mt){
  sapply(mt$acq, rep_mtu, col.nam = col.nam, mt = mt) %>% unlist()
}
# comme rep_mtm mais pour plusisieurs acquistions.

rep_mtu <- function(acq, col.nam, mt = mt){
  eph <- which(col.nam == colnames(mt$h5))
  rep(mt$h5[acq,eph], mt$h5$nbr_MS[acq])
}
# repete l'element de la colonne "col.nam" autant de fois que
# le nombre de spectre d'une acquisition.

sp_max <- function(ind = sp$indX, data){
  return(max(data[,ind]))
}
# retourne l'intensite maximale pour les index de "ind"

w.equal <- function(vec,nb){which(vec == nb)}
# retourne la position des elements de vec egaux a nb.

which.na <- function(x){which(is.na(x) == TRUE)}
# retourne la position des NA sur un vecteur.

which.not.na <- function(x){which(is.na(x) == FALSE)}
# retourne la position des non-NA sur un vecteur.

which.sup <- function(vec, threshold){
  return(which(vec > threshold))
}
# retourne la position des elements de vec superieur au seuil.

# Function inutile : ####
citation_list <- {list(
c("Il faut aller trop loin pour decouvrir les limites", "Joris Huguenin"),
c("Dieu, aie pitie de nous, nous sommes a la merci des ingenieurs !", 'Dr.Malcom, Jurassic Park'),
c("Grab a brush and put on a little make-up","System of a Down"),
c("The Sun Machine is Coming Down, and We're Gonna Have a Party", "David Bowie"),
c("I'm just a poor boy, I need no sympathy.", "Queen"),
c("Au village, sans pretention, J'ai mauvaise reputation","Georges Brassens"),
c("Debout les gars, reveillez-vous ! On va au bout du monde","Huges Aufray"),
c("A m'asseoir sur un banc cinq minutes avec toi, Et regarder les gens tant qu'y en a.","Renaud"),
c("Ready or not, here I come, you can't hide. Gonna find you and make you want me","The Fugees"),
c("Emancipate yourselves from mental slavery.","Bob Marley"),
c("Parce que c'est notre BROCHEEEET !!!.", "Manuel Macro"),
c("Hey DJ met nous donc du Funk, que je danse le MIA. Je danse le MIA.", "IAM"),
c("Doo, doo, doo, doo, doo, doo, doo, doo", "Lou Reed"),
c("L'obscurite ne peut pas chasser l'obscurite, seule la lumiere le peut. La haine ne peut pas chasser la haine, seul l'amour le peut.", "Martin Luther King"),
c("La vie, ce n'est pas d'attendre que les orages passent, c'est d'apprendre a danser sous la pluie.", "Seneque"),
c("Nos vies sont pleines de catastrophes qui n'ont jamais eu lieu.", "Auteur inconnu"),
c("S'il y a un probleme, il y a une solution. S'il n'y a pas de solution, alors ce n'est pas un probleme.", "Auteur inconnu"),
c("Si vous pouvez le rever, vous pouvez le faire.", "Walt Disney"),
c("Ils ne savaient pas que c'etait impossible, alors ils l'ont fait.", "Mark Twain"),
c("J'ai decide d'etre heureux parce que c'est bon pour la sante.", "Voltaire"),
c("Si vous pensez que l'aventure est dangereuse, essayez la routine, elle est mortelle.", "Paulo Coelho"),
c("Les gens les plus heureux n'ont pas tout ce qu'il y a de mieux. Ils font juste de leur mieux avec tout ce qu'ils ont.", "Auteur inconnu"),
c("Le veritable voyage ne consiste pas a chercher de nouveaux paysages, mais a avoir de nouveaux yeux.", "Marcel Proust"),
c("Avec trop on se perd. Avec moins on se trouve.", "Tchouang Tseu"),
c("N'aie pas peur d'avancer lentement. Aie peur de rester immobile.", "Proverbe chinois"),
c("Ne cherche pas le bonheur, cree-le.", "Auteur inconnu"),
c("Ne t'inquiete pas de l'echec. Inquiete-toi de ce que tu manques si tu n'essayes meme pas.", "Jack Canfield"),
c("Mieux vaut fait que parfait.", "Auteur inconnu"),
c("Lorsqu'on regarde dans la bonne direction, il ne reste plus qu'a avancer.", "Proverbe bouddhiste"),
c("Un objectif bien defini est a moitie atteint.", "Abraham Lincoln"),
c("Quand on ose, on se trompe souvent. Quand on n'ose pas, on se trompe toujours.", "Romain Rolland"),
c("La vie c'est comme une bicyclette, il faut avancer pour ne pas perdre l'equilibre.", "Albert Einstein"),
c("Il y a deux facons de penser. L'une est de croire que les miracles n'existent pas. L'autre est de croire que chaque chose est un miracle.", "Albert Einstein"),
c("Fais de ta vie un reve et d'un reve une realite.", "Antoine de St Exupery"),
c("Il y a plus de courage que de talent dans la plupart des reussites.", "Felix Leclerc"),
c("Ce que nous sommes est le resultat de ce que nous avons pense.", "Bouddha"),
c("Les gagnants cherchent des moyens, les perdants des excuses.", "Franklin Roosevelt"),
c("Un voyage de mille lieues commence toujours par un premier pas.", "Lao Tseu"),
c("Tous les jours a tous points de vue, je vais de mieux en mieux.", "Emile Coue"),
c("Il faut toujours viser la lune car meme en cas d'echec on atterrit dans les etoiles.", "Oscar Wilde"),
c("Ce n'est pas parce que les choses sont difficiles que nous n'osons pas les faire, c'est parce que nous n'osons pas les faire qu'elles sont difficiles.", "Seneque"),
c("N'attendez pas d'etre heureux pour sourire. Souriez plutot afin d'etre heureux.", "Edward L. Kramer"),
c("Si tu fais ce que tu as toujours fait, tu obtiendras ce que tu as toujours obtenu.", "Tony Robbins"))}

# ####
#
# ####
# Function de lancement des analysis : ####
PTR_MS_analysis <- function(sp, mt, citation_list){
  # Control Color ####
  if(mt$ctrl_col == TRUE){mt$h5$color[mt$acq] <- ctrl_color(mt)}
  # Control Plot ####
  if(mt$ctrl_plot == TRUE){ctrl_plot_fct(ls_h5, sp, mt)}

  # Control Alignement ####
  if(mt$ctrl_align == TRUE){
    if(("Alignement" %in% dir("Figures"))==FALSE){
      dir.create("Figures/Alignement")
    }
    align_print_ctr(mat = sp$MS_n_al, sp, suffixe = "not_align")
    align_print_ctr(mat = sp$MS, sp,  suffixe = "align")
  }

  # Realignement ####
  if(12 == 42){re.align(sp)}

  # Gestion du temps ####
  if(mt$T_para_initiaux == TRUE){
    sp <- re.init.T.para(sp)
  }

  if(mt$T_recalc_para == TRUE){
    sp <- re.calc.T.para(sp, mt)
  }

  # View Plot ####
  if((mt$view_plot[1] == TRUE)&(is.logical(mt$view_plot) == TRUE)){view_plot_fct(ls_h5, sp, mt)}

  if(is.numeric(mt$view_plot) == TRUE){
    sapply(mt$view_plot, view_one_plot_fct, ls_h5, sp, mt)
  }

  if(mt$view_dy_plot == TRUE){
    mt$review_dy_select <- dy.spectra(sp, mode = mt$view_dy_mode, select.sp = mt$view_dy_select)
  }

  if(mt$view_circular == TRUE){
    if(("Circular" %in% dir("Figures"))==FALSE){
      dir.create("Figures/Circular")
    }

    options(knitr.table.format = "html")
    sapply(mt$view_circular_num_acq, print_spectre_circle,
           sp = sp, mode = mt$view_circular_mode, thr = mt$view_circular_threshold)
    library(knitr)
  }

  # Soustraction de blanc ####
  if(mt$sub_blank == TRUE){
    sp <- blank_subtraction(sp, mt)
  }

  if(mt$sub_blank_reverse == TRUE){
    sp <- reverse_blank_subtraction(sp)
  }

  # Aire sous la courbe ####
  if(mt$AUC == TRUE){
    if(length(sp$MS.AUC) == 0){
      print("Calcul des AUC en cours")
      length(citation_list) %>% sample(1) %>% citation_list[[.]] %>% print()
      sp$MS.AUC <- peak_AUC(sp)
      print('Fin du calcul')
      print("Sauvegarde de l'environnement (.RData)")
      str_split(mt$wd,"/")[[1]] %>% length() %>% str_split(mt$wd,"/")[[1]][.] %>% paste0("_AUC.RData") %>% save.image()

    }
    monitor_plot2_AUC(sp, mt)

    if(mt$export.AUC != FALSE){
      export_AUC(sp, mt)
    }
  }

  # Calibration & validation (a finir) ####
  # Calibration

  if(mt$calibration == TRUE){
    if(length(sp$MS.AUC) == 0){
      print("Calcul des AUC en cours")
      length(citation_list) %>% sample(1) %>% citation_list[[.]] %>% print()
      sp$MS.AUC <- peak_AUC(sp)
      print('Fin du calcul')
      print("Sauvegarde de l'environnement (.RData)")
      str_split(wd,"/")[[1]] %>% length() %>% str_split(wd,"/")[[1]][.] %>% paste0("_AUC.RData") %>% save.image()
    }

    if(("Calibration" %in% dir())==FALSE){
      dir.create("Calibration")
    }

    if(mt$cal_plot_exp == TRUE){
      cal_exp <- "xy"
    }else{
      cal_exp <- ""
    }

    if((mt$M.conc == FALSE)||(length(mt$M.conc) == 1)){
      mt$M.conc <- length(mt$M.num) %>% rep(mt$M.conc,.)
    }else if(length(mt$M.conc) != length(mt$M.num)){
      mt$M.conc <- length(mt$M.num) %>% rep(1,.)
      print("WARNING ! Lengths of mt$M.conc and mt$M.nim are not equal. All M.conc are 1.")
    }

    coor <- coor_mt(mt$acq)
    C_cal <- rep_mtm("concentration", mt)
    Col_cal <- rep_mtm("color", mt)
    reg.coef <- matrix(NA, nrow = length(mt$M.num), 8,
                       dimnames = list(ch(mt$M.num),c("slope", "intercept", "std.error.slope", "std.error.inter",
                                                      "interval.confidence.slope.0.025","interval.confidence.slope.0.975",
                                                      "interval.confidence.intercept.0.025","interval.confidence.intercept.0.975")))
    p <- 0

    if(mt$cal_weight != FALSE){
      wgt <- rep_mtm(mt$cal_weight, mt)
    }else{
      wgt <- length(coor) %>% rep(1,.)
    }

    for(m in mt$M.num){
      p <- p + 1
      cal_conc <- data.frame(concentration = C_cal*mt$M.conc[p], AUC = sp$MS.AUC[coor,ch(m)], w4lm = wgt)
      cal_C <- lm(concentration~AUC, data=cal_conc, weights = w4lm)

      if(mt$cal_rapport == TRUE){
        titrepdf <- paste0("calibration_rapport_at_",m,"u")
        titrepdf <- paste0(wd,"/Calibration") %>% dir() %>% grep(titrepdf,.) %>% length() %>% add(1) %>%
          paste0(wd,"/Calibration/",titrepdf,"_",.,".pdf")
        rmarkdown::render(paste0(fct_dir,"/rapport_calibration.Rmd"), output_file = titrepdf)
      }

      reg.coef[p,c(2,1)] <- cal_C$coefficients

      eph <- summary(cal_C)
      reg.coef[p,c(4,3)] <- eph$coefficients[,2]

      reg.coef[p,c(7,5,8,6)] <- confint(cal_C)

      main.text <- str_c("Concentration = AUC * ", format(reg.coef[p,1], scientific=TRUE, digits = 3),
                         " + ", format(reg.coef[p,2], scientific=TRUE, digits = 3))

      x_data <- cal_conc$AUC
      y_data <- cal_conc$concentration

      if(mt$cal_plot_exp == TRUE){
        title.tiff <- paste0("Calibration/calibration_",m,"_lin_reg_exp.tiff")
        eph <- unique(c(w.equal(x_data,0), w.equal(y_data,0)))
        x_data[eph] <- NA
        y_data[eph] <- NA
      } else{
        title.tiff <- paste0("Calibration/calibration_",m,"_lin_reg.tiff")
      }

      tiff(file = title.tiff, width = 500, height = 450,units = "px")
      par(mar = c(5,5,2,1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0), xpd = FALSE)
      matplot(x_data, y_data, pch = 3, col = "white",
              xlab = "AUC", ylab = "concentration", log = cal_exp)
      points(cal_conc$AUC, cal_conc$concentration, pch = 3, col = Col_cal)
      abline(a = reg.coef[p,2], b = reg.coef[p,1], col = "red", lwd = 2, untf = TRUE)
      abline(a = reg.coef[p,7], b = reg.coef[p,5], col = "chocolate3", lwd = 2, untf = TRUE, lty= 2)
      abline(a = reg.coef[p,8], b = reg.coef[p,6], col = "chocolate3", lwd = 2, untf = TRUE, lty= 2)
      mtext(side = 3, cex = 1.5,line = 0.1,
            text = main.text)
      legend("bottomright", legend = paste("mass at",m,"(m/z)"), bty = "n", cex = 1.5)
      dev.off()
    }

    tiff(file = "Calibration/cal_slope.tiff", width = 500, height = 450,units = "px")
    par(mar = c(5,5,2,1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0), xpd = FALSE)
    matplot(mt$M.num, reg.coef[,1], pch = 16, ylim = c(0,0.02),
            xlab = "mass (m/z)", ylab = "slope")
    abline(h = 0, lty = 2)
    dev.off()

    tiff(file = "Calibration/cal_intercept.tiff", width = 500, height = 450,units = "px")
    par(mar = c(5,5,2,1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0), xpd = FALSE)
    matplot(mt$M.num, reg.coef[,2], pch = 16,
            xlab = "mass (m/z)", ylab = "y-intercept")
    abline(h = 0, lty = 2)
    dev.off()

    eph <- data.frame(Masse = mt$M.num, reg.coef)
    write.table(eph, file = "Calibration/cal_regression_coefficient.csv", row.names = FALSE, col.names = TRUE, sep = ";")
    write.table(cbind(mt$M.num, reg.coef[,1]), file = "Calibration/cal_slope_fityk.csv", row.names = FALSE, col.names = FALSE)
    write.table(cbind(mt$M.num, reg.coef[,1]), file = "Calibration/cal_intercept_fityk.csv", row.names = FALSE, col.names = FALSE)

    save(cal_C, file = "calibration.RData")
  }

  # Validation
  if(mt$validation == TRUE){
    coor <- coor_mt(mt$acq)
    val_conc <- data.frame(AUC = MS.AUC[coor,ch(mt$M.num)])
    predict(cal_C, newdata = val_conc)
    predict(cal_C, newdata = val_conc, interval="prediction")
    predict(cal_C, newdata = val_conc, interval="confidence")
  }

  # PCA ####
  if(mt$PCA == TRUE){

    if(("PCA" %in% dir())==FALSE){
      dir.create("PCA")
    }

    if(mt$PCA_AUC == TRUE){
      pca_mat <- sp$MS.AUC[coor_mt(mt$acq),] %>% as.matrix()
      eph <- apply(pca_mat,2, sd)
      sd_nn <- which(eph > 0)

      pca_mat <- pca_mat[,sd_nn]
      sp$xMS_PCA <- colnames(pca_mat) %>% as.numeric()
    }else{
      pca_mat <- sp$MS[coor_mt(mt$acq),]
      sp$xMS_PCA <- sp$xMS
    }

    MSpca <- pca(pca_mat, ncomp = mt$npc)

    MSpca$EV <- MSpca$explvarx[,3] %>% multiply_by(100) %>% round(2)
    MSpca$mt_col <- rep_mtm("color", mt)

    if(mt$PCA_AUC == TRUE){
    }else{
      options(warn=-1)
      MSpca$pk.ld <- t(MSpca$P) %>% pk.list3(sp) %>% lapply(pk.short)

      options(warn=0)
    }

    if(mt$PCA_report == TRUE){
      titrepdf <-  dir("PCA") %>% grep("rapport_pca",.) %>% length() %>%
        add(1) %>% paste0(wd,"/PCA/rapport_pca",.,".pdf")
      rmarkdown::render(paste0(fct_dir,"/rapport_pca.Rmd"), output_file = titrepdf)
    }

    if(mt$PCA_plot_scores == TRUE){

      npc_plot <- seq(1:mt$npc)

      for(i in 1:(mt$npc-1)){
        mt$PCx <- i
        for(j in (i+1):mt$npc){
          mt$PCy <- j
          PCA_plot_scores(mt, MSpca)
        }
      }
    }

    if(mt$PCA_plot_loading == TRUE){
      if(mt$PCA_plot_loading_PC[1] == "all"){
        mt$PCA_plot_loading_PC <- seq(1:mt$npc)
      }

      npc_plot <- mt$PCA_plot_loading_PC

      for(i in npc_plot){
        mt$PCA_plot_loading_PC <- i
        PCA_plot_loadings(sp, mt, MSpca)
      }
      mt$PCA_plot_loading_PC <- npc_plot
    }
  }
  # ICA ####
  if(mt$ICA == TRUE){

    if(("ICA" %in% dir())==FALSE){
      dir.create("ICA")
    }

    if((exists("nc_old") == FALSE)||((mt$nc != nc_old)|(mt$acq != acq_old))){
      MSica <- icajade(sp$MS.AUC[coor_mt(mt$acq),], nc = mt$nc)
      # MSica <- icajade(sp$MS[coor_mt(mt$acq),], nc = mt$nc)
      nc_old <- mt$nc
      acq_old <- mt$acq
    }

    dim(MSica$S) # 88 ; 8
    dim(MSica$M) # 6463 ; 8
    dim(MSica$Y) # 88 ; 8
    dim(MSica$Q) # 8 : 6463
    dim(MSica$W) # 8 ; 6463
    dim(MSica$R) # 8 ; 8
    length(MSica$vafs) # 8
    length(MSica$iter) # 1

    xAUC <- colnames(sp$MS.AUC) %>% as.numeric()
    layout(matrix(1:6,2,3, byrow = TRUE), TRUE)
    for(i in 1:6){
      hist(MSica$M[,i], breaks = 1000, xlim = c(-800,800))
      # plot(xAUC, MSica$M[,i], type = "l")
    }

    dev.off()
    icaplot(MSica)

    x <- pmin(3, pmax(-3, stats::rnorm(50)))
    y <- pmin(3, pmax(-3, stats::rnorm(50)))
    xhist <- hist(x, breaks = seq(-3,3,0.5), plot = FALSE)
    yhist <- hist(y, breaks = seq(-3,3,0.5), plot = FALSE)
    top <- max(c(xhist$counts, yhist$counts))
    xrange <- c(-3, 3)
    yrange <- c(-3, 3)
    layout.show(nf)

    par(mar = c(3,3,1,1))
    plot(x, y, xlim = xrange, ylim = yrange, xlab = "", ylab = "")
    par(mar = c(0,3,1,1))
    barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0)
    par(mar = c(3,0,1,1))
    barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0, horiz = TRUE)

    par(def.par)  #- reset to default
  }

  # end of analysis ####
}
# ####
#
# ####

########################
##### END OF CODE ######
########################

# Archive de fonctions

# reduction ####
#reduction <- function(ls_h5){
#   mat_h5 <- red_ls(ls_h5)
#   xMS_r <- colnames(mat_h5) %>% as.numeric()
#   seuil <- apply(mat_h5, 2, sd) %>% median() %>% multiply_by(20)
#   # Le seuil est defini en prenant 20 fois l'ecart-type du bruit de fond
#   x_coef <- which(mat_h5[1,]>=seuil)
#
#   ###
#   # nbr_pixel <- matrix(c(1:10,seq(15,100,by = 5)),2,28,byrow = TRUE)
#   #
#   # for(i in 1:28){
#   #   seuil <- apply(mat_h5, 2, sd) %>% median() %>% multiply_by(nbr_pixel[1,i])
#   #   nbr_pixel[2,i] <- which(mat_h5[1,]>=seuil) %>% length()
#   # }
#   #
#   # matplot(nbr_pixel[1,],nbr_pixel[2,], pch = 16,
#   #         xlab = "Ratio signal sur bruit", ylab = "Taille d'un spectre")
#   ###
#   for(i in 2:nrow(mat_h5)){
#     x_coef <- union(x_coef, which(mat_h5[i,]>seuil)) %>% sort()
#   }
#
#   x_zero <- c(1,x_coef[1]-1)
#   for(i in 2:(length(x_coef)-1)){
#     if( (x_coef[i] - x_coef[i-1]) != 1){
#       x_zero <- union(x_zero, x_coef[i]-1)
#     }
#     if( (x_coef[i] - x_coef[i+1]) != -1){
#       x_zero <- union(x_zero, x_coef[i]+1)
#     }
#   }
#   x_coef <- union(x_zero, x_coef) %>% sort()
#
#   MS <- matrix(NA,nrow(mat_h5),length(x_coef),dimnames = list(rownames(mat_h5),xMS_r[x_coef]))
#   eph <- mat_h5
#   eph[,x_zero] <- runif(length(x_zero)*nrow(eph)) %>% round(2) %>% matrix(ncol=length(x_zero))
#   MS <- eph[,x_coef]
#   xMS <- xMS_r[x_coef]
#   return(list(MS,xMS,seuil))
# }
####
#
# find_peak <- function(mass_mol, x_ax = xMS, y_mat = all_MS){
#   mass_mol <- add(mass_mol,1.00794) %>% round(4)
#   x_m <- mass_mol - 0.35
#   x_M <- mass_mol + 0.35
#   xbr_m <- det_c(x_m,x_ax)
#   xbr_M <- det_c(x_M,x_ax)
#   ind_calib <- max.col(y_mat[,xbr_m:xbr_M]) %>% add(xbr_m) %>% subtract(1)
#   return(ind_calib)
# }
# # trouve le maximum local pour une masse donnee (non protonnee)
#
# fit_xMS <- function(x_ind, x_vec = xMS_r, Mcal = mass_cal){
#   x <- x_vec[x_ind]
#   y <- add(Mcal,1.0079)
#   fit <- lm(y~x)
#   vec <- xMS_r*fit$coefficients[2] + fit$coefficients[1]
#   return(vec)
# }
# # fit une abscisse
#####

# alignement #####
# align_fct <- function(sp){
# for(i in 1:length(sp$indX)){
#   if(length(sp$indX[[i]])>1){
#     p.max <- apply(sp$MS[,sp$indX[[i]]],1,which.max)
#     shift.med <- median(p.max)
#     p.shift <- subtract(p.max,shift.med)
#     shift.wdth <- range(p.shift) %>% diff()
#     shift.start <- subtract(p.shift, max(p.shift)) %>% abs()
#     shift.end <- subtract(p.shift, min(p.shift))
#
#     MS.ctr <- matrix(NA, nrow = nrow(sp$MS), ncol = add(length(sp$indX[[i]]), shift.wdth))
#     for(j in 1:nrow(sp$MS)){
#       a <- runif(shift.start[j],0,1) %>% round(2)
#       b <- runif(shift.end[j],0,1) %>% round(2)
#       MS.ctr[j,] <- c(a, sp$MS[j, sp$indX[[i]]], b)
#     }
#
#     a <- max(p.shift) %>% add(1)
#     b <- length(sp$indX[[i]]) %>% add(a) %>% subtract(1)
#     sp$MS[,sp$indX[[i]]] <- MS.ctr[,a:b]
#   }
# }
# return(sp)
# }
# aligne les pics
#####

# detection de pics ####
#
# peak.max <- function(mat, sp = sp){
#   MS.max <- matrix(NA, nrow = nrow(mat), ncol = length(sp$indX),
#                    dimnames = list(rownames(mat),names(sp$indX)))
#   MS.max.ind <- MS.max
#   for (i in 1:nrow(mat)){
#     for (j in 1:length(sp$indX)){
#       MS.max[i,j] <- max(mat[i,sp$indX[[j]]])
#       MS.max.ind[i,j] <- which.max(mat[i,sp$indX[[j]]]) %>% sp$indX[[j]][.]
#     }
#   }
#   return(list("MS.max" = MS.max, "MS.max.ind" = MS.max.ind))
# }
## retourne une matrice mesurant le maximum de chaque masse detectee.
#
# pk.list <- function(mat, sp = sp, sh.l = FALSE, w.sub = 4){
#
#   mat_raw <- mat
#   if(min(mat)<0){
#     mat[which(mat <0)] <- 0
#   }
#
#   list.max <- peak.max(mat, sp)
#   mat.max <- list.max$MS.max
#
#   peak_mat <- apply(mat.max, 1, pk.mat, indX = sp$indX, sh.l = sh.l)
#   if(is.list(peak_mat)==FALSE){
#     eph <- nrow(peak_mat)/4
#     peak_mat <- matrix(peak_mat, 4, eph, byrow = FALSE, dimnames = list(c("W", "int", "pos", "coor"),ch(1:eph)))
#     colnames(peak_mat) <- ch(peak_mat[1,])
#     peak_mat <- list(peak_mat)
#   }
#   peak_mat <- lapply(peak_mat,pk.min.sub, win.sub = w.sub)
#
#   for(i in 1:length(peak_mat)){
#     pk.min.sub(peak_mat[[i]], win.sub = w.sub)
#   }
#
#   for(i in 1:length(peak_mat)){
#     if(length(peak_mat[[i]]) > 4){
#       eph <- sapply(peak_mat[[i]][1,], match, table = colnames(list.max$MS.max.ind))
#       peak_mat[[i]][4,] <- list.max$MS.max.ind[i,eph]
#     }else if(length(peak_mat[[i]]) == 4){
#       eph <- sapply(peak_mat[[i]][1], match, table = colnames(list.max$MS.max.ind))
#       peak_mat[[i]][4] <- list.max$MS.max.ind[i,eph]
#     }
#   }
#
#   if(min(mat_raw)<0){
#     peak_mat_p <- peak_mat
#     mat <- -mat_raw
#     mat[which(mat <0)] <- 0
#
#     list.max <- peak.max(mat, sp)
#     mat.max <- list.max$MS.max
#
#     peak_mat <- apply(mat.max, 1, pk.mat, indX = sp$indX, sh.l = sh.l)
#     peak_mat <- lapply(peak_mat,pk.min.sub, win.sub = w.sub)
#
#     for(i in 1:length(peak_mat)){
#       if(length(peak_mat[[i]]) > 4){
#         eph <- sapply(peak_mat[[i]][1,], match, table = colnames(list.max$MS.max.ind))
#         peak_mat[[i]][4,] <- list.max$MS.max.ind[i,eph]
#         peak_mat[[i]][2,] <- -peak_mat[[i]][2,]
#       }else if(length(peak_mat[[i]]) == 4){
#         eph <- sapply(peak_mat[[i]][1], match, table = colnames(list.max$MS.max.ind))
#         peak_mat[[i]][4] <- list.max$MS.max.ind[i,eph]
#         peak_mat[[i]][2] <- peak_mat[[i]][2]
#       }
#     }
#     peak_mat_n <- peak_mat
#
#     for(i in 1:length(peak_mat)){
#       pkmat <- cbind(peak_mat_p[[i]],peak_mat_n[[i]])
#       n.pk <- ncol(pkmat)
#       doublon <- 0
#       for(j in 1:(n.pk-1)){
#         dble <- match(pkmat[1,j],pkmat[1,(j+1):n.pk])
#         if(is.na(dble)==FALSE){
#           if(which.max(c(abs(pkmat[2,j]),abs(pkmat[2,j+dble])))==2){
#             doublon <- c(doublon,j)
#           } else{
#             doublon <- c(doublon,j+dble)
#           }
#         }
#       }
#
#       if(length(doublon)>1){
#         peak_mat[[i]] <- pkmat[,-doublon[-1]]
#       } else{
#         peak_mat[[i]] <- pkmat
#       }
#
#     }
#   }
#
#   return(peak_mat)
# }
# retourne une liste de peak pour chaque spectre de mat (adapte a MS, MS_m, loading,...)
#
# pk.mat <- function(matt, indX, sh.l = FALSE){
#   val <- boxplot.stats(matt)[[4]]
#   if(sh.l == TRUE){
#     val <- boxplot.stats(val)[[4]]
#   }
#   if(length(val)<2){
#     eph <- which.max(matt)
#     val <- matt[c(eph,which.max(matt[-eph]))]
#   }
#   xpk <- names(matt) %>% as.numeric()
#   idx <- which(matt %in% c(val))
#
#   cord <- lapply(indX[ch(xpk[idx])],median) %>% unlist() %>% round()
#
#   return(rbind("W"= xpk[idx],"int"= val, "pos"= idx, "coor" = cord))
# }
# # retourne une matrice des peak pour le spectre sp.m (utilise dans pk.list)
#
# pk.min.sub <- function(mats = peak_mat[[45]], win.sub = w.sub){
#   pk.min <- which.min(mats[2,])
#   for (j in 1:(ncol(mats)-1)){
#     eph <- which(abs(mats[1,] - mats[1,j]) <= win.sub)
#     if(length(eph)>1){
#       pk.min <- c(pk.min,eph[which(mats[2,eph] < mats[2,j], arr.ind = TRUE)])
#     }
#   }
#   pk.min <- unique(pk.min) %>% sort()
#
#   mats <- mats[,-pk.min]
#   if(length(which(mats[2,] < 0.05*max(mats[2,])))>0){
#     mats <- mats[,-which(mats[2,] < 0.05*max(mats[2,]))]
#   }
#   return(mats)
# }
## retourne une matrix sans les petits peaks (utilisee dans pk.list)
#####


