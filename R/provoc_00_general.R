#### Provoc_00_general ####

#' provoc : Perform a Rapid Overview for the Volatils Organic Compounds
#'
#' analyze data of VOC by PTR-ToF-MS Vocus
#'
#' @docType package
#' @name provoc
#'
#' @import dygraphs
#' @importFrom magrittr add
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr subtract
#' @importFrom magrittr %>%
#' @importFrom MALDIquant alignSpectra
#' @importFrom MALDIquant averageMassSpectra
#' @importFrom MALDIquant binPeaks
#' @importFrom MALDIquant createMassSpectrum
#' @importFrom MALDIquant detectPeaks
#' @importFrom MALDIquant filterPeaks
#' @importFrom MALDIquant intensityMatrix
#' @importFrom MALDIquant smoothIntensity
#' @importFrom rhdf5 H5Fclose
#' @importFrom rhdf5 H5Fopen
#' @importFrom rmarkdown render
#' @importFrom scales alpha
#' @importFrom stats density
#' @importFrom stats median
#' @importFrom stringr str_flatten
#' @importFrom stringr str_pad
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_split
#' @importFrom stringr str_sub
#' @importFrom usethis use_dev_package
#' @importFrom usethis use_package
#' @importFrom viridis viridis
#' @importFrom xts xts
NULL

usethis::use_package("dygraphs", min_version = "1.1.0.0")
usethis::use_package("magrittr", min_version = "2.0.0")
usethis::use_package("MALDIquant", min_version = "1.19.0")
usethis::use_dev_package("rhdf5", remote = "grimbough/rhdf5")
usethis::use_package("rmarkdown", min_version = "2.11")
usethis::use_package("scales", min_version = "1.1.0")
usethis::use_package("stringr", min_version = "1.4.0")
usethis::use_package("usethis", min_version = "2.0.0")
usethis::use_package("viridis", min_version = "0.6.0")
usethis::use_package("xts", min_version = "0.12.0")


#### Fonction inutile ####
citation.list <- {list(
  c("Il faut aller trop loin pour decouvrir les limites.", "Joris Huguenin"),
  c("Les trous dans les pantalons, c'est comme les enfants, ca n'arrete pas de grandir", "Joris Huguenin"),
  c("Dieu, aie pitie de nous, nous sommes a la merci des ingenieurs !", 'Dr.Malcom, Jurassic Park'),
  c("Grab a brush and put on a little make-up.","System of a Down"),
  c("The Sun Machine is Coming Down, and We're Gonna Have a Party.", "David Bowie"),
  c("I'm just a poor boy, I need no sympathy.", "Queen"),
  c("Au village, sans pretention, J'ai mauvaise reputation.","Georges Brassens"),
  c("Debout les gars, reveillez-vous ! On va au bout du monde.","Huges Aufray"),
  c("Tu dis qu'si les elections ca changeait vraiment la vie \n
     Y'a un bout d'temps, mon colon, qu'voter ca s'rait interdit","Renaud"),
  c("Ready or not, here I come, you can't hide. Gonna find you and make you want me.","The Fugees"),
  c("Emancipate yourselves from mental slavery.","Bob Marley"),
  c("Hey DJ met nous donc du Funk, que je danse le MIA. Je danse le MIA.", "IAM"),
  c("Doo, doo, doo, doo, doo, doo, doo, doo.", "Lou Reed"),
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
  c("Dieu existe-elle ?", "Patrick Sebastien"),
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
  c("Si tu fais ce que tu as toujours fait, tu obtiendras ce que tu as toujours obtenu.", "Tony Robbins"),
  c("Redemarrage de l'evaluation d'une promesse interrompue.", "R poetic warning message"),
  c("Je suis gentil avec tout le monde, celui qui dit le contraire je lui foutrai mon poing dans la gueule.", "Leo Ferre"),
  c("Le desespoir est une forme superieure de la critique.", "Leo Ferre"),
  c("Les diplomes sont faits pour les gens qui n'ont pas de talent.","Pierre Desproges"),
  c("Bal tragique a Colombey, un mort.","Hara Kiri"),
  c("Si la matiere grise etait plus rose, le monde aurait moins les idees noires.","Pierre Dac"),
  c("J'ai pris la decision de ne plus etre influencable. Qu'est-ce que vous en pensez ?","Patrick Sebastien"),
  c("Est-il indispensable d'etre cultive quand il suffit de fermer sa gueule pour briller en societe ?","Pierre Desproges"),
  c("On ne discute pas recettes de cuisine avec des anthropophages.", "Jean-Pierre Vernant"))}

#### Provoc_01_gest_meta ####
#### Gestion of time ####

# acq.time
#' Title
#'
#' @param ls.t
#'
#' @return
#'
#' @examples
acq.time <- function(ls.t = ls_h5[[1]]){
  oldw <- getOption("warn")
  options(warn = -1)

  eph <- ls.t$AcquisitionLog$Log$timestring[1] %>%
    as.POSIXct(format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")

  options(warn = oldw)
  # return(T_acq)
  return(eph)
}

# reinitialise les parametres date et time
#' Title
#'
#' @param L
#'
#' @return
#' @export
#'
#' @examples
re.init.T.para <- function(L = sp){
  L$Trecalc$date <- L$Tinit$date
  L$Trecalc$timing <- L$Tinit$timing
  L <- wf.update("re.init.T.para","sp", L)
  return(L)
}

# calcul les nouveaux date et temps
#' Title
#'
#' @param L
#'
#' @return
#' @export
#'
#' @examples
re.calc.T.para <- function(L = sp){
  vec_T0 <- L$mt$meta[,"acq_T0 (ID)"]
  vec_D <- L$mt$meta[,"delta_T (s)"]

  fmr <- c(which.na(vec_T0),which(vec_T0 == ""))
  if(length(fmr) > 0) vec_T0[fmr] <- fmr

  fmr <- c(which.na(vec_D), which(vec_D == ""))
  if(length(fmr) > 0) vec_D[fmr] <- 0

  # pour le temps
  vec_T0 <- as.numeric(vec_T0)
  Ti <- L$Tinit$timing
  Tr <- L$Trecalc$timing

  # pour la date
  # vec_D <- as.numeric(vec_D)
  Di <- L$Tinit$date
  # Dr <- L$Trecalc$date
  # subtract(Di[[17]][1], Di[[1]][1]) %>% as.numeric()

  for(i in 1:length(vec_T0)){
    # the time between acquisition at switch and Tref
    fmr <- difftime(Di[[i]][1], Di[[vec_T0[i]]][1], units = "secs")

    # apply the delta between T0 and Tref
    Tr[[i]] <- magrittr::add(Ti[[i]],fmr)

    # # Dref become D0 for acquisition
    # Dr[[i]] <- subtract(Di[[i]],fmr)
    #
    # # apply the delta in seconde
    # Dr[[i]] <- add(Dr[[i]],vec_D[i])
  }

  # L$Trecalc$date <- Dr
  L$Trecalc$date <- L$Tinit$date
  L$Trecalc$timing <- Tr
  L <- wf.update("re.calc.T.para","sp", L)
  return(L)
}

#### Shift of x mass ####

# mass shift
#' Title
#'
#' @param Li
#'
#' @return
#'
#' @examples
mass.shift <- function(Li){
  min_xMS <- lapply(Li, length.xMS) %>% unlist()

  if(length(unique(min_xMS)) == 1){
    return(Li[[1]]$xMS)
  }else{
    min_xMS <- min(min_xMS)

    mat_xMS <- Li[[1]]$xMS[1:min_xMS]
    for(i in 2:length(Li)) mat_xMS <- rbind(mat_xMS, Li[[i]]$xMS[1:min_xMS])

    fmr <- seq(50,500,10) %>% sapply(det.c, vec = mat_xMS[1,])
    diff_xMS <- t(mat_xMS[,fmr])-mat_xMS[1,fmr]
    ind_diff <- which(diff_xMS[1,] != 0)
    n_df <- length(ind_diff)
    sp_name <- sapply(Li,conc.lst, elem = 1)

    png("Figures/Control/files_shifted.png", width = 400, height = 350)
    par(mar = c(3,3,3,0.1), mgp = c(2,1,0), cex = 1.5)
    matplot(mat_xMS[1,fmr], diff_xMS[,ind_diff], type = "l", lty = 1,
            col = viridis(n_df, direction = -1), lwd = 2,
            ylab = "Mass shift", xlab = "m/z",
            main = paste("file(s) shifted"))
    legend("right",legend = ,sp_name[ind_diff], bty = "n", lty = 1, lwd = 2,
           col = viridis(n_df, direction = -1))
    dev.off()
    print.h("There is a shift in m/z. Look the figure control")
    return(mat_xMS[1,])
  }
}

#### Gestion of name ####

# List of names

#' Title
#'
#' @param f_h5
#'
#' @return
#'
#' @examples
nm.ls <- function(f_h5){
  nm_h5 <- str_remove_all(dir("h5")[f_h5],"_20......_......")
  nm_h5 <- str_remove_all(nm_h5,"20......_......_")
  nm_h5 <- str_remove_all(nm_h5,".h5")
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
# Make the name of samples. The date (20yymmdd_hhmmss.h5)
# is deleting and the acquisitions with the same name.

#### meta data ####

# Export a meta folder empty

#' Title
#'
#' @param L
#'
#' @return
#' @export
#'
#' @examples
empty.meta <- function(L = sp){
  nb_acq <- length(L$names)
  ne <- cumsum(L$nbr_sp)
  ns <- c(1,add(ne,1)[-nb_acq])

  header <- c("names","ID", "nbr_MS", "start", "end", "used", "blank (ID)", "color",
              "concentration","unit","acq_T0 (ID)", "delta_T (s)", "grp1", "grp2", "...")

  mt <- matrix("", nrow = nb_acq, ncol = length(header)-6) %>%
    cbind(L$names, 1:nb_acq, L$nbr_sp, ns, ne, rep(TRUE, nb_acq),.) %>%
    rbind(header,.)

  write.table(mt, file = "meta_empty.csv", sep = ";", dec = ",", row.names = FALSE, col.names = FALSE)

  colnames(mt) <- header
  mt <- mt[-1, -1]
  rownames(mt) <- L$names

  mt[, "color"] <- ctrl.color(mt[,"color"])

  L$mt <- list("name" = "import", "meta" = mt)
  return(L)
}

# Import of meta data

#' Title
#'
#' @param nm
#' @param L
#'
#' @return
#' @export
#'
#' @examples
import.meta <- function(nm = "meta_empty", L = sp){

  mt <- read.table(paste0(nm,".csv"), sep = ";", dec = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  colnames(mt)[1:11] <- c("ID", "nbr_MS", "start", "end", "used", "blank (ID)", "color",
                          "concentration","unit","acq_T0 (ID)", "delta_T (s)")
  mt <- as.matrix(mt)
  fmr <- as.logical(mt[,"used"])
  mt[is.na(mt)==TRUE] <- ""
  mt[,"used"] <- fmr
  mt[, "color"] <- ctrl.color(mt[, "color"])

  L$mt <- list("name" = nm, "meta" = mt)
  L <- wf.update("import.meta",nm, L)

  L$acq <- as.logical(mt[,"used"]) %>% which.equal(1)
  s_acq <- as.numeric(mt[,"start"])
  e_acq <- as.numeric(mt[,"end"])

  L$Sacq <- NULL
  for(i in L$acq) L$Sacq <- c(L$Sacq, seq(s_acq[i],e_acq[i]))

  return(L)
}

#### Provoc_02_gest_meta ####
#### Importation ####

# read.h5

#' Title
#'
#' @param num_fil
#' @param ll
#'
#' @return
#'
#' @examples
read.h5 <- function(num_fil=1, ll = f_h5){

  # find the name of file
  name_h5 <- nm.ls(num_fil)

  # files import
  act_h5 <- paste0("h5/",dir("h5")[num_fil]) %>% H5Fopen()

  # abscissa extraction [~160 000 pts]
  xMS <- act_h5$FullSpectra$MassAxis

  # intensity extraction
  all_MS <- act_h5$FullSpectra$TofData[,1,,]
  # intensities extraction [ acquisition number * 160 000 pts]
  fmr <- dim(all_MS)
  dim(all_MS) <- c(fmr[1], fmr[2] * fmr[3]) # for 2d array

  # timing extraction
  all_timing <- act_h5$TimingData$BufTimes
  dim(all_timing) <- c(fmr[2] * fmr[3])

  # date extraction
  all_date <- acq.time(act_h5) + all_timing

  #TPS2
  all_TPS2 <- act_h5$TPS2$TwData
  fmr <- dim(all_TPS2)
  dim(all_TPS2) <- c(fmr[1], fmr[2] * fmr[3])
  row.names(all_TPS2) <- act_h5$TPS2$TwInfo

  # File close
  H5Fclose(act_h5)

  # reduction
  fmr <- 1:det.c(xMS,50)
  xMS <- as.vector(xMS[-fmr])
  MS <- all_MS[-fmr,]

  # print the working progress and the time code
  print.h(paste0(name_h5, " # ",which(num_fil == ll), "/", length(ll)))

  # return
  list("name" = name_h5,
       "xMS" = xMS,
       "MS" = MS,
       "date" = all_date,
       "timing" = all_timing,
       "nbr_sp" = ncol(MS),
       "meta" = all_TPS2)
}

#' Title
#'
#' @param wdir
#'
#' @return
#' @export
#'
#' @examples
import.h5 <- function(wdir = wd){

  if(("Figures" %in% dir())==FALSE){
    dir.create("Figures")
    dir.create("Figures/Control")
  }

  if(("h5" %in% dir())==FALSE){
    cat("Sorry but the import can't continue. Create a \"h5\" folder with all .h5 fills that you
        want analyse.")
  }

  # data importation ####
  f_h5 <- dir("h5") %>% grep(".h5",.)     # localise h5 files

  length(citation.list) %>% sample(1) %>% citation.list[[.]] %>% cat()
  cat(" \n - - - - - - - - - - - - - - - \n")
  list_h5 <- lapply(f_h5, read.h5, ll = f_h5)

  # formating of sp list ####
  sp <- list()
  sp$names <- sapply(list_h5,conc.lst, elem = 1)
  sp$Tinit$date <- sapply(list_h5,conc.lst, elem = 4, simplify = FALSE)
  sp$Tinit$timing <- sapply(list_h5,conc.lst, elem = 5, simplify = FALSE)
  sp$nbr_sp <- sapply(list_h5,conc.lst, elem = 6)
  sp$meta <- sapply(list_h5,conc.lst, elem = 7, simplify = FALSE)

  sp$xMS <- mass.shift(list_h5)
  print.h("Concatene MS")

  sp$MS <- list()
  for(i in 1:length(list_h5)){
    sp$MS <- c(sp$MS, list(list_h5[[i]]$MS[1:length(sp$xMS),]))
    list_h5[[i]] <- 0
  }
  sp$MS <- do.call(cbind,sp$MS) %>% t()
  remove(list_h5)

  # size reduction ####
  print.h("Reduction")
  sp <- red.xMS(sp)

  # create MassSpectrum object ####
  print.h("Create MassSpectrum object")
  sp$MS <- apply(sp$MS,1, create_local_MS, xMS = sp$xMS)

  # smooth spectra ####
  print.h("Smooth spectra")
  oldw <- getOption("warn")
  options(warn = -1)
  sp$MS <- smoothIntensity(sp$MS,
                           method = "SavitzkyGolay",
                           halfWindowSize = 3)
  options(warn = oldw)

  # align spectra ####
  print.h("Align spectra")
  sp$names_acq <- prep.names(sp) %>% apply(2,names.samples)
  sp$MS <- alignSpectra(sp$MS, tolerance = 0.02)
  sp$MS <- averageMassSpectra(sp$MS, labels = convertStr2List(sp), method="mean")

  # peak detection ####
  print.h("Peak detection")

  sp$peaks <- detectPeaks(sp$MS, method="MAD", halfWindowSize=20, SNR=5)
  sp$peaks <- binPeaks(sp$peaks, tolerance=0.01)
  sp$peaks <- MALDIquant::filterPeaks(sp$peaks, minFrequency=0.25)
  sp$peaks <- intensityMatrix(sp$peaks, sp$MS)
  colnames(sp$peaks) <- colnames(sp$peaks) %>% as.numeric() %>% round(3)

  sp$xMS <- sapply(sp$MS, mass.spectra) %>% rowMeans() %>% round(3)
  sp$MS <- sapply(sp$MS, mat.spectra)

  rownames(sp$MS) <- sp$xMS
  colnames(sp$MS) <- unlist(sp$names_acq)
  rownames(sp$peaks) <- colnames(sp$MS)

  # export meta folder and finish ####
  sp$Trecalc <- sp$Tinit
  sp$workflow <- wdir
  names(sp$workflow)[[1]] <- "import.h5"
  sp$wd <- wdir
  sp$acq <- 1:length(sp$names)
  sp$Sacq <- 1:ncol(sp$MS)

  sp <- empty.meta(sp)
  sp <- list.order(sp)

  print.h("Import is completed")
  return(sp)
  # import function is finished ####
}

#' Title
#'
#' @param L
#'
#' @return
#'
#' @examples
red.xMS <- function(L=sp){
  maxMS <- apply(L$MS,2,max) # Imax for each mass
  fmr <- which(maxMS < 500)
  dMS <- density(maxMS[fmr], bw = 0.001)
  thr <- dMS$x[which.max(dMS$y)] %>% round(0) %>% multiply_by(2.5) %>% round(0)

  indinf <- which(maxMS < thr) # column where zero mass is superior at threshold
  nbi <- length(indinf)
  diffind <- subtract(indinf[-1], indinf[-nbi])
  fr <- which(diffind > 1)
  fr_t <- c(fr, fr+1) %>% sort()
  br <- sapply(-10:10, add, e2 = indinf[fr_t]) %>% as.vector() %>% sort() %>% unique()
  kp <- which(indinf %in% br)
  ind_del <- indinf[-kp]

  diffind <- subtract(ind_del[-1], ind_del[-length(ind_del)])
  fr <- which(diffind > 1)
  fr_t <- c(fr, fr+1) %>% sort()
  inull <- ind_del[-fr_t]
  idel <- ind_del[fr_t]

  L$MS[,idel] <- rep(0,nrow(L$MS))

  {
    # # graphe de vision
    # tiff(paste("Densité des masses maximum.tiff"))
    #  plot(dMS, main="Density",xlim = c(0, thr*2))
    #  abline(v = thr, lty = 2, lwd = 2, col = "red")
    #  legend("topright", bty = "n",
    #         legend = c(paste("threshold = ",round(thr,0)),
    #                    paste("nbr of mass deleted =", length(ind_del))))
    # dev.off()
    #
    # plot.threshold <- function(br, L= sp, z = c(0,300), ind_d = indinf, ind_k = indinf[kp],inu = inull){
    #  tiff(paste("Vue des masses supprimées de",br[1],"à",br[2],"Da.tiff"))
    #    a <- det.c(br[1],L$xMS):det.c(br[2],L$xMS)
    #    matplot(sp$xMS[a],maxMS[a], type = "l", ylim = z,
    #            xlab = "m/z (Da)", ylab = "intensité (u.a.)",
    #            main = "Spectre de l'intensité maximale de chaque masse")
    #    legend("topleft", bty = "n", lty = 1, col = c("black","blue","red"),
    #           legend = c("spectre max","masses supprimées", "masses gardées"))
    #    abline(h = thr, lty = 2, lwd = 0.8)
    #    b <- det.c(br[1],L$xMS[-inu]):det.c(br[2],L$xMS[-inu])
    #    matplot(L$xMS[-inu][b], t(L$MS[,-inu][,b]), type = "l", lty = 1,
    #            col = viridis(n = nrow(L$MS), alpha = 0.2), add = TRUE)
    #    matplot(L$xMS[ind_d],maxMS[ind_d], type = "l", add = TRUE, lwd = 2, col = "blue")
    #    matplot(L$xMS[ind_k],maxMS[ind_k], type = "l", lwd = 2, add = TRUE, col = alpha("red",0.5))
    #  dev.off()
    # }
    #
    # plot.threshold(L = L, br = c(50.9, 51.2))
    # plot.threshold(L = L, br = c(66.7, 67.4))
    # plot.threshold(L = L, br = c(64.1, 64.7), z = c(40,150))
    # plot.threshold(L = L, br = c(236, 238))
    # plot.threshold(L = L, br = c(340.5, 341.5))
    #
  }

  # mise en forme finale

  L$MS <- L$MS[,-inull]
  L$xMS <- L$xMS[-inull]

  return(L)
}

#### Provoc_03_gest_meta ####
#### Plot spectra ####
# Plot spectra dynamic

#' Title
#'
#' @param sel_sp
#' @param L
#' @param new_color
#'
#' @return
#' @export
#'
#' @examples
dy.spectra <- function(sel_sp = sp$mt$meta[sp$acq,"end"], L = sp, new_color = FALSE){
  if(is.character(sel_sp) == TRUE) sel_sp <- as.numeric(sel_sp)
  if(length(sel_sp) > 30) print("Caution /!\ The number of spectra is too big. Select less spectra.")
  if(length(sel_sp) <=30){

    sp_sel <- L$MS[,sel_sp]

    if(new_color == TRUE) dy_color <- viridis(length(sel_sp),alpha = 0.8)
    if(new_color == FALSE) dy_color <- rep.mtm("color", L, sel = "all")[sel_sp]

    if(length(sel_sp)==1){
      dysp <- sp_sel %>% cbind(L$xMS,.) %>% as.data.frame()
      colnames(dysp)[2] <- colnames(L$MS)[sel_sp]
      titre <- colnames(L$MS)[sel_sp]
    }else{
      dysp <- cbind(L$xMS,sp_sel) %>% as.data.frame()
      titre <- "sp_align"
    }
    rownames(dysp) <- L$xMS

    ftitre <- paste0(getwd(),"/Figures/") %>% dir() %>% grep(titre, .) %>% length() %>% add(1)
    ftitre <- paste0(getwd(),"/Figures/dy_",titre,"_",ftitre)
    fmr <- system.file("rmd", "print_dy_sp.Rmd", package = "provoc")
    rmarkdown::render(input = fmr, output_file = ftitre)
  }
}

# Plot spectra tiff

#' Title
#'
#' @param sel_sp
#' @param pkm
#' @param pkM
#' @param L
#' @param new_title
#' @param new_color
#' @param leg
#'
#' @return
#' @export
#'
#' @examples
fx.spectra <- function(sel_sp = sp$mt$meta[sp$acq,"end"], pkm = 59, pkM = 205,
                       L = sp, new_title = "fx_spectra", new_color = FALSE, leg = "r"){
  # check the selection
  if(is.character(sel_sp) == TRUE) sel_sp <- as.numeric(sel_sp)
  if(length(sel_sp) > 30) print("Caution /!\ The number of spectra is too big. Select less spectra.")

  # create variable for optimize zoom
  xmin <- pkm
  xmax <- pkM
  cmin <- det.c(pkm - 0.3, L$xMS)
  cmax <- det.c(pkM + 0.3, L$xMS)

  # selection of major peaks
  pk_short_list <- pk.short(pk_mat = L$peaks[c(sel_sp,sel_sp),])

  fmr <- pk_short_list[2,]
  ind_pk <- which.sup(fmr,mean(fmr))
  if(length(fmr)/2 < length(ind_pk)) ind_pk <- which.sup(fmr,median(fmr))
  pk_short_list <- pk_short_list[,ind_pk]

  # select color
  if(new_color == TRUE)  fx_color <- viridis(length(sel_sp),alpha = 0.8)
  if(new_color == FALSE) fx_color <- rep.mtm("color", L, sel = "all")[sel_sp]

  # define title
  if(length(sel_sp)==1) new_title <- paste0("fx_",colnames(L$MS)[sel_sp])

  ntitre <- paste0(getwd(),"/Figures/") %>% dir() %>% grep(new_title, .) %>% length() %>% add(1)
  ntitre <- paste0("Figures/",new_title,"_",ntitre,"_zoom_",xmin,"_to_",xmax,".tif")

  # calculate the max spectra
  if(length(sel_sp)==1) sp_max <- L$MS[cmin:cmax,sel_sp]
  if(length(sel_sp) >1) sp_max <- apply(L$MS[cmin:cmax,sel_sp],1,max)
  ymax <- max(sp_max)*1.1

  # define the legend postion
  leg_pos <- "topright"
  if(leg == "l") leg_pos <- "topleft"

  # the plot
  tiff(filename = ntitre, width = 1000, height = 580)
  par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

  # a plot with subline
  matplot(L$xMS[cmin:cmax], sp_max, type = "l",
          col = alpha("turquoise2",0.5), lwd = 5, ylim = c(0,ymax),
          xlab = "m/z", ylab = "Relative intensity (u.a.)",
          main = new_title, xaxt="n")

  # the MS plot
  matplot(L$xMS[cmin:cmax], L$MS[cmin:cmax,sel_sp],
          type = "l", col = fx_color, add = TRUE)

  # legend
  if((leg != "n")&(length(sel_sp) >1)){
    legend(leg_pos, bty = "n", col = fx_color,
           legend = colnames(L$MS)[sel_sp], lty = 1)
  }

  # axis ...
  axis(side=1,0:600, tcl=0,labels=FALSE)
  axis(side=2,(-ymax:ymax)*2,tcl=0,labels=FALSE)
  axis(side=3,0:600, tcl=0,labels=FALSE)
  axis(side=4,(-ymax:ymax)*2,tcl=0,labels=FALSE)

  # ... and ticks
  axis(1, at = seq(dizaine(xmin), dizaine(xmax), 10), lwd.ticks = 2, tck = -0.03)
  axis(1, at = seq(dizaine(xmin), dizaine(xmax) +10, 5), labels = FALSE, tck = -0.03)
  axis(1, at = seq(dizaine(xmin), dizaine(xmax) +10, 1), labels = FALSE, tck = -0.01)
  text(pk_short_list[1,], pk_short_list[2,], labels = pk_short_list[1,], cex = 0.8, pos = 3, offset = 0.5)
  # end of plot
  dev.off()
}


#### Plots kinetic ####

#' Title
#'
#' @param M_num
#' @param each_mass
#' @param group
#' @param graph_type
#' @param L
#' @param Y_exp
#' @param time_format
#'
#' @return
#' @export
#'
#' @examples
kinetic.plot <- function(M_num = M.Z(c(59, 137)), each_mass = TRUE,
                         group = FALSE, graph_type = "dy", L = sp,
                         Y_exp = FALSE, time_format = "date"){
  # Mise en forme :
  if(("kinetic" %in% dir("Figures"))==FALSE) dir.create("Figures/kinetic")
  vp <- list(exp = Y_exp, time = time_format, grp = group)

  # Graphe :

  # Etape 1 : chaque mass ?
  if(each_mass == TRUE){
    # OUI

    for(ma in M_num){ # ma = M_num[1]

      # Etape 2 : chaque groupe ?
      if(group != FALSE){
        # OUI
        grp <- L$mt$meta[,group][L$acq] %>% as.character() # on repere les groupes
        for(u in unique(grp)){ # u = grp[1]

          ind_PK <- L$acq[which(u == grp)]              # on repere les indices de chaques groupes

          # Etape 3 : plot statique ou dynamique ?
          if(graph_type == "fx"){
            # plot fixe

            titre <- c("Figures/kinetic/pk_at", ma,"of",u) %>%
              str_flatten(" ") %>% paste0("_",vp$time,".tiff")
            fx.kinetic.plot(L, titre, acq = ind_PK, MA = ma, VP = vp)

          }else if(graph_type == "dy"){
            # plot dynamique
            titre <- str_flatten(ma, collapse = " ") %>% paste("peak at",.,"of", u, vp$time)
            dy.kinetic.plot(L, titre, acq = ind_PK, MA = ma, VP = vp)
          }
        }

      }else{
        # NON

        # Etape 3 : plot dynamique ou statique ?
        if(graph_type == "fx"){
          # plot statique

          titre <- c("Figures/kinetic/pk_at", ma,"of",head(row.names(L$mt$meta)[L$acq])) %>%
            str_flatten(" ") %>% paste0("_",vp$time,".tiff")
          fx.kinetic.plot(L, titre, acq = L$acq, MA = ma, VP = vp)
        }else if(graph_type == "dy"){
          # plot dynamique
          titre <- str_flatten(ma, collapse = " ") %>% paste("peak at",.,"of all selected", vp$time)
          dy.kinetic.plot(L, titre, acq = L$acq, MA = ma, VP = vp)
        }
      }
    }
  }else{
    # NON
    ma <- M_num

    # Etape 2 : chaque groupe ?
    if(group != FALSE){
      # OUI

      grp <- L$mt$meta[,group][L$acq] %>% as.character() # on repere les groupes
      for(u in unique(grp)){# u = unique(grp)[1]

        ind_PK <- L$acq[which(u == grp)]              # on repere les indices de chaques groupes

        # Etape 3 : plot statique ou dynamique ?
        if(graph_type == "fx"){
          # plot fixe

          titre <- c("Figures/kinetic/pk_at", ma,"of",u) %>%
            str_flatten(" ") %>% paste0("_",vp$time,".tiff")
          fx.kinetic.plot(L, titre, acq = ind_PK, MA = ma, VP = vp)

        }else if(graph_type == "dy"){
          # plot dynamique
          titre <- str_flatten(ma, collapse = " ") %>% paste("peak at",.,"of", u, vp$time)
          dy.kinetic.plot(L, titre, acq = ind_PK, MA = ma, VP = vp)
        }
      }
    }else{
      # NON

      # Etape 3 : plot statique ou dynamique ?
      if(graph_type == "fx"){
        # plot statique

        titre <- c("Figures/kinetic/pk_at", ma,"of",head(row.names(L$mt$meta)[L$acq])) %>%
          str_flatten(" ") %>% paste0("_",vp$time,".tiff")
        fx.kinetic.plot(L, titre, acq = L$acq, MA = ma, VP = vp)

      }else if(graph_type == "dy"){
        # plot dynamique
        titre <- str_flatten(ma, collapse = " ") %>% paste("peak at",.,"of all selected", vp$time)
        dy.kinetic.plot(L, titre, acq = L$acq, MA = ma, VP = vp)
      }
    }
  }
}

#' Title
#'
#' @param L
#' @param titre
#' @param acq
#' @param MA
#' @param VP
#'
#' @return
#'
#' @examples
fx.kinetic.plot <- function(L, titre, acq = ind_PK, MA = ma, VP = vp){
  # index for peaks and acquisitions
  ind_pk <- which(colnames(L$peaks) %in% MA)
  ind_Sacq <- ind.acq(acq,L)

  # intensity
  Ibn <- c(0, max(L$peaks[ind_Sacq, ind_pk]))

  # time
  if(VP$time == "time"){
    Tbn <- c(0,0)
    for(i in acq) Tbn[2] <- max(Tbn[2], L$Trecalc$timing[[i]])

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
    Tbn <- Tbn/Tdiv

    Xlab <- paste0("Time (",unit,")")
    List_abs <- lapply(L$Trecalc$timing, divide_by, e2 = Tdiv)
  }

  # or date
  if(VP$time == "date"){
    Cdate <- as.POSIXct(Sys.time(), format="%m/%d/%Y %H:%M:%S")
    attr(Cdate, "tzone") <- "UTC"
    for(i in acq) Cdate <- c(Cdate, L$Trecalc$date[[i]])
    Cdate <- Cdate[-1]
    Tbn <- range(Cdate)

    Xlab <- paste0("Date")
    List_abs <- L$Trecalc$date
  }

  pk_col <- L$mt$meta[,"color"]
  if(VP$exp == FALSE) VP$exp <- ""
  if(VP$exp == TRUE){
    VP$exp <- "y"
    Ibn[1] <- 10
  }

  tiff(file = titre, width = 1200, height = 600,units = "px")
  par(mar = c(5,5,2,16),mgp = c(3.5,1.5,0),xpd = NA,
      cex.main=2, cex.lab = 2, cex.axis = 2)

  if(VP$time == "date"){
    matplot(Tbn, Ibn, type = "l", col = "white", log = VP$exp,
            xlab = Xlab, xaxt = "n", ylab = "Intensity (a.u.)")
    axis.POSIXct(1, x =  Cdate)
  }else{
    matplot(Tbn, Ibn, type = "l", col = "white", log = VP$exp,
            xlab = Xlab, ylab = "Intensity (a.u.)")
  }

  nq <- 0
  for(i in acq){ # i=1
    cl <- 15
    nq <- nq + 1
    for(j in ind_pk){ # j = ind_pk[2]
      coor <- ind.acq(i,L)

      if(length(coor)>1){
        fmr <- length(coor) + (1-nq)*round(length(coor)/30)
        matplot(List_abs[[i]], L$peaks[coor,j], type = "l", lwd = 2,
                col = pk_col[i], add = TRUE)
        matplot(List_abs[[i]], L$peaks[coor,j],
                pch = cl, col = pk_col[i], cex = 2, add = TRUE)
      }else{
        matplot(List_abs[[i]], L$peaks[coor,j],
                pch = cl, col = pk_col[i], add = TRUE)
      }
      cl <- cl + 1
    }
  }

  if(length(acq) <= 10){
    l.acq <- acq
  }else{
    fmr <- length(acq) %>% subtract(4)
    l.acq <- acq[c(1:5, fmr:length(acq))]
  }

  if(length(MA)<= 10){
    l.num <- MA
  }else{
    fmr <- length(MA) %>% subtract(4)
    l.num <- MA[c(1:5, fmr:length(MA))]
  }

  legend("topright", bty = "n", cex = 1.5, xpd = NA, inset = c(-0.26,0),
         legend = c("Sample(s) :", L$names[l.acq]," ","Masse(s) :", l.num),
         lty = c(NA,rep(1,length(l.acq)), NA, NA, rep(NA,length(l.num))), lwd = 2,
         pch = c(NA,rep(NA,length(l.acq)), NA, NA, 14 + seq(1,length(l.num))),
         col = c(NA,pk_col, NA, NA, rep("black", length(l.num))))

  dev.off()
}

#' Title
#'
#' @param L
#' @param titre
#' @param acq
#' @param MA
#' @param VP
#'
#' @return
#'
#' @examples
dy.kinetic.plot <- function(L, titre, acq = ind_PK, MA = ma, VP = vp){
  # index for peaks and acquisitions
  ind_pk <- which(colnames(L$peaks) %in% MA)
  ind_Sacq <- ind.acq(acq,L)

  # intensity
  Ibn <- c(0, max(L$peaks[ind_Sacq, ind_pk]))

  # time
  if(VP$time == "time"){
    Tbn <- c(0,0)
    for(i in acq) Tbn[2] <- max(Tbn[2], L$Trecalc$timing[[i]])

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
    Tbn <- Tbn/Tdiv

    Xlab <- paste0("Time (",unit,")")
    List_abs <- lapply(L$Trecalc$timing, divide_by, e2 = Tdiv)
  }

  # or date
  if(VP$time == "date"){
    Cdate <- as.POSIXct(Sys.time(), format="%m/%d/%Y %H:%M:%S")
    attr(Cdate, "tzone") <- "UTC"
    for(i in acq) Cdate <- c(Cdate, L$Trecalc$date[[i]])
    Cdate <- Cdate[-1]
    Tbn <- range(Cdate)

    Xlab <- paste0("Date")
    List_abs <- L$Trecalc$date
  }

  # pk_col <- L$mt$meta[,"color"]
  if(VP$exp == FALSE) VP$exp <- ""
  if(VP$exp == TRUE){
    VP$exp <- "y"
    Ibn[1] <- 10
  }

  # the list of data
  dy_mat <- sapply(acq, dy.mat.pk, ipk = ind_pk, La = List_abs, Li = L,
                   vp = VP, simplify = FALSE)

  # convert to data.frame
  ldy <- length(dy_mat)
  if(ldy == 1){
    dysp <- dy_mat[[1]] %>% as.data.frame()
  }else{
    dysp <- merge(dy_mat[[1]], dy_mat[[2]], all=TRUE)
    if(ldy >= 3){
      for(i in 3:ldy){
        dysp <- merge(dysp, dy_mat[[i]], all=TRUE)
      }
    }
  }
  #dysp <- dysp[order(dysp[,3], dysp[,1]),]
  dysp$ID <- NULL

  # and plot
  if(VP$time == "date"){
    dysp$xT <- as.POSIXct(dysp$xT, origin = "1970-01-01", tz = "GMT")
    dysp <- xts::xts(dysp[,-1], order.by = dysp$xT, tz="GMT")

    fmr <- system.file("rmd", "print_dypk_date.Rmd", package = "provoc")
    rmarkdown::render(input = fmr,
                      output_file = paste0(L$wd, "/Figures/kinetic/", titre))
  }else{
    fmr <- system.file("rmd", "print_dypk_time.Rmd", package = "provoc")
    rmarkdown::render(input = fmr,
                      output_file = paste0(L$wd, "/Figures/kinetic/", titre))
  }
}

#' Title
#'
#' @param ac
#' @param ipk
#' @param La
#' @param Li
#' @param vp
#'
#' @return
#'
#' @examples
dy.mat.pk <- function(ac = acq, ipk = ind_pk, La = List_abs, Li = L, vp = VP){

  if(vp$grp == FALSE) nid <- rownames(Li$mt$meta)[ac]
  if(vp$grp != FALSE) nid <- Li$mt$meta[ac,vp$grp]

  fmr <- rep(Li$mt$meta[ac,"ID"], length(La[[ac]])) %>% as.numeric()
  fmr <- cbind(La[[ac]],fmr, Li$peaks[ind.acq(ac,Li), ipk])
  colnames(fmr) <- c("xT","ID", paste(nid, colnames(Li$peaks)[ipk]))

  return(fmr)
}

#### Provoc_04_gest_meta ####
#### micro-functions ####

# print garbage collection

#' Title
#'
#' @return
#'
#' @examples
print.gc <- function(){
  fmr <- memory.size()
  gc()
  paste0("RAM : ",fmr," -> gc -> ", memory.size()) %>% print.h()
}

#concatenation

#' Title
#'
#' @param list_n
#' @param elem
#'
#' @return
#'
#' @examples
conc.lst <- function(list_n, elem = 1){
  list_n[[elem]]
}

#concatenation

#' Title
#'
#' @param list_n
#' @param elem
#'
#' @return
#'
#' @examples
dim.lst <- function(list_n, elem = 1){
  dim(list_n[[elem]])
}

#names of acquisition

#' Title
#'
#' @param L
#'
#' @return
#'
#' @examples
prep.names <- function(L){
  fmr <- log10(L$nbr_sp) %>% floor() %>% add(1)
  rbind(fmr, L$names, L$nbr_sp)
}

#' Title
#'
#' @param vec
#'
#' @return
#'
#' @examples
names.samples <- function(vec){str_pad(1:vec[3],vec[1], pad = "0") %>% paste(vec[2],.,sep = "_")}

#' Title
#'
#' @param L
#'
#' @return
#'
#' @examples
convertStr2List <- function(L){
  plip <- function(vec) return(vec)
  fmr <- unlist(L$names_acq) %>% lapply(plip)
  names(fmr) <- unlist(L$names_acq)
  return(fmr)
}

#control colo

#' Title
#'
#' @param vec_col
#'
#' @return
#'
#' @examples
ctrl.color <- function(vec_col = mt[,"color"]){

  fmr <- which.na(vec_col == "")
  if(length(fmr) > 0) vec_col[fmr] <- viridis(length(fmr)) %>% alpha(0.5)

  fmr <- which(vec_col == "")
  if(length(fmr) > 0) vec_col[fmr] <- viridis(length(fmr)) %>% alpha(0.5)

  return(vec_col)
}

# create Mass Spectrum objet
#' Title
#'
#' @param MS
#' @param xMS
#'
#' @return
#'
#' @examples
create_local_MS <- function(MS, xMS){createMassSpectrum(xMS,MS)}

# return to spectra

#' Title
#'
#' @param spobj
#'
#' @return
#'
#' @examples
mat.spectra <- function(spobj){spobj@intensity}

#' Title
#'
#' @param spobj
#'
#' @return
#'
#' @examples
mass.spectra <- function(spobj){spobj@mass}

# retourne l'index valable d'une borne le long d'un vecteur.

#' Title
#'
#' @param brn
#' @param vec
#'
#' @return
#' @export
#'
#' @examples
det.c <- function(brn,vec){subtract(vec,brn) %>% sapply(abs) %>% which.min()}

# detect length of x mass

#' Title
#'
#' @param splist
#'
#' @return
#'
#' @examples
length.xMS <- function(splist){length(splist$xMS)}

# print le texte suivi de l'heure

#' Title
#'
#' @param txt
#'
#' @return
#' @export
#'
#' @examples
print.h <- function(txt = "hello there"){heure() %>% paste0(txt,", ",.) %>% print()}

# trouve les pics proches

#' Title
#'
#' @param pk_x
#' @param mat
#' @param w.sub
#'
#' @return
#'
#' @examples
pk.red <- function(pk_x = pk_max[1,1], mat = pk_max, w.sub = 4){
  fmr <- subtract(mat[1,],pk_x) %>% abs() %>% multiply_by(-1) %>% which.sup(-(w.sub+1))
  return(fmr[which.max(mat[2,fmr])])
}

# list des peaks principaux

#' Title
#'
#' @param pk_mat
#'
#' @return
#'
#' @examples
pk.short <- function(pk_mat = L$peaks){
  pk_max <- colnames(pk_mat) %>% as.numeric() %>% rbind(apply(pk_mat,2,max))
  fmr <- sapply(pk_max[1,], pk.red, mat = pk_max, w.sub = 4) %>% unique() %>% sort()
  pk_mat <- pk_max[,fmr]
  rownames(pk_mat) <- c("m/z","int")
  return(pk_mat)
}

# donne l'heure
#' Title
#'
#' @return
#'
#' @examples
heure <- function(){str_split(Sys.time(),pattern = " ")[[1]][2]}

# order the list

#' Title
#'
#' @param L
#'
#' @return
#'
#' @examples
list.order <- function(L = sp){
  L <- list("MS" = L$MS,
            "peaks" = L$peaks,
            "xMS" = L$xMS,
            "names" = L$names,
            "wd" = L$wd,
            "acq" = L$acq,
            "Sacq" = L$Sacq,
            "nbr_sp" = L$nbr_sp,
            "names_acq" = L$names_acq,
            "Tinit" = L$Tinit,
            "Trecalc" = L$Trecalc,
            "workflow" = L$workflow,
            "mt" = L$mt,
            "meta" = L$meta)
}

# named workflow

#' Title
#'
#' @param nwf
#' @param L
#'
#' @return
#'
#' @examples
name.wf <- function(nwf = "randow", L = sp){
  fmr <- length(L$workflow)
  names(L$workflow)[[fmr]] <- nwf
  return(L)
}

# update workflow

#' Title
#'
#' @param nm_wf
#' @param obj_wf
#' @param L
#'
#' @return
#'
#' @examples
wf.update <- function(nm_wf, obj_wf, L = sp){
  L$workflow <- c(L$workflow, list(obj_wf))
  L <- name.wf(nm_wf, L)
  return(L)
}

# Repet meta parameter

#' Title
#'
#' @param col.nam
#' @param L
#' @param sel
#'
#' @return
#'
#' @examples
rep.mtm <- function(col.nam, L, sel = "acq"){
  fmr <- L$acq
  if(sel == "all") fmr <- as.numeric(L$mt$meta[,"ID"])
  sapply(fmr, rep.mtu, col.nam = col.nam, L = L, simplify = FALSE) %>% unlist()
}

# Repet meta parameter of a single aquisition

#' Title
#'
#' @param acq
#' @param col.nam
#' @param L
#'
#' @return
#'
#' @examples
rep.mtu <- function(acq, col.nam, L){
  fmr <- which(col.nam == colnames(L$mt$meta))
  rep(L$mt$meta[acq,fmr], L$nbr_sp[acq])
}

# arrondi a la dizaine inferieure.

#' Title
#'
#' @param x
#'
#' @return
#'
#' @examples
dizaine <- function(x){
  eph <- log(x,10) %>% floor() %>% multiply_by(10)
  divide_by(x,eph) %>% floor() %>% multiply_by(eph)
}

# search all peak in accord to a mass number

#' Title
#'
#' @param ma
#' @param L
#'
#' @return
#' @export
#'
#' @examples
M.Z <- function(ma,L=sp){
  vec_pk <- colnames(L$peaks) %>% as.numeric()
  fmr <- NULL
  for(maz in ma) fmr <- c(fmr, which((vec_pk < maz + 0.5)&(vec_pk > maz - 0.5)))
  return(vec_pk[fmr])
}

# return index of spectra for each acquistion

#' Title
#'
#' @param n_acq
#' @param L
#'
#' @return
#'
#' @examples
ind.acq <- function(n_acq,L){
  fmr <- NULL
  mat_mt <- cbind(as.numeric(L$mt$meta[,"start"]),
                  as.numeric(L$mt$meta[,"end"]))
  for(i in n_acq) fmr <- c(fmr, mat_mt[i,1]:mat_mt[i,2])
  return(fmr)
}

#### which pack ####

# equal :

#' Title
#'
#' @param vec
#' @param nb
#'
#' @return
#' @export
#'
#' @examples
which.equal <- function(vec,nb){which(vec == nb)}
# retourne la position des elements de vec egaux a nb.

# na :

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
which.na <- function(x){which(is.na(x) == TRUE)}
# retourne la position des NA sur un vecteur.

# not na :

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
which.not.na <- function(x){which(is.na(x) == FALSE)}
# retourne la position des non-NA sur un vecteur.

# superior :

#' Title
#'
#' @param vec
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
which.sup <- function(vec, threshold){return(which(vec > threshold))}
# retourne la position des elements de vec superieur au seuil.

#### End of Code ####
