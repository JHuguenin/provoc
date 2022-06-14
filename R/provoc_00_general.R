#### Provoc_00_general ####

#' provoc : Perform a Rapid Overview for the Volatils Organic Compounds
#'
#' analyze data of VOC by PTR-ToF-MS Vocus
#'
#' @docType package
#' @name provoc
#'
#' @import dygraphs
#' @import graphics
#' @import grDevices
#' @import utils
#' @importFrom ALS als
#' @importFrom baseline baseline.rollingBall
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
#' @importFrom rhdf5 h5closeAll
#' @importFrom rhdf5 h5createFile
#' @importFrom rhdf5 h5ls
#' @importFrom rhdf5 h5write
#' @importFrom rmarkdown render
#' @importFrom scales alpha
#' @importFrom stats approx
#' @importFrom stats density
#' @importFrom stats median
#' @importFrom stringr str_flatten
#' @importFrom stringr str_pad
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_split
#' @importFrom stringr str_sub
#' @importFrom viridis viridis
#' @importFrom xts xts
NULL

#### Fonction inutile ####
citation.list <- {list(
  c("Il faut aller trop loin pour decouvrir les limites.", "Joris Huguenin"),
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
  c("Error in chol.default(Winv) : le mineur dominant d'ordre 118 n'est pas defini positif", "rchemo error message"),
  c("Je suis gentil avec tout le monde, celui qui dit le contraire je lui foutrai mon poing dans la gueule.", "Leo Ferre"),
  c("Le desespoir est une forme superieure de critique.", "Leo Ferre"),
  c("Les diplomes sont faits pour les gens qui n'ont pas de talent.","Pierre Desproges"),
  c("Bal tragique a Colombey, un mort.","Hara Kiri"),
  c("Si la matiere grise etait plus rose, le monde aurait moins les idees noires.","Pierre Dac"),
  c("J'ai pris la decision de ne plus etre influencable. Qu'est-ce que vous en pensez ?","Patrick Sebastien"),
  c("Est-il indispensable d'etre cultive quand il suffit de fermer sa gueule pour briller en societe ?","Pierre Desproges"),
  c("On ne discute pas recettes de cuisine avec des anthropophages.", "Jean-Pierre Vernant"))}

#### Provoc_01_gest_meta ####
#### Gestion of time ####

#' Collecte the metadata in h5 file
#' @param ls.t a h5 file
#' @return an hour
#' @noRd
acq.time <- function(ls.t = ls_h5[[1]]){
  oldw <- getOption("warn")
  options(warn = -1)

  eph <- ls.t$AcquisitionLog$Log$timestring[1] %>%
    as.POSIXct(format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")

  options(warn = oldw)
  return(eph)
}

#' Reinitialize the parameters of time
#'
#' This function allows you to reset the parameters related to the acquisition
#' times. This makes it possible to modify the T0 of an acquisition group. The
#' modifications must be made by the meta file
#' @export
#' @param L sp
#' @return sp
#' @examples
#' # The first time configuration :
#' #
#' # sp <- import.meta("meta_1")
#' # sp <- re.calc.T.para(sp)
#' #
#' # The second time configuration :
#' #
#' # sp <- import.meta("meta_2")
#' # sp <- re.init.T.para(sp)
#' # sp <- re.calc.T.para(sp)
re.init.T.para <- function(L = sp){
  L$Trecalc$date <- L$Tinit$date
  L$Trecalc$timing <- L$Tinit$timing
  L <- wf.update("re.init.T.para","sp", L)
  return(L)
}

#' calculate the parameters of time
#'
#' Be careful. By default, the "time" option uses a relative T0 from the first
#' spectrum of each acquisition and the "date" option uses the actual date and
#' time of each spectrum.
#' Using the acq_T0 column with the "time" option allows different acquisitions
#' to be sequenced using the T0 of the specified acquisition. Using the acq_T0
#' column with the "date" option allows you to overlap acquisitions on the T0
#' of the specified acquisition.
#'
#' The delta_T column is used to add the specified time (in seconds) to the
#' acquisition.
#' @export
#' @param L sp
#' @return sp
#' @examples
#' # The first time configuration :
#' #
#' # sp <- import.meta("meta_1")
#' # sp <- re.calc.T.para(sp)
#' #
#' # The second time configuration :
#' #
#' # sp <- import.meta("meta_2")
#' # sp <- re.init.T.para(sp)
#' # sp <- re.calc.T.para(sp)
re.calc.T.para <- function(L = sp){
  vec_T0 <- L$mt$meta[,"acq_T0 (ID)"]
  vec_D <- L$mt$meta[,"delta_T (s)"]

  fmr <- c(which.na(vec_T0),which(vec_T0 == ""))
  if(length(fmr) > 0) vec_T0[fmr] <- fmr

  fmr <- c(which.na(vec_D), which(vec_D == ""))
  if(length(fmr) > 0) vec_D[fmr] <- 0

  # pour le temps
  vec_T0 <- as.numeric(vec_T0)
  vec_D <- as.numeric(vec_D)
  Ti <- L$Tinit$timing
  Tr <- L$Trecalc$timing

  # pour la date
  Di <- L$Tinit$date
  Dr <- L$Trecalc$date


  # subtract(Di[[17]][1], Di[[1]][1]) %>% as.numeric()

  for(i in 1:length(vec_T0)){ # i=100
    # the time between acquisition at switch and Tref
    fmr <- difftime(Di[[i]][1], Di[[vec_T0[i]]][1], units = "secs")
    fmr <- as.numeric(fmr)
    # apply the delta between T0 and Tref
    Tr[[i]] <- magrittr::add(Ti[[i]],fmr)
    Tr[[i]] <- magrittr::add(Tr[[i]],vec_D[i])

    # Dref become D0 for acquisition
    Dr[[i]] <- magrittr::subtract(Di[[i]],fmr)
    # apply the delta in seconde
    Dr[[i]] <- magrittr::add(Dr[[i]],vec_D[i])
  }

  # x <- (0:19)*12
  # y <- 7:12
  # xx <- lapply(x,add, e2=y) %>% do.call(c,.)
  # yy <- lapply(1:20, function(x) Tr[[x]]) %>% do.call(c,.)
  # plot(xx, yy, pch = 16)

  L$Trecalc$date <- Dr
  # L$Trecalc$date <- L$Tinit$date
  L$Trecalc$timing <- Tr
  L <- wf.update("re.calc.T.para","sp", L)
  return(L)
}

#### Shift of x mass ####

#' Shift abscissa for combine several acquisitions
#' @param Li sp
#' @return matrix
#' @noRd
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
    hprint("There is a shift in m/z. Look the figure control")
    return(mat_xMS[1,])
  }
}

#### Gestion of name ####

#' return the list of names
#' Make the name of samples. The date (20yymmdd_hhmmss.h5)
#'  is deleting and the acquisitions with the same name.
#' @noRd
#' @param f_h5 a string of character
#' @return another string of character
nm.ls <- function(f_h5, wd){
  nm_h5 <- str_remove_all(dir(paste0(wd,"/h5"))[f_h5],"_20......_......")
  nm_h5 <- str_remove_all(nm_h5,"20......_......_")
  nm_h5 <- str_remove_all(nm_h5,"\\.h5")
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

#### meta data ####

#' Export a meta file
#' Execute this function for create a initial file "meta_empty.csv".
#' @param L sp
#' @return sp
#' @export
#' @examples
#' # empty.meta()
empty.meta <- function(L = sp){
  nb_acq <- length(L$names)
  ne <- cumsum(L$nbr_sp)
  ns <- c(1,add(ne,1)[-nb_acq])

  header <- c("names","ID", "nbr_MS", "start", "end", "used", "blank (ID)", "color",
              "concentration","unit","acq_T0 (ID)", "delta_T (s)", "grp1", "grp2", "...")

  mt <- matrix("", nrow = nb_acq, ncol = length(header)-6) %>%
    cbind(L$names, 1:nb_acq, L$nbr_sp, ns, ne, rep(TRUE, nb_acq),.) %>%
    rbind(header,.)

  write.table(mt, file = paste0(L$wd,"/meta_empty.csv"), sep = ";", dec = ",", row.names = FALSE, col.names = FALSE)

  colnames(mt) <- header
  mt <- mt[-1, -1]
  rownames(mt) <- L$names

  mt[, "color"] <- ctrl.color(mt[,"color"])

  L$mt <- list("name" = "import", "meta" = mt)
  return(L)
}

#' Import a meta file
#' @param nm "meta_1" the name (whitout .csv) of a new meta data file.
#' @param L sp
#' @return sp
#' @export
#' @examples
#' # sp <- import.meta("meta_1")
import.meta <- function(nm = "meta_empty", L = sp){

  mt <- read.table(paste0(L$wd,"/",nm,".csv"), sep = ";", dec = ",",
                   header = TRUE, row.names = 1, stringsAsFactors = FALSE,
                   comment.char = "", check.names = FALSE)
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

#### Provoc_02_importation ####
#### Importation ####

#' Read one file h5
#' @param num_fil a number, e.g. 1
#' @param ll the f_h5 obj
#' @param sk skip, the number of non-imported spectra (starting with the first)
#' @return a temporary sp
#' @noRd
read.h5 <- function(num_fil=1, ll = f_h5, wd = wdir, sk = skip){

  # find the name of file
  name_h5 <- nm.ls(num_fil, wd)

  # files import
  act_h5 <- paste0(wd,"/h5/",dir(paste0(wd,"/h5"))[num_fil]) %>% H5Fopen()

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
  hprint(paste0(name_h5, " # ",which(num_fil == ll), "/", length(ll)))

  # return
  rd <- (sk+1):ncol(MS)
  list("name" = name_h5,
       "xMS" = xMS,
       "MS" = MS[,rd],
       "date" = all_date[rd],
       "timing" = all_timing[rd],
       "nbr_sp" = length(rd),
       "meta" = all_TPS2[,rd])
}

#' Import all file.h5 in the h5 folder of working directory.
#'
#' For use this functions, you must have a folder name "h5" whith acquisition inside.
#' @param wdir the working directory
#' @param pk_param a set of parameters for importation and detection peak
#' @param ctrl_peak logical. An option for create a visual detect peak checking plot
#' @param baseline_correction logical. If TRUE, the baseline is correted
#' @param skip, numeric. The number of non-imported spectra (starting with the first)
#' @return sp
#' @export
#' @examples
#' # wd <- "C:/Users/huguenin/Documents/R/provoc test/data test/miscalenous"
#' #
#' # /!\ Note : Your datas are store like that :
#' # "wd/h5/00_file_PTR_ToF_MS.h5"
#' # "wd/h5/01_file_PTR_ToF_MS.h5"
#' # "wd/h5/02_file_PTR_ToF_MS.h5"
#' #
#' # setwd(wd)
#' # sp <- import.h5(wd)
#' #
#' # For the parameters :
#' # pk_param =  c(NULL,"very hight", "hight", "medium", "low")
#' # This determine the sensibility of your peak target.
#' # list(method = "MAD", halfWindowSize = 5, SNR = 10, smooth = 6) #NULL
#' # list(method = "MAD", halfWindowSize = 2, SNR = 10, smooth = 6) # very hight
#' # list(method = "MAD", halfWindowSize = 2, SNR = 40, smooth = 6) # hight
#' # list(method = "MAD", halfWindowSize = 5, SNR = 40, smooth = 6) # medium
#' # list(method = "MAD", halfWindowSize = 10, SNR = 60, smooth = 6) # low
#' #
#' # enougth, you can directly choose parameters :
#' # pk_param = list(method = "MAD", halfWindowSize = 2, SNR = 40, smooth = 6)
#' # method = "MAD" or "SuperSmoother"
#' # halfwindSize, SNR and smooth are integers
import.h5 <- function(wdir = getwd(), pk_param = NULL, ctrl_peak = FALSE, baseline_correction = TRUE, skip = FALSE){
  # wdir <- "D:/Cao Li/spring_2022"

  if(("Figures" %in% dir(wdir))==FALSE){
    dir.create(paste0(wdir,"/Figures"))
    dir.create(paste0(wdir,"/Figures/Control"))
  }

  if(("h5" %in% dir(wdir))==FALSE){
    cat("Sorry but the import can't continue. Create a \"h5\" folder with all .h5 fills that you
        want analyse.")
  }

  # detect peak parameters
  if(is.null(pk_param) == TRUE) pk_param <- list(method = "MAD", halfWindowSize = 5, SNR = 10, smooth = 6)
  if(pk_param[[1]] == "very hight") pk_param <- list(method = "MAD", halfWindowSize = 2, SNR = 10, smooth = 6)
  if(pk_param[[1]] == "hight") pk_param <- list(method = "MAD", halfWindowSize = 2, SNR = 40, smooth = 6)
  if(pk_param[[1]] == "medium") pk_param <- list(method = "MAD", halfWindowSize = 5, SNR = 40, smooth = 6)
  if(pk_param[[1]] == "low") pk_param <- list(method = "MAD", halfWindowSize = 10, SNR = 60, smooth = 6)

  # data importation ####
  f_h5 <- dir(paste0(wdir,"/h5")) %>% grep("\\.h5",.)     # localise h5 files

  length(citation.list) %>% sample(1) %>% citation.list[[.]] %>% cat()
  cat(" \n - - - - - - - - - - - - - - - \n")
  list_h5 <- lapply(f_h5, read.h5, ll = f_h5, wd = wdir, sk = skip) # num_fil = f_h5

  # formating of sp list ####
  sp <- list()
  sp$names <- sapply(list_h5,conc.lst, elem = 1)
  sp$Tinit$date <- sapply(list_h5,conc.lst, elem = 4, simplify = FALSE)
  sp$Tinit$timing <- sapply(list_h5,conc.lst, elem = 5, simplify = FALSE)
  sp$nbr_sp <- sapply(list_h5,conc.lst, elem = 6)
  sp$meta <- sapply(list_h5,conc.lst, elem = 7, simplify = FALSE)

  sp$xMS <- mass.shift(list_h5)
  hprint("Concatene MS")

  sp$MS <- list()
  for(i in 1:length(list_h5)){
    sp$MS <- c(sp$MS, list(list_h5[[i]]$MS[1:length(sp$xMS),]))
    list_h5[[i]] <- 0
  }
  sp$MS <- do.call(cbind,sp$MS) %>% t()
  remove(list_h5)

  # size reduction ####
  hprint("Reduction")
  sp <- red.xMS(sp)

  # baseline correction ####
  if(baseline_correction == TRUE){
    hprint("Baseline correction")
    WM <- diff(sp$xMS) %>% mean() %>% divide_by(0.2,.) %>% ceiling()
    fmr <- baseline.rollingBall(sp$MS, wm = WM, ws = WM)
    sp$MS <- fmr$corrected
    remove(fmr)
  }

  # Mass spectrum objet ####
  hprint("Create MassSpectrum object")
  sp$MS <- apply(sp$MS,1, create_local_MS, xMS = sp$xMS)

  # smooth spectra ####
  hprint("Smooth spectra")
  oldw <- getOption("warn")
  options(warn = -1)
  sp$MS <- smoothIntensity(sp$MS,
                           method = "SavitzkyGolay",
                           halfWindowSize = pk_param$smooth)
  options(warn = oldw)

  # align spectra ####
  hprint("Align spectra")
  sp$names_acq <- prep.names(sp) %>% apply(2,names.samples)
  sp$MS <- alignSpectra(sp$MS, tolerance = 0.02)
  sp$MS <- averageMassSpectra(sp$MS, labels = convertStr2List(sp), method="mean")

  # peak detection ####
  hprint("Peak detection")

  sp$peaks <- detectPeaks(sp$MS,
                          method = pk_param$method,
                          halfWindowSize = pk_param$halfWindowSize,
                          SNR = pk_param$SNR)

  sp$peaks <- binPeaks(sp$peaks, tolerance=0.01)
  sp$peaks <- MALDIquant::filterPeaks(sp$peaks, minFrequency=0.005)
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

  sp <- empty.meta(L = sp)
  sp <- list.order(L = sp)

  hprint("Import is completed")

  if(ctrl_peak == TRUE){
    hprint("Check the peak detection")
    fmr <- colnames(sp$peaks) %>% as.numeric() %>% round(0) %>% unique()
    sapply(fmr, peak.ctrl, L=sp)
  }else if(ctrl_peak == FALSE){
    cat("\n Warning ! You don't control the quality of peak detection. You can do it whith :
        \n sapply(c(50:250), peak.ctrl)
        \n or
        \n peak.ctrl(137) \n")
  }

  return(sp)
  # import function is finished ####
}

#' Reduction of abscissa
#' @param L sp
#' @return a temporary sp
#' @noRd
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
    # tiff(paste("Densite des masses maximum.tiff"))
    #  plot(dMS, main="Density",xlim = c(0, thr*2))
    #  abline(v = thr, lty = 2, lwd = 2, col = "red")
    #  legend("topright", bty = "n",
    #         legend = c(paste("threshold = ",round(thr,0)),
    #                    paste("nbr of mass deleted =", length(ind_del))))
    # dev.off()
    #
    # plot.threshold <- function(br, L= sp, z = c(0,300), ind_d = indinf, ind_k = indinf[kp],inu = inull){
    #  tiff(paste("Vue des masses supprimees de",br[1],"a",br[2],"Da.tiff"))
    #    a <- det.c(br[1],L$xMS):det.c(br[2],L$xMS)
    #    matplot(sp$xMS[a],maxMS[a], type = "l", ylim = z,
    #            xlab = "m/z (Da)", ylab = "intensite (u.a.)",
    #            main = "Spectre de l'intensite maximale de chaque masse")
    #    legend("topleft", bty = "n", lty = 1, col = c("black","blue","red"),
    #           legend = c("spectre max","masses supprimees", "masses gardees"))
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

#### Provoc_03_export_figures ####
#### Plot spectra ####

#' Plot spectra dynamic.
#'
#' Export a html plot for an intuitive exploration of data. The number of spectra is limited at 60
#' for increase the comfort of viewing.
#' @param sel_sp a numeric vector with the selection of spectra.
#' @param L sp
#' @param new_color Logical TRUE/FALSE
#' @param name the title of plot and the name of the file.
#' @return a html plot
#' @export
#' @examples
#' # dy.spectra(sp$mt$meta[sp$acq,"end"], new_color = FALSE)
#' # dy.spectra(1:50)
#' # dy.spectra(c(1, 4, 8, 22, 30), new_color = TRUE)
dy.spectra <- function(sel_sp = sp$mt$meta[sp$acq,"end"], L = sp, new_color = FALSE, name = FALSE){
  if(is.character(sel_sp) == TRUE) sel_sp <- as.numeric(sel_sp)
  if(length(sel_sp) > 60) print("Caution /!\ The number of spectra is too big. Select less spectra.")
  if(length(sel_sp) <=60){

    sp_sel <- L$MS[,sel_sp]

    if(new_color == TRUE) dy_color <- viridis(length(sel_sp),alpha = 0.8)
    if(new_color == FALSE) dy_color <- rep.mtm("color", L, sel = "all")[sel_sp]

    if(length(sel_sp)==1){
      dysp <- sp_sel %>% cbind(L$xMS,.) %>% as.data.frame()
      colnames(dysp)[2] <- colnames(L$MS)[sel_sp]
      titre <- colnames(L$MS)[sel_sp]
      if(name != FALSE) titre <- name
    }else{
      dysp <- cbind(L$xMS,sp_sel) %>% as.data.frame()
      titre <- "sp_align"
      if(name != FALSE) titre <- name
    }
    rownames(dysp) <- L$xMS

    ftitre <- paste0(L$wd,"/Figures/") %>% dir() %>% grep(titre, .) %>% length() %>% add(1)
    ftitre <- str_pad(ftitre,width = 2,pad = "0")
    ftitre <- paste0(L$wd,"/Figures/dy_",titre,"_",ftitre)
    fmr <- system.file("rmd", "print_dy_sp.Rmd", package = "provoc")
    rmarkdown::render(input = fmr, output_file = ftitre)
  }
}

#' Plot spectra tiff
#'
#' Export a tiff plot for fixed idea, a presentation or an article.
#' @param sel_sp a numeric vector with the selection of spectra.
#' @param pkm the lower mass
#' @param pkM the upper mass
#' @param L sp
#' @param new_title an explicite title
#' @param new_color Logical TRUE/FALSE
#' @param leg legend on the rigth "r" or on the left "l"
#' @return a tiff plot
#' @export
#' @examples
#' # For a large spectre :
#' # fx.spectra(sel_sp = sp$mt$meta[sp$acq,"end"], pkm = 137, pkM = 137, leg = "l")
#' #
#' # for just one peak :
#' # fx.spectra(seq(1,100, by = 10), pkm = 59, pkM = 150)
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

  ntitre <- paste0(L$wd,"/Figures/") %>% dir() %>% grep(new_title, .) %>% length() %>% add(1)
  ntitre <- paste0(L$wd,"/Figures/",new_title,"_",ntitre,"_zoom_",xmin,"_to_",xmax,".tif")

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

#' a unique function for monitoring the kinectic of COV.
#'
#' This function allows to plot the intensity kinetics of several peaks.
#'  The graphs produced can group or not spectra belonging to the same class.
#'  It is also possible to observe a single peak or a group.
#' @param M_num a vector with exact masses. It's possible to specifique these masses like
#' a vector c(59.045, 137.121) or to be more evasive M.Z.max(c(59, 137)).
#' @param each_mass Logical TRUE of FALSE. If TRUE, a unique plot with all mass specifie in M_num.
#' Else if FALSE, a plot is create for each mass.
#' @param group FALSE or the name of a column of meta table. e.g. "grp1".
#' @param graph_type "dy" or " fx" for create a dynamic plot in html or a fixed plot in tiff.
#' @param L sp
#' @param Y_exp Logical TRUE or FALSE. The y axe is exponential ?
#' @param time_format "date" for abscissa in day hour minutes seconde with the real date of the acquisition spectra.
#' Or "time" for combine kinetic of several acquisition.
#' @return a plot
#' @export
kinetic.plot <- function(M_num = M.Z.max(c(59, 137)), each_mass = TRUE,
                         group = FALSE, graph_type = "dy", L = sp,
                         Y_exp = FALSE, time_format = "date"){
  # Mise en forme :
  tit_wd <- paste0(L$wd,"/Figures/kinetic")
  if(("kinetic" %in% dir(paste0(L$wd,"/Figures")))==FALSE) dir.create(tit_wd)
  vp <- list(exp = Y_exp, time = time_format, grp = group)

  tit_wd <- paste0(tit_wd,"/pk_at")

  # Graphe :

  # Etape 1 : chaque mass ?
  if(each_mass == TRUE){
    # OUI

    for(ma in M_num){ # ma = M_num[1]

      # Etape 2 : chaque groupe ?
      if(group != FALSE){
        # OUI
        grp <- L$mt$meta[,group][L$acq] %>% as.character() # on repere les groupes
        for(u in unique(grp)){ # u = unique(grp)[1]

          save_acq <- L$acq
          L$acq <- which(u == L$mt$meta[,group]) %>% intersect(L$acq)

          # Etape 3 : plot statique ou dynamique ?
          if(graph_type == "fx"){
            # plot fixe

            titre <- c(tit_wd, ma,"of",u) %>% str_flatten(" ") %>% paste0("_",vp$time,".tiff")
            fx.kinetic.plot(L, titre, acq = L$acq, MA = ma, VP = vp)

          }else if(graph_type == "dy"){
            # plot dynamique
            titre <- str_flatten(ma, collapse = " ") %>% paste("peak at",.,"of", u, vp$time)
            dy.kinetic.plot(L, titre, acq = L$acq, MA = ma, VP = vp)
          }
          L$acq <- save_acq
        }

      }else{
        # NON

        # Etape 3 : plot dynamique ou statique ?
        if(graph_type == "fx"){
          # plot statique

          titre <- c(tit_wd, ma,"of",head(row.names(L$mt$meta)[L$acq])) %>%
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

        save_acq <- L$acq
        L$acq <- which(u == L$mt$meta[,group]) %>% intersect(L$acq) # on repere les indices de chaques groupes

        # Etape 3 : plot statique ou dynamique ?
        if(graph_type == "fx"){
          # plot fixe

          titre <- c(tit_wd, ma,"of",u) %>%
            str_flatten(" ") %>% paste0("_",vp$time,".tiff")
          fx.kinetic.plot(L, titre, acq = L$acq, MA = ma, VP = vp)

        }else if(graph_type == "dy"){
          # plot dynamique
          titre <- str_flatten(ma, collapse = " ") %>% paste("peak at",.,"of", u, vp$time)
          dy.kinetic.plot(L, titre, acq = L$acq, MA = ma, VP = vp)
        }
        L$acq <- save_acq
      }
    }else{
      # NON

      # Etape 3 : plot statique ou dynamique ?
      if(graph_type == "fx"){
        # plot statique

        titre <- c(tit_wd, ma,"of",head(row.names(L$mt$meta)[L$acq])) %>%
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


#' a internal fonction. Use kinetic.plot()
#' @param L sp
#' @param titre a string of character
#' @param acq a group of spectra
#' @param MA a group of masse
#' @param VP a list of other data
#' @return a plot
#' @noRd
fx.kinetic.plot <- function(L, titre, acq = ind_PK, MA = ma, VP = vp){
    # index for peaks and acquisitions
    ind_pk <- ind.pk(MA,"exact",L)
    ind_Sacq <- ind.acq(acq,L)

    # intensity
    Ibn <- c(0, max(L$peaks[ind_Sacq, ind_pk]))

    # time
    if(VP$time == "time"){
      Tbn <- c(0,0)
      sapply(acq, function(X, Li) max(Li$Trecalc$timing[[X]]), Li = L) %>%
        max() -> Tbn[2]

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
        xx <- unlist(List_abs)[coor]
        if(length(coor)>1){
          fmr <- length(coor) + (1-nq)*round(length(coor)/30)
          matplot(xx, L$peaks[coor,j], type = "l", lwd = 2,
                  col = pk_col[i], add = TRUE)
          matplot(xx, L$peaks[coor,j],
                  pch = cl, col = pk_col[i], cex = 2, add = TRUE)
        }else{
          matplot(xx, L$peaks[coor,j],
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

#' a internal fonction. Use kinetic.plot()
#' @param L sp
#' @param titre a string of character
#' @param acq a group of spectra
#' @param MA a group of masse
#' @param VP a list of other data
#' @return a plot
#' @noRd
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
    l_acq <- apply(L$mt$meta, 1, function(mat) as.numeric(mat["start"]):as.numeric(mat["end"]))
    List_abs <- lapply(l_acq, function(liste, Tlist) unlist(Tlist)[liste], Tlist = List_abs)
  }

  # or date
  if(VP$time == "date"){
    Cdate <- as.POSIXct(Sys.time(), format="%m/%d/%Y %H:%M:%S")
    attr(Cdate, "tzone") <- "UTC"
    for(i in 1:nrow(L$mt$meta)) Cdate <- c(Cdate, L$Trecalc$date[[i]])
    Cdate <- Cdate[-1]
    Tbn <- range(Cdate)

    Xlab <- "Date"

    List_abs <- sapply(acq, function(acq) Cdate[ind.acq(acq,L)])
  }

  if(VP$exp == FALSE) VP$exp <- ""
  if(VP$exp == TRUE){
    VP$exp <- "y"
    Ibn[1] <- 10
  }

  # the list of data
  dy_mat <- sapply(acq, dy.mat.pk, ipk = ind_pk, La = List_abs, Li = L,
                   vp = VP, simplify = FALSE)

  if(VP$grp != FALSE) dy_mat <- lapply(dy_mat, dy.trans.ID)

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
  dysp$ID <- NULL
  dycol <- unname(L$mt$meta[acq,"color"]) %>% rep(each = length(MA))

  # and plot
  if(VP$time == "date"){
    dysp$xT <- as.POSIXct(Cdate[ind_Sacq], origin = "1970-01-01", tz = "")
    dysp <- xts::xts(dysp[,-1], order.by = dysp$xT, tz = "")

    fmr <- system.file("rmd", "print_dypk_date.Rmd", package = "provoc")
    rmarkdown::render(input = fmr,
                      output_file = paste0(L$wd, "/Figures/kinetic/", titre))
  }else{
    fmr <- system.file("rmd", "print_dypk_time.Rmd", package = "provoc")
    rmarkdown::render(input = fmr,
                      output_file = paste0(L$wd, "/Figures/kinetic/", titre))
  }
}

#' a internal fonction. Use kinetic.plot()
#' @param ac a group of spectra
#' @param ipk a group of mass
#' @param La a list of peaks
#' @param Li sp
#' @param vp a list of other data
#' @return a matrix
#' @noRd
dy.mat.pk <- function(ac = acq, ipk = ind_pk, La = List_abs, Li = L, vp = VP){

  eq_acq <- which(Li$acq == ac)

  # if(vp$grp == FALSE) nid <- rownames(Li$mt$meta)[ac]
  # if(vp$grp != FALSE) nid <- Li$mt$meta[ac,vp$grp]
  nid <- rownames(Li$mt$meta)[ac]

  fmr <- rep(Li$mt$meta[ac,"ID"], length(La[[eq_acq]])) %>% as.numeric()
  fmr <- cbind(La[[eq_acq]],fmr, Li$peaks[ind.acq(ac,Li), ipk])
  colnames(fmr) <- c("xT","ID", paste(nid, colnames(Li$peaks)[ipk]))

  return(fmr)
}


#' a internal fonction. Use kinetic.plot()
#' @param li a list
#' @return another list
#' @noRd
dy.trans.ID <- function(li){
  colnames(li)[-(1:2)] <- paste(colnames(li)[-(1:2)],li[1,2])
  return(li)
}

#### Control Peak et metadata ####

#' print a graph with MS and peaks detected.
#'
#' The graph will be centered on peak. It is not necessary to be precise for this number.
#' @param peak a number corresponding to a mass to be checked
#' @param L sp
#' @param suffixe a precision about this peak or this serie
#' @return a plot
#' @export
#' @examples
#' # for(i in 50:210) peak.ctrl(i)
peak.ctrl <- function(peak = 137, L = sp, suffixe = ""){

  pk_titre <- paste0("Peak_ctrl",suffixe)
  if(pk_titre %in% dir(paste0(L$wd,"/Figures/"))==FALSE){
    dir.create(paste0(L$wd,"/Figures/",pk_titre))
  }

  zm <- c(peak-.3,peak+.3)
  brm <- det.c(zm[1],L$xMS)
  brM <- det.c(zm[2],L$xMS)

  ntitre <- paste0("Figures/",pk_titre,"/MS_vs_peak_at_",peak,"_Da",suffixe,".tiff")
  pk <- colnames(L$peaks) %>% as.numeric()

  if(length(L$xMS[brm:brM]) >1){
    tiff(filename = ntitre, width = 1000, height = 580)
    par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

    matplot(L$xMS[brm:brM], L$MS[brm:brM,], type = "l",
            col = viridis(ncol(L$MS), alpha= 0.5),
            main = paste("MS (green) vs peaks (orange) at",peak,"Da",suffixe),
            xlab = "m/z (Da)", ylab = "intensity (a.u.)")
    matplot(pk, t(L$peaks), pch = 16, col = alpha("darkorange3",0.5), add = TRUE)
    dev.off()
  }
}

#' print a lot of graph for a rapid overview of metadata kinetic
#'
#' @param L sp
#' @param short a short selection of meta. FALSE for all graph
#' @return a plot
#' @export
#' @examples
#' # meta.ctrl()
meta.ctrl <- function(L = sp, short = TRUE){
  mtdt <- do.call(cbind, L$meta)
  max_mtdt <- apply(mtdt,1,max)
  mtdate <- do.call(c,L$Tinit$date)

  sel <- c(1,9,15,55,75,76,77,78,79,87,88,91,95,96,128,163,167,171,182,205,213,214)
  if(short == FALSE) sel <- which(max_mtdt != 0)

  for (nm in sel){ # nm = 75
    tiff(filename = paste0("Figures/Control/meta_data_",nm,"_",rownames(mtdt)[nm],".tiff"), width = 500, height = 300)
    par(mar = c(2.5,2.5,2,0.5), mgp = c(4,1,0))
    matplot(mtdate, mtdt[nm,], pch = 16, main = rownames(mtdt)[nm])
    dev.off()
  }
}

#### Provoc_04_others_functions ####
#### micro-functions ####

#' print garbage collection
#' @return
#' @noRd
print.gc <- function(){
  fmr <- memory.size()
  gc()
  paste0("RAM : ",fmr," -> gc -> ", memory.size()) %>% hprint()
}

#' concatenation
#' @param list_n a list
#' @param elem the index
#' @return just one element of the list
#' @noRd
conc.lst <- function(list_n, elem = 1){
  list_n[[elem]]
}

#' concatenation (dimension)
#' @param list_n a list
#' @param elem the index
#' @return the dimension
#' @noRd
dim.lst <- function(list_n, elem = 1){
  dim(list_n[[elem]])
}

#' names of acquisition
#' @param L sp
#' @return a name
#' @noRd
prep.names <- function(L){
  fmr <- log10(L$nbr_sp) %>% floor() %>% add(1)
  rbind(fmr, L$names, L$nbr_sp)
}

#' a internal function for preparation name
#' @param vec a section of meta table
#' @return a string character
#' @noRd
names.samples <- function(vec){str_pad(1:vec[3],vec[1], pad = "0") %>% paste(vec[2],.,sep = "_")}

#' convert string to list
#' @param L sp
#' @return a list
#' @noRd
convertStr2List <- function(L){
  plip <- function(vec) return(vec)
  fmr <- unlist(L$names_acq) %>% lapply(plip)
  names(fmr) <- unlist(L$names_acq)
  return(fmr)
}

#' control color for replace NA by a color
#' @param vec_col a string character with NA or color.
#' @return vec_col
#' @noRd
ctrl.color <- function(vec_col = mt[,"color"]){

  fmr <- which.na(vec_col == "")
  if(length(fmr) > 0) vec_col[fmr] <- viridis(length(fmr)) %>% alpha(0.5)

  fmr <- which(vec_col == "")
  if(length(fmr) > 0) vec_col[fmr] <- viridis(length(fmr)) %>% alpha(0.5)

  return(vec_col)
}

# ' create Mass Spectrum objet
#' @param MS a matrix
#' @param xMS a vector
#' @return a MassSpectrum objet
#' @noRd
create_local_MS <- function(MS, xMS){createMassSpectrum(xMS,MS)}

#' return to spectra
#' @param spobj a MassSpectrum obj
#' @return a vector of intensity
#' @noRd
mat.spectra <- function(spobj){spobj@intensity}

#' return to spectra
#' @param spobj a MassSpectrum obj
#' @return a vector of mass
#' @noRd
mass.spectra <- function(spobj){spobj@mass}

#' return the index more closed of brn in the numeric vector.
#' @param brn a number
#' @param vec a numeric vector
#' @return a number
#' @export
#' @examples
#' # vec <- seq(1,100,by = .2)
#' # a <- det.c(59.38, vec)
#' # a = 293. So vec[293] is the most closed of 59,38 than vec.
det.c <- function(brn,vec){
  subtract(vec,brn) %>% sapply(abs) %>% which.min()
}

#' detect length of x mass
#' @param splist sp
#' @return a number
#' @noRd
length.xMS <- function(splist){length(splist$xMS)}

#' print the text below by the hour
#' @param txt a character string
#' @return a character string more the hour
#' @export
#' @examples
#' # hprint(txt = "hello there")
hprint <- function(txt = "hello there"){
  heure() %>% paste0(txt,", ",.) %>% message()
}

#' return the closed peak
#' @param pk_x a peak
#' @param mat a matrix of peak intensity
#' @param w.sub the width of window
#' @return a vector of peaks closed
#' @noRd
pk.red <- function(pk_x = pk_max[1,1], mat = pk_max, w.sub = 4){
  fmr <- subtract(mat[1,],pk_x) %>% abs() %>% multiply_by(-1) %>% which.sup(-(w.sub+1))
  return(fmr[which.max(mat[2,fmr])])
}

#' list of principals peaks
#' @param pk_mat matrix of peaks
#' @return a matrix with principal peak
#' @noRd
pk.short <- function(pk_mat = L$peaks){
  pk_max <- colnames(pk_mat) %>% as.numeric() %>% rbind(apply(pk_mat,2,max))
  fmr <- sapply(pk_max[1,], pk.red, mat = pk_max, w.sub = 4) %>% unique() %>% sort()
  pk_mat <- pk_max[,fmr]
  rownames(pk_mat) <- c("m/z","int")
  return(pk_mat)
}

#' give the time in the format that I want
#' @return an hour
#' @noRd
heure <- function(){str_split(Sys.time(),pattern = " ")[[1]][2]}

#' order the list
#' @param L sp
#' @return sp
#' @noRd
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

#' named workflow
#' @param nwf I forget
#' @param L sp
#' @return sp
#' @noRd
name.wf <- function(nwf = "randow", L = sp){
  fmr <- length(L$workflow)
  names(L$workflow)[[fmr]] <- nwf
  return(L)
}

#' update workflow
#' @param nm_wf name of operation
#' @param obj_wf param of that
#' @param L sp
#' @return sp
#' @noRd
wf.update <- function(nm_wf, obj_wf, L = sp){
  L$workflow <- c(L$workflow, list(obj_wf))
  L <- name.wf(nm_wf, L)
  return(L)
}

#' Repet meta parameter
#' @param col.nam a name of column
#' @param L sp
#' @param sel "acq" or "all"
#' @return a matrix
#' @noRd
rep.mtm <- function(col.nam, L, sel = "acq"){
  fmr <- L$acq
  if(sel == "all") fmr <- as.numeric(L$mt$meta[,"ID"])
  sapply(fmr, rep.mtu, col.nam = col.nam, L = L, simplify = FALSE) %>% unlist()
}

#' Repet meta parameter of a single aquisition
#' @param acq a number
#' @param col.nam a name of column
#' @param L sp
#' @return a vector
#' @noRd
rep.mtu <- function(acq, col.nam, L){
  fmr <- which(col.nam == colnames(L$mt$meta))
  rep(L$mt$meta[acq,fmr], L$nbr_sp[acq])
}

#' round to the lower ten
#' @param x a number
#' @return a number
#' @noRd
dizaine <- function(x){
  eph <- log(x,10) %>% floor() %>% multiply_by(10)
  divide_by(x,eph) %>% floor() %>% multiply_by(eph)
}

#' search all peak in accord to a mass number
#' @param ma a number of mass
#' @param L sp
#' @return a numeric vector with all the exact mass closed to the mass number
#' @export
#' @examples
#' # M.Z(c(59, 137))
#' # 58.873  59.045  59.233  59.267  59.320  59.405 137.037 137.121
M.Z <- function(ma,L=sp){
  vec_pk <- colnames(L$peaks) %>% as.numeric()
  fmr <- NULL
  for(maz in ma) fmr <- c(fmr, which((vec_pk < maz + 0.5)&(vec_pk > maz - 0.5)))
  return(vec_pk[fmr])
}

#' search the highest peak in accord to a mass number
#' @param ma a number of mass
#' @param L sp
#' @return a numeric vector with all the exact mass closed to the mass number
#' @export
#' @examples
#' # M.Z.max(c(59, 137))
#' # 59.233 137.121
M.Z.max <- function(ma, L = sp){
  j <- 0
  for(i in ma){
    j <- j + 1
    fmr <- match(M.Z(i),colnames(L$peaks))
    if(length(fmr) == 0){
      ma[j] <- NA
    }else if(length(fmr)==1){
      ma[j] <- colnames(sp$peaks)[fmr] %>% as.numeric()
    }else{
      mInd <- apply(sp$peaks[,fmr],2, max) %>% which.max()
      ma[j] <- colnames(sp$peaks)[fmr[mInd]] %>% as.numeric()
    }
  }
  return(ma)
}

#' return index of spectra for each acquistion
#' @param n_acq a number
#' @param L sp
#' @return a numeric vector
#' @export
#' @examples
#' # For just one acquisition :
#' # ind.acq(1,sp)
#' # Or for more :
#' # sapply(1:10,ind.acq,L=sp,simplify = TRUE)
ind.acq <- function(n_acq,L){
  fmr <- NULL
  mat_mt <- cbind(as.numeric(L$mt$meta[,"start"]),
                  as.numeric(L$mt$meta[,"end"]))
  for(i in n_acq) fmr <- c(fmr, mat_mt[i,1]:mat_mt[i,2])
  return(fmr)
}

#' return index of peaks
#' @param ms the mass of the peak
#' @param mode a character. Specifies whether the returned index is the exact peak (exact), the maximum peak (max) or all peaks (all) related to a unit mass.
#' @param L sp
#' @return a numeric vector
#' @export
#' @examples
#' # For get the index for one peak :
#' # ind.pk(137)
#' # [1] 318
#' # Or for more :
#' # sapply(c(81,137,205),ind.pk)
ind.pk <- function(ms, mode = c("max","exact","all"), L = sp){
  ms <- as.numeric(ms)
  if(mode[1] == "all") ind <- match(M.Z(ms),colnames(L$peaks))
  if(mode[1] == "max") ind <- match(M.Z.max(ms),colnames(L$peaks))
  if(mode[1] == "exact") ind <- match(ms,colnames(L$peaks))
  if((mode[1] == "exact")&(is.na(ind[1])==TRUE)) print("Warning! There is no peak for this mass.")
  return(ind)
}


#### which pack ####

#' return the index of elements egals at nb.
#' @param vec a numeric vector
#' @param nb a number
#' @return index
#' @noRd
which.equal <- function(vec,nb){which(vec == nb)}

#' return the index of NA.
#' @param x a vector
#' @return index
#' @noRd
which.na <- function(x){which(is.na(x) == TRUE)}

#' return the index of not NA elements.
#' @param x a vector
#' @return index
#' @noRd
which.not.na <- function(x){which(is.na(x) == FALSE)}

#' return the index of elements superior of theshold.
#' @param vec a numeric vector
#' @param threshold a number
#' @return index
#' @noRd
which.sup <- function(vec, threshold){return(which(vec > threshold))}

#### Provoc_05_modify H5 ####
#### Manip h5 files ####


#' extrat info of h5 file
#' @param num a numeric vector
#' @param w_dir working directory
#' @return a list with informations
#' @noRd
export.info <- function(num = 1 , w_dir = w_d){

  # files import
  act_h5 <- dir(paste0(w_dir,"/h5"))[num] %>% paste0(w_dir,"/h5/",.) %>% H5Fopen()
  all_MS <- act_h5$FullSpectra$TofData[,1,,]
  all_timing <- act_h5$TimingData$BufTimes
  H5Fclose(act_h5)

  # title
  nm_h5 <- dir(paste0(w_dir,"/h5"))[num]
  nm_h5 <- str_remove_all(nm_h5,paste0(w_dir,"/h5/"))
  nm_h5 <- str_remove_all(nm_h5,"_20......_......")
  nm_h5 <- str_remove_all(nm_h5,"20......_......_")
  nm_h5 <- str_remove_all(nm_h5,"\\.h5")

  # intensities extraction [ acquisition number * 160 000 pts]
  fmr <- dim(all_MS)
  dim(all_MS) <- c(fmr[1], fmr[2] * fmr[3]) # for 2d array

  # timing extraction
  dim(all_timing) <- c(fmr[2] * fmr[3])
  timing <-  diff(all_timing) %>% mean() %>% round(2)

  pres <- c(timing,dim(all_MS)) %>% as.matrix()
  colnames(pres) <- nm_h5
  rownames(pres) <- c("acq time (s)", "length abs", "nbr spectra")

  # combine
  mat <- rbind(apply(all_MS,2,mean),
               apply(all_MS,2,median),
               apply(all_MS,2,max))
  rownames(mat) <- c("mean", " median", "max")

  fmr <- dim(all_MS)[2] %>% dizaine() %>% log(10) %>% round() %>% add(1)
  fmr <- str_pad(1:dim(all_MS)[2], width = fmr, pad = "0")
  colnames(mat) <- paste0(nm_h5, "_", fmr)

  # rep <- list("name" = nm_h5, "pres" = c(timing,dim(all_MS)), "mat" =  mat)
  rep <- list("pres" = pres, "mat" =  mat)
  return(rep)
}

#' summarise the h5 files
#'
#' the function return a list of two objects, one matrix and one list. The summary show all acquisitions in column
#' and others informations in row. Theses informations are the ID, the acquisition time (in second),
#' the length of x mass, the number of spectra, index of start and end, and
#' the smallest value among the maximum intensities of each spectrum of the acquisition.
#'
#' the other list is order by ID. Each table includes the average, median and
#' maximum intensities of each spectrum of the acquisition.
#' @param w_d working directory
#' @return a list with informations
#' @export
#' @examples
#' # l_info <- info.h5()
#' # l_info$summary
#'
#' # A2_20  A2_30  A2_40  A3_20  A3_30  A3_40
#' # ID                1      2      3      4      5      6
#' # acq time (s)     10     10     10     10     10      0
#' # length abs   158768 158768 158768 158768 158768 158768
#' # nbr spectra     180    180    180    180    180    176
#' # start             1    181    361    541    721    901
#' # end             180    360    540    720    900   1076
#' # min Imax     184083 112189 251471  47338 137394      0
#'
#' # We can view than the A3_40, ID 6, have an intensity at zero. And more, all
#' # acquisitions are 180 spectra and acquisition during 10 seconds except the ID 6.
#'
#' # l_info$$byID[[6]]
#' #         A3_40_174    A3_40_175     A3_40_176
#' # mean    407.22151    404.94106             0
#' # median  13.36089     13.32512              0
#' # max     792628.75    789557.87             0
#'
#' # The last spectra is corrompted. We should deleted the number 176.
#'
info.h5 <- function(w_d = getwd()){
  all_fil <- grep("\\.h5", dir(paste0(w_d,"/h5")))

  oldw <- getOption("warn")
  options(warn = -1)
  l_info <- lapply(all_fil, export.info, w_dir = w_d)
  options(warn = oldw)

  l_mat <-  lapply(l_info, function(liste) liste$mat) %>% as.data.frame()
  l_pres <- lapply(l_info, function(liste) liste$pres) %>% as.data.frame()
  fmr <- lapply(colnames(l_pres), grep, x = colnames(l_mat)) %>% sapply(range)
  colnames(fmr) <- colnames(l_pres)
  l_pres <- rbind(1:ncol(l_pres),l_pres,fmr,rep(0,ncol(l_pres)))
  for(i in 1:ncol(l_pres)) l_pres[7,i] <- round(min(l_mat[3,l_pres[5,i]:l_pres[6,i]]))
  rownames(l_pres)[c(1,5,6,7)] <- c("ID","start","end","min Imax")

  l_mat <-  lapply(l_info, function(liste) liste$mat)

  l_info <- list("byID" = l_mat, "summary" = l_pres)
  return(l_info)
}

#' delete a couple of specrta in a h5 file
#'
#' sometimes, an experience is stopped brutaly. The last spectra is corrupted.
#' The file can be reading by the import.h5() function. info.h5 and delete.spectra.h5()
#' can be used for resolving the problem.
#' @param ID_h5 the ID give by info.h5()
#' @param num the number correspoding at the spectra corrupted
#' @param w_d working directory
#' @return a modifiy h5 file.
#' @export
#' @examples
#' # Use the info.h5() function before use that. And read the help for understand why
#' # we deleted the spectra number 176 of ID6.
#'
#' # delete.spectra.h5(6,176)
#'
#' # After that, you have a new h5 file (mod_A3_40.h5). You can exhile the original file.
delete.spectra.h5 <- function(ID_h5 = 1, num = 1, w_d = getwd()){

  #mod option
  oldw <- getOption("warn")
  options(warn = -1)

  # File opened
  act_h5 <- dir(paste0(w_d,"/h5"))[ID_h5] %>% paste0(w_d,"/h5/",.) %>% H5Fopen()

  # copy in temporary object
  H5 <- list()
  H5$FullSpectra <- act_h5$FullSpectra
  H5$FullSpectra$SaturationWarning <- NULL
  H5$FullSpectra$SumSpectrum <- NULL
  H5$TimingData$BufTimes <- act_h5$TimingData$BufTimes
  H5$TPS2$TwData <- act_h5$TPS2$TwData
  H5$TPS2$TwInfo <- act_h5$TPS2$TwInfo
  H5$AcquisitionLog <- act_h5$AcquisitionLog

  # File closed
  H5Fclose(act_h5)

  # modified H5
  fmr <- dim(H5$FullSpectra$TofData[,1,,])
  d3_num <- divide_by(num, fmr[2]) %>% floor()
  H5$FullSpectra$TofData <- H5$FullSpectra$TofData[,1,,-d3_num]
  dim(H5$FullSpectra$TofData) <- c(dim(H5$FullSpectra$TofData)[1],1,dim(H5$FullSpectra$TofData)[2:3])
  H5$TimingData$BufTimes <- H5$TimingData$BufTimes[,-d3_num]
  H5$TPS2$TwData <- H5$TPS2$TwData[,,-d3_num]

  #create .h5 file
  nom <- dir(paste0(w_d,"/h5"))[ID_h5] %>% paste0(w_d,"/h5/mod_",.)
  h5createFile(nom)
  h5write(H5$FullSpectra, nom, "FullSpectra")
  h5write(H5$TPS2, nom, "TPS2")
  h5write(H5$TimingData, nom, "TimingData")
  h5write(H5$AcquisitionLog, nom, "AcquisitionLog")
  h5closeAll()

  # option re-init
  options(warn = oldw)
}


#### Provoc_06_MCR ####

#' prepar a matrix for the preprocess step
#' @param indM an integer
#' @param selec a vector of selected peaks, or "all"
#' @param s_T a charater string "date" or "time"
#' @param L the list with spectra (sp)
#'
#' @return conform matrix organised in a list
#' @noRd
sort_mat <- function(indM = 14, selec = "all", s_T = "date", L){
  if(selec[1] == "all") selec <- 1:ncol(L$peaks)
  Xmat <- L$peaks[ind.acq(indM,L),selec] # list of select peak
  if(s_T == "date") row.names(Xmat) <- unlist(L$Trecalc$date)[ind.acq(indM,L)] %>% round(0)
  if(s_T == "time") row.names(Xmat) <- unlist(L$Trecalc$timing)[ind.acq(indM,L)] %>% round(0)
  return(Xmat)
}

#' gestion of metadata organised in fonction of the specified group
#' @param name_col a character string of the group 's column name
#' @param L the list with spectra (sp)
#'
#' @return a list with some information
#' @noRd
gest.gr.mcr <- function(name_col = "samples", L = Li){
  vec_ech <- factor(L$mt$meta[,name_col])
  u_ech <- levels(vec_ech)
  list_ech <- lapply(u_ech, grep, vec_ech)
  levels(vec_ech) <- length(u_ech) %>% viridis(begin = 0.15, end = 0.85, option = "turbo")
  return(list("name_gr" = name_col,
              "nbr" = length(u_ech),
              "grp" = u_ech,
              "col" = as.character(vec_ech),
              "lvl" = levels(vec_ech),
              "list" = list_ech))
}

#' for print three differents figures relating to MCR
#' @param col_mt a character string of the groupe
#' @param nc an integer for the composant number
#' @param Lg the list with spectra (sp)
#' @param s_T the time format (date or time)
#' @param res_als a list with the result of MCR
#'
#' @return some figures
#' @noRd
mcr.print <- function(col_mt = "samples", nc = ncMCR, Lg = Li, s_T = "date", res_als = tea.als){

  # la gestion de la couleur
  gr_mcr <- gest.gr.mcr(name_col = col_mt, Lg)

  col_nc <- viridis(nc, begin = 0.2, option = "plasma")

  # l'amplitude temporelle des datas :
  if(s_T == "date") x_zoom <- lapply(res_als$CList, rownames) %>% range() %>% as.numeric() %>% as.POSIXct(origin = "1970-01-01", tz = "GMT")
  if(s_T == "time") x_zoom <- lapply(res_als$CList, rownames) %>% unlist() %>% as.numeric() %>% range()

  # la gestion du repertoire

  fol <- paste0("MCR_",gr_mcr$name_gr,"_nc",nc)

  if(fol %in% dir(paste0(Lg$wd,"/Figures/"))==FALSE){
    dir.create(paste0(Lg$wd,"/Figures/",fol))
  }

  # les graphes

  for(comp in 1:nc){ # comp = 1

    ####
    # figure des spectres purs
    ####

    tiff(filename = paste0("Figures/",fol,"/MCR_spectre_comp",comp,"_on_",nc,".tiff"), width = 1000, height = 580)
    par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

    xpMS <- rownames(res_als$S) %>% as.numeric() # le premier echantillon
    matplot(xpMS, res_als$S[,comp], type = "l", col = col_nc[comp],
            main = paste("MCR spectrum n", comp, "(model with",nc,"cp)"),xlim = c(50,210),
            xlab = "Mass (m/z)", ylab = "Intensity (a.u.)")
    dev.off()

    # le range d'intensite
    y_zoom <- lapply(res_als$CList, function(mat,cp) range(mat[-1,cp]), cp = comp) %>% range()

    ####
    # figure des concentrations purs
    ####

    tiff(filename = paste0("Figures/",fol,"/MCR_score_comp",comp,"_on_",nc,".tiff"), width = 1000, height = 580)
    par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

    xt <- as.numeric(rownames(res_als$CList[[1]])) %>% as.POSIXct(origin = "1970-01-01", tz = "GMT")

    matplot(xt , res_als$CList[[1]][,comp], pch = 16, col = gr_mcr$lvl[1],
            ylim = y_zoom, main = paste("MCR score n", comp, "(model with",nc,"cp)"), xlim = x_zoom,
            xlab = "Time (s)", ylab = "Intensity (a.u.)")

    for(i in 1:length(res_als$CList)){
      xt <- as.numeric(rownames(res_als$CList[[i]]))
      matplot(xt, res_als$CList[[i]][,comp], pch = 16, col = gr_mcr$col[i], add= TRUE)
      matplot(xt , res_als$CList[[i]][,comp], lty = 1, type = "l", col = gr_mcr$col[i], add = TRUE)
    }
    legend("topright", legend = gr_mcr$grp, bty = "n", pch = 16, col = gr_mcr$lvl)
    dev.off()
  }

  for(f in 1:gr_mcr$nbr){ # f = 1
    # l'index de chaque voie :
    indL <- unlist(gr_mcr$list[f])

    # le range d'intensite par voie :
    y_zoom <- lapply(indL, function(ii,mat) range(mat[[ii]][-1,]), mat = res_als$CList) %>% range()

    # graphe des concentrations par voie :
    tiff(filename = paste0("Figures/",fol,"/",gr_mcr$name_gr,"_",gr_mcr$grp[f],"_",nc,"_cp.tiff"), width = 1000, height = 580)
    par(mar = c(5,5,5,0.1), cex.main=2, cex.lab = 2, cex.axis = 2, mgp = c(3.5,1.5,0))

    # les dates :
    xt <- as.numeric(rownames(res_als$CList[[indL[1]]])) %>% as.POSIXct(origin = "1970-01-01", tz = "GMT")

    # premiere acquisition :
    matplot(xt , res_als$CList[[indL[1]]], pch = 16, col = col_nc,
            ylim = y_zoom, main = gr_mcr$grp[f], xlim = x_zoom,
            xlab = "Time (min)", ylab = "Intensity (a.u.)")
    matplot(xt ,res_als$CList[[indL[1]]], lty = 1, type = "l", col = col_nc, add = TRUE)

    # les autres acquisitons :
    for(i in indL[-1]){
      xt <- as.numeric(rownames(res_als$CList[[i]]))
      matplot(xt , res_als$CList[[i]], pch = 16, col = col_nc, add= TRUE)
      matplot(xt , res_als$CList[[i]], lty = 1, type = "l", col = col_nc, add = TRUE)
    }

    legend("topright", legend = paste("Comp", 1:nc), bty = "n", pch = 16, col =  col_nc)
    dev.off()
  }
}

#' a fork of the function 'preprocess2' from the alsace package
#' @param X a conform matrix
#' @param dim1 a vector with the time points
#' @param dim2 a vector with the MS peak
#'
#' @return a conformed matrix for ALS
#' @noRd
preprocess2 <- function (X = Xraw[[1]], dim1 = tpoints, dim2 = lambdas){
  dX <- list(h = 1:nrow(X), v = 1:ncol(X))
  X <- apply(X, 2, function(xx) stats::approx(dX$h, xx, dX$h)$y)
  X <- apply(X, 1, function(xx) stats::approx(dX$v, xx, dX$v)$y) %>% t()
  if (min(X) < 0) X <- X - min(X)
  dimnames(X) <- list(dX$h, dX$v)
  return(X)
}

#' a fork of the function 'opa' from the alsace package. "Finding the most
#' dissimilar variables in a data matrix: the Orthogonal Projection Approach"
#' @param x the tea object
#' @param ncomp integer number of componant
#'
#' @return matrix with the result of the Orthogonal projection approach
#' @noRd
opa2 <- function (x = tea, ncomp = ncMCR){
  if (is.list(x)) x <- do.call("rbind", x)
  lambdas <- colnames(x)
  x <- t(apply(x, 1, function(xx) xx/rep(sqrt(crossprod(xx)), length(xx))))
  Xref <- matrix(0, ncomp, ncol(x))
  huhn <- colMeans(x)
  Xref[1, ] <- huhn/rep(sqrt(crossprod(huhn)), length(huhn))
  for (i in 1:ncomp) {
    Xs <- lapply(1:nrow(x), function(ii, xx, xref) rbind(xref, xx[ii, ]), x, Xref[1:(i - 1), ])
    dissims <- sapply(Xs, function(xx) det(tcrossprod(xx)))
    Xref[i, ] <- x[which.max(dissims), ]
  }
  colnames(Xref) <- lambdas
  t(Xref)
}

#' a fork of the function 'doALS' from the alsace package.
#' "Wrapper function for als, plus some support functions"
#' @param Xl the output of preprocess2
#' @param PureS the output of opa2
#'
#' @return an object with MCR ALS result
#' @noRd
doALS2 <- function (Xl = tea, PureS = tea.opa){

  rd <- sapply(Xl,dim)[1,]
  if(max(rd) != min(rd)){
    dL <- abs(subtract(rd, max(rd)))
    for(i in 1:length(rd)) Xl[[i]] <- rbind(Xl[[i]], matrix(0,dL[i],ncol(Xl[[i]])))
  }

  Cini <- lapply(Xl, function(xl) xl[, 1:ncol(PureS)])

  result <- ALS::als(PsiList = Xl, CList = Cini, S = PureS, normS = 0.5,  optS1st = FALSE, )

  if(max(rd) != min(rd)){
    for(i in 1:length(rd)){
      result$CList[[i]] <- result$CList[[i]][1:rd[i],]
      result$resid[[i]] <- result$resid[[i]][1:rd[i],]
      Xl[[i]] <- Xl[[i]][1:rd[i],]
    }
  }

  colnames(result$S) <- paste("Component", 1:ncol(PureS))

  for (i in 1:length(result$CList)) colnames(result$CList[[i]]) <- colnames(result$S)
  predicted.values <- lapply(1:length(Xl), function(ii) Xl[[ii]] - result$resid[[ii]])

  predvals2 <- sum(unlist(predicted.values)^2)

  npoints <- sapply(rd,multiply_by,e2= nrow(result$S)) %>% sum()

  result$summ.stats <- list(lof = 100 * sqrt(result$rss/predvals2),
                            rms = sqrt(result$rss/npoints), r2 = 1 - result$rss/predvals2)
  class(result) <- "ALS"
  return(result)
}

#' mcr.voc is a large function to easily do an MCR analysis for the dataset.
#' As output, you can retrieve figures and MCR results.
#'
#' @param ncMCR integer; number of componant of MCR
#' @param grp a character string of the group 's column name
#' @param pk_sel a vector of selected peaks, or "all"
#' @param time_format a charater string "date" or "time"
#' @param Li the list with spectra (sp)
#'
#' @return return figures and results
#' @export
#'
#' @examples
#' # mcr.result <- mcr.voc(ncMCR = 3, grp = "modality")
mcr.voc <- function(ncMCR = 3, grp = NULL, pk_sel = "all", time_format = c("date","time"), Li = sp){
  Xraw <- lapply(Li$acq, sort_mat, selec = pk_sel, s_T = time_format[1], L = Li) # liste des spectres en fct du temps trier pour les 4 chambres
  tea <- lapply(Xraw, preprocess2) # preprocess2 is forked from alsace::preprocess
  tea.opa <- opa2(tea, ncMCR) # opa2 is forked from alsace::opa
  tea.als <- doALS2(tea, tea.opa) # this step can be long ; doALS2 is forked from alsace::doALS

  Scomp <- tea.als$S # spectres pures
  Ccomp <- do.call(rbind, tea.als$CList)  # concentration pures
  rownames(Ccomp) <- sp$names_acq[sp$Sacq]

  if(is.null(grp)){
    Li$mt$meta <- cbind(Li$mt$meta,"no_grp" = rownames(Li$mt$meta))
    grp <- "no_grp"
  }

  mcr.print(col_mt = grp, nc = ncMCR,
            Lg = Li, s_T = time_format[1], res_als = tea.als)

  return(list("Scomp" = Scomp, "Ccomp" = Ccomp,
              "ncMCR" = ncMCR, "group" = grp))
}

#### End of Code ####

# library(baseline)
# library(dygraphs)
# library(magrittr)
# library(MALDIquant)
# library(rhdf5)
# library(rmarkdown)
# library(scales)
# library(stringr)
# library(viridis)
# library(xts)
