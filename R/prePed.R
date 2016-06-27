

"prePed" <- function(Pedig, keep=NULL, thisBreed=NA, lastNative=NA, addNum=FALSE, I=0){
    if(!is.null(keep)){keep <- setdiff(keep, c(NA))}
    colnames(Pedig)[1:3]<-c("Indiv", "Sire", "Dam")
    rownames(Pedig)<-Pedig$Indiv
    if(!("Sex" %in% colnames(Pedig))){Pedig$Sex<-NA}
    if(!("Breed" %in% colnames(Pedig)) & !is.na(thisBreed)){Pedig$Breed<-thisBreed}
    withBreed <- ("Breed" %in% colnames(Pedig))
    withBorn  <- ("Born"  %in% colnames(Pedig))
    Pedig[,"Indiv"] <- as.character(Pedig[,"Indiv"])
    Pedig[,"Sire"]  <- as.character(Pedig[,"Sire"])
    Pedig[,"Dam"]   <- as.character(Pedig[,"Dam"])
    if(withBreed){Pedig[,"Breed"] <- as.character(Pedig[,"Breed"])}
    
    Pedig[Pedig[, "Sire"] %in% "0", "Sire"] <- NA
    Pedig[Pedig[, "Dam"] %in% "0",   "Dam"] <- NA
    
    ######### Cut Pedigree loops #########
    suppressWarnings(ord<-pedigree::orderPed(Pedig[,1:3]))
    if(sum(ord==-1)>0){
      cat("Pedigree loops were detected. We recommend to correct them manually before\n")
      cat("using prePed(). The parents of the following individuals are set to unknown\n")
      cat("to remove the loops.\n")
      print(Pedig[ord==-1, 2:3])
      cat("\n")
    
      Pedig[ord==-1,"Sire"]<-NA
      Pedig[ord==-1,"Dam"]<-NA
      Pedig[ord==-1,"Breed"]<-"Pedigree Error"
    }
 
    ####### Add imaginary ancestors #######
    if(!is.na(lastNative)){
      rownames(Pedig)<-Pedig[,"Indiv"]
      ID <-  Pedig[is.na(Pedig[,"Sire"]) & !is.na(Pedig[,"Dam"]), "Indiv"]
      Pedig[ID, "Sire"]<- paste("S", ID, sep="")
      ID <-  Pedig[!is.na(Pedig[,"Sire"]) & is.na(Pedig[,"Dam"]), "Indiv"]
      Pedig[ID, "Dam"]<- paste("D", ID, sep="")
    }
    
    #### Add lines for ancestors, sort pedigree ####
    Pedig <- nadiv::prepPed(Pedig)  
    Pedig[,"Indiv"] <- as.character(Pedig[,"Indiv"])
    rownames(Pedig)<- Pedig[,"Indiv"]
    
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    ### code sexes as 1 (males) and 2 (females) ####
    if(sum(!is.na(Pedig[,"Sex"]))>0){
      sexes<-names(table(Pedig[,"Sex"]))
      if(length(sexes)>2){
        cat("Warning: The sex has more than 2 levels. Please correct it.\n")
      }
      if(length(sexes)==1){sexes<-c(sexes, "dummysex")}
      Mval <- Mode(Pedig[Pedig[,"Indiv"] %in% Pedig[,"Sire"], "Sex"])
      Fval <- Mode(Pedig[Pedig[,"Indiv"] %in% Pedig[,"Dam"],  "Sex"])
      if(is.na(Mval)){Mval<-setdiff(sexes, Fval)}
      if(is.na(Fval)){Fval<-setdiff(sexes, Mval)}
      if(!is.na(Mval)&!is.na(Fval)&(Mval!=Fval)){
        MF <- c(1,2)
        names(MF) <- c(Mval, Fval)
        Pedig[,"Sex"] <- MF[as.character(Pedig[,"Sex"])]
      }else{
      cat("Meaning of sex labels cannot be determined from pedigree structure.\n")
      }
    }

    #### determine sexes from pedigree structure ####
    wrongMale   <- Pedig[,"Indiv"] %in% Pedig[, "Sire"] & !(Pedig[,"Sex"] %in% c(1, NA))
    wrongFemale <- Pedig[,"Indiv"] %in% Pedig[, "Dam"]  & !(Pedig[,"Sex"] %in% c(2, NA))
    if(sum(wrongMale)+sum(wrongFemale)>0){
      cat("The sex of the following animals was not compatible with the pedigree, so\n")
      cat("it was modified:\n")
      print(Pedig[wrongMale|wrongFemale, 2:3])
      cat("\n")
    }
    Pedig[Pedig[,"Indiv"] %in% Pedig[, "Sire"],"Sex"] <- 1
    Pedig[Pedig[,"Indiv"] %in% Pedig[, "Dam"], "Sex"] <- 2

    ####             Estimate missing breeds              ###
    ###    Animals with missing breeds are assumed to     ###
    ### be from the same breed as most of their offspring ###
    if(withBreed){
      ID  <- Pedig[is.na(Pedig[,"Breed"]), "Indiv"]
      Tab <- Pedig[Pedig[,"Sire"]%in% ID | Pedig[,"Dam"] %in% ID & !is.na(Pedig[,"Breed"]),c("Sire", "Dam", "Breed")] 
      Tab <- data.frame(ID=c(Tab[,1], Tab[,2]), Breed=c(Tab[,3], Tab[,3]))
      Pedig[ID, "Breed"] <- tapply(as.character(Tab[,2]),list(Tab[,1]),Mode)[ID]
    }
    
    ####   convert breed name of founders   ###
    ####  born after lastNative to unknown  ###
    if(withBreed & withBorn & !is.na(lastNative)){
      isFounder <- is.na(Pedig[,"Sire"]) & is.na(Pedig[,"Dam"]) & (Pedig[,"Breed"]==thisBreed | is.na(Pedig[,"Breed"]))
      notNative <- isFounder & Pedig[,"Born"]>lastNative
      names(notNative)<-Pedig[,"Indiv"]
      ID  <- Pedig[is.na(notNative), "Indiv"]
      Tab <- Pedig[Pedig[,"Sire"]%in% ID | Pedig[,"Dam"] %in% ID & !is.na(Pedig[,"Born"]),c("Sire", "Dam", "Born")] 
      Tab <- data.frame(ID=c(Tab[,1], Tab[,2]), Born=c(Tab[,3], Tab[,3]))
      geb <- tapply(Tab[,2],list(Tab[,1]),min)-I
      notNative[ID] <- geb[ID]>lastNative
      Pedig[is.na(notNative) | notNative, "Breed"] <- "unknown"
    }
    
    
    ######  prune Pedigree   #######
    if(!is.null(keep)){
      Pedig<-nadiv::prunePed(Pedig, phenotyped=keep)
      }
    ######    Add numeric IDs    #######
    if(addNum){
      nP<-nadiv::numPed(Pedig[,1:3])
      nP[nP==-998]<-0
      Pedig$numIndiv <- nP[,1]
      Pedig$numSire <- nP[,2]
      Pedig$numDam <- nP[,3]
    }
    
    cols <- c("Indiv", "Sire", "Dam", "Sex")
    if(withBreed){cols <- c(cols, "Breed")}
    if(withBorn){cols <- c(cols, "Born")}
    if(addNum){cols <- c(cols, c("numIndiv", "numSire", "numDam"))}
    cols  <- c(cols, setdiff(colnames(Pedig),  cols))
    Pedig <- Pedig[, cols]

    Pedig
}