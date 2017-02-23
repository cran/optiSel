


"pedBreedComp"<-function(Pedig, thisBreed){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  Indiv<-1; Sire<-2; Dam<-3; Sex<-4; Breed<-5
  for(i in c(Indiv, Sire, Dam, Breed)){Pedig[,i]<-as.character(Pedig[,i])}

  Rassen<-setdiff(names(table(Pedig[,Breed])),c(thisBreed))
  if("unknown"%in% Rassen){Rassen<-c(setdiff(Rassen, c("unknown")),"unknown")}
  
  Pedig[Pedig[,Breed] %in% Rassen, Sire]<- Pedig[Pedig[,Breed] %in% Rassen, Breed]
  Pedig[Pedig[,Breed] %in% Rassen, Dam] <- Pedig[Pedig[,Breed] %in% Rassen, Breed]

  Breeds <- data.frame(Rassen,"Migrant","Migrant",0,"Dummy",stringsAsFactors=FALSE)
  Origin <- data.frame(c("Migrant","0"),"?","?", 0, "Dummy",stringsAsFactors=FALSE)
  colnames(Breeds)<-colnames(Pedig)[1:5]
  colnames(Origin)<-colnames(Pedig)[1:5]
  Groups<-rbind.data.frame(Origin, Breeds)
  Pedig <-rbind.data.frame(Groups, Pedig[,1:5])
  Pedig[is.na(Pedig[,Sire]),Sire]<-"0"
  Pedig[is.na(Pedig[,Dam]),  Dam]<-"0"
  n <- nrow(Groups)
  Cont <- genecont(Pedig[,Indiv], Pedig[,Sire], Pedig[,Dam], NAncestors=n)[- (1:n),3:n]
  Cont <- Cont[,rev(order(colMeans(Cont)))]
  Cont <- data.frame(Indiv=rownames(Cont), native=1-rowSums(Cont), Cont, stringsAsFactors = FALSE)
  rownames(Cont)<-Cont$Indiv
  if(PedigAsDataTable){
    setDT(Cont)
    }
  Cont
}


