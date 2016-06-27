

"getNe"<-function(Year, Div,df=3.5,I){
  IEst   <- c(predict(sm.spline(Year, I, df=df), Year, 0))  
	for(j in 2:length(Div)){if(Div[j]>Div[j-1])Div[j]<-Div[j-1]}
	DivEst <- c(predict(sm.spline(Year,Div,df=df),Year,0))
	DivAbl <- c(predict(sm.spline(Year,Div,df=df),Year,1))
	deltaf <- -IEst*DivAbl/DivEst
	
	Ne       <-1/(2*(deltaf))
	names(Ne) <-Year
	Ne
}

