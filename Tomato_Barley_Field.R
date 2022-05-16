## Mean probing duration (Nonviruliferous BLH)
NPro_B=939/86400
NPro_T=1039/86400
NPro=NPro_B+NPro_T

## Mean probing duration (Viruliferous BLH)
VPro_B=473/86400
VPro_T=1295/86400
VPro=VPro_B+VPro_T

## Relative proging preference (Nonviruliferous BLH)
NBT_T=0.4
NBT_B=0.6

## Relative probing preference (Viruliferous BLH)
VBT_T=0.53
VBT_B=0.47

## Preference for infection status of plants
VtoN=VtoV=NtoV=NtoN=0.5

# Transmission rate
prob=0.38

# Tomato composition (10 to 100%)
composition=seq(0.1,1,0.05)

# Initial number of viruliferous BLH
BLH=10

# Initial number of plants
ini_plant=100

# Running time
Time=100

##################################################################################
## Manipulated
F_V_T=F_N_T=NULL
Time_V=F_Time_V=NULL
Final_BT_T=NULL


for(k in composition){
	# Parameter initialize
	Time_V=NULL

	N_T=ini_plant*k
	N_B=ini_plant*(1-k)
	V_T=V_B=0

	for(i in 1:Time){
		V_T[i+1]=V_T[i]+(VBT_T*VtoN*N_T[i]*VPro_T/((VBT_T*VtoN*N_T[i]*VPro_T+VBT_T*VtoV*V_T[i]*VPro_T)+(VBT_B*VtoN*VPro_B*(ini_plant*(1-k))))*prob*BLH)	# Number of infected tomato plants
		N_T[i+1]=N_T[1]-V_T[i+1]	# Number of non-infected tomato plants
	}

	Time_V=sum(V_T<=0.2*ini_plant*k)	# Time to 20% infection

	F_V_T=rbind(F_V_T,V_T)
	F_N_T=rbind(F_N_T,N_T)

	F_Time_V=rbind(F_Time_V,Time_V)

	Final_BT_T=rbind(Final_BT_T,V_T/(k*100))	# Proportion of infected tomato plants

}

###############################################################################
## Innate (without manipulated preference)

# Parameter define
F_V_T2=F_N_T2=NULL
Time_V2=F_Time_V2=NULL
Final_BT_T2=NULL

for(k in composition){

	# Parameter initialize
	Time_V2=NULL
	N_T2=ini_plant*k
	N_B2=ini_plant*(1-k)
	V_T2=V_B2=0

	for(i in 1:Time){
		V_T2[i+1]=V_T2[i]+(NBT_T*NtoN*N_T2[i]*NPro_T/((NBT_T*VtoN*N_T2[i]*NPro_T+NBT_T*VtoV*V_T2[i]*NPro_T)+(NBT_B*NtoN*NPro_B*(ini_plant*(1-k))))*prob*BLH)	# Number of infected tomato plants
		N_T2[i+1]=N_T2[1]-V_T2[i+1]	# Number of non-infected tomato plants
	}

	Time_V2=c(Time_V2,sum(V_T2<=0.2*ini_plant*k))	# Time to 20% infection

	F_V_T2=rbind(F_V_T2,V_T2)
	F_N_T2=rbind(F_N_T2,N_T2)
	F_Time_V2=rbind(F_Time_V2,Time_V2)

	Final_BT_T2=rbind(Final_BT_T2,V_T2/(k*100))	# Proportion of infected tomato plants

}

####################################################################
###### Graph
####################################################################
## 1. BCTV-infected tomato plants(%)

Graph_BT_T=apply(Final_BT_T,2,mean)
Graph_BT_T2=apply(Final_BT_T2,2,mean)

plot(Graph_BT_T*100,xlab="Time",ylab='Infected tomato plants (%)',type='l',lwd=3,cex.axis=2,cex.lab=2,ylim=range(0,100),col=2)
lines(Graph_BT_T2*100,lwd=3,col=2,lty=2)


#########################################
## 2. Time to 20% infected tomato

infection=0.2

Graph_BT_Time=NULL
for(i in 1:length(composition)){
	n=max(which(Final_BT_T[i,]<0.2))
	a=Final_BT_T[i,n+1]-Final_BT_T[i,n]
	b=infection-Final_BT_T[i,n]
	Graph_BT_Time[i]=n+(b/a)
}

Graph_BT_Time2=NULL
for(i in 1:length(composition)){
	n=max(which(Final_BT_T2[i,]<0.2))
	a=Final_BT_T2[i,n+1]-Final_BT_T2[i,n]
	b=infection-Final_BT_T2[i,n]
	Graph_BT_Time2[i]=n+(b/a)
}

plot(composition*100,Graph_BT_Time,ylim=range(0,10),xlab='Tomato composition (%)',ylab='Time to 20% infection of tomato',type='l',lwd=3,cex.axis=2,cex.lab=2,col=2)
lines(composition*100,Graph_BT_Time2,lty=2,lwd=3,col=2)




