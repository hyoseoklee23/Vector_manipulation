## Mean probing duration (Nonviruliferous BLH)
NPro_P=1053/86400
NPro_T=1039/86400
NPro=NPro_P+NPro_T

## Mean probing duration (Viruliferous BLH)
VPro_P=1201/86400
VPro_T=1295/86400
VPro=VPro_P+VPro_T

## Relative proging preference (Nonviruliferous BLH)
NPT_T=0.37
NPT_P=0.63

## Relative probing preference (Viruliferous BLH)
VPT_T=0.47
VPT_P=0.53

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

# Parameter define
F_V_T=F_N_T=NULL
Time_V=F_Time_V=NULL
Final_PT_T=NULL


for(k in composition){

	# Parameter initialize
	Time_V=NULL
	N_T=ini_plant*k
	N_P=ini_plant*(1-k)
	V_T=V_P=0

	for(i in 1:Time){
		V_T[i+1]=V_T[i]+(VPT_T*VtoN*N_T[i]*VPro_T/((VPT_T*VtoN*N_T[i]*VPro_T+VPT_T*VtoV*V_T[i]*VPro_T)+(VPT_P*VtoN*N_P[i]*VPro_P+VPT_P*VtoV*V_P[i]*VPro_P))*prob*BLH)	# Number of infected tomato plants
		N_T[i+1]=N_T[1]-V_T[i+1]	# Number of non-infected tomato plants

		V_P[i+1]=V_P[i]+(VPT_P*VtoN*N_P[i]*VPro_P/((VPT_T*VtoN*N_T[i]*VPro_T+VPT_T*VtoV*V_T[i]*VPro_T)+(VPT_P*VtoN*N_P[i]*VPro_P+VPT_P*VtoV*V_P[i]*VPro_P))*prob*BLH)	# Number of infected ribwort plantain
		N_P[i+1]=N_P[1]-V_P[i+1]	# Number of non-infected ribwort plantain
	}

	Time_V=sum(V_T<=0.2*ini_plant*k)	# Time to 20% infection

	F_V_T=rbind(F_V_T,V_T)
	F_N_T=rbind(F_N_T,N_T)
	F_Time_V=rbind(F_Time_V,Time_V)

	Final_PT_T=rbind(Final_PT_T,V_T/(k*100))	# Proportion of infected tomato plants

}

###############################################################################
## Innate (without manipulated preference)

# Parameter define
F_V_T2=F_N_T2=NULL
Time_V2=F_Time_V2=NULL
Final_PT_T2=NULL

for(k in composition){

	# Parameter initialize
	Time_V2=NULL
	N_T2=ini_plant*k
	N_P2=ini_plant*(1-k)
	V_T2=V_P2=0

	for(i in 1:Time){
		V_T2[i+1]=V_T2[i]+(NPT_T*NtoN*N_T2[i]*NPro_T/((NPT_T*NtoN*N_T2[i]*NPro_T+NPT_T*NtoV*V_T2[i]*NPro_T)+(NPT_P*NtoN*N_P2[i]*NPro_P+NPT_P*NtoV*V_P2[i]*NPro_P))*prob*BLH)	# Number of infected tomato plants
		N_T2[i+1]=N_T2[1]-V_T2[i+1]	# Number of non-infected tomato plants

		V_P2[i+1]=V_P2[i]+(NPT_P*NtoN*N_P2[i]*NPro_P/((NPT_T*NtoN*N_T2[i]*NPro_T+NPT_T*NtoV*V_T2[i]*NPro_T)+(NPT_P*NtoN*N_P2[i]*NPro_P+NPT_P*NtoV*V_P2[i]*NPro_P))*prob*BLH)	# Number of infected ribwort plantain
		N_P2[i+1]=N_P2[1]-V_P2[i+1]	# Number of non-infected ribwort plantain
	}
	Time_V2=sum(V_T2<=0.2*ini_plant*k)	# Time to 20% infection

	F_V_T2=rbind(F_V_T2,V_T2)
	F_N_T2=rbind(F_N_T2,N_T2)
	F_Time_V2=rbind(F_Time_V2,Time_V2)

	Final_PT_T2=rbind(Final_PT_T2,V_T2/(k*100))	# Proportion of infected tomato plants
}

####################################################################
###### Graph
####################################################################
## 1. BCTV-infected tomato plants(%)

Graph_PT_T=apply(Final_PT_T,2,mean)
Graph_PT_T2=apply(Final_PT_T2,2,mean)

plot(Graph_PT_T*100,xlab="Time",ylab='Infected tomato plants (%)',type='l',lwd=3,cex.axis=2,cex.lab=2,ylim=range(0,100))
lines(Graph_PT_T2*100,lwd=3,col=1,lty=2)


#########################################
## 2. Time to 20% infected tomato

infection=0.2

Graph_PT_Time=NULL
for(i in 1:length(composition)){
	n=max(which(Final_PT_T[i,]<0.2))
	a=Final_PT_T[i,n+1]-Final_PT_T[i,n]
	b=infection-Final_PT_T[i,n]
	Graph_PT_Time[i]=n+(b/a)
}

Graph_PT_Time2=NULL
for(i in 1:length(composition)){
	n=max(which(Final_PT_T2[i,]<0.2))
	a=Final_PT_T2[i,n+1]-Final_PT_T2[i,n]
	b=infection-Final_PT_T2[i,n]
	Graph_PT_Time2[i]=n+(b/a)
}

plot(composition*100,Graph_PT_Time,ylim=range(0,10),xlab='Tomato composition (%)',ylab='Time to 20% infection of tomato',type='l',lwd=3,cex.axis=2,cex.lab=2)
lines(composition*100,Graph_PT_Time2,lty=2,lwd=3)




