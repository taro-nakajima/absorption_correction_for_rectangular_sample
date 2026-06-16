import math
import sys

args=sys.argv

hkl=[]
FH1 = open(args[1],"r")
for line in FH1:
	h_temp, k_temp, l_temp = line.split()
	hkl.append([float(h_temp),float(k_temp),float(l_temp)])

FH3 = open("result.txt","w")

N2=20	# mesh for positions in the sample.

a_star=1.263		# A-1
c_star=0.1729		# A-1 (reciprocal axis perpendicular to the plate)

Ei=34.05
ki=math.sqrt(Ei/2.072)		

sigma_abs=0.058+0.070	# absorption cross section. calculated by https://www.ncnr.nist.gov/resources/activation/

L=0.35      # width of the sample (unit: cm)
t=0.04     # thickness of the sample (unit: cm)

FH3.write("# Absorption correction for a rectangular shaped sample \n")
FH3.write("# Incident energy : {0} \n".format(Ei))
FH3.write("# Absorption cross section : {0} \n".format(sigma_abs))
FH3.write("# Sample shape : L= {0} cm, t= {1}\n".format(L,t))

FH3.write("# H  K  L  AbsFactor  meanPath  IncAngle  SctAngle \n")

for HKLindex in range(len(hkl)):

	Qx=c_star*hkl[HKLindex][2]
	Qy=a_star*hkl[HKLindex][0]
	Q_len=math.sqrt(Qx**2.0+Qy**2.0)

	alpha=math.atan2(Qx,Qy)/math.pi*180.0		# angle between the scattering vector and c* axis. 

	TwoTheta=(math.asin(Q_len/(2.0*ki))*2.0)/math.pi*180.0		# degrees

	Omega=(TwoTheta/2.0-alpha)
	
	IncAngle=Omega
	SctAngle=IncAngle+180.0-abs(TwoTheta)
	if(IncAngle<0):
		IncAngle+=360.0
	if(SctAngle<0):
		SctAngle+=360.0

	tanIncAngle=math.tan(IncAngle/180.0*math.pi)
	tanSctAngle=math.tan(SctAngle/180.0*math.pi)

	AbsFactor=0.0
	meanPath=0.0
	for p in range(N2):
		for q in range(N2):
			dx=t/float(N2)
			dy=L/float(N2)
			x=t/2.0-dx*(float(p)+0.5)
			y=L/2.0-dy*(float(q)+0.5)
			path=0.0
			if(IncAngle>=0.0) and (IncAngle<90.0):
				path+=math.sqrt(min((t/2.0-x)**2.0+((t/2.0-x)*tanIncAngle)**2.0,((L/2.0-y)/tanIncAngle)**2.0+((L/2.0-y))**2.0 ))
			elif(IncAngle>=90.0) and (IncAngle<180.0):
				path+=math.sqrt(min((-t/2.0-x)**2.0+((-t/2.0-x)*tanIncAngle)**2.0,((L/2.0-y)/tanIncAngle)**2.0+((L/2.0-y))**2.0 ))
			elif(IncAngle>=180.0) and (IncAngle<270.0):
				path+=math.sqrt(min((-t/2.0-x)**2.0+((-t/2.0-x)*tanIncAngle)**2.0,((-L/2.0-y)/tanIncAngle)**2.0+((-L/2.0-y))**2.0 ))
			elif(IncAngle>=270.0) and (IncAngle<360.0):
				path+=math.sqrt(min((t/2.0-x)**2.0+((t/2.0-x)*tanIncAngle)**2.0,((-L/2.0-y)/tanIncAngle)**2.0+((-L/2.0-y))**2.0 ))

			if(SctAngle>=0.0) and (SctAngle<90.0):
				path+=math.sqrt(min((t/2.0-x)**2.0+((t/2.0-x)*tanSctAngle)**2.0,((L/2.0-y)/tanSctAngle)**2.0+((L/2.0-y))**2.0 ))
			elif(SctAngle>=90.0) and (SctAngle<180.0):
				path+=math.sqrt(min((-t/2.0-x)**2.0+((-t/2.0-x)*tanSctAngle)**2.0,((L/2.0-y)/tanSctAngle)**2.0+((L/2.0-y))**2.0 ))
			elif(SctAngle>=180.0) and (SctAngle<270.0):
				path+=math.sqrt(min((-t/2.0-x)**2.0+((-t/2.0-x)*tanSctAngle)**2.0,((-L/2.0-y)/tanSctAngle)**2.0+((-L/2.0-y))**2.0 ))
			elif(SctAngle>=270.0) and (SctAngle<360.0):
				path+=math.sqrt(min((t/2.0-x)**2.0+((t/2.0-x)*tanSctAngle)**2.0,((-L/2.0-y)/tanSctAngle)**2.0+((-L/2.0-y))**2.0 ))

			meanPath+=path
			
			AbsFactor+=math.exp(-sigma_abs*path)

	AbsFactor=AbsFactor/float(N2*N2)
	meanPath=meanPath/float(N2*N2)

	FH3.write("{0}  {1}  {2}  {3}  {4}  {5}  {6}\n".format(hkl[HKLindex][0],hkl[HKLindex][1],hkl[HKLindex][2],AbsFactor,meanPath,IncAngle,SctAngle))

FH3.close()
