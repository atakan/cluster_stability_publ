import matplotlib.pyplot as plt

fo = open('r_u_plane.txt', 'r')  ## data file

datas = fo.readlines()

radiis = [] ## radii
u = [] # radial velocity



i = 0
for j in range(1,16): # 15 pieces!
	while(i < len(datas)-1):
		if(float(datas[i]) == 0):
			break
		radiis.append(float(datas[i]))
		u.append(float(datas[i+1]))
		i += 2
	
	fig = plt.figure(j)
	lines = plt.plot(radiis,u)
	plt.setp(lines, color='r' , linestyle='None' , marker='.' , markeredgewidth=0.1)
	#plt.xscale('log')
	plt.axis([0, 10, -1, 1]) ## graph interval
	#plt.title("r-u state plane for t = " + str((j-1)*1.5) + " crossing-time")
	plt.xlabel('R')
	plt.ylabel('U')
	#plt.axis("off")
	#plt.rcParams['figure.facecolor'] = 'black'
	#fig1 = plt.gcf()
	plt.show()
	#fig1.savefig('test1.png', dpi=75)
	#plt.grid(True)
	
	radiis = []
	u = []
	i += 1
	
fo.close();
