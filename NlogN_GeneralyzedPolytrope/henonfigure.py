import matplotlib.pyplot as plt


unstable_n = [0.5 ,0.5, 0.5, 0.5, 0.5, 0.5,0.6,0.6,0.6,0.6,0.6,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.8,0.8,0.8,0.8,0.8,0.8,0.9,0.9,0.9,0.9,1,1,1,1.1,1.2]
unstable_m = [-0.5,-0.4,-0.3,-0.2,-0.1,0,-0.5,-0.4,-0.3,-0.2,-0.1,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.7,-0.6,-0.5,-0.4,-0.7,-0.6,-0.5,-0.7,-0.7]
stable_n = [0.5,0.6,0.6,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,1,1,1,1,1,1,1.1,1.1,1.1,1.1,1.1,1.1,1.2,1.2,1.2,1.2,1.2,1.2]
stable_m = [0.1,0,0.1,0,0.1,-0.1,0,0.1,-0.1,0,0.1,-0.4,-0.3,-0.2,-0.1,0,0.1,-0.4,-0.3,-0.2,-0.1,0,0.1,-0.4,-0.3,-0.2,-0.1,0,0.1]
undetermined_n = [0.9,0.9,1.1,1.1,1.2,1.2]
undetermined_m = [-0.3,-0.2,-0.6,-0.5,-0.6,-0.5]
	

lines1 = plt.plot(unstable_n,unstable_m , 'kx' , label = 'unstable')
lines2 = plt.plot(stable_n , stable_m , 'ko' , label = 'stable')
lines3 = plt.plot( undetermined_n , undetermined_m , 'c*' , label = 'undetermined')
plt.title("Stability of Generalyzed Polytropes")
plt.xlabel('n')
plt.ylabel('m')
plt.axis([0.4, 1.3, -0.8, 0.2]) ## graph interval
#plt.legend()

#plt.grid(True)
plt.show()