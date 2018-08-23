import matplotlib.pyplot as plt

try:
    star_number=int(raw_input('Input(0 to 9):'))
except ValueError:
    print "Not a number"

fo = open('rminmax.txt', 'r')  ## data file

datas = fo.readlines()

seconds = [] ## time
star = [] ## radiis of star

i = 0
t = 0
while(i < len(datas)):
	
	seconds.append(t)
	t += 0.004  # will be determined by main code
	star.append(float(datas[i + star_number]))
	i += 10
	

lines = plt.plot(seconds,star)
plt.setp(lines, color='r', linewidth=0.9)
#plt.yscale('log')
plt.title("Radii of Star vs. Time")
plt.xlabel('T')
plt.ylabel('R')
#plt.axis([0, 2000, 0, 3]) ## graph interval
plt.grid(True)
plt.show()

fo.close();
