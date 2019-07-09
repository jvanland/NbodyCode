#!/usr/bin/python

from pylab import *
import sys

#if len(sys.argv) != 3:
#    print "usage: ./elements.py bin_num iter_num"
#    sys.exit(2)

#n = int(sys.argv[1])
#m = int(sys.argv[2])    

#filename = "kozai5k_%d_%d.dat" % (n,m)
filename = "elem2.dat"
f1 = open(filename, 'r')
data = [map(float, line.split()) for line in f1] 
time = array([bit[0] for bit in data]) 
abin = array([bit[1] for bit in data]) 
ebin = array([bit[2] for bit in data])
pbin = array([bit[3] for bit in data])
acm = array([bit[4] for bit in data])
ecm = array([bit[5] for bit in data])
inc = array([bit[6] for bit in data])#*(180/3.1415926)
f1.close

inc_max = [141]*len(time)
inc_min = [39]*len(time)
inc_good = [90]*len(time)
cons_1 = sqrt(1-ebin*ebin)
cons_2 = cos((3.1415926/180.0)*inc)
close = min(pbin)
emax = max(ebin)
ctime = 94.872*pbin.argmin()
clabel = "closest approach %f AU" % (close)
print emax

figure()
subplot(411)
plot(time/1000,abin,label="abin (AU)",linewidth=2)
ylabel(r'$a_1(\mathrm{AU})$',fontsize=16)
tick_params(axis='x',labelbottom='off')
#ylim(-2,1.5)
subplot(414)
plot(time/1000, pbin,linewidth=2)
ylabel(r'$q_1(\mathrm{AU})$',fontsize=16)
xlabel(r'$\mathrm{Time} (10^{3}\mathrm{yrs})$',fontsize=16)
#ylim(0,4)
subplot(413)
plot(time/1000, (1-ebin*ebin),label="ebin",linewidth=2)
yscale('log')
ylabel(r'$1-e_{1}^{2}$',fontsize=16)
#xlabel('Time ($10^{3}$yrs)')
tick_params(axis='x',labelbottom='off')
#ylim(0,1)
subplot(412)
plot(time/1000,inc,label="inc (rad)",linewidth=2)
ylabel(r'$i_{tot}(\mathrm{Deg})$',fontsize=16)
tick_params(axis='x',labelbottom='off')
tight_layout()
plot(time/1000, inc_max, 'k--')
plot(time/1000, inc_min, 'k--')
plot(time/1000, inc_good, 'k--')
#plot(ctime,close,'*',label=clabel) 
#title(filename,fontsize=22)
#xticks([0,0.4,0.8,1.2,1.6, 2.0, 2.4],fontsize=20),yticks(fontsize=20)
#ylim(0,5)
#xlim(0,160000)
#legend(loc=6,prop={'size':20})
tight_layout()
savefig("exampleIns.eps")

#figure()
#subplot(211)
#plot(time/1000, cons_1*cons_2)
#subplot(212)
#plot(time/1000, cons_2)
#subplot(122)
#plot(time/100000,abin,label="abin (AU)",linewidth=2)
#plot(time/100000, pbin,label="pericenter (AU)",linewidth=2)
#plot(time/100000, ebin,label="ebin",linewidth=2)
#ylim(0,1.1), xlim(0.82,0.84)
#xticks(fontsize=20),yticks(fontsize=20),xlabel('time (10^5 yrs)',fontsize=20)

#figure()
#plot(time,acm*(1-ecm),label="c.o.m. pericenter")
#plot(time,acm,label="c.o.m. a")
#xlabel('Time (yrs)')
#ylabel('Distance (AU)')
#title('Binary 5k_3')
#legend(loc=3)
#ylim(0,1000)
#xlim(0,160000)

#f2 = open('superelem9995.dat', 'r')
#data2 = [map(float, line.split()) for line in f2]
#acm2 = array([bit[1] for bit in data2])
#pcm = array([bit[4] for bit in  data2])
#f2.close

#figure()
#plot(time,acm2)
#plot(time,pcm)
#xlim(0,10000)

show()
