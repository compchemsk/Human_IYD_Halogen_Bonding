for i in {1..100}
do
cd Pull-${i}
echo "
parm hIYD-I-Tyr.parm7
trajin pull-${i}.nc 
distance D1 @7305 @1898 out d1.dat
distance D2 @7303 @1437 out d2.dat
distance D3 @7306 @1503 out d3.dat
distance D4 @7303 @7270 out d4.dat
distance D5 @7312 @7267 out d5.dat
distance D6 @7324 @7262 out d6.dat
run
quit " > cp.in
cpptraj cp.in
paste d1.dat d2.dat d3.dat d4.dat d5.dat d6.dat > dist.dat
awk '{print($1*0.02),$2,$4,$6,$8,$10,$12}' dist.dat > seq-event.dat
sed -i '1d' seq-event.dat
rm d1.dat d2.dat d3.dat d4.dat d5.dat d6.dat dist.dat
echo "
import numpy as np
import matplotlib.pyplot as plt
data=np.loadtxt('seq-event.dat')
fig = plt.figure()
fig.set_size_inches(10, 10)
plt.plot(data[:,0],data[:,1],color='black',label='K182')
plt.plot(data[:,0],data[:,2],color='red',label='E157')
plt.plot(data[:,0],data[:,3],color='green',label='Y161')
plt.plot(data[:,0],data[:,4],color='blue',label='Fl1')
plt.plot(data[:,0],data[:,5],color='orange',label='Fl2')
plt.plot(data[:,0],data[:,6],color='magenta',label='phenolate')
plt.xlabel('Time (ns)',fontsize=20)
plt.ylabel('Distance (\u212B)',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim(0,27)
plt.ylim(0,50)
plt.legend(loc='upper right',frameon=False)
plt.show() " > pl.py
python pl.py
cd ../
done
