for i in {1..100}
do
cd Pull-${i}
echo "
parm hIYD-I-Tyr.parm7
trajin pull-${i}.nc 
vector v0 center @7314,7308,7311,7313,7309,7310 out com-I-Tyr-traj.dat
run
quit " > cp.in
cpptraj cp.in
sed -i '1d' com-I-Tyr-traj.dat 
awk '{print$2,$3,$4}' com-I-Tyr-traj.dat > com-I-Tyr.dat
rm com-I-Tyr-traj.dat
rm cp.in
echo "
import numpy as np
data=np.loadtxt('com-I-Tyr.dat')
list_int_d=[]
for i in range(len(data)):
  x = data[i,0]
  y = data[i,1]
  z = data[i,2]
  d = (89.898041*x - 68.303435*y - 84.832676*z + 3530.988309)/np.sqrt((89.898041**2) + (68.303435**2) + (84.832676**2))
  list_int_d.append(d)

def calculate_angle(p1, p2, p3):
  v1 = np.array(p1) - np.array(p2)
  v2 = np.array(p3) - np.array(p2)
  cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
  angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))  # Clip to handle numerical precision issues
  angle_deg = np.degrees(angle_rad)
  return angle_deg

list_int_ang=[]
for i in range(len(data)):
  x = data[i,0]
  y = data[i,1]
  z = data[i,2]
  p3 = (x, y, z)
  p1 = (26.5691, 46.0421, 32.7075)
  p2 = (33.6849, 54.9862, 33.0468)
  list_int_ang.append(calculate_angle(p1, p2, p3))

with open('internal_coord_kxxx_${i}.dat', 'w') as f:
  for i in range(len(list_int_d)):
    f.write(str(list_int_d[i]) + '\t' + str(list_int_ang[i]) + '\n')

" > int_coord.py
python int_coord.py
rm int_coord.py
cd ..
done
