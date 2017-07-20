import numpy as np
import cmath
import sys

system = sys.argv[1]
version = sys.argv[2]

coeffxyz_ref = [line.strip() for line in open("%s/base-line-%s/run-coeff-1.xyz" % (system, version), 'r')] 
coeffxyz = [line.strip() for line in open("new/run-coeff-1.xyz", 'r')]

#coeffxyz_ref = [line.strip() for line in open("base_line/run-TRIMER_VACUO/run-coeff-1.xyz", 'r')] 
j=0
h=0
t = 0
t_prev = 0
eps = 10e-7
for i in range(0,len(coeffxyz)):
    line = coeffxyz[i].split()
    if (line[0] == 'i' and i != 0 and t != int(line[2].replace(',','')) ):
        t_prev = t
        t = int(line[2].replace(',',''))
        #print("now t is", t, "and t_prev is", t_prev)
    if (line[0] != 'Psi;' and line[0] != 'i'):
        if (t != t_prev):
            line_ref = coeffxyz_ref[i].split()
            #print("entered new time-step loop")
            index = int(line[0])-1
            #print ("assigning exp_coeffs for site %d", index)
            exp_coeff = complex(float(line[2]), float(line[3]))
            exp_coeff_ref = complex(float(line_ref[2]), float(line_ref[3]))
            diff = exp_coeff - exp_coeff_ref
#            print(diff)
            diff = float(abs(diff))
#            print(diff)
            if (diff > eps):
                print("diff is significant:", diff, "in line", i, ", time step", t)
		#print 'FAIL'
		#exit(2)
            else:
                print("values identical in line", i, "time step", t)
#print 'PASS'
