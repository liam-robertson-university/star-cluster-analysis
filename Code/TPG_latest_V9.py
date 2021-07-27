# Import modules
import numpy as np
import pylab
import matplotlib.pyplot as plt
import ClustModsV9 as cm

# Import input file
data = np.loadtxt("Pleides_edit.csv", delimiter=",")

# Initialising B-V, U-B, V lists
B_minus_Vlist = []
U_minus_Blist = []
V_list = [] 

# Adding data points to B-V, U-B, V lists from input file
for i in range(len(data)):
	B_minus_Vlist.append(data[i][0])
	U_minus_Blist.append(data[i][1])
	V_list.append(data[i][2])

# Converting the lists to arrays
B_minus_V = np.array(B_minus_Vlist)
U_minus_B = np.array(U_minus_Blist)
V = np.array(V_list)

# Transformation from B-V, U-B, V arrays to G-R, U-G, G arrays
# Transformation 1 from karaali 2005
g_minus_r = 1.023*B_minus_V + 0.016*U_minus_B

# Transformation 2 from bilir 2005
g_minus_r_2 = 1.124*B_minus_V - 0.252

# g uses transformation 2
g = V + 0.634*B_minus_V - 0.108

# Transformation 1 from karaali 2005
u_minus_g = 0.779*U_minus_B + 0.755*B_minus_V + 0.801

# Setting arrays equal to variables
x = g_minus_r
y = g_minus_r_2
z = g
q = u_minus_g
x_copy = x[:]
q_copy = q[:]

# Creates a series of bins and creates an array showing which bins each star belong in
bins = np.linspace(-0.15, 2, num = 40)
digitized = np.digitize(x, bins)

# Removes stars with y values greater or less than one sigma away from the mean in each bin. Repeats this process multiple times
for j in range(4):
        ylist = np.array([])
        xlist = np.array([])
        diglist = np.array([])
        for i in range(len(bins)):
                ylist = np.hstack((ylist, cm.reject(q_copy[digitized==i], q_copy[digitized==i])))
                xlist = np.hstack((xlist, cm.reject(x_copy[digitized==i], q_copy[digitized==i])))
                diglist = np.hstack((diglist, cm.reject(digitized[digitized==i], q_copy[digitized==i])))
        q_copy = np.copy(ylist)
        x_copy = np.copy(xlist)
        digitized = np.copy(diglist)

# Creating polynomial line of best fit for the x against q scatter plot
l = np.polyfit(x_copy, q_copy, 10)
f = np.poly1d(l)

# Calculating x and y points for line of best fit
x_new = np.linspace(x_copy[0], x_copy[-1], 50)
y_new = f(x_new)
plt.plot(x_new, y_new)

# Dereddening U-G and G-R
Au = cm.dered(0.3543, "u")
Ag = cm.dered(0.477, "Ag")
Ar = cm.dered(0.6231, "Ar")
ug_true = q - (Au - Ag)
gr_true = x - (Ag - Ar)

# Plotting colour colour graphs
cm.graph(x_copy, q_copy, "", "", "", "magenta")
cm.graph(x, z, "g-r", "g", "g vs g-r w/ transformation 1", "blue")
cm.graph(y, z, "g-r", "g", "g vs g-r w/ transformation 2", "blue")
cm.graph(x, q, "g-r", "u-g", "u-g vs g-r w/ transformation 1", "green")
