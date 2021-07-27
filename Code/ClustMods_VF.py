import numpy as np
import pylab
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import warnings

warnings.filterwarnings("ignore")
"""
Function used to find magnitude of extinction A for different filter bands. It does so by finding y in different filters where y = 1/x - 1.82 and x is the wavelength of the filter. Then it finds a(x) and b(x) using relations from Cardelli et al. Then it uses the relation A/Av = a(x) + b(x)/R to find A where R = 3.1
"""
def dered(x, Av):
    y = (1 / x) - 1.82
    a = 1 + (0.17699 * y) - (0.50447 * y**2) - (0.02427 * y**3) + (0.72085 * y**4) + (0.01979 * y**5) - (0.77530 * y**6) + (0.32999 * y**7)
    b = (1.41338 * y) + (2.28305 * y**2) + (1.07233 * y**3) - (5.38434 * y**4) - (0.62251 * y**5) + (5.30260 * y**6) - (2.09002 * y**7)
    A = Av * (a + (b / 3.1))
    return A

# Function used to plot graphs from data points
def graph(y, z, xlab, ylab, title, color):
    plt.scatter(y, z, c=color)
    plt.gca().invert_yaxis()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)

def rdf(Asep, length, maxbin):
    dr = length / maxbin
    bins = np.linspace(0, length, maxbin)
    hist = np.array([])
    dig = np.digitize(Asep, bins)
    for i in range(int(maxbin)):
        hist = np.hstack((hist, len(Asep[dig==i])))
        if i > 0:
            hist[i] = hist[i] /((4 *np.pi) * ((dr * i)**2 - (dr * (i-1))**2 ))
    hist = hist / np.sum(hist)
    return bins, hist

def imf(arr, maxbin):
    abs_mag = arr - (5 * np.log10(2300/10))
    mag_poly = np.array([-10.24, -9.99, -9.31, -8.79, -8.33, -7.72, -7.04, -5.76, -4.64, -3.45, -2.55, -2.0, -1.51, -0.89, -0.19, 0.42, 0.89, 1.21, 1.44, 1.88])
    mass_poly = np.array([90., 60., 37., 30., 23., 23.3, 17.5, 14.2, 10.9, 7.6, 5.9, 5.2, 4.5, 3.8, 3.35, 2.9, 2.72, 2.54, 2.36, 2.0])
    poly = np.polyfit(mag_poly, mass_poly, 15)
    line = np.poly1d(poly)
    mass = line(abs_mag)
    bins = np.logspace(np.log10(2.5), np.log10(14.5), maxbin)
    hist = np.array([])
    dig = np.digitize(mass, bins)
    for i in range(maxbin):
        hist = np.hstack((hist, len(mass[dig==i])))
    return bins, hist



# Importing data from spreadsheet and appending it to variables
data = np.loadtxt("Pleides_edit.csv", delimiter=",", skiprows=1)
B_V, U_B, V = data.T

#Transformations from B-V, U-B, V arrays to G-R, U-G, G arrays. Transformations from karaali 2005 and bilir 2005
g_r = 1.023 * B_V + 0.016 * U_B
g_r_2 = 1.124 * B_V - 0.252
g = V + 0.634 * B_V - 0.108
u_g = 0.779 * U_B + 0.755 * B_V + 0.801

# Setting arrays equal to variables
x, y, z, q =  g_r, g_r_2, g, u_g
x_copy, q_copy, z_copy, g_r_copy, g_copy = x[:], q[:], z[:], x[:], z[:]

# Creates a series of bins and creates an array showing which bins each star belong in
bins = np.linspace(-0.1, 2, num = 40)
digitized = np.digitize(x, bins)
#######
for j in range(4):
    xlist, ylist, zlist, diglist = (np.array([]) for i in range(4))
    for i in range(len(bins)):
        mask_out = [abs(q_copy[digitized==i] - np.nanmean(q_copy[digitized==i])) < np.nanstd(q_copy[digitized==i])]
        ylist = np.hstack((ylist, q_copy[digitized==i][mask_out]))
        xlist = np.hstack((xlist, x_copy[digitized==i][mask_out]))
        zlist = np.hstack((zlist, z_copy[digitized==i][mask_out]))
        diglist = np.hstack((diglist, digitized[digitized==i][mask_out]))
    q_copy, x_copy, z_copy, digitized = ylist, xlist, zlist, diglist

#outlier elim. for better polynomial fit of HR diagram
bins2 = np.linspace(-0.1, 2, num = 40)
digitized2 = np.digitize(x, bins2)
for k in range(1):
    xlist2, ylist2, diglist2 = (np.array([]) for i in range(3))
    for i in range(len(bins2)):
        mask_out = [abs(g_copy[digitized2==i] - np.nanmean(g_copy[digitized2==i]))< np.nanstd(g_copy[digitized2==i])]
        xlist2 = np.hstack((xlist2, g_r_copy[digitized2==i][mask_out]))
        ylist2 = np.hstack((ylist2, g_copy[digitized2==i][mask_out]))
        diglist2 = np.hstack((diglist, digitized2[digitized2==i][mask_out]))
    g_r_copy, g_copy, digitized2 = xlist2, ylist2, diglist2

# Creating polynomial line of best fit for the g-r against u-g scatter plot
l = np.polyfit(x_copy, q_copy, 15)
poly_line = np.poly1d(l)

#Calculating x and y points for line of best fit
x_lit = np.linspace(-0.3, 1.5, 50)
y_lit = poly_line(x_lit)
#plt.plot(x_lit, y_lit)
print("hello")
