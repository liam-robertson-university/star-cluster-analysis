import numpy as np
import pylab
import scipy as sp
import matplotlib.pyplot as plt
import math
import ClustMods_VF     as cm

# Import data from spreadsheets and appending it to variables
data = np.loadtxt("g_r_u_pos_error_clusterdata.csv", delimiter=",", skiprows=1)
g_r_data = np.loadtxt("g_r_clusterdata.csv", delimiter=",", skiprows=1)
data_quad = np.loadtxt("g_quadrant_clusterdata.csv", delimiter=',', skiprows=1)
g, r, u, RA, DEC, g_err, r_err, u_err = data.T
g_full, r_full = g_r_data.T
g_quad, RA_quad, DEC_quad = data_quad.T

# Unaltered data for g-r against g plot
g_copy = np.copy(g)
r_copy = np.copy(r)
g_r_copy = g_copy - r_copy

# RA and Dec origin are the coordinates at the centre of the cluster in RA and Dec
# Rdist_max is the maximum distance stars can be from the centre of the cluster before they are filtered
RA_origin = 34.775
DEC_origin = 57.15
Rdist_max = 0.7649434474

# Calculating the error in u-g and the angular separation between each star and the centre of the cluster
u_g_error = np.sqrt(np.square(g_err) + np.square(u_err))
Angsep = np.sqrt((np.square((np.absolute(DEC - DEC_origin))) + np.square((15 * np.cos(((DEC + DEC_origin) * 0.5 * math.pi / 180)) * np.absolute(RA - RA_origin)))))
Angsep_quad = np.sqrt((np.square((np.absolute(DEC_quad - DEC_origin))) + np.square((15 * np.cos(((DEC_quad + DEC_origin) * 0.5 * math.pi / 180)) * np.absolute(RA_quad - RA_origin)))))

# Filtering out stars which are a specified distance from the centre of the cluster (Rdist_max)
var_matrix = np.array([g, r, u, u_g_error])
Rdist_matrix = var_matrix[:, Angsep < Rdist_max]
g_quad1 = g_quad[Angsep_quad < (2 * Rdist_max)]

# Filtering out fainter stars in g magnitude
g = Rdist_matrix[0,:]
g_mask_matrix = Rdist_matrix[:, g < 12]
g, r, u, u_g_error = g_mask_matrix
# Creating and implementing a mask to filter out stars which are towards the red end of the colour spectrum
g_r = g - r
u_g = u - g
filter_matrix = np.array([g_r, u_g, u_g_error])
colour_matrix = filter_matrix[:, g_r < 0.3]
g_r, u_g, u_g_error = colour_matrix
colour_matrix = colour_matrix[:, g_r > 0.125]
g_r, u_g, u_g_error = colour_matrix

# Setting up an array of possible Av values and finding dereddening coefficients Au-Ag, Ag-Ar
Av = np.linspace(0.1, 2.0, 100)
Au_Ag = cm.dered(0.3543, Av) - cm.dered(0.477, Av)
Ag_Ar = cm.dered(0.477, Av) - cm.dered(0.6231, Av)

#turns the array with size (37) into matrices of size (100, 37)
matrix_u_g = np.tile(u_g,(len(Au_Ag),1))
matrix_g_r = np.tile(g_r,(len(Ag_Ar),1))
matrix_u_g_error = np.tile(u_g_error,(len(Au_Ag),1))

#transposing the extinction vectors so they can be taken away from the colour matrices
dered_Au_Ag = Au_Ag.reshape((-1, 1))
dered_Ag_Ar = Ag_Ar.reshape((-1, 1))

#subtracting the extinction vectors so they can be taken away from the colour matrices
dered_u_g = np.subtract(matrix_u_g, dered_Au_Ag)
dered_g_r = np.subtract(matrix_g_r, dered_Ag_Ar)

#finding the chi squared and Calculating the Av value from this
u_g_model = cm.poly_line(dered_g_r)
chi_sq_matrix = np.square((dered_u_g - u_g_model)/(matrix_u_g_error))
chi_sum = np.sum(chi_sq_matrix, axis=1)
min_chi_sum = np.amin(chi_sum)
min_Av = Av[chi_sum.argmin()]

print(min_Av, min_chi_sum)

#finding extinction magnitudes in each filter from the Av found above
Ag_true = cm.dered(0.477, min_Av)
Ar_true = cm.dered(0.6231, min_Av)
Au_true = cm.dered(0.3543, min_Av)

#correcting the cluster data to remove reddening
g_full_cor = g_full - Ag_true
r_full_cor = r_full - Ar_true
g_r_cor = g_r -Ag_true + Ar_true
u_g_cor = u_g -Au_true + Ag_true
g_r_full_cor = g_full_cor - r_full_cor

#finding line of best fit for Pleiades HR diagram
lit_g_fit = np.polyfit(cm.g_r_copy, cm.g_copy, 6)
lit_g_fit_line = np.poly1d(lit_g_fit)
g_r_lit = np.linspace(cm.g_r_copy[0], 1.08, 50)
g_lit = lit_g_fit_line(g_r_lit)

#binning the g-r data for the corrected Cluster NGC869 & taking mean to be used in the dist mod calculation
bins = np.linspace(0.12, 0.39, num = 20)
digitized = np.digitize(g_r_full_cor, bins)
bin_mean_y = [np.nanmean(g_full_cor[digitized == i]) for i in range(len(bins))]
bin_stdev_y = [np.nanstd(g_full_cor[digitized == i]) for i in range(len(bins))]

#Calculating the distance modulus
poly_lit_fit_y = lit_g_fit_line(bins)
g_shift = poly_lit_fit_y - bin_mean_y
g_shift_mean = np.mean(g_shift)
print g_shift_mean
g_shift_err = np.std(g_shift) / len(g_shift)
print g_shift_err

print max(g_full_cor)

#rdf and imf codes
rdf_x, rdf_y = cm.rdf(Angsep_quad, 0.8, 20)
imf_x, imf_y = cm.imf(g_quad1, 10)
imf_x = np.log10(imf_x)
imf_y = np.log10(imf_y)
imf_y_err = np.sqrt(imf_y)
imf_x1 = imf_x[imf_x < 1.]
imf_y1 = imf_y[imf_x < 1.]
imf_x2 = imf_x1[imf_x1 > 0.45]
imf_y2 = imf_y1[imf_x1 > 0.45]

power_fit = np.polyfit(imf_x2, imf_y2, 1)
power_fit_line = np.poly1d(power_fit)
gradient = (power_fit_line(imf_x2)[-1] - power_fit_line(imf_x2)[0]) / (imf_x2[-1] - imf_x2[0])

#
# plt.figure(1)
# plt.plot(g_r_full_cor, g_full_cor,  color='b', marker='.', linestyle='None', label='NGC869 Cluster')
# #plt.plot(cm.x, cm.z, color='#DC143C', marker='.', linestyle='None', label='Pleiades Cluster')
# #plt.plot(g_r_lit, g_lit, color='k', marker='.', linestyle='-', linewidth=2, label='Pleiades CLuster Polynomial Fit')
# plt.title("HR diagram for h-persei Cluster", fontsize = 16, y=1.025)
# plt.xlabel("g-r Colour", fontsize = 16)
# plt.ylabel("g Magnitude", fontsize = 16)
# #plt.legend(loc='upper right')
# #plt.annotate(s='', xy=(0.04,7.25), xytext=(0.04,12.5), arrowprops=dict(facecolor='black', arrowstyle='<->'))
# plt.xlim(xmax=1.0)
# plt.gca().invert_yaxis()
# plt.grid(True)


# plt.figure(2)
# plt.errorbar(bins, bin_mean_y, yerr=bin_stdev_y, color='b', marker='.', ms=10, linestyle='None', label='Binned NGC869 Cluster Data')
# plt.plot(g_r_lit, g_lit, color='k', marker='.', linestyle='-', linewidth=2, label='Pleiades CLuster Polynomial Fit')
# plt.title("HR Diagram showing Pleiades Polynomial Fit \n with Binned NGC869 Cluster Data", fontsize = 16, y = 1.)
# plt.xlabel("g-r Colour", fontsize = 16)
# plt.ylabel("g Magnitude", fontsize = 16)
# plt.gca().invert_yaxis()
# plt.grid(True)
# plt.legend()
#
# plt.figure(3)
# plt.plot(Av, chi_sum, color='#7F007F', linestyle='-')
# plt.title("Graph of $\chi^2$ against $A_v$", y=1.0)
# plt.xlabel("$A_v$", fontsize = 20)
# plt.ylabel("$\chi^2$", fontsize = 20)
# plt.ylim(ymin=-1E+8)
# plt.annotate('$A_v$ value where $\chi^2$ in minimised', xy=(1.215, 2971329), xytext=(0.75, -9E+7),fontsize = 16, arrowprops=dict(facecolor='black', shrink=0.1, width=2))
# plt.grid(True)
#
#
plt.figure(4)
#plt.plot(g_r, u_g, color='b', marker='.', linestyle='None', label='NGC869 Cluster')
plt.plot(cm.x, cm.q, color='r', marker='.', linestyle='None', label='Pleiades Cluster')
plt.plot(cm.x_lit, cm.y_lit, color='k', marker='.', linestyle='-', linewidth=2, label='Pleiades Cluster Polynomial Fit')
#plt.plot(g_r_cor, u_g_cor, color='b', marker='.',linestyle='None')
plt.gca().invert_yaxis()
plt.title("Colour-Colour Plot for Pleiades")
plt.xlabel("g-r Colour", fontsize=16)
plt.ylabel("u-g Colour", fontsize=16)
plt.ylim(ymin=3.2, ymax=-0.5)
plt.legend(fontsize=12)
plt.grid(True)
#
# plt.figure(5)
# plt.plot(g_r_cor, u_g_cor, color='b', marker='.', markeredgecolor='k',linestyle='None', label='NGC869 Cluster')
# plt.plot(cm.x, cm.q, color='r', marker='.', linestyle='None', label='Pleiades Cluster')
# plt.plot(cm.x_lit, cm.y_lit, color='k', marker='.', linestyle='-', linewidth=2, label='Pleiades Cluster Polynomial Fit')
# plt.gca().invert_yaxis()
# plt.title("Colour-Colour Plot \n including De-Reddened NGC869 Cluster")
# plt.xlabel("g-r Colour", fontsize=16)
# plt.ylabel("u-g Colour", fontsize=16)
# plt.ylim(ymin=3.2, ymax=-0.5)
# plt.legend(fontsize=12)
# plt.grid(True)
#
#
# plt.figure(6)
# plt.errorbar(bins, g_shift, yerr=bin_stdev_y, ecolor='#8B008B', capsize=5, capthick=2, color='b', linestyle='None', label = 'Distance Modulus from each NGC869 data bin')
# plt.plot([min(bins), max(bins)], [g_shift_mean, g_shift_mean], color = 'r', linestyle = '--', linewidth = 2, label = 'Overall Average Distance Modulus')
# plt.title("Distance Modulus Graph showing plot of $\Delta$g against g-r Colour")
# plt.xlabel("g-r colour", fontsize = 16)
# plt.ylabel("Distance Modulus ($\Delta$g)", fontsize = 16)
# plt.grid(True)
# plt.legend(fontsize=12)
#
# plt.figure(7)
# plt.plot(rdf_x, rdf_y)
# plt.xlabel(r"Angular Distance from Cluster Centre /  $\theta$ ")
# plt.ylabel("Star Number Density")
# plt.title("Radial Distribution Function for h-persei Cluster")
# plt.grid(True)
#
# plt.figure(8)
# plt.errorbar(imf_x, imf_y, yerr = np.log10(imf_y_err))
# plt.plot(imf_x2, power_fit_line(imf_x2), color='r', linestyle='-', linewidth=2)
# plt.xlabel("log(Mass / $M_\odot$)")
# plt.ylabel("log(Mass Function / $\phi (m) \Delta m $")
# plt.title("Initial Mass Function for h-persei Cluster")
# plt.grid(True)

plt.show()
