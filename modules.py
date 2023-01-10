from pathlib import Path
from os import chdir
import json
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from regions import RectangleSkyRegion
from astropy.coordinates import SkyCoord
from gammapy.datasets import MapDataset
from gammapy.maps import Map
from gammapy.modeling.models import Models

# function to add to JSON
def write_json(new_data, filename='sample.json'):
    with open(filename, 'r+') as file:
        # First we load existing data into a dict.
        file_data = json.load(file)
        # Join new_data with file_data inside emp_details
        file_data.update(new_data)
        # Sets file's current position at offset.
        file.seek(0)
        # convert back to json.
        json.dump(file_data, file, indent=4)

# Reading the JSON file
def read_json(filename='results_sliding_window.json'):
    with open(filename) as f:
        data = json.load(f)
    # Create a list of lists where to store useful results
    n_val = len(data[list(data.keys())[0]])
    lists = [[] for _ in range(n_val)]
    # Fill the list
    for key, val in data.items():
        lists[0].append(key)
        ls = list(val.values())
        for i in range(1, n_val):
            lists[i].append(ls[i])

    return lists

# Plot the results
def plot_results(lists, spectra=1, spatial=1):
    fig, ax = plt.subplots(len(lists)-2, 1, sharex='col',
                           sharey='row', figsize=(15, 7))
    fig.subplots_adjust(hspace=1, wspace=1)

    ax[0].semilogy(lists[0], lists[2])
    ax[0].set_ylim(1.e-16, 1.e-12)
    if spectra == 1:
        ax[0].set_title("amplitude")
    else:
        ax[0].set_title("const")
    ax[0].set_ylabel("1 / (cm2 s TeV)")

    ax[1].plot(lists[0], lists[3])
    ax[1].set_ylim(-4, 4)
    ax[1].set_title("lat Position")
    ax[1].set_ylabel("deg")

    ax[2].plot(lists[0], lists[4])
    ax[2].set_ylim(0, 5)
    if spatial == 0:
        ax[2].set_title("Sigma parameter")
    else:
        ax[2].set_title("r_0 parameter")
    if spatial != 0 and spatial !=1:
        ax[3].plot(lists[0], lists[5])
        ax[3].set_ylim(0, 1)
        ax[3].set_title("eta")
        ax[4].plot(lists[0], lists[6])
        ax[4].set_ylim(0, 1)
        ax[4].set_title("e")
        ax[5].plot(lists[0], lists[7])
        ax[5].set_ylim(0, 180)
        ax[5].set_title("phi parameter")
    else:
        ax[3].plot(lists[0], lists[5])
        ax[3].set_ylim(0, 1)
        ax[3].set_title("e")
        ax[4].plot(lists[0], lists[6])
        ax[4].set_ylim(0, 180)
        ax[4].set_title("phi parameter")

    # Combine all the operations and display
    plt.show()

# If the fit failed in one "window", create a fake "zero" result for that window
def addMissingElement(list1, list2): #list1 is the one with the missing window, list2 is the one complete (the nullhypothesis)
    for i, item in enumerate(list1[0]): 
        if item not in list2[0]:
            list2[0].insert(i, list1[0][i])
            list2[1].insert(i, [False, 'Optimization failed. Estimated distance to minimum too large.', 0])
            for j in range(2,len(list2)):
                list2[j].insert(i, 0)
    return list2
    

#Calculate TS using total_stat output from fitting
def computeTS(nullhyp, newhyp):
    H0 = [i[2] for i in nullhyp[1]] #store total_stat in list
    H1 = [i[2] for i in newhyp[1]]
    #Compute TS as difference H0 - H1. H0 is null hypothesis
    TS = [element1 - element2 for (element1, element2) in zip(H0, H1)]
    return TS

#Check if the Optimization terminated successfully    
def IsSuccessful(list):
    succes = []
    for i in list[1]:
        if i[1] == "Optimization terminated successfully.":
            succes.append(1)
        else:
            succes.append(0)
    return succes

#Plot the TS values    
def plotTS(nullhyp, newhyp, max=None, min=None):
    s_TS, f_TS, s_lon, f_lon = [], [], [], []
    list1, list2 = read_json(nullhyp), read_json(newhyp)
    new_list = addMissingElement(list1, list2)
    TS = computeTS(list1, new_list)
    succes = IsSuccessful(new_list)
    for i, ele in enumerate(succes): #Two different arrays for the TS where IsSuccessful is True or False
        if ele > 0:
            s_TS.append(TS[i])
            s_lon.append(float(new_list[0][i]))
        else:
            f_TS.append(TS[i])
            f_lon.append(float(new_list[0][i]))
    avg = np.mean(s_TS)
    fig, ax1 = plt.subplots(figsize=(10, 3))
    ax1.scatter(s_lon, s_TS, color="green")
    ax1.scatter(f_lon, f_TS, color="red")
    ax1.axhline(y = avg, color = 'w', linestyle = '-') #The average TS value
    if max or min:
        ax1.set_ylim([min, max])
    ax1.set_ylabel(r'TS', fontweight='bold')
    ax1.set_xlabel(r'Longitude', fontweight='bold')
    fig.tight_layout()
    plt.show()

#Show the residual for a fit
def plotRes(dataset, model):
    cwd = Path.cwd()
    fdatasets = cwd / 'datasets' / dataset
    fmodels = cwd / model
    res = read_json(model+'/results_sliding_window.json')

    for i in res[0]:
        analysis = MapDataset.read(filename= fdatasets / f"{i}-dataset.fits.gz")
        chdir("..")
        models_read = Models.read(fmodels / f"model_{i}.yaml")
        analysis.models = models_read
        analysis.plot_residuals_spatial(method="diff/sqrt(model)", vmin=-0.5, vmax=0.5, add_cbar=True)
        analysis.models.plot_regions(color="white")
        plt.title(f'Position {i}')
        chdir(cwd)
        plt.show()

#Create a movie to show the position of each window
def showWindow(dataset):
    fig = plt.figure(figsize=(15, 5))
    ims = []
    excess = Map.read("excess_map.fits")
    steps = np.arange(68, 79, 0.5)
    if dataset in ["analysis_1", "analysis_2", "analysis_3", "analysis_4"]: 
        width, height = 10*u.deg, 8*u.deg
    else:
        width, height = 10*u.deg, 12*u.deg
    for i in steps:
        ax = excess.smooth(2).plot(add_cbar=True, animated=True)
        ax.set_title(f"l = {i}")
        region=RectangleSkyRegion(SkyCoord(l=f'{i}', b='2', unit='deg', frame='galactic'), width=width, 
                                    height=height)
        pix_reg = region.to_pixel(wcs=ax.wcs)
        pix_reg.plot(ax=ax, facecolor="none", edgecolor="white")
        ims.append([ax])
    plt.close()
    ani = animation.ArtistAnimation(fig, ims, interval = 500)
    ani.save("movie.mp4")