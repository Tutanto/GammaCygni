import json
import matplotlib.pyplot as plt

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
    ax[0].set_ylim(1.e-15, 1.e-10)
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
