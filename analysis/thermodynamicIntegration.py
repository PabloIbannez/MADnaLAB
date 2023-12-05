import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import scipy

# Process when ctrl+c is pressed

import signal

stop = False
def signal_handler(sig, frame):
    print('You pressed Ctrl+C!')
    global stop
    stop = True
    sys.exit(0)


signal.signal(signal.SIGINT, signal_handler)

kT = 0.59
units = 'KcalMol_A'

kcal2kj = 4.184

skip = 0.5

tiFileName = "ti.dat"

with open(sys.argv[1], 'r') as f:
    VLMPsession = json.load(f)

# Get the VLMPsession folder fropm sys.argv[1]
simFolder = sys.argv[1].split("/")
simFolder = "/".join(simFolder[:-1])

allTogether = False
if len(sys.argv) > 2:
        title = sys.argv[2]
        allTogether = True

for sim in VLMPsession['simulations']:
    name,s,r,sinfo = sim

    # Check units
    simUnits = sinfo['units'][0]['type']
    if simUnits != units:
        print("Units are not correct, expected: %s, got: %s" % (units, simUnits))
        exit()

    tiDataFilePath = simFolder + "/" + s + "/" + tiFileName
    outputFilePath = simFolder + "/" + s + "/" + "tiProcessed.dat"

    print("Processing simulation: ", name, ", at ", s)

    try:

        lambda2Data   = {}
        currentLambda = -1.0
        with open(tiDataFilePath, 'r') as f:
            for line in f:
                if "#" in line:
                    currentLambda = float(line.split()[1])
                else:
                    lambda2Data.setdefault(currentLambda, []).append(float(line))

        #Apply skip
        for l in lambda2Data.keys():
            lambda2Data[l] = lambda2Data[l][int(len(lambda2Data[l])*skip):]

        lambdas = []
        mean    = []
        std     = []
        err     = []
        for l in lambda2Data.keys():
            #plt.plot(lambda2Data[l], label="lambda = %s" % l)
            #plt.legend()
            #plt.show()
            m = np.mean(lambda2Data[l])
            s = np.std(lambda2Data[l])
            e  = s / np.sqrt(len(lambda2Data[l]))
            lambdas.append(l)
            mean.append(m)
            std.append(s)
            err.append(e)

        #Interpolate
        f = scipy.interpolate.interp1d(lambdas, mean, kind='cubic')

        #xmin = np.min(lambdas)
        #xmax = np.max(lambdas)

        xmin = 0.0
        xmax = 1.0

        x = np.linspace(xmin, xmax, 100)

        with open(outputFilePath, 'w') as fout:
            for i in range(len(lambdas)):
                fout.write("%s %s %s %s\n" % (lambdas[i], mean[i], std[i], err[i]))

        if not allTogether:

            plt.errorbar(lambdas, mean, yerr=err, fmt='.', label=name)
            plt.plot(x, f(x))
            energy =scipy.integrate.quad(f, xmin, xmax)[0]

            plt.title("Free Energy: %s (kcal/mol), %s (kj/mol) %s (kT)" % (round(energy,2),
                                                                           round(energy*kcal2kj,2),
                                                                           round(energy/kT,2)))

            plt.xlabel(r'$\lambda$', fontsize=20)
            plt.ylabel(r'$\left< \frac{dU}{d\lambda} \right>_{\lambda}$', fontsize=20)

            plt.legend()
            plt.show()
        else:

            #Take the 3 characters after "eps"
            name_p = name.split("eps")[1][:3].split("_")
            name_p = float(".".join(name_p))
            name_p = "$\epsilon_{aminoacids} = %s$" % name_p

            energy =scipy.integrate.quad(f, xmin, xmax)[0]
            plt.errorbar(lambdas, mean, yerr=err, fmt='.', label=r"%s: %s (kcal/mol), %s (kj/mol) %s (kT)" % (name_p,
                                                                                                              round(energy,2),
                                                                                                              round(energy*kcal2kj,2),
                                                                                                              round(energy/kT,2)))

            plt.plot(x, f(x))

    except:
        print("Error processing simulation: ", name)
        print("Error: ", sys.exc_info()[0])
        print("Error message: ", sys.exc_info()[1])
        plt.close()
        if stop:
            exit()
        else:
            continue

if allTogether:

    # Increase title, axis and legend size
    plt.rcParams.update({'font.size': 20})

    plt.xlim([0.0, 1.05])
    plt.ylim([ymin, ymax])

    plt.xlabel(r'$\lambda$', fontsize=20)
    plt.ylabel(r'$\left< \frac{dU}{d\lambda} \right>_{\lambda}$', fontsize=20)
    plt.title(title)
    plt.legend()
    plt.show()



