import numpy
import random
import sys

name1 = sys.argv[1]
name2 = sys.argv[2]

#Zamenjaj '_interfacearea.tsv' s primerno končnico
volumes1 = numpy.loadtxt(name1 + '_interfacearea.tsv',delimiter='\t').transpose()
volumes2 = numpy.loadtxt(name2 + '_interfacearea.tsv',delimiter='\t').transpose()

B=10000
averages = numpy.zeros(B)
stdevs = numpy.zeros(B)
covs = numpy.zeros(B)

for x in range(B):
    #Array mora imeti toliko naključnih števil, kot je neodvisnih simulacij (oz podenot v neodvisnih simulacijah, če gre za lastnosti, odvisne le od ene podenote) - v tem primeru 5 za 5 neodvisnih simulacij dimera,
    #za volumen hidrofobnega žepa pa 10 - za 5 neodvisnih simulacij in vsako podenoto v dimeru posebej
    resampled_simulations_index1 = numpy.array([random.randint(0,4), random.randint(0,4), random.randint(0,4), random.randint(0,4), random.randint(0,4)])
    resampled_simulations_index2 = numpy.array([random.randint(0,4), random.randint(0,4), random.randint(0,4), random.randint(0,4), random.randint(0,4)])

    if x%100==0:
        print("Number", x)
    data1 = numpy.zeros(100000)
    data2 = numpy.zeros(100000)

    for i in range(5):
        for j in range(20000):
            data1[i*1000 + j] = volumes1[resampled_simulations_index1[i]][random.randint(0,19999)]
            data2[i*1000 + j] = volumes2[resampled_simulations_index2[i]][random.randint(0,19999)]
    averages[x] = numpy.mean(data1) - numpy.mean(data2)
    stdevs[x] = numpy.std(data1, ddof=1) - numpy.std(data2, ddof=1)
    covs[x] = (numpy.std(data1, ddof=1) / numpy.mean(data1)) - (numpy.std(data2, ddof=1) / numpy.mean(data2))

#numpy.savetxt(name + "_bootstrap_onelevel_averages.csv", averages, delimiter=",")
#numpy.savetxt(name + "_bootstrap_onelevel_stdevs.csv", stdevs, delimiter=",")
#numpy.savetxt(name + "_bootstrap_onelevel_covs.csv", covs, delimiter=",")

print("\n")
print("Diff averages p-value ", 2*(averages<0).mean()," or ", 2*(averages>0).mean(), " CI 95%: (", numpy.quantile(averages, 0.025), "-", numpy.quantile(averages, 0.975), "). 99.5%: (", numpy.quantile(averages, 0.0025), "-", numpy.quantile(averages, 0.9975), ").")
print("Diff stdev p-value ", 2*(stdevs<0).mean()," ", 2*(stdevs>0).mean(), " CI 95%: (", numpy.quantile(stdevs, 0.025), "-", numpy.quantile(stdevs, 0.975), "). 99.5%: (", numpy.quantile(stdevs, 0.0025), "-", numpy.quantile(stdevs, 0.9975), ").")
print("Diff cov p-value ", 2*(covs<0).mean()," ", 2*(covs>0).mean(), "CI 95%: (", numpy.quantile(covs, 0.025), "-", numpy.quantile(covs, 0.975), "). 99.5%: (", numpy.quantile(covs, 0.0025), "-", numpy.quantile(covs, 0.9975), ").")
