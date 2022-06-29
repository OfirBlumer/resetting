import numpy as np

def sample(func,params,rate):
    searching = True
    time = 0
    while searching:
        nextPassage = func(*params)
        nextRate = np.random.exponential(1/rate)
        if nextRate>nextPassage:
            searching = False
            time += nextPassage
        else:
            time += nextRate
    return time

samples = np.array([sample(np.random.gamma,[0.25,1e5],10) for i in range(10000)])
print(samples,samples.mean(),samples.std())