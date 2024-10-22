'''
Creating over 20 unique legend colors using matplotlib
'''
import matplotlib.pyplot as plt
import numpy as np

NCOL = 20



# mcoll = []
# cm = plt.get_cmap('gist_rainbow')
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_prop_cycle(color=[cm(1.*i/NCOL) for i in range(NCOL)])
# for i in range(NCOL):
#     lin, = ax.plot(np.arange(10)*(i+1))
#     mcoll.append(lin.get_color())

# # mcoll = np.reshape(mcoll, (16))

    
# # fig.savefig('moreColors.png')

# plt.show()
 
from matplotlib.pyplot import cm
import numpy as np

#variable n below should be number of curves to plot

n = 20
fig = plt.figure()
color = cm.rainbow(np.linspace(0, 1, n))
for i, c in enumerate(color):
   lin, = plt.plot(np.arange(10)*(i+1), c=c)

#or version 2:
fig = plt.figure()
color = iter(cm.rainbow(np.linspace(0, 1, n)))
for i in range(n):
   c = next(color)
   lin, = plt.plot(np.arange(10)*(i+1), c=c)


#or version 3:
fig = plt.figure()
color = cm.rainbow(np.linspace(0, 1, n))
for i in range(n):
   c = color[i]
   lin, = plt.plot(np.arange(10)*(i+1), c=c)


plt.show()
