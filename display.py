from pylab import *

def Display(posteriorFile='posterior_sample.txt',DataFile='Data/dummy_data.dat',plotfilebase=None,showfullsample=True):
  data = loadtxt(DataFile)
  posterior_sample = atleast_2d(loadtxt(posteriorFile))

  ion()
  figure(4)
  if showfullsample:
    for i in xrange(0, posterior_sample.shape[0]):
      hold(False)
      plot(data[:,0], data[:,1], 'bo')
      hold(True)
      plot(data[:,0], posterior_sample[i, -data.shape[0]:], 'r-')
      ylim([0, 1.1*data[:,1].max()])
      draw()
  else:
    hold(False)
    plot(data[:,0], data[:,1], 'bo')
    hold(True)
    plot(data[:,0], posterior_sample[-1, -data.shape[0]:], 'r-')
    ylim([0, 1.1*data[:,1].max()])
    draw()

#ioff()
  show()
  if plotfilebase != None:
    savefig('fig.lightcurve_fit.%s.png'%plotfilebase)
  figure(5)
  hist(posterior_sample[:,9], 20)
  xlabel('Number of Bursts')
  show()
  if plotfilebase != None:
    savefig('fig.burst_hist.%s.png'%plotfilebase)

  figure(6)
  pos = posterior_sample[:, 10:110]
  pos = pos[pos != 0.]
  hist(pos, 1000)
  xlabel('Time')
  title('Positions of Bursts')
  ioff()
  show()
  if plotfilebase != None:
    savefig('fig.burst_positions.%s.png'%plotfilebase)

