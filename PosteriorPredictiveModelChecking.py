from pylab import *
import numpy as np

def calcT(T,pred,data,times=None,err=5*np.sqrt(2)):
    if T == 'T0':
        return (pred-data)**2
    elif T == 'T1':
        return (pred-data)**2/data
    elif T == 'T2':
        return (pred-data)**2/err**2
    elif ((T == 'ACF')|(T=='ACFT')):
        tau=pred
        #print tau
        tshift = np.reshape(np.ones(len(tau)),((len(tau),1)))*times+np.reshape(tau,((len(tau),1)))
        corr = np.zeros(len(tau))
        for i in range(0,len(tau)):
            datashift = np.interp(times,tshift[i],data[i],left=0,right=0)
            corr[i] = np.correlate(data[i],datashift)*1./len(datashift[datashift>0])
        return corr

def PPMC(T,posteriorFile='posterior_sample.txt',DataFile='Data/dummy_data.dat',ntrials=10000,doplot=True,plotfile=None,plotfile2=None,figN=1,testplot=False,errstd=1.,truetau=None,errtype=None):
    data = np.loadtxt(DataFile)
    times = data[:,0]
    #derrs = data[:,1]
    #print errstd
    posterior_sample = atleast_2d(loadtxt(posteriorFile))
    taus = posterior_sample[:,1]
    goodtaus = np.where(np.abs(taus)<times[-2]-times[0])[0]
    posterior_sample,taus = posterior_sample[goodtaus],taus[goodtaus]
    post_size = posterior_sample.shape[0]
    post_sample_inds = np.random.randint(0,post_size,ntrials)
    data_template = np.reshape(np.ones(post_size),((post_size,1)))*data[:,1]
    post_data = posterior_sample[:,-data.shape[0]:]
    if errtype == None:
        derrs = np.ones(len(times))*errstd
        post_sample_rnoise = np.random.normal(0,np.reshape(np.ones(ntrials),((ntrials,1)))*derrs)
    elif errtype == 'Poisson':
        derrs = errstd*np.sqrt(post_data)
        post_sample_rnoise = np.random.normal(0,derrs[post_sample_inds])
    #data_template = np.reshape(np.ones(ntrials),((ntrials,1)))*data[:,0]
    pred_data = post_data[post_sample_inds]+post_sample_rnoise
    if T == 'ACF':
        T_arr = calcT(T,taus[post_sample_inds],pred_data[post_sample_inds],times=times)
        Td_arr = calcT(T,taus,data_template,times=times)
    elif T == 'ACFT':
        ttimes = np.arange(1,times[-5]-times[0],1)
        T_arr,Td_arr = np.zeros((len(post_sample_inds),len(ttimes))),np.zeros((len(taus),len(ttimes)))
        for ti in range(0,len(ttimes)):
            T_arr[:,ti] = calcT(T,ttimes[ti]*np.ones(len(post_sample_inds)),pred_data[post_sample_inds],times=times)
            Td_arr[:,ti] = calcT(T,ttimes[ti]*np.ones(len(taus)),data_template,times=times)
        print T_arr,Td_arr
    else:
        T_arr = np.sum(calcT(T,pred_data,post_data[post_sample_inds]),axis=1)
        Td_arr = np.sum(calcT(T,data_template,post_data),axis=1)
    if doplot:
        figure(figN)
        clf()
        if T == 'ACFT':
            for ti in range(0,ntrials):
                plot(T_arr[ti],color='blue')
            for ti in range(0,len(taus)):
                plot(Td_arr[ti],color='red')
            gm = np.argsort(np.abs(taus-np.median(taus)))[0]
            plot(Td_arr[gm],color='cyan')
            xlabel('Time Lag (days)')
            ylabel('ACF')
        else:
            hist(T_arr,50,color='blue')
            hist(Td_arr,10,color='r')
            xlabel('T(d)')
            ylabel('Num. of trials')
        show()
        if plotfile != None: savefig(plotfile)
        if plotfile2 != None:
            figure(figN+1)
            clf()
            if T == 'ACFT':
                Td = Td_arr[gm]
                Td_dum = np.reshape(np.ones(ntrials),((ntrials,1)))*Td
                Tvtau_arr = np.sum(calcT('T0',T_arr,Td_dum),axis=1)
                scatter(np.abs(taus[post_sample_inds]),Tvtau_arr,color='blue',s=3)
                axhline(0,ls='dashed',color='black')
                if truetau != None: scatter(np.array([truetau]),np.array([0]),marker='x',color='red',s=50)
            else:
                scatter(np.abs(taus[post_sample_inds]),T_arr,color='blue')
                scatter(np.abs(taus),Td_arr,color='r')
            ylabel('T(d)')
            xlabel('Time Delay')
            show()
            if plotfile2 != None: savefig(plotfile2)
    if testplot:
        for k in range(0,10):
            figure(8)
            clf()
            plot(times,post_data[k])
            scatter(times,pred_data[k])
            savefig('/home/rumbaugh/Fermi/TimeBombs/plots/PPMC_testplot_%i.png'%k)
