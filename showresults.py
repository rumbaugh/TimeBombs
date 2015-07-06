def ShowResults(sampleFile="sample.txt",sample_infoFile="sample_info.txt", levelsFile="levels.txt",posteriorFile="posterior_sample.txt",DataFile='Data/dummy_data.dat',weightFile='weights.txt',plotfilebase=None,clearall=True,showfullsample=True):
    import postprocess
    import display
    if clearall:
        for p in range(0,6): 
            figure(p+1)
            clf()
    postprocess.postprocess(sampleFile=sampleFile,sample_infoFile=sample_infoFile, levelsFile=levelsFile,posteriorFile=posteriorFile,weightFile=weightFile)
    display.Display(posteriorFile=posteriorFile,DataFile=DataFile,plotfilebase=plotfilebase,showfullsample=showfullsample)

