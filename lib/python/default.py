def SelectAntennas(energy, direction, position, topo): #GeV, meters

    import numpy as np
    from scipy import interpolate

    #adaptation
    energy=energy*1e9 #from GeV to eV
    theta=np.arccos(direction[2])
    phi=np.arctan2(direction[1], direction[0])
    xt=position[0] #position of tau decay
    yt=position[1]
    zt=position[2]
    theta=theta*180/np.pi #degs
    phi=phi*180/np.pi

    gamma=0.47*np.pi/180*np.log10(energy/1e17)+0.9*np.pi/180 #angle of the cone large (ag)
    gamma=0.42*np.pi/180*np.log10(energy/1e17)+0.45*np.pi/180 #angle of the cone small (cons)
    zcmin=14000 #metres
    zcmax=39000*np.log10(energy/1e17)+55000 #large (ag)
    zcmax=27000*np.log10(energy/1e17)+22000 #small (cons)

    deltar=200
    antxmin=-30000
    antxmax=30000
    antymin=-30000
    antymax=30000

    def TopoToShower(x,y,z,xt,yt,zt,theta,phi): #from coordinates in the topography frame to coordinates in the shower frame
        xp=x-xt
        yp=y-yt
        zp=z-zt
        alpha=theta*np.pi/180 #around y=theta
        beta=phi*np.pi/180 #around z=phi
        #(RzRy)-1
        xpp=np.cos(beta)*np.cos(alpha)*xp + np.sin(beta)*np.cos(alpha)*yp - np.sin(alpha)*zp
        ypp=-np.sin(beta)*xp + np.cos(beta)*yp
        zpp=np.cos(beta)*np.sin(alpha)*xp + np.sin(beta)*np.sin(alpha)*yp + np.cos(alpha)*zp
        return xpp,ypp,zpp



    def ShowerToTopo(x,y,z,xt,yt,zt,theta,phi): #from coordinates in the shower frame to coordinates in the topography frame
        alpha=theta*np.pi/180 #around y=theta
        beta=phi*np.pi/180 #around z=phi
        #RzRy
        xp=np.cos(beta)*np.cos(alpha)*x - np.sin(beta)*y + np.cos(beta)*np.sin(alpha)*z
        yp=np.sin(beta)*np.cos(alpha)*x + np.cos(beta)*y + np.sin(beta)*np.sin(alpha)*z
        zp=-np.sin(alpha)*x +  np.cos(alpha)*z
        xpp=xp+xt
        ypp=yp+yt
        zpp=zp+zt
        return xpp,ypp,zpp



    def GetTopography(antxmin,antxmax,antymin,antymax,deltar,topo):
        xbef=np.arange(antxmin,antxmax,deltar)
        ybef=np.arange(antymin,antymax,deltar)
        nx=len(xbef)
        ny=len(ybef)
        z=np.zeros((nx,ny))
        x=np.zeros((nx,ny))
        y=np.zeros((nx,ny))
        for i in range(nx):
            for j in range(ny):
                z[i,j]=topo.ground_altitude(xbef[i],ybef[j])
        for i in range(nx):
            x[i,:]=xbef[i]
        for i in range(ny):
            y[:,i]=ybef[i]
        return x,y,z #2D arrays




    x,y,z=GetTopography(antxmin,antxmax,antymin,antymax,deltar,topo)
    z += 3.

    xpp,ypp,zpp=TopoToShower(x,y,z,xt,yt,zt,theta,phi)

    #test which ground points (antennas) are in the cone
    radius2=xpp**2+ypp**2
    radius2max=((zpp+zcmin)*np.tan(gamma))**2
    testradius=np.sqrt(radius2)-np.sqrt(radius2max)
    indr=np.nonzero(testradius<=0)
    testedgemin=zcmin+zpp
    indemin=np.nonzero(testedgemin<=0)
    testedgemax=zcmax+zpp
    indemax=np.nonzero(testedgemax>=0)
    ind=np.nonzero((testradius<=0) & (testedgemin<=0) & (testedgemax>=0))

    #check if shower crash into a mountain early (before xmax)
    zinshower=np.arange(0,-zcmin-deltar,-deltar)
    xinshower=np.zeros(len(zinshower))
    yinshower=np.zeros(len(zinshower))
    xintopo,yintopo,zintopo=ShowerToTopo(xinshower,yinshower,zinshower,xt,yt,zt,theta,phi)
    indinterp=np.nonzero((x<xintopo.max()+1000) & (x>xintopo.min()-1000) & (y<yintopo.max()+1000) & (y>yintopo.min()-1000))
    realztopo=interpolate.RectBivariateSpline(np.unique(x[indinterp]), np.unique(y[indinterp]), np.reshape(z[indinterp],(len(np.unique(x[indinterp])),len(np.unique(y[indinterp])) ) ) ) (xintopo,yintopo,grid=False)
    if len(np.nonzero(realztopo>zintopo)[0]) >0:
        ind=[]

    #shadowing
    keep=np.zeros(len(x[ind]))+1
    thetaincone=np.arccos( (zpp[ind]+zcmin)/np.sqrt((xpp[ind])**2+(ypp[ind])**2+(zpp[ind]+zcmin)**2) ) #rad, ref cone
    cosphiincone=xpp[ind]/(np.sqrt((xpp[ind])**2+(ypp[ind])**2+(zpp[ind]+zcmin)**2)*np.sin(thetaincone))
    phiincone=np.arccos( np.round(cosphiincone)  ) #rad, ref cone
    for i in range(len(x[ind])):
        deltarcount=np.arange(0,(zpp[ind][i]+zcmin)/np.cos(thetaincone[i]),deltar)
        xinshower=deltarcount*np.sin(thetaincone[i])*np.cos(phiincone[i]) #ref shower
        yinshower=deltarcount*np.sin(thetaincone[i])*np.sin(phiincone[i]) #ref shower
        zinshower=deltarcount*np.cos(thetaincone[i])-zcmin #ref shower
        xintopo,yintopo,zintopo=ShowerToTopo(xinshower,yinshower,zinshower,xt,yt,zt,theta,phi)
        indinterp=np.nonzero((x<xintopo.max()+1000) & (x>xintopo.min()-1000) & (y<yintopo.max()+1000) & (y>yintopo.min()-1000))
        realztopo=interpolate.RectBivariateSpline(np.unique(x[indinterp]), np.unique(y[indinterp]), np.reshape(z[indinterp],(len(np.unique(x[indinterp])),len(np.unique(y[indinterp])) ) ) ) (xintopo,yintopo,grid=False)
        if len(np.nonzero(realztopo>zintopo)[0]) >0:
            keep[i]=-1

    return ind[0][keep==1] * xpp.shape[1] + ind[1][keep==1]
